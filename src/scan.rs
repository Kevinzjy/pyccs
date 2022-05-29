use std::cmp::{min, max};
use std::collections::{HashMap, BTreeMap};

use anyhow::{Result, anyhow};
use bio::alignment::{pairwise, AlignmentOperation};
use rust_spoa::poa_consensus;

use crate::kmer::{mbleven, split_kmers};
use crate::math::{call_peak, editing_distance};
// use crate::spoa::poa_consensus;


pub fn find_concensus(
    seq: &[u8],
) -> Result<(Vec<(usize, usize)>, String)> {
    let seq_len = seq.len();

    if seq_len <= 50 {
        return Err(anyhow!("Reads too short"));
    }

    let cand_param = vec![(11u8, false), (8u8, false), (11u8, true), (8u8, true)];
    for (k, is_hpc) in &cand_param {
        let tmp_chain = match circular_finder(seq, &k, &is_hpc) {
            Ok(x) => x,
            Err(_) => continue,
        };

        let mut seqs: Vec<Vec<u8>>= Vec::with_capacity(tmp_chain.len());
        for (s, e) in &tmp_chain {
            let mut tmp_seq = String::from_utf8(seq[*s..*e].to_vec())?;
            tmp_seq.push_str("\0");
            seqs.push((tmp_seq.into_bytes()).to_vec());
        }
        let consensus = String::from(poa_consensus(&seqs, 2, 10, -4, -8, -2, -24, -1));

        // Check segment similarity
        let tail = seqs[seqs.len()-1].to_vec();
        let ccs_tail = (&consensus[0..min(tail.len(), consensus.len())].as_bytes()).to_vec();
        let dis_tail = editing_distance(&tail, &ccs_tail)? as f32 / tail.len() as f32;

        let mut dis_body = 0f32;
        if seqs.len() == 2 {
            let seq1 = (&seqs[0][0..min(tail.len(), seqs[0].len())]).to_vec();
            dis_body = editing_distance(&seq1, &ccs_tail)? as f32 / tail.len() as f32;
        } else {
            let seq1 = consensus.as_bytes().to_vec();
            for i in 0..seqs.len()-1 {
                let seq2 = seqs[i].to_vec();
                let tmp_dis = editing_distance(&seq1, &seq2)? as f32 / consensus.len() as f32;
                if tmp_dis > dis_body {
                    dis_body = tmp_dis;
                }
            }
        }

        if dis_body > 0.2 || dis_tail > 0.35 {
            return Err(anyhow!("Filter by distance threshold."));
        }

        return Ok((tmp_chain, consensus));
    }
    
    Err(anyhow!("No circular sequence found."))
}

fn circular_finder(
    seq: &[u8],
    k: &u8,
    is_hpc: &bool
) -> Result<Vec<(usize, usize)>> {
    let p_match = 0.85;
    let p_indel = 0.1;
    let d_min = 40;
    let support_min = 2i32;

    // Split into kmers
    let tmp_ret = split_kmers(&seq, &k, &is_hpc);
    let kmer_idx: BTreeMap::<usize, u64> = tmp_ret.0;
    let kmer_occ: HashMap::<u64, Vec<usize>> = tmp_ret.1;

    // Estimate distance
    let d_estimate = estimate_distance(&kmer_occ, &p_indel, &d_min)?;
    let d_mean = d_estimate.0;
    let d_delta = d_estimate.1;
    let d_support = d_estimate.2;

    if d_support < support_min {
        return Err(anyhow!("No enough supportting reads."));
    }

    // Find optimal hits
    let hits = circular_hits(&kmer_idx, &k, &d_mean, &d_delta, &p_match, &is_hpc)?;
    if hits.len() == 0 {
        return Err(anyhow!("No optimal hits."));
    }
    
    let chain = optimal_chains(&hits)?;
    
    let segments = split_sequence(&chain, seq)?;

    if segments[0].0 > 100 || segments[segments.len()-1].1 < seq.len() - 100 {
        return Err(anyhow!("Clip sequence too long."));
    }

    if segments.len() < 2 {
        return Err(anyhow!("No full length sequence."));
    }

    if segments.len() == 2 && segments[segments.len()-1].1 - segments[segments.len()-1].0 < 30 {
        return Err(anyhow!("Less than 30bp overlap."));
    }

    Ok(segments)
}

/// Calculate approximate length of reptitive sequence
/// Return:
///   d_mean, d_delta, d_support

fn estimate_distance(
    kmer_occ: &HashMap::<u64, Vec<usize>>,
    p_indel: &f32,
    d_min: &usize,
) -> Result<(i32, i32, i32)> {
    let mut tuple_dis: Vec<i32> = Vec::new();
    for (_, pos) in kmer_occ {
        if pos.len() == 1 {
            continue;
        }
        for i in 0 .. pos.len() - 1 {
            let x1 = pos[i];
            let x2 = pos[i+1];
            if x2 - x1 >= *d_min {
                tuple_dis.push((x2 - x1) as i32);
            }
        }
    }
    tuple_dis.sort_unstable();
    if tuple_dis.len() <= 2 {
        return Ok((0, 0, 0));
    }

    let mut dis_uniq = tuple_dis.clone();
    dis_uniq.dedup();
    let d_mean: i32;
    if dis_uniq.len() == 1 {
        d_mean = dis_uniq[0];
    } else {
        match call_peak(&tuple_dis, "scott") {
            Ok(x) => d_mean = x,
            Err(_) => return Ok((0, 0, 0)),
        }
    }
    let mut d_support = 0i32;
    for i in &tuple_dis {
        if *i == d_mean {
            d_support += 1i32;
        }
    }
    assert!(d_support > 0i32);

    // Random walk model
    let d_delta = (2.3 * (p_indel * (d_mean as f32)).sqrt()) as i32;

    Ok((d_mean, d_delta, d_support))
}

/// Find circular hits
/// Return:
///   vec![x, y, y-x, distance]

fn circular_hits(
    all_kmers: &BTreeMap<usize, u64>,
    k: &u8, 
    d_mean: &i32,
    d_delta: &i32,
    p_match: &f32,
    is_hpc: &bool,
) -> Result<BTreeMap::<usize, (usize, usize, usize, usize)>> {
    let mut hits = BTreeMap::<usize, (usize, usize, usize, usize)>::default();
    for (pos, seed) in all_kmers {
        let s: usize = (*pos as i32 + *d_mean - *d_delta) as usize;
        let e: usize = (*pos as i32 + *d_mean + *d_delta) as usize;
        let ret = optimal_hit(all_kmers, seed, k, &s, &e, &is_hpc)?;
        let j = ret.0;
        let score = ret.2;
        if score > ((*k as f32) * (1.0 - *p_match)) as usize {
            continue;
        }
        hits.entry(*pos).or_insert((*pos, j, j - *pos, score));
    }
    Ok(hits)
}

/// Find best hit in [s, e] region
/// Return:
///   vec![pos, kmer, distance, shift]

fn optimal_hit(
    all_kmers: &BTreeMap<usize, u64>,
    seed: &u64,
    k: &u8,
    s: &usize,
    e: &usize,
    is_hpc: &bool,
) -> Result<(usize, u64, usize, usize)> {
    let mut hits: Vec<(usize, u64, usize, usize)> = Vec::with_capacity(e - s + 1);
    for i in *s..(*e+1) {
        if !all_kmers.contains_key(&i) {
            continue;
        }
        let tmp = all_kmers.get(&i).unwrap();       
        let tmp_dis = mbleven(seed, tmp, k)?;
        hits.push((i, *tmp, tmp_dis, max(i.wrapping_sub(*s), e.wrapping_sub(i))));
    }

    if hits.len() == 0 {
        return Ok((0, 0, *k as usize, 0));
    }
    if *is_hpc && hits.len() > 1 {
        hits = collapse_hits(&hits)?;
    }

    let mut min_dis = *k as usize;
    let mut min_shift = *e + *s;
    let mut best_hit = hits[0];

    for x in &hits {
        if x.2 < min_dis || (x.2 == min_dis && x.3 < min_shift) {
            best_hit = *x;
            min_dis = x.2;
            min_shift = x.3;
        }
    }
    Ok(best_hit)
}

pub fn collapse_hits(
    hits: &Vec<(usize, u64, usize, usize)>,
) -> Result<Vec<(usize, u64, usize, usize)>> {
    let mut collapsed: Vec<(usize, u64, usize, usize)> = Vec::with_capacity(hits.len());
    for i in 0..hits.len()-1 {
        if hits[i].0 == hits[i+1].0 + 1 && hits[i].1 == hits[i+1].1 {
            continue;
        }
        if i == 0 {
            collapsed.push(hits[i]);
        }
        collapsed.push(hits[i+1]);
    }
    Ok(collapsed)
}

/// Intersect hits into longest chains
/// Return:
///   vec![x, y, distance]

fn optimal_chains(
    hits: &BTreeMap::<usize, (usize, usize, usize, usize)>,
) -> Result<Vec<(usize, usize, usize)>>{
    let mut init_hit: (usize, usize, usize, usize) = (0, 0, 0, 0);
    for (_, j) in hits {
        init_hit = *j;
        break;
    }

    let mut chains: Vec<Vec<(usize, usize, usize)>> = Vec::new();
    chains.push(vec![(init_hit.0, init_hit.1, init_hit.3), ]);
    for (_, j) in hits {
        let tmp_len1 = chains.len() - 1;
        let tmp_len2 = chains[tmp_len1].len() - 1;

        if chains[tmp_len1][tmp_len2].0 < j.0 && j.0 < chains[tmp_len1][tmp_len2].1 && chains[tmp_len1][tmp_len2].1 < j.1 {
            chains[tmp_len1].push((j.0, j.1, j.3));
        } else if chains[tmp_len1][tmp_len2].1 < j.0 {
            chains.push(vec![(j.0, j.1, j.3)]);
        } else {
            continue;
        }
    }

    let mut best_chain = &chains[0];
    let mut best_len = best_chain[best_chain.len().wrapping_sub(1)].1 - best_chain[0].0;
    for x in &chains {
        let x_len = x[x.len().wrapping_sub(1)].1 - x[0].0; 
        if x_len > best_len {
            best_chain = x;
            best_len = x_len;
        }
    }

    Ok(best_chain.to_vec())
}

/// Split reads using best chain
/// Return:
///   segments

fn split_sequence(
    chain: &Vec<(usize, usize, usize)>,
    seq: &[u8],
) -> Result<Vec<(usize, usize)>> {
    let m_score = 10i32;
    let n_score = -4i32;
    let g_score = -8i32;
    let e_score = -2i32;
    let score_func = |a: u8, b: u8| if a == b { m_score} else { n_score };

    // Find first perfect match
    let mut boundaries: Vec<usize> = Vec::with_capacity(chain.len());
    for (s, e, score) in chain {
        if *score == 0 {
            boundaries.push(*s);
            boundaries.push(*e);
            break; 
        }
    }
    if boundaries.len() == 0 {
        return Err(anyhow!("No perfect match for chain"));
    }

    // Init hits
    let mut hits = HashMap::<usize, (usize, usize, usize)>::default();
    let mut hit_pos: Vec<usize> = Vec::with_capacity(chain.len());
    for hit in chain {
        hits.entry(hit.0).or_insert(*hit);
        hit_pos.push(hit.0);
    }

    let final_x = chain[chain.len()-1].1;
    let mut idx = 0usize;
    let mut last_x: usize = boundaries[boundaries.len()-1];
    
    while last_x < final_x {
        if hits.contains_key(&last_x) {
            idx = hit_pos.iter().position(|&x| x == last_x).unwrap();
            last_x = hits.get(&last_x).unwrap().1;
            boundaries.push(last_x);
        } else {
            while hit_pos[idx] < last_x {
                if idx >= hit_pos.len() - 1 || hit_pos[idx+1] >= last_x {
                    break;
                }
                idx += 1;
            }
            if idx >= hit_pos.len() - 1 {
                break;
            }

            let (q_s, r_s, _) = hits.get(&(hit_pos[idx])).unwrap();
            let (q_e, r_e, _) = hits.get(&(hit_pos[idx+1])).unwrap();

            let q_seq: Vec<u8> = seq[*q_s..*q_e+1].to_vec();
            let r_seq: Vec<u8> = seq[*r_s..*r_e+1].to_vec();

            let mut aligner = pairwise::Aligner::with_capacity(
                q_seq.len(), r_seq.len(), g_score, e_score, score_func
            );
            let alignment = aligner.global(&q_seq, &r_seq);

            let mut q_shift: usize = 0;
            let mut r_shift: usize = 0;

            for x in alignment.operations {
                match x {
                    AlignmentOperation::Match => {q_shift += 1; r_shift += 1;},
                    AlignmentOperation::Subst => {q_shift += 1; r_shift += 1;},
                    AlignmentOperation::Ins => {q_shift += 1;},
                    AlignmentOperation::Del => {r_shift += 1;},
                    AlignmentOperation::Xclip(n) => {q_shift += n;},
                    AlignmentOperation::Yclip(n) => {r_shift += n;},
                }
                if q_shift == last_x - q_s {
                    break;
                }
            }
            last_x = r_s + r_shift;
            boundaries.push(last_x);
        }
    }
    
    let mut segments: Vec<(usize, usize)> = Vec::with_capacity(boundaries.len());
    for i in 0..boundaries.len()-1 {
        segments.push((boundaries[i], boundaries[i+1]));
    }
    if segments[segments.len()-1].1 < final_x {
        segments.push((segments[segments.len()-1].1, final_x));
    }
    
    Ok(segments)
}