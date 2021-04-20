use anyhow::{Result};
use std::collections::{HashMap, BTreeMap};

pub fn kmer_to_vec(kmer: &u64, k: &u8) -> Result<Vec<u8>> {
    let mask = 3;
    let mut marker = *kmer;
    let mut bases: Vec<u8> = Vec::with_capacity((*k).into());
    for _ in 0..*k {
        let c: u8 = match marker & mask {
            0 => 65,
            1 => 67,
            2 => 71,
            3 => 84,
            _ => continue,
        };
        // bases.push((marker * mask) as u8);
        bases.push(c);
        marker = marker >> 2;
    }
    bases.reverse();
    Ok(bases)
}

pub fn kmer_to_str(kmer: &u64, k: &u8) -> Result<String> {
    let bases = kmer_to_vec(kmer, k)?;
    let kstr = String::from_utf8(bases)?;
    Ok(kstr)
}

pub fn iter_kmers(seq: &[u8], k: &u8, is_hpc: &bool) -> Vec<(usize, u64)> {
    let mask = (1 << 2 * *k) - 1;

    let mut kmers: Vec<(usize, u64)> = Vec::with_capacity(seq.len());
    let mut marker = 0u64;
    let mut last_c = 4u64;
    let mut last_size = 0usize;
    for i in 0..seq.len() {
        let c = match seq[i] {
            65 => 0, // A
            67 => 1, // C
            71 => 2, // G
            84 => 3, // T
            _ => continue,
        };
        if *is_hpc && last_c == c {
            continue;
        }
        last_size += 1;
        marker = (marker << 2 | c ) & mask;
        last_c = c;
        if i < (*k - 1) as usize || last_size < *k as usize{
            continue;
        }
        kmers.push((i, marker));
    }
    kmers
}

/// Split sequence into kmers.
/// # Arguments
///
/// * `seq` - sequence
/// * `k` - K-mer size
///
/// # Returns
/// * returns all kmers as HashSet and occurence of each kmer as HashMap

pub fn split_kmers(
    seq: &[u8], 
    k: &u8,
    is_hpc: &bool,
) -> (BTreeMap<usize, u64>, HashMap::<u64, Vec<usize>>) {
    let mut kmer_idx = BTreeMap::<usize, u64>::default();
    let mut kmer_occ = HashMap::<u64, Vec<usize>>::default();

    let kmers: Vec<(usize, u64)> = iter_kmers(seq, k, is_hpc);
    for (x, kmer) in &kmers {
        kmer_idx.entry(*x).or_insert(*kmer);
        kmer_occ.entry(*kmer)
            .or_default()
            .push(*x);
    }

    (kmer_idx, kmer_occ)
}

// Editing model
// Replace -> 0, INSERT -> 1, DELETION -> -1

pub static MODELS: [[i8; 2]; 3] = [
    [1, -1], [-1, 1], [0, 0], 
];

// Mbleven algorithm
// From: https://github.com/fujimotos/mbleven/blob/master/mbleven.py

pub fn mbleven(
    kmer1: &u64,
    kmer2: &u64,
    k: &u8,
) -> Result<usize> {
    let mut res = 3usize;
    for model in &MODELS {
        let cost = check_model(&kmer1, &kmer2, &k, &model);
        if cost < res {
            res = cost;
        }
    }
    Ok(res)
}

pub fn check_model(
    kmer1: &u64,
    kmer2: &u64,
    k: &u8,
    model: &[i8; 2],
) -> usize {
    let mask = 3u8;
    let mut idx1 = 0u8;
    let mut idx2 = 0u8;
    let mut cost = 0u8;
    let mut pad = 0u8;

    while idx1 < *k && idx2 < *k {
        let c1 = kmer1 >> 2*(*k-1-idx1) & mask as u64;
        let c2 = kmer2 >> 2*(*k-1-(idx2-pad)) & mask as u64;
        if c1 != c2 {
            cost += 1;
            if 2 < cost {
                return cost as usize;
            }
            let option = model[cost as usize-1];
            match option {
                -1i8 => idx1 += 1,
                1i8 => idx2 += 1,
                0i8 => {
                    idx1 += 1;
                    idx2 += 1;
                    pad = 0;
                },
                _ => continue,
            };
        } else {
            idx1 += 1;
            idx2 += 1;
            pad = 0;
        }

    }
    (cost + *k - idx1 + *k - idx2) as usize
}


/// Hash function from minimap2
/// https://github.com/lh3/minimap2/blob/master/sketch.c

// fn hash64(kmer: u64, mask: &u64) -> u64 {
//     let mut key = kmer;
//     key = (!key).wrapping_add(key<<21) & mask;
//     key = key ^ (key >> 24);
//     key = ((key + (key<<3)) + (key<<8)) & mask;
//     key = key ^ (key >> 14);
//     key = ((key + (key<<2)) + (key<<4)) & mask;
//     key = key ^ (key >> 28);
//     key = (key + (key<<31)) & mask;
//     key
// }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer() {
        let is_hpc = true;
        let seq = "AAGGTCGATCGAAGCTGATCGATCGATCGTGCTACGTGATGATGCTAGCCTGACTGATCGTAGCAGC";
        let _ = iter_kmers(&seq.as_bytes(), &8, &is_hpc);
        // for (pos, kmer) in &kmers {
        //     let _kmer = kmer_to_str(kmer, &8);
        //     println!("{:?}", _kmer);
        // }
        println!("{:?}", split_kmers(&seq.as_bytes(), &8, &is_hpc));
    }

    #[test]
    fn test_mbeleven() {
        let kmer1: u64 = 65028; // "TTTGAACA"
        let kmer2: u64 = 63506; // "TTAGAACA"
        assert_eq!(2usize, mbleven(&kmer1, &kmer2, &8).unwrap());
    }
}

