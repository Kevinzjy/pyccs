use anyhow::{Result};
use std::collections::{HashMap, BTreeMap};

pub fn kmer_to_vec(kmer: &u64, k: &u8) -> Result<Vec<u8>> {
    let mask = 3;
    let mut marker = *kmer;
    let mut bases: Vec<u8> = Vec::with_capacity((*k).into());
    for i in 0..*k {
        let c: u8 = match marker & mask {
            0 => 65,
            1 => 67,
            2 => 71,
            3 => 84,
            _ => continue,
        };
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

pub fn iter_kmers(seq: &[u8], k: &u8) -> Vec<(usize, u64)> {
    let mask = (1 << 2 * *k) - 1;

    let mut kmers: Vec<(usize, u64)> = Vec::with_capacity(seq.len());
    let mut marker = 0u64;

    for i in 0..seq.len() {
        let c = match seq[i] {
            65 => 0, // A
            67 => 1, // C
            71 => 2, // G
            84 => 3, // T
            _ => continue,
        };
        marker = (marker << 2 | c ) & mask;
        if i < (*k - 1) as usize {
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
) -> (BTreeMap<usize, u64>, HashMap::<u64, Vec<usize>>) {
    let mut kmer_idx = BTreeMap::<usize, u64>::default();
    let mut kmer_occ = HashMap::<u64, Vec<usize>>::default();

    let kmers: Vec<(usize, u64)> = iter_kmers(seq, k);
    for (x, kmer) in &kmers {
        kmer_idx.entry(*x).or_insert(*kmer);
        kmer_occ.entry(*kmer)
            .or_default()
            .push(*x);
    }

    (kmer_idx, kmer_occ)
}


/// Hash function from minimap2
/// https://github.com/lh3/minimap2/blob/master/sketch.c

fn hash64(kmer: u64, mask: &u64) -> u64 {
    let mut key = kmer;
    key = (!key).wrapping_add(key<<21) & mask;
    key = key ^ (key >> 24);
    key = ((key + (key<<3)) + (key<<8)) & mask;
    key = key ^ (key >> 14);
    key = ((key + (key<<2)) + (key<<4)) & mask;
    key = key ^ (key >> 28);
    key = (key + (key<<31)) & mask;
    key
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer() {
        let seq = "AAGTCGATCGAAGCTGATCGATCGATCGTGCTACGTGATGATGCTAGCCTGACTGATCGTAGCAGC";
        let kmers = iter_kmers(&seq.as_bytes(), &8);
        for (pos, kmer) in &kmers {
            let _kmer = kmer_to_str(kmer, &8);
            println!("{:?}", _kmer);
        }
        println!("{:?}", split_kmers(&seq.as_bytes(), &8));
    }
}

