
use libc::c_char;
use std::ffi::CStr;
use std::str;

extern "C" {
    fn poa(
        seqs: *const *const u8,
        num_seqs: i32,
        alignment_type: i32, // 0 = local, 1 = global, 2 = gapped
        match_score: i32,
        mismatch_score: i32,
        gap_open: i32,
        gap_extend: i32,
        gap2_open: i32,
        gap2_extend: i32,
    ) -> *const c_char;
}

pub fn poa_consensus(
    seqs: &Vec<Vec<u8>>,
    alignment_type: i32,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
    gap2_open: i32,
    gap2_extend: i32,
) -> String {
    let num_seqs = seqs.len() as i32;
    let mut seq_ptrs: Vec<*const u8> = Vec::with_capacity(seqs.len());
    for seq in seqs {
        seq_ptrs.push(seq.as_ptr());
    }

    let c_buf: *const c_char = unsafe {
        poa(
            seq_ptrs.as_ptr(),
            num_seqs,
            // consensus.as_ptr(),
            alignment_type,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            gap2_open,
            gap2_extend,
        )
    };
    let c_str: &CStr = unsafe { CStr::from_ptr(c_buf) };
    let consensus: &str = c_str.to_str().unwrap();

    String::from(consensus)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_poa() {
        let seqs = vec!["ATTGCCCGTT",
        "AATGCCGTT",
        "AATGCCCGAT",
        "AACGCCCGTC",
        "AGTGCTCGTT",
        "AATGCTCGTT"];
        let mut cseqs: Vec<Vec<u8>> = Vec::with_capacity(seqs.len()); 
        for seq in seqs {
            let mut tmp_seq = String::from(seq);
            tmp_seq.push_str("\0");
            cseqs.push((tmp_seq.into_bytes()).to_vec());
        }
        
        let consensus = poa_consensus(&cseqs, 1, 5, -4, -3, -1, -3, -1);

        let expected = "AATGCCCGTT";
        assert_eq!(consensus, expected);
    }

    #[test]
    fn test_dna_consensus() {
        let mut seqs = vec![];

        // generated each string by adding small tweaks to the expected consensus "AATGCCCGTT"
        for seq in ["ATTGCCCGTT\0",
            "AATGCCGTT\0",
            "AATGCCCGAT\0",
            "AACGCCCGTC\0",
            "AGTGCTCGTT\0",
            "AATGCTCGTT\0"].iter() {
            seqs.push((*seq).bytes().map(|x|{x as u8}).collect::<Vec<u8>>());
        }

        let consensus = poa_consensus(&seqs, 1, 5, -4, -3, -1, -3, -1);

        let expected = "AATGCCCGTT";
        assert_eq!(consensus, expected);
    }


    #[test]
    fn test_protein_consensus() {
        let mut seqs = vec![];
        // expect consensus "FNLKPSWDDCQ"
        for seq in ["FNLKESWDDCQ\0".to_string(),
            "FNLKPSWDCQ\0".to_string(),
            "FNLKSPSWDDCQ\0".to_string(),
            "FNLKASWCQ\0".to_string(),
            "FLKPSWDDCQ\0".to_string(),
            "FNLKPSWDADCQ\0".to_string()].iter() {
            seqs.push(seq.chars().into_iter().map(|x|{x as u8}).collect::<Vec<u8>>());
        }

        let consensus = poa_consensus(&seqs, 1, 5, -4, -3, -1, -3, -1);
        eprintln!("{:?}", &consensus);

        let expected = "FNLKPSWDDCQ";
        assert_eq!(consensus, expected);

    }

    // #[test]
    // // #[should_panic]
    // fn test_not_null_terminated() {
    //     let mut seqs = vec![];

    //     // generated each string by adding small tweaks to the expected consensus "AATGCCCGTT"
    //     for seq in ["ATTGCCCGTT",
    //         "AATGCCGTT",
    //         "AATGCCCGAT",
    //         "AACGCCCGTC",
    //         "AGTGCTCGTT",
    //         "AATGCTCGTT"].iter() {
    //         seqs.push((*seq).bytes().map(|x|{x as u8}).collect::<Vec<u8>>());
    //     }

    //     poa_consensus(&seqs, 1, 5, -4, -3, -1, -3, -1);
    // }
}
