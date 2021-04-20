// Deal with Nanosim reads id
// Example: 
// ["ENSMUST00000144554", "ENSMUSG00000028546", "4:110251281-110251521", "-", "241_72_unaligned_0_R_0_491_0"]
// 
// To explain the information in the header, we have two examples:
// 
// ```fasta
// >ref|NC-001137|-[chromosome=V]_468529_unaligned_0_F_0_3236_0
// ```
// 
// - All information before the first _ are chromosome information. 
// - 468529 is the start position.
// - unaligned suggesting it should be unaligned to the reference.
// - The first 0 is the sequence index.
// - F represents a forward strand.
// - 0_3236_0 means that sequence length extracted from the reference is 3236 bases.
// 
// ```fasta
// >ref|NC-001143|-[chromosome=XI]_115406_aligned_16565_R_92_12710_2
// ```
// 
// - This is an aligned read coming from chromosome XI at position 115406. 
// - 16565 is the sequence index. 
// - R represents a reverse complement strand.
// - 92_12710_2 means that this read has 92-base head region (cannot be aligned),
//   followed by 12710 bases of middle region, and then 2-base tail region.

pub fn filter_read_id(read_id: &str) -> bool {
    let id_vec: Vec<&str> = read_id.split("|").collect();
    let _tscp_id = id_vec[0];
    let _gene_id = id_vec[1];
    let _circ_id = id_vec[2];
    let _strand = id_vec[3];
    let attr = id_vec[4];
    let attr_vec: Vec<&str> = attr.split("_").collect();
    let circ_len: u32 = attr_vec[0].parse::<u32>().unwrap();
    let _st: u32 = attr_vec[1].parse::<u32>().unwrap();
    let is_aligned = attr_vec[2];
    let _idx = attr_vec[3];
    let _direction = attr_vec[4];
    let head: u32 = attr_vec[5].parse::<u32>().unwrap();
    let middle: u32 = attr_vec[6].parse::<u32>().unwrap();
    let tail: u32 = attr_vec[7].parse::<u32>().unwrap();


    // if read_id != "ENSMUST00000021181|ENSMUSG00000020831|11:70236877-70237625|-|197_1126_aligned_43558_F_52_821_60" {
    //     return false;
    // }

    // aligned to circRNA
    if is_aligned == "unaligned" {
        return false;
    }

    // head and tail region
    if head > 100u32 || tail > 100u32 {
        return false;
    }

    // more than one lap
    if (middle as f32) < (2.0 * circ_len as f32) {
        return false
    }

    return true
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_filter() {
        let read1 = "ENSMUST00000144554|ENSMUSG00000028546|4:110251281-110251521|-|241_72_unaligned_0_R_0_491_0";
        assert_eq!(filter_read_id(&read1), false);

        let read2 = "ENSMUST00000063690|ENSMUSG00000027068|2:69397616-69397779|+|164_1007_aligned_43557_R_29_405_9";
        assert_eq!(filter_read_id(&read2), true);
    }

}