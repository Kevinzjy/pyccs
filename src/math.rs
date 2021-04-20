
use std::f32::consts::PI;
use std::cmp::{min, max};
use std::collections::HashMap;

use anyhow::{Result, anyhow};

use crate::fastpdf_data::{FAST_MAX, FAST_H, FASTPDFX, FASTPDFY};
use crate::kmer::kmer_to_vec;

pub struct GaussianKDE {
    pub samples: Vec<i32>,
    pub bandwidth: f32,
}

pub trait Density {
    fn density(&self, x: i32) -> f32;
    fn fast_density(&self, x: i32) -> f32;
}

impl Density for GaussianKDE {
    fn density(&self, x: i32) -> f32 {
        let n = self.samples.len() as f32;
        let mut s = 0f32;
        for i in &self.samples {
            if *i == x{
                s += 1.0;
                continue; 
            }
            match *i > x {
                true => {
                    let r = (i - x) as f32 / self.bandwidth;
                    s += (-0.5 * r.powi(2)).exp();
                },
                false => {
                    let r = (x - i) as f32 / self.bandwidth;
                    s += (-0.5 * r.powi(2)).exp();
                }
            }            
        }
        let sqrt_2pi = (2.0 * PI).sqrt();
        s / n / self.bandwidth / sqrt_2pi
    }

    fn fast_density(&self, x: i32) -> f32 {
        let n = self.samples.len() as f32;
        let mut s = 0f32;
        for i in &self.samples {
            if *i == x {
                s += 1.0;
                continue
            }
            match *i > x {
                true => {
                    let r = (i - x) as f32 / self.bandwidth;
                    s += norm_pdf(&r);
                },
                false => {
                    let r = (x - i) as f32 / self.bandwidth;
                    s += norm_pdf(&r);
                }
            }
        }
        let sqrt_2pi = (2.0 * PI).sqrt();
        s / n / self.bandwidth / sqrt_2pi
    }
}

// Scott's rule for choosing bandwidth

pub fn scott_bw(x: &[i32]) -> f32 {
    let n = x.len() as f32;
    let d = 1f32; // 1d-kde
    n.powf(-1.0 / (d + 4 as f32))
}

// Silverman’s rule for choosing bandwidth

pub fn silverman_bw(x: &[i32]) -> f32 {
    let n = x.len() as f32;
    let d = 1f32;
    (n * (d + 2.) / 4.).powf(-1. / (d + 4.))
}

pub fn gaussian_kde(tuple_dis: &[i32], bw: f32) -> Box<dyn Density> {
    let n = tuple_dis.len();
    assert!(n > 0);
    assert!(bw > 0.0);

    Box::new(GaussianKDE {
        samples: tuple_dis.to_vec(),
        bandwidth: bw,
    })
}

pub fn fast_peak(x: &[i32]) -> Result<i32> {
    // Get frequency
    let mut tmp = HashMap::<i32, i32>::default();
    for i in x {
        *tmp.entry(*i).or_insert(0) += 1;
    }
    let mut freq: Vec<_> = tmp.clone().into_iter().collect();
    freq.sort_by(|x,y| x.0.cmp(&y.0));

    let mut init_freq = freq.clone();
    init_freq.sort_by(|x,y| y.1.cmp(&x.1));
    let mut init_peak = init_freq[0].0;

    // Merge intervals
    let mut delta = 5i32;
    loop {
        let mut last_s = 0i32;
        let mut last_e = 0i32;
        let mut last_f = 0i32;
        let mut merged_freq: Vec<_> = Vec::with_capacity(freq.len());
        for (i, j) in &freq {
            if *i > last_e + delta {
                if last_e != 0 {
                    merged_freq.push((last_s, last_e, last_f));
                }
                last_s = *i;
                last_e = *i;
                last_f = *j;
            } else {
                last_e = max(last_e, *i);
                last_f += j;
            }
        }
        merged_freq.push((last_s, last_e, last_f));
        merged_freq.sort_by(|x,y| y.2.cmp(&x.2));

        let tmp_center: i32 = (merged_freq[0].0 + merged_freq[0].1) / 2;
        let mut tmp_peak = merged_freq[0].0;
        let mut tmp_support = 0;
        for (i, j) in &freq {
            if *i < merged_freq[0].0 || *i > merged_freq[0].1 {
                continue;
            }
            if *j > tmp_support || (*j == tmp_support && *i-tmp_center < tmp_peak-tmp_center) {
                tmp_peak = *i;
                tmp_support = *j;
            }
        }

        if (tmp_peak - init_peak).abs() <= 2 {
            init_peak = tmp_peak;
            delta = (2.3 * (0.1 * tmp_peak as f32).sqrt()) as i32;
        } else {
            break;
        }
    }
    Ok(init_peak)
}

pub fn call_peak(x: &[i32], method: &str) -> Result<i32> {
    if method == "fast" {
        return fast_peak(x);
    }

    let bw: f32 = match method {
        "silverman" => silverman_bw(x),
        "scott" => scott_bw(x),
        _ => 0f32,
    };

    if bw == 0f32 {
        return Err(anyhow!("Wrong bin width parameter!"));
    }

    let kernel = gaussian_kde(x, bw);

    let mut x_max = &x[0];
    let mut p_max = 0f32;
    for i in x {
        let pdf = kernel.density(*i);
        if pdf > p_max {
            x_max = i;
            p_max = pdf;
        }
    }
    Ok(*x_max)
}

// Fast estimation of normal pdf
// Inspired from: https://github.com/yixuan/fastncdf

pub fn norm_pdf(x: &f32) -> f32{
    if *x >= FAST_MAX {
        return 0f32;
    }
    let i: f32 = match *x > 0.0 {
        true => *x,
        false => *x * -1.0,
    };

    let x_est = (i * FAST_H) as usize;
    let w = (i - FASTPDFX[x_est]) * FAST_H;
    let y_est = w * FASTPDFY[x_est+1] + (1.0 - w) * FASTPDFY[x_est];
    y_est
}

// Editing distance

// pub fn editing_distance(seq1: &Vec<u8>, seq2: &Vec<u8>) -> Result<usize> {
//     // Create 2-d array
//     let m = seq1.len();
//     let n = seq2.len();
//     let mut dis_raw = vec![0usize; (m+1) * (n+1)];
//     let mut dis_base: Vec<_> = dis_raw.as_mut_slice().chunks_mut(n+1).collect();
//     let dis_mtx: &mut[&mut[_]] = dis_base.as_mut_slice();

//     for i in 0..n+1 {
//         dis_mtx[0][i] = i;
//     }
//     for i in 0..m+1 {
//         dis_mtx[i][0] = i;
//     }

//     for i in 1..m+1 {
//         for j in 1..n+1 {
//             let left = dis_mtx[i-1][j] + 1;
//             let down = dis_mtx[i][j-1] + 1;
//             let mut left_down = dis_mtx[i-1][j-1];
//             if seq1[i-1] != seq2[j-1] {
//                 left_down += 1;
//             }
//             dis_mtx[i][j] = min(left, min(down, left_down));
//         }
//     }

//     Ok(dis_mtx[m][n])
// }

pub fn editing_distance(seq1: &Vec<u8>, seq2: &Vec<u8>) -> Result<usize> {
    let len_a = seq1.len();
    let len_b = seq2.len();
    if len_a < len_b {
        return editing_distance(seq2, seq1);
    }
    // handle special case of 0 length
    if len_a == 0 {
        return Ok(len_b);
    } else if len_b == 0 {
        return Ok(len_a);
    }

    let len_b = len_b + 1;

    let mut pre;
    let mut tmp;
    let mut cur = vec![0; len_b];

    // initialize string b
    for i in 1..len_b {
        cur[i] = i;
    }

    // calculate edit distance
    for i in 0..seq1.len() {
        let ca = seq1[i];
        pre = cur[0];
        cur[0] = i + 1;
        for j in 0..seq2.len() {
            let cb = seq2[j];
            tmp = cur[j + 1];
            cur[j + 1] = min(
                // deletion
                tmp + 1,
                min(
                    // insertion
                    cur[j] + 1,
                    // match or substitution
                    pre + if ca == cb { 0 } else { 1 },
                ),
            );
            pre = tmp;
        }
    }
    Ok(cur[len_b - 1])
}

// Kmer distance

pub fn kmer_distance(kmer1: &u64, kmer2: &u64, k: &u8) -> Result<usize> {
    assert!(*k > 0);
    let kvec1: Vec<u8> = kmer_to_vec(kmer1, k)?;
    let kvec2: Vec<u8> = kmer_to_vec(kmer2, k)?;
    editing_distance(&kvec1, &kvec2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_call_peak() {
        let x: Vec<i32> = [47, 48, 48, 49, 319, 319, 319, 324].to_vec();
        assert_eq!(319i32, call_peak(&x, "scott").unwrap());
        assert_eq!(319i32, call_peak(&x, "silverman").unwrap());
    }

    #[test]
    fn test_density() {
        // let x: Vec<i32> = [47, 48, 49, 52, 60, 300, 310, 319, 319, 321, 324].to_vec();
        let x: Vec<i32> = [78, 78, 78, 78, 78, 78, 81, 81, 81, 81, 81, 82, 82, 82, 82, 82].to_vec();
        let bw: f32 = scott_bw(&x);

        let kernel = gaussian_kde(&x, bw);
        let mut _x_max1 = &x[0];
        let mut p_max1 = 0f32;
        let mut _x_max2 = &x[0];
        let mut p_max2 = 0f32;
        for i in &x {
            let pdf1 = kernel.fast_density(*i);
            let pdf2 = kernel.density(*i);

            // println!("{:?} {:?}", pdf1, pdf2);

            if pdf1 > p_max1 {
                _x_max1 = i;
                p_max1 = pdf1;
            }

            if pdf2 > p_max2 {
                _x_max2 = i;
                p_max2 = pdf2;
            }
        }
        // assert_eq!(x_max1, x_max2);
    }

    #[test]
    fn test_distance() {
        let kmer1: u64 = 65028; // "TTTGAACA"
        let kmer2: u64 = 63506; // "TTGAACAG"
        let dis = kmer_distance(&kmer1, &kmer2, &8).unwrap();
        assert_eq!(2usize, dis);
    }
    
    // #[test]
    // fn test_fast_peak() {
    //     let x: Vec<i32> = [47, 48, 48, 49, 319, 319, 319, 319, 324].to_vec();
    //     assert_eq!(319i32, call_peak(&x, "fast").unwrap());
    // }
}