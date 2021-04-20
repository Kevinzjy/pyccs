// use std::fs::File;
// use std::path::Path;
// use std::io::{BufRead, BufReader};
// use flate2::read::MultiGzDecoder;
// use anyhow::{Context, Result, anyhow, bail};

// // Reading Fastq
// pub fn read_fastx (file_name: &Path) -> Result<dyn BufRead> {
//     if file_name.extension().unwrap() == "gz" {
//         _read_fastx_gz(file_name)
//     } else {
//         _read_fastx(file_name)
//     }
// }

// fn _read_fastx_gz (file_name: &Path) -> Result<dyn BufRead> {
//     let fastq_fn = File::open(file_name).unwrap();
//     let gz = MultiGzDecoder::new(fastq_fn);
//     let buf_reader = BufReader::with_capacity(32*1024, gz);
//     Ok(buf_reader)
// }

// fn _read_fastx (file_name: &Path) -> Result<dyn BufRead> {
//     let fastq_fn = File::open(file_name).unwrap();
//     let buf_reader = BufReader::with_capacity(32*1024, fastq_fn);
//     Ok(buf_reader)
// }
