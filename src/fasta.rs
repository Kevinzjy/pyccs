use std::fs::File;
use std::path::Path;
use std::io::{Write, BufWriter, BufReader};

use indicatif::{ProgressBar, ProgressStyle, ProgressDrawTarget};
use anyhow::{Context, Result, anyhow};
use flate2::read::MultiGzDecoder;
use seq_io::fasta::Record;
use seq_io::parallel::read_parallel;

use crate::scan::find_concensus;
use crate::logger::info;

fn check_fasta(
    infile: &str,
    flag_t: &str,
) -> Result<u64> {
    let n_threads = flag_t.parse::<u32>()
        .with_context(|| format!("Unvalid -t parameters"))?;

    // Input file
    let in_file = Path::new(&infile);
    if !in_file.exists() {
        return Err(anyhow!("Input {} does not exist!", &infile));
    }

    // Open input
    let fasta_fn = File::open(in_file).unwrap();
    let buf_reader = BufReader::with_capacity(1<<22, fasta_fn);
    let reader = seq_io::fasta::Reader::new(buf_reader);

    let mut total_n = 0u64;
    read_parallel(reader, n_threads, n_threads as usize, |record_set| {
        let mut chunk_n = 0u64;
        for _ in record_set.into_iter().enumerate() {
            chunk_n += 1;
        }
        chunk_n
    }, |record_sets| {
        while let Some(result) = record_sets.next() {
            let (_, tmp_ret) = result.unwrap();
            total_n += tmp_ret;
        }
    });

    Ok(total_n)
}

fn check_fasta_gz(
    infile: &str,
    flag_t: &str,
) -> Result<u64> {
    let n_threads = flag_t.parse::<u32>()
        .with_context(|| format!("Unvalid -t parameters"))?;

    // Input file
    let in_file = Path::new(&infile);
    if !in_file.exists() {
        return Err(anyhow!("Input {} does not exist!", &infile));
    }

    // Open input
    let fasta_fn = File::open(in_file).unwrap();
    let gz = MultiGzDecoder::new(fasta_fn);
    let buf_reader = BufReader::with_capacity(1<<22, gz);
    let reader = seq_io::fasta::Reader::new(buf_reader);

    let mut total_n = 0u64;
    read_parallel(reader, n_threads, n_threads as usize, |record_set| {
        let mut chunk_n = 0u64;
        for _ in record_set.into_iter().enumerate() {
            chunk_n += 1;
        }
        chunk_n
    }, |record_sets| {
        while let Some(result) = record_sets.next() {
            let (_, tmp_ret) = result.unwrap();
            total_n += tmp_ret;
        }
    });

    Ok(total_n)
}

pub fn scan_fasta(
    infile: &str,
    ccsfile: &str,
    rawfile: &str,
    flag_t: &str,
) -> Result<()> {
    let n_reads = check_fasta(&infile, &flag_t)?;
    info(format!("{} total reads", n_reads));

    let n_threads = flag_t.parse::<u32>()?;
    let in_file = Path::new(&infile);

    // CCS output
    let ccs_f = File::create(ccsfile)
        .with_context(|| format!("Unvalid output CCS file"))?;
    let mut ccs_writer = BufWriter::new(ccs_f);

    // Raw output
    let raw_f = File::create(rawfile)
        .with_context(|| format!("Unvalid output CCS file"))?;
    let mut raw_writer = BufWriter::new(raw_f);

    // Process
    info("Step 1 - Scanning raw reads to find circRNA candidates");
    let fastx_fn = File::open(in_file).unwrap();
    let buf_reader = BufReader::with_capacity(1<<22, fastx_fn);
    let reader = seq_io::fasta::Reader::new(buf_reader);

    let bar = ProgressBar::new(n_reads);
    bar.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] [{bar:40./}] {pos:>7}/{len:7} ({eta})")
        .progress_chars("#>-"));
    bar.set_draw_target(ProgressDrawTarget::stderr());

    let mut total_ccs = 0u64;
    read_parallel(reader, n_threads, n_threads as usize, |record_set| {
        let mut chunk_n = 0u64;
        let mut chunk_ret: Vec<(usize, String, String)> = Vec::with_capacity(100);
        for (i, record) in record_set.into_iter().enumerate() {
            chunk_n += 1;
            let seq = record.seq();
            let _read_id = record.id().unwrap();

            // let is_circ = filter_read_id(&read_id);
            // if !is_circ {
            //     continue;
            // }

            match find_concensus(seq) {
                Err(_) => continue,
                Ok(x) => {
                    let ccs_seq = x.1;
                    let mut segments: Vec<String> = Vec::default();
                    for (s, e) in x.0 {
                        segments.push(String::from(vec![s.to_string(), e.to_string()].join("-")));
                    }
                    let segment_str = segments.join(";");
                    chunk_ret.push((i, segment_str, ccs_seq));
                }
            };
        }
        (chunk_n, chunk_ret)
    }, |record_sets| {
        while let Some(result) = record_sets.next() {
            let (record_set, tmp_ret) = result.unwrap();
            let chunk_n = tmp_ret.0;
            let chunk_ret: Vec<_> = tmp_ret.1;
            let chunk_reads: Vec<_> = record_set.into_iter().collect();
            
            for (i, segment_str, ccs_seq) in chunk_ret {
                write!(
                    raw_writer, 
                    ">{}\n{}\n",
                    String::from_utf8(chunk_reads[i].head().to_vec()).unwrap(),
                    String::from_utf8(chunk_reads[i].seq().to_vec()).unwrap()
                ).unwrap();
                write!(
                    ccs_writer,
                    ">{}\t{}\t{}\n{}\n",
                    chunk_reads[i].id().unwrap(),
                    segment_str,
                    ccs_seq.len(), 
                    ccs_seq
                ).unwrap();
                total_ccs += 1;
            }
            bar.inc(chunk_n);
        };
    });
    bar.finish();

    info(format!("{}/{} circular candidate reads", total_ccs, n_reads));

    Ok(())
}

pub fn scan_fasta_gz(
    infile: &str,
    ccsfile: &str,
    rawfile: &str,
    flag_t: &str,
) -> Result<()> {
    let n_reads = check_fasta_gz(&infile, &flag_t)?;
    info(format!("{} total reads", n_reads));

    let n_threads = flag_t.parse::<u32>()?;
    let in_file = Path::new(&infile);

    // CCS output
    let ccs_f = File::create(ccsfile)
        .with_context(|| format!("Unvalid output CCS file"))?;
    let mut ccs_writer = BufWriter::new(ccs_f);

    // Raw output
    let raw_f = File::create(rawfile)
        .with_context(|| format!("Unvalid output CCS file"))?;
    let mut raw_writer = BufWriter::new(raw_f);

    // Process
    info("Step 1 - Scanning raw reads to find circRNA candidates");
    let fastx_fn = File::open(in_file).unwrap();
    let gz = MultiGzDecoder::new(fastx_fn);
    let buf_reader = BufReader::with_capacity(1<<22, gz);
    let reader = seq_io::fasta::Reader::new(buf_reader);

    let bar = ProgressBar::new(n_reads);
    bar.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] [{bar:40./}] {pos:>7}/{len:7} ({eta})")
        .progress_chars("#>-"));
    bar.set_draw_target(ProgressDrawTarget::stderr());

    let mut total_ccs = 0u64;
    read_parallel(reader, n_threads, n_threads as usize, |record_set| {
        let mut chunk_n = 0u64;
        let mut chunk_ret: Vec<(usize, String, String)> = Vec::with_capacity(100);
        for (i, record) in record_set.into_iter().enumerate() {
            chunk_n += 1;
            let seq = record.seq();
            let _read_id = record.id().unwrap();

            // let is_circ = filter_read_id(&read_id);
            // if !is_circ {
            //     continue;
            // }

            match find_concensus(seq) {
                Err(_) => continue,
                Ok(x) => {
                    let ccs_seq = x.1;
                    let mut segments: Vec<String> = Vec::default();
                    for (s, e) in x.0 {
                        segments.push(String::from(vec![s.to_string(), e.to_string()].join("-")));
                    }
                    let segment_str = segments.join(";");
                    chunk_ret.push((i, segment_str, ccs_seq));
                }
            };
        }
        (chunk_n, chunk_ret)
    }, |record_sets| {
        while let Some(result) = record_sets.next() {
            let (record_set, tmp_ret) = result.unwrap();
            let chunk_n = tmp_ret.0;
            let chunk_ret: Vec<_> = tmp_ret.1;
            let chunk_reads: Vec<_> = record_set.into_iter().collect();
            
            for (i, segment_str, ccs_seq) in chunk_ret {
                write!(
                    raw_writer, 
                    ">{}\n{}\n",
                    String::from_utf8(chunk_reads[i].head().to_vec()).unwrap(),
                    String::from_utf8(chunk_reads[i].seq().to_vec()).unwrap()
                ).unwrap();
                write!(
                    ccs_writer,
                    ">{}\t{}\t{}\n{}\n",
                    chunk_reads[i].id().unwrap(),
                    segment_str,
                    ccs_seq.len(), 
                    ccs_seq
                ).unwrap();
                total_ccs += 1;
            }
            bar.inc(chunk_n);
        };
    });
    bar.finish();

    info(format!("{}/{} circular candidate reads", total_ccs, n_reads));

    Ok(())
}