use std::path::Path;

use docopt::Docopt;
use serde::Deserialize;
use anyhow::{Result, anyhow};

use pyccs::{fasta, fastq, logger};

const PKG_NAME: &'static str = env!("CARGO_PKG_NAME");
const PKG_VERSION: &'static str = env!("CARGO_PKG_VERSION");
const USAGE: &'static str = "
ccs (scan for circular consensus sequence)

Usage:
  ccs [options] -i FASTX -o FASTX -r FASTX
  ccs (-h | --help)
  ccs --version

Options:
  -i FILE       Input file (fasta/fastq)
  -o FILE       Output CCS sequences
  -r FILE       Output raw sequences
  -t INT        threads number [default: 1]
  -h --help     Show this screen
  --version     Show version
";

#[derive(Clone, Debug, Deserialize)]
struct Args {
    flag_i: String,
    flag_o: String,
    flag_r: String,
    flag_t: String,
    flag_version: bool,
    flag_help: bool,
}

fn main() -> Result<()> {
    let args: Args = Docopt::new(USAGE)
        .and_then(|d| d.deserialize())
        .unwrap_or_else(|e| e.exit());

    if args.flag_version {
        print! ("{} {}", PKG_NAME, PKG_VERSION);
        return Ok(());
    }

    logger::info("Start running");

    let in_file = Path::new(&args.flag_i);
    if !in_file.exists() {
        return Err(anyhow!("Input {} does not exist!", &args.flag_i));
    }

    if args.flag_i.ends_with(".fa") || args.flag_i.ends_with(".fasta") {
        let _ = fasta::scan_fasta(&args.flag_i, &args.flag_o, &args.flag_r, &args.flag_t)?;
    } else if args.flag_i.ends_with(".fa.gz") || args.flag_i.ends_with(".fasta.gz") {
        let _ = fasta::scan_fasta_gz(&args.flag_i, &args.flag_o, &args.flag_r, &args.flag_t)?;
    } else if args.flag_i.ends_with(".fq") || args.flag_i.ends_with(".fastq") {
        let _ = fastq::scan_fastq(&args.flag_i, &args.flag_o, &args.flag_r, &args.flag_t)?;
    } else if args.flag_i.ends_with(".fq.gz") || args.flag_i.ends_with(".fastq.gz") {
        let _ = fastq::scan_fastq_gz(&args.flag_i, &args.flag_o, &args.flag_r, &args.flag_t)?;
    } else {
        return Err(anyhow!("Can't recognize file extension of {}, must be either fa(.gz)/fasta(.gz)/fq(.gz)/fastq(.gz)", &args.flag_i));
    }

    Ok(())
}