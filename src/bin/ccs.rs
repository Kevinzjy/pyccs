use docopt::Docopt;
use serde::Deserialize;
use anyhow::Result;

use circtools::{scan, logger};

const PKG_NAME: &'static str = env!("CARGO_PKG_NAME");
const PKG_VERSION: &'static str = env!("CARGO_PKG_VERSION");
const USAGE: &'static str = "
ccs (scan for circular consensus sequence)

Usage:
  circtools [options] -i FASTX -o FASTX -r FASTX
  circtools (-h | --help)
  circtools --version

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

    let _ = scan::scan_ccs_reads(
        &args.flag_i,
        &args.flag_o,
        &args.flag_r,
        &args.flag_t
    )?;

    Ok(())
}