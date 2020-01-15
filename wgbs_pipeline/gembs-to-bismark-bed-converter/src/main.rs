use csv::{ReaderBuilder, StringRecord, WriterBuilder};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io;
use structopt::StructOpt;

#[derive(StructOpt)]
struct Cli {
    #[structopt(parse(from_os_str))]
    gembs_bed_file: std::path::PathBuf,
    #[structopt(parse(from_os_str))]
    bismark_bed_outfile: std::path::PathBuf,
}

#[derive(Debug, Deserialize, Serialize)]
struct Row {
    contig: String,
    // Positions could be as large as 2.4e8 (chr1 length)
    start: u32,
    stop: u32,
    #[serde(skip_deserializing)]
    methylation_percentage: f64,
    #[serde(skip)]
    reference_call: String,
    #[serde(skip)]
    genotype_call: String,
    #[serde(skip)]
    flags: String,
    #[serde(skip)]
    methylation_estimate: String,
    // We could probably get away with u8 here, but use u16 to be safe
    non_converted_count: u16,
    converted_count: u16,
    #[serde(skip)]
    total_bases_supporting_call: u16,
    // I thought this would be a useful field, but it's more accurate to sum the
    // converted and non-converted counts, since the gemBS filtering will be consistent.
    #[serde(skip)]
    total_bases: u16,
}

fn main() -> io::Result<()> {
    let args = Cli::from_args();
    let gembs_bed_file = File::open(&args.gembs_bed_file)?;
    let bismark_bed_file = File::create(&args.bismark_bed_outfile)?;

    let mut gembs_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(&gembs_bed_file);

    // Not happy about this boilerplate, but gemBS header names depend on the sample
    // barcode, so we need to rename them manually (can't fit to one of serde's renaming
    // patterns)
    gembs_reader.set_headers(StringRecord::from(vec![
        "contig",
        "start",
        "stop",
        "reference_call",
        "genotype_call",
        "flags",
        "methylation_estimate",
        "non_converted_count",
        "converted_count",
        "total_bases_supporting_call",
        "total_bases",
    ]));

    let mut bismark_writer = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(&bismark_bed_file);

    let mut records = gembs_reader.deserialize::<Row>();
    // The first record consists of the old headers, so we need to skip it. Weird quirk
    // the crate API, ideally there should be a way to overwrite them in place.
    records.next();

    for i in records {
        let mut record = i?;
        let total = record.non_converted_count + record.converted_count;
        // gemBS will output either 0 or 1 for its methylation model, but its not
        // immediately obvious why from the gemBS C plugin code. Here we'll consistently
        // output 0 if the total count is 0 to avoid NaNs in output.
        if total == 0 {
            record.methylation_percentage = 0.;
        } else {
            record.methylation_percentage = record.non_converted_count as f64 / total as f64;
        }
        bismark_writer.serialize(record)?;
    }

    Ok(())
}
