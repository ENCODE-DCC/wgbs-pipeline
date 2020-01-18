use csv::{Reader, ReaderBuilder, StringRecord, Writer, WriterBuilder};
use serde::{Deserialize, Serialize};
use std::error::Error;
use std::fs::File;
use std::io;
use structopt::StructOpt;

fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::from_args();
    let gembs_bed_file = File::open(&args.gembs_bed_file)?;
    let bismark_bed_file = File::create(&args.bismark_bed_outfile)?;
    let gembs_reader = get_reader(&gembs_bed_file);
    let mut bismark_writer = get_writer(&bismark_bed_file);
    process(gembs_reader, &mut bismark_writer)?;
    Ok(())
}

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

fn get_reader<R: io::Read>(rdr: R) -> Reader<R> {
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_reader(rdr);
    // Not happy about this boilerplate, but gemBS header names depend on the sample
    // barcode, so we need to rename them manually (can't fit to one of serde's renaming
    // patterns)
    reader.set_headers(StringRecord::from(vec![
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
    reader
}

fn get_writer<W: io::Write>(wtr: W) -> Writer<W> {
    WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(wtr)
}

fn process<R: io::Read, W: io::Write>(
    mut rdr: Reader<R>,
    wtr: &mut Writer<W>,
) -> Result<(), io::Error> {
    let mut records = rdr.deserialize::<Row>();
    // Need to get rid of the old header row, csv doesn't overwrite them in place.
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
        wtr.serialize(record)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_get_args() {
        let args = Cli::from_iter(vec!["prog", "infile", "outfile"]);
        assert_eq!(args.gembs_bed_file, std::path::PathBuf::from("infile"));
        assert_eq!(
            args.bismark_bed_outfile,
            std::path::PathBuf::from("outfile")
        );
    }

    #[test]
    fn test_get_reader() -> Result<(), Box<dyn Error>> {
        let data = "\
        Contig\tPos0\tPos1\tRef\tsample_E:Call\tsample_E:Flags\tsample_E:Meth\tsample_E:non_conv\tsample_E:conv\tsample_E:support_call\tsample_E:total\n\
        chr1\t10649\t10650\tC\tC\tGQ=11;MQ=22\t1.000\t4\t0\t5\t5\n\
        ";
        let mut rdr = get_reader(data.as_bytes());
        let records = rdr
            .records()
            .collect::<Result<Vec<StringRecord>, csv::Error>>()?;
        assert_eq!(
            records[1],
            vec![
                "chr1",
                "10649",
                "10650",
                "C",
                "C",
                "GQ=11;MQ=22",
                "1.000",
                "4",
                "0",
                "5",
                "5"
            ]
        );
        Ok(())
    }

    #[test]
    fn test_get_writer() -> Result<(), Box<dyn Error>> {
        let mut wtr = get_writer(vec![]);
        wtr.write_record(&["a", "b", "c"])?;
        let data = String::from_utf8(wtr.into_inner()?)?;
        assert_eq!(data, "a\tb\tc\n");
        Ok(())
    }

    #[test]
    fn test_process() -> Result<(), Box<dyn Error>> {
        let data = "\
        Contig\tPos0\tPos1\tRef\tsample_E:Call\tsample_E:Flags\tsample_E:Meth\tsample_E:non_conv\tsample_E:conv\tsample_E:support_call\tsample_E:total\n\
        chr1\t10649\t10650\tC\tC\tGQ=11;MQ=22\t1.000\t4\t0\t5\t5\n\
        chr1\t10650\t10651\tC\tC\tGQ=11;MQ=22\t1.000\t0\t0\t5\t5\n\
        ";
        let rdr = get_reader(data.as_bytes());
        let mut wtr = get_writer(vec![]);
        process(rdr, &mut wtr)?;
        let data = String::from_utf8(wtr.into_inner()?)?;
        assert_eq!(
            data,
            "chr1\t10649\t10650\t1.0\t4\t0\nchr1\t10650\t10651\t0.0\t0\t0\n"
        );
        Ok(())
    }
}
