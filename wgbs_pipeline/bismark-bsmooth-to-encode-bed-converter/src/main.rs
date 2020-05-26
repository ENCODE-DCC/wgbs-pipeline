use csv::{Reader, ReaderBuilder, Writer, WriterBuilder};
use palette::rgb::Rgb;
use palette::Hsv;
use serde::{Deserialize, Serialize};
use std::error::Error;
use std::fs::File;
use std::io;
use std::{cmp, fmt};
use structopt::StructOpt;

fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::from_args();

    let bismark_bed_file = File::open(&args.bismark_bed_file)?;
    let smoothed_methylation = File::open(&args.smoothed_methylation_tsv)?;
    let encode_bed_file = File::create(&args.encode_bed_outfile)?;

    let bismark_reader = get_reader(&bismark_bed_file);
    let smoothed_methylation_reader = get_reader(&smoothed_methylation);
    let mut encode_writer = get_writer(&encode_bed_file);

    process(
        bismark_reader,
        smoothed_methylation_reader,
        &mut encode_writer,
    )?;

    Ok(())
}

// Bsmooth will sometimes output invalid floats (NA in the tsv), for these rows we skip
// outputting them in the bed
fn process<R: io::Read, W: io::Write>(
    mut bismark_reader: Reader<R>,
    mut smoothed_methylation_reader: Reader<R>,
    encode_writer: &mut Writer<W>,
) -> Result<(), io::Error> {
    let records = bismark_reader.deserialize::<Row>();
    let smoothed_methylation_records =
        smoothed_methylation_reader.deserialize::<SmoothedMethylationRow>();

    for (i, smoothed_methylation_record) in records.zip(smoothed_methylation_records) {
        let mut record = i?;
        match smoothed_methylation_record?.smoothed_methylation_percentage {
            SmoothedMethylationPercentage::Valid(valid) => {
                record.smoothed_methylation_percentage = valid * 100.
            }
            SmoothedMethylationPercentage::Nan(_) => continue,
        };
        record.item_name = ".".to_string();
        record.coverage = record.converted_count + record.non_converted_count;
        record.score = cmp::min(record.coverage, ENCODE_SCORE_CAP);
        record.strandedness = UNKNOWN_STRANDEDNESS;
        record.start_thick_display = record.start;
        record.stop_thick_display = record.stop;
        record.color_value =
            rgb_from_methylation(record.smoothed_methylation_percentage).to_string();
        encode_writer.serialize(record)?;
    }
    Ok(())
}

// We interpolate in HSV color space so that at 50% between red and green we obtain
// the yellow RGB value (255, 255, 0), interpolating in RGB space results in a dark
// yellow (127, 127, 0) instead. Green = 0% methylation, red = 100% methylation
fn rgb_from_methylation(methylation: f64) -> RgbWithDisplay {
    let rgb: Rgb = Hsv::new((1. - (methylation / 100.) as f32) * ENCODE_HSV_MAX_HUE, 1., 1.).into();
    RgbWithDisplay::from(rgb.into_format::<u8>().into_components())
}

fn get_reader<R: io::Read>(rdr: R) -> Reader<R> {
    ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(rdr)
}

fn get_writer<W: io::Write>(wtr: W) -> Writer<W> {
    WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(wtr)
}

const ENCODE_SCORE_CAP: u16 = 1000;
const ENCODE_HSV_MAX_HUE: f32 = 120.;
const UNKNOWN_STRANDEDNESS: char = '.';

#[derive(StructOpt)]
struct Cli {
    #[structopt(parse(from_os_str))]
    bismark_bed_file: std::path::PathBuf,
    #[structopt(parse(from_os_str))]
    smoothed_methylation_tsv: std::path::PathBuf,
    #[structopt(parse(from_os_str))]
    encode_bed_outfile: std::path::PathBuf,
}

struct RgbWithDisplay {
    red: u8,
    green: u8,
    blue: u8,
}

impl fmt::Display for RgbWithDisplay {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{},{},{}", self.red, self.green, self.blue)
    }
}

type RgbTriplet = (u8, u8, u8);

impl From<RgbTriplet> for RgbWithDisplay {
    fn from(rgb: RgbTriplet) -> Self {
        let (red, green, blue) = rgb;
        RgbWithDisplay { red, green, blue }
    }
}

#[derive(Debug, Deserialize, Serialize)]
#[serde(untagged)]
enum SmoothedMethylationPercentage {
    Valid(f64),
    Nan(String),
}

#[derive(Debug, Deserialize, Serialize)]
struct Row {
    contig: String,
    // Positions could be as large as 2.4e8 (chr1 length)
    start: u32,
    stop: u32,
    // For pipeline purposes, will just be "." (no name), but could change in future
    #[serde(skip_deserializing)]
    item_name: String,
    #[serde(skip_serializing)]
    methylation_percentage: f64,
    // We could probably get away with u8 here, but use u16 to be safe
    #[serde(skip_serializing)]
    non_converted_count: u16,
    #[serde(skip_serializing)]
    converted_count: u16,
    // Coverage capped at 1000: min(coverage, 1000)
    // See https://www.encodeproject.org/data-standards/wgbs/
    #[serde(skip_deserializing)]
    score: u16,
    #[serde(skip_deserializing)]
    strandedness: char,
    #[serde(skip_deserializing)]
    start_thick_display: u32,
    #[serde(skip_deserializing)]
    stop_thick_display: u32,
    #[serde(skip_deserializing)]
    color_value: String,
    #[serde(skip_deserializing)]
    coverage: u16,
    #[serde(skip_deserializing)]
    smoothed_methylation_percentage: f64,
}

#[derive(Debug, Deserialize, Serialize)]
struct SmoothedMethylationRow {
    smoothed_methylation_percentage: SmoothedMethylationPercentage,
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_get_args() {
        let args = Cli::from_iter(vec!["prog", "infile", "methvec", "outfile"]);
        assert_eq!(args.bismark_bed_file, std::path::PathBuf::from("infile"));
        assert_eq!(
            args.smoothed_methylation_tsv,
            std::path::PathBuf::from("methvec")
        );
        assert_eq!(args.encode_bed_outfile, std::path::PathBuf::from("outfile"));
    }

    #[test]
    fn test_rgb_from_methylation() {
        let RgbWithDisplay { red, green, blue } = rgb_from_methylation(50.);
        assert_eq!(red, 255);
        assert_eq!(green, 255);
        assert_eq!(blue, 0);
    }
    #[test]
    fn test_rgbwithdisplay_from_rgbtriplet() {
        let rgbtriplet: RgbTriplet = (100, 200, 3);
        let RgbWithDisplay { red, green, blue } = RgbWithDisplay::from(rgbtriplet);
        assert_eq!(red, 100);
        assert_eq!(green, 200);
        assert_eq!(blue, 3);
    }

    #[test]
    fn test_rgbwithdislay_to_string() {
        let x = RgbWithDisplay {
            red: 200,
            green: 200,
            blue: 200,
        };
        assert_eq!(x.to_string(), String::from("200,200,200"));
    }

    #[test]
    fn test_get_reader() -> Result<(), Box<dyn Error>> {
        let data = "\
                    chr1\t10649\t10650\t1.0\t4\t0\n\
                    ";
        let mut rdr = get_reader(data.as_bytes());
        let records = rdr
            .records()
            .collect::<Result<Vec<csv::StringRecord>, csv::Error>>()?;
        assert_eq!(records[0], vec!["chr1", "10649", "10650", "1.0", "4", "0",]);
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
        let bismark_data = "\
                            chr1\t10649\t10650\t1.0\t4\t0\n\
                            chr1\t10650\t10651\t1.0\t4\t0\n\
                            ";
        let bsmooth_data = "\
        0.356\n
        NA\n
        ";
        let bismark_rdr = get_reader(bismark_data.as_bytes());
        let bsmooth_rdr = get_reader(bsmooth_data.as_bytes());
        let mut wtr = get_writer(vec![]);
        process(bismark_rdr, bsmooth_rdr, &mut wtr)?;
        let result = String::from_utf8(wtr.into_inner()?)?;
        assert_eq!(
            result,
            "chr1\t10649\t10650\t.\t4\t.\t10649\t10650\t219,255,0\t4\t35.6\n"
        );
        Ok(())
    }
}
