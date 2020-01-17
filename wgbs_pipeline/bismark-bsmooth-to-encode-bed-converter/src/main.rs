use csv::{ReaderBuilder, WriterBuilder};
use palette::rgb::Rgb;
use palette::Hsv;
use serde::{Deserialize, Serialize};
use std::error::Error;
use std::fs::File;
use std::{cmp, fmt};
use structopt::StructOpt;

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
    smoothed_methylation_percentage: f64,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::from_args();
    let bismark_bed_file = File::open(&args.bismark_bed_file)?;
    let smoothed_methylation = File::open(&args.smoothed_methylation_tsv)?;
    let encode_bed_file = File::create(&args.encode_bed_outfile)?;

    let mut bismark_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(&bismark_bed_file);

    let mut smoothed_methylation_reader = ReaderBuilder::new()
        .has_headers(false)
        .from_reader(&smoothed_methylation);

    let mut encode_writer = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(&encode_bed_file);

    let records = bismark_reader.deserialize::<Row>();
    let smoothed_methylation_records =
        smoothed_methylation_reader.deserialize::<SmoothedMethylationRow>();

    for (i, smoothed_methylation_record) in records.zip(smoothed_methylation_records) {
        let mut record = i?;
        record.item_name = ".".to_string();
        record.coverage = record.converted_count + record.non_converted_count;
        record.score = cmp::min(record.coverage, ENCODE_SCORE_CAP);
        record.strandedness = UNKNOWN_STRANDEDNESS;
        record.start_thick_display = record.start;
        record.stop_thick_display = record.stop;
        record.smoothed_methylation_percentage =
            smoothed_methylation_record?.smoothed_methylation_percentage;
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
    let rgb: Rgb = Hsv::new((1. - methylation as f32) * ENCODE_HSV_MAX_HUE, 1., 1.).into();
    RgbWithDisplay::from(rgb.into_format::<u8>().into_components())
}