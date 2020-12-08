# ENCODE WGBS Pipeline Reference

## Table of Contents

* [Inputs](reference.md#inputs)
* * [Input Descriptions](reference.md#input-descriptions)
* [Outputs](reference.md#outputs)

## Inputs

A input JSON to run the full pipeline, including building `gemBS` indexes, will look like this. In this example we have two biological replicates, the first with two libraries and the second with one library, with each library sequenced paired-ended:

```json
{
  "wgbs.fastqs": [
    [
      [
        "rep1_lib1_R1.fastq.gz",
        "rep1_lib1_R2.fastq.gz"
      ],
      [
        "rep1_lib2_R1.fastq.gz",
        "rep1_lib2_R2.fastq.gz"
      ]
    ],
    [
      [
        "rep2_lib1_R1.fastq.gz",
        "rep2_lib1_R2.fastq.gz"
      ]
    ]
  ],
  "wgbs.reference": "hg38.fasta.gz",
  "wgbs.extra_reference": "lambda.fasta.gz",
  "wgbs.underconversion_sequence_name": "chrL"
}
```

If you plan on running the pipeline multiple times it is recommended to create the `gemBS` index file once, then use that as input to the runs instead of the raw sequence `fasta` files. To build the references, your input should look like the following:

```json
{
  "wgbs.extra_reference": "lambda.fasta.gz",
  "wgbs.index_only": true,
  "wgbs.reference": "hg38.fasta.gz",
  "wgbs.underconversion_sequence_name": "chrL"
}
```

This will generate a combined contig sizes for the reference and extra reference as well as the `gemBS` index file. Then, to use the built references, specfiy `indexed_contig_sizes` and `indexed_reference` in the input JSON, like so:

```json
{
  "wgbs.fastqs": [
    [
      [
        "rep1_lib1_R1.fastq.gz",
        "rep1_lib1_R2.fastq.gz"
      ]
    ]
  ],
  "wgbs.indexed_contig_sizes": "my.contig.sizes",
  "wgbs.indexed_reference": "indexes.tar.gz",
  "wgbs.reference": "hg38.fasta.gz",
  "wgbs.extra_reference": "lambda.fasta.gz",
  "wgbs.underconversion_sequence_name": "chrL"
}
```

### Input Descriptions

`wgbs.fastqs` is a 3-D array of `fastq` file URIs (e.g. local file paths, Google Cloud Storage URIs, HTTP links) indexed by `[replicate][library][read pair]`. Replicates are processed independently. Libraries are merged together within each replicate to produce a single bam file and set of downstream methylation files for each replicate. Note that mixed paired and single ended libraries are not allowed.

`wgbs.reference` is a `fasta` file corresponding to the genome of interest, e.g. https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/ for GRCh38 (human).

`wgbs.extra_reference` is a `fasta` file containing extra control sequences used for QC purposes. For example, lambda DNA, which is unmethylated by nature, is commonly added to WGBS libraries as a control to verify that unmethylated cytosines are converted to uracil at a satisfactory rate (> 99% per [ENCODE standards](https://www.encodeproject.org/data-standards/wgbs/)). As such, for a WGBS experiment using lambda control, a `fasta` containing the sequence of lambda genome would be used, e.g. https://www.encodeproject.org/files/lambda.fa/

`wgbs.underconversion_sequence_name` is a string corresponding to the name of the sequence in the `wgbs.extra_reference` `fasta` that contains the control sequence to be used for calculating the bisulfite conversion rate. For lambda, this would typically be `chrL`.

`wgbs.indexed_reference` is a `.tar.gz` archive containing the following files produced from `gemBS index` during this pipeline's index generation step. Here `${prefix}` is the basename of the fasta file used to produce it, without the extension:
* `${prefix}.gemBS.contig_md5`
* `${prefix}.gemBS.ref.gzi`
* `${prefix}.BS.gem`
* `${prefix}.gemBS.ref`
* `${prefix}.contig.sizes`
* `${prefix}.gemBS.ref.fai`

`wgbs.indexed_contig_sizes` is a chromosome (chrom) sizes-style file, which is simply a two-column `tsv` where the rows contain the chromosome name and size, in that order. This file is produced by this pipeline during index generation from the input reference `fasta` files.

Every task also has customizable resources (number of cpus, RAM, and disk size). There
are too many to list individually here, but they all take the same form in the input
JSON file: `wgbs.[TASK_NAME]_[RESOURCE] = [INTEGER]`, where `RESOURCE` is one of
`num_cpus`, `ram_gb`, and `disk_size_gb`. As an example, if you need to increase RAM to
80 GB for the `map` task then in the input JSON you would put the following line:

```json
"wgbs.map_ram_gb": 80
```

## Outputs

The pipeline produces more outputs than listed here, only the most useful ones are listed here.

#### Task `wgbs.index`

* `bam`: a `BAM` file containing mapping results for the replicate, with all libraries for that replicate merged together.
* `csi`: a `BAM` index in `csi` format for the above `BAM`, useful for visualization and subsetting, among other purposes. See [specifications](http://samtools.github.io/hts-specs/CSIv1.pdf) for details.

#### Task `wgbs.calculate_average_coverage`

* `average_coverage_qc`: a `JSON` file containing the average coverage of its input BAM file and `samtools stats` quality metrics. Average coverage is defined as `(aligned read counts * read length ) / total genome size`, where read length and aligned read counts are taken from the `samtools stats` output and `total genome size` is the sum of chromosome sizes in the `index_contig_sizes` file.

#### Task `wgbs.extract`

Unless otherwise specified, `bigBed` files are in ENCODE `bedMethyl` (`bed9+2`) format, see [ENCODE standards](https://www.encodeproject.org/data-standards/wgbs/) for details. One of each of these files is produced per replicate.

* `plus_strand_bw`: a `bigWig` file containing the methylation percentage signal on the plus strand
* `minus_strand_bw`: a `bigWig` file containing the methylation percentage signal on the minus strand
* `chg_bb`: a `bigBed` file containing the methylation at reference CHG sites
* `chh_bb`: a `bigBed` file containing the methylation at reference CHH sites
* `cpg_bb`: a `bigBed` file containing the methylation at reference CpG sites
* `chg_bed`: a `bed` file containing the methylation at reference CHG sites
* `chh_bed`: a `bed` file containing the methylation at reference CHH sites
* `cpg_bed`: a `bed` file containing the methylation at reference CpG sites
* `cpg_txt`: a `gemBS`-style `bed` file containing the methylation at genotyped CpG sites (i.e. sites that could confidently be called by `gemBS` as CpGs). For details on this format, see the [gemBS documentation](http://statgen.cnag.cat/GEMBS/UserGuide/_build/html/pipelineExtract.html#gembs-style-output-files)

#### Task `wgbs.make_coverage_bigwig`

* `coverage_bigwig`: a `bigWig` file containing the coverage genome-wide, one file per replicate.

#### Task `wgbs.calculate_bed_pearson_correlation`

* `bed_pearson_correlation_qc`: a `JSON` file containing the Pearson correlation between replicates of methylation percentage at reference CpG sites with 10 or more reads. Note that this file is only produced when there are exactly two replicates.

#### Task `wgbs.bsmooth`

* `smoothed_cpg_bed`: a `bed` file in ENCODE `bedMethyl` format containing the smoothed methylation percentage at **genotyped** CpG sites for the input replicate, as calculated using `Bsmooth` from the [bsseq R package](https://www.bioconductor.org/packages/release/bioc/html/bsseq.html)
* `smoothed_cpg_bigbed`: the above `bed` file converted to `bigBed` format

#### Task `wgbs.qc_report`

* `portal_map_qc_json`: a `JSON` file containing mapping quality metrics for the replicate. These values are scraped out of the `gemBS` HTML QC reports. See the [ENCODE schema](https://www.encodeproject.org/profiles/gembs_alignment_quality_metric) for details on which metrics are provided. Further descriptions of the metrics can be found in the [gemBS documentation](http://statgen.cnag.cat/GEMBS/UserGuide/_build/html/qualityControl.html)
* `map_qc_insert_size_plot_png`: a `PNG` file containing a plot of insert size distribution across reads
* `map_qc_mapq_plot_png`: a `PNG` file containing a plot of mapping quality (MAPQ) distribution across reads
