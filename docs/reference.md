## Inputs

A typical input JSON for the pipeline will look like this:

```json
{
  "wgbs.extra_reference": "tests/data/conversion_control.fa.gz",
  "wgbs.fastqs": [
    [
      [
        "tests/data/sample5_data_1_200000.fastq.gz",
        "tests/data/sample5_data_2_200000.fastq.gz"
      ]
    ]
  ],
  "wgbs.include_conf_file": "/software/conf/IHEC_standard.conf",
  "wgbs.reference": "tests/data/sacCer3.fa.gz",
  "wgbs.run_bsmooth": false,
  "wgbs.underconversion_sequence_name": "NC_001416.1"
}
```

### gemBS index

The gemBS index is a `.tar.gz` archive containing the following files produced by running gemBS index, where `${prefix}` is the name of the fasta file used to produce it (without file extension):
* `${prefix}.gemBS.contig_md5`
* `${prefix}.gemBS.ref.gzi`
* `${prefix}.BS.gem`
* `${prefix}.gemBS.ref`
* `${prefix}.contig.sizes`
* `${prefix}.gemBS.ref.fai`


## Outputs
