import "../../../wgbs-pipeline.wdl" as wgbs

workflow test_make_coverage_bigwig {
    File encode_bed
    File chrom_sizes

    call wgbs.make_coverage_bigwig { input:
        encode_bed = encode_bed,
        chrom_sizes = chrom_sizes,
    }
}
