import "../../../wgbs-pipeline.wdl" as wgbs

workflow test_calculate_average_coverage {
    File bam
    File chrom_sizes
    Int ncpus

    call wgbs.calculate_average_coverage { input:
        bam = bam,
        chrom_sizes = chrom_sizes,
        ncpus = ncpus,
    }
}
