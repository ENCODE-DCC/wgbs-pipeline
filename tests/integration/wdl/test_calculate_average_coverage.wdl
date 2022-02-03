import "../../../wgbs-pipeline.wdl" as wgbs

workflow test_calculate_average_coverage {
    File bam
    File chrom_sizes
    Int num_cpus
    String docker

    call wgbs.calculate_average_coverage { input:
        bam = bam,
        chrom_sizes = chrom_sizes,
        num_cpus = num_cpus,
        docker = docker,
    }
}
