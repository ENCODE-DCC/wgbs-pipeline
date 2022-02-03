import "../../../wgbs-pipeline.wdl" as wgbs

workflow test_calculate_bed_pearson_correlation {
    File bed1
    File bed2
    String docker

    call wgbs.calculate_bed_pearson_correlation { input:
        bed1 = bed1,
        bed2 = bed2,
        docker = docker,
    }
}
