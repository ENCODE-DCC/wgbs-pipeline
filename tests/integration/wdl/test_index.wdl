import "../../../wgbs-pipeline.wdl" as wgbs

workflow test_index {
    File configuration_file
    File reference
    File extra_reference
    String docker

    call wgbs.index { input:
        configuration_file = configuration_file,
        reference = reference,
        extra_reference = extra_reference,
        docker = docker,
    }
}
