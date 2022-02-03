import "../../../wgbs-pipeline.wdl" as wgbs

workflow test_prepare {
    File configuration_file
    File metadata_file
    File index
    String reference
    String extra_reference
    String docker

    call wgbs.prepare { input:
        configuration_file = configuration_file,
        metadata_file = metadata_file,
        index = index,
        reference = reference,
        extra_reference = extra_reference,
        docker = docker,
    }
}
