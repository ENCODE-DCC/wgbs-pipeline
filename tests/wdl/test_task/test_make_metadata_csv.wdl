import "../../../wgbs-pipeline.wdl" as wgbs

workflow test_make_metadata_csv {
    Array[String] sample_names
    Array[Array[File]] fastqs

    call wgbs.make_metadata_csv { input:
        sample_names = sample_names,
        fastqs = fastqs
    }
}