import "../../../wgbs-pipeline.wdl" as wgbs

workflow test_make_metadata_csv {
    Array[String] sample_names
    Array[Array[Array[File]]] fastqs
    String barcode_prefix

    call wgbs.make_metadata_csv { input:
        sample_names = sample_names,
        fastqs = write_json(fastqs),
        barcode_prefix = barcode_prefix,
    }
}
