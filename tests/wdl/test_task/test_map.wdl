import "../../../wgbs-pipeline.wdl" as wgbs

workflow test_map {
    File reference
    File index
    File gemBS_json
    Array[File] fastqs
    String sample_name
    String sample_barcode = "sample_sample5"

    call wgbs.map { input:
        reference = reference,
        index = index,
        gemBS_json = gemBS_json,
        fastqs = fastqs,
        sample_barcode = sample_barcode,
        sample_name = sample_name
    }
}
