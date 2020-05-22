import "../../../wgbs-pipeline.wdl" as wgbs

workflow test_bscaller {
    File reference
    File gemBS_json
    File bam
    File csi
    File index
    String sample_name
    String sample_barcode

    call wgbs.bscaller { input:
        reference = reference,
        gemBS_json = gemBS_json,
        bam = bam,
        csi = csi,
        index = index,
        sample_barcode = sample_barcode,
        sample_name = sample_name
    }
}
