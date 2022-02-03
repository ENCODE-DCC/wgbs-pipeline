import "../../../wgbs-pipeline.wdl" as wgbs

workflow test_extract {
    File reference
    File gemBS_json
    File contig_sizes
    File bcf
    File bcf_csi
    String sample_barcode
    String docker

    call wgbs.extract { input:
        reference = reference,
        gemBS_json = gemBS_json,
        contig_sizes = contig_sizes,
        bcf = bcf,
        bcf_csi = bcf_csi,
        sample_barcode = sample_barcode,
        docker = docker,
    }
}
