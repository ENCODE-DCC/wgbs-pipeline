import "../../../wgbs-pipeline.wdl" as wgbs

workflow test_bsmooth {
    File gembs_cpg_bed
    File chrom_sizes
    Int num_threads
    Int num_workers

    call wgbs.bsmooth { input:
        gembs_cpg_bed = gembs_cpg_bed,
        chrom_sizes = chrom_sizes,
        num_threads = num_threads,
        num_workers = num_workers
    }
}
