import "../../../wgbs-pipeline.wdl" as wgbs

workflow test_make_conf {
    Int num_threads
    Int num_jobs
    String reference
    File extra_reference
    Boolean benchmark_mode = false
    String? include_file
    String? underconversion_sequence_name
    String docker

    call wgbs.make_conf { input:
        num_threads = num_threads,
        num_jobs = num_jobs,
        reference = reference,
        extra_reference = extra_reference,
        include_file = include_file,
        underconversion_sequence_name = underconversion_sequence_name,
        benchmark_mode = benchmark_mode,
        docker = docker,
    }
}
