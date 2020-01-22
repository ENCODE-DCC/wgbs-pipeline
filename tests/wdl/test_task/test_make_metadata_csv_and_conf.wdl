import "../../../wgbs-pipeline.wdl" as wgbs

workflow test_make_metadata_csv_and_conf {
    Array[String] sample_names
    Array[Array[File]] fastqs
    String barcode_prefix

	Int num_threads
	Int num_jobs
	String reference
	String extra_reference
	String? include_file
	String? underconversion_sequence_name

    call wgbs.make_metadata_csv_and_conf { input:
        sample_names = sample_names,
        fastqs = write_tsv(fastqs),
        barcode_prefix = barcode_prefix,
		num_threads = num_threads,
		num_jobs = num_jobs,
		reference = reference,
		extra_reference = extra_reference,
		include_file = include_file,
		underconversion_sequence_name = underconversion_sequence_name
    }
}
