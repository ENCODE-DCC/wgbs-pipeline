workflow wgbs {
	File reference
	File? indexed_reference
	File? indexed_contig_sizes
	File? extra_reference
	Array[Array[File]] fastqs
	Array[String] sample_names
	Array[String] sample_barcodes

	Array[File] fastqs_ = flatten(fastqs)

	call flatten_ { input:
		fastqs = fastqs_
	}

}

task flatten_ {
	Array[File] fastqs

	command { 
		mkdir mapping
		cat ${write_lines(fastqs)} | xargs -I % ln -s % mapping
		ls mapping
	}

	output {
		
	}

}
