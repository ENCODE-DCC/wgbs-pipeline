workflow wgbs {
	File reference
	File? indexed_reference
	File? indexed_contig_sizes
	File? extra_reference
	Array[Array[File]] fastqs
	Array[String] sample_names
	Array[String] sample_barcodes

	call flatten { input:
		fastqs = fastqs
	}

}

task flatten {
	Array[Array[File]] fastqs

	command { 
		mkdir mapping
		flatten.py --tsv=${write_tsv(fastqs)} | xargs -I % ln -s % mapping
		ls mapping
	}

	output {
		
	}

}
