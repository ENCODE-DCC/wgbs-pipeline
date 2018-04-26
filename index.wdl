workflow index {
	File reference_fasta

	call index_job { input:
		reference_fasta = reference_fasta
	}

	call pyglob as fetch_gem_file { input:
		pattern = "./**/*.BS.gem",
	}

	call pyglob as fetch_info_file { input:
		pattern = "./**/*.BS.info",
	}

	output {
		File reference_gem = fetch_gem_file.files[0]
		File reference_info = fetch_info_file.files[0]
	}
}

task index_job {

	File reference_fasta

	command {
		gemBS index -i ${reference_fasta}
	}

	output {
		String placeholder = read_string(stdout())
	}
}

task pyglob {

	String pattern
	String? nearest_node

	command {
		python3 /software/helpers/glob_helper.py \
			${"--pattern '" + pattern + "'"} 
	}

	output {
		Array[File] files = read_lines(stdout())
	}
}