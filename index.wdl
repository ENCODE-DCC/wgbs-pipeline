workflow index {

	File reference_fasta

	call index_job { input:
		reference_fasta = reference_fasta
	}

		
	call pyglob as gem_file { input:
		pattern = "*.BS.gem",
		job_rc = index_job.return_code
	}

	call pyglob as info_file { input:
		pattern = "*.BS.info",
		job_rc = index_job.return_code
	}

	call pyglob as fasta_file { input:
		pattern = "*.fa",
		job_rc = index_job.return_code
	}

	output {
		#File reference_gem = if defined(gem_file.files) then gem_file.files[0] else ""
		#File reference_info = if defined(info_file.files) then info_file.files[0] else ""

		File reference_gem = gem_file.matches[0]
		File reference_info = info_file.matches[0]
		File original_fasta = fasta_file.matches[0]

	}
}

task index_job {

	File reference_fasta

	command {
		gemBS index -i ${reference_fasta}
		cd .. && cd .. && ls -lR
	}

	output {
		Int return_code = read_int("rc")
		Array[String] ls = read_lines(stdout())
	}
}

task pyglob {

	String pattern
	Int? nearness
	Int? nearness_value = if defined(nearness) then nearness else 5
	Int? job_rc

	command {
		python3 /software/helpers/glob_helper.py \
			${"--pattern '" + pattern + "'"} \
			${"--nearness " + nearness_value}
	}

	output {
		Array[File] matches = read_lines('matching_files.txt')
		Array[String] ls = read_lines(stdout())	
	}
}