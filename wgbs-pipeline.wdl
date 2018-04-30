# ENCODE WGBS Pipeline running https://github.com/heathsc/gemBS
# Maintainer: Ulugbek Baymuradov
# import 'index.wdl' as step1

workflow wgbs {

	File reference_fasta
	File metadata 
	File sequence
	
	call index_job { input:
		reference_fasta = reference_fasta,
	}

	output {
		File reference_info = index_job.reference_info
		#File reference_info = index.reference_info
		#File test_file = index.original_fasta
	}
}


task index_job {

	File reference_fasta

	command {
		gemBS index -i ${reference_fasta}
		python3 /software/helpers/glob_helper.py \
			${"--pattern '*.BS.info'"} \
			${"--nearness " + 1} \
			${"--matched-files-name info"}
	}

	output {
		File reference_info = glob(read_lines("info.txt")[0])[0]
	}
}