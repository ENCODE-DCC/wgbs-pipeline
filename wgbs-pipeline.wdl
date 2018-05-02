# ENCODE WGBS Pipeline running https://github.com/heathsc/gemBS
# Maintainer: Ulugbek Baymuradov
# import 'index.wdl' as step1

workflow wgbs {

	File reference_fasta
	File metadata 
	File sequence
	
	call index { input:
		reference_fasta = reference_fasta,
	}

	call prepare_config { input:
		metadata_csv = metadata
	}

	output {
		File reference_info = index.reference_info
		File reference_gem = index.reference_gem
		File metadata_json = prepare_config.metadata_json
	}
}


task index {

	File reference_fasta

	command {
		gemBS index -i ${reference_fasta}
		pyglob.py \
			${"--pattern '*.BS.info'"} \
			${"--nearness " + 1} \
			${"--matched-files-name info"}
		pyglob.py \
			${"--pattern '*.BS.gem'"} \
			${"--nearness " + 1} \
			${"--matched-files-name gem"}
		mkdir index_out
		cat info.txt | xargs -I % ln -s % index_out
		cat gem.txt | xargs -I % ln -s % index_out

	}

	output {
		File reference_info = glob("index_out/*.BS.info")[0]
		File reference_gem = glob("index_out/*.BS.gem")[0]
	}
}

task prepare_config {
	File metadata_csv

	command {
		gemBS prepare-config -t ${metadata_csv} -j metadata.json
	}

	output {
		File metadata_json = glob("metadata.json")[0]
	}
}