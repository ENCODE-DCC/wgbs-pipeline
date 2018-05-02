# ENCODE WGBS Pipeline running https://github.com/heathsc/gemBS
# Maintainer: Ulugbek Baymuradov
# import 'index.wdl' as step1

workflow wgbs {

	File reference_fasta
	File metadata 
	Array[File] fastq_files
	
	call index { input:
		reference_fasta = reference_fasta
	}

	call prepare_config { input:
		metadata_csv = metadata
	}

	call generate_mapping_commands { input:
		reference_gem = index.reference_gem,
		metadata_json = prepare_config.metadata_json
	}

	call mapping { input:
		reference_gem = index.reference_gem,
		metadata_json = prepare_config.metadata_json,
		fastq_files = fastq_files,
		commands = generate_mapping_commands.mapping_commands
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

task generate_mapping_commands {

	File reference_gem
	File metadata_json

	command {
		mkdir fastqs
		mkdir -p data/mappings 
		mkdir tmp
		mkdir reference
		ln -s ${reference_gem} reference/
		ln -s ${metadata_json} .
		gemBS mapping-commands -I reference/$(basename ${reference_gem}) -j metadata.json -i ./fastqs/ -o ./data/mappings/ -d tmp/ -t 8 -p
	}

	output {
		Array[String] mapping_commands = read_lines(stdout())
	}
}

task mapping {
	File reference_gem
	File metadata_json
	Array[File] fastq_files
	Array[String] commands

	command {
		mkdir fastqs
		mkdir -p data/mappings 
		mkdir tmp
		mkdir reference
		ln -s ${reference_gem} reference/
		ln -s ${metadata_json} .
		cat ${write_lines(fastq_files)} | xargs -I % ln -s % fastqs
		${commands[0]}
	}

	output {

	}
}











