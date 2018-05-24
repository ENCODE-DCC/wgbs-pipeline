# ENCODE WGBS Pipeline running https://github.com/heathsc/gemBS
# Maintainer: Ulugbek Baymuradov
import 'bs-call.wdl' as bscaller

workflow wgbs {
	String organism
	File reference_fasta
	File? indexed_reference_gem
	File? indexed_reference_info
	File metadata 
	Int pyglob_nearness
	Array[File] fastq_files
	Array[String] chromosomes


	if (!defined(indexed_reference_gem)) {
		call index { input:
			reference_fasta = reference_fasta,
			nearness = pyglob_nearness
		}
	}

	File indexed_gem = select_first([indexed_reference_gem, index.reference_gem])
	File indexed_info = select_first([indexed_reference_info, index.reference_info])

	call prepare_config { input:
		metadata_csv = metadata
	}

	call get_sample_names { input:
		metadata_json = prepare_config.metadata_json
	}

	call generate_mapping_commands { input:
		reference_gem = indexed_gem,
		metadata_json = prepare_config.metadata_json
	}

	call mapping { input:
		reference_gem = indexed_gem,
		metadata_json = prepare_config.metadata_json,
		fastq_files = fastq_files,
		commands = generate_mapping_commands.mapping_commands
	}


	scatter(name in get_sample_names.names) {
		call merging_sample { input:
			mapping_outputs = mapping.all_outputs,
			metadata_json = prepare_config.metadata_json,
			sample = name
		}

		call bscaller.bscall { input:
			organism = organism,
			reference_fasta = reference_fasta,
			bam = merging_sample.bam,
			bai = merging_sample.bai,
			sample = name,
			chromosomes = chromosomes
		}

		call bscall_concatenate { input:
			sample = name,
			bcf_files = bscall.bcf_files
		}

		call methylation_filtering {input:
			merged_call_file = bscall_concatenate.merged_file
		}

	}


	output {
		File reference_info = indexed_info
		File reference_gem = indexed_gem
		File metadata_json = prepare_config.metadata_json
		Array[File] mapping_step_outputs = mapping.all_outputs
		Array[File] merged_bam = merging_sample.bam
		Array[File] merged_bai = merging_sample.bai
		Array[File] bscall_concatenated_file = bscall_concatenate.merged_file
		Array[File] methylation_filtered_file = methylation_filtering.filtered_meth_file
	}
}


task index {

	File reference_fasta
	Int nearness

	Int? memory_gb
	Int? cpu
	String? disks


	command {
		gemBS index -i ${reference_fasta}
		pyglob.py \
			${"--pattern '*.BS.info'"} \
			${"--nearness " + nearness} \
			${"--matched-files-name info"}
		pyglob.py \
			${"--pattern '*.BS.gem'"} \
			${"--nearness " + nearness} \
			${"--matched-files-name gem"}
		mkdir index_out
		cat info.txt | xargs -I % ln -s % index_out
		cat gem.txt | xargs -I % ln -s % index_out

	}

	output {
		File reference_info = glob("index_out/*.BS.info")[0]
		File reference_gem = glob("index_out/*.BS.gem")[0]
	}

	runtime {
		cpu: select_first([cpu,16])
		memory : "${select_first([memory_gb,'52'])} GB"
		disks : select_first([disks,"local-disk 100 HDD"])
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

task get_sample_names {
	File metadata_json

	command <<<
		python3 <<CODE
		import json
		with open('${metadata_json}') as file:
			metadata_json = json.load(file)
		names = []
		for key, value in metadata_json.items():
			names.append(value['sample_barcode'])
		for name in set(names):
			print(name)
		CODE
	>>>

	output {
		Array[String] names = read_lines(stdout())
	}
}

task generate_mapping_commands {

	String reference_gem
	File metadata_json

	command {
		mkdir fastqs
		mkdir -p data/mappings 
		mkdir tmp
		mkdir reference
		ln -s ${reference_gem} reference/
		ln -s ${metadata_json} .
		gemBS mapping-commands -I reference/$(basename ${reference_gem}) -j $(basename ${metadata_json}) -i fastqs/ -o data/mappings/ -d tmp/ -t 16 -p
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

	Int? memory_gb
	Int? cpu
	String? disks


	command {
		mkdir tmp
		mkdir fastqs
		mkdir -p data/mappings 
		mkdir reference
		ln -s ${reference_gem} reference/
		ln -s ${metadata_json} .
		cat ${write_lines(fastq_files)} | xargs -I % ln -s % fastqs
		${commands[0]}
	}

	output {
		# May need individual outputs for testing purposes
		# File mapping_bai_file = glob("data/mappings/*.bai")[0]
		# File mapping_bam_file = glob("data/mappings/*.bam")[0]
		# File mapping_json_file = glob("data/mappings/*.json")[0]
		Array[File] all_outputs = glob("data/mappings/*")
	}

	runtime {
		cpu: select_first([cpu,16])
		memory : "${select_first([memory_gb,'104'])} GB"
		disks : select_first([disks,"local-disk 100 HDD"])
	}
}

task merging_sample {
	Array[File] mapping_outputs
	File metadata_json
	String sample 

	command {
		mkdir temp
		mkdir -p data/mappings
		cat ${write_lines(mapping_outputs)} | xargs -I % ln % data/mappings
		gemBS merging-sample -i data/mappings -j ${metadata_json} -s ${sample} -t 8 -o data/sample_mappings -d tmp/
	}

	output {
		File bam = glob("data/sample_mappings/*.bam")[0]
		File bai = glob("data/sample_mappings/*.bai")[0]
	}
}

task bscall_concatenate {
	String sample
	Array[File] bcf_files

	command {
		mkdir -p data/merged_calls
		gemBS bscall-concatenate -s ${sample} -l ${sep=" " bcf_files} -o data/merged_calls
	}

	output {
		File merged_file = glob("data/merged_calls/*.raw.bcf")[0]
	}
}

task methylation_filtering {
	File merged_call_file
	command {
		mkdir -p data/filtered_meth_calls
		gemBS methylation-filtering -b ${merged_call_file} -o data/filtered_meth_calls
	}

	output {
		File filtered_meth_file = glob("data/filtered_meth_calls/*.txt.gz")[0]
	}
}





































