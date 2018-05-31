# ENCODE WGBS Pipeline running https://github.com/heathsc/gemBS
# Maintainer: Ulugbek Baymuradov
import 'bs-call.wdl' as bscaller
import 'mapping.wdl' as mapper

workflow wgbs {
	String organism
	File reference_fasta
	File? indexed_reference_gem
	File? indexed_reference_info
	File metadata 
	Int pyglob_nearness
	Array[Pair[String, Array[Pair[String,Array[File]]]]] fastq_files
	Array[File]? mapping_outputs
	Array[String] chromosomes
	Array[String] samples


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


	call generate_mapping_commands { input:
		reference_gem = indexed_gem,
		metadata_json = prepare_config.metadata_json
	}

	call flowcell_to_commands { input:
		metadata_json = prepare_config.metadata_json,
		commands = generate_mapping_commands.mapping_commands
	}

	scatter(sample_files in fastq_files) {
		
		call mapper.mapping as map { input:
			reference_gem = indexed_gem,
			metadata_json = prepare_config.metadata_json,
			sample_files = sample_files,
			commands = flowcell_to_commands.paired_commands
		}

		call merging_sample { input:
			mapping_outputs = map.mapping_outputs,
			metadata_json = prepare_config.metadata_json,
			sample = sample_files.left
		}

		call bscaller.bscall { input:
			organism = organism,
			reference_fasta = reference_fasta,
			bam = merging_sample.bam,
			bai = merging_sample.bai,
			sample = sample_files.left,
			chromosomes = chromosomes
		}

		call bscall_concatenate { input:
			sample = sample_files.left,
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
		Array[Array[Array[File]]] mapping_step_outputs = map.mapping_outputs
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
		memory : "${select_first([memory_gb,'60'])} GB"
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

task flowcell_to_commands {
	File metadata_json
	Array[String] commands

	command <<<
		python3 <<CODE
		import json
		with open('${metadata_json}') as file:
			metadata_json = json.load(file)
		lane_names = []
		for key, value in metadata_json.items():
			lane_names.append("{}_{}_{}".format(value['flowcell_name'], 
											  value['index_name'], 
											  value['lane_number']))
		for name in set(lane_names):
			for command_string in '${sep="<>" commands}'.split('<>'):
				if name in command_string:
					print("{}\t{}".format(name, command_string))
		CODE
	>>>

	output {
		Map[String,String] paired_commands = read_map(stdout())
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
		gemBS mapping-commands -I reference/$(basename ${reference_gem}) -j $(basename ${metadata_json}) -i fastqs/ -o data/mappings/ -d tmp/ -t 14 -p
	}

	output {
		Array[String] mapping_commands = read_lines(stdout())
	}
}


task merging_sample {
	Array[Array[File]] mapping_outputs
	File metadata_json
	String sample 

	Int? memory_gb
	Int? cpu
	String? disks


	command {
		mkdir temp
		mkdir -p data/mappings
		cat ${write_lines(mapping_outputs[0])} | xargs -I % ln % data/mappings
		gemBS merging-sample -i data/mappings -j ${metadata_json} -s ${sample} -t 14 -o data/sample_mappings -d tmp/
	}

	output {
		File bam = glob("data/sample_mappings/*.bam")[0]
		File bai = glob("data/sample_mappings/*.bai")[0]
	}

	runtime {
		cpu: select_first([cpu,16])
		memory : "${select_first([memory_gb,'60'])} GB"
		disks : select_first([disks,"local-disk 500 HDD"])
	}
}

task bscall_concatenate {
	String sample
	Array[File] bcf_files

	Int? memory_gb
	Int? cpu
	String? disks


	command {
		mkdir -p data/merged_calls
		gemBS bscall-concatenate -s ${sample} -l ${sep=" " bcf_files} -o data/merged_calls
	}

	output {
		File merged_file = glob("data/merged_calls/*.raw.bcf")[0]
	}

	runtime {
		cpu: select_first([cpu,2])
		memory : "${select_first([memory_gb,'7'])} GB"
		disks : select_first([disks,"local-disk 100 HDD"])
	}
}

task methylation_filtering {
	File merged_call_file

	Int? memory_gb
	Int? cpu
	String? disks

	command {
		mkdir -p data/filtered_meth_calls
		gemBS methylation-filtering -b ${merged_call_file} -o data/filtered_meth_calls
	}

	output {
		File filtered_meth_file = glob("data/filtered_meth_calls/*.txt.gz")[0]
	}

	runtime {
		cpu: select_first([cpu,16])
		memory : "${select_first([memory_gb,'60'])} GB"
		disks : select_first([disks,"local-disk 500 HDD"])
	}
}