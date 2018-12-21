workflow mapping {
	File reference_gem
	File metadata_json
	Map[String, String] commands

	Pair[String, Array[Pair[String,Array[File]]]] sample_files

	scatter(lane_index in sample_files.right) {
		call mapping_job { input:
			fastq_files = lane_index.right,
			reference_gem = reference_gem,
			metadata_json = metadata_json,
			command_string = commands[lane_index.left]
		}
	}

	output {
		Array[Array[File]] mapping_outputs = mapping_job.all_outputs
	}
}

task mapping_job {
	File reference_gem
	File metadata_json
	String command_string
	Array[File] fastq_files

	Int? memory_gb
	Int? cpu
	String? disks


	command {
		mkdir tmp
		mkdir fastqs
		mkdir -p data/mappings 
		mkdir reference
		ln ${reference_gem} reference/
		ln ${metadata_json} .
		cat ${write_lines(fastq_files)} | xargs -I % ln % fastqs
		${command_string}
	}

	runtime {
		cpu: select_first([cpu,16])
		memory : "${select_first([memory_gb,'60'])} GB"
		disks : select_first([disks,"local-disk 700 HDD"])
		preemptible: 0
	}

	output {
		# May need individual outputs for testing purposes
		# File mapping_bai_file = glob("data/mappings/*.bai")[0]
		# File mapping_bam_file = glob("data/mappings/*.bam")[0]
		# File mapping_json_file = glob("data/mappings/*.json")[0]
		Array[File] all_outputs = glob("data/mappings/*")
	}
}