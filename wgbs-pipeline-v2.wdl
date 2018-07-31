workflow wgbs {
	File reference
	File? extra_reference
	Array[Array[File]] fastqs
	Array[String] samples

	call prepare {input:
		reference = reference,
		extra_reference = extra_reference
	}

	call index {input:
		reference = reference,
		extra_reference = extra_reference,
		gemBS_json = prepare.gemBS_json
		}
	
}

task prepare {
	File configuration_file
	File metadata_file
	String reference
	String? extra_reference


	command {
		mkdir reference && touch reference/$(basename ${reference})
		touch reference/$(basename ${extra_reference})
		gemBS prepare -c ${configuration_file} \
					  -t ${metadata_file} \
					  -o gemBS.json \
					  --no-db
	}

	output {
		File gemBS_json = glob("gemBS.json")[0]
	}
}

task index {
	File reference
	File? extra_reference
	File gemBS_json

	command {
		mkdir reference && ln -s ${reference} reference && ln ${extra_reference} reference
		gemBS -j ${gemBS_json} index
	}

	output {
		File BS_gem = glob("indexes/*.BS.gem")[0]
		File BS_info = glob("indexes/*.BS.info")[0]
		File contig_sizes = glob("indexes/*.contig.sizes")[0]
	}
}