workflow wgbs {
	File reference
	File? indexed_reference
	File? indexed_contig_sizes
	File? extra_reference
	Array[Array[File]] fastqs
	Array[String] sample_names
	Array[String] sample_barcodes

	call prepare { input:
		reference = reference,
		extra_reference = extra_reference
	}

	if (!defined(indexed_reference)) {
		call index as index_reference { input:
			reference = reference,
			extra_reference = extra_reference,
			gemBS_json = prepare.gemBS_json
		}
	}

	File index = select_first([indexed_reference, index_reference.BS_gem])
	File contig_sizes = select_first([indexed_contig_sizes, index_reference.contig_sizes])

	scatter(i in range(length(sample_names))) {
		call map { input:
			fastqs = fastqs[i],
			sample_barcode = sample_barcodes[i],
			sample_name = sample_names[i],
			index = index,
			reference = reference,
			gemBS_json = prepare.gemBS_json
		}

		call bscaller { input:
			reference = reference,
			gemBS_json = prepare.gemBS_json,
			sample_barcode = sample_barcodes[i],
			sample_name = sample_names[i],
			bam = map.bam,
			bai = map.bai,
			contig_sizes = contig_sizes
		}

		call extract { input:
			gemBS_json = prepare.gemBS_json,
			reference = reference,
			sample_barcode = sample_barcodes[i],
			sample_name = sample_names[i],
			bcf = bscaller.bcf,
			bcf_csi = bscaller.bcf_csi,
			contig_sizes = contig_sizes
		}
	}

	Array[File] map_qc_json_ = flatten(map.qc_json)
	Array[File] bscaller_qc_json_ = flatten(bscaller.qc_json)

	call qc_report { input:
		map_qc_json = map_qc_json_,
		bscaller_qc_json = bscaller_qc_json_,
		gemBS_json = prepare.gemBS_json,
		sample_names = sample_names,
		sample_barcodes = sample_barcodes,
		reference = reference
	}


	output {
		File index_used = index
		Array[File] bams = map.bam
		Array[File] bais = map.bai
		Array[File] bam_md5s = map.bam_md5
		Array[File] bscaller_bcf = bscaller.bcf
		Array[File] bscaller_bcf_csi = bscaller.bcf_csi
		Array[File] bscaller_bcf_md5 = bscaller.bcf_md5
		Array[File] extracted_bw = extract.bw
		Array[File] extracted_chg_bb = extract.chg_bb
		Array[File] extracted_chh_bb = extract.chh_bb
		Array[File] extracted_cpg_bb = extract.cpg_bb
		Array[File] extracted_chg_bed = extract.chg_bed
		Array[File] extracted_chh_bed = extract.chh_bed
		Array[File] extracted_cpg_bed = extract.cpg_bed
		Array[File] extracted_cpg_txt = extract.cpg_txt
		Array[File] extracted_cpg_txt_tbi = extract.cpg_txt_tbi
		Array[File] extracted_non_cpg_txt = extract.non_cpg_txt
		Array[File] extracted_non_cpg_txt_tbi = extract.non_cpg_txt_tbi
	}
}

task prepare {
	File configuration_file
	File metadata_file
	String reference
	String? extra_reference

	command {
		mkdir reference
		touch reference/$(basename ${reference})
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
		mkdir reference
		ln -s ${reference} reference && ln -s ${extra_reference} reference
		gemBS -j ${gemBS_json} index
	}

	output {
		File BS_gem = glob("indexes/*.BS.gem")[0]
		File BS_info = glob("indexes/*.BS.info")[0]
		File contig_sizes = glob("indexes/*.contig.sizes")[0]
	}
}

task map {
	File reference
	File index
	File gemBS_json
	Array[File] fastqs
	String sample_barcode
	String sample_name

	command {
		mkdir reference && ln ${reference} reference
		mkdir indexes && ln ${index} indexes
		mkdir -p fastq/${sample_name}
		cat ${write_lines(fastqs)} | xargs -I % ln % fastq/${sample_name}
		mkdir -p mapping/${sample_barcode}
		gemBS -j ${gemBS_json} map -b ${sample_barcode} --ignore-db
	}

	output {
		File bam = glob("mapping/**/*.bam")[0]
		File bai = glob("mapping/**/*.bai")[0]
		File bam_md5 = glob("mapping/**/*.bam.md5")[0]
		Array[File] qc_json = glob("mapping/**/*.json")
	}

}

task bscaller {
	File reference
	File gemBS_json
	File bam
	File bai
	File contig_sizes
	String sample_barcode 
	String sample_name

	command {
		mkdir reference && ln ${reference} reference
		mkdir indexes && ln ${contig_sizes} indexes
		mkdir -p mapping/${sample_barcode}
		ln -s ${bam} mapping/${sample_barcode}
		ln -s ${bai} mapping/${sample_barcode}
		gemBS -j ${gemBS_json} call --ignore-db --ignore-dep
	}

	output {
		File bcf = glob("calls/**/*.bcf")[0]
		File bcf_csi = glob("calls/**/*.bcf.csi")[0]
		File bcf_md5 = glob("calls/**/*.bcf.md5")[0]
		Array[File] qc_json = glob("calls/**/*.json")
	}
}

task extract {
	File gemBS_json
	File reference
	File contig_sizes
	File bcf
	File bcf_csi
	String sample_barcode 
	String sample_name

	command {
		mkdir reference && ln -s ${reference} reference
		mkdir indexes && ln -s ${contig_sizes} indexes
		mkdir -p calls/${sample_barcode}
		mkdir -p extract/${sample_barcode}
		ln -s ${bcf} calls/${sample_barcode} 
		ln -s ${bcf_csi} calls/${sample_barcode}
		gemBS -j ${gemBS_json} extract -w -B --ignore-db --ignore-dep
	}

	output {
		File bw = glob("extract/**/*.bw")[0]
		File chg_bb = glob("extract/**/*_chg.bb")[0]
		File chh_bb = glob("extract/**/*_chh.bb")[0]
		File cpg_bb = glob("extract/**/*_cpg.bb")[0]
		File chg_bed = glob("extract/**/*_chg.bed.gz")[0]
		File chh_bed = glob("extract/**/*_chh.bed.gz")[0]
		File cpg_bed = glob("extract/**/*_cpg.bed.gz")[0]
		File cpg_txt = glob("extract/**/*_cpg.txt.gz")[0]
		File cpg_txt_tbi = glob("extract/**/*_cpg.txt.gz.tbi")[0]
		File non_cpg_txt = glob("extract/**/*_non_cpg.txt.gz")[0]
		File non_cpg_txt_tbi = glob("extract/**/*_non_cpg.txt.gz.tbi")[0]
	}
}

task qc_report {
	Array[File] map_qc_json 
	Array[File] bscaller_qc_json
	Array[String] sample_names
	Array[String] sample_barcodes
	File reference
	File gemBS_json


	command {
		mkdir reference && mkdir mapping_reports && mkdir calls_reports
		ln -s ${reference} reference
		cat ${write_lines(sample_barcodes)} | xargs -I % mkdir -p mapping/%
		cat ${write_lines(sample_barcodes)} | xargs -I % mkdir -p calls/%
		cat ${write_lines(map_qc_json)} | xargs -I % ln -s % mapping/$(jq --raw-output '.ReadGroup | split("\t")[3] | split("BC:")[1]' %)
		cat ${write_lines(bscaller_qc_json)} | xargs -I % ln -s % calls/$(jq --raw-output '.ReadGroup | split("\t")[3] | split("BC:")[1]' %)
		gemBS -j ${gemBS_json} map-report -p ENCODE -o mapping_reports
		gemBS -j ${gemBS_json} call-report -p ENCODE -o calls_reports
	}

	output {
		File map_html = glob("mapping_reports/*.html")[0]
		File bscaller_html = glob("calls_reports/variant_calling/*.html")[0]
	}
}
