#CAPER docker quay.io/encode-dcc/wgbs-pipeline:0.1.0
#CAPER singularity docker://quay.io/encode-dcc/wgbs-pipeline:0.1.0

workflow wgbs {
	File configuration_file
	File reference
	File? indexed_reference
	File? indexed_contig_sizes
	File? extra_reference
	Array[Array[File]] fastqs
	Array[String] sample_names

	String barcode_prefix = "sample_"
	Array[String] sample_barcodes = prefix(barcode_prefix, sample_names)

	call make_metadata_csv { input:
		sample_names = sample_names,
		fastqs = write_tsv(fastqs),  # don't need the file contents, so avoid localizing
		barcode_prefix = barcode_prefix
	}

	if (!defined(indexed_reference)) {
		call index as index_reference { input:
			configuration_file = configuration_file,
			metadata_file = make_metadata_csv.metadata_csv,
			reference = reference,
			extra_reference = extra_reference
		}
	}

	File index = select_first([indexed_reference, index_reference.BS_gem])
	File contig_sizes = select_first([indexed_contig_sizes, index_reference.contig_sizes])

	if (defined(indexed_reference) && defined(indexed_contig_sizes)) {
		call prepare { input:
			configuration_file = configuration_file,
			metadata_file = make_metadata_csv.metadata_csv,
			contig_sizes = indexed_contig_sizes,
			reference = reference,
			index = index,
			extra_reference = extra_reference
		}
	}

	File gemBS_json = select_first([prepare.gemBS_json, index_reference.gemBS_json])

	scatter(i in range(length(sample_names))) {
		call map { input:
			fastqs = fastqs[i],
			sample_barcode = sample_barcodes[i],
			sample_name = sample_names[i],
			index = index,
			reference = reference,
			gemBS_json = gemBS_json
		}

		call bscaller { input:
			reference = reference,
			gemBS_json = gemBS_json,
			sample_barcode = sample_barcodes[i],
			sample_name = sample_names[i],
			bam = map.bam,
			bai = map.bai,
			contig_sizes = contig_sizes
		}

		call extract { input:
			gemBS_json = gemBS_json,
			reference = reference,
			sample_barcode = sample_barcodes[i],
			sample_name = sample_names[i],
			bcf = bscaller.bcf,
			bcf_csi = bscaller.bcf_csi,
			contig_sizes = contig_sizes
		}

		call bsmooth { input:
			cpg_gembs_bed = extract.cpg_txt,
			chrom_size = contig_sizes
		}
	}

	Array[File] map_qc_json_ = flatten(map.qc_json)
	Array[File] bscaller_qc_json_ = flatten(bscaller.qc_json)

	call qc_report { input:
		map_qc_json = map_qc_json_,
		bscaller_qc_json = bscaller_qc_json_,
		gemBS_json = gemBS_json,
		reference = reference,
		contig_sizes = contig_sizes
	}
}

task make_metadata_csv {
	Array[String] sample_names
	File fastqs
	String barcode_prefix

	command {
		set -euo pipefail
		python3 $(which make_metadata_csv.py) \
			-n "${sep=' ' sample_names}" \
			--files "${fastqs}" \
			-b "${barcode_prefix}"
	}

	output {
		File metadata_csv = glob("*_metadata.csv")[0]
	}
}

task prepare {
	File configuration_file
	File metadata_file
	File contig_sizes
	String reference
	String index
	String? extra_reference

	command {
		set -euo pipefail
		mkdir reference && mkdir indexes
		touch reference/$(basename ${reference})
		touch reference/$(basename ${extra_reference})
		touch indexes/$(basename ${index})
		ln ${contig_sizes} indexes/
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
	File configuration_file
	File metadata_file
	File reference
	File? extra_reference

	command {
		set -euo pipefail
		mkdir reference
		ln -s ${reference} reference && ln -s ${extra_reference} reference
		gemBS prepare -c ${configuration_file} \
					  -t ${metadata_file} \
					  -o gemBS.json \
					  --no-db
		gemBS -j gemBS.json index
	}

	output {
		File BS_gem = glob("indexes/*.BS.gem")[0]
		File BS_info = glob("indexes/*.BS.info")[0]
		File contig_sizes = glob("indexes/*.contig.sizes")[0]
		File gemBS_json = glob("gemBS.json")[0]
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
		set -euo pipefail
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
		set -euo pipefail
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
		set -euo pipefail
		mkdir reference && ln ${reference} reference
		mkdir indexes && ln ${contig_sizes} indexes
		mkdir -p calls/${sample_barcode}
		mkdir -p extract/${sample_barcode}
		ln ${bcf} calls/${sample_barcode}
		ln ${bcf_csi} calls/${sample_barcode}
		gemBS -j ${gemBS_json} extract -q 0 -l 0 -w -B --ignore-db --ignore-dep
	}

	output {
		File bw = glob("extract/**/*.bw")[0]
		File chg_bb = glob("extract/**/*_chg.bb")[0]
		File chh_bb = glob("extract/**/*_chh.bb")[0]
		File cpg_bb = glob("extract/**/*_cpg.bb")[0]
		File chg_bed = glob("extract/**/*_chg.bed.gz")[0]
		File chh_bed = glob("extract/**/*_chh.bed.gz")[0]
		File cpg_bed = glob("extract/**/*_cpg.bed.gz")[0]
		File cpg_txt = glob("extract/**/*_cpg.txt.gz")[0]  # gemBS-style output bed file
		File cpg_txt_tbi = glob("extract/**/*_cpg.txt.gz.tbi")[0]
		File non_cpg_txt = glob("extract/**/*_non_cpg.txt.gz")[0]
		File non_cpg_txt_tbi = glob("extract/**/*_non_cpg.txt.gz.tbi")[0]
	}
}

task bsmooth {
	File cpg_gembs_bed
	File chrom_sizes

	command {
		set -euo pipefail
		gembs-to-bismark-bed-converter ${cpg_gembs_bed} converted_bismark.bed
		Rscript bsmooth.R converted_bismark.bed smoothed.tsv
		bismark-bsmooth-to-encode-bed-converter converted_bismark.bed smoothed.tsv smoothed_encode.bed
		bedToBigBed smoothed_encode.bed ${chrom_sizes} smoothed_encode.bb -type=bed9+2 -tab
		gzip -n smoothed_encode.bed
	}

	output {
		File smoothed_cpg_bed = glob("smoothed_encode.bed.gz")[0]
		File smoothed_cpg_bigbed = glob("smoothed_encode.bb")[0]
	}
}

task qc_report {
	Array[File] map_qc_json
	Array[File] bscaller_qc_json
	File reference
	File gemBS_json
	File contig_sizes


	command {
		set -euo pipefail
		mkdir reference && mkdir mapping_reports && mkdir calls_reports && mkdir indexes
		ln -s ${contig_sizes} indexes
		ln -s ${reference} reference
		cat ${write_lines(map_qc_json)} | while read line
		do
			barcode=$(jq --raw-output '.ReadGroup | split("\t")[3] | split("BC:")[1]' $line)
			mkdir -p mapping/$barcode
			ln $line mapping/$barcode
			touch mapping/$barcode/"$barcode.bam"
		done
		cat ${write_lines(bscaller_qc_json)} | while read line
		do
			barcode=$(cut -d '_' -f1 <<< $(basename $line))
			mkdir -p calls/$barcode
			ln $line calls/$barcode
			touch calls/$barcode/"$barcode.bcf"
		done
		gemBS -j ${gemBS_json} map-report -p ENCODE -o mapping_reports
		gemBS -j ${gemBS_json} call-report -p ENCODE -o calls_reports
	}

	output {
		Array[File] map_html_assets = glob("mapping_reports/mapping/*")
		Array[File] bscaller_html_assets = glob("calls_reports/variant_calling/*")
	}
}
