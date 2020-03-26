#CAPER docker quay.io/encode-dcc/wgbs-pipeline:0.1.0
#CAPER singularity docker://quay.io/encode-dcc/wgbs-pipeline:0.1.0
#CROO out_def https://raw.githubusercontent.com/ENCODE-DCC/wgbs-pipeline/dev/croo_out_def.json

workflow wgbs {
	File reference
	File? indexed_reference
	File? indexed_contig_sizes
	File extra_reference
	# biological replicate[technical replicate[one (single ended) or two paired fastqs]]]
	Array[Array[Array[File]]] fastqs
	Array[String] sample_names

	Int num_gembs_threads = 8
	Int num_gembs_jobs = 3
	String? underconversion_sequence_name
	String? include_conf_file

	String barcode_prefix = "sample_"
	Array[String] sample_barcodes = prefix(barcode_prefix, sample_names)

	Int bsmooth_num_workers = 8
	Int bsmooth_num_threads = 2
	Boolean run_bsmooth = true
	Boolean benchmark_mode = false

	# Optional params to tweak gemBS extract phred threshold and mininimum informative coverage for genotyped cytosines, otherwise will use gembs defaults
	Int? gembs_extract_phred_threshold
	Int? gembs_extract_min_inform

	call make_metadata_csv_and_conf { input:
		sample_names = sample_names,
		fastqs = write_json(fastqs),  # don't need the file contents, so avoid localizing
		barcode_prefix = barcode_prefix,
		num_threads = num_gembs_threads,
		num_jobs = num_gembs_jobs,
		reference = reference,
		extra_reference = extra_reference,
		include_file = include_conf_file,
		underconversion_sequence_name = underconversion_sequence_name,
		benchmark_mode = benchmark_mode
	}

	if (!defined(indexed_reference)) {
		call index as index_reference { input:
			configuration_file = make_metadata_csv_and_conf.gembs_conf,
			metadata_file = make_metadata_csv_and_conf.metadata_csv,
			reference = reference,
			extra_reference = extra_reference
		}
	}

	File index = select_first([indexed_reference, index_reference.gembs_indexes])
	File contig_sizes = select_first([indexed_contig_sizes, index_reference.contig_sizes])

	if (defined(indexed_reference) && defined(indexed_contig_sizes)) {
		call prepare { input:
			configuration_file = make_metadata_csv_and_conf.gembs_conf,
			metadata_file = make_metadata_csv_and_conf.metadata_csv,
			reference = reference,
			index = index,
			extra_reference = extra_reference
		}
	}

	File gemBS_json = select_first([prepare.gemBS_json, index_reference.gemBS_json])

	scatter(i in range(length(sample_names))) {
		call map { input:
			fastqs = flatten(fastqs[i]),
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
			csi = map.csi,
			index = index,
		}

		call extract { input:
			gemBS_json = gemBS_json,
			reference = reference,
			sample_barcode = sample_barcodes[i],
			sample_name = sample_names[i],
			bcf = bscaller.bcf,
			bcf_csi = bscaller.bcf_csi,
			contig_sizes = contig_sizes,
			phred_threshold = gembs_extract_phred_threshold,
			min_inform = gembs_extract_min_inform,
		}

		if (run_bsmooth) {
			call bsmooth { input:
				gembs_cpg_bed = extract.cpg_txt,
				chrom_sizes = contig_sizes,
				num_workers = bsmooth_num_workers,
				num_threads = bsmooth_num_threads
			}
		}
		call qc_report { input:
			map_qc_json = map.qc_json,
			gemBS_json = gemBS_json,
			reference = reference,
			contig_sizes = contig_sizes,
			sample_barcode = sample_barcodes[i]
		}
	}
}

task make_metadata_csv_and_conf {
	Array[String] sample_names
	File fastqs
	String barcode_prefix

	Boolean benchmark_mode
	Int num_threads
	Int num_jobs
	String reference
	File extra_reference
	String? include_file
	String? underconversion_sequence_name

	command {
		set -euo pipefail
		python3 $(which make_metadata_csv.py) \
			-n ${sep=' ' sample_names} \
			--files "${fastqs}" \
			-b "${barcode_prefix}" \
			-o "gembs_metadata.csv"

		python3 $(which make_conf.py) \
			-t "${num_threads}" \
			-j "${num_jobs}" \
			-r "${reference}" \
			-e "${extra_reference}" \
			${if defined(underconversion_sequence_name) then ("-u " + underconversion_sequence_name) else ""} \
			${if defined(include_file) then ("-i " + include_file) else ""} \
			${if benchmark_mode then ("--benchmark-mode") else ""} \
			-o "gembs.conf"
	}

	output {
		File metadata_csv = glob("gembs_metadata.csv")[0]
		File gembs_conf = glob("gembs.conf")[0]
	}

	runtime {
		cpu: 1
		memory: "2 GB"
		disks: "local-disk 10 SSD"
	}
}

task prepare {
	File configuration_file
	File metadata_file
	File index
	String reference
	String extra_reference

	command {
		set -euo pipefail
		mkdir reference && mkdir indexes
		touch reference/$(basename ${reference})
		touch reference/$(basename ${extra_reference})
		tar xf ${index} -C indexes --strip-components 1
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
	File extra_reference

	command {
		set -euo pipefail
		mkdir reference
		ln -s ${reference} reference && ln -s ${extra_reference} reference
		gemBS prepare -c ${configuration_file} \
					  -t ${metadata_file} \
					  -o gemBS.json \
					  --no-db
		gemBS -j gemBS.json index
		# See https://stackoverflow.com/a/54908072 . Want to make tar idempotent
		find indexes -type f -not -path '*.err' -not -path '*.info' -print0 | LC_ALL=C sort -z | tar --owner=0 --group=0 --numeric-owner --mtime='2019-01-01 00:00Z' --pax-option=exthdr.name=%d/PaxHeaders/%f,delete=atime,delete=ctime --no-recursion --null -T - -cf indexes.tar
		gzip -nc indexes.tar > indexes.tar.gz
	}

	output {
		File gembs_indexes = glob("indexes.tar.gz")[0]
		File contig_sizes = glob("indexes/*.contig.sizes")[0]
		File gemBS_json = glob("gemBS.json")[0]
	}

	runtime {
		memory: "64 GB"
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
		mkdir indexes && tar xf ${index} -C indexes --strip-components 1
		mkdir -p fastq/${sample_name}
		cat ${write_lines(fastqs)} | xargs -I % ln % fastq/${sample_name}
		mkdir -p mapping/${sample_barcode}
		gemBS -j ${gemBS_json} map -b ${sample_barcode} --ignore-db
	}

	output {
		File bam = glob("mapping/**/*.bam")[0]
		File csi = glob("mapping/**/*.csi")[0]
		File bam_md5 = glob("mapping/**/*.bam.md5")[0]
		Array[File] qc_json = glob("mapping/**/*.json")
	}

	runtime {
		memory: "64 GB"
	}
}

task bscaller {
	File reference
	File gemBS_json
	File bam
	File csi
	File index
	String sample_barcode
	String sample_name

	command {
		set -euo pipefail
		mkdir reference && ln ${reference} reference
		mkdir indexes && tar xf ${index} -C indexes --strip-components 1
		mkdir -p mapping/${sample_barcode}
		ln -s ${bam} mapping/${sample_barcode}
		ln -s ${csi} mapping/${sample_barcode}
		gemBS -j ${gemBS_json} call --ignore-db --ignore-dep
	}

	output {
		File bcf = glob("calls/**/*.bcf")[0]
		File bcf_csi = glob("calls/**/*.bcf.csi")[0]
		File bcf_md5 = glob("calls/**/*.bcf.md5")[0]
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
	Int? phred_threshold
	Int? min_inform

	command {
		set -euo pipefail
		mkdir reference && ln ${reference} reference
		mkdir indexes && ln ${contig_sizes} indexes
		mkdir -p calls/${sample_barcode}
		mkdir -p extract/${sample_barcode}
		ln ${bcf} calls/${sample_barcode}
		ln ${bcf_csi} calls/${sample_barcode}
		# htslib can complain the index is older than the bcf, so touch it to update timestamp
		touch ${bcf_csi}
		gemBS -j ${gemBS_json} extract \
			${if defined(phred_threshold) then ("-q " + phred_threshold) else ""} \
			${if defined(min_inform) then ("-l " + min_inform) else ""} \
			-B --ignore-db --ignore-dep
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

	runtime {
		disks: "local-disk 500 SSD"
		memory: "96 GB"
	}
}

task bsmooth {
	File gembs_cpg_bed
	File chrom_sizes
	Int num_workers
	Int num_threads

	command {
		set -euo pipefail
		gzip -cdf ${gembs_cpg_bed} > gembs_cpg.bed
		gembs-to-bismark-bed-converter gembs_cpg.bed converted_bismark.bed
		Rscript $(which bsmooth.R) -i converted_bismark.bed -o smoothed.tsv -w ${num_workers} -t ${num_threads}
		bismark-bsmooth-to-encode-bed-converter converted_bismark.bed smoothed.tsv smoothed_encode.bed
		bedToBigBed smoothed_encode.bed ${chrom_sizes} smoothed_encode.bb -type=bed9+2 -tab
		gzip -n smoothed_encode.bed
	}

	output {
		File smoothed_cpg_bed = glob("smoothed_encode.bed.gz")[0]
		File smoothed_cpg_bigbed = glob("smoothed_encode.bb")[0]
	}

	runtime {
		cpu: 16
		disks: "local-disk 500 SSD"
		memory: "128 GB"
	}
}

task qc_report {
	Array[File] map_qc_json
	File reference
	File gemBS_json
	File contig_sizes
	String sample_barcode

	command {
		set -euo pipefail
		mkdir reference && mkdir mapping_reports && mkdir calls_reports && mkdir indexes
		ln -s ${contig_sizes} indexes
		ln -s ${reference} reference
		mkdir -p "mapping/${sample_barcode}"
		cat ${write_lines(map_qc_json)} | xargs -I % ln % "mapping/${sample_barcode}"
		touch mapping/${sample_barcode}/"${sample_barcode}.bam"
		gemBS -j ${gemBS_json} map-report -p ENCODE -o mapping_reports
		python3 $(which parse_map_qc_html.py) -i mapping_reports/mapping/${sample_barcode}.html -o gembs_map_qc.json
	}

	output {
		Array[File] map_html_assets = glob("mapping_reports/mapping/*")
		File portal_map_qc_json = glob("gembs_map_qc.json")[0]
		File map_qc_insert_size_plot_png = glob("mapping_reports/mapping/${sample_barcode}.isize.png")[0]
		File map_qc_mapq_plot_png = glob("mapping_reports/mapping/${sample_barcode}.mapq.png")[0]
	}
}
