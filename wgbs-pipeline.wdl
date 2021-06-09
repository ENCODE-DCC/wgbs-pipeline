workflow wgbs {
    meta {
        version: "1.1.5"
        caper_docker: "encodedcc/wgbs-pipeline:1.1.5"
        caper_singularity: "docker://encodedcc/wgbs-pipeline:1.1.5"
        croo_out_def: "https://raw.githubusercontent.com/ENCODE-DCC/wgbs-pipeline/dev/croo_out_def.json"
    }

    File reference
    File? indexed_reference
    File? indexed_contig_sizes
    File? extra_reference
    # biological replicate[technical replicate[one (single ended) or two paired fastqs]]]
    Array[Array[Array[File]]]? fastqs
    Array[String]? sample_names

    Int num_gembs_threads = 8
    Int num_gembs_jobs = 3
    String? underconversion_sequence_name
    String pipeline_type = "wgbs"
    String include_conf_file = if pipeline_type=="wgbs" then "/software/conf/IHEC_standard.conf" else "/software/conf/IHEC_RRBS.conf"

    String barcode_prefix = "sample_"

    # Flag to only generate index tarball
    Boolean index_only = false
    Boolean benchmark_mode = false

    # Optional params to tweak gemBS extract phred threshold and mininimum informative coverage for genotyped cytosines, otherwise will use gembs defaults
    Int? gembs_extract_phred_threshold
    Int? gembs_extract_min_inform

    # Per task resource parameters
    Int? bscaller_disk_size_gb
    Int? bscaller_num_cpus
    Int? bscaller_ram_gb
    Int? calculate_average_coverage_disk_size_gb
    Int? calculate_average_coverage_num_cpus
    Int? calculate_average_coverage_ram_gb
    Int? calculate_bed_pearson_correlation_disk_size_gb
    Int? calculate_bed_pearson_correlation_num_cpus
    Int? calculate_bed_pearson_correlation_ram_gb
    Int? extract_disk_size_gb
    Int? extract_num_cpus
    Int? extract_ram_gb
    Int? index_disk_size_gb
    Int? index_num_cpus
    Int? index_ram_gb
    Int? make_conf_disk_size_gb
    Int? make_conf_num_cpus
    Int? make_conf_ram_gb
    Int? make_coverage_bigwig_disk_size_gb
    Int? make_coverage_bigwig_num_cpus
    Int? make_coverage_bigwig_ram_gb
    Int? make_metadata_csv_disk_size_gb
    Int? make_metadata_csv_num_cpus
    Int? make_metadata_csv_ram_gb
    Int? map_disk_size_gb
    Int? map_num_cpus
    Int? map_ram_gb
    String? map_sort_memory
    Int? map_sort_threads
    Int? prepare_disk_size_gb
    Int? prepare_num_cpus
    Int? prepare_ram_gb
    Int? qc_report_disk_size_gb
    Int? qc_report_num_cpus
    Int? qc_report_ram_gb

    # Don't need metadata csv to create indexes
    if (!index_only) {
        call make_metadata_csv { input:
            sample_names = select_first([sample_names, range(length(select_first([fastqs])))]),
            fastqs = write_json(fastqs),  # don't need the file contents, so avoid localizing
            barcode_prefix = barcode_prefix,
            num_cpus = make_metadata_csv_num_cpus,
            ram_gb = make_metadata_csv_ram_gb,
            disk_size_gb = make_metadata_csv_disk_size_gb,
        }
    }

    call make_conf { input:
        num_threads = num_gembs_threads,
        num_jobs = num_gembs_jobs,
        reference = reference,
        extra_reference = extra_reference,
        include_file = include_conf_file,
        underconversion_sequence_name = underconversion_sequence_name,
        benchmark_mode = benchmark_mode,
        num_cpus = make_conf_num_cpus,
        ram_gb = make_conf_ram_gb,
        disk_size_gb = make_conf_disk_size_gb,
    }

    if (!defined(indexed_reference)) {
        call index as index_reference { input:
            configuration_file = make_conf.gembs_conf,
            reference = reference,
            extra_reference = extra_reference,
            num_cpus = index_num_cpus,
            ram_gb = index_ram_gb,
            disk_size_gb = index_disk_size_gb,
        }
    }

    File index = select_first([indexed_reference, index_reference.gembs_indexes])
    File contig_sizes = select_first([indexed_contig_sizes, index_reference.contig_sizes])

    if (!index_only) {
        call prepare { input:
            configuration_file = make_conf.gembs_conf,
            metadata_file = select_first([make_metadata_csv.metadata_csv]),
            reference = reference,
            index = index,
            extra_reference = select_first([extra_reference, reference]),
            num_cpus = prepare_num_cpus,
            ram_gb = prepare_ram_gb,
            disk_size_gb = prepare_disk_size_gb,
        }
    }

    if (!index_only && defined(fastqs)) {
        Array[Array[Array[File]]] fastqs_ = select_first([fastqs])
        File gemBS_json = select_first([prepare.gemBS_json])
        Array[String] sample_names_ = select_first([sample_names, range(length(fastqs_))])
        Array[String] sample_barcodes = prefix(barcode_prefix, sample_names_)

        scatter(i in range(length(sample_names_))) {
            call map { input:
                fastqs = flatten(fastqs_[i]),
                sample_barcode = sample_barcodes[i],
                sample_name = sample_names_[i],
                index = index,
                reference = reference,
                gemBS_json = gemBS_json,
                num_cpus = map_num_cpus,
                ram_gb = map_ram_gb,
                disk_size_gb = map_disk_size_gb,
                sort_threads = map_sort_threads,
                sort_memory = map_sort_memory,
            }

            call bscaller { input:
                reference = reference,
                gemBS_json = gemBS_json,
                sample_barcode = sample_barcodes[i],
                sample_name = sample_names_[i],
                bam = map.bam,
                csi = map.csi,
                index = index,
                num_cpus = bscaller_num_cpus,
                ram_gb = bscaller_ram_gb,
                disk_size_gb = bscaller_disk_size_gb,
            }

            call calculate_average_coverage { input:
                bam = map.bam,
                chrom_sizes = contig_sizes,
                num_cpus = calculate_average_coverage_num_cpus,
                ram_gb = calculate_average_coverage_ram_gb,
                disk_size_gb = calculate_average_coverage_disk_size_gb,
            }

            call extract { input:
                gemBS_json = gemBS_json,
                reference = reference,
                sample_barcode = sample_barcodes[i],
                bcf = bscaller.bcf,
                bcf_csi = bscaller.bcf_csi,
                contig_sizes = contig_sizes,
                phred_threshold = gembs_extract_phred_threshold,
                min_inform = gembs_extract_min_inform,
                num_cpus = extract_num_cpus,
                ram_gb = extract_ram_gb,
                disk_size_gb = extract_disk_size_gb,
            }

            call make_coverage_bigwig { input:
                encode_bed = extract.cpg_bed,
                chrom_sizes = contig_sizes,
                num_cpus = make_coverage_bigwig_num_cpus,
                ram_gb = make_coverage_bigwig_ram_gb,
                disk_size_gb = make_coverage_bigwig_disk_size_gb,
            }

            call qc_report { input:
                map_qc_json = map.qc_json,
                gemBS_json = gemBS_json,
                reference = reference,
                contig_sizes = contig_sizes,
                sample_barcode = sample_barcodes[i],
                num_cpus = qc_report_num_cpus,
                ram_gb = qc_report_ram_gb,
                disk_size_gb = qc_report_disk_size_gb,
            }
        }

        Array[File] bedmethyls = extract.cpg_bed
        if (length(bedmethyls) == 2) {
            call calculate_bed_pearson_correlation { input:
                bed1 = bedmethyls[0],
                bed2 = bedmethyls[1],
                num_cpus = calculate_bed_pearson_correlation_num_cpus,
                ram_gb = calculate_bed_pearson_correlation_ram_gb,
                disk_size_gb = calculate_bed_pearson_correlation_disk_size_gb,
            }
        }
    }
}

task make_metadata_csv {
    Array[String] sample_names
    File fastqs
    String barcode_prefix
    Int? num_cpus = 1
    Int? ram_gb = 2
    Int? disk_size_gb = 10

    command {
        set -euo pipefail
        python3 "$(which make_metadata_csv.py)" \
            -n ${sep=' ' sample_names} \
            --files "${fastqs}" \
            -b "${barcode_prefix}" \
            -o "gembs_metadata.csv"
    }

    output {
        File metadata_csv = "gembs_metadata.csv"
    }

    runtime {
        cpu: "${num_cpus}"
        memory: "${ram_gb} GB"
        disks: "local-disk ${disk_size_gb} SSD"
    }
}

task make_conf {
    Boolean benchmark_mode
    Int num_threads
    Int num_jobs
    String reference
    File? extra_reference
    String? include_file
    String? underconversion_sequence_name
    Int? num_cpus = 1
    Int? ram_gb = 2
    Int? disk_size_gb = 10

    command {
        set -euo pipefail
        python3 "$(which make_conf.py)" \
            -t "${num_threads}" \
            -j "${num_jobs}" \
            -r "${reference}" \
            ${if defined(extra_reference) then ("-e " + extra_reference) else ""} \
            ${if defined(underconversion_sequence_name) then ("-u " + underconversion_sequence_name) else ""} \
            ${if defined(include_file) then ("-i " + include_file) else ""} \
            ${if benchmark_mode then ("--benchmark-mode") else ""} \
            -o "gembs.conf"
    }

    output {
        File gembs_conf = "gembs.conf"
    }

    runtime {
        cpu: "${num_cpus}"
        memory: "${ram_gb} GB"
        disks: "local-disk ${disk_size_gb} SSD"
    }
}

task prepare {
    File configuration_file
    File metadata_file
    File index
    String reference
    String extra_reference
    Int? num_cpus = 8
    Int? ram_gb = 32
    Int? disk_size_gb = 500

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
        File gemBS_json = "gemBS.json"
    }

    runtime {
        cpu: "${num_cpus}"
        memory: "${ram_gb} GB"
        disks: "local-disk ${disk_size_gb} HDD"
    }
}

task index {
    File configuration_file
    File reference
    File? extra_reference
    Int? num_cpus = 8
    Int? ram_gb = 64
    Int? disk_size_gb = 500

    command {
        set -euo pipefail
        mkdir reference
        ln -s ${reference} reference
        ${if defined(extra_reference) then ("ln -s " + extra_reference + " reference") else ""}
        # We need a gemBS JSON to create the index, and we need a metadata csv to create
        # the gemBS JSON. Create a dummy one, we don't use the config later so it
        # doesn't matter.
        cat << EOF > metadata.csv
        Barcode,Name,Dataset,File1,File2
        1,2,3,1.fastq.gz,2.fastq.gz
        EOF
        gemBS prepare \
            -c ${configuration_file} \
            -t metadata.csv \
            -o gemBS.json \
            --no-db
        gemBS -j gemBS.json index
        # See https://stackoverflow.com/a/54908072 . Want to make tar idempotent
        find indexes -type f -not -path '*.err' -not -path '*.info' -print0 | LC_ALL=C sort -z | tar --owner=0 --group=0 --numeric-owner --mtime='2019-01-01 00:00Z' --pax-option=exthdr.name=%d/PaxHeaders/%f,delete=atime,delete=ctime --no-recursion --null -T - -cf indexes.tar
        gzip -nc indexes.tar > indexes.tar.gz
    }

    output {
        File gembs_indexes = "indexes.tar.gz"
        File contig_sizes = glob("indexes/*.contig.sizes")[0]
    }

    runtime {
        cpu: "${num_cpus}"
        memory: "${ram_gb} GB"
        disks: "local-disk ${disk_size_gb} HDD"
    }
}

task map {
    File reference
    File index
    File gemBS_json
    Array[File] fastqs
    String sample_barcode
    String sample_name
    Int? num_cpus = 8
    Int? ram_gb = 64
    Int? disk_size_gb = 500
    Int? sort_threads
    String? sort_memory

    command {
        set -euo pipefail
        mkdir reference && ln ${reference} reference
        mkdir indexes && tar xf ${index} -C indexes --strip-components 1
        mkdir -p fastq/${sample_name}
        cat ${write_lines(fastqs)} | xargs -I % ln % fastq/${sample_name}
        mkdir -p mapping/${sample_barcode}
        gemBS \
            -j ${gemBS_json} \
            map \
            -b ${sample_barcode} \
            ${if defined(sort_threads) then ("--sort-threads " + sort_threads) else ""} \
            ${if defined(sort_memory) then ("--sort-memory " + sort_memory) else ""} \
            --ignore-db
    }

    output {
        File bam = glob("mapping/**/*.bam")[0]
        File csi = glob("mapping/**/*.csi")[0]
        File bam_md5 = glob("mapping/**/*.bam.md5")[0]
        Array[File] qc_json = glob("mapping/**/*.json")
    }

    runtime {
        cpu: "${num_cpus}"
        memory: "${ram_gb} GB"
        disks: "local-disk ${disk_size_gb} HDD"
    }
}

task calculate_average_coverage {
    File bam
    File chrom_sizes
    Int? num_cpus = 1
    Int? ram_gb = 8
    Int? disk_size_gb = 200

    command {
        python3 "$(which calculate_average_coverage.py)" \
            --bam ${bam} \
            --chrom-sizes ${chrom_sizes} \
            --threads ${num_cpus} \
            --outfile average_coverage_qc.json
    }

    output {
        File average_coverage_qc = "average_coverage_qc.json"
    }

    runtime {
        cpu: "${num_cpus}"
        memory: "${ram_gb} GB"
        disks: "local-disk ${disk_size_gb} SSD"
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
    Int? num_cpus = 8
    Int? ram_gb = 32
    Int? disk_size_gb = 500

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

    runtime {
        cpu: "${num_cpus}"
        memory: "${ram_gb} GB"
        disks: "local-disk ${disk_size_gb} HDD"
    }
}

task extract {
    File gemBS_json
    File reference
    File contig_sizes
    File bcf
    File bcf_csi
    String sample_barcode
    Int? phred_threshold
    Int? min_inform
    Int? num_cpus = 8
    Int? ram_gb = 96
    Int? disk_size_gb = 500

    command {
        set -euo pipefail
        CHROM_SIZES_FILENAME=$(jq -r '.config.DEFAULT.reference |
            split("/")[-1] |
            rtrimstr(".gz") |
            rtrimstr(".fa") |
            rtrimstr(".fasta") +
            ".contig.sizes"' \
            ${gemBS_json})
        mkdir reference && ln ${reference} reference
        mkdir indexes && mv ${contig_sizes} indexes/"$CHROM_SIZES_FILENAME"
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

        # UCSC validation code does not work with BED track definition line.
        for i in extract/**/*.bed.gz; do
            FILENAME=$(basename "$i" .bed.gz)
            gzip -dck "$i" |
                tail -n +2 |
                gzip -n > "$FILENAME"_no_header.bed.gz
        done
    }

    output {
        File plus_strand_bw = glob("extract/**/*_pos.bw")[0]
        File minus_strand_bw = glob("extract/**/*_neg.bw")[0]
        File chg_bb = glob("extract/**/*_chg.bb")[0]
        File chh_bb = glob("extract/**/*_chh.bb")[0]
        File cpg_bb = glob("extract/**/*_cpg.bb")[0]
        File chg_bed = glob("extract/**/*_chg.bed.gz")[0]
        File chh_bed = glob("extract/**/*_chh.bed.gz")[0]
        File cpg_bed = glob("extract/**/*_cpg.bed.gz")[0]
        File chg_bed_no_header = glob("*_chg_no_header.bed.gz")[0]
        File chh_bed_no_header = glob("*_chh_no_header.bed.gz")[0]
        File cpg_bed_no_header = glob("*_cpg_no_header.bed.gz")[0]
        File cpg_txt = glob("extract/**/*_cpg.txt.gz")[0]  # gemBS-style output bed file
        File cpg_txt_tbi = glob("extract/**/*_cpg.txt.gz.tbi")[0]
        File non_cpg_txt = glob("extract/**/*_non_cpg.txt.gz")[0]
        File non_cpg_txt_tbi = glob("extract/**/*_non_cpg.txt.gz.tbi")[0]
    }

    runtime {
        cpu: "${num_cpus}"
        memory: "${ram_gb} GB"
        disks: "local-disk ${disk_size_gb} SSD"
    }
}

task make_coverage_bigwig {
    File encode_bed
    File chrom_sizes
    Int? num_cpus = 4
    Int? ram_gb = 8
    Int? disk_size_gb = 50

    command {
        set -euo pipefail
        BEDGRAPH_FILE="coverage.bedGraph"
        echo "track type=bedGraph" > "$BEDGRAPH_FILE"
        gzip -dc ${encode_bed} |
            tail -n +2 |
            xsv select -d '\t' -n 1,2,3,10 |
            xsv fmt -t '\t' >> "$BEDGRAPH_FILE"
        bedGraphToBigWig "$BEDGRAPH_FILE" ${chrom_sizes} coverage.bw
    }

    output {
        File coverage_bigwig = "coverage.bw"
    }

    runtime {
        cpu: "${num_cpus}"
        memory: "${ram_gb} GB"
        disks: "local-disk ${disk_size_gb} SSD"
    }
}

task calculate_bed_pearson_correlation {
    File bed1
    File bed2
    Int? num_cpus = 1
    Int? ram_gb = 8
    Int? disk_size_gb = 50

    command {
        python3 "$(which calculate_bed_pearson_correlation.py)" --bedmethyls ${bed1} ${bed2} --outfile bed_pearson_correlation_qc.json
    }

    output {
        File bed_pearson_correlation_qc = "bed_pearson_correlation_qc.json"
    }

    runtime {
        cpu: "${num_cpus}"
        memory: "${ram_gb} GB"
        disks: "local-disk ${disk_size_gb} SSD"
    }
}

task qc_report {
    Array[File] map_qc_json
    File reference
    File gemBS_json
    File contig_sizes
    String sample_barcode
    Int? num_cpus = 1
    Int? ram_gb = 4
    Int? disk_size_gb = 50

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
        File portal_map_qc_json = "gembs_map_qc.json"
        Array[File] map_qc_insert_size_plot_png = glob("mapping_reports/mapping/${sample_barcode}.isize.png")
        File map_qc_mapq_plot_png = "mapping_reports/mapping/${sample_barcode}.mapq.png"
    }

    runtime {
        cpu: "${num_cpus}"
        memory: "${ram_gb} GB"
        disks: "local-disk ${disk_size_gb} HDD"
    }
}
