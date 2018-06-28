workflow bscall {
	File reference_fasta
	File bam
	File bai
	String organism
	String sample
	Array[String] chromosomes

	scatter(chromosome in chromosomes) {
		call bscall_job { input:
			organism = organism,
			reference_fasta = reference_fasta,
			bam = bam,
			bai = bai,
			sample = sample,
			chromosome = chromosome
		}
	}

	output {
		Array[File] bcf_files = bscall_job.bcf
		Array[File] json_files = bscall_job.json
	}

}

task bscall_job {
	File reference_fasta
	File bam
	File bai
	String organism
	String sample
	String chromosome

	Int? memory_gb
	Int? cpu
	String? disks

	
	command {
		mkdir -p data/chr_snp_calls
		mkdir -p data/sample_mappings
		ln -s ${bam} data/sample_mappings/
		ln -s ${bai} data/sample_mappings/
		gemBS bscall \
			-r ${reference_fasta} \
			-e ${organism} \
			-s ${sample} \
			-c ${chromosome} \
			-i data/sample_mappings/$(basename ${bam}) \
			-o data/chr_snp_calls
	}

	output {
		File json = glob("data/chr_snp_calls/*.json")[0]
		File bcf = glob("data/chr_snp_calls/*.bcf")[0]
	}

	runtime {
		cpu: select_first([cpu,8])
		memory : "${select_first([memory_gb,'30'])} GB"
		disks : select_first([disks,"local-disk 200 HDD"])
		preemptible: 0
	}
}
