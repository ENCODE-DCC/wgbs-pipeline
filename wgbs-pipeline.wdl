# ENCODE WGBS Pipeline running https://github.com/heathsc/gemBS
# Maintainer: Ulugbek Baymuradov

workflow wgbs {

	File reference_fasta
	File metadata 
	File sequence
	
	call index { input:
		reference_fasta = reference_fasta
	} 
}

task index {

	File reference_fasta

	command {
		gemBS index -i reference
	}

	output {
		File reference_gem
		File reference_info
	}

	runtime {
		docker : "quay.io/encode-dcc/wgbs-pipeline:latest"
	}
}