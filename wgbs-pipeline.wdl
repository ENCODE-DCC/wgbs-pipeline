# ENCODE WGBS Pipeline running https://github.com/heathsc/gemBS
# Maintainer: Ulugbek Baymuradov
import 'index.wdl' as step1

workflow wgbs {

	File reference_fasta
	File metadata 
	File sequence
	
	call step1.index { input:
		reference_fasta = reference_fasta
	}

	output {
		File reference_gem = index.reference_gem
		File reference_info = index.reference_info
	}
}