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
		gemBS index -i ${reference_fasta}
	}

	output {
		#Array[String] reference_bs = read_lines(stdout())
		File reference_gem = glob("./*.BS.gem")[0]
		File reference_info = glob("./*.BS.info")[0]
	}
}