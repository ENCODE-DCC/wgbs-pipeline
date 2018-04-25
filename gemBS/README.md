GEMBS
=====

GEMBS is a bioinformatic pipeline designed for hightroughput analysis of DNA methylation from whole genome bisulfites sequencing data (WGBS). It implements GEM3, a high performance read aligner and BScall, a variant caller specifically for bisulfite sequencing data.


The manuscript describing the pipeline is available [here](https://www.biorxiv.org/content/early/2017/10/11/201988)

---------   
Licensing
---------

GEMBS is licensed under GPL.

--------
Download
--------

Use ``git clone --recursive`` to retrieve the complete source code including the code from external projects such as ``bs_call`` and ``gem3-mapper``.

    git clone --recursive https://github.com/heathsc/gemBS.git

------------
Installation
------------

1) Before starting the installation of GEMBS, please check if your system has the GSL library already installed.

    If your system does not have GSL library then you can download it from [GSL](https://www.gnu.org/software/gsl/) and follow the installation steps. 

    Once GSL is already available on your system then you can compile and install GEMBS.

2) Change GSL library paths. In order to compile bscall you must specify the GSL headers and library directories. 
   To do that, edit file ./tools/bs_call/Gsl.mk file with the proper paths. Just modify two lines starting with GSL_LIB and GSL_INC.

    ./tools/bs_call/Gsl.mk:  

        #1. MODIFY HERE THE GSL LIBRARY LOCATION. FOR Example: GSL_LIB = -L/path/to/gsl/lib
        GSL_LIB = -L/path/to/GSL/lib/
        #2. MODIFY HERE THE GSL HEADERS LOCATION. FOR Example: GSL_LIB = -L/path/to/gsl/include
        GSL_INC = -I/path/to/GSL/include/ 

3) Use python install command:

    ``python setup.py install --user``

-----------------------
Check your installation
-----------------------

For checking your installation follow the commands found at [README_TEST.md](test/README_TEST.md).


-------------
Documentation
-------------

Documentation can be found at [GEMBS](http://statgen.cnag.cat/GEMBS/)

----------
Changelog:
----------
    2.1.0 Add new command gemBS direct-mapping
    2.1.0 Fixed Makefile clean in tools
    2.0.2 Fixed bug related with the percentage of High Quality Variant in Variants summary report.
    2.0.2 Check temporary directory existence.
    2.0.2 Fixed QualityNonRefCpg sample name in png image.
    2.0.2 Fixed mapper issues related with aligning performace.
    2.0.2 Fixed arguments for Under/Over Conversion sequence name in gem3-mapper
    2.0.1 On bscall repository, fixed argument -k about discarded reads that do not form proper pairs.
    2.0 Check tmp folder before starting mapping process.
    2.0 Added Left and Right Trimming optional arguments to gemBS bscall.
    2.0 Added GC Coverage correlation value to BS Call Stats Summary.
    2.0 Fixed error when reporting complete path to not found bam files.
    2.0 Fixed iteration over sampleBams dictionary in MergeAll method.
    2.0 Updated: Avoid redo indexing when merging just one file.
    2.0 Changed conversion formula.
    2.0 Added parameter for dbSNP.
    2.0 Added threads to bscall.
    2.0 Removed CpGs reports. Already done from bscall report.
    2.0 Fixed bs_call makefile for the gcc to be used.
    2.0 New bscall version. Generates JSON report.
    2.0 Removed gemBS options snp-stats,cpg-report,cpg-stats.
    2.0 Added summary report from the bs_call json stats
    2.0 New BSCall Report. From bscall son file generates three types of reports:
        Mapping and Coverage Report
        Bs-Genotypes Calls Report
        Methylation Statistics report
    1.7 Added non stranded read conversion parameter
    1.7 Fixed SE crash when estimating overlapped bases.
    1.7 Fixed gem-index (gem3) to follow fastq and SAM specifications. 
        Modified gem3-mapper repository external module.
        New external module https://github.com/heathsc/gem3-mapper.git
    1.7 Fixed threads parameter to samtools merge
    1.7 Fixed threads parameter to gem-mapper
    1.7 Removed Indels Field on Variants Report.
    1.7 Added Overlapping Bases at Mapping Report
    1.7 Modified Base Counts Overall, removed Base Counts general and Base Counts Overall
    1.7 New Dinucleotide CpGs Report
        New table dinucleotide stats
        New plots for Informative Reads and CpGs
        Methylation levels plots for different types of CpGs
        Summary Table
    1.7 New Readme file to inform about report test
    1.7 New basic statis table for Variants Report
    1.7 Removed parameter -r (reference length) parameter for mapping reports command (gemBS bsMap).
    1.6 New CpGs Density plot, include box plos, bar plot and fitting curve
    1.6 Change name at CpG report:
        "Heterozygous" for "Alternative CX"
        "De Novo CpGs Methylation Status" for "Non Reference CpGs"
        "CpGs with SNP" for "SNPs (CX) at Reference CpGs"
    1.6 CpGs Report Simplified to Q>20
    1.6 BigWig Default parameters for filtering CpG per a given quality and a total number of supported informative reads   
    1.5 Initial Release  


----------
Developers
----------
 
 GEMBS:
 * Marcos Fernandez-Callejo - marcos.fernandez@cnag.crg.eu
 
 gem mapper:
 * Santiago Marco-Sola - santiago.marco@cnag.crg.eu

 bisulfite caller and filtering:
 * Simon Heath - simon.heath@cnag.crg.eu


