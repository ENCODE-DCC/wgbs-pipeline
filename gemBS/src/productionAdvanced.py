#!/usr/bin/env python
"""Production pipelines"""
import os
import logging
import src
from utils import CommandException
import production

class DirectMapping(production.BasicPipeline):
    title = "Bisulphite mapping from a given file or standard input"
    description = """Maps single end or paired end bisulfite sequence using the gem mapper from a given file or standard input. 
                      
                  Exemple:
                      1) Given Interleaved FASTQ:
                      gemBS direct-mapping -I ref.BS.gem --fli flowcellname_lanename_indexname --json myfile.json -i interleaved.fq.gz --output-dir OUTPUTPATH --tmp-dir $TMPDIR --threads 8 -p
                      2) Standard Input FASTQ:
                      zcat interleaved.fq.gz | gemBS direct-mapping -I ref.BS.gem --fli flowcellname_lanename_indexname --json myfile.json --output-dir OUTPUTPATH --tmp-dir $TMPDIR --threads 8 -p
                      3) Single End FASTQ:
                      gemBS direct-mapping -I ref.BS.gem --fli flowcellname_lanename_indexname --json myfile.json -i singlend.fq.gz --output-dir OUTPUTPATH --tmp-dir $TMPDIR --threads 8 
                  """   
                                    
    def register(self,parser):
        ## required parameters
        parser.add_argument('-I', '--index', dest="index", metavar="index_file.BS.gem", help='Path to the Bisulfite Index Reference file.', required=True)
        parser.add_argument('-f', '--fli', dest="fli", metavar="FLOWCELL_LANE_INDEX", help='Lane to be mapped, format must be FLOWCELL_LANE_INDEX.', required=True)
        parser.add_argument('-j', '--json', dest="json_file", metavar="JSON_FILE", help='JSON file configuration.', required=True)
        parser.add_argument('-i', '--input-file', dest="input_file", metavar="FILE", help='Input file to be aligned.(FASTA/FASTQ/BAM).Default: Standard Input.', required=False)
        parser.add_argument('-1', '--i1', dest="pair_one", metavar="FILE", help='Pair One FASTQ file.', required=False)
        parser.add_argument('-2', '--i2', dest="pair_two", metavar="FILE", help='Pair Two FASTQ file.', required=False)
        parser.add_argument('-o', '--output-dir', dest="output_dir", metavar="PATH",default=".", help='Directory to store Bisulfite mapping results. Default: %s' %self.output_dir)
        parser.add_argument('-d', '--tmp-dir', dest="tmp_dir", metavar="PATH", default="/tmp/", help='Temporary folder to perform sorting operations. Default: %s' %self.tmp_dir)      
        parser.add_argument('-t', '--threads', dest="threads",default="1", help='Number of threads to perform sorting operations. Default %s' %self.threads)
        parser.add_argument('-b', '--is-bam', dest="is_bam", action="store_true", default=False, help="Input data (File or STDIN) is a bam or sam file.")        
        parser.add_argument('-p', '--paired-end', dest="paired_end", action="store_true", default=False, help="Input data is Paired End")
        parser.add_argument('-s', '--read-non-stranded', dest="read_non_stranded", action="store_true",default=False, 
                            help='Automatically selects the proper C->T and G->A read conversions based on the level of Cs and Gs on the read.')     
        parser.add_argument('-n', '--underconversion-sequence', dest="underconversion_sequence", metavar="SEQUENCE", help='Name of Lambda Sequence used to control unmethylated cytosines which fails to be\
                             deaminated and thus appears to be Methylated.', default=None,required=False)
        parser.add_argument('-v', '--overconversion-sequence', dest="overconversion_sequence", metavar="SEQUENCE", help='Name of Lambda Sequence used to control methylated cytosines which are\
                             deaminated and thus appears to be Unmethylated.', default=None,required=False)
                             
                             
    def run(self, args):     
        #Flowcell Lane Index
        self.name = args.fli
        #Flowcell Lane Index Info
        self.fliInfo = production.FLIdata(args.json_file).sampleData[self.name] 
        #Index
        self.index = args.index
        #Output Directory
        self.output_dir = args.output_dir
        #Paired
        self.paired = args.paired_end
        #Read Non Standard
        self.read_non_stranded = args.read_non_stranded
        #TMP
        self.tmp_dir = args.tmp_dir 
        #Threads
        self.threads = args.threads
        #Inputs
        self.pair_one = args.pair_one
        self.pair_two = args.pair_two
        self.input_file = args.input_file
        self.is_bam = args.is_bam
                
        #Check for input data
        if args.paired_end:
            if args.pair_one and args.pair_two and args.input_file:
                raise CommandException("Sorry!! Pair One and Pair Two arguments can not be used with Input file Argument." %(args.pair_one))
            if (args.pair_one and not args.pair_two) or ( not args.pair_one and args.pair_two):
                raise CommandException("Sorry!! Pair One and Pair Two arguments must be specified together." %(args.pair_one))  
            if args.pair_one:
                if not os.path.isfile(args.pair_one):
                    raise CommandException("Sorry!! Pair one file %s was not found." %(args.pair_one))
            if args.pair_two:
                if not os.path.isfile(args.pair_two):
                    raise CommandException("Sorry!! Pair two file %s was not found." %(args.pair_two))
            if args.input_file:
                if not os.path.isfile(args.input_file):
                    raise CommandException("Sorry!! File %s was not found." %(args.input_file))
        else:
            if args.pair_one or args.pair_two:
                raise CommandException("Sorry!! Pair One and Pair Two arguments can not be used in single end mode.")
            if args.input_file:
                if not os.path.isfile(args.input_file):
                    raise CommandException("Sorry!! File %s was not found." %(args.input_file))
        
        if args.is_bam and (args.pair_one or args.pair_two):
            raise CommandException("Sorry!! Bam argument can not be present if Pair one or Pair Two FASTQ files are determined.") 
            
        #Check Bisulfite Conversion
        self.underconversion_sequence = ""
        if args.underconversion_sequence is not None:
            self.underconversion_sequence = args.underconversion_sequence
            
        self.overconversion_sequence = ""
        if args.overconversion_sequence is not None: 
            self.overconversion_sequence = args.overconversion_sequence
            
        #Check Temp Directory
        if not os.path.isdir(self.tmp_dir):
            raise CommandException("Temporary directory %s does not exists or is not a directory." %(self.tmp_dir))  

        self.log_parameter()
        logging.gemBS.gt("Bisulfite Direct Mapping...")      

        ret = src.direct_mapping(name=self.name,index=self.index,fliInfo=self.fliInfo,paired=self.paired,threads=self.threads,
                   file_pe_one=self.pair_one,file_pe_two=self.pair_two,file_input=self.input_file,is_bam=self.is_bam,
                   read_non_stranded=self.read_non_stranded,
                   outputDir=self.output_dir,tmpDir=self.tmp_dir,
                   under_conversion=self.underconversion_sequence,over_conversion=self.overconversion_sequence)

        if ret:
            logging.gemBS.gt("Bisulfite Mapping done!! Output File: %s" %(ret))

    def extra_log(self):
        """Extra Parameters to be printed"""
        #Virtual methos, to be define in child class
        printer = logging.gemBS.gt
        
        printer("------------ Mappings Parameter ------------")
        printer("Name             : %s", self.name)
        printer("Index            : %s", self.index)
        printer("Paired           : %s", self.paired)
        printer("Read non stranded: %s", self.read_non_stranded)
        printer("BAM/SAM input    : %s", self.is_bam)
    
        if self.paired:
            if self.pair_one:
                printer("Input Pair One   : %s", self.pair_one)
                printer("Input Pair Two   : %s", self.pair_two)
            elif self.input_file:
                printer("Input Interleaved: %s", self.pair_one)
            else:
                printer("Input Paired End : STDIN")
        else:
            if self.input_file:
                printer("Input Single End : %s", self.pair_one)
            else:
                printer("Input Single End : STDIN")
                        
        printer("")





         