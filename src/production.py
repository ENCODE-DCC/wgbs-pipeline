#!/usr/bin/env python
"""Production pipelines"""
import os
import logging
import json
import sys
from sys import exit
import subprocess

from utils import Command, CommandException

import src
import report
import sphinx
import bsCallReports


class Fli(object):
    
    def __init__(self):
        #fli Members
        self.sample_barcode = None
        self.library = None
        self.flowcell = None
        self.lane = None
        self.index = None
        
    def getFli(self):
        #Get FLI (flowcell lane index)
        return "%s_%s_%s" %(self.flowcell,self.lane,self.index)


class FLIdata(object):
     #Class to manage the flowcell lane index information of the project
     def __init__(self,json_file=None):
        self.json_file = json_file
        self.sampleData = {}
   
        with open(self.json_file, 'r') as fileJson:
            data = json.load(fileJson)      
        
            for fli in data:
                fliCommands = Fli()            
                
                for key, value in data[fli].iteritems():
                    if key == "sample_barcode":
                        fliCommands.sample_barcode = value    
                    elif key == "library_barcode":
                        fliCommands.library = value    
                    elif key == "flowcell_name":
                        fliCommands.flowcell = value
                    elif key == "lane_number":
                        fliCommands.lane = value
                    elif key == "index_name":
                        fliCommands.index = value
                
                self.sampleData[fliCommands.getFli()] = fliCommands

class BasicPipeline(Command):
    """General mapping pipeline class."""

    def __init__(self):
        # general parameter
        self.input = None # input files
        self.output = None #Output files
        self.output_dir = "."
        self.tmp_dir = "/tmp/"
        self.threads = "1"
        
        self.membersInitiation()
        
    def membersInitiation(self):
        #To fullfill in the child class
        pass
        
    def log_parameter(self):
        """Print selected parameters"""
        printer = logging.gemBS.gt
      
        printer("------------ Input Parameter ------------")
        printer("Input File(s)    : %s", self.input)
        printer("Output File(s)   : %s", self.output)
        printer("Output Directory : %s", self.output_dir)
        printer("TMP Directory    : %s", self.tmp_dir)
        printer("Threads          : %s", self.threads)
        printer("")
       
        self.extra_log()
        
    def extra_log(self):
        """Extra Parameters to be printed"""
        #Virtual methos, to be define in child class
        pass


class PrepareConfiguration(Command):
    title = "Prepare Configuration"
    description = """ Creates a json configuration file to perform the different steps of the Bisulfite pipeline.
                      
                      Input Files:
                          Option 1: Simple text file, comma separated with 5 columns.
                          FORMAT: sample_id,library,flowcell,lane,index
                          
                          Option 2: CNAG Lims Subproject json file
                    
                      If you are managing CNAG bisulfite sequencing data and you have acces to the CNAG lims then Option 2 is the most user friendly.
                      Otherwise Option 1.
                  """
                
    def register(self, parser):
        ## required parameters
        parser.add_argument('-t', '--text-metadata', dest="text_metadata", help="""Metadata configured in 5 columns text file. Comma separated values (CSV) file.
                                                                                   FORMAT:  sample_id,library,flowcell,lane,index                                                                           
                                                                                """,default=None)
        parser.add_argument('-l', '--lims-cnag-json', dest="lims_cnag_json", help="""Lims Cnag subproject json file.""",default=None)
        parser.add_argument('-j', '--json', dest="json", help='JSON ouput file',required=True)
        
    def run(self,args):        
        #Try text metadata file
        if args.text_metadata is not None:
            if os.path.isfile(args.text_metadata):
                src.prepareConfiguration(text_metadata=args.text_metadata,jsonOutput=args.json)
            else:
                raise CommandException("Sorry!! File %s not found!" %(args.text_metadata))
        elif args.lims_cnag_json is not None:
            if os.path.isfile(args.lims_cnag_json):
                src.prepareConfiguration(lims_cnag_json=args.lims_cnag_json,jsonOutput=args.json)
            else:
                raise CommandException("Sorry!! File %s not found!" %(args.lims_cnag_json))
        else:
            raise CommandException("No input file inserted!!")
                    
     
class Index(BasicPipeline):
    title = "Index genomes"
    description = """Reference indexing for Bisulfite GEM mapping 
                     Generates a file called reference.BS.gem (Index) and 
                     reference.BS.info (Information about the index process) 
    """

    def register(self, parser):
        ## required parameters
        parser.add_argument('-i', '--input', dest="input", help='Path to a single fasta reference genome file.', required=True)
        parser.add_argument('-t', '--threads', dest="threads", help='Number of threads. By default GEM indexer will use the maximum available on the system.',default=None)
        parser.add_argument('-d','--list-dbSNP-files',dest="list_db_snp_files",nargs="+",metavar="FILES",
                            help="List of dbSNP files (can be compressed) to create an index to later use it at the bscall step. The bed files should have the name of the SNP in column 4.",default=[],required=False)
        parser.add_argument('-x', '--dbsnp-index', dest="dbsnp_index", help='dbSNP output index file name.',default="",required=False)

    def run(self, args):
        self.input = args.input
        self.threads = args.threads
        self.list_dbSNP_files = args.list_db_snp_files
        self.dbsnp_index = args.dbsnp_index
        if len(self.list_dbSNP_files)>0:     
            if self.dbsnp_index == "":
                raise CommandException("dbSNP Index file must be specified through --dbsnp-index parameter.")
        
        if not os.path.exists(self.input):
            raise CommandException("Input file not found : %s" % self.input)
        
        if self.input.endswith(".fasta.gz"):
            self.output = "%s.BS" %(self.input[:-9])
        elif self.input.endswith(".fa.gz"):
            self.output = "%s.BS" %(self.input[:-6])
        elif self.input.endswith(".fasta"):
            self.output = "%s.BS" %(self.input[:-6])    
        elif self.input.endswith(".fa"):
            self.output = "%s.BS" %(self.input[:-3])
        else:
            raise CommandException("Sorry!! Input file %s should be a Fasta file with one of these suffixes: .fa .fasta .fa.gz fasta.gz.")
        
        self.log_parameter()
        logging.gemBS.gt("Creating index")
        ret = src.index(self.input, self.output, threads=self.threads,list_dbSNP_files=self.list_dbSNP_files,dbsnp_index=self.dbsnp_index)
        if ret:
            logging.gemBS.gt("Index done: %s.gem" %(ret))
            
       
class MappingCommnads(BasicPipeline):
    title = "Show Mapping commands"
    description = """ From a json input file, generates the set of mapping commands to run for mapping all Bisulfite data involved in a Project """

    def register(self,parser):
        ## required parameters
        parser.add_argument('-I', '--index', dest="index", metavar="index_file.BS.gem", help='Path to the Bisulfite Index Reference file.', required=True)
        parser.add_argument('-j', '--json', dest="json_file", metavar="JSON_FILE", help='JSON file configuration.', required=True)
        parser.add_argument('-i', '--input-dir', dest="input_dir", metavar="PATH", help='Directory where is located input data. FASTQ or BAM format.', required=True)
        parser.add_argument('-o', '--output-dir', dest="ouput_dir", metavar="PATH",default=".", help='Directory to store Bisulfite mapping results. Default: %s' %self.output_dir)
        parser.add_argument('-d', '--tmp-dir', dest="tmp_dir", metavar="PATH", default="/tmp/", help='Temporary folder to perform sorting operations. Default: %s' %self.tmp_dir)      
        parser.add_argument('-t', '--threads', dest="threads",default="1", help='Number of threads to perform sorting operations. Default: %s' %self.threads)
        parser.add_argument('-p', '--paired-end', dest="paired_end", action="store_true", default=None, help="Input data is Paired End")
        parser.add_argument('-s', '--read-non-stranded', dest="read_non_stranded", action="store_true",default=False, 
                            help='Automatically selects the proper C->T and G->A read conversions based on the level of Cs and Gs on the read.') 
        parser.add_argument('-n', '--underconversion-sequence', dest="underconversion_sequence", metavar="SEQUENCE", help='Name of Lambda Sequence used to control unmethylated cytosines which fails to be\
                             deaminated and thus appears to be Methylated.', default=None, required=False)
        parser.add_argument('-v', '--overconversion-sequence', dest="overconversion_sequence", metavar="SEQUENCE", help='Name of Lambda Sequence used to control methylated cytosines which are\
                             deaminated and thus appears to be Unmethylated.', default=None, required=False)
        
        
    def run(self,args):
        ## All Flowcell Lane Index        
        for k,v in FLIdata(args.json_file).sampleData.iteritems():
            ##Non Stranded
            non_stranded = ""
            if args.read_non_stranded:
                non_stranded += "-s"
            ## Conversion Parameters
            conversion_parameters = ""            
            if args.underconversion_sequence is not None:  
                conversion_parameters += " -n %s" %(args.underconversion_sequence)
            if args.overconversion_sequence is not None:  
                conversion_parameters += " -v %s" %(args.overconversion_sequence)
            
            
            if args.paired_end:
                print "gemBS mapping -I %s -f %s -j %s -i %s -o %s -d %s -t %s -p %s %s"\
                       %(args.index,k,args.json_file,args.input_dir,args.ouput_dir,args.tmp_dir,str(args.threads),non_stranded,conversion_parameters)
            else:
                print "gemBS mapping -I %s -f %s -j %s -i %s -o %s -d %s -t %s %s %s"\
                      %(args.index,k,args.json_file,args.input_dir,args.ouput_dir,args.tmp_dir,str(args.threads),non_stranded,conversion_parameters)
            
        
class Mapping(BasicPipeline):
    title = "Bisulphite mapping"
    description = """Maps a single end or paired end bisulfite sequence using the gem mapper. 
    
                  Each time a mapping is called a fastq file is mapped. (Two Paired fastq files in case of paired end).
                  Files must be located in an input directory in the form: FLOWCELL_LANE_INDEX.suffix
                  
                  Suffix could be: _1.fq _2.fq _1.fastq _2.fastq _1.fq.gz _2.fq.gz _1.fastq.gz _2.fastq.gz for paired end.
                                   .fq .fq.gz .fastq .fastq.gz for single end or interleaved paired end files.

                  Suffix could also be .bam in case of aligned files.
                  
                  Exemple:
                      gemBS mapping -I ref.BS.gem --fli flowcellname_lanename_indexname --json myfile.json --input-dir INPUTPATH --output-dir OUTPUTPATH --tmp-dir $TMPDIR --threads 8 -p
        
                  """   
 
    def register(self,parser):
        ## required parameters
        parser.add_argument('-I', '--index', dest="index", metavar="index_file.BS.gem", help='Path to the Bisulfite Index Reference file.', required=True)
        parser.add_argument('-f', '--fli', dest="fli", metavar="FLOWCELL_LANE_INDEX", help='Lane to be mapped, format must be FLOWCELL_LANE_INDEX.', required=True)
        parser.add_argument('-j', '--json', dest="json_file", metavar="JSON_FILE", help='JSON file configuration.', required=True)
        parser.add_argument('-i', '--input-dir', dest="input_dir", metavar="PATH", help='Directory where is located input data. FASTQ or BAM format.', required=True)
        parser.add_argument('-o', '--output-dir', dest="output_dir", metavar="PATH",default=".", help='Directory to store Bisulfite mapping results. Default: %s' %self.output_dir)
        parser.add_argument('-d', '--tmp-dir', dest="tmp_dir", metavar="PATH", default="/tmp/", help='Temporary folder to perform sorting operations. Default: %s' %self.tmp_dir)      
        parser.add_argument('-t', '--threads', dest="threads",default="1", help='Number of threads to perform sorting operations. Default %s' %self.threads)
        parser.add_argument('-p', '--paired-end', dest="paired_end", action="store_true", default=None, help="Input data is Paired End")
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
        self.fliInfo = FLIdata(args.json_file).sampleData[self.name] 
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
        
        #Input Data
        self.input_pair_one = None
        self.input_pair_two = None
        self.input_interleaved = None        
        self.input_se = None
        self.input_bam = None        
        
        self.inputPath = args.input_dir + "/" + self.fliInfo.getFli()
        
        #Check for input data
        if args.paired_end:
            #Pair One
            if os.path.isfile(self.inputPath + "_1.fastq"):
                self.input_pair_one = self.inputPath + "_1.fastq"
            elif os.path.isfile(self.inputPath + "_1.fastq.gz"):
                self.input_pair_one = self.inputPath + "_1.fastq.gz"
            elif os.path.isfile(self.inputPath + "_1.fq"):
                self.input_pair_one = self.inputPath + "_1.fq"
            elif os.path.isfile(self.inputPath + "_1.fq.gz"):
                self.input_pair_one = self.inputPath + "_1.fq.gz"
            
            #Pair Two
            if os.path.isfile(self.inputPath + "_2.fastq"):
                self.input_pair_two = self.inputPath + "_2.fastq"
            elif os.path.isfile(self.inputPath + "_2.fastq.gz"):
                self.input_pair_two = self.inputPath + "_2.fastq.gz"
            elif os.path.isfile(self.inputPath + "_2.fq"):
                self.input_pair_two = self.inputPath + "_2.fq"
            elif os.path.isfile(self.inputPath + "_2.fq.gz"):
                self.input_pair_two = self.inputPath + "_2.fq.gz"
                
            #Interleaved
            if self.input_pair_one is None and self.input_pair_two is None:
                if os.path.isfile(self.inputPath + ".fastq"):
                    self.input_interleaved = self.inputPath + ".fastq"
                elif os.path.isfile(self.inputPath + ".fastq.gz"):
                    self.input_interleaved = self.inputPath + ".fastq.gz"
                elif os.path.isfile(self.inputPath + ".fq"):
                    self.input_interleaved = self.inputPath + ".fq"
                elif os.path.isfile(self.inputPath + ".fq.gz"):
                    self.input_interleaved = self.inputPath + ".fq.gz"
        else:
            #Single End
                if os.path.isfile(self.inputPath + ".fastq"):
                    self.input_se = self.inputPath + ".fastq"
                elif os.path.isfile(self.inputPath + ".fastq.gz"):
                    self.input_se = self.inputPath + ".fastq.gz"
                elif os.path.isfile(self.inputPath + ".fq"):
                    self.input_se = self.inputPath + ".fq"
                elif os.path.isfile(self.inputPath + ".fq.gz"):
                    self.input_se =self.inputPath + ".fq.gz"  
                    
        #Check for BAM input data
        if self.input_pair_one is None and self.input_pair_two is None and self.input_interleaved is None and self.input_se is None:
            if os.path.isfile(self.inputPath + ".bam"):
                self.input_bam = self.inputPath + ".bam"
                
        #Check for input existance
        if self.input_pair_one is None and self.input_pair_two is None and self.input_interleaved is None and self.input_se is None and self.input_bam is None: 
            raise CommandException("No input files where found in %s directory." %(args.input_dir))
            
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
        logging.gemBS.gt("Bisulfite Mapping...")
        ret = src.mapping(name=self.name,index=self.index,fliInfo=self.fliInfo,
                    file_pe_one=self.input_pair_one,file_pe_two=self.input_pair_two,
                    file_interleaved = self.input_interleaved,
                    file_se = self.input_se,
                    read_non_stranded = self.read_non_stranded,
                    file_bam = self.input_bam,
                    outputDir=self.output_dir,paired = self.paired,
                    tmpDir=self.tmp_dir,threads=self.threads,
                    under_conversion=self.underconversion_sequence,
                    over_conversion=self.overconversion_sequence)
                    
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
        printer("Input Pair One   : %s", self.input_pair_one)
        printer("Input Pair Two   : %s", self.input_pair_two)
        printer("Input Interleaved: %s", self.input_interleaved)
        printer("Input Single End : %s", self.input_se)        
        printer("Input BAM        : %s", self.input_bam)
        
        
        printer("")
        
        

class MergingAll(BasicPipeline):
    title = "Merging bams"
    description = """Merges all bam alignments involved in a given Bisulfite project.
                     Each bam alignment file belonging to a sample should be merged to perform the methylation calling."""
                     
    def register(self,parser):
        ## required parameters                     
        parser.add_argument('-i', '--input-dir', dest="input_dir",metavar="PATH", help='Path where are located the BAM aligned files.', required=True)
        parser.add_argument('-j', '--json', dest="json_file", metavar="JSON_FILE", help='JSON file configuration.', required=True)
        parser.add_argument('-t', '--threads', dest="threads", metavar="THREADS", default="1", help='Number of threads, Default: %s' %self.threads)
        parser.add_argument('-o', '--output-dir', dest="output_dir", metavar="PATH",help='Output directory to store merged results.',required=True)
        parser.add_argument('-d', '--tmp-dir', dest="tmp_dir", metavar="PATH", default="/tmp/", help='Temporary folder to perform sorting operations. Default: %s' %self.tmp_dir)
        
    def run(self, args):
        #Threads
        self.threads = args.threads
        #Output Directory
        self.output_dir = args.output_dir
        #TMP DIR
        self.tmp_dir = args.tmp_dir
        
        #Create Dictionary of samples and bam file        
        self.samplesBams = {}  
        self.records = 0
        for k,v in FLIdata(args.json_file).sampleData.iteritems():
            fileBam = "%s/%s.bam" %(args.input_dir,v.getFli())
            self.records = self.records + 1
            if os.path.isfile(fileBam):
                if v.sample_barcode not in self.samplesBams:
                    self.samplesBams[v.sample_barcode] = [fileBam]
                else:
                    self.samplesBams[v.sample_barcode].append(fileBam) 
                    
        #Check list of file
        self.totalFiles = 0
        for sample,listBams in self.samplesBams.iteritems():
            self.totalFiles += len(listBams)
              
        if self.totalFiles < 1:
            raise CommandException("Sorry not bam files were found!!")
        elif self.totalFiles != self.records:
            raise CommandException("Sorry not all bam files were created!!")
            
        self.log_parameter()
        logging.gemBS.gt("Merging process started...")
        ret = src.merging(inputs=self.samplesBams,threads=self.threads,output_dir=self.output_dir,tmpDir=self.tmp_dir)
         
        if ret:
            logging.gemBS.gt("Merging process done!! Output files generated:")
            for sample,outputBam  in ret.iteritems():
                logging.gemBS.gt("%s: %s" %(sample, outputBam))
                    

class MergingSample(BasicPipeline):   
    title = "Merging bams"
    description = """Merges all bam alignments for a given Bisulfite Sample."""

    def register(self,parser):
        ## required parameters                     
        parser.add_argument('-i', '--input-dir', dest="input_dir",metavar="PATH", help='Path were are located the BAM aligned files.', required=True)
        parser.add_argument('-j', '--json', dest="json_file", metavar="JSON_FILE", help='JSON file configuration.', required=True)
        parser.add_argument('-s', '--sample-id',dest="sample_id",metavar="SAMPLE",help="Sample unique identificator") 
        parser.add_argument('-t', '--threads', dest="threads", metavar="THREADS", default="1", help='Number of threads, Default: %s' %self.threads)
        parser.add_argument('-o', '--output-dir', dest="output_dir", metavar="PATH",help='Output directory to store merged results.',required=True)
        parser.add_argument('-d', '--tmp-dir', dest="tmp_dir", metavar="PATH", default="/tmp/", help='Temporary folder to perform sorting operations. Default: %s' %self.tmp_dir) 
        
    def run(self, args):
        #Threads
        self.threads = args.threads
        #Output Directory
        self.output_dir = args.output_dir
        #sample
        self.sample_id = args.sample_id
        #TMP DIR
        self.tmp_dir = args.tmp_dir
        
        #Create Dictionary of samples and bam file
        self.samplesBams = {} 
        self.records = 0
        
        for k,v in FLIdata(args.json_file).sampleData.iteritems():
            if self.sample_id == v.sample_barcode:            
                fileBam = "%s/%s.bam" %(args.input_dir,v.getFli())
                self.records = self.records + 1
                if os.path.isfile(fileBam):
                    if v.sample_barcode not in self.samplesBams:
                        self.samplesBams[v.sample_barcode] = [fileBam]
                    else:
                        self.samplesBams[v.sample_barcode].append(fileBam) 
                        
        #Check list of file
        self.totalFiles = 0
        for sample,listBams in self.samplesBams.iteritems():
            self.totalFiles += len(listBams)
              
        if self.totalFiles < 1:
            raise CommandException("Sorry not bam files were found!!")
        elif self.totalFiles != self.records:
            raise CommandException("Sorry not all bam files were created!!")
            
        self.log_parameter()
        logging.gemBS.gt("Merging process started...")
        ret = src.merging(inputs=self.samplesBams,threads=self.threads,output_dir=self.output_dir,tmpDir=self.tmp_dir)
         
        if ret:
            logging.gemBS.gt("Merging process done!! Output files generated:")
            for sample,outputBam  in ret.iteritems():
                logging.gemBS.gt("%s: %s" %(sample, outputBam))
        
        
    def extra_log(self):
        """Extra Parameters to be printed"""
        #Virtual methos, to be define in child class
        printer = logging.gemBS.gt
        
        printer("------- Merging Sample Parameter -------")
        printer("Sample ID           : %s", self.sample_id)
        printer("")
                     


        
class MethylationCall(BasicPipeline):
    title = "Methylation Calling"
    description = """Performs a methylation calling from a bam aligned file.
                     This process is performed over a list of chromosomes in a sequentially way.
                     If you prefer to run the methylation calls in parallel you should consider bscall
                     command.
                  """
    def membersInitiation(self):
        self.species = "HomoSapiens"
        self.chroms = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

                                   
    def register(self, parser):
        ## required parameters
        parser.add_argument('-r','--fasta-reference',dest="fasta_reference",metavar="PATH",help="Path to the fasta reference file.",required=True)
        parser.add_argument('-e','--species',dest="species",metavar="SPECIES",default="HomoSapiens",help="Sample species name. Default: %s" %self.species)
        parser.add_argument('-j','--json',dest="json_file",metavar="JSON_FILE",help='JSON file configuration.')
        parser.add_argument('-p','--path-bam',dest="path_bam",metavar="PATH_BAM",help='Path where are stored sample BAM files.',default=None)
        parser.add_argument('-g','--right-trim', dest="right_trim", metavar="BASES",type=int, default=0, help='Bases to trim from right of read pair, Default: 0')
        parser.add_argument('-f','--left-trim', dest="left_trim", metavar="BASES",type=int, default=5, help='Bases to trim from left of read pair, Default: 5')        
        parser.add_argument('-o','--output-dir',dest="output_dir",metavar="PATH",help='Output directory to store the results.',default=None)
        parser.add_argument('-d','--paired-end', dest="paired_end", action="store_true", default=False, help="Input data is Paired End")
        parser.add_argument('-t','--threads', dest="threads", metavar="THREADS", default="1", help='Number of threads, Default: %s' %self.threads)
        parser.add_argument('-k','--keep-unmatched', dest="keep_unmatched", action="store_true", default=False, help="Do not discard reads that do not form proper pairs.")
        parser.add_argument('-u','--keep-duplicates', dest="keep_duplicates", action="store_true", default=False, help="Do not merge duplicate reads.")    
        parser.add_argument('-b','--dbSNP-index-file', dest="dbSNP_index_file", metavar="FILE", help="dbSNP index file.",required=False,default="")
        parser.add_argument('-l','--list-chroms',dest="list_chroms",nargs="+",metavar="CHROMS",help="""List of chromosomes to perform the methylation pipeline.
                                                                                                       Can be a file where every line is a chromosome contig. 
                                                                                                       By default human chromosomes: %s """ %self.chroms,
                            default=["chr1","chr2","chr3","chr4","chr5","chr6","chr7",
                                     "chr8","chr9","chr10","chr11","chr12","chr13",
                                     "chr14","chr15","chr16","chr17","chr18","chr19",
                                     "chr20","chr21","chr22","chrX","chrY"])
       
    def run(self,args):
        self.threads = args.threads
        self.fasta_reference = args.fasta_reference 
        self.species = args.species
        self.input_dir = args.path_bam
        self.right_trim = args.right_trim
        self.left_trim = args.left_trim
        self.json_file = args.json_file
        self.output_dir = args.output_dir  
        self.paired = args.paired_end
        self.keep_unmatched = args.keep_unmatched
        self.keep_duplicates = args.keep_duplicates
        self.dbSNP_index_file = args.dbSNP_index_file
        
        self.list_chroms = []
    
        if len(args.list_chroms) > 1:
            self.list_chroms = args.list_chroms
        elif os.path.isfile(args.list_chroms[0]):
            #Check if List_chroms is a file or just a list of chromosomes
            #Parse file to extract chromosme list 
            with open(args.list_chroms[0] , 'r') as chromFile:
                for line in chromFile:
                    self.list_chroms.append(line.rstrip())
        else:
            self.list_chroms = args.list_chroms          
        
        #Check fasta existance
        if not os.path.isfile(args.fasta_reference):
            raise CommandException("Sorry path %s was not found!!" %(args.fasta_reference))
        
        #Check input bam existance
        self.sampleBam = {}
        for k,v in FLIdata(args.json_file).sampleData.iteritems():
            fileBam = "%s/%s.bam" %(self.input_dir,v.sample_barcode)

            if not os.path.isfile(fileBam):
                raise CommandException("Sorry path %s was not found!!" %(fileBam))

            if v.sample_barcode not in self.sampleBam:
                self.sampleBam[v.sample_barcode] = fileBam
        
        #Call for everything
        self.log_parameter()
        logging.gemBS.gt("Methylation Calling...")
        if len(args.list_chroms) > 0:
            ret = src.methylationCalling(reference=self.fasta_reference,species=self.species,
                                         right_trim=self.right_trim, left_trim=self.left_trim,
                                         sample_bam=self.sampleBam,chrom_list=self.list_chroms,
                                         output_dir=self.output_dir,paired_end=self.paired,keep_unmatched=self.keep_unmatched,
                                         keep_duplicates=self.keep_duplicates,dbSNP_index_file=self.dbSNP_index_file,threads=self.threads)   
                                   
            if ret:
                logging.gemBS.gt("Methylation call done, samples performed: %s" %(ret))
                
                
    def extra_log(self):
        """Extra Parameters to be printed"""
        #Virtual methos, to be define in child class
        printer = logging.gemBS.gt
        
        printer("----------- Methylation Calling --------")
        printer("Reference       : %s", self.fasta_reference)
        printer("Species         : %s", self.species)
        printer("Right Trim      : %i", self.right_trim)
        printer("Left Trim       : %i", self.left_trim)
        printer("Chromosomes     : %s", self.list_chroms)
        printer("json File       : %s", self.json_file)
        printer("Threads         : %s", self.threads)
        if self.dbSNP_index_file != "":
            printer("dbSNP File      : %s", self.dbSNP_index_file)
        for sample,input_bam in self.sampleBam.iteritems():
            printer("Sample: %s    Bam: %s" %(sample,input_bam))
        printer("")
                                  
class MethylationFiltering(BasicPipeline):
    title = "Filtering of the output generated by the Methylation Calling."
    description = """ Filters all sites called as homozygous CC or GG with a 
                      probability of genotyping error <= 0.01
                      
                      Subset of dinucleotides called as CC/GG
                  """
                  
    def register(self,parser):
        ## required parameters
        parser.add_argument('-b','--bcf',dest='bcf_file',metavar="PATH",help="bcf Methylation call file", required=True)
        parser.add_argument('-o','--output-dir',dest="output_dir",metavar="PATH",help='Output directory to store the results.',default=None)
        
    def run(self,args):
        self.output_dir = args.output_dir
        self.bcf_file  = args.bcf_file       
        
        #Check bcf file existance
        if not os.path.isfile(args.bcf_file):
            raise CommandException("Sorry path %s was not found!!" %(args.bcf_file))
            
        #Call methylation filtering
        self.log_parameter()
        logging.gemBS.gt("Methylation Filtering...")
        ret = src.methylationFiltering(bcfFile=self.bcf_file,output_dir=self.output_dir)
        if ret:
            logging.gemBS.gt("Methylation filtering done, results located at: %s" %(ret))
            
            
    def extra_log(self):
        """Extra Parameters to be printed"""
        #Virtual methos, to be define in child class
        printer = logging.gemBS.gt
        
        printer("----------- Methylation Filtering ---------")
        printer("bcfFile         : %s", self.bcf_file)       
        printer("")
        
class BsCall(BasicPipeline):
    title = "Bisulfite calling for sample and chromosome."
    description = """ Tool useful for a cluster application manager. Methylation
                      calls for a given sample and chromosome.
                  """
     
    def membersInitiation(self):
        self.species = "HomoSapiens"
             
    def register(self,parser):
        ## required parameters
        parser.add_argument('-r','--fasta-reference',dest="fasta_reference",metavar="PATH",help="Path to the fasta reference file.",required=True)
        parser.add_argument('-e','--species',dest="species",metavar="SPECIES",default="HomoSapiens",help="Sample species name. Default: %s" %self.species)
        parser.add_argument('-s','--sample-id',dest="sample_id",metavar="SAMPLE",help="Sample unique identificator")  
        parser.add_argument('-c','--chrom',dest="chrom",metavar="CHROMOSOME",help="Chromosome name where is going to perform the methylation call")  
        parser.add_argument('-i','--input-bam',dest="input_bam",metavar="INPUT_BAM",help='Input BAM aligned file.',default=None)
        parser.add_argument('-g','--right-trim', dest="right_trim", metavar="BASES",type=int, default=0, help='Bases to trim from right of read pair, Default: 0')
        parser.add_argument('-f','--left-trim', dest="left_trim", metavar="BASES", type=int, default=5, help='Bases to trim from left of read pair, Default: 5')
        parser.add_argument('-o','--output-dir',dest="output_dir",metavar="PATH",help='Output directory to store the results.',default=None)
        parser.add_argument('-p','--paired-end', dest="paired_end", action="store_true", default=False, help="Input data is Paired End") 
        parser.add_argument('-t','--threads', dest="threads", metavar="THREADS", default="1", help='Number of threads, Default: %s' %self.threads)     
        parser.add_argument('-k','--keep-unmatched', dest="keep_unmatched", action="store_true", default=False, help="Do not discard reads that do not form proper pairs.")
        parser.add_argument('-u','--keep-duplicates', dest="keep_duplicates", action="store_true", default=False, help="Do not merge duplicate reads.")
        parser.add_argument('-d','--dbSNP-index-file', dest="dbSNP_index_file", metavar="FILE", help="dbSNP index file.",required=False,default="")

    def run(self,args):
        self.threads = args.threads
        self.reference = args.fasta_reference
        self.species = args.species
        self.input = args.input_bam 
        self.right_trim = args.right_trim
        self.left_trim = args.left_trim        
        self.chrom = args.chrom
        self.sample_id = args.sample_id
        self.output_dir = args.output_dir
        self.paired = args.paired_end
        self.keep_unmatched = args.keep_unmatched
        self.keep_duplicates = args.keep_duplicates
        self.dbSNP_index_file = args.dbSNP_index_file
        
        #Check fasta existance
        if not os.path.isfile(args.fasta_reference):
            raise CommandException("Sorry path %s was not found!!" %(args.fasta_reference))
        
        #Check input bam existance 
        if not os.path.isfile(args.input_bam):
            raise CommandException("Sorry path %s was not found!!" %(args.input_bam))
                        
        #Bs Calling per chromosome
        self.log_parameter()
        logging.gemBS.gt("BsCall per sample and chromosome...")
        
        ret = src.bsCalling (reference=self.reference,species=self.species,input_bam=self.input,chrom=self.chrom,
                             right_trim=self.right_trim, left_trim=self.left_trim,
                             sample_id=self.sample_id,output_dir=self.output_dir,
                             paired_end=self.paired,keep_unmatched=self.keep_unmatched,
                             keep_duplicates=self.keep_duplicates,dbSNP_index_file=self.dbSNP_index_file,threads=self.threads)
        if ret:
            logging.gemBS.gt("Bisulfite calling done: %s" %(ret)) 
       
    def extra_log(self):
        """Extra Parameters to be printed"""
        #Virtual methos, to be define in child class
        printer = logging.gemBS.gt
        
        printer("-------------- BS Call ------------")
        printer("Reference       : %s", self.reference)
        printer("Species         : %s", self.species) 
        printer("Chromosomes     : %s", self.chrom)
        printer("Sample ID       : %s", self.sample_id)
        printer("Threads         : %s", self.threads)
        printer("Right Trim      : %i", self.right_trim)
        printer("Left Trim       : %i", self.left_trim)
                
        if self.dbSNP_index_file != "":
            printer("dbSNP File      : %s", self.dbSNP_index_file)
        printer("")       
            
                  
class BsCallConcatenate(BasicPipeline):
    title = "Concatenation of methylation calls for different chromosomes for a given sample."  
    description = """ Concatenates bcf files comming from different methylation calls of 
                      different chromosomes.
                  """
    
    def register(self,parser):
        ## required parameters
        parser.add_argument('-s','--sample-id',dest="sample_id",metavar="SAMPLE",help="Sample unique identificator",required=True)
        parser.add_argument('-l','--list-bcfs',nargs="+",dest="list_bcfs",metavar="BCFLIST",help="List of bcfs to be concatenated.",required=True)
        parser.add_argument('-o','--output-dir',dest="output_dir",metavar="PATH",help='Output directory to store the results.',default=None)
        
    
    def run(self,args):
        self.list_bcf = args.list_bcfs  
        self.sample_id = args.sample_id
        self.output_dir = args.output_dir
        
        #Check bcf files to concatenate
        if len(args.list_bcfs) < 1:
            raise CommandException("No bcf files to concatenate.")
            
        for bcfFile in args.list_bcfs:
            #Check bcf existance 
            if not os.path.isfile(bcfFile):
                raise CommandException("Sorry path %s was not found!!" %(bcfFile))
                
        #Bs Calling Concatenate
        self.log_parameter()
        logging.gemBS.gt("BCF concatenate files...")
        args.list_bcfs.sort(key=lambda x: '{0:0>8}'.format(x).lower())        
        ret = src.bsConcat(list_bcfs=self.list_bcf,sample=self.sample_id,output_dir=self.output_dir)
        if ret:
            logging.gemBS.gt("BCF Concatenation Done: %s" %(ret))
            
    def extra_log(self):
        """Extra Parameters to be printed"""
        #Virtual methos, to be define in child class
        printer = logging.gemBS.gt
        
        printer("------- BS Call Concatenate ----------")
        printer("List BCF        : %s", self.list_bcf)
        printer("Sample ID       : %s", self.sample_id)
        printer("")       
            
      

      
class MappingReports(BasicPipeline):
    title = "Bisulfite Mapping reports. Builds a HTML and SPHINX report per lane and Sample."
    description = """ From json files lane stats, builds a HTML and SPHINX report per lane and sample """
    
    def register(self,parser):
        ## Mapping report stats parameters
        parser.add_argument('-j','--json',dest="json_file",metavar="JSON_FILE",help='JSON file configuration.',required=True)
        parser.add_argument('-i', '--input-dir', dest="input_dir",metavar="PATH", help='Path were are located the JSON stat files.', required=True)
        parser.add_argument('-n', '--name', dest="name", metavar="NAME", help='Output basic name',required=True)
        parser.add_argument('-o', '--output-dir', dest="output_dir", metavar="PATH",help='Output directory to store html report.',required=True)
         
         
    def run(self, args):
        self.name = args.name
        self.output_dir = args.output_dir
        
        #Recover json files from input-dir according to json file
        self.sample_lane_files = {}   
      
        self.records = 0
        for k,v in FLIdata(args.json_file).sampleData.iteritems():
            self.records = self.records + 1
            fileJson = "%s/%s.json" %(args.input_dir,v.getFli())
            if os.path.isfile(fileJson):
                if v.sample_barcode not in self.sample_lane_files: 
                   newFli = {}
                   newFli[v.getFli()] = [fileJson]
                   self.sample_lane_files[v.sample_barcode] = newFli
                elif v.getFli() not in self.sample_lane_files[v.sample_barcode]:
                   newFli = {}
                   newFli[v.getFli()] = [fileJson]
                   self.sample_lane_files[v.sample_barcode].update(newFli)
                elif v.getFli() in self.sample_lane_files[v.sample_barcode]:
                    self.sample_lane_files[v.sample_barcode][v.getFli()].append(fileJson)     
               
        #Check list of file
        if len(self.sample_lane_files) < 1:
            raise CommandException("Sorry not json files were found!!")

        self.log_parameter()
        logging.gemBS.gt("Building html reports...")
        report.buildReport(inputs=self.sample_lane_files,output_dir=self.output_dir,name=self.name)
        logging.gemBS.gt("Building sphinx reports...")
        sphinx.buildReport(inputs=self.sample_lane_files,output_dir="%s/SPHINX/" %(self.output_dir),name=self.name)
        logging.gemBS.gt("Report Done.")
         
    def extra_log(self):
        """Extra Parameters to be printed"""
        #Virtual methos, to be define in child class
        printer = logging.gemBS.gt
        
        printer("------- Mapping Report ----------")
        printer("Name            : %s", self.name)
        printer("")             
        
class VariantsReports(BasicPipeline):
    title = "BS Calls reports. Builds a HTML and SPHINX report per Sample."
    description = """ From chromosome stats json files, builds a HTML and SPHINX report per Sample """

    def membersInitiation(self):
        self.chroms = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

    def register(self,parser):
        ## variants reports stats parameters
        parser.add_argument('-j','--json',dest="json_file",metavar="JSON_FILE",help='JSON file configuration.',required=True)
        parser.add_argument('-i', '--input-dir', dest="input_dir",metavar="PATH", help='Path were are located the JSON variants stats files.', required=True)
        parser.add_argument('-n', '--name', dest="name", metavar="NAME", help='Output basic name',required=True)
        parser.add_argument('-o', '--output-dir', dest="output_dir", metavar="PATH",help='Output directory to store html and Sphinx Variants report.',required=True)
        
        parser.add_argument('-l','--list-chroms',dest="list_chroms",nargs="+",metavar="CHROMS",help="""List of chromosomes to perform the methylation pipeline.
                                                                                                       Can be a file where every line is a chromosome contig. 
                                                                                                       By default human chromosomes: %s """ %self.chroms,
                            default=["chr1","chr2","chr3","chr4","chr5","chr6","chr7",
                                     "chr8","chr9","chr10","chr11","chr12","chr13",
                                     "chr14","chr15","chr16","chr17","chr18","chr19",
                                     "chr20","chr21","chr22","chrX","chrY"])
                                     

    def run(self, args):
        self.name = args.name
        self.output_dir = args.output_dir
        self.json_file = args.json_file
        self.list_chroms = []
       
        if len(args.list_chroms) > 1:
            self.list_chroms = args.list_chroms
        elif os.path.isfile(args.list_chroms[0]):
            #Check if List_chroms is a file or just a list of chromosomes
            #Parse file to extract chromosme list 
            with open(args.list_chroms[0] , 'r') as chromFile:
                for line in chromFile:
                    self.list_chroms.append(line.rstrip())
        else:
            self.list_chroms = args.list_chroms            
       
        #Recover json files from input-dir according to json file
        self.sample_chr_files = {}
        self.sample_list = {}
        
        for k,v in FLIdata(args.json_file).sampleData.iteritems():
            self.sample_list[v.sample_barcode] = 0

        for sample,num in self.sample_list.iteritems():
            for chrom in self.list_chroms:
                fileJson = "%s/%s_%s.json" %(args.input_dir,sample,chrom)
                if os.path.isfile(fileJson):
                    if sample not in self.sample_chr_files:
                        self.sample_chr_files[sample] = [fileJson]
                    else:
                        self.sample_chr_files[sample].append(fileJson)
            
                
        self.log_parameter()
        logging.gemBS.gt("Building Bs Calls html and sphinx reports...")
        bsCallReports.buildBscallReports(inputs=self.sample_chr_files,output_dir=self.output_dir,name=self.name)        
        logging.gemBS.gt("Report Done.")                         

    def extra_log(self):
        """Extra Parameters to be printed"""
        #Virtual methos, to be define in child class
        printer = logging.gemBS.gt
        
        printer("------- Variants Reports ----------")
        printer("Name            : %s", self.name)
        printer("Json            : %s", self.json_file)
        printer("")   
        
                
class CpgBigwig(BasicPipeline):
    title = "Build BigWig files."
    description = """ Creates BigWig files to show pipeline results in Genome Browsers.
                      
                      Input Files:
                          CpG: Compressed CpG File to be transformed to Methylation and Coverage BigWig files.
                          Chromosome Sizes: File of chromosome lengths
                    
                      Creates methylation and coverage BigWig files.
                  """
                
    def register(self, parser):
        ## required parameters
        parser.add_argument('-c','--cpg-file', dest="cpg_file", help="""CpG gzipped Compressed File.""",required=True,default=None)
        parser.add_argument('-l','--chrom-length', dest="chrom_length", help="""Chromosome Length Text File.
                                                                                 Format: Two Columns: <chromosome name> <size in bases>""",required=True,default=None)
        parser.add_argument('-n','--name', dest="name", metavar="NAME", help='Output basic name',required=True)
        parser.add_argument('-q', '--quality', dest="quality", metavar="QUAL", help='Quality filtering criteria for the CpGs. By default 20.',required=False,default="20")
        parser.add_argument('-i', '--informative-reads', dest="informative_reads", metavar="READS", help='Total number of informative reads to filter CpGs.By default 5.',required=False,default="5")   
        parser.add_argument('-o','--output-dir',dest="output_dir",metavar="PATH",help='Output directory to store the results.',required=True,default=None)
                                                                                 

    def run(self,args):        
        self.name = args.name
        self.output_dir = args.output_dir
        self.cpg_file = args.cpg_file
        self.chrom_length = args.chrom_length
        self.quality = args.quality
        self.informative_reads = args.informative_reads
                
        #Check CpG gzipped compressed file
        if not os.path.isfile(args.cpg_file):
            raise CommandException("Sorry path %s was not found!!" %(args.cpg_file))    
            
        #Check chromosome length file exitance
        if not os.path.isfile(args.chrom_length):
            raise CommandException("Sorry path %s was not found!!" %(args.chrom_length)) 
        
      
                
        #Bs Calling Concatenate
        self.log_parameter()
        logging.gemBS.gt("CpG BigWig Conversion...")
                
        ret = src.cpgBigWigConversion(name=self.name,output_dir=self.output_dir,cpg_file=self.cpg_file,
                                      chr_len=self.chrom_length,quality=self.quality,informative_reads=self.informative_reads)
        if ret:
            logging.gemBS.gt("CpG Bigwig Conversion Done: %s" %(ret))
            
    def extra_log(self):
        """Extra Parameters to be printed"""
        #Virtual methos, to be define in child class
        printer = logging.gemBS.gt
        
        printer("------- CpG BigWig Conversion ----------")
        printer("Name            : %s", self.name)
        printer("CpG File        : %s", self.cpg_file)
        printer("Chrom Length    : %s", self.chrom_length)
        printer("Quality         : %s", self.quality)
        printer("Info. Reads     : %s", self.informative_reads)
        printer("")       
            
       

