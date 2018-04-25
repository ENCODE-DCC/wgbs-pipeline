# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 12:11:17 2017

@author: marcos
"""

from sphinx import ConfigSphinx

class SphinxBsCallReport(object):
    """ Basic Methods for BsCall Reports """
    
    def __init__(self,sphinx_file_name=None,currentName=""):
        """ Initializates Sphinx BsCall Basic Memebers """
        self.contents = [] #Vector of Sphinx Tags
        self.sphinx_file_name = sphinx_file_name
        self.currentName = currentName
        
    def addLevelSection(self,ident=20,title=None,char="="):
        """ Write Level Section
        
            ident - Characters to ident
            title - Section title
            char - Type of character
        """ 
        #1.Title Line
        titleLine = ""       
        for i in range(ident):
            titleLine += " "
        titleLine += title
        self.contents.append("%s\n" %(titleLine))
        
        #2.Section Line
        sectionLine = ""
        for i in range(ident):
            sectionLine += " "
        for i in range(len(title)):
            sectionLine += char
        
        #3. Append Section Line
        self.contents.append("%s\n" %(sectionLine))
        
    def addTopSection(self,ident=20,title=None):
        """ Write Top Section
        
            ident - Characters to ident
            title - Section title
            char - Type of character
        """ 
        self.addLevelSection(ident=ident,title=title,char="=")
    
    def addSubSection(self,ident=20,title=None):
        """ Write Sub Section
        
            ident - Characters to ident
            title - Section title
            char - Type of character
        """ 
        self.addLevelSection(ident=ident,title=title,char="-")
        
    def addSubSubSection(self,ident=20,title=None):
        """ Write Sub sub section
        
            ident - Characters to ident
            title - Section title
            char - Type of character
        """
        self.addLevelSection(ident=ident,title=title,char="^")
        
    def addLevelHighlight(self,ident=20,title=None,char="*"):
        """ Write Level Section
        
            ident - Characters to ident
            title - Section title
            char - Type of character
        """
        #1. Over Line
        overLine = ""
        for i in range(ident):
            overLine += " "          
        for i in range(len(title)):
            overLine += char
        self.contents.append("%s\n" %(overLine))          
        #2.Title Line
        titleLine = ""              
        for i in range(ident):
            titleLine += " "
        titleLine += title
        self.contents.append("%s\n" %(titleLine))    
        #3.Section Line
        sectionLine = ""    
        for i in range(ident):
            sectionLine += " "
        for i in range(len(title)):
            sectionLine += char       
        #3. Append Section Line
        self.contents.append("%s\n" %(sectionLine))
        
    def addPartLine(self,ident=20,title=None):
        """ Add Part Line (Main Highline)
        
            ident - Characters to ident
            title - Section title
            char - Type of character
        """
        self.addLevelHighlight(ident=ident,title=title,char="#")
        
    def addChapterLine(self,ident=20,title=None):
        """ Add Chapter Line (Chapter Highline)
        
            ident - Characters to ident
            title - Section title
            char - Type of character
        """
        self.addLevelHighlight(ident=ident,title=title,char="*")
               
    def addTocTree(self,vectorDocuments=None):
        """ Add Toc Tree Documents

            vectorDocuments - List of documents to write at the index
        """
        content =   ".. toctree::\n"
        content +=  "   :maxdepth: 3\n"
        content +=  "\n"
        for document in vectorDocuments:
            content += "   %s\n" %(document)
        self.contents.append(content)        
                
    def addLine(self,ident=20,cells=0,lenCell=30,char="-"):
        """ Write Line
        
            ident - Characters to ident
            cells - Number of cell in the rows
            char - Type of character
        """       
        underLine = []
        for i in range(ident):
            underLine.append(" ")   
            
        for cell in range(cells):
            underLine.append("+")
            for chars in range(lenCell-1):
                underLine.append(char)
            if cell == (cells-1):
                underLine[-1] = "+"
            
        underLine.append("\n")
        self.contents.append(''.join(underLine))    
     
    
    def addSimpleLine(self,ident=20,cells=0,lenCell=30):
        """ Underline
        
            ident - Characters to ident
            cells - Number of cell in the rows
            isLeft - length of each cell
        """       
        self.addLine(ident=ident,cells=cells,lenCell=lenCell,char="-")
        
    def addHeaderLine(self,ident=20,cells=0,lenCell=30):
        """ Underline
        
            ident - Characters to ident
            cells - Number of cell in the rows
            isLeft - length of each cell
        """       
        self.addLine(ident=ident,cells=cells,lenCell=lenCell,char="=")

    def addCell(self,text="",lenCell=30,isLeft=True,isRight=False):
        """ Generate Cell
        
            ident - Characters to ident
            text - Cell Contents
            lenCell - Cell length
            isLeft - True if is the first cell
            isRight - True if is the last cell
            return cell text
        """
        cellText = []
        #Add Contents
        cellText.append("|")
        cellText.append(text)
        #Add with spaces
        for i in range(lenCell - len(''.join(cellText))):
            cellText.append(" ")
        if isRight:
            cellText[-1] = "|"
            cellText.append("\n")
        
        return ''.join(cellText)
        
        
    def addRow(self,ident=20,lenCell=30,values=None): 
        """ Add row sphinx
            
            ident - Characters to ident
            lenCell - Length of the Cell in each row
            values - Vector of values to be printed in a row
        """
        cells = 0
        cellContents = ""
        #Add ident chars
        for i in range(ident):
            cellContents += " " 
 
        lenFields = len(values)
 
        for field in values:
            isLeft = False
            isRight = False
            if cells == 0:
                isLeft = True
            elif cells == (lenFields-1):
                isRight = True            
            
            cellContents += self.addCell(text="%s" %(field),lenCell=lenCell,isLeft=isLeft,isRight=isRight)
            cells += 1
            
        self.contents.append(cellContents) 
            
        #Underline
        self.addSimpleLine(ident=ident,cells=cells,lenCell=lenCell)
        
    def createTable(self,ident=20,lenCell=30,values=[]):
        """ Create General Variants Table

           ident - Characters to ident 
           lenCell - Cell length
           values - Array of values
        """
        cellsHeader = 0
        
        #1.Table header
        cellsHeader = len(values[0])
        #2. Write overline
        self.addSimpleLine(ident=ident,cells=cellsHeader,lenCell=lenCell)
            
        #3.Stats table
        cellContents = ""
        
         #Add ident chars
        for i in range(ident):
            cellContents += " "  
        
        i = 0
        total = len(values[0])
                
        for conceptHeader in values[0]:
            isLeft = False
            isRight = False
            if i == 0:
                isLeft = True
            elif i == total-1:
                isRight = True
                    
            cellContents += self.addCell(text=conceptHeader,lenCell=lenCell,isLeft=isLeft,isRight=isRight)
            i = i + 1
        
        self.contents.append(cellContents)
            
        #Header Line 
        self.addHeaderLine(ident=ident,cells=cellsHeader,lenCell=lenCell)
                
        #2. Table Contents
        for i in range(1,len(values)):        
            self.addRow(ident=ident,lenCell=lenCell,values=values[i])    
            
            
    def buildImage(self,pathImage=None):
        """ Build image
  
            pathImage - PNG file to be shown at sphinx document      
        """
        #1. SPHINX TAG BUILDING 
        imageLine = "    .. image:: %s\n" %(pathImage)
        self.contents.append(imageLine)


            
    def saveReport(self):
        '''Save the report to a given name file'''
        with open(self.sphinx_file_name, 'w') as fileDocument:
            for line in self.contents:
                fileDocument.write(line)

    def createPage(self):
        '''To Be Defined in the Child Class'''
        pass
    
    def configureStats(self,stats_vector=None):
        '''To Be Defined in the Child Class'''
        pass
    
    
class SphinxMappingCoverage(SphinxBsCallReport):
    """Class to manage Mapping Coverage Report for Sphinx Format"""
    
    
    def __init__(self,sphinx_file_name=None,current_name=None):
        """  Class constructor
        
             sphinx_file_name -- SPHINX file document to be build
             current_name -- Current main name
        """
        #Call Parent class constructor
        SphinxBsCallReport.__init__(self,sphinx_file_name=sphinx_file_name,currentName=current_name)
        
        
    def configureStats(self,stats_vector=None):
        """ Configure Stats Methods
        
            stats_vector - Vector of BsCallStats to build HtmlMappingCoverage  [readLevelStats,baseLevelStats,allCoverage,gcCoverage,qualityAll,nonCpGReadProfile]
        """
        self.read_level_stats = stats_vector[0]
        self.base_level_stats = stats_vector[1]
        self.all_coverage = stats_vector[2]
        self.gc_coverage = stats_vector[3]
        self.quality_all = stats_vector[4]
        self.non_cpg_read_profile = stats_vector[5]
        
    def createPage(self):
        '''Define SPHINX  documents'''
        ident = 0
        lenCell = 30
        #0. Mapping and Coverage Statistics
        self.addTopSection(ident=0,title="Alignments and Coverage Statistics for %s" %(self.currentName))        
        self.contents.append("\n")         
        #1. Read Level Counts
        self.addSubSection(ident=ident,title="Read Level Counts")
        self.contents.append("\n")
        self.createTable(ident=ident,lenCell=lenCell,values=self.read_level_stats.getTable())        
        self.contents.append("\n")        
        #2. Base Level Counts
        self.addSubSection(ident=ident,title="Base Level Counts")
        self.contents.append("\n")
        self.createTable(ident=ident,lenCell=lenCell,values=self.base_level_stats.getTable())        
        self.contents.append("\n")        
        #3. Coverage All Bases and Quality All Bases
        self.addSubSection(ident=ident,title="Coverage Distribution")
        self.contents.append("\n")
        self.buildImage(pathImage=self.all_coverage.relativeSphinxPathImage)
        self.contents.append("\n")
        self.addSubSection(ident=ident,title="Quality Distribution")
        self.contents.append("\n")
        self.buildImage(pathImage=self.quality_all.relativeSphinxPathImage)
        self.contents.append("\n")
        #4. GC Coverage and Non CpG Read Profile
        self.addSubSection(ident=ident,title="GC/Coverage Heatmap")
        self.contents.append("\n")
        self.buildImage(pathImage=self.gc_coverage.relativeSphinxPathImage)
        self.contents.append("\n")
        self.addSubSection(ident=ident,title="Non-CpG Read Profile")
        self.contents.append("\n")
        self.buildImage(pathImage=self.non_cpg_read_profile.relativeSphinxPathImage)
        self.contents.append("\n")
        #5. Save HTML Document
        self.saveReport()
 

class SphinxBsGenotypeCalls(SphinxBsCallReport):
    """Class to manage Bisulfite Genotype Calling Report for Sphinx Format"""    
    
    def __init__(self,sphinx_file_name=None,current_name=None):
        """  Class constructor
        
             sphinx_file_name -- SPHINX file document to be build
             current_name -- Current main name
        """
        #Call Parent class constructor
        SphinxBsCallReport.__init__(self,sphinx_file_name=sphinx_file_name,currentName=current_name)
           
    def configureStats(self,stats_vector=None):
        """ Configure Stats Methods
        
            stats_vector - Vector of BsCallStats to build HtmlBsGenotypeCalls 
                           [totalStats,vcfFilterStats,variantCoverage,dbSnpCoverage,qualityVariant,
                            [fsVariant,qdVariant,qdNonVariant,rmsmqVariant,rmsmqNonVariant,gofVariant,gofNonVariant],
                            mutations]
        """
        self.total_stats = stats_vector[0]
        self.vcf_filter_stats = stats_vector[1]
        self.variant_coverage = stats_vector[2]
        self.dbSNP_coverage = stats_vector[3]
        self.quality_variant = stats_vector[4]

        self.fsVariant = stats_vector[5][0]
        self.qdVariant = stats_vector[5][1]
        self.qdNonVariant = stats_vector[5][2]
        self.rmsmqVariant = stats_vector[5][3]
        self.rmsmqNonVariant = stats_vector[5][4]
        self.gofVariant = stats_vector[5][5]
        self.gofNonVariant = stats_vector[5][6]        
                
        self.mutations = stats_vector[6]
        
    def createPage(self):
        '''Define SPHINX documents'''
        ident = 0
        lenCell = 30
        #0. Bisulfite Calling Statistics
        self.addTopSection(ident=0,title="Bisulfite Calling for %s" %(self.currentName))        
        self.contents.append("\n")                 
        #1. Variant Counts
        self.addSubSection(ident=ident,title="Variant Counts")
        self.contents.append("\n")
        self.createTable(ident=ident,lenCell=lenCell,values=self.total_stats.getTableVariants())        
        self.contents.append("\n")             
        #2.VCF Filtering Stats
        self.addSubSection(ident=ident,title="VCF Filtering Stats")
        self.contents.append("\n")
        self.createTable(ident=ident,lenCell=lenCell,values=self.vcf_filter_stats.getTable())        
        self.contents.append("\n")  
        #3. Coverage SNP Bases and Quality Variants
        self.addSubSection(ident=ident,title="Coverage Variants")
        self.contents.append("\n")
        self.buildImage(pathImage=self.variant_coverage.relativeSphinxPathImage)
        self.contents.append("\n")
        self.addSubSection(ident=ident,title="Quality Variants")
        self.contents.append("\n")
        self.buildImage(pathImage=self.quality_variant.relativeSphinxPathImage)
        self.contents.append("\n")        
        #4. Coverage DbSNPS
        if self.dbSNP_coverage.getHaveCounts():
            self.addSubSection(ident=ident,title="Coverage DbSNPs Variant Sites")
            self.contents.append("\n")
            self.buildImage(pathImage=self.dbSNP_coverage.relativeSphinxPathImage)
            self.contents.append("\n")        
        #5.QC Distribution
        self.addSubSection(ident=ident,title="Filtering Criteria Distribution")
        self.contents.append("\n")
        #5.1 FisherStrand Variant
        self.addSubSubSection(ident=ident,title="Strand bias estimated using Fisher's Exact Test.")
        self.contents.append("\n")
        self.buildImage(pathImage=self.fsVariant.relativeSphinxPathImage)
        self.contents.append("\n")
        #5.2 QualityByDepth Variant NonVariant
        self.addSubSubSection(ident=ident,title="Allele-specific call confidence normalized by depth of sample reads supporting the allele. Variants.")
        self.contents.append("\n")
        self.buildImage(pathImage=self.qdVariant.relativeSphinxPathImage)
        self.contents.append("\n")
        self.addSubSubSection(ident=ident,title="Allele-specific call confidence normalized by depth of sample reads supporting the allele. Non-Variants.")
        self.contents.append("\n")
        self.buildImage(pathImage=self.qdNonVariant.relativeSphinxPathImage)
        self.contents.append("\n")              
        #5.3 RMSMappingQuality Variant NonVariant
        self.addSubSubSection(ident=ident,title="Allele-specific Root Mean Square of the mapping quality of reads across all samples. Variants.")
        self.contents.append("\n")
        self.buildImage(pathImage=self.rmsmqVariant.relativeSphinxPathImage)
        self.contents.append("\n") 
        self.addSubSubSection(ident=ident,title="Allele-specific Root Mean Square of the mapping quality of reads across all samples. Non-Variants.")
        self.contents.append("\n")
        self.buildImage(pathImage=self.rmsmqNonVariant.relativeSphinxPathImage)
        self.contents.append("\n") 
        #5.4 GoodnessOfFit Variant NonVariant
        self.addSubSubSection(ident=ident,title="Phred scaled goodness of fit. Variants.")
        self.contents.append("\n")
        self.buildImage(pathImage=self.gofVariant.relativeSphinxPathImage)
        self.contents.append("\n") 
        self.addSubSubSection(ident=ident,title="Phred scaled goodness of fit. Non-Variants.")
        self.contents.append("\n")
        self.buildImage(pathImage=self.gofNonVariant.relativeSphinxPathImage)
        self.contents.append("\n")         
        #6. Mutation All 
        self.addSubSection(ident=ident,title="Mutation All")
        self.contents.append("\n")
        self.createTable(ident=ident,lenCell=lenCell,values=self.mutations.getTableMutationProfile())        
        self.contents.append("\n")  
        #7. TiTv
        self.addSubSection(ident=ident,title="Ti/Tv Ratio")
        self.contents.append("\n")
        self.createTable(ident=ident,lenCell=lenCell,values=self.mutations.getTiTvTable())        
        self.contents.append("\n")          
        #8. Save SPHINX Document
        self.saveReport()
             
        
class SphinxMethylation(SphinxBsCallReport):
    """Class to manage Methylation Report for Sphinx Format"""
 
    def __init__(self,sphinx_file_name=None,current_name=None):
        """  Class constructor
        
             sphinx_file_name -- SPHINX file document to be build
             current_name -- Current main name
        """
        #Call Parent class constructor
        SphinxBsCallReport.__init__(self,sphinx_file_name=sphinx_file_name,currentName=current_name)   
        
    def configureStats(self,stats_vector=None):
        """ Configure Stats Methods
        
            stats_vector - Vector of BsCallStats to build HtmlMethylation
                           [totalStats,refCpGcoverage,refCpGInfCoverage,nonRefCpGcoverage,nonRefCpGinfCoverage,
                            qualityRefCpG,qualityNonRefCpG,plotMethylation,summaryMethylation]                            
        """          
        self.total_stats = stats_vector[0]
        self.ref_cpg_coverage = stats_vector[1]
        self.ref_cpg_inf_coverage = stats_vector[2]
        self.non_ref_cpg_coverage = stats_vector[3]
        self.non_ref_cpg_inf_coverage = stats_vector[4]
        self.quality_ref_cpg = stats_vector[5]
        self.quality_non_ref_cpg = stats_vector[6]
        self.methylation_plot_cpg = stats_vector[7]
        self.summary_methylation = stats_vector[8]
        
    def createPage(self):
        '''Define SPHINX documents'''
        #0. Methylation Statistics
        ident = 0
        lenCell = 30
        self.addTopSection(ident=0,title="Methylation statistics for %s" %(self.currentName))        
        self.contents.append("\n")           
        #1. Total Stats Ref CpGs/NonRefCpG
        self.addSubSection(ident=ident,title="CpG Counts")
        self.contents.append("\n")
        self.createTable(ident=ident,lenCell=lenCell,values=self.total_stats.getTableMethylation())        
        self.contents.append("\n")
        #2. Coverage CpG and Coverage CpG Inf
        self.addSubSection(ident=ident,title="Coverage Reference CpGs")
        self.contents.append("\n")
        self.buildImage(pathImage=self.ref_cpg_coverage.relativeSphinxPathImage)
        self.contents.append("\n")
        self.addSubSection(ident=ident,title="Coverage Information Reads Reference CpGs")
        self.contents.append("\n")
        self.buildImage(pathImage=self.ref_cpg_inf_coverage.relativeSphinxPathImage)
        self.contents.append("\n")   
        #3. Non Coverage CpG and Non Coverage CpG Inf
        self.addSubSection(ident=ident,title="Coverage Non-Reference CpGs")
        self.contents.append("\n")
        self.buildImage(pathImage=self.non_ref_cpg_coverage.relativeSphinxPathImage)
        self.contents.append("\n")
        self.addSubSection(ident=ident,title="Coverage Information Reads Non-Reference CpGs")
        self.contents.append("\n")
        self.buildImage(pathImage=self.non_ref_cpg_inf_coverage.relativeSphinxPathImage)
        self.contents.append("\n")   
        #4. Quality RefCpG and Quality Non Reference CpG
        self.addSubSection(ident=ident,title="Quality Reference CpGs")
        self.contents.append("\n")
        self.buildImage(pathImage=self.quality_ref_cpg.relativeSphinxPathImage)
        self.contents.append("\n")
        self.addSubSection(ident=ident,title="Quality Non-Reference CpGs")
        self.contents.append("\n")
        self.buildImage(pathImage=self.quality_non_ref_cpg.relativeSphinxPathImage)
        self.contents.append("\n")
        #5. Methylation distribution at CpGs
        self.addSubSection(ident=ident,title="Methylation distribution at CpGs")
        self.contents.append("\n")
        self.buildImage(pathImage=self.methylation_plot_cpg.relativeSphinxPathImage)
        self.contents.append("\n")
        #6. Summary Methylation Profile
        self.addSubSection(ident=ident,title="CpGs Methylation Profile (From Individual Strands)")
        self.contents.append("\n")
        self.createTable(ident=ident,lenCell=lenCell,values=self.summary_methylation.getTable())        
        self.contents.append("\n")
        #8. Save SPHINX Document
        self.saveReport()        
        
class SphinxIndexBsCall(SphinxBsCallReport):
    """ Creates and Manages the Index SPHINX Document for BSCall Reports """
    
    def __init__(self,output_dir=None,name_project=None,samples_stats=None):
        """  Class constructor
        
             output_dir -- Output directory to store HTML reports
             name_project -- Project name
             samples_stats -- Dictionary of samples and list of stats objects. [sample][ListStatsObject]
        """
        #Call Parent class constructor
        SphinxBsCallReport.__init__(self,sphinx_file_name="%s/%s.rst" %(output_dir,name_project),currentName=name_project)
        
        self.output_dir = output_dir
        self.samples_stats = samples_stats
                
    def createPage(self):
        """ Run Index for BsCall""" 
        #1. Add SPHINX Report Header
        #Scaped name
        scaped_name = ""
        for character in self.currentName:
            if character == "_":
                scaped_name += "\\%s" %(character)
            else:
                scaped_name += "%s" %(character)        
        
        self.addPartLine(ident=0,title="BScall Reports Project %s"%(scaped_name))

        vectorDocuments = []
          
        #2. Create Sample and lanes Report 
        for sample in self.samples_stats:
            #2.1 Create Mapping and Coverage Statistics
            mappingCoverageSphinx = "%s/%s_mapping_coverage.rst" %(self.output_dir,sample) 
            sample_mapping_coverage = SphinxMappingCoverage(sphinx_file_name=mappingCoverageSphinx,current_name=sample)
            #2.1.2 Setup mapping coverage report
            sample_mapping_coverage.configureStats(stats_vector=self.samples_stats[sample]["mappingCoverage"])
            #2.1.3 Create Document
            sample_mapping_coverage.createPage()

            #2.2 Create Bs-Genotype Calls Report
            variantsSphinx = "%s/%s_variants.rst" %(self.output_dir,sample) 
            sample_variants = SphinxBsGenotypeCalls(sphinx_file_name=variantsSphinx,current_name=sample)
            #2.2.1 Setup Variants Report
            sample_variants.configureStats(stats_vector=self.samples_stats[sample]["calls"])
            #2.2.2 Create Document
            sample_variants.createPage()

            #2.3 Create Methylation Statitics
            methylationSphinx = "%s/%s_methylation.rst" %(self.output_dir,sample)
            sample_methylation = SphinxMethylation(sphinx_file_name=methylationSphinx,current_name=sample)
            #2.3.1 Setup Methylation Report
            sample_methylation.configureStats(stats_vector=self.samples_stats[sample]["methylation"])
            #2.3.2 Create Document
            sample_methylation.createPage()
            
            #2.4 Vector of links to create index Table
            vectorDocuments.extend(["%s_mapping_coverage"%(sample),"%s_variants"%(sample),"%s_methylation"%(sample)])
                     
        #3.Toc Table
        self.contents.append("\n")
        self.addTocTree(vectorDocuments=vectorDocuments)
        #3.1 Toc Table Save Report
        self.saveReport()
