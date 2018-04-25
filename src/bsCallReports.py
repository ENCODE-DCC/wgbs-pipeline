# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 12:02:40 2017

@author: marcos
"""

from reportStats import RunBasicStats
from report import BasicHtml
from bsCallStats import *
from bsCallSphinxReports import *

import os
import json


class HtmlBsCallReport(object):
    """ Basic Methods for BsCall Reports """
    
    def __init__(self,html_file_name=None,parentName="",currentName="",parentDocument=""):
        """ Initializates Html BsCall Basic Memebers """
        self.tags = [] #Vector of html tags
        self.html_file_name = html_file_name
        self.parentName = parentName
        self.currentName = currentName
        self.parentDocument = parentDocument
        
    def start(self):
        ''' Add HTML Report Start tags to a given vector'''
        self.tags.append("<HTML>\n")

        self.tags.append(" <HEAD>\n")
        self.tags.append(" <STYLE TYPE=\"text/css\">\n")
        self.tags.append("  <!--\n")
        self.tags.append("   @import url(\"style.css\"); \n")
        self.tags.append("  -->\n")
        self.tags.append(" </STYLE>\n")
        self.tags.append(" </HEAD>\n")
        self.tags.append(" <BODY>\n")
        
    def stop(self):
        ''' Close HTML header to a given vector object '''  
        self.tags.append(" </BODY>\n")
        self.tags.append("</HTML>\n")

    def space(self):
        ''' Add HTML Report Header to a given vector'''
        self.tags.append("<BR><BR><BR>\n")
        
    def header(self):
        ''' Add HTML Report Header to a given vector'''
        self.start()

        self.tags.append("\n <P id='path'> /%s/%s </P> \n" %(self.parentName,self.currentName))
        self.tags.append("\n <a class=\"link\" href=\"%s\"><B>BACK</B></a> <br>\n" %(os.path.basename(self.parentDocument)))   
        self.tags.append("  <H1 id=\"title\"> <U> SAMPLE %s </U> </H1>\n" %(self.currentName))
        
    def headerIndex(self,title=""):
        ''' Add HTML Report Header for Indexing Report '''
        self.start()
        self.tags.append('<H1 id="title"> <U> %s </U> </H1> \n' %(title))        
        
    def sectionTitle(self,title=""):
        ''' Adds Section Title '''
        self.tags.append('<H3 id="section"> %s </H3>\n' %(title))
              
    def addTableHeader(self,values=None):
        """ Add Row html tags to a table
            
            values - vector of row values to be printed
        """        
        self.tags.append("   <TR>\n")
        for concept in values:        
            self.tags.append("    <TH scope=\"col\">%s</TH>\n"%(concept))
            
        self.tags.append("   </TR>\n")        

    def emptyFields(self,values):
        """ Return True if oll fields in values are empty """
        totalText = ""
        for field in values:
            if field != "":
                totalText += "%s" %(field)
        if totalText == "":
            return True
        else:
            return False
     
    def addRow(self,values=None,isOdd=True):
        """ Add Row html tags to a table
            
            values - vector of row values to be printed
            isOdd - True/False for the proper row color
        """
        empty = ""        
        if self.emptyFields(values):
            empty += "empty"

        if isOdd == True: 
            self.tags.append("   <TR class=\"odd %s\">\n"%(empty))
        else:
            classEmpty = ""
            if empty != "":
                classEmpty = 'class="empty"'                
            self.tags.append("   <TR %s>\n" %(classEmpty))
        
        for field in values:        
            self.tags.append("   <TD> %s </TD> \n" %(field))    

        self.tags.append("   </TR>\n") 
        
    def tableDefinition(self,color=None):
        """ HTML Table tag definition
        
        color - Table color, could be green or blue
        """
        if color == "green":
            self.tags.append("  <TABLE id=\"green\">\n")
        else:
            self.tags.append("  <TABLE id=\"hor-zebra\">\n")
            
    def createTable(self,color=None,values=None,tableTitle=None):
        """Create Statistics Table
        
           color - Table color, could be green or blue  
           values - Vector of values to be printed, first row corrresponds to the header (Titles)
           tableTitle - Table Title
        """
        #0.Table Title
        self.sectionTitle(title=tableTitle)
     
        #1.Table Definition
        self.tableDefinition(color=color)

        #2.Table header
        self.addTableHeader(values[0])        
        
        #3.Fullfill table
        isOdd = True
        totalRows = len(values)
        for index in range(1,totalRows):    
            self.addRow(values=values[index],isOdd=isOdd)
            isOdd = not isOdd
        
        self.tags.append(" </TABLE>\n")
                
    def createTablePlots(self,color=None,plotsTitle=[],imagePaths=[]):
        """Build Plots and locate on html report
        
          color - Color for the table were the plot is going to be located
          plotsTitle - List of title of plots
          imagePaths - List of path to the images to be shown
        """
        self.tableDefinition(color=color)

        self.tags.append('   <TR class="odd">\n')
        for title in plotsTitle:
            self.tags.append("    <TH scope=\"col\"> %s </TH> \n" %(title))
        self.tags.append('   </TR>\n')
            
        self.tags.append('   <TR>\n')    
        for plotPath in imagePaths:       
            self.tags.append('   <TD> <img src="%s" alt="%s"> </TD> \n' %(plotPath,plotPath))
        self.tags.append('   </TR>\n')
        
        self.tags.append("  </TABLE>\n")
        
    def createPlot(self,stats=None,color="blue",title=""):
        """ Adds a table Document of Plots
        
            stats -- BsCall Stats Instances to get the plots
            titles -- Title plot
        """
        self.space()
        if stats.getHaveCounts():
            stats.plot()
            self.createTablePlots(color=color,plotsTitle=[title],imagePaths=[stats.relativePathImage])
            
    def createSetPlots(self,stats_objects=[],color="blue",titles=[]):
        """ Adds a table Document of Plots
        
            stats -- List of BsCall Stats instances to get the plots
            titles -- List Of titlesTitle plot
        """
        self.space()
        title_list = []
        image_paths = []
        for stats,title in zip(stats_objects,titles):
            if stats.getHaveCounts:
                stats.plot()
                image_paths.append(stats.relativePathImage)
                title_list.append(title)
        
        if len(image_paths) > 0 or len(title_list) > 0:
            self.createTablePlots(color=color,plotsTitle=title_list,imagePaths=image_paths)
               
    def createLinksTable(self,samples_links=[],summary_field=None,color=""):
        """ Create links Table
        
            samples_links - Vector of links to be generate per each sample. sampleName,linkHtmlMappingCoverage,linkHtmlGenotype,linkHtmlMethylation
            color - Color for the table
            summary_field - [sample]{"alignments":ListAligmentAndCoverage,"variants":ListVariantsAndCoverage,"methylation":ListMethylation}             
        """
        #1. Table Tags and color definition
        self.tableDefinition(color=color)
        #2. Define link color
        link_class = ''	
        if color == "green":
            link_class += '"linkGreen"'
        else:
            link_class += '"link"'
        #3. Headers
        header = ["Sample","Aligned","Uniquely Aligned","Used for Calling","GC Coverage Correlation"]
        header.extend(["Variants","Variants Passed","Mean Coverage Variants Passed","Ti/Tv Ratio"])
        header.extend(["Mean CpG Methylation","Mean CpG Coverage","Passed CpGs"])
        header.append("Reports")
        #3.1. HTML Header 
        self.tags.append('   <TR>\n')
        for row in header:
            self.tags.append('      <TH scope=\"col\"> %s </TH> \n'%(row))
        self.tags.append('   </TR>\n')
        #4. Fullfill table
        isOdd = True
        totalRows = len(samples_links)
        for index in range(1,totalRows):
            #4.0 Links to reports
            sample_name = samples_links[index][0]
            alignments_link = samples_links[index][1]
            variants_link = samples_links[index][2]
            methylation_link = samples_links[index][3]
            #4.1 Row Definition
            class_styles = "" 
            if isOdd:
                class_styles += 'class="odd"'
            #4.2.1 Row for Headers of Alignents and coverage   
            rowValues = [sample_name]
            rowValues.extend(summary_field[sample_name]['alignments'])
            rowValues.extend(summary_field[sample_name]['variants'])
            rowValues.extend(summary_field[sample_name]['methylation'])
            #4.2.2 Report links
            links_report = '<A class=%s href="%s"> &#187 Alignments & Coverage </A> <BR> \n'%(link_class,alignments_link)
            links_report += '<A class=%s href="%s"> &#187 Variants </A> <BR> \n'%(link_class,variants_link)
            links_report += '<A class=%s href="%s"> &#187 Methylation </A> <BR> \n'%(link_class,methylation_link)
            rowValues.append(links_report)
            #4.2.3 Contents 
            self.tags.append('    <TR %s>\n' %(class_styles))
            for row in rowValues:
                self.tags.append('      <TD>%s</TD>\n'%(row))
            #4.2.4 Close row            
            self.tags.append('   </TR>\n')
            #4.3 Row closing            
            isOdd = not isOdd
            
        #5.Close Table
        self.tags.append(" </TABLE>\n")

            
    def saveReport(self):
        '''Save the report to a given name file'''
        with open(self.html_file_name, 'w') as fileDocument:
            for line in self.tags:
                fileDocument.write(line)

    def createPage(self):
        '''To Be Defined in the Child Class'''
        pass
    
    def configureStats(self,stats_vector=None):
        '''To Be Defined in the Child Class'''
        pass
    
class HtmlMappingCoverage(HtmlBsCallReport):
    """ Class which Mapping Coverage """
    
    def __init__(self,html_file_name=None,current_name=None,parent_name=None,parent_document=None):
        """  Class constructor
        
             html_file_name -- HTML file document to be build
             current_name -- Current main name
             parent_name -- Parent name of the parent document
             parent_document -- Parent document
        """
        #Call Parent class constructor
        HtmlBsCallReport.__init__(self,html_file_name=html_file_name,parentName=parent_name,currentName=current_name,parentDocument=parent_document)

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
        '''Define HTML web documents'''
        #0. Mapping and Coverage Statistics
        self.header()
        #1. Read Level Counts
        self.space()
        self.createTable(color="blue",values=self.read_level_stats.getTable(),tableTitle="Read Level Counts")
        #2. Base Level Counts
        self.space()
        self.createTable(color="green",values=self.base_level_stats.getTable(),tableTitle="Base Level Counts")
        #3. Coverage All Bases and Quality All Bases
        self.sectionTitle(title="Coverage and Quality")
        self.createSetPlots(stats_objects=[self.all_coverage,self.quality_all],color="blue",titles=["Coverage Distribution","Quality Distribution"])
        #4. GC Coverage and Non CpG Read Profile
        self.createSetPlots(stats_objects=[self.gc_coverage,self.non_cpg_read_profile],color="green",titles=["GC/Coverage Heatmap","Non-CpG Read Profile"])
        #5. Save HTML Document
        self.stop()
        self.saveReport()
        
        
class HtmlBsGenotypeCalls(HtmlBsCallReport):
    """ Class which BS-Genotype Calls """
    
    def __init__(self,html_file_name=None,current_name=None,parent_name=None,parent_document=None):
        """  Class constructor
        
             html_file_name -- HTML file document to be build
             current_name -- Current main name
             parent_name -- Parent name of the parent document
             parent_document -- Parent document
        """
        #Call Parent class constructor
        HtmlBsCallReport.__init__(self,html_file_name=html_file_name,parentName=parent_name,currentName=current_name,parentDocument=parent_document)

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
        
        #Set GoF Axis X Equal for both
        number_locations_gof = self.gofVariant.selectTotalNumberLocations(self.gofNonVariant)
        self.gofVariant.setNumberOfLocations(locationsToRecover=number_locations_gof)
        self.gofNonVariant.setNumberOfLocations(locationsToRecover=number_locations_gof)
                
        self.mutations = stats_vector[6]
        
    def createPage(self):
        '''Define HTML web documents'''
        #0. Bisulfite Calling Statistics
        self.header()
        #1. Total Stats
        self.space()
        self.createTable(color="blue",values=self.total_stats.getTableVariants(),tableTitle="Variant counts")
        #3.VCF Filtering Stats
        self.space()
        self.createTable(color="green",values=self.vcf_filter_stats.getTable(),tableTitle="VCF Filtering Stats")
        #3. Coverage SNP Bases and Quality Variants        
        self.sectionTitle(title="Coverage and Quality")
        self.createSetPlots(stats_objects=[self.variant_coverage,self.quality_variant],color="green",titles=["Coverage Variants","Quality Variants"])        
        
        #4. Coverage DbSNPS
        self.createPlot(stats=self.dbSNP_coverage,color="blue",title="Coverage DbSNPs Variant Sites")
        
        #5.QC Distribution
        self.sectionTitle(title="Filtering Criteria Distribution")
        #5.1 FisherStrand Variant
        self.createPlot(stats=self.fsVariant,color="green",title="Phred scale strand bias estimated using Fisher's Exact Test.")
        #5.2 QualityByDepth Variant NonVariant
        self.createSetPlots(stats_objects=[self.qdVariant,self.qdNonVariant],color="blue",
                            titles=["Allele-specific call confidence normalized by depth of sample reads supporting the allele. Variants.",
                                    "Allele-specific call confidence normalized by depth of sample reads supporting the allele. Non-Variants."])        
        #5.3 RMSMappingQuality Variant NonVariant
        self.createSetPlots(stats_objects=[self.rmsmqVariant,self.rmsmqNonVariant],color="green",
                            titles=["Root Mean Square of the mapping quality of reads. Variants.",
                                    "Root Mean Square of the mapping quality of reads. Non-Variants."])      
        #5.4 GoodnessOfFit Variant NonVariant
        self.createSetPlots(stats_objects=[self.gofVariant,self.gofNonVariant],color="blue",
                            titles=["Phred scaled goodness of fit to the diploid model. Variants.",
                                    "Phred scaled goodness of fit to the diploid model. Non-Variants."])        
        
        #6. Mutation All 
        self.space()
        self.createTable(color="green",values=self.mutations.getTableMutationProfile(),tableTitle="Mutations")        
        #7. TiTv
        self.space()
        self.createTable(color="blue",values=self.mutations.getTiTvTable(),tableTitle="Ti/Tv Ratio")
      
        #8. Save HTML Document
        self.stop()
        self.saveReport()
        
class HtmlMethylation(HtmlBsCallReport):
    """ Class which BS-Genotype Calls """
    
    def __init__(self,html_file_name=None,current_name=None,parent_name=None,parent_document=None):
        """  Class constructor
        
             html_file_name -- HTML file document to be build
             current_name -- Current main name
             parent_name -- Parent name of the parent document
             parent_document -- Parent document
        """
        #Call Parent class constructor
        HtmlBsCallReport.__init__(self,html_file_name=html_file_name,parentName=parent_name,currentName=current_name,parentDocument=parent_document)

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
        '''Define HTML web documents'''
        #0. Mapping and Coverage Statistics
        self.header()
        #1. Total Stats Ref CpGs/NonRefCpG
        self.space()
        self.createTable(color="blue",values=self.total_stats.getTableMethylation(),tableTitle="CpG counts")
        #2. Coverage CpG and Coverage CpG Inf
        self.sectionTitle(title="Coverage CpGs")
        self.createSetPlots(stats_objects=[self.ref_cpg_coverage,self.ref_cpg_inf_coverage],color="green",titles=["Coverage Reference CpGs","Coverage Information Reads Reference CpGs"])        
        #3. Non Coverage CpG and Non Coverage CpG Inf
        self.createSetPlots(stats_objects=[self.non_ref_cpg_coverage,self.non_ref_cpg_inf_coverage],color="blue",titles=["Coverage Non Reference CpGs","Coverage Information Reads Non-Reference CpGs"])        
        #4. Quality RefCpG and Quality Non Reference CpG
        self.sectionTitle(title="Quality CpGs")
        self.createSetPlots(stats_objects=[self.quality_ref_cpg,self.quality_non_ref_cpg],color="green",titles=["Genotype Quality Phred scaled Reference CpGs","Genotype Quality Phred scaled Non-Reference CpGs"])        
        #5. Methylation distribution at CpGs
        self.sectionTitle(title="Methylation CpGs")
        self.createSetPlots(stats_objects=[self.methylation_plot_cpg],color="blue",titles=["CpGs Methylation Distribution"])    
        #6. Summary Methylation Profile
        self.createTable(color="green",values=self.summary_methylation.getTable(),tableTitle="CpG Methylation Profile (From Individual Strands) ")
        #7. Save HTML Document
        self.stop()
        self.saveReport()
        

class HtmlIndexBsCall(HtmlBsCallReport):
    """ Creates and Manages the Index HTML Document for BSCall Reports """
    

    def __init__(self,output_dir=None,name_project=None,samples_stats=None,samples_summary=None):
        """ Construct class for Creting Index HTML Document and BsCall Reports 
        
             output_dir -- Output directory to store HTML reports
             name_project -- Project name
             samples_stats -- Dictionary of samples and list of stats objects. [sample][ListStatsObject]
             samples_summary -- Dictionar of samples and summary fields. [sample][dictFieldSummary]   
                                                                         [dictFieldsSummary]{"alignments":ListAligmentAndCoverage,"variants":ListVariantsAndCoverage,"methylation":ListMethylation}
        """
        self.index_html_document = "%s/%s.html" %(output_dir,name_project)
        HtmlBsCallReport.__init__(self,html_file_name=self.index_html_document,parentName="",currentName="",parentDocument="")
        self.output_dir = output_dir
        self.name_project = name_project
        self.samples_stats = samples_stats
        self.samples_summary = samples_summary
        
    def createPage(self):
        #0. Index Header HTML  Document
        self.headerIndex(title="BScall Reports Project %s" %(self.name_project))
        #1. Sample Reports Building
        #1.0 Vector of Samples and its links to the reports
        samples_links = [["SAMPLE NAME","ALIGNMENTS AND COVERAGE REPORT","VARIANTS REPORT","METHYLATION REPORT"]]
        for sample in self.samples_stats:
            #1.1 Create Mapping and Coverage Statistics
            mappingCoverageHtml = "%s/%s_mapping_coverage.html" %(self.output_dir,sample) 
            parent_name = "%s/mapping_coverage" %(self.name_project)
            sample_mapping_coverage = HtmlMappingCoverage(html_file_name=mappingCoverageHtml,current_name=sample,parent_name=parent_name,parent_document=self.index_html_document)
            #1.1.2 Setup mapping coverage report
            sample_mapping_coverage.configureStats(stats_vector=self.samples_stats[sample]["mappingCoverage"])
            #1.1.3 Create Document
            sample_mapping_coverage.createPage()

            #1.2 Create Bs-Genotype Calls Report
            variantsHtml = "%s/%s_variants.html" %(self.output_dir,sample) 
            parent_variants = "%s/variants" %(self.name_project)
            sample_variants = HtmlBsGenotypeCalls(html_file_name=variantsHtml,current_name=sample,parent_name=parent_variants,parent_document=self.index_html_document)
            #1.2.1 Setup Variants Report
            sample_variants.configureStats(stats_vector=self.samples_stats[sample]["calls"])
            #1.2.2 Create Document
            sample_variants.createPage()

            #1.3 Create Methylation Statitics
            methylationHtml = "%s/%s_methylation.html" %(self.output_dir,sample)
            parent_methylation = "%s/methylation" %(self.name_project)           
            sample_methylation = HtmlMethylation(html_file_name=methylationHtml,current_name=sample,parent_name=parent_methylation,parent_document=self.index_html_document)
            #1.3.1 Setup Methylation Report
            sample_methylation.configureStats(stats_vector=self.samples_stats[sample]["methylation"])
            #1.3.2 Create Document
            sample_methylation.createPage()
            
            #1.4 Vector of links to create index Table
            samples_links.append([sample,os.path.basename(mappingCoverageHtml),os.path.basename(variantsHtml),os.path.basename(methylationHtml)])
        
        #2.Create Links Table
        self.createLinksTable(samples_links=samples_links,summary_field=self.samples_summary,color="blue")      
        
        #3. Save HTML Document
        self.stop()
        self.saveReport()
    
def buildBscallReports(inputs=None,output_dir=None,name=None):
    """ Build variant report.
    
        inputs -- Dictionary of samples and list of files. [sample] [barcode_chrN.json,barcode_chrN+1.json]
        output_dir -- Output directory to store reports.
        name --  Name basic to build output results.
    """
    #1. Check output directory
    if not os.path.exists("%s/IMG/" %(output_dir)):
        os.makedirs("%s/IMG/" %(output_dir)) 
        
    if not os.path.exists("%s/SPHINX/" %(output_dir)):
        os.makedirs("%s/SPHINX/" %(output_dir))

    #Proces list chromosome files
    dict_samples = {}
    dict_samples_summaries = {}
       
    for sample,chrom_json_files in inputs.iteritems():
        #Parsing json file
        readLevelStats = ReadsAndBases()
        baseLevelStats = BaseLevel()
        #Coverages
        allCoverage = Coverage(concept="All",yLabel="#Sites",pngFile="%s/IMG/%s_coverage_all.png" %(output_dir,sample) )
        variantCoverage = Coverage(concept="Variants",yLabel="#SNPs",pngFile="%s/IMG/%s_coverage_variants.png" %(output_dir,sample))
        dbSnpCoverage = Coverage(concept="dbSnp",yLabel="#SNPs",pngFile="%s/IMG/%s_coverage_dbsnp.png" %(output_dir,sample))
        refCpGcoverage = Coverage(concept="RefCpG",yLabel="#CpGs",pngFile="%s/IMG/%s_coverage_refCpG.png" %(output_dir,sample))
        refCpGInfCoverage = Coverage(concept="RefCpGInf",yLabel="#CpGs",pngFile="%s/IMG/%s_coverage_refCpGInf.png" %(output_dir,sample))
        nonRefCpGcoverage = Coverage(concept="NonRefCpG",yLabel="#CpGs",pngFile="%s/IMG/%s_coverage_nonRefCpG.png" %(output_dir,sample))
        nonRefCpGinfCoverage = Coverage(concept="NonRefCpGInf",yLabel="#CpGs",pngFile="%s/IMG/%s_coverage_nonRefCpGinf.png" %(output_dir,sample))
        #GC Coverage
        gcCoverage = GCcoverage("%s/IMG/%s_gc_coverage.png" %(output_dir,sample))
        #TotalStats
        totalStats = TotalStats()
        #Quality
        qualityAll = Quality(concept="All",yLabel="#Sites",pngFile="%s/IMG/%s_quality_all.png" %(output_dir,sample))
        qualityVariant = Quality(concept="Variants",yLabel="#SNPs",pngFile="%s/IMG/%s_quality_variant.png" %(output_dir,sample))
        qualityRefCpG = Quality(concept="RefCpG",yLabel="#CpGs",pngFile="%s/IMG/%s_quality_refcpg.png" %(output_dir,sample))
        qualityNonRefCpG = Quality(concept="NonRefCpG",yLabel="#CpGs",pngFile="%s/IMG/%s_quality_nonRefCpg.png" %(output_dir,sample))
        #QC Distributions
        fsVariant = QCDistribution(concept="FisherStrandVariant",typeDistribution="FisherStrand",typeBaseLocation="",pngFile="%s/IMG/%s_fs_variant.png" %(output_dir,sample))
        fsVariant.setAxisXLabel(newLabel="Fisher Strand Phred scale probability")
        qdVariant = QCDistribution(concept="QualityByDepthVariant",typeDistribution="QualityByDepth",typeBaseLocation="Variant",pngFile="%s/IMG/%s_qd_variant.png" %(output_dir,sample))
        qdNonVariant = QCDistribution(concept="QualityByDepthNonVariant",typeDistribution="QualityByDepth",typeBaseLocation="NonVariant",pngFile="%s/IMG/%s_qd_nonvariant.png" %(output_dir,sample))
        rmsmqVariant = QCDistribution(concept="RMSMappingQualityVariant",typeDistribution="RMSMappingQuality",typeBaseLocation="Variant",pngFile="%s/IMG/%s_rmsmq_variant.png" %(output_dir,sample))
        rmsmqNonVariant = QCDistribution(concept="RMSMappingQualityNonVariant",typeDistribution="RMSMappingQuality",typeBaseLocation="NonVariant",pngFile="%s/IMG/%s_rmsmq_nonvariant.png" %(output_dir,sample))
        gofVariant = QCDistribution(concept="GoodnessOfFitVariant",typeDistribution="GoodnessOfFit",typeBaseLocation="Variant",pngFile="%s/IMG/%s_gof_variant.png" %(output_dir,sample))
        gofNonVariant = QCDistribution(concept="GoodnessOfFitNonVariant",typeDistribution="GoodnessOfFit",typeBaseLocation="NonVariant",pngFile="%s/IMG/%s_gof_nonvariant.png" %(output_dir,sample))        
        #VCF Filter Stats
        vcfFilterStats = VCFFilterStats()
        #Mutations
        mutationsStats = Mutations("Mutations")
        #Methylation
        methylationAllRefCpG = Methylation("AllRefCpG")
        methylationPassRefCpG = Methylation("PassedRefCpG")
        methylationNonRefCpG = Methylation("AllNonRefCpG")
        methylationPassNonRefCpG = Methylation("PassedNonRefCpG")
        #NonCpGReadProfile
        nonCpGReadProfile = NonCpGReadProfile("%s/IMG/%s_nonCpgReadProfile.png" %(output_dir,sample))
        #SummaryMethylation
        summaryMethylation = SummaryMethylation() 
        
        #Load al json files
        for json_file in chrom_json_files: 
            with open(json_file, 'r') as file_json:
                try:
                    data = json.load(file_json)
        
                    readLevelStats.add(data["filterStats"]["ReadLevel"])
                    baseLevelStats.add(data["filterStats"]["BaseLevel"])
                    #Coverages
                    allCoverage.add(data["totalStats"]["coverage"]["All"])
                    variantCoverage.add(data["totalStats"]["coverage"]["Variant"])
                    if "dbSNP" in data["totalStats"]["coverage"]:
                        dbSnpCoverage.add(data["totalStats"]["coverage"]["dbSNP"])
                    refCpGcoverage.add(data["totalStats"]["coverage"]["RefCpG"])
                    refCpGInfCoverage.add(data["totalStats"]["coverage"]["RefCpGInf"])
                    nonRefCpGcoverage.add(data["totalStats"]["coverage"]["NonRefCpG"])
                    nonRefCpGinfCoverage.add(data["totalStats"]["coverage"]["NonRefCpGInf"])
                    #GC Coverage
                    gcCoverage.add(data["totalStats"]["coverage"]["GC"])
                    #TotalStats
                    totalStats.add(data["totalStats"])
                    #Quality
                    qualityAll.add(data["totalStats"]["quality"]["All"])
                    qualityVariant.add(data["totalStats"]["quality"]["Variant"])
                    qualityRefCpG.add(data["totalStats"]["quality"]["RefCpG"])
                    qualityNonRefCpG.add(data["totalStats"]["quality"]["NonRefCpG"])
                    #QC Distributions
                    fsVariant.add(data["totalStats"]["QCDistributions"]["FisherStrand"])
                    qdVariant.add(data["totalStats"]["QCDistributions"]["QualityByDepth"])
                    qdNonVariant.add(data["totalStats"]["QCDistributions"]["QualityByDepth"])
                    rmsmqVariant.add(data["totalStats"]["QCDistributions"]["RMSMappingQuality"])
                    rmsmqNonVariant.add(data["totalStats"]["QCDistributions"]["RMSMappingQuality"])
                    gofVariant.add(data["totalStats"]["QCDistributions"]["GoodnessOfFit"])
                    gofNonVariant.add(data["totalStats"]["QCDistributions"]["GoodnessOfFit"])  
                    #VCF Filter Stats
                    vcfFilterStats.add(data["totalStats"]["VCFFilterStats"])                
                    #Mutations
                    mutationsStats.add(data["totalStats"]["mutations"])              
                    #Methylation
                    methylationAllRefCpG.add(data["totalStats"]["methylation"]["AllRefCpg"])
                    methylationPassRefCpG.add(data["totalStats"]["methylation"]["PassedRefCpg"])
                    methylationNonRefCpG.add(data["totalStats"]["methylation"]["AllNonRefCpg"])
                    methylationPassNonRefCpG.add(data["totalStats"]["methylation"]["PassedNonRefCpg"])
                    #NonCpGReadProfile
                    nonCpGReadProfile.add(data["totalStats"]["methylation"]["NonCpGreadProfile"])

                except ValueError, e:
                    pass # invalid json

        

        #Prepare plot for Methylation levels
        plotMethylation = PlotMethylationLevels(concept="Methylation Levels",pngFile="%s/IMG/%s_methylation_levels.png" %(output_dir,sample),
                              meth_list=[methylationAllRefCpG,methylationPassRefCpG,methylationNonRefCpG,methylationPassNonRefCpG])

        #Get Table Summary
        summaryMethylation.setData(concept = "AllRefCpg",values=methylationAllRefCpG.methylation_cpgs)
        summaryMethylation.setData(concept = "PassedRefCpg",values=methylationPassRefCpG.methylation_cpgs)
        summaryMethylation.setData(concept = "AllNonRefCpg",values=methylationNonRefCpG.methylation_cpgs)
        summaryMethylation.setData(concept = "PassedNonRefCpg",values=methylationPassNonRefCpG.methylation_cpgs)
        
        #Sample Summary 
        #Prepare GC Values for plotting and getting GC Correlation
        gcCoverage.selectDataToPlot()                                              
        samples_summary = SummarySample(sampleName=sample,readLevelStats=readLevelStats,baseLevelStats=baseLevelStats,gcCoverage=gcCoverage,totalStats=totalStats,variantCoverage=variantCoverage,
                                       mutationStats=mutationsStats,methylationPassRefCpg=methylationPassRefCpG,refCpgCoverage=refCpGcoverage)  
                                       
        dict_samples_summaries[sample] = samples_summary.getTable()
        
        #DataSet Per Samples
        dict_samples[sample] = {"mappingCoverage": [readLevelStats,baseLevelStats,allCoverage,gcCoverage,qualityAll,nonCpGReadProfile],
                                "calls": [totalStats,vcfFilterStats,variantCoverage,dbSnpCoverage,qualityVariant,
                                           [fsVariant,qdVariant,qdNonVariant,rmsmqVariant,rmsmqNonVariant,gofVariant,gofNonVariant],
                                          mutationsStats
                                         ],
                                "methylation": [totalStats,refCpGcoverage,refCpGInfCoverage,nonRefCpGcoverage,nonRefCpGinfCoverage,
                                                qualityRefCpG,qualityNonRefCpG,plotMethylation,summaryMethylation]
                                }
                                                
                                            
                                               
            

    #Create Index BSCall Report
    htmlIndexBsCall = HtmlIndexBsCall(output_dir=output_dir,name_project=name,samples_stats=dict_samples,samples_summary=dict_samples_summaries)
    htmlIndexBsCall.createPage()
    
    #CSS object
    cssBuilder = BasicHtml()
    vector_css = []
    cssBuilder.buildStyleSheet(vector_css)
    RunBasicStats.saveDocument(file_name="%s/style.css" %(output_dir),vectorContent=vector_css) 
    
    #Create Index BSCall Report
    sphinxIndexBsCall = SphinxIndexBsCall(output_dir="%s/SPHINX/"%(output_dir),name_project=name,samples_stats=dict_samples)
    sphinxIndexBsCall.createPage()
    
    #Config python file
    cfgFile = ConfigSphinx(path_config_file="%s/SPHINX/conf.py" %(output_dir),path_makefile_file="%s/SPHINX/Makefile" %(output_dir),master_file=name,project_name=name,main_title='BSCALL REPORT')
    cfgFile.run()