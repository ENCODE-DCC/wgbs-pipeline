# -*- coding: utf-8 -*-
#!/usr/bin/env python

import os
import json

from reportStats import NucleotideStats,LaneStats,SampleStats,RunBasicStats

"""gemBS parsers JSON files to build an HTML report"""
class BasicHtml(RunBasicStats):
    """ Class responsable of basic html functions """

    def __init__(self,mapping_stats=None,png_mapq_histogram=None,png_insert_size_histogram=None):
        #Call Parent class constructor
        #super(RunBasicStats, self).__init__() 
        RunBasicStats.__init__(self,mapping_stats=mapping_stats,png_mapq_histogram=png_mapq_histogram,png_insert_size_histogram=png_insert_size_histogram)
        
    def run(self,vectorHtml=None):
        '''To Be Defined in the Child Class'''
        pass


    def addHtmlReportHeader(self,vectorHtml=None):
        ''' Add HTML Report Header to a given vector'''
        vectorHtml.append("<HTML>\n")

        vectorHtml.append(" <HEAD>\n")
        vectorHtml.append(" <STYLE TYPE=\"text/css\">\n")
        vectorHtml.append("  <!--\n")
        vectorHtml.append("   @import url(\"style.css\"); \n")
        vectorHtml.append("  -->\n")
        vectorHtml.append(" </STYLE>\n")
        vectorHtml.append(" </HEAD>\n")
        vectorHtml.append(" <BODY>\n")
        
    def closeHtmlReport(self,vectorHtml=None):
        ''' Close HTML header to a given vector object '''  
        vectorHtml.append(" </BODY>\n")
        vectorHtml.append("</HTML>\n")
        
    def addRow(self,vectorHtml=None,statValue=None,name_concept=None,isOdd=True):
        """ Add row html tags
            
            vectorHtml - Vector of tags to add
            statValue - Value to be built
            name_concept - Name of the concept
        """
        statValue.sumPair()
        
        if isOdd == True: 
            vectorHtml.append("   <TR class=\"odd\">\n")
        else:
            vectorHtml.append("   <TR>\n")
        if self.is_paired:
            vectorHtml.append("   <TD> %s </TD> <TD> %s </TD> <TD> %.2f %% </TD> <TD> %s </TD> <TD> %.2f %% </TD> <TD> %s </TD> <TD> %.2f %% </TD> \n" \
                              %(name_concept,statValue.total,statValue.total_percentage,statValue.pair_one,statValue.pair_one_percentage,statValue.pair_two,statValue.pair_two_percentage))
                               
        else:        
            vectorHtml.append("   <TD> %s </TD> <TD> %s </TD> <TD> %.2f %%</TD> \n" %(name_concept,statValue.total,statValue.total_percentage))
        vectorHtml.append("   </TR>\n")
        
    def addRowSimpleValue(self,vectorHtml=None,value=None,name_concept=None,isOdd=True):
        """ Add Simple Value Row html tags
            
            vectorHtml - Vector of tags to add
            value - stats value
            name_concept - Name of the concept
        """
        if isOdd == True: 
            vectorHtml.append("   <TR class=\"odd\">\n")
        else:
            vectorHtml.append("   <TR>\n")
        
        vectorHtml.append("   <TD> %s </TD> <TD> %s </TD> \n" %(name_concept,value))
        vectorHtml.append("   </TR>\n")
       
 
    def tableDefinition(self,vectorHtml=None,color=None):
        """ HTML Table tag definition
        
        vectorHtml - List of html to add table tags
        color - Table color, could be green or blue
        """
        if color == "green":
            vectorHtml.append("  <TABLE id=\"green\">\n")
        else:
            vectorHtml.append("  <TABLE id=\"hor-zebra\">\n")
            
    def createStatsTable(self,color=None):
        """Create Statistics Table
        
           color - Table color, could be green or blue  
           returns vector of HTML tags
        """
        vTableHtml = []   
     
        #1.Table Definition
        self.tableDefinition(vectorHtml=vTableHtml,color=color)

        #2.Table header
        vTableHtml.append("   <TR>\n")
        vTableHtml.append("    <TH scope=\"col\">Concept</TH> <TH scope=\"col\">Total Reads</TH> <TH scope=\"col\">%</TH>\n")
        if self.is_paired:
            vTableHtml.append("    <TH scope=\"col\">Pair One Reads</TH> <TH scope=\"col\">%</TH> <TH scope=\"col\">Pair Two Reads</TH> <TH scope=\"col\">%</TH> \n")        
        vTableHtml.append("   </TR>\n")
        
        #3.Stats table
        self.createStatsForTableReport()
        #Sequenced Reads
        self.addRow(vTableHtml,self.mapping_stats.sequenced_reads,"Sequenced Reads",True)
        #General reads 
        self.addRow(vTableHtml,self.mapping_stats.reads,"General Reads",False)
        #Reads in Control sequences
        self.addRow(vTableHtml,self.mapping_stats.reads_sequencing_control,"Reads in Control sequences",True)  
        #Reads under conversion control  
        self.addRow(vTableHtml,self.mapping_stats.reads_under_conversion_control,"Reads under conversion control ",False)
        #Reads over conversion control
        self.addRow(vTableHtml,self.mapping_stats.reads_over_conversion_control,"Reads over conversion control",True)                  
        #Unmapped reads
        self.addRow(vTableHtml,self.mapping_stats.reads_unmapped,"Unmapped reads",False) 
        #Bisulfite_reads
        self.addRow(vTableHtml,self.mapping_stats.reads_reference_C2T,"Bisulfite_reads C2T",True) 
        self.addRow(vTableHtml,self.mapping_stats.reads_reference_G2A,"Bisulfite_reads G2A",False)
        
        vTableHtml.append(" </TABLE>\n")
        return vTableHtml
        
    def createUniqueFragmentsTable(self,color=None,unique_fragments=None,average_unique=None):
        """Create Unique and Fragments Statistics Table
        
           color - Table color, could be green or blue  
           unique_fragments - Fragments that were considere as unique
           average_unique - Average unique fragments
           returns vector of HTML tags
        """
        vTableHtml = []   

        #1.Table Definition
        self.tableDefinition(vectorHtml=vTableHtml,color=color)

        #2.Table header
        vTableHtml.append("   <TR>\n")
        vTableHtml.append("    <TH scope=\"col\">Concept</TH> <TH scope=\"col\">Value</TH>\n")
        vTableHtml.append("   </TR>\n")
        
        #3.Stats table
        #Uniquely mapped reads
        vTableHtml.append('  <TR>')
        vTableHtml.append("   <TD> %s </TD> <TD> %s </TD> \n" %("Unique Fragments",unique_fragments))
        vTableHtml.append('  </TR>')
        vTableHtml.append('  <TR class="odd">')
        vTableHtml.append("   <TD> %s </TD> <TD> %.2f %% </TD> \n" %("Average Unique",average_unique))
        vTableHtml.append("  </TR>\n")  

        vTableHtml.append(" </TABLE>\n")
        return vTableHtml
        
        
        
    def createBasesStatsTable(self,color=None):
        """Create Bases Statistics Table
        
           color - Table color, could be green or blue  
           returns vector of HTML tags
        """
        vTableHtml = []   
        
        #1.Table Definition
        self.tableDefinition(vectorHtml=vTableHtml,color=color)

        #2.Table header
        vTableHtml.append("   <TR>\n")
        vTableHtml.append("    <TH scope=\"col\">Concept</TH> <TH scope=\"col\">Total Bases</TH> <TH scope=\"col\">%</TH>\n")
        if self.is_paired:
            vTableHtml.append("    <TH scope=\"col\">Pair One Bases</TH> <TH scope=\"col\">%</TH> <TH scope=\"col\">Pair Two Bases</TH> <TH scope=\"col\">%</TH> \n")        
        vTableHtml.append("   </TR>\n")
        
        #Table Counts
        self.createBasesStats() 
        
        #Base counts Overall
        self.addRow(vTableHtml,self.mapping_stats.reads_overall_A,"Base Counts Overall A",True)
        self.addRow(vTableHtml,self.mapping_stats.reads_overall_C,"Base Counts Overall C",False)
        self.addRow(vTableHtml,self.mapping_stats.reads_overall_G,"Base Counts Overall G",True)
        self.addRow(vTableHtml,self.mapping_stats.reads_overall_T,"Base Counts Overall T",False)
        self.addRow(vTableHtml,self.mapping_stats.reads_overall_N,"Base Counts Overall N",True)

        #Base counts GeneralC2T
        #self.addRow(vTableHtml,self.mapping_stats.reads_general_C2T_A,"Base Counts GeneralC2T A",False)
        #self.addRow(vTableHtml,self.mapping_stats.reads_general_C2T_C,"Base Counts GeneralC2T C",True)
        #self.addRow(vTableHtml,self.mapping_stats.reads_general_C2T_G,"Base Counts GeneralC2T G",False)
        #self.addRow(vTableHtml,self.mapping_stats.reads_general_C2T_T,"Base Counts GeneralC2T T",True)
        #self.addRow(vTableHtml,self.mapping_stats.reads_general_C2T_N,"Base Counts GeneralC2T N",False)

        #Base counts GeneralG2A
        #self.addRow(vTableHtml,self.mapping_stats.reads_general_G2A_A,"Base Counts GeneralG2A A",True)
        #self.addRow(vTableHtml,self.mapping_stats.reads_general_G2A_C,"Base Counts GeneralG2A C",False)
        #self.addRow(vTableHtml,self.mapping_stats.reads_general_G2A_G,"Base Counts GeneralG2A G",True)
        #self.addRow(vTableHtml,self.mapping_stats.reads_general_G2A_T,"Base Counts GeneralG2A T",False)
        #self.addRow(vTableHtml,self.mapping_stats.reads_general_G2A_N,"Base Counts GeneralG2A N",True)

        #Base Counts UnderConversionControlC2T
        #self.addRow(vTableHtml,self.mapping_stats.reads_under_conversion_control_C2T_A,"Base Counts UnderConversionControlC2T A",False)
        #self.addRow(vTableHtml,self.mapping_stats.reads_under_conversion_control_C2T_C,"Base Counts UnderConversionControlC2T C",True)
        #self.addRow(vTableHtml,self.mapping_stats.reads_under_conversion_control_C2T_G,"Base Counts UnderConversionControlC2T G",False)
        #self.addRow(vTableHtml,self.mapping_stats.reads_under_conversion_control_C2T_T,"Base Counts UnderConversionControlC2T T",True)
        #self.addRow(vTableHtml,self.mapping_stats.reads_under_conversion_control_C2T_N,"Base Counts UnderConversionControlC2T N",False)

        #Base Counts UnderConversionControlG2A
        #self.addRow(vTableHtml,self.mapping_stats.reads_under_conversion_control_G2A_A,"Base Counts UnderConversionControlG2A A",True)
        #self.addRow(vTableHtml,self.mapping_stats.reads_under_conversion_control_G2A_C,"Base Counts UnderConversionControlG2A C",False)
        #self.addRow(vTableHtml,self.mapping_stats.reads_under_conversion_control_G2A_G,"Base Counts UnderConversionControlG2A G",True)
        #self.addRow(vTableHtml,self.mapping_stats.reads_under_conversion_control_G2A_T,"Base Counts UnderConversionControlG2A T",False)
        #self.addRow(vTableHtml,self.mapping_stats.reads_under_conversion_control_G2A_N,"Base Counts UnderConversionControlG2A N",True)

        #Base Counts OverConversionControl_C2T
        #self.addRow(vTableHtml,self.mapping_stats.reads_over_conversion_control_C2T_A,"Base Counts OverConversionControlC2T A",False)
        #self.addRow(vTableHtml,self.mapping_stats.reads_over_conversion_control_C2T_C,"Base Counts OverConversionControlC2T C",True)
        #self.addRow(vTableHtml,self.mapping_stats.reads_over_conversion_control_C2T_G,"Base Counts OverConversionControlC2T G",False)
        #self.addRow(vTableHtml,self.mapping_stats.reads_over_conversion_control_C2T_T,"Base Counts OverConversionControlC2T T",True)
        #self.addRow(vTableHtml,self.mapping_stats.reads_over_conversion_control_C2T_N,"Base Counts OverConversionControlC2T N",False)
        
        #Base Counts OverConversionControl_G2A
        #self.addRow(vTableHtml,self.mapping_stats.reads_over_conversion_control_G2A_A,"Base Counts OverConversionControlG2A A",True)
        #self.addRow(vTableHtml,self.mapping_stats.reads_over_conversion_control_G2A_C,"Base Counts OverConversionControlG2A C",False)
        #self.addRow(vTableHtml,self.mapping_stats.reads_over_conversion_control_G2A_G,"Base Counts OverConversionControlG2A G",True)
        #self.addRow(vTableHtml,self.mapping_stats.reads_over_conversion_control_G2A_T,"Base Counts OverConversionControlG2A T",False)
        #self.addRow(vTableHtml,self.mapping_stats.reads_over_conversion_control_G2A_N,"Base Counts OverConversionControlG2A N",True)
        
        vTableHtml.append(" </TABLE>\n")
        return vTableHtml



    def createOverlappedBasesTable(self,color='green',bases_overlapped=0,total_bases=0,average_overlapped_bases=0):
        """ Create Overlapped Bases Table
        
            color - Table color, coould be green or blue
            bases_overlapped -- Total estimated overlapping bases
            total_bases -- Total number of bases 
            average_overlapped_bases -- average overlapping bases
            returns vector of HTML tags
        """
        vTableHtml = []
        
        #1.Table Definition
        self.tableDefinition(vectorHtml=vTableHtml,color=color)
        
        #2. Table header
        vTableHtml.append("   <TR>\n")
        vTableHtml.append("    <TH scope=\"col\">Overlapped Bases (Properly Aligned) </TH> <TH scope=\"col\">Total Bases (Properly Aligned)</TH> <TH scope=\"col\">% Overlapping Bases (Properly Aligned)</TH>\n")
        vTableHtml.append("   </TR>\n")
        
        #3.Stats table
        #Read Values  
        vTableHtml.append('  <TR>')
        vTableHtml.append("   <TD> %i </TD>  <TD> %i </TD> <TD> %.2f %%</TD> \n" %(bases_overlapped,total_bases,average_overlapped_bases))
        vTableHtml.append('  </TR>')

        vTableHtml.append(" </TABLE>\n")
        
        return vTableHtml  
        
        
    def createBisulfiteConversionRate(self,color='blue'):
        """Create Bisulfite Conversion Rate

           color - Table color, could be green or blue
           returns vector of HTML tags
        """
        vTableHtml = []   
     
        #1.Table Definition
        self.tableDefinition(vectorHtml=vTableHtml,color=color)

        #2.Table header
        vTableHtml.append("   <TR>\n")
        vTableHtml.append("    <TH scope=\"col\">Bisulfite Conversion Type</TH> <TH scope=\"col\">Conversion Rate</TH>\n")
        vTableHtml.append("   </TR>\n")
        
        #3.Stats table
        #Read Values  
        vTableHtml.append('  <TR>')
        vTableHtml.append("   <TD> %s </TD> <TD> %s </TD> \n" %("Under Conversion Rate",self.mapping_stats.getUnderConversionRate()))
        vTableHtml.append('  </TR>')
        vTableHtml.append('  <TR class="odd">')
        vTableHtml.append("   <TD> %s </TD> <TD> %s </TD> \n" %("Over Conversion Rate",self.mapping_stats.getOverConversionRate()))
        vTableHtml.append("  </TR>\n")
        vTableHtml.append(" </TABLE>\n")
        
        return vTableHtml  
        

    def createSingleValueTable(self,color=None,name=None,value=None):
        """Create Single Value Table
        
           color - Table color, could be green or blue  
           returns vector of HTML tags
        """
        vTableHtml = []   
     
        #1.Table Definition
        self.tableDefinition(vectorHtml=vTableHtml,color=color)

        #2.Table header
        vTableHtml.append("   <TR>\n")
        vTableHtml.append("    <TH scope=\"col\">Concept</TH> <TH scope=\"col\">Total Reads</TH>\n")
        vTableHtml.append("   </TR>\n")
        
        #3.Stats table
        #Read Values  
        vTableHtml.append('  <TR>')
        vTableHtml.append("   <TD> %s </TD> <TD> %s </TD> \n" %(name,value))
        vTableHtml.append("   </TR>\n")
        vTableHtml.append(" </TABLE>\n")
        
        return vTableHtml    
        
 
    def createLinksTables(self,vector_html=None,vector_links=None,title=None,color=None):
        ''' Create a set of link tables '''
        self.tableDefinition(vectorHtml=vector_html,color=color)

        vector_html.append("   <TR> <TH scope=\"col\"> %s </TH> </TR> \n" %(title))

        linksClass = ''	
        if color == "green":
            linksClass += '"linkGreen"'
        else:
            linksClass += '"link"'
  
        isOdd = True
        for sampleLink in vector_links:
            if isOdd == True: 
                vector_html.append('   <TR class="odd"> <TD> <a class=%s href="%s"> %s </TD> </TR> \n' %(linksClass,sampleLink[1],sampleLink[0]))
            else:
                vector_html.append('   <TR> <TD> <a class=%s href="%s"> %s </TD> </TR> \n' %(linksClass,sampleLink[1],sampleLink[0]))

            isOdd = not isOdd
            
        vector_html.append(" </TABLE>\n")
        
        
    def createLinksSumupSampleTables(self,vector_html=None,vector_links=None,samplesStats=None,color=None):
        ''' Create a set of link tables '''
        self.tableDefinition(vectorHtml=vector_html,color=color)
        
        linksClass = ''	
        if color == "green":
            linksClass += '"linkGreen"'
        else:
            linksClass += '"link"'

        vector_html.append("   <TR> <TH scope=\"col\"> Sample </TH>  \
                                    <TH scope=\"col\"> Yield Reads </TH>  \
                                    <TH scope=\"col\"> Uniquely Mapped Fragments </TH>  \
                                    <TH scope=\"col\"> Average Unique </TH>  \
                                    <TH scope=\"col\"> Under Conversion Rate </TH>  \
                                    <TH scope=\"col\"> Over Conversion Rate </TH>  \
                                    <TH scope=\"col\"> Report </TH>  \
                                    </TR> \n")
	
        isOdd = True
        for sampleLink,sampleStats in zip(vector_links,samplesStats):
            
            if isOdd == True: 
                vector_html.append('   <TR class="odd"> \n')
            else:
                vector_html.append('   <TR>\n')

            vector_html.append("    <TD> %s </TD>  \
                                    <TD> %s </TD>  \
                                    <TD> %s </TD>  \
                                    <TD> %.2f %% </TD>  \
                                    <TD> %s </TD>  \
                                    <TD> %s </TD>  \
                                    <TD> <a class=%s href=\"%s\"> &#187 %s </TD>  \
                                    </TR> \n" %(sampleStats.name,sampleStats.sequenced_reads.total,sampleStats.totalSampleUniqueReads,
                                                 sampleStats.averageSampleUniqueReads,sampleStats.getUnderConversionRate(),
                                                 sampleStats.getOverConversionRate(),
                                                 linksClass,sampleLink,sampleStats.name
                                                 ))
            
            isOdd = not isOdd
            
        vector_html.append(" </TABLE>\n")

    def addSpaceSection(self,vectorHtml=None):
        ''' Add HTML Report Header to a given vector'''
        vectorHtml.append("<BR><BR><BR>\n")
        
    def addSectionTitle(self,vectorHtml=None,title=None):
        ''' Add HTML Report Header to a given vector'''
        vectorHtml.append('<H1 id="section"> %s </H1>\n' %(title))

    def buildStyleSheet(self,vector_css=None):
        ''' Add HTML header to a given vector '''
        vector_css.append(" #title\n")
        vector_css.append(" {\n")
        vector_css.append("   font-family: \"Sans-Serif\";\n")
        vector_css.append("   font-size: 20px;\n")
        vector_css.append("   text-align: center;\n")
        vector_css.append("   color: #039;\n")
        vector_css.append("   border-collapse: collapse;\n")
        vector_css.append(" }\n")
        vector_css.append(" #section\n")
        vector_css.append(" {\n")
        vector_css.append("   font-family: \"Sans-Serif\";\n")
        vector_css.append("   font-size: 18px;\n")
        vector_css.append("   text-align: center;\n")
        vector_css.append("   color: #039;\n")
        vector_css.append("   border-collapse: collapse;\n")
        vector_css.append(" }\n")
        #LINKS TABLE
        vector_css.append("#linksTable-b\n")
        vector_css.append("{\n")
        vector_css.append("   font-family: \"Sans-Serif\";\n")
        vector_css.append("	font-size: 12px;\n")
        vector_css.append("	background: #fff;\n")
        vector_css.append("	margin: 5px;\n")
        vector_css.append("	width: 200px;\n")
        vector_css.append("	border-collapse: collapse;\n")
        vector_css.append("	text-align: left;\n")
        vector_css.append("}\n")
        vector_css.append("#linksTable-b th\n")
        vector_css.append("{\n")
        vector_css.append("	font-size: 14px;\n")
        vector_css.append("	font-weight: normal;\n")
        vector_css.append("	color: #039;\n")
        vector_css.append("	padding: 10px 8px;\n")
        vector_css.append("	border-bottom: 2px solid #6678b1;\n")
        vector_css.append("}\n")
        vector_css.append("#linksTable-b td\n")
        vector_css.append("{\n")
        vector_css.append("	border-bottom: 1px solid #ccc;\n")
        vector_css.append("	color: #669;\n")
        vector_css.append("	padding: 6px 8px;\n")
        vector_css.append("}\n")
        vector_css.append("linksTable-b tbody tr:hover td\n")
        vector_css.append("{\n")
        vector_css.append("	color: #009;\n")
        vector_css.append("}\n")
        #BLUE TABLE
        vector_css.append(" #hor-zebra\n")
        vector_css.append(" {\n")
        vector_css.append("   font-family: \"Sans-Serif\";\n")
        vector_css.append("   font-size: 12px;\n")
        vector_css.append("   margin: 0px;\n")
        vector_css.append("   width: 100%;\n")
        vector_css.append("   text-align: left;\n")
        vector_css.append("   border-collapse: collapse;\n")
        vector_css.append(" }\n")
        vector_css.append(" #hor-zebra th\n")
        vector_css.append(" {\n")
        vector_css.append("   font-size: 14px;\n")
        vector_css.append("   font-weight: normal;\n")
        vector_css.append("   padding: 10px 8px;\n")
        vector_css.append("   color: #039;\n")
        vector_css.append("   border-bottom: 2px solid #6678b1;\n")
        vector_css.append("   border-right: 1px solid #6678b1; \n")
        vector_css.append("	  border-left: 1px solid #6678b1;\n")
        vector_css.append(" }\n")
        vector_css.append(" #hor-zebra td\n")
        vector_css.append(" {\n")
        vector_css.append("   padding: 8px;\n")
        vector_css.append("   color: #669;\n")
        vector_css.append("   border-right: 1px solid #6678b1; \n")
        vector_css.append("	  border-left: 1px solid #6678b1;\n")
        vector_css.append(" }\n")
        vector_css.append(" #hor-zebra .odd\n")
        vector_css.append(" {\n")
        vector_css.append("   background: #e8edff;\n")
        vector_css.append("   border-right: 1px solid #6678b1; \n")
        vector_css.append("	  border-left: 1px solid #6678b1;\n")
        vector_css.append(" }\n")
        vector_css.append(" #hor-zebra .empty\n")
        vector_css.append(" {\n")
        vector_css.append("   border-bottom: 2px solid #6678b1;\n")
        vector_css.append(" }\n") 
        vector_css.append(" #hor-zebra .oddvalue\n")
        vector_css.append(" {\n")
        vector_css.append("   background: #f5f7ff;\n")
        vector_css.append("   border-right: 1px solid #6678b1; \n")
        vector_css.append("	  border-left: 1px solid #6678b1;\n")
        vector_css.append(" }\n")
        vector_css.append(" #hor-zebra .pairheader\n")
        vector_css.append(" {\n")
        vector_css.append("   background: #eeffe8;\n")
        vector_css.append("   border-right: 1px solid #66b17b; \n")
        vector_css.append("	  border-left: 1px solid #66b17b;\n")
        vector_css.append(" }\n")
        vector_css.append(" #hor-zebra .pairvalue\n")
        vector_css.append(" {\n")
        vector_css.append("   background: #f6fff4;\n")
        vector_css.append("   border-right: 1px solid #66b17b; \n")
        vector_css.append("	  border-left: 1px solid #66b17b;\n")
        vector_css.append(" }\n")
        #GREEN TABLE
        vector_css.append(" #green \n")
        vector_css.append(" {\n")
        vector_css.append("   font-family: \"Sans-Serif\";\n")
        vector_css.append("   font-size: 12px;\n")
        vector_css.append("   margin: 0px;\n")
        vector_css.append("   width: 100%;\n")
        vector_css.append("   text-align: left;\n")
        vector_css.append("   border-collapse: collapse;\n")
        vector_css.append(" }\n")
        vector_css.append(" #green th\n")
        vector_css.append(" {\n")
        vector_css.append("   font-size: 14px;\n")
        vector_css.append("   font-weight: normal;\n")
        vector_css.append("   padding: 10px 8px;\n")
        vector_css.append("   color: #2b9900;\n")
        vector_css.append("   border-bottom: 2px solid #66b16f;\n")
        vector_css.append("   border-right: 1px solid #66b16f; \n")
        vector_css.append("	  border-left: 1px solid #66b16f;\n")
        vector_css.append(" }\n")
        vector_css.append(" #green td\n")
        vector_css.append(" {\n")
        vector_css.append("   padding: 8px;\n")
        vector_css.append("   color: #578055;\n")
        vector_css.append("   border-right: 1px solid #66b16f; \n")
        vector_css.append("	  border-left: 1px solid #66b16f;\n")
        vector_css.append(" }\n")
        vector_css.append(" #green .odd\n")
        vector_css.append(" {\n")
        vector_css.append("   background: #eaffe8;\n")
        vector_css.append("   border-right: 1px solid #66b16f; \n")
        vector_css.append("	  border-left: 1px solid #66b16f;\n")
        vector_css.append(" }\n")
        vector_css.append(" #green .empty\n")
        vector_css.append(" {\n")
        vector_css.append("   border-bottom: 2px solid #66b16f;\n")
        vector_css.append(" }\n")          
        #LINKS
        vector_css.append("a.link:link {font-family: \"Sans-Serif\";font-size: 12px;color: #039;text-decoration:none;}")
        vector_css.append("a.link:visited {font-family: \"Sans-Serif\";font-size: 12px;color: #039;text-decoration:none;}")
        vector_css.append("a.link:hover {font-family: \"Sans-Serif\";font-size: 12px;color: #039;text-decoration:underline;}")
        #LINKS GREEN
        vector_css.append("a.linkGreen:link {font-family: \"Sans-Serif\";font-size: 12px;color: #2b9900;text-decoration:none;}")
        vector_css.append("a.linkGreen:visited {font-family: \"Sans-Serif\";font-size: 12px;color: #2b9900;text-decoration:none;}")
        vector_css.append("a.linkGreen:hover {font-family: \"Sans-Serif\";font-size: 12px;color: #2b9900;text-decoration:underline;}")
        #P BLUE
        vector_css.append(" #descriptionBlue \n")
        vector_css.append(" {\n")
        vector_css.append("   font-family: \"Sans-Serif\";\n")
        vector_css.append("   font-size: 14px;\n")
        vector_css.append("   text-align: left;\n")
        vector_css.append("   color: #039;\n")
        vector_css.append("   border-collapse: collapse;\n")
        vector_css.append(" }\n")
        #P GREEN
        vector_css.append(" #descriptionGreen \n")
        vector_css.append(" {\n")
        vector_css.append("   font-family: \"Sans-Serif\";\n")
        vector_css.append("   font-size: 14px;\n")
        vector_css.append("   text-align: left;\n")
        vector_css.append("   color: #009900;\n")
        vector_css.append("   border-collapse: collapse;\n")
        vector_css.append(" }\n")
        #P SUBTITLE BLUE
        vector_css.append(" #subtitleBlue \n")
        vector_css.append(" {\n")
        vector_css.append("   font-family: \"Sans-Serif\";\n")
        vector_css.append("   font-size: 16px;\n")
        vector_css.append("   font-weight: bold;\n")
        vector_css.append("   text-decoration:underline;\n")
        vector_css.append("   text-align: left;\n")
        vector_css.append("   color: #039;\n")
        vector_css.append("   border-collapse: collapse;\n")
        vector_css.append(" }\n")
        #P SUBTITLE GREEN
        vector_css.append(" #subtitleGreen \n")
        vector_css.append(" {\n")
        vector_css.append("   font-family: \"Sans-Serif\";\n")
        vector_css.append("   font-size: 16px;\n")
        vector_css.append("   font-weight: bold;\n")
        vector_css.append("   text-decoration:underline;\n")
        vector_css.append("   text-align: left;\n")
        vector_css.append("   color: #009900;\n")
        vector_css.append("   border-collapse: collapse;\n")
        vector_css.append(" }\n")
        #P SOLID
        vector_css.append("#path { \n")
        vector_css.append("   border-style: solid; \n")
        vector_css.append("   border-width: 2px; \n")
        vector_css.append("   border-color: #039; \n")
        vector_css.append('   font-family: "Sans-Serif"; \n')
        vector_css.append("   color: #039;\n")
        vector_css.append("   font-size: 12px;\n")
        vector_css.append(" }\n")
                

class LaneHtml(BasicHtml):
    """ Class responsable for printing html report per Lane """
    
    def __init__(self,project_name=None,sample_name=None,lane_stats=None,png_mapq_histogram=None,png_insert_size_histogram=None,html_parent_path=None):
        """ Lane Html manage the html report building
        
            project_name - Project name
            sample_name - Sample Name to which belongs the lane
            lane_stats - LaneStats object to build HTML report
            png_mapq_histogram - File png histogram for the MAPQ plot
            png_insert_size_histogram - File png histogram for insert size histogram plot
            html_parent_path - HTML parent document
        """
        #Call Parent class constructor
        #super(LaneHtml, self).__init__()
        BasicHtml.__init__(self,mapping_stats=lane_stats,png_mapq_histogram=png_mapq_histogram,png_insert_size_histogram=png_insert_size_histogram)
        
        self.project_name = project_name
        self.sample_name = sample_name
        self.html_parent_path = html_parent_path
        
    def run(self,vectorHtml=None):
        """ Run Lane HTML Documentation Building
        
            vectorHtml - vector of html tags
        """        
        #1. Add HTML Report Header
        self.addHtmlStartReport(vectorHtml=vectorHtml)
  
        self.addSectionTitle(vectorHtml=vectorHtml,title="Mapping Stats (Reads)")
        vectorHtml.extend(self.createStatsTable(color='blue'))

        self.addSpaceSection(vectorHtml=vectorHtml) 
        self.addSectionTitle(vectorHtml=vectorHtml,title="Uniqueness (Fragments)")
        vectorHtml.extend(self.createUniqueFragmentsTable(color='green',unique_fragments=self.mapping_stats.getUniqueMappedReads(),
                          average_unique=self.mapping_stats.getAverageUniqueMappedReads()))
        
        self.addSpaceSection(vectorHtml=vectorHtml)
        self.addSectionTitle(vectorHtml=vectorHtml,title="Mapping Stats (Bases)")
        vectorHtml.extend(self.createBasesStatsTable(color='blue'))
        
        self.addSpaceSection(vectorHtml=vectorHtml)
        self.addSectionTitle(vectorHtml=vectorHtml,title="Bisulfite Conversion Rate")
        vectorHtml.extend(self.createBisulfiteConversionRate(color='blue'))

        if self.is_paired:
            self.addSpaceSection(vectorHtml=vectorHtml) 
            self.addSectionTitle(vectorHtml=vectorHtml,title="Correct Pairs")
            vectorHtml.extend(self.createSingleValueTable(color='green',name='Correct Pairs',value=self.mapping_stats.correct_pairs))
  
        self.addSpaceSection(vectorHtml=vectorHtml) 
        self.addSectionTitle(vectorHtml=vectorHtml,title="Mapping Quality")
        vectorHtml.extend(self.buildMapqHistogram(color='blue'))
  
        self.addSpaceSection(vectorHtml=vectorHtml) 
        self.addSectionTitle(vectorHtml=vectorHtml,title="Read Length")
        vectorHtml.extend(self.buildReadsLengthTable(color='green'))
        
        if self.is_paired:
            self.addSpaceSection(vectorHtml=vectorHtml) 
            self.addSectionTitle(vectorHtml=vectorHtml,title="Insert Size Plot")
            vectorHtml.extend(self.buildInsertSizePlot(color='blue'))
            self.closeHtmlReport(vectorHtml=vectorHtml)
                       
    def addHtmlStartReport(self,vectorHtml=None):
        ''' Add HTML Report Header to a given vector'''
        self.addHtmlReportHeader(vectorHtml)

        vectorHtml.append("\n <P id='path'> /%s/%s/%s </P> \n" %(self.project_name,self.sample_name,self.mapping_stats.name))
        vectorHtml.append("\n <a class=\"link\" href=\"%s\"><B>BACK</B></a> <br>\n" %(os.path.basename(self.html_parent_path)))        
        vectorHtml.append("  <H1 id=\"title\"> <U> SAMPLE %s LANE %s </U> </H1>\n" %(self.sample_name,self.mapping_stats.name))
               
    def buildMapqHistogram(self,color=None):
        """ From matplot lib plots a Mappping qualty histogram
        
            color - Color for the table were the plot is goind to be located
        """
        #1. PLOT BUILDING
        self.drawMapqHistogram()
                
        #2. HTML TAG BUILDING 
        vTableHtml = []
   
        self.tableDefinition(vectorHtml=vTableHtml,color=color)

        vTableHtml.append("   <TR class=\"odd\"> <TH scope=\"col\"> Mapping Quality Histogram</TH> </TR> \n")
        vTableHtml.append('   <TR> <TD> <img src="%s" alt="%s"> </TD> </TR> \n' %(os.path.basename(self.png_mapq_histogram),os.path.basename(self.png_mapq_histogram)))
        vTableHtml.append("  </TABLE>\n")
        
        return vTableHtml
        
    def buildReadsLengthTable(self,color=None):
        """ Build Read Length Table
        
            color - Color for the table were the plot is going to be located
        """
        
        vTableHtml = []
        #1 Table Definition
        self.tableDefinition(vectorHtml=vTableHtml,color=color)
        
        vTableHtml.append("   <TR class=\"odd\"> <TH scope=\"col\"> Read Length </TH> <TH scope=\"col\"> Reads </TH> </TR> \n")
        
        isOdd = False
        
        for dictValue in self.mapping_stats.read_length_histogram:
            for readLength, reads in dictValue.iteritems():
                self.addRowSimpleValue(vectorHtml=vTableHtml,value=reads,name_concept=readLength,isOdd=isOdd)
                isOdd = not isOdd
                
        vTableHtml.append("  </TABLE>\n")
         
        return vTableHtml
        
    def buildInsertSizePlot(self,color=None):
        """Build Insert Size plot and locate on html report
        
          color - Color for the table were the plot is going to be located
        """
        vTableHtml = []
        
        #1. PLOT BUILDING
        self.drawInsertSizePlot()
        
        #2. HTML TAG BUILDING 
        vTableHtml = []
   
        self.tableDefinition(vectorHtml=vTableHtml,color=color)

        vTableHtml.append("   <TR class=\"odd\"> <TH scope=\"col\"> Insert Size Histogram</TH> </TR> \n")
        vTableHtml.append('   <TR> <TD> <img src="%s" alt="%s"> </TD> </TR> \n' %(os.path.basename(self.png_insert_size_histogram),os.path.basename(self.png_insert_size_histogram)))
        vTableHtml.append("  </TABLE>\n")
        
        return vTableHtml        
        
        
class SampleHtml(BasicHtml):
    """ Class to manage Sample report, builds a html report per sample and call for Lane html report building """

    def __init__(self,project_name=None,sample_stats=None,html_parent_path=None,html_sample=None,png_insert_size_histogram=None,png_mapq_histogram=None):
        """ SampleHtml constructor
        
            project_name - Project name
            sample_stats - SampleStats object from which create html report
            html_parent_path - HTML parent document
            html_sample - HTML sample document name
            png_insert_size_histogram - File png histogram for insert size histogram plot
        """
        #Call Parent class constructor
        #super(SampleHtml, self).__init__()    
        BasicHtml.__init__(self,mapping_stats=sample_stats,png_mapq_histogram=png_mapq_histogram,png_insert_size_histogram=png_insert_size_histogram)        
        
        self.project_name = project_name
        
        self.is_paired = sample_stats.is_paired
                
        self.html_parent_path = html_parent_path
        self.html_sample = html_sample
        
    def run(self,vectorHtml=None):
        """ Run Lane HTML Documentation Building
        
            vectorHtml - vector of html tags
        """    
        
         #1. Add HTML Report Header
        self.addHtmlStartReport(vectorHtml=vectorHtml)
        self.addSpaceSection(vectorHtml=vectorHtml) 
        self.addSectionTitle(vectorHtml=vectorHtml,title="Mapping Stats (Reads)")
        vectorHtml.extend(self.createStatsTable(color='blue'))
        
        self.addSpaceSection(vectorHtml=vectorHtml) 
        self.addSectionTitle(vectorHtml=vectorHtml,title="Uniqueness (Fragments)")
        vectorHtml.extend(self.createUniqueFragmentsTable(color='green',unique_fragments=self.mapping_stats.totalSampleUniqueReads,
                          average_unique=self.mapping_stats.averageSampleUniqueReads))
        
        self.addSpaceSection(vectorHtml=vectorHtml)
        self.addSectionTitle(vectorHtml=vectorHtml,title="Mapping Stats (Bases)")
        vectorHtml.extend(self.createBasesStatsTable(color='blue'))

        self.addSpaceSection(vectorHtml=vectorHtml)
        self.addSectionTitle(vectorHtml=vectorHtml,title="Bisulfite Conversion Rate")
        vectorHtml.extend(self.createBisulfiteConversionRate(color='blue'))        
             
        if self.is_paired:
            self.addSpaceSection(vectorHtml=vectorHtml) 
            self.addSectionTitle(vectorHtml=vectorHtml,title="Correct Pairs")
            vectorHtml.extend(self.createSingleValueTable(color='green',name='Correct Pairs',value=self.mapping_stats.correct_pairs))
        
        #2. Lane Stats and links per lane
        vLaneLinksVector = []
        for laneStats in self.mapping_stats.list_lane_stats:
            htmlLane = "%s/%s.html" %(os.path.dirname(self.html_parent_path),laneStats.name)
            mapqHistogram = "%s/%s.mapq.png" %(os.path.dirname(self.html_parent_path),laneStats.name)
            isizeHistogram = "%s/%s.isize.png" %(os.path.dirname(self.html_parent_path),laneStats.name)
            vLaneHtml = []
            laneHtml = LaneHtml(project_name=self.project_name,sample_name=self.mapping_stats.name,lane_stats=laneStats,png_mapq_histogram=mapqHistogram,png_insert_size_histogram=isizeHistogram,html_parent_path=self.html_sample)
            laneHtml.run(vLaneHtml)
            RunBasicStats.saveDocument(file_name=htmlLane,vectorContent=vLaneHtml)
            vLaneLinksVector.append([laneStats.name,os.path.basename(htmlLane)])
          
        #3.Mapping Quality
        self.addSpaceSection(vectorHtml=vectorHtml) 
        self.addSectionTitle(vectorHtml=vectorHtml,title="Mapping Quality")
        vectorHtml.extend(self.buildMapqHistogram(color='blue'))
     
        #4. Insert Size Plot    
        if self.is_paired:
            self.addSpaceSection(vectorHtml=vectorHtml) 
            self.addSectionTitle(vectorHtml=vectorHtml,title="Insert Size Plot")
            vectorHtml.extend(self.buildInsertSizePlot(color='green'))          
          
        #5.Links table
        self.addSpaceSection(vectorHtml=vectorHtml) 
        self.addSectionTitle(vectorHtml=vectorHtml,title="Mapping Lanes Reports")
        self.createLinksTables(vectorHtml,vLaneLinksVector,title="LANE REPORTS",color="blue")
        #6. Close html
        self.closeHtmlReport(vectorHtml=vectorHtml)
            
            
    def addHtmlStartReport(self,vectorHtml=None):
        ''' Add HTML Report Header to a given vector'''
        self.addHtmlReportHeader(vectorHtml)
        
        vectorHtml.append("\n <P id='path'> /%s/%s </P> \n" %(self.project_name,self.mapping_stats.name))
        vectorHtml.append("\n <a class=\"link\" href=\"%s\"><B>BACK</B></a> <br>\n" %(os.path.basename(self.html_parent_path)))        
        vectorHtml.append("  <H1 id=\"title\"> <U> SAMPLE %s </U> </H1>\n" %(self.mapping_stats.name))
   
    def buildMapqHistogram(self,color=None):
        """ From matplot lib plots a Mappping qualty histogram
        
            color - Color for the table were the plot is goind to be located
        """
        #1. PLOT BUILDING
        self.drawMultipleMapqHistogram()
                
        #2. HTML TAG BUILDING 
        vTableHtml = []
   
        self.tableDefinition(vectorHtml=vTableHtml,color=color)

        vTableHtml.append("   <TR class=\"odd\"> <TH scope=\"col\"> Mapping Quality Histogram</TH> </TR> \n")
        vTableHtml.append('   <TR> <TD> <img src="%s" alt="%s"> </TD> </TR> \n' %(os.path.basename(self.png_mapq_histogram),os.path.basename(self.png_mapq_histogram)))
        vTableHtml.append("  </TABLE>\n")
        
        return vTableHtml     
        
    def buildInsertSizePlot(self,color=None):
        """Build Insert Size plot and locate on html report
        
          color - Color for the table were the plot is going to be located
        """
        #1. PLOT BUILDING
        self.drawMultipleInsertSizePlot()
        
        #2. HTML TAG BUILDING 
        vTableHtml = []
   
        self.tableDefinition(vectorHtml=vTableHtml,color=color)

        vTableHtml.append("   <TR class=\"odd\"> <TH scope=\"col\"> Insert Size Histogram</TH> </TR> \n")
        vTableHtml.append('   <TR> <TD> <img src="%s" alt="%s"> </TD> </TR> \n' %(os.path.basename(self.png_insert_size_histogram),os.path.basename(self.png_insert_size_histogram)))
        vTableHtml.append("  </TABLE>\n")
        
        return vTableHtml 
        
                
class IndexHtml(BasicHtml):
    """ Class which defines Index Sample Report """
    
    def __init__(self,output_dir=None,name_project=None,vector_samples=None):
        """  Class constructor
        
             output_dir -- Output directory to store HTML reports
             name_project -- Project name
        """
        #Call Parent class constructor
        BasicHtml.__init__(self)        
        
        self.output_dir = output_dir
        self.name_project = name_project
        self.project_html_document = "%s/%s.html" %(self.output_dir,self.name_project)
        self.vector_samples = vector_samples
        
        
    def run(self,vectorHtml=None):
        """ Run Lane HTML Documentation Building
        
            vectorHtml - vector of html tags
        """ 
        #1. Add HTML Report Header
        self.addHtmlReportHeader(vectorHtml)
        
        vectorHtml.append('<H1 id="title"> <U> Methylation Pipeline Report Project %s </U> </H1>' %(self.name_project))
        
        vector_sample_links = []
        
        for sampleStats in self.vector_samples:
            sampleHtml = "%s/%s.html" %(self.output_dir,sampleStats.name)
            isizeHistogram = "%s/%s.isize.png" %(self.output_dir,sampleStats.name)
            png_mapq_histogram = "%s/%s.mapq.png" %(self.output_dir,sampleStats.name)
            sampleReport = SampleHtml(project_name=self.name_project,sample_stats=sampleStats,html_parent_path=self.project_html_document,\
                                      html_sample=sampleHtml,png_insert_size_histogram=isizeHistogram,png_mapq_histogram=png_mapq_histogram)
            vSampleHtml = []
            sampleReport.run(vSampleHtml)
            RunBasicStats.saveDocument(file_name=sampleHtml,vectorContent=vSampleHtml)
            vector_sample_links.append(os.path.basename(sampleHtml))
            
            
        #3.Links table
        self.createLinksSumupSampleTables(vector_html=vectorHtml,vector_links=vector_sample_links,samplesStats=self.vector_samples)        
        
        #4.Close Html
        self.closeHtmlReport(vectorHtml=vectorHtml)
        
        
def buildReport(inputs=None,output_dir=None,name=None):
    """ Build report per lane and sample.
    
        inputs -- Dictionary of samples and lanes [sample][fli]json_file
        output_dir -- Output directory to store html documents.
        name --  Name basic to build output results.
    """
      
    #Check output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)      
      
    #Proces list Lane files
    vector_samples = []
    for sample,fli_json in inputs.iteritems():
        list_stats_lanes = []
        for fli,json_files in fli_json.iteritems():  
            for json_file in json_files:
                lane = LaneStats(name=fli,json_file=json_file)
                list_stats_lanes.append(lane)
            
        vector_samples.append(SampleStats(name=sample,list_lane_stats=list_stats_lanes))    

    #IndexHtml object
    vector_index_html = []
    indexHtml = IndexHtml(output_dir=output_dir,name_project=name,vector_samples=vector_samples)
    indexHtml.run(vectorHtml=vector_index_html)
    RunBasicStats.saveDocument(file_name=indexHtml.project_html_document,vectorContent=vector_index_html)
    
    #CSS object
    cssBuilder = BasicHtml()
    vector_css = []
    cssBuilder.buildStyleSheet(vector_css)
    RunBasicStats.saveDocument(file_name="%s/style.css" %(output_dir),vectorContent=vector_css)
    
    
    

            
    
    
