# -*- coding: utf-8 -*-
#!/usr/bin/env python

import os
from reportStats import NucleotideStats,LaneStats,SampleStats,RunBasicStats

class BasicSphinx(RunBasicStats):
    """ Class responsable of basic Sphinx functions """

    def __init__(self,mapping_stats=None,png_mapq_histogram=None,png_insert_size_histogram=None):
        #Call Parent class constructor
        RunBasicStats.__init__(self,mapping_stats=mapping_stats,png_mapq_histogram=png_mapq_histogram,png_insert_size_histogram=png_insert_size_histogram)
        
    def run(self,vectorSphinx=None):
        '''To Be Defined in the Child Class'''
        pass
    
    def addLine(self,ident=20,vectorSphinx=None,cells=0,lenCell=30,char="-"):
        """ Write Line
        
            ident - Characters to ident
            vectorSphinx - Vector of Sphinx Contents
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
        vectorSphinx.append(''.join(underLine))    
     
    
    def addSimpleLine(self,ident=20,vectorSphinx=None,cells=0,lenCell=30):
        """ Underline
        
            ident - Characters to ident
            vectorSphinx - Vector of Sphinx Contents
            cells - Number of cell in the rows
            isLeft - length of each cell
        """       
        self.addLine(ident=ident,vectorSphinx=vectorSphinx,cells=cells,lenCell=lenCell,char="-")
        
    def addHeaderLine(self,ident=20,vectorSphinx=None,cells=0,lenCell=30):
        """ Underline
        
            ident - Characters to ident
            vectorSphinx - Vector of Sphinx Contents
            cells - Number of cell in the rows
            isLeft - length of each cell
        """       
        self.addLine(ident=ident,vectorSphinx=vectorSphinx,cells=cells,lenCell=lenCell,char="=")

    def addCell(self,text="",lenCell=30,isLeft=True,isRight=False):
        """ Generate Cell
        
            ident - Characters to ident
            vectorSphinx - Vector of Sphinx Contents
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
                    

    def addRow(self,ident=20,lenCell=30,vectorSphinx=None,statValue=None,name_concept=None):
        """ Add row sphinx
            
            ident - Characters to ident
            lenCell - Length of the Cell in each row
            vectorSphinx - Vector of tags to add
            statValue - Value to be built
            name_concept - Name of the concept
        """
        statValue.sumPair()
        cells = 0
        cellContents = ""
        #Add ident chars
        for i in range(ident):
            cellContents += " "        
        
        if self.is_paired:
            cellContents += self.addCell(text="%s" %(name_concept),lenCell=lenCell,isLeft=True,isRight=False)
            cells += 1
            cellContents += self.addCell(text="%s" %(statValue.total),lenCell=lenCell,isLeft=False,isRight=False)
            cells += 1
            cellContents += self.addCell(text="%.2f %%" %(statValue.total_percentage),lenCell=lenCell,isLeft=False,isRight=False)
            cells += 1
            cellContents += self.addCell(text="%s" %(statValue.pair_one),lenCell=lenCell,isLeft=False,isRight=False)
            cells += 1
            cellContents += self.addCell(text="%.2f %%" %(statValue.pair_one_percentage),lenCell=lenCell,isLeft=False,isRight=False)
            cells += 1
            cellContents += self.addCell(text="%s" %(statValue.pair_two),lenCell=lenCell,isLeft=False,isRight=False)
            cells += 1
            cellContents += self.addCell(text="%.2f %%" %(statValue.pair_two_percentage),lenCell=lenCell,isLeft=False,isRight=True)
            cells += 1   
            vectorSphinx.append(cellContents)                           
        else: 
            cellContents += self.addCell(text="%s" %(name_concept),lenCell=lenCell,isLeft=True,isRight=False)
            cells += 1
            cellContents += self.addCell(text="%s" %(statValue.total),lenCell=lenCell,isLeft=False,isRight=False)
            cells += 1
            cellContents += self.addCell(text="%.2f %%" %(statValue.total_percentage),lenCell=lenCell,isLeft=False,isRight=True)
            cells += 1   
            vectorSphinx.append(cellContents) 
            
        #Underline
        self.addSimpleLine(ident=ident,vectorSphinx=vectorSphinx,cells=cells,lenCell=lenCell)


    def addRowSimpleValue(self,ident=20,lenCell=30,vectorSphinx=None,name_concept=None,value=None):
        """ Add Simple Value Row Sphinx tags
            
            ident - Characters to ident
            lenCell - Length of the Cell in each row
            vectorSphinx - Vector of tags to add
            statValue - Value to be built
            name_concept - Name of the concept
        """
        cells = 0
        cellContents = ""
        #Add ident chars
        for i in range(ident):
            cellContents += " "
            
        cellContents += self.addCell(text="%s" %(name_concept),lenCell=lenCell,isLeft=True,isRight=False)
        cells += 1
        cellContents += self.addCell(text="%s" %(value),lenCell=lenCell,isLeft=False,isRight=True)
        cells += 1   
        vectorSphinx.append("%s" %(cellContents))

        #Underline
        self.addSimpleLine(ident=ident,vectorSphinx=vectorSphinx,cells=cells,lenCell=lenCell)
        
            
    def createStatsTable(self,ident=20,lenCell=30):
        """Create Statistics Table
        
           ident - Characters to ident 
           lenCell - Cell length
           returns vector of Sphinx tags
        """
        vTableSphinx = []   
        cellsHeader = 0
        
        #1.Table header
        if self.is_paired:
            cellsHeader = 7            
            #2. Write overline
            self.addSimpleLine(ident=ident,vectorSphinx=vTableSphinx,cells=cellsHeader,lenCell=lenCell)

            #3.Stats table
            cellContents = ""
            cellContents += self.addCell(text="Concept",lenCell=lenCell,isLeft=True,isRight=False)
            cellContents += self.addCell(text="Total Reads",lenCell=lenCell,isLeft=False,isRight=False)
            cellContents += self.addCell(text="%",lenCell=lenCell,isLeft=False,isRight=False)
            cellContents += self.addCell(text="P1 Reads",lenCell=lenCell,isLeft=False,isRight=False)
            cellContents += self.addCell(text="%",lenCell=lenCell,isLeft=False,isRight=False)
            cellContents += self.addCell(text="P2 Reads",lenCell=lenCell,isLeft=False,isRight=False)
            cellContents += self.addCell(text="%",lenCell=lenCell,isLeft=False,isRight=True)
            vTableSphinx.append(cellContents)                
        else:
            cellsHeader = 3
            #2. Write overline
            self.addSimpleLine(ident=ident,vectorSphinx=vTableSphinx,cells=cellsHeader,lenCell=lenCell)
            
            #3.Stats table
            cellContents = ""
            cellContents += self.addCell(text="Concept",lenCell=lenCell,isLeft=True,isRight=False)
            cellContents += self.addCell(text="Total Reads",lenCell=lenCell,isLeft=False,isRight=False)
            cellContents += self.addCell(text="%",lenCell=lenCell,isLeft=False,isRight=False)
            vTableSphinx.append(cellContents)
            
        #Header Line 
        self.addHeaderLine(ident=ident,vectorSphinx=vTableSphinx,cells=cellsHeader,lenCell=lenCell)   

        #2. Table Contents
        self.createStatsForTableReport()
        #Sequenced Reads
        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.sequenced_reads,name_concept="Sequenced Reads")
        #General reads 
        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads,name_concept="General Reads")
        #Reads in Control sequences
        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_sequencing_control,name_concept="Reads in Control sequences")
        #Reads under conversion control  
        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_under_conversion_control,name_concept="Reads under conversion control")
        #Reads over conversion control
        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_over_conversion_control,name_concept="Reads over conversion control")
        #Unmapped reads
        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_unmapped,name_concept="Unmapped reads")
        #Bisulfite_reads
        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_reference_C2T,name_concept="Bisulfite_reads C2T")
        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_reference_G2A,name_concept="Bisulfite_reads G2A")

        return vTableSphinx
        
    def createUniqueFragmentsTable(self,ident=20,lenCell=30,unique_fragments=None,average_unique=None):
        """Create Unique and Fragments Statistics Table
        
           ident - Characters to ident 
           lenCell - Cell length  
           unique_fragments - Fragments that were considere as unique
           average_unique - Average unique fragments
           returns vector of Sphinx tags
        """
        vTableSphinx = []

        #0. OverLine
        self.addSimpleLine(ident=ident,vectorSphinx=vTableSphinx,cells=2,lenCell=lenCell)        
        
        #1.Table header
        cellContents = ""
        cellContents += self.addCell(text="Concept",lenCell=lenCell,isLeft=True,isRight=False)
        cellContents += self.addCell(text="Value",lenCell=lenCell,isLeft=False,isRight=True)
        vTableSphinx.append(cellContents)
        
        #Header Line 
        self.addHeaderLine(ident=ident,vectorSphinx=vTableSphinx,cells=2,lenCell=lenCell)
        
        #2.Stats table
        #Uniquely mapped reads
        self.addRowSimpleValue(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,name_concept="Unique Fragments",value=unique_fragments)
        self.addRowSimpleValue(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,name_concept="Average Unique",value="%.2f %%" %(average_unique))

        return vTableSphinx
        
        
    def createBasesStatsTable(self,ident=20,lenCell=30):
        """Create Bases Statistics Table
        
           ident - Characters to ident 
           lenCell - Cell length  
           returns vector of Sphinx tags
        """
        vTableSphinx = []   

        #2.Table header
        cellContents = ""
        cells = 0
        cellContents += self.addCell(text="Concept",lenCell=lenCell,isLeft=True,isRight=False)
        cells += 1
        cellContents += self.addCell(text="Total Bases",lenCell=lenCell,isLeft=False,isRight=False)
        cells += 1
        cellContents += self.addCell(text="%",lenCell=lenCell,isLeft=False,isRight=not self.is_paired)
        cells += 1
        if self.is_paired:
            cellContents += self.addCell(text="P1 Bases",lenCell=lenCell,isLeft=False,isRight=False)
            cells += 1
            cellContents += self.addCell(text="%",lenCell=lenCell,isLeft=False,isRight=False)
            cells += 1
            cellContents += self.addCell(text="P2 Bases",lenCell=lenCell,isLeft=False,isRight=False)
            cells += 1
            cellContents += self.addCell(text="%",lenCell=lenCell,isLeft=False,isRight=True)
            cells += 1
        
        #Overline
        self.addSimpleLine(ident=ident,vectorSphinx=vTableSphinx,cells=cells,lenCell=lenCell)
        vTableSphinx.append(cellContents)
        #Header Line 
        self.addHeaderLine(ident=ident,vectorSphinx=vTableSphinx,cells=cells,lenCell=lenCell)
                
        #Table Counts
        self.createBasesStats()        
        
        #Base counts Overall
        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_overall_A,name_concept="Base Counts Overall A")
        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_overall_C,name_concept="Base Counts Overall C")
        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_overall_G,name_concept="Base Counts Overall G")
        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_overall_T,name_concept="Base Counts Overall T")
        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_overall_N,name_concept="Base Counts Overall N")

        #Base counts GeneralC2T
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_general_C2T_A,name_concept="Base Counts GeneralC2T A")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_general_C2T_C,name_concept="Base Counts GeneralC2T C")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_general_C2T_G,name_concept="Base Counts GeneralC2T G")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_general_C2T_T,name_concept="Base Counts GeneralC2T T")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_general_C2T_N,name_concept="Base Counts GeneralC2T N")

        #Base counts GeneralG2A
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_general_G2A_A,name_concept="Base Counts GeneralG2A A")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_general_G2A_C,name_concept="Base Counts GeneralG2A C")          
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_general_G2A_G,name_concept="Base Counts GeneralG2A G")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_general_G2A_T,name_concept="Base Counts GeneralG2A T")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_general_G2A_N,name_concept="Base Counts GeneralG2A N")
        
        #Base Counts UnderConversionControlC2T
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_under_conversion_control_C2T_A,name_concept="Base Counts UnderConversionControlC2T A")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_under_conversion_control_C2T_C,name_concept="Base Counts UnderConversionControlC2T C")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_under_conversion_control_C2T_G,name_concept="Base Counts UnderConversionControlC2T G")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_under_conversion_control_C2T_T,name_concept="Base Counts UnderConversionControlC2T T")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_under_conversion_control_C2T_N,name_concept="Base Counts UnderConversionControlC2T N")

        #Base Counts UnderConversionControlG2A
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_under_conversion_control_G2A_A,name_concept="Base Counts UnderConversionControlG2A A")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_under_conversion_control_G2A_C,name_concept="Base Counts UnderConversionControlG2A C")   
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_under_conversion_control_G2A_G,name_concept="Base Counts UnderConversionControlG2A G")   
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_under_conversion_control_G2A_T,name_concept="Base Counts UnderConversionControlG2A T")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_under_conversion_control_G2A_N,name_concept="Base Counts UnderConversionControlG2A N")
        
        #Base Counts OverConversionControl_C2T
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_over_conversion_control_C2T_A,name_concept="Base Counts OverConversionControlC2T A")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_over_conversion_control_C2T_C,name_concept="Base Counts OverConversionControlC2T C")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_over_conversion_control_C2T_G,name_concept="Base Counts OverConversionControlC2T G")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_over_conversion_control_C2T_T,name_concept="Base Counts OverConversionControlC2T T")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_over_conversion_control_C2T_N,name_concept="Base Counts OverConversionControlC2T N")
                
        #Base Counts OverConversionControl_G2A
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_over_conversion_control_G2A_A,name_concept="Base Counts OverConversionControlG2A A")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_over_conversion_control_G2A_C,name_concept="Base Counts OverConversionControlG2A C")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_over_conversion_control_G2A_G,name_concept="Base Counts OverConversionControlG2A G")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_over_conversion_control_G2A_T,name_concept="Base Counts OverConversionControlG2A T")
#        self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=self.mapping_stats.reads_over_conversion_control_G2A_N,name_concept="Base Counts OverConversionControlG2A N")
        
        return vTableSphinx
                
    def createBisulfiteConversionRate(self,ident=20,lenCell=30):
        """Create Bisulfite Conversion Rate

           ident - Characters to ident 
           lenCell - Cell length  
           returns vector of Sphinx tags
        """
        vTableSphinx = []   
     
        #0. Overline
        self.addSimpleLine(ident=ident,vectorSphinx=vTableSphinx,cells=2,lenCell=lenCell)
        
        #1.Table header
        cellContents = ""
        cellContents += self.addCell(text="Bisulfite Conversion Type",lenCell=lenCell,isLeft=True,isRight=False)
        cellContents += self.addCell(text="Conversion Rate",lenCell=lenCell,isLeft=False,isRight=True)
        vTableSphinx.append(cellContents)
              
        #1. Header Line 
        self.addHeaderLine(ident=ident,vectorSphinx=vTableSphinx,cells=2,lenCell=lenCell)

        #2.Stats table
        #Read Values  
        self.addRowSimpleValue(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,name_concept="Under Conversion Rate",value=self.mapping_stats.getUnderConversionRate())
        self.addRowSimpleValue(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,name_concept="Over Conversion Rate",value=self.mapping_stats.getOverConversionRate())
       
        return vTableSphinx  

    def createSingleValueTable(self,ident=20,lenCell=30,name=None,value=None):
        """Create Single Value Table
        
           ident - Characters to ident 
           lenCell - Cell length  
           returns vector of Sphinx tags
        """
        vTableSphinx = []   
     
        #0. Overline
        self.addSimpleLine(ident=ident,vectorSphinx=vTableSphinx,cells=2,lenCell=lenCell)
        
        #1.Table header
        cellContents = ""
        cellContents += self.addCell(text="Concept",lenCell=lenCell,isLeft=True,isRight=False)
        cellContents += self.addCell(text="Total Reads",lenCell=lenCell,isLeft=False,isRight=True)
        vTableSphinx.append(cellContents)
        
        #1. Header Line 
        self.addHeaderLine(ident=ident,vectorSphinx=vTableSphinx,cells=2,lenCell=lenCell)
        
        #2.Stats table
        #Read Values  
        self.addRowSimpleValue(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,name_concept=name,value=value)

        return vTableSphinx    
        

    def createSumupSampleTables(self,ident=20,lenCell=30,vector_sphinx=None,samplesStats=None):
        ''' Create a set of link tables '''


        #1.Table header
        cellContents = ""
        
        #Add ident chars
        for i in range(ident):
            cellContents += " "
            
        numCells = 0
        cellContents += self.addCell(text="Sample",lenCell=lenCell,isLeft=True,isRight=False)
        numCells += 1
        cellContents += self.addCell(text="Yield Reads",lenCell=lenCell,isLeft=False,isRight=False)
        numCells += 1
        cellContents += self.addCell(text="Uniquely Mapped Fragments",lenCell=lenCell,isLeft=False,isRight=False)
        numCells += 1
        cellContents += self.addCell(text="Average Unique",lenCell=lenCell,isLeft=False,isRight=False)
        numCells += 1
        cellContents += self.addCell(text="Under Conversion Rate",lenCell=lenCell,isLeft=False,isRight=False)
        numCells += 1
        cellContents += self.addCell(text="Over Conversion Rate",lenCell=lenCell,isLeft=False,isRight=True)
        numCells += 1
        
        #Header Line 
        self.addSimpleLine(ident=ident,vectorSphinx=vector_sphinx,cells=numCells,lenCell=lenCell)      
        vector_sphinx.append("%s"%(cellContents))
        self.addHeaderLine(ident=ident,vectorSphinx=vector_sphinx,cells=numCells,lenCell=lenCell)

        #3. Table Contents
        for sampleStats in samplesStats:
            cellContents = ""
            
            #Add ident chars
            for i in range(ident):
                cellContents += " "
            
            numCells = 0
            cellContents += self.addCell(text="%s" %(sampleStats.name),lenCell=lenCell,isLeft=True,isRight=False)
            numCells += 1
            cellContents += self.addCell(text="%s" %(sampleStats.sequenced_reads.total),lenCell=lenCell,isLeft=False,isRight=False)
            numCells += 1
            cellContents += self.addCell(text="%s" %(sampleStats.totalSampleUniqueReads),lenCell=lenCell,isLeft=False,isRight=False)
            numCells += 1
            cellContents += self.addCell(text="%s" %(format(sampleStats.averageSampleUniqueReads, '.2f')),lenCell=lenCell,isLeft=False,isRight=False)
            numCells += 1
            cellContents += self.addCell(text="%s" %(sampleStats.getUnderConversionRate()),lenCell=lenCell,isLeft=False,isRight=False)
            numCells += 1
            cellContents += self.addCell(text="%s" %(sampleStats.getOverConversionRate()),lenCell=lenCell,isLeft=False,isRight=True)
            numCells += 1
            vector_sphinx.append("%s"%(cellContents))
            self.addSimpleLine(ident=ident,vectorSphinx=vector_sphinx,cells=numCells,lenCell=lenCell)
            

        
    def addLevelSection(self,ident=20,vectorSphinx=None,title=None,char="="):
        """ Write Level Section
        
            ident - Characters to ident
            vectorSphinx - Vector of Sphinx Contents
            title - Section title
            char - Type of character
        """ 
        #1.Title Line
        titleLine = ""       
         
        for i in range(ident):
            titleLine += " "
        titleLine += title
        vectorSphinx.append("%s\n" %(titleLine))
        
        #2.Section Line
        sectionLine = ""
        
        for i in range(ident):
            sectionLine += " "
    
        for i in range(len(title)):
            sectionLine += char
        
        #3. Append Section Line
        vectorSphinx.append("%s\n" %(sectionLine))
        
    def addTopSection(self,ident=20,vectorSphinx=None,title=None):
        """ Write Top Section
        
            ident - Characters to ident
            vectorSphinx - Vector of Sphinx Contents
            title - Section title
            char - Type of character
        """ 
        self.addLevelSection(ident=ident,vectorSphinx=vectorSphinx,title=title,char="=")
    
    def addSubSection(self,ident=20,vectorSphinx=None,title=None):
        """ Write Sub Section
        
            ident - Characters to ident
            vectorSphinx - Vector of Sphinx Contents
            title - Section title
            char - Type of character
        """ 
        self.addLevelSection(ident=ident,vectorSphinx=vectorSphinx,title=title,char="-")
        
    def addSubSubSection(self,ident=20,vectorSphinx=None,title=None):
        """ Write Sub sub section
        
            ident - Characters to ident
            vectorSphinx - Vector of Sphinx Contents
            title - Section title
            char - Type of character
        """
        self.addLevelSection(ident=ident,vectorSphinx=vectorSphinx,title=title,char="^")
        
    def addLevelHighlight(self,ident=20,vectorSphinx=None,title=None,char="*"):
        """ Write Level Section
        
            ident - Characters to ident
            vectorSphinx - Vector of Sphinx Contents
            title - Section title
            char - Type of character
        """
        #1. Over Line
        overLine = ""
        for i in range(ident):
            overLine += " "          
        for i in range(len(title)):
            overLine += char
        vectorSphinx.append("%s\n" %(overLine))
              
        #2.Title Line
        titleLine = ""              
        for i in range(ident):
            titleLine += " "
        titleLine += title
        vectorSphinx.append("%s\n" %(titleLine))
        
        #3.Section Line
        sectionLine = ""    
        for i in range(ident):
            sectionLine += " "
        for i in range(len(title)):
            sectionLine += char
        
        #3. Append Section Line
        vectorSphinx.append("%s\n" %(sectionLine))
        
    def addPartLine(self,ident=20,vectorSphinx=None,title=None):
        """ Add Part Line (Main Highline)
        
            ident - Characters to ident
            vectorSphinx - Vector of Sphinx Contents
            title - Section title
            char - Type of character
        """
        self.addLevelHighlight(ident=ident,vectorSphinx=vectorSphinx,title=title,char="#")
        
    def addChapterLine(self,ident=20,vectorSphinx=None,title=None):
        """ Add Chapter Line (Chapter Highline)
        
            ident - Characters to ident
            vectorSphinx - Vector of Sphinx Contents
            title - Section title
            char - Type of character
        """
        self.addLevelHighlight(ident=ident,vectorSphinx=vectorSphinx,title=title,char="*")
               
    def addTocTree(self,vectorSphinx=None,vectorDocuments=None):
        """ Add Toc Tree Documents

            vectorSphinx - Vector to add Sphinx tags
            vectorDocuments - List of documents to write at the index
        """
        content =   ".. toctree::\n"
        content +=  "   :maxdepth: 3\n"
        content +=  "\n"
        for document in vectorDocuments:
            content += "   %s\n" %(document)
        vectorSphinx.append(content)


class LaneSphinx(BasicSphinx):
    """ Class responsable for printing Sphinx report per Lane """
    
    def __init__(self,project_name=None,sample_name=None,lane_stats=None,png_mapq_histogram=None,png_insert_size_histogram=None):
        """ Lane Sphinx manage the html report building
        
            project_name - Project name
            sample_name - Sample Name to which belongs the lane
            lane_stats - LaneStats object to build HTML report
            png_mapq_histogram - File png histogram for the MAPQ plot
            png_insert_size_histogram - File png histogram for insert size histogram plot
        """
        #Call Parent class constructor
        BasicSphinx.__init__(self,mapping_stats=lane_stats,png_mapq_histogram=png_mapq_histogram,png_insert_size_histogram=png_insert_size_histogram)        
        
        self.project_name = project_name
        self.sample_name = sample_name
        
    def run(self,vectorSphinx=None,ident=20,lenCell=30):
        """ Run Lane SPHINX Documentation Building
        
            ident - Characters to ident
            lenCell - Length of table cells in characters
            vectorSphinx - vector of Sphinx tags
        """  
        #0. Add Chapter Line       
        self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Sample %s Lane %s" %(self.sample_name,self.mapping_stats.name))        
        vectorSphinx.append("\n")        
        
        #1. Add SPHINX Report Header
        self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Mapping Stats (Reads)")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.createStatsTable(ident=ident,lenCell=lenCell))
        vectorSphinx.append("\n")
        
        self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Uniqueness (Fragments)")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.createUniqueFragmentsTable(ident=ident,lenCell=lenCell,unique_fragments=self.mapping_stats.getUniqueMappedReads(),
                                                          average_unique=self.mapping_stats.getAverageUniqueMappedReads()))
        vectorSphinx.append("\n")                                                          

        self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Mapping Stats (Bases)")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.createBasesStatsTable(ident=ident,lenCell=lenCell))
        vectorSphinx.append("\n")

        self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Bisulfite Conversion Rate")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.createBisulfiteConversionRate(ident=ident,lenCell=lenCell))
        vectorSphinx.append("\n")

        if self.is_paired:
            self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Correct Pairs") 
            vectorSphinx.append("\n")
            vectorSphinx.extend(self.createSingleValueTable(ident=ident,lenCell=lenCell,name='Correct Pairs',value=self.mapping_stats.correct_pairs))
            vectorSphinx.append("\n")
          
        self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Mapping Quality Histogram")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.buildMapqHistogram(ident=ident))
        vectorSphinx.append("\n")
 
        self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Read Length")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.buildReadsLengthTable(ident=ident,lenCell=lenCell))
        vectorSphinx.append("\n")
        
        if self.is_paired:
            self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Insert Size Histogram")
            vectorSphinx.append("\n")
            vectorSphinx.extend(self.buildInsertSizePlot(ident=ident,lenCell=lenCell))
            vectorSphinx.append("\n")
                       
    def addHtmlStartReport(self,ident=20,vectorSphinx=None):
        ''' Add HTML Report Header to a given vector'''
        self.addChapterLine(ident=ident,vectorSphinx=vectorSphinx,title="SAMPLE %s LANE %s"%(self.sample_name,self.mapping_stats.name))
               
    def buildMapqHistogram(self,ident=20):
        """ From matplot lib plots a Mappping qualty histogram
        
        """
        #1. PLOT BUILDING
        self.drawMapqHistogram()
                
        #2. SPHINX TAG BUILDING 
        vectorSphinx = []
   
        imageLine = ".. image:: %s\n" %(os.path.basename(self.png_mapq_histogram))
        vectorSphinx.append(imageLine)

        return vectorSphinx

    def buildReadsLengthTable(self,ident=20,lenCell=30):
        """ Build Read Length Table

            ident - Characters to ident
            lenCell - Length of table cells in characters
            returns vector of Sphinx tags
        """
        vectorSphinx = []
        #1 Table Definition
        self.addSimpleLine(ident=ident,vectorSphinx=vectorSphinx,cells=2,lenCell=lenCell)        
 
        #1.Table header
        cellContents = ""
        cellContents += self.addCell(text="Read Length",lenCell=lenCell,isLeft=True,isRight=False)
        cellContents += self.addCell(text="Reads",lenCell=lenCell,isLeft=False,isRight=True)
        vectorSphinx.append(cellContents)

        self.addHeaderLine(ident=ident,vectorSphinx=vectorSphinx,cells=2,lenCell=lenCell)        
 
        #2.Contents
        for dictValue in self.mapping_stats.read_length_histogram:
            for readLength, reads in dictValue.iteritems():
                self.addRowSimpleValue(ident=ident,lenCell=lenCell,vectorSphinx=vectorSphinx,name_concept=readLength,value=reads)
    
        return vectorSphinx
        
    def buildInsertSizePlot(self,ident=20,lenCell=30):
        """Build Insert Size plot and locate on html report

            ident - Characters to ident
            lenCell - Length of table cells in characters
            returns vector of Sphinx tags
        """
        vectorSphinx = []

        #1. DRAW INSERT SIZE PLOT
        self.drawInsertSizePlot()
        
        #2. HTML TAG BUILDING 
        imageLine = ".. image:: %s\n" %(os.path.basename(self.png_insert_size_histogram))
        vectorSphinx.append(imageLine)

        return vectorSphinx        
        
 
class SampleSphinx(BasicSphinx):
    """ Class to manage Sample report, builds a Sphinx report per sample and call for Lane html report building """

    def __init__(self,project_name=None,sample_stats=None,sphinx_sample=None,png_insert_size_histogram=None,png_mapq_histogram=None):
        """ SampleHtml constructor
        
            project_name - Project name
            sample_stats - SampleStats object from which create html report
            sphinx_sample - HTML sample document name
            png_insert_size_histogram - File png histogram for insert size histogram plot
            png_mapq_histogram - PNG file to store mapq histogram
        """
        #Call Parent class constructor
        BasicSphinx.__init__(self,mapping_stats=sample_stats,png_mapq_histogram=png_mapq_histogram,png_insert_size_histogram=png_insert_size_histogram)        
        
        self.project_name = project_name
        self.is_paired = sample_stats.is_paired
        self.sphinx_sample = sphinx_sample
        
    def run(self,vectorSphinx=None,ident=20,lenCell=30):
        """ Run Sample Sphinx Documentation Building
        
            ident - Characters to ident
            lenCell - Length of table cells in characters
            vectorSphinx - vector of Sphinx tags
            return vectorDocuments to create Toc Tree Index File
        """    
        #0. Add Chapter Line       
        self.addTopSection(ident=ident,vectorSphinx=vectorSphinx,title="Sample %s" %(self.mapping_stats.name))        
        vectorSphinx.append("\n")
        
        #1. Add Sphinx Report Header
        self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Mapping Stats (Reads)")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.createStatsTable(ident=ident,lenCell=lenCell))
        vectorSphinx.append("\n")
        
        #2. Uniqueness Fragments                
        self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Uniqueness (Fragments)")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.createUniqueFragmentsTable(ident=ident,lenCell=lenCell,unique_fragments=self.mapping_stats.totalSampleUniqueReads,
                                                          average_unique=self.mapping_stats.averageSampleUniqueReads))
        vectorSphinx.append("\n")                                                          
                                
        #3. Mapping Stats                                                                                                                                                   
        self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Mapping Stats (Bases)")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.createBasesStatsTable(ident=ident,lenCell=lenCell))
        vectorSphinx.append("\n")
        
        #4. Mapping Quality Histogram
        self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Mapping Quality Histogram")
        vectorSphinx.append("\n") 
        vectorSphinx.extend(self.buildMultipleMapqHistogram(ident=ident))
        vectorSphinx.append("\n")
        
        #6. Bisulfite Conversion Rate        
        self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Bisulfite Conversion Rate")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.createBisulfiteConversionRate(ident=ident,lenCell=lenCell))
        vectorSphinx.append("\n")
                            
        #7. Insert Size Plot   
        if self.is_paired:
            self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Insert Size Histogram")
            vectorSphinx.append("\n")
            vectorSphinx.extend(self.buildMultipleInsertSizePlot(ident=ident)) 
            vectorSphinx.append("\n")
        
        #8. Correct Pairs
        if self.is_paired:
            self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Correct Pairs") 
            vectorSphinx.append("\n")
            vectorSphinx.extend(self.createSingleValueTable(ident=ident,lenCell=lenCell,name='Correct Pairs',value=self.mapping_stats.correct_pairs))
            vectorSphinx.append("\n")
               
        #9. Lane Stats and links per lane
        for laneStats in self.mapping_stats.list_lane_stats:
            mapqHistogram = "%s/%s_mapq.png" %(os.path.dirname(self.sphinx_sample),laneStats.name)
            isizeHistogram = "%s/%s_isize.png" %(os.path.dirname(self.sphinx_sample),laneStats.name)
            laneSphinx = LaneSphinx(project_name=self.project_name,sample_name=self.mapping_stats.name,lane_stats=laneStats,png_mapq_histogram=mapqHistogram,png_insert_size_histogram=isizeHistogram)
            laneSphinx.run(vectorSphinx,ident=0,lenCell=45)
            
                                             
    def addSphinxStartReport(self,vectorSphinx=None,ident=20):
        ''' Add SPHINX Report Header to a given vector'''      
        self.addChapterLine(ident=ident,vectorSphinx=vectorSphinx,title="SAMPLE %s" %(self.mapping_stats.name))
   
    def buildMultipleMapqHistogram(self,ident=20):
        """ From matplot lib plots a Mappping qualty histogram
        
            ident - Characters to ident
        """
        #1. PLOT BUILDING
        self.drawMultipleMapqHistogram()
                
        #2. HTML TAG BUILDING 
        vTableSphinx = []
        imageLine = ".. image:: %s\n" %(os.path.basename(self.png_mapq_histogram))
        vTableSphinx.append(imageLine)   
        
        return vTableSphinx     
        
    def buildMultipleInsertSizePlot(self,ident=20):
        """Build Multiple Insert Size plots and locate on html report
        
          ident - Characters to ident
        """
        #1. PLOT BUILDING
        self.drawMultipleInsertSizePlot()
        
        #2. HTML TAG BUILDING 
        vectorSphinx = []
        imageLine = ".. image:: %s\n" %(os.path.basename(self.png_insert_size_histogram))
        vectorSphinx.append(imageLine)

        return vectorSphinx           
   
class SumupSphinx(BasicSphinx):
    """ Class which defines Index Sample Report """
    
    def __init__(self,output_dir=None,name_project=None,vector_samples=None):
        """  Class constructor
        
             output_dir -- Output directory to store HTML reports
             name_project -- Project name
             vector_sample -- Vector of samples
        """
        #Call Parent class constructor
        BasicSphinx.__init__(self)
        
        self.output_dir = output_dir
        self.name_project = name_project
        self.project_sphinx_document = "%s/%s.rst" %(self.output_dir,self.name_project)
        self.vector_samples = vector_samples
                
    def run(self,vectorSphinx=None):
        """ Run Sumup Sphinx Documentation Building
        
            vectorSphinx - vector of html tags
        """ 
        #1. Add SPHINX Report Header
        self.addPartLine(ident=0,vectorSphinx=vectorSphinx,title="Methylation Pipeline Report Project %s"%(self.name_project))

        vectorDocuments = []
        vectorDocuments.append("SUMUP")        
          
        #2. Create Sample and lanes Report 
        for sampleStats in self.vector_samples:
            sampleSphinx = "%s/%s.rst" %(self.output_dir,sampleStats.name)
            vectorDocuments.append("%s" %(sampleStats.name))
            isizeHistogram = "%s/%s_isize.png" %(os.path.dirname(self.output_dir),sampleStats.name)
            png_mapq_histogram = "%s/%s_mapq.png" %(os.path.dirname(self.output_dir),sampleStats.name)
            sampleReport = SampleSphinx(project_name=self.name_project,sample_stats=sampleStats,sphinx_sample=sampleSphinx,\
                                        png_insert_size_histogram=isizeHistogram,png_mapq_histogram=png_mapq_histogram)
            vSampleSphinx = []
            sampleReport.run(vSampleSphinx,ident=0,lenCell=45)
            RunBasicStats.saveDocument(file_name=sampleSphinx,vectorContent=vSampleSphinx)

        #3. Sumup Table  
        sumupTable = "%s/SUMUP.rst" %(self.output_dir)  
        vectorSummary = []
        #2.1 Create Summary Table 
        self.addTopSection(ident=0,vectorSphinx=vectorSummary,title="Sample Sumup")
        vectorSummary.append("\n")
        self.createSumupSampleTables(ident=0,lenCell=30,vector_sphinx=vectorSummary,samplesStats=self.vector_samples)
        #2.2 Save Sumary Table 
        RunBasicStats.saveDocument(file_name=sumupTable,vectorContent=vectorSummary)

         
        #3.Toc Table
        vectorSphinx.append("\n")
        self.addTocTree(vectorSphinx=vectorSphinx,vectorDocuments=vectorDocuments)
        
class ConfigSphinx(object):
    """Sphinx configuration python file generation"""
    def __init__(self,path_config_file=None,path_makefile_file=None,master_file=None,project_name=None,main_title='BISULFITE PIPELINE'):
        """Constructor Class for config file generation
        
            path_config_file - Configuration python file for sphinx compiling
            master_file - Master file where it is found main toc tree index
        """
        self.config_file = path_config_file
        self.master_file = master_file
        self.project_name = project_name
        self.path_makefile_file = path_makefile_file
        
        #Scaped name
        self.main_title = ""
        for character in main_title:
            if character == "_":
                self.main_title += "\\%s" %(character)
            else:
                self.main_title += "%s" %(character)  
                
        

        main_title
        
    def run(self):
        """Creates config.py file"""
        vectorConfig = []
        vectorConfig.append("import sys, os")
        vectorConfig.append("extensions = ['sphinx.ext.autodoc', 'sphinx.ext.doctest', 'sphinx.ext.mathjax']")
        vectorConfig.append("templates_path = ['_templates']")
        vectorConfig.append("source_suffix = '.rst'")
        vectorConfig.append("master_doc = '%s'" %(self.master_file))
        vectorConfig.append("project = u'Bisulfite Pipeline'")
        vectorConfig.append("exclude_patterns = []")
        vectorConfig.append("pygments_style = 'sphinx'")
        vectorConfig.append("import sphinx_rtd_theme")
        vectorConfig.append('html_theme = "sphinx_rtd_theme"')
        vectorConfig.append("html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]")
        vectorConfig.append("html_static_path = ['_static']")
        vectorConfig.append("htmlhelp_basename = 'BISULFITEPIPELINEdoc'")
        vectorConfig.append("#Latex Elements")
        
        vectorConfig.append("latex_elements = {")
        vectorConfig.append("# The paper size ('letterpaper' or 'a4paper').")
        vectorConfig.append("#'papersize': 'letterpaper',")
        vectorConfig.append("# The font size ('10pt', '11pt' or '12pt').")
        vectorConfig.append("'pointsize': '10pt',")
        vectorConfig.append("# Additional stuff for the LaTeX preamble.")
        vectorConfig.append("#'preamble': '',")
        vectorConfig.append("'classoptions': ',oneside',")
        vectorConfig.append("'babel': '\\usepackage[english]{babel}',")

        #Scaped name
        scaped_name = ""
        for character in self.project_name:
            if character == "_":
                scaped_name += "\\%s" %(character)
            else:
                scaped_name += "%s" %(character)         
        
        
        vectorConfig.append("'releasename':'%s'" %(scaped_name))
        vectorConfig.append("}")
        
        vectorConfig.append("latex_documents = [")
        vectorConfig.append("('%s', '%s.bisulfite.tex', u'%s'," %(self.master_file,self.master_file,self.main_title))
        vectorConfig.append("u'gemBS', 'manual'),")
        vectorConfig.append("]")

        with open(self.config_file, 'w') as fileDocument:
            for line in vectorConfig:
                fileDocument.write("%s\n" %(line))
    
        self.runMakefile()
        
    def runMakefile(self):
        """Creates makefile sphinx configuration"""
        vectorConfig = []
        
        vectorConfig.append("SPHINXOPTS    =")
        vectorConfig.append("SPHINXBUILD   = sphinx-build")
        vectorConfig.append("PAPER         =")
        vectorConfig.append("BUILDDIR      = build")        
        
        vectorConfig.append("PAPEROPT_a4     = -D latex_paper_size=a4")
        vectorConfig.append("PAPEROPT_letter = -D latex_paper_size=letter")
        vectorConfig.append("ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .")
        vectorConfig.append("I18NSPHINXOPTS  = $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .")
        vectorConfig.append("\n")
        vectorConfig.append(".PHONY: help clean html dirhtml singlehtml htmlhelp latex latexpdf")
        vectorConfig.append("\n")
        vectorConfig.append("help:")
        vectorConfig.append('\t@echo \"Please use make <target> where <target> is one of"')
        vectorConfig.append('\t@echo "  html       to make standalone HTML files"')
        vectorConfig.append('\t@echo "  dirhtml    to make HTML files named index.html in directories"')
        vectorConfig.append('\t@echo "  singlehtml to make a single large HTML file"')
        vectorConfig.append('\t@echo "  htmlhelp   to make HTML files and a HTML help project"')
        vectorConfig.append('\t@echo "  latex      to make LaTeX files, you can set PAPER=a4 or PAPER=letter"')
        vectorConfig.append('\t@echo "  latexpdf   to make LaTeX files and run them through pdflatex"')
        vectorConfig.append("\n")
   
        vectorConfig.append('clean:')
        vectorConfig.append('\t-rm -rf $(BUILDDIR)/*')
        vectorConfig.append("\n")

        vectorConfig.append('html:')
        vectorConfig.append('\t$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html')
        vectorConfig.append('\t@echo')
        vectorConfig.append('\t@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."')
        vectorConfig.append("\n")

        vectorConfig.append('dirhtml:')
        vectorConfig.append('\t$(SPHINXBUILD) -b dirhtml $(ALLSPHINXOPTS) $(BUILDDIR)/dirhtml')
        vectorConfig.append('\t@echo')
        vectorConfig.append('\t@echo "Build finished. The HTML pages are in $(BUILDDIR)/dirhtml."')
        vectorConfig.append("\n")

        vectorConfig.append('singlehtml:')
        vectorConfig.append('\t$(SPHINXBUILD) -b singlehtml $(ALLSPHINXOPTS) $(BUILDDIR)/singlehtml')
        vectorConfig.append('\t@echo')
        vectorConfig.append('\t@echo "Build finished. The HTML page is in $(BUILDDIR)/singlehtml."')
        vectorConfig.append("\n")
 
        vectorConfig.append('latex:')
        vectorConfig.append('\t$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex')
        vectorConfig.append('\t@echo')
        vectorConfig.append('\t@echo "Build finished; the LaTeX files are in $(BUILDDIR)/latex."')
        vectorConfig.append('\t@echo "Run make in that directory to run these through (pdf)latex"')
        vectorConfig.append('\t@echo "(use make latexpdf here to do that automatically)."')
        vectorConfig.append("\n")

        vectorConfig.append('latexpdf:')
        vectorConfig.append('\t$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex')
        vectorConfig.append('\t@echo "Running LaTeX files through pdflatex..."')
        vectorConfig.append('\t$(MAKE) -C $(BUILDDIR)/latex all-pdf')
        vectorConfig.append('\t@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."')
        vectorConfig.append("\n")

        with open(self.path_makefile_file, 'w') as fileDocument:
            for line in vectorConfig:
                fileDocument.write("%s\n" %(line))

           
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

    #SumupSphinx object
    vector_sumup_sphinx = []
    sumupSphinx = SumupSphinx(output_dir=output_dir,name_project=name,vector_samples=vector_samples)
    sumupSphinx.run(vectorSphinx=vector_sumup_sphinx)
    RunBasicStats.saveDocument(file_name=sumupSphinx.project_sphinx_document,vectorContent=vector_sumup_sphinx)
    
    #Config python file
    cfgFile = ConfigSphinx(path_config_file="%s/conf.py" %(output_dir),path_makefile_file="%s/Makefile" %(output_dir),master_file=name,project_name=name)
    cfgFile.run()
    