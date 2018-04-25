# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 15:35:46 2017

@author: marcos
"""

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

import numpy as np
import os
import math


class StatsMother(object):
    """Common methods to all type of class stats"""
    
    def __init__(self):
        """ Defines a Dictionary of Values"""
        self.data = {}
    
    def sumValues(self,keyDict,dictionary):
        """ Sums value in conter if the  key is present in a given dictionary
        
            counter - member to add the new value
            keyDict - Key concept
            dictionary  - Dictionary of concept and values
        """
        if keyDict in dictionary:
            if keyDict in self.data:
                self.data[keyDict] += int(dictionary[keyDict])
            else:
                self.data[keyDict] = int(dictionary[keyDict])
            
    def getTotalSumatory(self,vector):
        """ Gets the total number from the sum of all elements of a vector 
        
            returns -- Total summatory
        """
        total = 0
        for i in vector:
            total += i
        return total
    
    def getPercentage(self,subtotal,total):
        """Calculates the percentage for a given value"""
        if total != 0.0:
            return (float(subtotal)/float(total))*100
        else:
            return 0.0


class PlotMother(object):
    """Common methods for plotting"""
    
    def __init__(self,pngFile,concept,relativeImgDirectory="./IMG/",relativeSphinxPathImage="../IMG/"):
        """Plot Mother Constructor"""
        self.pngFile = pngFile
        self.relativePathImage = "%s/%s" %(relativeImgDirectory,os.path.basename(pngFile))
        self.relativeSphinxPathImage = "%s/%s" %(relativeSphinxPathImage,os.path.basename(pngFile))
        self.haveCounts = False
        self.concept = concept
        self.ax = None
        self.figure = None
        
    def setHaveCounts(self,value):
        """Set have counts Status"""
        self.haveCounts = value
        
    def getHaveCounts(self):
        """Get Have Counts Status"""
        return self.haveCounts
        
    def initSubplots(self,figure_size=(5.7,5),yText="",xText=""):
        """Initializate subplots"""
        self.figure, self.ax = plt.subplots(1,1,figsize=figure_size)           
        self.ax.set_ylabel(yText)
        self.ax.set_xlabel(xText)
        self.ax.set_title(self.concept)        
                
    def setYAxisExponential(self):
        """ Force Y axis to be exponential """
        #self.ax.ticklabel_format(style='sci', axis='y')
        #self.ax.yaxis.major.formatter.set_powerlimits((-1,1))        
        self.move_sn_y(offs=-.05,dig=1,side='left')
                
    def move_sn_y(self,offs=0, dig=0, side='left', omit_last=False):
        """Method copied from: https://werthmuller.org/blog/2014/move-scientific-notation/
           Move scientific notation exponent from top to the side.
           Additionally, one can set the number of digits after the comma for the y-ticks, hence if it should state 1, 1.0, 1.00 and so forth.

           offs -      float, optional; <0> Horizontal movement additional to default.
           dig  -      int, optional; <0> Number of decimals after the comma.
           side -      string, optional; {<'left'>, 'right'} To choose the side of the y-axis notation.
           omit_last - bool, optional; <False> If True, the top y-axis-label is omitted.
     
           return locs - list List of y-tick locations.
        """
        #Get the ticks
        locs, _ = plt.yticks()
        #Put the last entry into a string, ensuring it is in scientific notation
        llocs = '%.3e' % locs[-1]
        #Get the magnitude, hence the number after the 'e'. E.g: '1.235e+08' => 8
        yoff = int(str(llocs).split('e')[1])
        #If omit_last, remove last entry
        if omit_last:
            slocs = locs[:-1]
        else:
            slocs = locs
        #Set ticks to the requested precision
        form = r'$%.'+str(dig)+'f$'
        plt.yticks(locs, list(map(lambda x: form % x, slocs/(10**yoff))))
        #Define offset depending on the side
        if side == 'left':
            offs = -.18 - offs # Default left: -0.18
        elif side == 'right':
            offs = 1 + offs    # Default right: 1.0
        #Plot the exponent
        self.ax.text(offs, .98, r'$\times10^{%i}$' % yoff,transform=plt.gca().transAxes, verticalalignment='top')
        #Return the locs
        return locs        

    def saveAndClose(self):
        """Save the figure and close it"""
        pylab.savefig(self.pngFile)
        plt.close(self.figure)

    def add(self):
        """Virtual Method to be defined in Child Class"""
        pass        
        
    def plot(self):
        """Virtual Method to be Defined in Child Class"""
        pass

class SimpleStats(StatsMother):
    """ Simple stats Variants """
    
    def __init__(self,concept=None):
        """ Initializates basic members 
       
            concept - Name of the concept        
        """
        StatsMother.__init__(self)
        self.concept = concept
        self.data["All"] = 0
        self.data["Passed"] = 0
        self.has_counts = False
        
    def add(self,values={}):
        """
            Add (Sum) new Values
            
            values - Dictionary of values for Simple Stats
        """
        self.sumValues("All",values)
        self.sumValues("Passed",values)
        self.has_counts = True
        
    def hasCounts(self):
        """
            Returns True If the Simple Stats must be reported
        """
        return self.has_counts


    def getList(self):
        """ 
            Get List to be printed in a row
            
            return list to be printed in a row table
        """
        percentage = self.getPercentage(self.data["Passed"],self.data["All"])
        return [self.concept,"%i"%(self.data["All"]),"%i"%(self.data["Passed"]),"%.2f %%"%(percentage)]            
        
class TotalStats(object):
    """Manage Total Stats Table"""
    
    def __init__(self):
        """ Initialize Total Stats Counter"""
        self.snps = SimpleStats("SNPs")
        self.indels = SimpleStats("Indels")
        self.multiallelic = SimpleStats("Multiallelic")
        self.dbSNPsites = SimpleStats("dbSNPsites")
        self.dbSNPVariantSites = SimpleStats("dbSNPVariantSites")
        self.refCpGs = SimpleStats("RefCpG")
        self.nonRefCpGs = SimpleStats("NonRefCpG")
        
    def addSimpleStats(self,simpleStat,concept,value):
        """ Add count simple stats for a given concept and value
        
            simpleStat -- simple stat to be processed
            concept -- key to be extracted from values
            value -- dictionary to be checked
        """
        if concept in value:
            simpleStat.add(value[concept])
    
    def add(self,values):
        """ Add dictionary of values to the total Stats counters
        
            values -- Dictionary of values
        """
        self.addSimpleStats(self.snps,"SNPS",values)
        self.addSimpleStats(self.indels,"Indels",values)
        self.addSimpleStats(self.multiallelic,"Multiallelic",values)
        self.addSimpleStats(self.dbSNPsites,"dbSNPsites",values)
        self.addSimpleStats(self.dbSNPVariantSites,"dbSNPVariantSites",values)
        self.addSimpleStats(self.refCpGs,"RefCpG",values)
        self.addSimpleStats(self.nonRefCpGs,"NonRefCpG",values)
    
    def setListTable(self,sampleStats,vector):
        """Add sampleStats to the table"""
        if sampleStats.hasCounts():
            vector.append(sampleStats.getList())
    
    def getTableVariants(self):
        """Return a Table for the Variants Report"""
        table_vector = [["Type","Total","Pass","%"]]
        #1.SNPs
        self.setListTable(self.snps,table_vector)
        #2.Multiallelic
        self.setListTable(self.multiallelic,table_vector)
        #3.dbSNPsites
        self.setListTable(self.dbSNPsites,table_vector)
        #4.dbSNPVariantSites
        self.setListTable(self.dbSNPVariantSites,table_vector)
        return table_vector
        
    def getTableMethylation(self):
        """Return a table for the Methylation Report"""
        table_vector = [["Type","Total","Pass","%"]]
        #1.Reference CpGs
        self.setListTable(self.refCpGs,table_vector)
        #2.Non Reference CpGs
        self.setListTable(self.nonRefCpGs,table_vector)
        return table_vector
        
class ReadLevelStats(StatsMother):
    """ Read level Stats"""
    
    def __init__(self,concept=None):
        """ Initializes Read Level Stats
        
            concept - Name of the concept
        """
        StatsMother.__init__(self)
        self.concept = concept
        self.data["Reads"] = 0
        self.data["Bases"] = 0
        
    def add(self,values={}):
        """
            Add (Sum) new Values
            
            values - Dictionary of values for Read Level Stats
        """
        self.sumValues("Reads",values)
        self.sumValues("Bases",values)
        
class ReadsAndBases(StatsMother):
    """ Manages The set of statistic grouped by ReadLevel Stats"""
    
    def __init__(self):
        """
            Class constructor to Initialize ReadLevelStats memebers
        """
        StatsMother.__init__(self)
        self.passed = ReadLevelStats("Passed")
        self.unmapped = ReadLevelStats("Unmapped")        
        self.mateUnmapped = ReadLevelStats("MateUnmapped")
        self.duplicate = ReadLevelStats("Duplicate")
        self.badOrientation = ReadLevelStats("BadOrientation")
        self.largeInsertSize = ReadLevelStats("LargeInsertSize")        
        self.lowMAPQ = ReadLevelStats("LowMAPQ")
        self.notCorrectlyAligned = ReadLevelStats("NotCorrectlyAligned")
    
    def addValues(self,feature=None,featureName="",values={}):
        """
            Call read level stats checking the exitance in the dictionary
            
            feature - Feature instance to add values
            featureName - Concept of the feature
            values - Dictionary of values
        """
        if featureName in values:
            feature.add(values[featureName])
    
    def add(self,values={}):
        """
            Updates ReadLevel Registers From a dictionary
        """
        self.addValues(feature=self.passed,featureName="Passed",values=values)            
        self.addValues(feature=self.unmapped,featureName="Unmapped",values=values)            
        self.addValues(feature=self.mateUnmapped,featureName="MateUnmapped",values=values)            
        self.addValues(feature=self.duplicate,featureName="Duplicate",values=values)            
        self.addValues(feature=self.badOrientation,featureName="BadOrientation",values=values)            
        self.addValues(feature=self.largeInsertSize,featureName="LargeInsertSize",values=values)
        self.addValues(feature=self.lowMAPQ,featureName="LowMAPQ",values=values)
        self.addValues(feature=self.notCorrectlyAligned,featureName="NotCorrectlyAligned",values=values)

    def getTotalBases(self):
        """
            Get total number of bases
        """
        total_mapped_bases = self.getTotalSumatory([self.passed.data["Bases"],self.duplicate.data["Bases"],self.badOrientation.data["Bases"],
                                                    self.largeInsertSize.data["Bases"],self.lowMAPQ.data["Bases"],self.notCorrectlyAligned.data["Bases"]])
        total_bases = total_mapped_bases + self.unmapped.data["Bases"]
        return total_bases
    
    def getTable(self):        
        """ Get vector of rows to be printed in a table report"""
        vector_table = []
        vector_table.append(["Type","#reads","%","#bases","%"])
        #1.Get Total Mapped Reads
        total_mapped_reads = self.getTotalSumatory([self.passed.data["Reads"],self.duplicate.data["Reads"],self.badOrientation.data["Reads"],
                                                    self.largeInsertSize.data["Reads"],self.lowMAPQ.data["Reads"],self.notCorrectlyAligned.data["Reads"]])
        total_reads = total_mapped_reads + self.unmapped.data["Reads"]
        #2.Get Total Mapped Bases
        total_mapped_bases = self.getTotalSumatory([self.passed.data["Bases"],self.duplicate.data["Bases"],self.badOrientation.data["Bases"],
                                                    self.largeInsertSize.data["Bases"],self.lowMAPQ.data["Bases"],self.notCorrectlyAligned.data["Bases"]])
        total_bases = total_mapped_bases + self.unmapped.data["Bases"]
                
        #3.Row Total
        vector_table.append(["Total","%i"%(total_reads),"100%","%i"%(total_bases),"100%"])
        #4.Row Mapped
        vector_table.append(["  Mapped","%i"%(total_mapped_reads),"%.2f %%"%(self.getPercentage(total_mapped_reads,total_reads)),
                                        "%i"%(total_mapped_bases),"%.2f %%"%(self.getPercentage(total_mapped_bases,total_bases))])
        #5.Row Passed
        vector_table.append(["  Passed","%i"%(self.passed.data["Reads"]),"%.2f %%"%(self.getPercentage(self.passed.data["Reads"],total_reads)),
                                        "%i"%(self.passed.data["Bases"]),"%.2f %%"%(self.getPercentage(self.passed.data["Bases"],total_bases))])
        #6.Row Not Passed
        total_not_passed_reads = total_mapped_reads - self.passed.data["Reads"]
        total_not_passed_bases = total_mapped_bases - self.passed.data["Bases"]
        vector_table.append(["  NotPassed","%i"%(total_not_passed_reads),"%.2f %%"%(self.getPercentage(total_not_passed_reads,total_reads)),
                                           "%i"%(total_not_passed_bases),"%.2f %%"%(self.getPercentage(total_not_passed_bases,total_bases))])
        #7.Empty Row
        vector_table.append(["","","","",""])
        #8.LowMAPQ Row
        vector_table.append(["    LowMAPQ","%i"%(self.lowMAPQ.data["Reads"]),"%.2f %%"%(self.getPercentage(self.lowMAPQ.data["Reads"],total_reads)),
                                           "%i"%(self.lowMAPQ.data["Bases"]),"%.2f %%"%(self.getPercentage(self.lowMAPQ.data["Bases"],total_bases))])
        #9.MateUnmapped Row
        vector_table.append(["    MateUnmapped","%i"%(self.mateUnmapped.data["Reads"]),"%.2f %%"%(self.getPercentage(self.mateUnmapped.data["Reads"],total_reads)),
                                           "%i"%(self.mateUnmapped.data["Bases"]),"%.2f %%"%(self.getPercentage(self.mateUnmapped.data["Bases"],total_bases))])
        #10.Duplicate Row
        vector_table.append(["    Duplicate","%i"%(self.duplicate.data["Reads"]),"%.2f %%"%(self.getPercentage(self.duplicate.data["Reads"],total_reads)),
                                             "%i"%(self.duplicate.data["Bases"]),"%.2f %%"%(self.getPercentage(self.duplicate.data["Bases"],total_bases))])
        #11.BadOrientation Row
        vector_table.append(["    BadOrientation","%i"%(self.badOrientation.data["Reads"]),"%.2f %%"%(self.getPercentage(self.badOrientation.data["Reads"],total_reads)),
                                             "%i"%(self.badOrientation.data["Bases"]),"%.2f %%"%(self.getPercentage(self.badOrientation.data["Bases"],total_bases))])                
        #12.LargeInsertSize Row
        vector_table.append(["    LargeInsertSize","%i"%(self.largeInsertSize.data["Reads"]),"%.2f %%"%(self.getPercentage(self.largeInsertSize.data["Reads"],total_reads)),
                                             "%i"%(self.largeInsertSize.data["Bases"]),"%.2f %%"%(self.getPercentage(self.largeInsertSize.data["Bases"],total_bases))])   
                                             
        #13.Empty Row
        vector_table.append(["","","","",""])

        #14.Not Correctly Aligned Row
        vector_table.append(["    NotCorrectlyAligned","%i"%(self.notCorrectlyAligned.data["Reads"]),"%.2f %%"%(self.getPercentage(self.notCorrectlyAligned.data["Reads"],total_reads)),
                                             "%i"%(self.notCorrectlyAligned.data["Bases"]),"%.2f %%"%(self.getPercentage(self.notCorrectlyAligned.data["Bases"],total_bases))])
        #15.Unmapped Row
        vector_table.append(["    Unmapped","%i"%(self.unmapped.data["Reads"]),"%.2f %%"%(self.getPercentage(self.unmapped.data["Reads"],total_reads)),
                                             "%i"%(self.unmapped.data["Bases"]),"%.2f %%"%(self.getPercentage(self.unmapped.data["Bases"],total_bases))])  
        return vector_table
        
    
class BaseLevel(StatsMother):
    """Base Level Statistics"""
    
    def __init__(self):
        """ Initializes BaseLevel class members """
        StatsMother.__init__(self)
        self.data["Passed"] = 0
        self.data["Trimmed"] = 0
        self.data["Clipped"] = 0
        self.data["Overlapping"] = 0
        self.data["LowQuality"] = 0
        
    def add(self, values={}):
        """
            Add (Sum) new Values
            
            values - Dictionary of values for Base Level Stats
        """   
        self.sumValues("Passed",values)
        self.sumValues("Trimmed",values)
        self.sumValues("Clipped",values)
        self.sumValues("Overlapping",values)
        self.sumValues("LowQuality",values)
            
    def getTable(self):
        """ Get vector of rows to be printed in a table report"""
        vector_table = []
        vector_table.append(["Bases","#","%"])
        total_bases = self.getTotalSumatory([self.data["Passed"],self.data["Trimmed"],self.data["Clipped"],self.data["Overlapping"],self.data["LowQuality"]])
        vector_table.append(["Total",total_bases,"100%"])
        vector_table.append(["Passed","%i" %(self.data["Passed"]),"%.2f %%" %(self.getPercentage(self.data["Passed"],total_bases)) ])
        vector_table.append(["Trimmed","%i" %(self.data["Trimmed"]),"%.2f %%" %(self.getPercentage(self.data["Trimmed"],total_bases))])
        vector_table.append(["Clipped","%i" %(self.data["Clipped"]),"%.2f %%" %(self.getPercentage(self.data["Clipped"],total_bases))])
        vector_table.append(["Overlapping","%i" %(self.data["Overlapping"]),"%.2f %%" %(self.getPercentage(self.data["Overlapping"],total_bases))])
        vector_table.append(["LowQuality","%i" %(self.data["LowQuality"]),"%.2f %%" %(self.getPercentage(self.data["LowQuality"],total_bases))])
        return vector_table
        
class DistributionPlot(PlotMother,StatsMother):
    """Manages Distribution Plots"""
    
    def __init__(self,concept,pngFile):
        """
            Initializes coverage class        
        
            concept - Coverage concept name
            pngFile - Png File name for the coverge plot 
        """
        PlotMother.__init__(self,pngFile=pngFile,concept=concept)
        self.percentage_limit_tail = 95
        self.value_locations = {}
                
    def getVectorToPlot(self,cleanTail=False):
        """From Dictionary Get Vector of Y values to be bar plotted 
        
           cleanTail - if true cleans all registers that does not represents the 95% of the data        
        """
        vector_plot = []
        
        maximumX = 0
         
        if cleanTail:
             #1. Get Total Number of Sites
             totalSites = 0
             for value in self.value_locations:
                 totalSites += self.value_locations[value]
             #2.Get Maximum X 95%
             currentSites = 0
             for value in self.value_locations:
                 currentSites += self.value_locations[value]
                 maximumX = value
                                                  
                 if self.getPercentage(currentSites,totalSites) > self.percentage_limit_tail:
                     if maximumX == 0:
                         maximumX = 1
                     break
        else:
            list_keys = self.value_locations.keys()
            list_keys.sort()        
            maximumX = list_keys[-1]
                   
        for x in range(0,(maximumX+1)):
            if x in self.value_locations:
                vector_plot.append(self.value_locations[x])
            else:
                vector_plot.append(0)
        return vector_plot
        
class Coverage(DistributionPlot):
    """ Class Coverage """

    def __init__(self,concept=None,yLabel="#SNPs",pngFile=None):
        """
            Initializes coverage class        
        
            concept - Coverage concept name
            pngFile - Png File name for the coverge plot 
            yLabel - Y label for the plot
        """
        
        DistributionPlot.__init__(self,pngFile=pngFile,concept=concept)
        self.yLabel = yLabel
        self.percentage_limit_tail = 98
        
    def add(self,values):
        """
            Add dictionary of new coverage values
            
            values - Dictionary of coverage and bases
        """
        for key in values:
            coverage = int(key)
            if coverage in self.value_locations:
                self.value_locations[coverage] += values[key]
            else:
                self.value_locations[coverage] = values[key]
        
        self.setHaveCounts(True)        
                
    def plot(self):
        """Builds a Coverage Plot and print it to a file"""
        #Number of Bins
        widhtBin = 1
        
        #Create vector of Y values
        yValues = self.getVectorToPlot(True)
                
        #Create Vector of X ticks
        xValues = range(len(yValues))
        
        #Build Plot           
        self.initSubplots(figure_size=(5.7,5),yText=self.yLabel,xText="Coverage")           
                     
        width = widhtBin     
        self.ax.bar(xValues, yValues, width, color="blue")
        self.ax.set_xlim(0, len(xValues))

        self.setYAxisExponential()
        self.saveAndClose()
        
    def getMean(self):
        """
            Get Mean Value from value_locations
        """
        total_number_of_coverage_events = 0
        total_number_of_bases = 0
        for coverage in self.value_locations:
            total_number_of_coverage_events += coverage * self.value_locations[coverage]
            total_number_of_bases += self.value_locations[coverage]
            
        return float(total_number_of_coverage_events)/float(total_number_of_bases)
        
    def getTotalMinimumCoverage(self, minimum_coverage = 0):
        """
            Get total bases from a minimum coverage
            
            minimum_coverage -- Minimum coverage to count total number of bases
        """
        total_counted_bases = 0
        for coverage in self.value_locations:
            if coverage >= minimum_coverage:
                total_counted_bases += self.value_locations[coverage]
        return total_counted_bases
        
class Quality(DistributionPlot):
    """ Class quality """
    
    def __init__(self,concept=None,yLabel="#SNPs",pngFile=None):
        """
            Initializes quality class        
        
            concept - Quality concept name
            pngFile - Png File name for the coverge plot 
        """
        DistributionPlot.__init__(self,pngFile=pngFile,concept=concept)
        self.quality_sites = []
        self.yLabel = yLabel
        
    def add(self,values):
        """
            Add dictionary of new coverage values
            
            values - Dictionary of coverage and bases
        """
        quality = 0
        for sites in values:
            if quality in self.value_locations:
                self.value_locations[quality] += sites
            else:
                self.value_locations[quality] = sites
            quality = quality + 1
        
        self.setHaveCounts(True)  
        
    def plot(self):
        """Builds a Coverage Plot and print it to a file"""
        #Number of Bins
        bins = 25
        
        #Create vector of Y values
        yValues = [0 for x in range(bins)]
        #Quality Sites
        quality_sites = self.getVectorToPlot(False)
        total_quality = len(quality_sites)        
        
        for quality in range(0,total_quality):
            index_for_quality = int ((float(bins)/total_quality)*(quality)) 
            yValues[index_for_quality] = yValues[index_for_quality] + quality_sites[quality]
            
        #Create Vector of X ticks
        xValues = [(total_quality/bins) * x for x in range(bins)]
        
        #Build Plot           
        self.initSubplots(figure_size=(5.7,5),yText=self.yLabel,xText="Quality")           
                     
        width = int (float(total_quality)/bins)
        self.ax.bar(xValues, yValues, width, color="#624ea7")

        self.setYAxisExponential()
        self.saveAndClose()
                        
class Mutations(StatsMother,PlotMother):
    """ Class to Manage Mutation Rate and Ti/Tv Ratio """
    
    def __init__(self,concept):
        """
            Initialize Mutation Concept Name
            
            concept - Mutations Concept Name
        """
        PlotMother.__init__(self,pngFile="",concept=concept)
        self.mutations_hits = {}
        #Status of mutations
        self.status = ["All","Passed","dbSNPAll","dbSNPPassed"]
        #Mutation types definition
        self.mutations_types = [
                                 ("Transition","A>G"),("Transition","G>A"),("Transition","T>C"),("Transition","C>T"),
                                 ("Transversion","A>C"),("Transversion","C>A"),("Transversion","T>G"),("Transversion","G>T"),
                                 ("Transversion","A>T"),("Transversion","T>A"),("Transversion","C>G"),("Transversion","G>C")
                               ]   
        
    def add(self,values):
        """
            Set of values to be updated at the vector of mutations
        
            values - dictionary of values for a given set of mutation
        """
        for mutation in values:
            if mutation in self.mutations_hits:
                for status in values[mutation]:
                    if status in self.mutations_hits[mutation]:
                        self.mutations_hits[mutation][status] += values[mutation][status]
                    else:
                        self.mutations_hits[mutation][status] = values[mutation][status]
            else:
                self.mutations_hits[mutation] = values[mutation]
        self.setHaveCounts(True)
        
    def getTotalStatus(self,status):
        """ Get total number of mutations for a given Status """
        total_status = 0        
        for mutation in self.mutations_hits:
            if status in self.mutations_hits[mutation]:
                total_status += self.mutations_hits[mutation][status]
        return total_status        
                
    def getTiTv(self,status):
        """ Calculates Transition Transversion Ratio 
        
            status - Status ot the mutation        
            return list of Ts/Tv, Transitions,Transversions        
        """
        transition = 0
        transversion = 0
        
        for mutation in self.mutations_hits:
            if status in self.mutations_hits[mutation]:
                if mutation == "A>G" or mutation == "T>C" or mutation == "C>T" or mutation == "G>A":
                    transition += self.mutations_hits[mutation][status]
                elif mutation == "A>C" or mutation == "T>G" or mutation == "A>T" or mutation == "T>A" or mutation == "C>A" or mutation == "G>T" or mutation == "C>G" or mutation == "G>C" :
                    transversion += self.mutations_hits[mutation][status]
          
        if transversion == 0:
            return ["0","%i" %(transition),"%i" %(transversion)]
        else:
            return ["%.2f" %(float(transition) / float(transversion)),"%i" %(transition),"%i" %(transversion)]
        
    def getTableMutationProfile(self):
        """
            Get a list of vector to build a Mutation Ratio
        """
        vector_table = []
        vector_table.append(["Type","Mutation","Status","#","%"])
        
        totalStatus = {}
        for mutation_status in self.status:
            totalStatus[mutation_status] = self.getTotalStatus(mutation_status)
                
        for mutation_status in self.status:
            for mutation_type in self.mutations_types:
                if mutation_status in self.mutations_hits[mutation_type[1]]:
                    vector_table.append([
                                        mutation_type[0],mutation_type[1],mutation_status,
                                        "%i" %(self.mutations_hits[mutation_type[1]][mutation_status]),
                                        "%.2f %%" %(self.getPercentage(self.mutations_hits[mutation_type[1]][mutation_status],totalStatus[mutation_status]))
                                       ])
            #Separator Row beetween each Mutation Status                           
            vector_table.append(["","","","",""])                              
          
        return vector_table
        
    def getTiTvTable(self):
        """
            Get a table of for the transition transversion Rate
        """
        vector_table = []
        vector_table.append(["Status","Ratio","Transitions","Transversions"])
        
        for mutation_status in self.status:
            ti_tv_fields = self.getTiTv(mutation_status)
            if ti_tv_fields[-2] != 0 and ti_tv_fields[-1] != 0:
                ti_tv_fields.insert(0, mutation_status)
                vector_table.append(ti_tv_fields)
                
        return vector_table

class Methylation(PlotMother):
    """ Class Methylation on CpGs """
    
    def __init__(self,concept="",pngFile=""):
        """
            Initializes quality class        
        
            concept - Quality concept name
            pngFile - Png File name for the coverge plot 
        """
        PlotMother.__init__(self,pngFile=pngFile,concept=concept)
        self.methylation_cpgs = []
        
    def add(self,values):
        """
            Add vector of new quality values
            
            values - Vector of qualities and bases
        """
        if len(self.methylation_cpgs) == 0:
            self.methylation_cpgs = values
        else:
            meth = 0
            for cpgs in values:
                self.methylation_cpgs[meth] += cpgs
                meth += 1
        self.setHaveCounts(True)
        
    def getPercentageCpGsVector(self):
        """
            Modify CpGs vector on PerCentages
        """
        #1.Get total CpGs
        total_cpgs = 0
        for cpgs in self.methylation_cpgs:
            total_cpgs += cpgs
        if total_cpgs == 0:
            return [0.0 for i in range(len(self.methylation_cpgs))]
        #2.Calculate CpGs Percentage
        percentage_vector = []
        for cpgs in self.methylation_cpgs:
            percentage_vector.append((float(cpgs)/total_cpgs)*100)
        return percentage_vector
        
    def getMean(self):
        """
            Get Mean Value from value_locations
        """
        total_number_of_methylation_events = 0
        total_number_of_cpgs = 0
        length_meth_cpgs = len(self.methylation_cpgs)
        for methylation in range(0,length_meth_cpgs):
            total_number_of_methylation_events += methylation * self.methylation_cpgs[methylation]
            total_number_of_cpgs += self.methylation_cpgs[methylation]

        if  total_number_of_cpgs == 0 or length_meth_cpgs == 0:
            return 0.0
            
        return (float(total_number_of_methylation_events)/float(total_number_of_cpgs))/float(length_meth_cpgs)
        
       
        
class PlotMethylationLevels(PlotMother):
    """ Class for plotting four types of methylation levels """
    
    def __init__(self,concept,pngFile,meth_list=[]):
        """
            Initializes plot of methylation levels        
        
            concept - Quality concept name
            pngFile - Png File name for the coverge plot
            meth_list - List of Methylation object to representated at the same plot
        """
        PlotMother.__init__(self,pngFile=pngFile,concept=concept)
        self.allRefCpgsMethylation = meth_list[0]
        self.passedRefCpgsMethylation = meth_list[1]
        self.allNonRefCpgsMethylation = meth_list[2]
        self.passedNonRefCpgsMethylation = meth_list[3]

    def plot(self):
        """Builds a Quality Plot and print it to a file
        
           Based on ideas found at: http://www.southampton.ac.uk/~fangohr/training/python/notebooks/Matplotlib.html
        
        """
        #1. Plot Configuration
        self.initSubplots(figure_size=(5.7,5),yText="% CpGs",xText="% Methylation")
        #2. Plot building
        self.ax.set_title("Methylation Status")
        x_scale = [x for x in range(101)]          
        self.ax.plot(x_scale,self.allRefCpgsMethylation.getPercentageCpGsVector(),color="black",lw=1,ls='-',label=self.allRefCpgsMethylation.concept)
        self.ax.plot(x_scale,self.passedRefCpgsMethylation.getPercentageCpGsVector(),color="red",lw=1,linestyle='-',label=self.passedRefCpgsMethylation.concept)
        self.ax.plot(x_scale,self.allNonRefCpgsMethylation.getPercentageCpGsVector(),color="green",lw=1,ls='-',label=self.allNonRefCpgsMethylation.concept)
        self.ax.plot(x_scale,self.passedNonRefCpgsMethylation.getPercentageCpGsVector(),color="blue",lw=1,linestyle='-',label=self.passedNonRefCpgsMethylation.concept)
        self.ax.legend(loc=0)
        self.ax.set_xlim(0, 100)
        #3.Plot Save And Close
        self.saveAndClose()        
 
class SummaryMethylation(StatsMother):
    """Base Level Statistics"""
    
    def __init__(self):
        """ Initializes BaseLevel class members """
        StatsMother.__init__(self)

    def setData(self, concept = "", values=[]):
        """
            Sets Data For plotting a Summary Methylation
            
            concept - Concept associated: AllRefCpG | PassedRefCpG | AllNonRefCpG | PassedNonRefCpG
            values - List of values for each type of genotyped CpG 
        """   
        self.data[concept] = values

    def getMedian(self,concept="",totalCpGs=None):
        """
            Get Median Value for a Given concept of Genotyped CpG

            concept - Concept associated: AllRefCpg | PassedRefCpg | AllNonRefCpg | PassedNonRefCpg
            totalCpGs - Total Number of CpGs
        """
        #1. Get Median CpGs
        medianCpG = []
        if totalCpGs % 2 == 0:
            medianCpG.append(int(totalCpGs/2))
            medianCpG.append(medianCpG[0] + 1)
        else:
            medianCpG.append(int(totalCpGs/2))
            
        #1. Get Methylation Value
        counted_cpgs = 0
        median_methylation = []
             
        meth_value = 0
        for cpgs in self.data[concept]:
            counted_cpgs += cpgs
            if len(medianCpG) == 2:
                if counted_cpgs >= medianCpG[0]:
                    if len(median_methylation) == 0:
                        median_methylation.append(meth_value)
                elif counted_cpgs >= medianCpG[1]:
                    return (float(median_methylation[0] + meth_value) / 2)/100
                    break
            else:
                if counted_cpgs >= medianCpG[0]:
                    return float(meth_value)/100
            meth_value = meth_value + 1
        return 0.0
            
        
    def getMethylationStatusList(self,concept=""):
        """ Get a list of Methylation Status 

            concept - Concept associated: AllRefCpg | PassedRefCpg | AllNonRefCpg | PassedNonRefCpg
        """  
        unmethylated = 0
        intermediate_methylated = 0        
        methylated = 0
        for i in range(len(self.data[concept])):
            if i>= 0 and i<31:
                unmethylated += self.data[concept][i]
            elif i >= 31 and i < 70:
                intermediate_methylated += self.data[concept][i]
            else:
                methylated += self.data[concept][i]
        return [unmethylated,intermediate_methylated,methylated]
    
    def getTable(self):
        """ Get vector of rows to be printed in a table report"""
        vector_table = []
        vector_table.append(["Type","#AllReferenceCpG","%","#PassedReferenceCpG","%","#AllNonReferenceCpG","%","#PassedNonReferenceCpG","%"])

        total_all_reference_cpg = self.getTotalSumatory(self.data["AllRefCpg"])
        total_passed_reference_cpg = self.getTotalSumatory(self.data["PassedRefCpg"])
        total_all_non_reference_cpg = self.getTotalSumatory(self.data["AllNonRefCpg"])
        total_passed_non_reference_cpg = self.getTotalSumatory(self.data["PassedNonRefCpg"])

        list_all_reference_cpg = self.getMethylationStatusList(concept="AllRefCpg")
        list_passed_reference_cpg = self.getMethylationStatusList(concept="PassedRefCpg")
        list_all_non_reference_cpg = self.getMethylationStatusList(concept="AllNonRefCpg")
        list_passed_non_reference_cpg = self.getMethylationStatusList(concept="PassedNonRefCpg")

        #Unmethylated
        vector_table.append(
                             [
                               "Unmethylated (0-0.3)",
                               "%i" %(list_all_reference_cpg[0]),"%.2f %%" %(self.getPercentage(list_all_reference_cpg[0],total_all_reference_cpg)),
                               "%i" %(list_passed_reference_cpg[0]),"%.2f %%" %(self.getPercentage(list_passed_reference_cpg[0],total_passed_reference_cpg)),
                               "%i" %(list_all_non_reference_cpg[0]),"%.2f %%" %(self.getPercentage(list_all_non_reference_cpg[0],total_all_non_reference_cpg)),
                               "%i" %(list_passed_non_reference_cpg[0]),"%.2f %%" %(self.getPercentage(list_passed_non_reference_cpg[0],total_passed_non_reference_cpg))
                             ]
                           )

        #Intermediate
        vector_table.append(
                             [
                               "Intermediate (0.3-0.7)",
                               "%i" %(list_all_reference_cpg[1]),"%.2f %%" %(self.getPercentage(list_all_reference_cpg[1],total_all_reference_cpg)),
                               "%i" %(list_passed_reference_cpg[1]),"%.2f %%" %(self.getPercentage(list_passed_reference_cpg[1],total_passed_reference_cpg)),
                               "%i" %(list_all_non_reference_cpg[1]),"%.2f %%" %(self.getPercentage(list_all_non_reference_cpg[1],total_all_non_reference_cpg)),
                               "%i" %(list_passed_non_reference_cpg[1]),"%.2f %%" %(self.getPercentage(list_passed_non_reference_cpg[1],total_passed_non_reference_cpg))
                             ]
                           )

        #Methylated
        vector_table.append(
                             [
                               "Methylated (0.7-1)",
                               "%i" %(list_all_reference_cpg[2]),"%.2f %%" %(self.getPercentage(list_all_reference_cpg[2],total_all_reference_cpg)),
                               "%i" %(list_passed_reference_cpg[2]),"%.2f %%" %(self.getPercentage(list_passed_reference_cpg[2],total_passed_reference_cpg)),
                               "%i" %(list_all_non_reference_cpg[2]),"%.2f %%" %(self.getPercentage(list_all_non_reference_cpg[2],total_all_non_reference_cpg)),
                               "%i" %(list_passed_non_reference_cpg[2]),"%.2f %%" %(self.getPercentage(list_passed_non_reference_cpg[2],total_passed_non_reference_cpg))
                             ]
                           )
                           
        #separator
        vector_table.append(
                             [
                               "",
                               "","",
                               "","",
                               "","",
                               "",""
                             ]
                           )

        #Median
        vector_table.append(
                             [
                               "Median Methylation",
                               "%.3f" %(self.getMedian(concept="AllRefCpg",totalCpGs=total_all_reference_cpg)),"",
                               "%.3f" %(self.getMedian(concept="PassedRefCpg",totalCpGs=total_passed_reference_cpg)),"",
                               "%.3f" %(self.getMedian(concept="AllNonRefCpg",totalCpGs=total_all_non_reference_cpg)),"",
                               "%.3f" %(self.getMedian(concept="PassedNonRefCpg",totalCpGs=total_passed_non_reference_cpg)),""
                             ]
                           )

        return vector_table
              
class NonCpGReadProfile(PlotMother):
    """Manage Non CpG Read Profile"""
    
    def __init__(self,pngFile):
        """
           Initialize vector of values
        """
        PlotMother.__init__(self,pngFile=pngFile,concept="")
        self.vector_values = []        
        
    def add(self, vector_values):
        """
            Updates vector of values with new values
            
            vector_values - New vector of values
        """
        if len(self.vector_values) == 0:
            self.vector_values = vector_values
        else:
            position = 0
            for readPosition in vector_values:
                for i in range(len(readPosition)):
                    self.vector_values[position][i] += readPosition[i]
                position += 1
        self.setHaveCounts(True)
                                
    def getNonConversionRatio(self):
        """
            Get Vector For Non Conversion

        """
        vectorNonConversion = []
        for readInformation in self.vector_values:
            percentage_non_coversion = 0
            total_non_conversion_info = readInformation[0]+readInformation[1]
            if total_non_conversion_info != 0:
                percentage_non_coversion = float(float(readInformation[0])/float(total_non_conversion_info))*100
            vectorNonConversion.append(percentage_non_coversion)
        
        return vectorNonConversion
                          
    def plot(self):
        """
            Builds a NonCpGReadProfile plot
        """
        #1. Line Space  
        total_values = len(self.vector_values) 
        x = np.linspace(0,total_values,num=total_values)
        #2. Plot Configuration
        self.initSubplots(figure_size=(5.7,5),yText="% Non-Conversion",xText="Read Position")
        #2.1 Plot Building        
        self.ax.set_title("% Non-Conversion at Non-CpG Sites")
        plt.plot(x,self.getNonConversionRatio(), '--', linewidth=2,color="magenta")
        #3. Axis Limit
        self.ax.set_xlim(0, total_values)
        #4. Save and close
        self.saveAndClose()     
  
class GCcoverage(StatsMother,PlotMother):
    """
        Manages the plot for GC /Coverage
    """
    
    def __init__(self,pngFile):
        """ Initialize GCCoverage Values """
        PlotMother.__init__(self,pngFile=pngFile,concept="")
        self.coverage_gc_bases = {}
        self.x_vector = []
        self.y_vector = []
        self.z_vector = []
        
    def add(self,values):
        """
            Add Values
            
            values - GC/Coverage Value
        """
        for coverage in values:
            int_coverage = int(coverage)
            if int_coverage in self.coverage_gc_bases:
                for i in range(len(values[coverage])):
                    self.coverage_gc_bases[int_coverage][i] += values[coverage][i]
            else:
                self.coverage_gc_bases[int_coverage] = values[coverage]
        self.setHaveCounts(True)
                
    def cleanNotSignificantTail(self):
        """
           Clean from GC/Coverage values the tail of coverages with low number of bases, less than 5%
        """
        #1. Get Total Number of bases
        totalBases = 0
        for coverage in self.coverage_gc_bases:
            totalBases += self.getTotalSumatory(self.coverage_gc_bases[coverage])

        #2. Get Until 95 and Cut
        new_coverage_gc_bases = {}
        totalProcessedCoverages =  0                
        for coverage in sorted(self.coverage_gc_bases):
            if self.getPercentage(totalProcessedCoverages,totalBases) < 95:
                new_coverage_gc_bases [coverage] = self.coverage_gc_bases[coverage]
                totalProcessedCoverages += self.getTotalSumatory(new_coverage_gc_bases[coverage])
            else:
                break                       
                        
        self.coverage_gc_bases = new_coverage_gc_bases
        
    def selectDataToPlot(self):
        """
            Selectd the data to be plotted
        """
        #0. Clean not significant tail of values
        self.cleanNotSignificantTail()
        #Inspired in: https://stackoverflow.com/questions/22712219/heat-map-using-matplotlib
        #1.Vector To Represent X axis. GC Coverage From 0% to 100%
        x_gc_percentage = range(101)
        #2.Vector to Represent Y axis. Coverage From 0 to len keys on coverage_gc_bases
        list_keys = self.coverage_gc_bases.keys()
        list_keys.sort()        
        maximumCoverage = list_keys[-1]
        y_coverage = range(maximumCoverage + 1)
        #3.Form matrix for the Heatmap
        z = []
        for coverage in y_coverage:
            if coverage in self.coverage_gc_bases:
                z.append(self.coverage_gc_bases[coverage])
            else:
                z.append([0 for x in range(101)])
                
        self.x_vector = x_gc_percentage
        self.y_vector = y_coverage
        self.z_vector = z
                
    def plot(self):
        """
            GC/Coverage Heatmap plot
        """
        #1. Create Plot
        self.initSubplots(figure_size=(5.7,5),yText="Coverage",xText="GC%")        
        
        self.ax.set_title("GC/Coverage Heatmap")
        
        im = self.ax.pcolormesh(self.x_vector,self.y_vector, self.z_vector)
        cbar = self.figure.colorbar(im)
        cbar.ax.set_ylabel("# Sites")
        
        #4.1 Axis Configuration        
        self.ax.axis('tight')

        self.saveAndClose()         
         

    def getCorrelationCoeficient(self):
        """
            Get Correlation Coeficient From Heatmap plotted values
            
            returns pearson correlation or -2 if it was not possible
        """
        #Estimate Sx, Sx2, Sy, Sy2,n
        yLen = len(self.y_vector)
        xLen = len(self.x_vector)
        sx = 0
        sx2 = 0
        sy = 0
        sy2 = 0
        sxy = 0
        n = 0
        for y in range(0,yLen):
            for x in range(0,xLen):
                z = self.z_vector[y][x]               
                sx += x*z
                sx2 += (x*x)*z
                sy += y*z
                sy2 += (y*y)*z
                sxy += z*x*y
                n += z
        #Estimate correlation
        if n != 0:
            numerator = float(sxy) - ((float(sx)*float(sy))/n)
            denominator = (float(sx2) - ((float(sx)*float(sx))/n)) * (float(sy2) - ((float(sy)*float(sy))/n))
            return float(numerator) / math.sqrt(denominator)
        else:
            return -2
        
            
class QCDistribution(DistributionPlot,StatsMother):
    """Class to plot QCDistributions

       Fisher Strand: Strand bias estimated using Fisher's Exact Test.
                      Strand bias is a type of sequencing bias in which one DNA strand is favored over the other, 
                      which can result in incorrect evaluation of the amount of evidence observed for one allele vs. the other.
                      The FisherStrand annotation is one of several methods that aims to evaluate whether there is strand bias in the data. 
                      It uses Fisher's Exact Test to determine if there is strand bias between forward and reverse strands for the reference or alternate allele.
                      The output is a Phred-scaled p-value. The higher the output value, the more likely there is to be bias. More bias is indicative of false positive calls.
                      
       Quality By Depth: Allele-specific call confidence normalized by depth of sample reads supporting the allele.
                         This annotation puts the variant confidence QUAL score into perspective by normalizing for the amount of coverage available. 
                         Because each read contributes a little to the QUAL score, variants in regions with deep coverage can have artificially inflated QUAL scores, 
                         giving the impression that the call is supported by more evidence than it really is. To compensate for this, we normalize the variant confidence by depth, 
                         which gives us a more objective picture of how well supported the call is.
      
       RMS Mapping Quality: Allele-specific Root Mean Square of the mapping quality of reads across all samples.
                            This annotation provides an estimation of the mapping quality of reads supporting each alternate allele in a variant call. 
                            Depending on the tool it is called from, it produces either raw data (sum of squared MQs) or the calculated root mean square.
                            The raw data is used to accurately calculate the root mean square when combining more than one sample. 
                            
       Goodness of Fit: Phred scaled goodness of fit.Likelihood-ratio test (LR) for best genotype call against best call with free allele frequency parameter.
                        The goodness of fit of a statistical model describes how well it fits a set of observations. 
                        Measures of goodness of fit typically summarize the discrepancy between observed values and the values expected under the model in question. 
    """

    def __init__(self,concept="",typeDistribution="FisherStrand",typeBaseLocation="Variant",pngFile=""):
        """ Initialize GCCoverage Values 
            
            pngFile -- png plot file
            concept -- Plot Concept
            typeDistribution -- Can be: "FisherStrand" | "QualityByDepth" | "RMSMappingQuality" | "GoodnessOfFit"
            typeBaseLocation -- Can be: "" | "Variant" | "NonVariant"
        """
        DistributionPlot.__init__(self,pngFile=pngFile,concept=concept)
        self.type_distribution = typeDistribution
        self.type_base_location = typeBaseLocation
        self.x_legend = self.type_distribution
        self.y_legend = ""
        self.percentage_limit_tail = 95
        if typeBaseLocation == "":
            self.y_legend += "Variant Sites"
        else:
            self.y_legend += "%s Sites" %(typeBaseLocation)
            
        #Number of Locations to be recovered in order to have the same axis for variants and non variants
        self.gof_locations = 0
        
    def setAxisXLabel(self,newLabel=""):
        """
            Updates X Label
            
            newLabel -- New Label
        """
        self.x_legend = newLabel        
                
    def add(self,values):
        """
            Add Values
            
            values - QCDistributions Value
        """
        for value in values:
            int_value = int(value)
            if int_value in self.value_locations:
                if self.type_base_location == "":
                    self.value_locations[int_value] += values[value] 
                else:
                    self.value_locations[int_value] += values[value][self.type_base_location]
            else:
                if self.type_base_location == "":
                    self.value_locations[int_value] = values[value]
                else:
                    self.value_locations[int_value] = values[value][self.type_base_location]
                       
        self.setHaveCounts(True)
        
    def getBarColor(self):
        """Get color according to the type of distribution"""
        if self.type_distribution == "FisherStrand":
            return "#ef8f8f"
        elif self.type_distribution == "QualityByDepth":
            return "#9bef8f"
        elif self.type_distribution == "RMSMappingQuality":
            return "#28e7bb"
        elif self.type_distribution == "GoodnessOfFit":
            return "#dbb028"
        return "#ef8f8f"

    def selectTotalNumberLocations(self,other):
        """
            Compare two different QCDistributions objects and get the largest one to unify X axis of two plots
            This method is designed to have the same x scale on GoodnessOfFit plots
            
            other - QCDistributions object to compare with
            returns total number of locations to be plotted
        """
        self.percentage_limit_tail = 98
        own_locations = len(self.getVectorToPlot(True))
        other_locations = len(other.getVectorToPlot(True))
        if own_locations >= other_locations:
            return own_locations
        else:
            return other_locations
            
    def setNumberOfLocations(self,locationsToRecover):
        """
            Updates the number of locations to be recovered
            
            locationsToRecover - New Locations Value
        """            
        self.gof_locations = locationsToRecover
                 
    def getUnifiedVectorToPlot(self,locationsToRecover):
        """From Dictionary Get Vector of Y values to be bar plotted
           Method designed to have the same bar for both plots of GoodnessOfFit Variants and NonVariants
        
           locationsToRecover - Total number of locations to recover from the value_locations
        """
        vector_plot = []
        for x in range(0,locationsToRecover):
            if x in self.value_locations:
                vector_plot.append(self.value_locations[x])
            else:
                vector_plot.append(0)
        return vector_plot       
                       
    def plot(self):
        """Builds a QCDistribution Plot and print it to a file"""
        #Number of Bins
        widhtBin = 1
        
        #Create vector of Y values
        yValues = []
        
        if self.type_distribution == "FisherStrand" or self.type_distribution == "QualityByDepth":
            yValues = self.getVectorToPlot(cleanTail=True)
        elif self.type_distribution == "GoodnessOfFit":
            yValues = self.getUnifiedVectorToPlot(locationsToRecover=self.gof_locations)
        else: 
            yValues = self.getVectorToPlot(cleanTail=False)
                
        #Create Vector of X ticks
        xValues = range(len(yValues))
        
        #Build Plot           
        self.initSubplots(figure_size=(5.7,5),yText=self.y_legend,xText=self.x_legend)           
                     
        width = widhtBin     
        self.ax.bar(xValues, yValues, width, color=self.getBarColor())

        self.setYAxisExponential()
        
        self.saveAndClose()
          
class VCFFilterStats(StatsMother):
    """VCF Filter Statistics"""
    
    def __init__(self):
        """ Initializes BaseLevel class members """
        StatsMother.__init__(self)
        
    def add(self, values={}):
        """
            Add (Sum) new Values
            
            values - Dictionary of values for Base Level Stats
        """ 
        for filtering_criteria in values:
            if filtering_criteria in self.data:
                self.data[filtering_criteria]["Variant"] += values[filtering_criteria]["Variant"]
                self.data[filtering_criteria]["NonVariant"] += values[filtering_criteria]["NonVariant"]
            else:
                self.data[filtering_criteria] = {"Variant": values[filtering_criteria]["Variant"],"NonVariant": values[filtering_criteria]["NonVariant"]}
                   
    def getTable(self):
        """ Get vector of rows to be printed in a table report"""
        vector_table = []

        #0.Get Total Values
        total_passed_variants = 0
        total_passed_non_variants = 0
        total_filtered_variants = 0
        total_filtered_non_variants = 0
        #0.1. Loop over values
        for filtering_criteria in self.data:
            if filtering_criteria == "PASS":
                total_passed_variants += self.data[filtering_criteria]["Variant"]
                total_passed_non_variants += self.data[filtering_criteria]["NonVariant"]
            else:
                total_filtered_variants += self.data[filtering_criteria]["Variant"]
                total_filtered_non_variants += self.data[filtering_criteria]["NonVariant"]
                
        total_non_variants = total_passed_non_variants + total_filtered_non_variants
        total_variants = total_passed_variants + total_filtered_variants
        total_passed = total_passed_non_variants + total_passed_variants
        total_filtered = total_filtered_non_variants + total_filtered_variants
        total_sites = total_variants + total_non_variants
        
        #1.Table Header
        vector_table.append(["Type","#Sites","%","#Non-Variant Sites","%","#Variant Sites","%"])
        #2.All Values
        vector_table.append(["All","%i"%(total_sites),"100%",
                             "%i"%(total_non_variants),"%.2f %%"%(self.getPercentage(total_non_variants,total_sites)),
                             "%i"%(total_variants),"%.2f %%"%(self.getPercentage(total_variants,total_sites))])
        #3.Empty Row (Separator Row)
        vector_table.append(["","","","","","",""])
        #4.PASS AND FILTERED
        vector_table.append(["Passed","%i"%(total_passed),"%.2f %%"%(self.getPercentage(total_passed,total_sites)),
                             "%i"%(total_passed_non_variants),"%.2f %%"%(self.getPercentage(total_passed_non_variants,total_non_variants)),
                             "%i"%(total_passed_variants),"%.2f %%"%(self.getPercentage(total_passed_variants,total_passed))])
        vector_table.append(["Filtered","%i"%(total_filtered),"%.2f %%"%(self.getPercentage(total_filtered,total_sites)),
                             "%i"%(total_filtered_non_variants),"%.2f %%"%(self.getPercentage(total_filtered_non_variants,total_non_variants)),
                             "%i"%(total_filtered_variants),"%.2f %%"%(self.getPercentage(total_filtered_variants,total_passed))])
        #5.Empty Row (Separator Row)
        vector_table.append(["","","","","","",""])
        #6.ROW FILTERED MOTIVES
        filtered_citerias = ["q20","qd2","q20,qd2", "fs60","q20,fs60","qd2,fs60","q20,qd2,fs60","mq40","q20,mq40","qd2,mq40","q20,qd2,mq40","fs60,mq40","q20,fs60,mq40",
                              "qd2,fs60,mq40","q20,qd2,fs60,mq40","gof20","q20,gof20","qd2,gof20","q20,qd2,gof20","fs60,gof20","q20,fs60,gof20","qd2,fs60,gof20","q20,qd2,fs60,gof20",
                              "mq40,gof20","q20,mq40,gof20","qd2,mq40,gof20","q20,qd2,mq40,gof20","fs60,mq40,gof20","q20,fs60,mq40,gof20","qd2,fs60,mq40,gof20","q20,qd2,fs60,mq40,gof20"]
         
        tuples_vector = []                     
        for filtered_criteria in filtered_citerias:
            criteria_non_variant = self.data[filtered_criteria]["NonVariant"]
            criteria_variant = self.data[filtered_criteria]["Variant"]
            criteria_total = criteria_non_variant + criteria_variant
            tuples_vector.append([filtered_criteria,criteria_total])
            
        def getKey(item):
            return item[1]
             
        for criteria_pair in sorted(tuples_vector, key=getKey,reverse=True): 
            criteria_non_variant = self.data[criteria_pair[0]]["NonVariant"]
            criteria_variant = self.data[criteria_pair[0]]["Variant"]
            criteria_total = criteria_non_variant + criteria_variant   
            
            vector_table.append([criteria_pair[0],"%i"%(criteria_total),"%.2f %%"%(self.getPercentage(criteria_total,total_filtered)),
                             "%i"%(criteria_non_variant),"%.2f %%"%(self.getPercentage(criteria_non_variant,total_filtered_non_variants)),
                             "%i"%(criteria_variant),"%.2f %%"%(self.getPercentage(criteria_variant,total_filtered_variants))])
                             
        return vector_table

class SummarySample(StatsMother):
    
    def __init__(self,sampleName="",readLevelStats=None,baseLevelStats=None,gcCoverage=None,totalStats=None,variantCoverage=None,
                 mutationStats=None,methylationPassRefCpg=None,refCpgCoverage=None):
        """ Summary Sample Constructor """
        self.sampleName = sampleName
        self.readLevelStats = readLevelStats
        self.baseLevelStats = baseLevelStats
        self.gcCoverage = gcCoverage
        self.totalStats = totalStats
        self.variantCoverage = variantCoverage
        self.mutationStats = mutationStats
        self.methylationPassRefCpg = methylationPassRefCpg
        self.refCpgCoverage = refCpgCoverage
  
    def fromBasesToEasyToRead(self,bases=0):
        """ 
            From a given number of bases, get the easy reduction to be read.
            
            bases  -  Number of bases to be transformed
            returns string with the simplified number and its units
        """
        gigaBase = 1000000000
        megaBase = 1000000
        kiloBase = 1000
        if bases >= gigaBase:
            return "%.2f Gb"%(float(bases)/float(gigaBase))
        elif bases >= megaBase:
            return "%.2f Mb"%(float(bases)/float(megaBase))
        else:
            return "%.2f Kb"%(float(bases)/float(kiloBase))
        
    def fromGenomeEventsToEasyToRead(self,events=0):
        """
           From a given number of genome events, get the easy reduction to be read.
           
           events  -  Number of events to be transformed
           returns string with the simplified number and its units
        """
        million = 1000000
        kilo = 1000
        if events >= million:
            return "%.2f M"%(float(events)/float(million))
        elif events >= kilo:
            return "%.2f K"%(float(events)/float(kilo))
        else:
            return "%i"%(events)
      
    def getTable(self):
        """Get Table For Printing Summary Sample Table """
        #1.Define Table Alignment And Coverage
        aligmentAndCoverage = []
        #1.1 Bases Aligned
        totalBases = self.readLevelStats.getTotalBases()
        aligmentAndCoverage.append("%s"%(self.fromBasesToEasyToRead(bases=totalBases)))
        #1.2 Bases Uniquely Aligned
        uniquely = self.readLevelStats.passed.data["Bases"]
        aligmentAndCoverage.append("%s (%.2f %%)"%(self.fromBasesToEasyToRead(bases=uniquely),self.getPercentage(uniquely,totalBases)))
        #1.3 Bases Used for Calling
        passed = self.baseLevelStats.data["Passed"]
        aligmentAndCoverage.append("%s (%.2f %%)"%(self.fromBasesToEasyToRead(bases=passed),self.getPercentage(passed,totalBases)))
        #1.4 Get Correlation coeficient for GC Coverage
        correlation = self.gcCoverage.getCorrelationCoeficient()
        aligmentAndCoverage.append("%.2f"%(correlation))
        #2. Variants
        variantsAndCoverage = []
        #2.1 Variants called
        allSnps = self.totalStats.snps.data["All"]
        variantsAndCoverage.append("%s"%(self.fromGenomeEventsToEasyToRead(events=allSnps)))
        #2.2 Variants High Quality
        snps_passed = self.totalStats.snps.data["Passed"]
        variantsAndCoverage.append("%s (%.2f %%)"%(self.fromGenomeEventsToEasyToRead(events=snps_passed),self.getPercentage(snps_passed,allSnps)))
        #2.3 Mean Coverage Variants High Quality
        variantsAndCoverage.append("%.2f X" %(self.variantCoverage.getMean()))
        #2.4 Ti/Tv Ration Passed
        variantsAndCoverage.append("%s" %(self.mutationStats.getTiTv("Passed")[0]))
        #3. Methylation
        methylation = []
        #3.1 Mean CpG Methylation
        methylation.append("%.2f"%(self.methylationPassRefCpg.getMean()))
        #3.2 Mean CpG Coverage
        methylation.append("%.2f X"%(self.refCpgCoverage.getMean()))
        #3.3. Total Stats Passed
        cpgs_passed = self.totalStats.refCpGs.data["Passed"]
        methylation.append("%s"%(self.fromGenomeEventsToEasyToRead(events=cpgs_passed)))

        #4.Define return value
        return {"alignments":aligmentAndCoverage,"variants":variantsAndCoverage,"methylation":methylation} 

        
        