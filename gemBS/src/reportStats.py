# -*- coding: utf-8 -*-
#!/usr/bin/env python
import json

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

from matplotlib import cm
import matplotlib.colors as colors

#from mpl_toolkits.mplot3d import Axes3D

import importlib
importlib.import_module('mpl_toolkits.mplot3d').Axes3D

import math

class NucleotideStats(object):
    """ Gets percentage of nucleotide statistics """
    
    def __init__(self,a_bp,c_bp,t_bp,g_bp,n_bp):
        """ From A,C,G,T total nucleotides get percentage statistics per each nucleotide """

        a_bp.sumPair()
        c_bp.sumPair()
        g_bp.sumPair()
        t_bp.sumPair()
        n_bp.sumPair()

        self.all_total = a_bp.total + c_bp.total + g_bp.total + t_bp.total + n_bp.total     
        self.pair_one_total = a_bp.pair_one + c_bp.pair_one + g_bp.pair_one + t_bp.pair_one + n_bp.pair_one
        self.pair_two_total = a_bp.pair_two + c_bp.pair_two + g_bp.pair_two + t_bp.pair_two + n_bp.pair_two

        self.all_a_percentage_bp = self.getPercentage(self.all_total,a_bp.total)
        self.pair_one_a_percentage_bp = self.getPercentage(self.pair_one_total,a_bp.pair_one)
        self.pair_two_a_percentage_bp = self.getPercentage(self.pair_two_total,a_bp.pair_two)

        self.all_c_percentage_bp = self.getPercentage(self.all_total,c_bp.total)
        self.pair_one_c_percentage_bp = self.getPercentage(self.pair_one_total,c_bp.pair_one)
        self.pair_two_c_percentage_bp = self.getPercentage(self.pair_two_total,c_bp.pair_two)

        self.all_g_percentage_bp = self.getPercentage(self.all_total,g_bp.total) 
        self.pair_one_g_percentage_bp = self.getPercentage(self.pair_one_total,g_bp.pair_one)
        self.pair_two_g_percentage_bp = self.getPercentage(self.pair_two_total,g_bp.pair_two)

        self.all_t_percentage_bp = self.getPercentage(self.all_total,t_bp.total)
        self.pair_one_t_percentage_bp = self.getPercentage(self.pair_one_total,t_bp.pair_one)
        self.pair_two_t_percentage_bp = self.getPercentage(self.pair_two_total,t_bp.pair_two)

        self.all_n_percentage_bp = self.getPercentage(self.all_total,n_bp.total)
        self.pair_one_n_percentage_bp = self.getPercentage(self.pair_one_total,n_bp.pair_one)
        self.pair_two_n_percentage_bp = self.getPercentage(self.pair_two_total,n_bp.pair_two)
        
    def getPercentage(self,total,subtotal):
        if total != 0.0:
            return (float(subtotal)/float(total))*100
        else:
            return 0.0
        

class Value(object):
    """ Class value to define a given metric, could be total or per pairs """
    def __init__(self):
        self.total = 0
        self.pair_one = 0
        self.pair_two = 0
        self.total_percentage = 0.0
        self.pair_one_percentage = 0.0
        self.pair_two_percentage = 0.0
        
    def sumPair(self):
        self.total = self.pair_one + self.pair_two
        
    def sumValues(self,other):
        self.total = self.total + other.total
        self.pair_one = self.pair_one + other.pair_one
        self.pair_two = self.pair_two + other.pair_two

    def getPercentage(self,total,subtotal):
        if total != 0.0:
            return (float(subtotal)/float(total))*100
        else:
            return 0.0
        
    def setPercentages(self,total,total_pair_one,total_pair_two):
        self.sumPair()
        self.total_percentage = self.getPercentage(total,self.total)
        self.pair_one_percentage = self.getPercentage(total_pair_one,self.pair_one)
        self.pair_two_percentage = self.getPercentage(total_pair_two,self.pair_two)
        
    def setPercentagesNucleotides(self,percent_all,percent_one,percent_two):
        self.total_percentage = percent_all
        self.pair_one_percentage = percent_one
        self.pair_two_percentage = percent_two
        

class BsStats(object):
    
    def __init__(self):
        #Is paired or single end
        self.is_paired = True
        #Sequenced reads       
        self.sequenced_reads = Value()
        #General reads       
        self.reads = Value()
        #Reads in Control sequences
        self.reads_sequencing_control = Value()
        #Reads under conversion control        
        self.reads_under_conversion_control = Value()        
        #Reads over conversion control
        self.reads_over_conversion_control = Value()
        #Unmapped reads
        self.reads_unmapped = Value()
        #Bisulfite_reads
        self.reads_reference_C2T = Value()
        self.reads_reference_G2A = Value()
        #Correct pairs
        self.correct_pairs = 0
        #Base counts
        #Overall
        self.reads_overall_A = Value()
        self.reads_overall_C = Value()
        self.reads_overall_G = Value()
        self.reads_overall_T = Value()
        self.reads_overall_N = Value()
        #GeneralC2T
        self.reads_general_C2T_A = Value()
        self.reads_general_C2T_C = Value()
        self.reads_general_C2T_G = Value()
        self.reads_general_C2T_T = Value()
        self.reads_general_C2T_N = Value()
        #GeneralG2A       
        self.reads_general_G2A_A = Value()
        self.reads_general_G2A_C = Value()
        self.reads_general_G2A_G = Value()
        self.reads_general_G2A_T = Value()
        self.reads_general_G2A_N = Value()
        #UnderConversionControlC2T
        self.reads_under_conversion_control_C2T_A = Value()  
        self.reads_under_conversion_control_C2T_C = Value()  
        self.reads_under_conversion_control_C2T_G = Value()  
        self.reads_under_conversion_control_C2T_T = Value()  
        self.reads_under_conversion_control_C2T_N = Value()  
        #UnderConversionControlG2A
        self.reads_under_conversion_control_G2A_A = Value()  
        self.reads_under_conversion_control_G2A_C = Value()  
        self.reads_under_conversion_control_G2A_G = Value()  
        self.reads_under_conversion_control_G2A_T = Value()  
        self.reads_under_conversion_control_G2A_N = Value()  
        #OverConversionControl_C2T
        self.reads_over_conversion_control_C2T_A = Value()
        self.reads_over_conversion_control_C2T_C = Value()
        self.reads_over_conversion_control_C2T_G = Value()
        self.reads_over_conversion_control_C2T_T = Value()
        self.reads_over_conversion_control_C2T_N = Value()
        #OverConversionControl_G2A
        self.reads_over_conversion_control_G2A_A = Value()
        self.reads_over_conversion_control_G2A_C = Value()
        self.reads_over_conversion_control_G2A_G = Value()
        self.reads_over_conversion_control_G2A_T = Value()
        self.reads_over_conversion_control_G2A_N = Value()
        
        #Vector History Mapq (From Mapping Quality 0 to Mapping Quality 60)
        self.mapping_quality_reads = []   
        self.lowestUniquelyQuality = 20  #Lowes mapping quality to consider a read as unique      
        
        #read_length_histogram
        self.read_length_histogram = {}
        self.read_insert_size_histogram = {}
        
    def assign_value(self,stat_value,vector_value):
        #Assign a json value to a stats value
        if self.is_paired:
            stat_value.pair_one = int(vector_value[0])
            stat_value.pair_two = int(vector_value[1])
        else:
            stat_value.pair_one = int(vector_value[0])
            
            
    def sum_values(self, other):
        """Sum all values"""
        #Sequenced Reads
        self.sequenced_reads.sumValues(other.sequenced_reads)   
        #General reads       
        self.reads.sumValues(other.reads)
        #Reads in Control sequences
        self.reads_sequencing_control.sumValues(other.reads_sequencing_control)
        #Reads under conversion control        
        self.reads_under_conversion_control.sumValues(other.reads_under_conversion_control)        
        #Reads over conversion control
        self.reads_over_conversion_control.sumValues(other.reads_over_conversion_control)
        #Unmapped reads
        self.reads_unmapped.sumValues(other.reads_unmapped)
        #Bisulfite_reads
        self.reads_reference_C2T.sumValues(other.reads_reference_C2T)
        self.reads_reference_G2A.sumValues(other.reads_reference_G2A)
        #Correct pairs
        self.correct_pairs = self.correct_pairs + other.correct_pairs
        #Base counts
        #Overall
        self.reads_overall_A.sumValues(other.reads_overall_A)
        self.reads_overall_C.sumValues(other.reads_overall_C)
        self.reads_overall_G.sumValues(other.reads_overall_G)
        self.reads_overall_T.sumValues(other.reads_overall_T)
        self.reads_overall_N.sumValues(other.reads_overall_N)
        #GeneralC2T
        self.reads_general_C2T_A.sumValues(other.reads_general_C2T_A)
        self.reads_general_C2T_C.sumValues(other.reads_general_C2T_C)
        self.reads_general_C2T_G.sumValues(other.reads_general_C2T_G)
        self.reads_general_C2T_T.sumValues(other.reads_general_C2T_T)
        self.reads_general_C2T_N.sumValues(other.reads_general_C2T_N)
        #GeneralG2A       
        self.reads_general_G2A_A.sumValues(other.reads_general_G2A_A)
        self.reads_general_G2A_C.sumValues(other.reads_general_G2A_C)
        self.reads_general_G2A_G.sumValues(other.reads_general_G2A_G)
        self.reads_general_G2A_T.sumValues(other.reads_general_G2A_T)
        self.reads_general_G2A_N.sumValues(other.reads_general_G2A_N)
        #UnderConversionControlC2T
        self.reads_under_conversion_control_C2T_A.sumValues(other.reads_under_conversion_control_C2T_A) 
        self.reads_under_conversion_control_C2T_C.sumValues(other.reads_under_conversion_control_C2T_C)  
        self.reads_under_conversion_control_C2T_G.sumValues(other.reads_under_conversion_control_C2T_G)  
        self.reads_under_conversion_control_C2T_T.sumValues(other.reads_under_conversion_control_C2T_T) 
        self.reads_under_conversion_control_C2T_N.sumValues(other.reads_under_conversion_control_C2T_N)  
        #UnderConversionControlG2A
        self.reads_under_conversion_control_G2A_A.sumValues(other.reads_under_conversion_control_G2A_A)
        self.reads_under_conversion_control_G2A_C.sumValues(other.reads_under_conversion_control_G2A_C)
        self.reads_under_conversion_control_G2A_G.sumValues(other.reads_under_conversion_control_G2A_G)
        self.reads_under_conversion_control_G2A_T.sumValues(other.reads_under_conversion_control_G2A_T)
        self.reads_under_conversion_control_G2A_N.sumValues(other.reads_under_conversion_control_G2A_N)
        #OverConversionControl_C2T
        self.reads_over_conversion_control_C2T_A.sumValues(other.reads_over_conversion_control_C2T_A)
        self.reads_over_conversion_control_C2T_C.sumValues(other.reads_over_conversion_control_C2T_C)
        self.reads_over_conversion_control_C2T_G.sumValues(other.reads_over_conversion_control_C2T_G)
        self.reads_over_conversion_control_C2T_T.sumValues(other.reads_over_conversion_control_C2T_T)
        self.reads_over_conversion_control_C2T_N.sumValues(other.reads_over_conversion_control_C2T_N)
        #OverConversionControl_G2A
        self.reads_over_conversion_control_G2A_A.sumValues(other.reads_over_conversion_control_G2A_A)
        self.reads_over_conversion_control_G2A_C.sumValues(other.reads_over_conversion_control_G2A_C)
        self.reads_over_conversion_control_G2A_G.sumValues(other.reads_over_conversion_control_G2A_G)
        self.reads_over_conversion_control_G2A_T.sumValues(other.reads_over_conversion_control_G2A_T)
        self.reads_over_conversion_control_G2A_N.sumValues(other.reads_over_conversion_control_G2A_N) 
        
        
    def getUnderConversionRate(self):   
        """Get Under Conversion Rate"""
        #A
        a_bp_pair_one = self.reads_under_conversion_control_C2T_A.pair_one + self.reads_under_conversion_control_G2A_A.pair_one
        a_bp_pair_two = self.reads_under_conversion_control_C2T_A.pair_two + self.reads_under_conversion_control_G2A_A.pair_two
        #C
        c_bp_pair_one = self.reads_under_conversion_control_C2T_C.pair_one + self.reads_under_conversion_control_G2A_C.pair_one
        c_bp_pair_two = self.reads_under_conversion_control_C2T_C.pair_two + self.reads_under_conversion_control_G2A_C.pair_two
        #G
        g_bp_pair_one = self.reads_under_conversion_control_C2T_G.pair_one + self.reads_under_conversion_control_G2A_G.pair_one
        g_bp_pair_two = self.reads_under_conversion_control_C2T_G.pair_two + self.reads_under_conversion_control_G2A_G.pair_two
        #T
        t_bp_pair_one = self.reads_under_conversion_control_C2T_T.pair_one + self.reads_under_conversion_control_G2A_T.pair_one
        t_bp_pair_two = self.reads_under_conversion_control_C2T_T.pair_two + self.reads_under_conversion_control_G2A_T.pair_two
        
        return self.getConversionRate(a_bp_pair_one,a_bp_pair_two,c_bp_pair_one,c_bp_pair_two,g_bp_pair_one,g_bp_pair_two,t_bp_pair_one,t_bp_pair_two)
        
    def getOverConversionRate(self):
        """Get Over Conversion Rate"""
        #A
        a_bp_pair_one = self.reads_over_conversion_control_C2T_A.pair_one + self.reads_over_conversion_control_G2A_A.pair_one
        a_bp_pair_two = self.reads_over_conversion_control_C2T_A.pair_two + self.reads_over_conversion_control_G2A_A.pair_two
        #C
        c_bp_pair_one = self.reads_over_conversion_control_C2T_C.pair_one + self.reads_over_conversion_control_G2A_C.pair_one
        c_bp_pair_two = self.reads_over_conversion_control_C2T_C.pair_two + self.reads_over_conversion_control_G2A_C.pair_two
        #G
        g_bp_pair_one = self.reads_over_conversion_control_C2T_G.pair_one + self.reads_over_conversion_control_G2A_G.pair_one
        g_bp_pair_two = self.reads_over_conversion_control_C2T_G.pair_two + self.reads_over_conversion_control_G2A_G.pair_two
        #T
        t_bp_pair_one = self.reads_over_conversion_control_C2T_T.pair_one + self.reads_over_conversion_control_G2A_T.pair_one
        t_bp_pair_two = self.reads_over_conversion_control_C2T_T.pair_two + self.reads_over_conversion_control_G2A_T.pair_two
       
        return self.getConversionRate(a_bp_pair_one,a_bp_pair_two,c_bp_pair_one,c_bp_pair_two,g_bp_pair_one,g_bp_pair_two,t_bp_pair_one,t_bp_pair_two)
        
    def getConversionRate(self,a_bp_pair_one=0,a_bp_pair_two=0,c_bp_pair_one=0,c_bp_pair_two=0,g_bp_pair_one=0,g_bp_pair_two=0,t_bp_pair_one=0,t_bp_pair_two=0):
        """Get Over Conversion Rate"""
        if float(a_bp_pair_one + g_bp_pair_one + c_bp_pair_two + t_bp_pair_two) == 0:
            return 0.0
            
        z = float(a_bp_pair_one + t_bp_pair_two) / float(a_bp_pair_one + g_bp_pair_one + c_bp_pair_two + t_bp_pair_two)
        a = (t_bp_pair_one + a_bp_pair_two) * float(1.0 - z) - (c_bp_pair_one + g_bp_pair_two) * z 
        b = (t_bp_pair_one + a_bp_pair_two + c_bp_pair_one + g_bp_pair_two) * float(1.0 - z)

        return float(a)/float(b)        
        
    def getUniqueMappedReads(self):
        """ Get total Uniquely mapped reads
            Uniquely mapped reads are those set of reads that have a mapping quality higher than 20 """
        
        totalUniquelyReads = 0        
        highestQuality = len(self.mapping_quality_reads)
        
        for qual in range(self.lowestUniquelyQuality,highestQuality):
            totalUniquelyReads += self.mapping_quality_reads[qual]
            
        return totalUniquelyReads
        
    def getAverageUniqueMappedReads(self):
        """ Get Average Uniquely mapped reads
            Uniquely mapped reads are those set of reads that have a mapping quality higher than 20 """
        
        totalUniquelyReads = 0        
        totalMappedReads = 0
        highestQuality = len(self.mapping_quality_reads)
        
        for qual in range(0,highestQuality):
            if qual >= self.lowestUniquelyQuality:
                totalUniquelyReads += self.mapping_quality_reads[qual]
            totalMappedReads += self.mapping_quality_reads[qual]
            
            
        return (float(totalUniquelyReads)/float(totalMappedReads))*100
    
    def getTotalMappedReads(self):
        """ Get Total Mapped Reads """      
        totalMappedReads = 0
        highestQuality = len(self.mapping_quality_reads)
        for qual in range(0,highestQuality):
            totalMappedReads += self.mapping_quality_reads[qual]
        return totalMappedReads
        
class LaneStats(BsStats):
    """Statistics per lane """
    
    def __init__(self,name,json_file=None):
        
        #Call Parent class constructor
        super(LaneStats, self).__init__()
        
        self.name = name  
        
        #Parsing json file
        with open(json_file, 'r') as file_json:
            data = json.load(file_json) 
            #Process data
            if data["MapperType"] != "Paired":
                self.is_paired = False
            
            #Sequenced Reads
            sequenced_reads = []
            
            seq_reads_one  = int(data["Reads"]["General"][0])+int(data["Reads"]["Unmapped"][0]) 
                
            if "SequencingControl" in data["Reads"]:
                seq_reads_one  += int(data["Reads"]["SequencingControl"][0])  
            if "UnderConversionControl" in data["Reads"]:
                seq_reads_one  += int(data["Reads"]["UnderConversionControl"][0])  
            if "OverConversionControl" in data["Reads"]:
                seq_reads_one  += int(data["Reads"]["OverConversionControl"][0])
                
            sequenced_reads.append(seq_reads_one)

            if self.is_paired:  
                seq_reads_two  = int(data["Reads"]["General"][1])+int(data["Reads"]["Unmapped"][1]) 
                
                if "SequencingControl" in data["Reads"]:
                    seq_reads_two  += int(data["Reads"]["SequencingControl"][1])  
                if "UnderConversionControl" in data["Reads"]:
                    seq_reads_two  += int(data["Reads"]["UnderConversionControl"][1])  
                if "OverConversionControl" in data["Reads"]:
                    seq_reads_two  += int(data["Reads"]["OverConversionControl"][1])
                    
                sequenced_reads.append(seq_reads_two)
                        
            self.assign_value(self.sequenced_reads,sequenced_reads)
            #General reads       
            self.assign_value(self.reads,data["Reads"]["General"])
            #Reads in Control sequences
            if "SequencingControl" in data["Reads"]:
                self.assign_value(self.reads_sequencing_control,data["Reads"]["SequencingControl"])  
            #Reads under conversion control      
            if "UnderConversionControl" in data["Reads"]:
                self.assign_value(self.reads_under_conversion_control,data["Reads"]["UnderConversionControl"])
            #Reads over conversion control
            if "OverConversionControl" in data["Reads"]:
                self.assign_value(self.reads_over_conversion_control,data["Reads"]["OverConversionControl"])                  
            #Unmapped reads
            self.assign_value(self.reads_unmapped,data["Reads"]["Unmapped"]) 
            #Bisulfite_reads
            if "NumReadsBS" in data:
                self.assign_value(self.reads_reference_C2T,data["NumReadsBS"]["C2T"]) 
                self.assign_value(self.reads_reference_G2A,data["NumReadsBS"]["G2A"])
            #Base counts Overall
            self.assign_value(self.reads_overall_A,data["BaseCounts"]["Overall"]["A"])
            self.assign_value(self.reads_overall_C,data["BaseCounts"]["Overall"]["C"])
            self.assign_value(self.reads_overall_G,data["BaseCounts"]["Overall"]["G"])
            self.assign_value(self.reads_overall_T,data["BaseCounts"]["Overall"]["T"])
            self.assign_value(self.reads_overall_N,data["BaseCounts"]["Overall"]["N"])
            #Base counts GeneralC2T
            if "GeneralC2T" in data["BaseCounts"]:
                self.assign_value(self.reads_general_C2T_A,data["BaseCounts"]["GeneralC2T"]["A"])
                self.assign_value(self.reads_general_C2T_C,data["BaseCounts"]["GeneralC2T"]["C"])
                self.assign_value(self.reads_general_C2T_G,data["BaseCounts"]["GeneralC2T"]["G"])
                self.assign_value(self.reads_general_C2T_T,data["BaseCounts"]["GeneralC2T"]["T"])
                self.assign_value(self.reads_general_C2T_N,data["BaseCounts"]["GeneralC2T"]["N"])
            #Base counts GeneralG2A
            if "GeneralG2A" in data["BaseCounts"]:
                self.assign_value(self.reads_general_G2A_A,data["BaseCounts"]["GeneralG2A"]["A"])
                self.assign_value(self.reads_general_G2A_C,data["BaseCounts"]["GeneralG2A"]["C"])
                self.assign_value(self.reads_general_G2A_G,data["BaseCounts"]["GeneralG2A"]["G"])
                self.assign_value(self.reads_general_G2A_T,data["BaseCounts"]["GeneralG2A"]["T"])
                self.assign_value(self.reads_general_G2A_N,data["BaseCounts"]["GeneralG2A"]["N"])
            #Base Counts UnderConversionControlC2T
            if "UnderConversionControlC2T" in data["BaseCounts"]:
                self.assign_value(self.reads_under_conversion_control_C2T_A,data["BaseCounts"]["UnderConversionControlC2T"]["A"])
                self.assign_value(self.reads_under_conversion_control_C2T_C,data["BaseCounts"]["UnderConversionControlC2T"]["C"])
                self.assign_value(self.reads_under_conversion_control_C2T_G,data["BaseCounts"]["UnderConversionControlC2T"]["G"])
                self.assign_value(self.reads_under_conversion_control_C2T_T,data["BaseCounts"]["UnderConversionControlC2T"]["T"])
                self.assign_value(self.reads_under_conversion_control_C2T_N,data["BaseCounts"]["UnderConversionControlC2T"]["N"])
            #Base Counts UnderConversionControlG2A
            if "UnderConversionControlG2A" in data["BaseCounts"]:            
                self.assign_value(self.reads_under_conversion_control_G2A_A,data["BaseCounts"]["UnderConversionControlG2A"]["A"])
                self.assign_value(self.reads_under_conversion_control_G2A_C,data["BaseCounts"]["UnderConversionControlG2A"]["C"])
                self.assign_value(self.reads_under_conversion_control_G2A_G,data["BaseCounts"]["UnderConversionControlG2A"]["G"])
                self.assign_value(self.reads_under_conversion_control_G2A_T,data["BaseCounts"]["UnderConversionControlG2A"]["T"])
                self.assign_value(self.reads_under_conversion_control_G2A_N,data["BaseCounts"]["UnderConversionControlG2A"]["N"])
            #Base Counts OverConversionControl_C2T
            if "OverConversionControlC2T" in data["BaseCounts"]:            
                self.assign_value(self.reads_over_conversion_control_C2T_A,data["BaseCounts"]["OverConversionControlC2T"]["A"])
                self.assign_value(self.reads_over_conversion_control_C2T_C,data["BaseCounts"]["OverConversionControlC2T"]["C"])
                self.assign_value(self.reads_over_conversion_control_C2T_G,data["BaseCounts"]["OverConversionControlC2T"]["G"])
                self.assign_value(self.reads_over_conversion_control_C2T_T,data["BaseCounts"]["OverConversionControlC2T"]["T"])
                self.assign_value(self.reads_over_conversion_control_C2T_N,data["BaseCounts"]["OverConversionControlC2T"]["N"])
            #Base Counts OverConversionControl_G2A
            if "OverConversionControlG2A" in data["BaseCounts"]:
                self.assign_value(self.reads_over_conversion_control_G2A_A,data["BaseCounts"]["OverConversionControlG2A"]["A"])
                self.assign_value(self.reads_over_conversion_control_G2A_C,data["BaseCounts"]["OverConversionControlG2A"]["C"])
                self.assign_value(self.reads_over_conversion_control_G2A_G,data["BaseCounts"]["OverConversionControlG2A"]["G"])
                self.assign_value(self.reads_over_conversion_control_G2A_T,data["BaseCounts"]["OverConversionControlG2A"]["T"])
                self.assign_value(self.reads_over_conversion_control_G2A_N,data["BaseCounts"]["OverConversionControlG2A"]["N"])
            
            #Vector History Mapq (From Mapping Quality 0 to Mapping Quality 60)
            self.mapping_quality_reads = data["HistMapq"]
            
            #Read Length histogram
            self.read_length_histogram = data["HistReadLen"]
            
            if self.is_paired:
                 #Correct Pairs                 
                 self.correct_pairs = int (data["CorrectPairs"])
                 #Histrogram Insert Sizes            
                 self.read_insert_size_histogram = data["HistTemplateLen"]
                 
                 
class SampleStats(BsStats):
    """ Statistics per sample. Sum of all metrics per lane """

    def __init__(self,name=None,list_lane_stats=None):

        #Call Parent class constructor
        super(SampleStats, self).__init__()        
        
        # Global Stats from all set of lanes belonging to a given sample
        self.name = name
        
        #Sum of all statistics
        self.totalSampleUniqueReads = 0
        self.totalSampleReads = 0
        self.averageSampleUniqueReads = 0
        
        for lane_stats in list_lane_stats:
            self.sum_values(lane_stats)
            self.totalSampleUniqueReads += lane_stats.getUniqueMappedReads()
            self.totalSampleReads += lane_stats.getTotalMappedReads()
           
        self.averageSampleUniqueReads = float(self.totalSampleUniqueReads)/float(self.totalSampleReads) * 100
           
        #List of lanes
        self.list_lane_stats = list_lane_stats
        #Get First Lane Paired Status
        self.is_paired = list_lane_stats[0].is_paired

class RunBasicStats(object):
    """ Class responsable of basic functions """

    def __init__(self,mapping_stats=None,png_mapq_histogram=None,png_insert_size_histogram=None):
        self.mapping_stats = mapping_stats
        self.png_mapq_histogram = png_mapq_histogram
        self.png_insert_size_histogram = png_insert_size_histogram
        self.is_paired = True
        
    def createStatsForTableReport(self):
        """Create Stats for Table Report either HTML or Sphinx"""
        #Sequenced Reads
        self.mapping_stats.sequenced_reads.sumPair()
        self.mapping_stats.sequenced_reads.setPercentages(self.mapping_stats.sequenced_reads.total,self.mapping_stats.sequenced_reads.pair_one,self.mapping_stats.sequenced_reads.pair_two)
        #General reads 
        self.mapping_stats.reads.setPercentages(self.mapping_stats.sequenced_reads.total,self.mapping_stats.sequenced_reads.pair_one,self.mapping_stats.sequenced_reads.pair_two)
        #Reads in Control sequences
        self.mapping_stats.reads_sequencing_control.setPercentages(self.mapping_stats.sequenced_reads.total,self.mapping_stats.sequenced_reads.pair_one,self.mapping_stats.sequenced_reads.pair_two)       
        #Reads under conversion control  
        self.mapping_stats.reads_under_conversion_control.setPercentages(self.mapping_stats.sequenced_reads.total,self.mapping_stats.sequenced_reads.pair_one,self.mapping_stats.sequenced_reads.pair_two)
        #Reads over conversion control
        self.mapping_stats.reads_over_conversion_control.setPercentages(self.mapping_stats.sequenced_reads.total,self.mapping_stats.sequenced_reads.pair_one,self.mapping_stats.sequenced_reads.pair_two)
        #Unmapped reads
        self.mapping_stats.reads_unmapped.setPercentages(self.mapping_stats.sequenced_reads.total,self.mapping_stats.sequenced_reads.pair_one,self.mapping_stats.sequenced_reads.pair_two)
        #Bisulfite_reads
        self.mapping_stats.reads_reference_C2T.setPercentages(self.mapping_stats.sequenced_reads.total,self.mapping_stats.sequenced_reads.pair_one,self.mapping_stats.sequenced_reads.pair_two)
        self.mapping_stats.reads_reference_G2A.setPercentages(self.mapping_stats.sequenced_reads.total,self.mapping_stats.sequenced_reads.pair_one,self.mapping_stats.sequenced_reads.pair_two)

    def createBasesStats(self):
        """Create Bases Statistics"""  
                
        self.is_paired = self.mapping_stats.is_paired 

        #Base counts Overall
        base_counts_nucleotide_stats = NucleotideStats(a_bp=self.mapping_stats.reads_overall_A,c_bp=self.mapping_stats.reads_overall_C,\
                                                       t_bp=self.mapping_stats.reads_overall_T,g_bp=self.mapping_stats.reads_overall_G,n_bp=self.mapping_stats.reads_overall_N)
                                                       
        self.mapping_stats.reads_overall_A.setPercentagesNucleotides(base_counts_nucleotide_stats.all_a_percentage_bp,\
                                                                     base_counts_nucleotide_stats.pair_one_a_percentage_bp,base_counts_nucleotide_stats.pair_two_a_percentage_bp) 
        self.mapping_stats.reads_overall_C.setPercentagesNucleotides(base_counts_nucleotide_stats.all_c_percentage_bp,base_counts_nucleotide_stats.pair_one_c_percentage_bp,\
                                                                     base_counts_nucleotide_stats.pair_two_c_percentage_bp)       
        self.mapping_stats.reads_overall_G.setPercentagesNucleotides(base_counts_nucleotide_stats.all_g_percentage_bp,base_counts_nucleotide_stats.pair_one_g_percentage_bp,\
                                                                     base_counts_nucleotide_stats.pair_two_g_percentage_bp)   
        self.mapping_stats.reads_overall_T.setPercentagesNucleotides(base_counts_nucleotide_stats.all_t_percentage_bp,base_counts_nucleotide_stats.pair_one_t_percentage_bp,\
                                                                     base_counts_nucleotide_stats.pair_two_t_percentage_bp)  
        self.mapping_stats.reads_overall_N.setPercentagesNucleotides(base_counts_nucleotide_stats.all_n_percentage_bp,base_counts_nucleotide_stats.pair_one_n_percentage_bp,\
                                                                     base_counts_nucleotide_stats.pair_two_n_percentage_bp)   
                                                                     
        #Base counts GeneralC2T
        general_C2T_nucleotide_stats = NucleotideStats(a_bp=self.mapping_stats.reads_general_C2T_A,c_bp=self.mapping_stats.reads_general_C2T_C,t_bp=self.mapping_stats.reads_general_C2T_T,\
                                                       g_bp=self.mapping_stats.reads_general_C2T_G,n_bp=self.mapping_stats.reads_general_C2T_N)
                
        self.mapping_stats.reads_general_C2T_A.setPercentagesNucleotides(general_C2T_nucleotide_stats.all_a_percentage_bp,general_C2T_nucleotide_stats.pair_one_a_percentage_bp,\
                                                                         general_C2T_nucleotide_stats.pair_two_a_percentage_bp)
        self.mapping_stats.reads_general_C2T_C.setPercentagesNucleotides(general_C2T_nucleotide_stats.all_c_percentage_bp,general_C2T_nucleotide_stats.pair_one_c_percentage_bp,\
                                                                         general_C2T_nucleotide_stats.pair_two_c_percentage_bp)    
        self.mapping_stats.reads_general_C2T_G.setPercentagesNucleotides(general_C2T_nucleotide_stats.all_g_percentage_bp,general_C2T_nucleotide_stats.pair_one_g_percentage_bp,\
                                                                         general_C2T_nucleotide_stats.pair_two_g_percentage_bp)    
        self.mapping_stats.reads_general_C2T_T.setPercentagesNucleotides(general_C2T_nucleotide_stats.all_t_percentage_bp,general_C2T_nucleotide_stats.pair_one_t_percentage_bp,\
                                                                         general_C2T_nucleotide_stats.pair_two_t_percentage_bp)
        self.mapping_stats.reads_general_C2T_N.setPercentagesNucleotides(general_C2T_nucleotide_stats.all_n_percentage_bp,general_C2T_nucleotide_stats.pair_one_n_percentage_bp,\
                                                                         general_C2T_nucleotide_stats.pair_two_n_percentage_bp)                                                                         
                                                                     
        #Base counts GeneralG2A
        general_G2A_nucleotide_stats = NucleotideStats(a_bp=self.mapping_stats.reads_general_G2A_A,c_bp=self.mapping_stats.reads_general_G2A_C,t_bp=self.mapping_stats.reads_general_G2A_T,\
                                                       g_bp=self.mapping_stats.reads_general_G2A_G,n_bp=self.mapping_stats.reads_general_G2A_N)
        
        self.mapping_stats.reads_general_G2A_A.setPercentagesNucleotides(general_G2A_nucleotide_stats.all_a_percentage_bp,general_G2A_nucleotide_stats.pair_one_a_percentage_bp,\
                                                                         general_G2A_nucleotide_stats.pair_two_a_percentage_bp)
        self.mapping_stats.reads_general_G2A_C.setPercentagesNucleotides(general_G2A_nucleotide_stats.all_c_percentage_bp,general_G2A_nucleotide_stats.pair_one_c_percentage_bp,\
                                                                         general_G2A_nucleotide_stats.pair_two_c_percentage_bp)
        self.mapping_stats.reads_general_G2A_G.setPercentagesNucleotides(general_G2A_nucleotide_stats.all_g_percentage_bp,general_G2A_nucleotide_stats.pair_one_g_percentage_bp,\
                                                                         general_G2A_nucleotide_stats.pair_two_g_percentage_bp)     
        self.mapping_stats.reads_general_G2A_T.setPercentagesNucleotides(general_G2A_nucleotide_stats.all_t_percentage_bp,general_G2A_nucleotide_stats.pair_one_t_percentage_bp,\
                                                                         general_G2A_nucleotide_stats.pair_two_t_percentage_bp)                                                                         
        self.mapping_stats.reads_general_G2A_N.setPercentagesNucleotides(general_G2A_nucleotide_stats.all_n_percentage_bp,general_G2A_nucleotide_stats.pair_one_n_percentage_bp,\
                                                                         general_G2A_nucleotide_stats.pair_two_n_percentage_bp)           
                                                                         
        #Base Counts UnderConversionControlC2T
        under_conversion_control_C2T_stats = NucleotideStats(a_bp=self.mapping_stats.reads_under_conversion_control_C2T_A,c_bp=self.mapping_stats.reads_under_conversion_control_C2T_C,\
                                                             t_bp=self.mapping_stats.reads_under_conversion_control_C2T_T,g_bp=self.mapping_stats.reads_under_conversion_control_C2T_G,\
                                                             n_bp=self.mapping_stats.reads_under_conversion_control_C2T_N)
                
        self.mapping_stats.reads_under_conversion_control_C2T_A.setPercentagesNucleotides(under_conversion_control_C2T_stats.all_a_percentage_bp,under_conversion_control_C2T_stats.pair_one_a_percentage_bp,\
                                                                                          under_conversion_control_C2T_stats.pair_two_a_percentage_bp)
        self.mapping_stats.reads_under_conversion_control_C2T_C.setPercentagesNucleotides(under_conversion_control_C2T_stats.all_c_percentage_bp,under_conversion_control_C2T_stats.pair_one_c_percentage_bp,\
                                                                                          under_conversion_control_C2T_stats.pair_two_c_percentage_bp)
        self.mapping_stats.reads_under_conversion_control_C2T_G.setPercentagesNucleotides(under_conversion_control_C2T_stats.all_g_percentage_bp,under_conversion_control_C2T_stats.pair_one_g_percentage_bp,\
                                                                                          under_conversion_control_C2T_stats.pair_two_g_percentage_bp)                                                                                          
        self.mapping_stats.reads_under_conversion_control_C2T_T.setPercentagesNucleotides(under_conversion_control_C2T_stats.all_t_percentage_bp,under_conversion_control_C2T_stats.pair_one_t_percentage_bp,\
                                                                                          under_conversion_control_C2T_stats.pair_two_t_percentage_bp)
        self.mapping_stats.reads_under_conversion_control_C2T_N.setPercentagesNucleotides(under_conversion_control_C2T_stats.all_n_percentage_bp,under_conversion_control_C2T_stats.pair_one_n_percentage_bp,\
                                                                                          under_conversion_control_C2T_stats.pair_two_n_percentage_bp)  
                                                                                          
        #Base Counts UnderConversionControlG2A
        under_conversion_control_G2A_stats = NucleotideStats(a_bp=self.mapping_stats.reads_under_conversion_control_G2A_A,c_bp=self.mapping_stats.reads_under_conversion_control_G2A_C,\
                                                             t_bp=self.mapping_stats.reads_under_conversion_control_G2A_T,g_bp=self.mapping_stats.reads_under_conversion_control_G2A_G,\
                                                             n_bp=self.mapping_stats.reads_under_conversion_control_G2A_N)
                
        self.mapping_stats.reads_under_conversion_control_G2A_A.setPercentagesNucleotides(under_conversion_control_G2A_stats.all_a_percentage_bp,under_conversion_control_G2A_stats.pair_one_a_percentage_bp,\
                                                                                          under_conversion_control_G2A_stats.pair_two_a_percentage_bp)   
        self.mapping_stats.reads_under_conversion_control_G2A_C.setPercentagesNucleotides(under_conversion_control_G2A_stats.all_c_percentage_bp,under_conversion_control_G2A_stats.pair_one_c_percentage_bp,\
                                                                                          under_conversion_control_G2A_stats.pair_two_c_percentage_bp)    
        self.mapping_stats.reads_under_conversion_control_G2A_G.setPercentagesNucleotides(under_conversion_control_G2A_stats.all_g_percentage_bp,under_conversion_control_G2A_stats.pair_one_g_percentage_bp,\
                                                                                          under_conversion_control_G2A_stats.pair_two_g_percentage_bp) 
        self.mapping_stats.reads_under_conversion_control_G2A_T.setPercentagesNucleotides(under_conversion_control_G2A_stats.all_t_percentage_bp,under_conversion_control_G2A_stats.pair_one_t_percentage_bp,\
                                                                                          under_conversion_control_G2A_stats.pair_two_t_percentage_bp)
        self.mapping_stats.reads_under_conversion_control_G2A_N.setPercentagesNucleotides(under_conversion_control_G2A_stats.all_n_percentage_bp,under_conversion_control_G2A_stats.pair_one_n_percentage_bp,\
                                                                                          under_conversion_control_G2A_stats.pair_two_n_percentage_bp)
                                                                                          
        #Base Counts OverConversionControl_C2T
        over_conversion_control_C2T_stats = NucleotideStats(a_bp=self.mapping_stats.reads_over_conversion_control_C2T_A,c_bp=self.mapping_stats.reads_over_conversion_control_C2T_C,\
                                                       t_bp=self.mapping_stats.reads_over_conversion_control_C2T_T,g_bp=self.mapping_stats.reads_over_conversion_control_C2T_G,\
                                                       n_bp=self.mapping_stats.reads_over_conversion_control_C2T_N)
                
        self.mapping_stats.reads_over_conversion_control_C2T_A.setPercentagesNucleotides(over_conversion_control_C2T_stats.all_a_percentage_bp,over_conversion_control_C2T_stats.pair_one_a_percentage_bp,\
                                                                                         over_conversion_control_C2T_stats.pair_two_a_percentage_bp)        
        self.mapping_stats.reads_over_conversion_control_C2T_C.setPercentagesNucleotides(over_conversion_control_C2T_stats.all_c_percentage_bp,over_conversion_control_C2T_stats.pair_one_c_percentage_bp,\
                                                                                         over_conversion_control_C2T_stats.pair_two_c_percentage_bp)    
        self.mapping_stats.reads_over_conversion_control_C2T_G.setPercentagesNucleotides(over_conversion_control_C2T_stats.all_g_percentage_bp,over_conversion_control_C2T_stats.pair_one_g_percentage_bp,\
                                                                                         over_conversion_control_C2T_stats.pair_two_g_percentage_bp) 
        self.mapping_stats.reads_over_conversion_control_C2T_T.setPercentagesNucleotides(over_conversion_control_C2T_stats.all_t_percentage_bp,over_conversion_control_C2T_stats.pair_one_t_percentage_bp,\
                                                                                         over_conversion_control_C2T_stats.pair_two_t_percentage_bp)   
        self.mapping_stats.reads_over_conversion_control_C2T_N.setPercentagesNucleotides(over_conversion_control_C2T_stats.all_n_percentage_bp,over_conversion_control_C2T_stats.pair_one_n_percentage_bp,\
                                                                                         over_conversion_control_C2T_stats.pair_two_n_percentage_bp) 
                                                                                         
        #Base Counts OverConversionControl_G2A
        over_conversion_control_G2A_stats = NucleotideStats(a_bp=self.mapping_stats.reads_over_conversion_control_G2A_A,c_bp=self.mapping_stats.reads_over_conversion_control_G2A_C,\
                                                       t_bp=self.mapping_stats.reads_over_conversion_control_G2A_T,g_bp=self.mapping_stats.reads_over_conversion_control_G2A_G,\
                                                       n_bp=self.mapping_stats.reads_over_conversion_control_G2A_N)
        
        self.mapping_stats.reads_over_conversion_control_G2A_A.setPercentagesNucleotides(over_conversion_control_G2A_stats.all_a_percentage_bp,over_conversion_control_G2A_stats.pair_one_a_percentage_bp,\
                                                                                         over_conversion_control_G2A_stats.pair_two_a_percentage_bp)
        self.mapping_stats.reads_over_conversion_control_G2A_C.setPercentagesNucleotides(over_conversion_control_G2A_stats.all_c_percentage_bp,over_conversion_control_G2A_stats.pair_one_c_percentage_bp,\
                                                                                         over_conversion_control_G2A_stats.pair_two_c_percentage_bp)
        self.mapping_stats.reads_over_conversion_control_G2A_G.setPercentagesNucleotides(over_conversion_control_G2A_stats.all_g_percentage_bp,over_conversion_control_G2A_stats.pair_one_g_percentage_bp,\
                                                                                         over_conversion_control_G2A_stats.pair_two_g_percentage_bp)
        self.mapping_stats.reads_over_conversion_control_G2A_T.setPercentagesNucleotides(over_conversion_control_G2A_stats.all_t_percentage_bp,over_conversion_control_G2A_stats.pair_one_t_percentage_bp,\
                                                                                         over_conversion_control_G2A_stats.pair_two_t_percentage_bp)
        self.mapping_stats.reads_over_conversion_control_G2A_N.setPercentagesNucleotides(over_conversion_control_G2A_stats.all_n_percentage_bp,over_conversion_control_G2A_stats.pair_one_n_percentage_bp,\
                                                                                         over_conversion_control_G2A_stats.pair_two_n_percentage_bp)

    @staticmethod
    def saveDocument(file_name=None,vectorContent=None):
        '''Save the document to a given name file
        
        file_name - File path to store the document
        vectorContent - Vector of contents to write in a file.
        '''
        with open(file_name, 'w') as fileDocument:
            for line in vectorContent:
                fileDocument.write(line)
                                              
    def drawMapqHistogram(self):
        """ From matplot lib plots a Mappping qualty histogram
        """
        #1. PLOT BUILDING
        readsMapq = self.mapping_stats.mapping_quality_reads
       
        mapqList = []
        for mapq in range(61):
            mapqList.append(mapq)
     
        matplotlib.pyplot.ioff()
        figure = plt.figure()
        plt.bar(mapqList,readsMapq,width=1,align='center',facecolor='blue', alpha=0.75)
        plt.xlabel('MapQ')
        plt.ylabel('Fragments')
        plt.title('MapQ Histogram')
        plt.axis([0, 60,min(readsMapq), max(readsMapq)])
        plt.grid(True)
        
        pylab.savefig(self.png_mapq_histogram)
        
        plt.close(figure)
        
        
    def drawInsertSizePlot(self):
        """ From matplot lib plots a Insert Size Plot
        """
        iSizeList = []
        readsList = []
        
        #1. PLOT BUILDING
        histogram_template_len = self.mapping_stats.read_insert_size_histogram
        
        for insert_size_length,reads in histogram_template_len.iteritems():
            iSizeList.append(int(insert_size_length))
            readsList.append(int(reads))
            
        matplotlib.pyplot.ioff()            
        figure = plt.figure()

        plt.plot(iSizeList, readsList, '.',color="r")
        plt.xlabel('Insert Size (bp)')
        plt.ylabel('Reads')
        plt.title('Insert Size Histogram')
        plt.axis([min(iSizeList), 800,min(readsList), max(readsList)])
        plt.grid(True)
        
        pylab.savefig(self.png_insert_size_histogram)
        
        plt.close(figure)
                           
    def drawMultipleInsertSizePlot(self):
        """ From matplot lib plots a Insert Size Plot
        """ 
        iSizeList = []
        readsList = []
        readsYaxis = []
        sizeXaxis = []
        lanesNames = []
        
        #1. PLOT BUILDING
        for laneStats in self.mapping_stats.list_lane_stats:
            lanesNames.append(laneStats.name)
            histogram_template_len = laneStats.read_insert_size_histogram
            sizeList  = []
            readList = []
            for insert_size_length,reads in histogram_template_len.iteritems():
                sizeList.append(int(insert_size_length))
                readList.append(int(reads))
                
            readsYaxis.extend(readList)
            sizeXaxis.extend(sizeList)
            iSizeList.append(sizeList)
            readsList.append(readList)
            
        matplotlib.pyplot.ioff()            
        figure = plt.figure()

        for iSize, readList in zip(iSizeList, readsList):        
            plt.plot(iSize, readList,'.')
        
        plt.xlabel('Insert Size (bp)')
        plt.ylabel('Reads')
        plt.title('Insert Size Histogram Per Lane')
        plt.axis([min(sizeXaxis), 800,min(readsYaxis), max(readsYaxis)])
        plt.legend(lanesNames,loc='upper right')        
        
        plt.grid(True)
       
        pylab.savefig(self.png_insert_size_histogram)
        
        plt.close(figure)                                         
                                                                     
    def drawMultipleMapqHistogram(self):
        """ From matplot lib plots a Mappping qualty histogram
        """  
        mapqFragmentsLanes = []
              
        for laneStats in self.mapping_stats.list_lane_stats:
            mapqFragmentsLanes.append(laneStats.mapping_quality_reads)
    
        matplotlib.pyplot.ioff()
        figure = plt.figure()
        ax = figure.add_subplot(111,projection='3d')
                       
        lane = 0
        qualityRange = list(range(61))        
        for fragmentsQuality in mapqFragmentsLanes:
            ax.bar(qualityRange,fragmentsQuality, zs=lane, zdir='y', color='b', alpha=0.8)
            lane = lane + 1
            
        ax.set_xlabel('MapQ')
        ax.set_ylabel('Lanes')
        ax.set_zlabel('Fragments')        
                
        #http://people.duke.edu/~ccc14/pcfb/numpympl/MatplotlibBarPlots.html
        pylab.savefig(self.png_mapq_histogram)
        
        plt.close(figure)                                                                  
                                                                     
                                                                     
                                                                     
                                                                     
                                                                     
                                                                     
                                                                     
                                                                     
                                                                     
                                                                     
                                                                     
                                                                     
                                                                     
                                                                     
                                                                     
                                                                     
                                                                     
                                                                     
                                                                     
                                                                     
                                                                     
                                                                     
                                                       
    