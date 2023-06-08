#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 14:27:53 2023

@authors: jasonsahl, grantpemberton
"""

from __future__ import division
import sys
import optparse
from optparse import OptionParser
import os
import collections
import time

try:
    from Bio import SeqIO
except:
    print("Biopython is not installed but needs to be...exiting")
    sys.exit()

#logPrint stuff
OUTSTREAM = sys.stdout
ERRSTREAM = sys.stderr
DEBUG = False

def logPrint(msg, stream=None):
    if stream is None:
        stream = OUTSTREAM
    stream.write('LOG: %s - %s\n' % (timestamp(), removeRecursiveMsg(msg)))
    stream.flush()

def errorPrint(msg, stream=None):
    if stream is None:
        stream = ERRSTREAM

    stream.write('ERROR: %s - %s\n' % (timestamp(), removeRecursiveMsg(msg)))
    stream.flush()

def debugPrint(fmsg, stream=None):
    """In this case msg is a function, so the work is only done if debugging is one"""
    if DEBUG:
        if stream is None:
            stream = ERRSTREAM

        stream.write('DEBUG: %s - %s\n' % (timestamp(), removeRecursiveMsg(fmsg())))
        stream.flush()

def timestamp():
    return time.strftime('%Y/%m/%d %H:%M:%S')

def removeRecursiveMsg(msg):
    """
    This takes a message and if it starts with something that looks like
    a message generated with these tools it chops it off.  Useful if using
    one of these logging functions to print output from a program using
    the same logging functions
    """
    if msg.startswith('ERROR: ') or msg.startswith('DEBUG: ') or msg.startswith('LOG: '):
        return msg.split(' - ', 1)[1]
    else:
        return msg

#This function tests that required files are actually present
def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s cannot be opened' % option)
        sys.exit()

def parse_fasta_by_coverage(in_fasta,min_cov):
    passing_records = []
    failed_records = []
    with open(in_fasta) as my_fasta:
        for record in SeqIO.parse(my_fasta,"fasta"):
            #The coverage will be the 3rd field
            header_fields = record.id.split("_")
            if int(header_fields[2])<min_cov:
                failed_records.append(record.id)
            else:
                passing_records.append(record.id)
    if len(failed_records)>0:
        print("%s records were below %sX and will be filtered" % (len(failed_records),str(min_cov)))
    return passing_records

def parse_zygosity(in_fasta,passing_records,proportion):
    passing = []
    failed = []
    sample_dict = {}
    #Parse through the FASTA again and create a dictionary
    with open(in_fasta) as my_fasta:
        for record in SeqIO.parse(my_fasta,"fasta"):
            #The sample name will be in the first field
            header_fields = record.id.split("_")
            #Only look at sequences that have passed the first filter
            if record.id in passing_records:
                #This will create a dictionary of sample: all_associated_coverages
                try:
                    sample_dict[header_fields[0]+"_"+header_fields[1]].append(int(header_fields[2]))
                except KeyError:
                    sample_dict[header_fields[0]+"_"+header_fields[1]] = [int(header_fields[2])]
    #Iterate through the dictionary
    for k,v in sample_dict.items():
        values = sorted(v,reverse=True)
        #I assume this is unlikely, but if a sample is only present once, keep it
        if len(values) == 1:
            passing.append(k+"_"+str(v))
        else:
            #This only keeps the two largest values in a list
            kept_values = values[:2]
            #Here I divide the large value by the small value
            if float(kept_values[0]/kept_values[1])>=proportion:
                #I want to keep the entire sample name for subsequent filtering
                passing.append(k+"_"+str(kept_values[0]))
                #Here I will also keep the second read as it has passed the proportion
                passing.append(k+"_"+str(kept_values[1]))
            #Now we need to analyze the case where the proportion filter is failed
            elif float(kept_values[0]/kept_values[1])<proportion:
                #In this case, I just want to keep the most common one
                passing.append(k+"_"+str(kept_values[0]))
    #I will now find differences between the two lists and report how many were filtered
    diffs = set(passing_records).difference(set(passing))
    if len(diffs)>0:
        logPrint("%s samples failed the proportion filter and will be removed" % len(diffs))
    return passing
       

def get_alleles(in_fasta):
    sequences = []
    with open(in_fasta) as my_fasta:
        for record in SeqIO.parse(my_fasta,"fasta"):
            sequences.append(str(record.seq))
    #This creates a dictionary of sequence and count_number
    frequency = collections.Counter(sequences)
    #Sort the dictionary by value in decending order: most common element is listed first
    sorted_frequency = sorted(frequency.items(), key=lambda x:x[1], reverse=True)
    #returns this dictionary to be used later
    return sorted_frequency

def assign_alleles(fasta,allele_file,allele_list,passing_records):
    #This covers the situation where there are no previous alleles
    #This list will include the previous alleles where necessary
    previous_allele_numbers = []
    file_name = os.path.basename(fasta).strip(".fasta")
    if "NULL" in allele_file:
        #We will need this dictionary for renaming the alleles
        allele_dict = {}
        #First, get the name of the allele, remove fasta extension
        #Now I will write the alleles in terms of order
        alleles_out = open("%s_alleles.tsv" % file_name, "w")
        counter = (i for i in range(len(allele_list)))
        #The second counter is needed just for the dictionary
        second_counter = (i for i in range(len(allele_list)))
        for item in allele_list:
            alleles_out.write(item[0]+"\t"+str(next(counter)+1)+"\n")
            allele_dict.update({item[0]:str(next(second_counter)+1)})
            #alleles_out.write(item[0]+"\t"+str(autoIncrement())+"\n")
        #close the file
        alleles_out.close()
    else:
        #This is the case where you already have allele numbers assigned
        allele_dict = {}
        with open(allele_file) as my_known_alleles:
            for line in my_known_alleles:
                fields = line.split()
                allele_dict.update({fields[0]:fields[1]})
                previous_allele_numbers.append(fields[1])
    #Now I will rename the FASTA file to include the information that I need
    with open(fasta) as my_fasta:
        genotype_out = open("%s.genotyped.fasta" % file_name,"w")
        for record in SeqIO.parse(my_fasta,"fasta"):
            #Check to make sure that I want to process this one
            if record.id in passing_records:
                #split the name for renaming purposes
                name_fields = record.id.split("_")
                #Sample name = field[0], count = field[2]
                #Now I will parse through the allele list and find a match
                if record.seq in allele_dict:
                    #This will pull out the allele number from the dictionary
                    allele = allele_dict.get(record.seq)
                else:
                    #I do this so I know where to start renumbering
                    last_allele = previous_allele_numbers[-1]
                    allele = int(last_allele)+1
                    #I also need to update the allele list
                    previous_allele_numbers.append(allele)
                    #I need to update the dictionary to include the new allele
                    allele_dict.update({record.seq:str(allele)})
                    #This indicates that the allele has not been seen before
                    #I will need to overwrite the existing alleles file
                new_header = name_fields[0]+"_"+file_name+"_"+name_fields[2]+"_"+str(allele)
                #From above, the locus should be the file_name
                #Write the new file
                genotype_out.write(">"+str(new_header)+"\n"+str(record.seq)+"\n")
                #Close the file
        #I will overwrite the alleles file to make sure I capture any new diversity
        if "NULL" not in allele_file:
            newout = open("%s_alleles.tsv" % file_name, "w")
            for k,v in allele_dict.items():
                newout.write(str(k)+"\t"+str(v)+"\n")
            newout.close()
        genotype_out.close()
               
def main(fasta,min_cov,proportion,alleles):
    passing_records = parse_fasta_by_coverage(fasta,min_cov)
    passing_records2 = parse_zygosity(fasta,passing_records,proportion)
    allele_list = get_alleles(fasta)
    assign_alleles(fasta,alleles,allele_list,passing_records2)

if __name__ == "__main__":
    #TODO: bump version when significant changes are made
    parser = OptionParser(usage="usage: %prog [options]",version="%prog 0.0.2")
    #First step: make sure it works on a single file
    parser.add_option("-f","--fasta",dest="fasta",
                     help="input FASTA file to filter [REQUIRED]",
                     type="string",action="callback",callback=test_file)
    parser.add_option("-c","--min_cov",dest="min_cov",
                     help="filter under minimum coverage; defaults to 2",
                     type="int",action="store",default=2)
    parser.add_option("-p","--proportion",dest="proportion",
                     help="reads under this proportion will be filtered; defaults to 0.3",
                     type="float",action="store",default=0.3)
    parser.add_option("-a","--alleles",dest="alleles",
                     help="TSV file of seuqence'\tprevious_allele_number",
                     type="string",action="store",default="NULL")
    options, args = parser.parse_args()
    mandatories = ["fasta"]
    for m in mandatories:
        if not getattr(options,m,None):
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)
    main(options.fasta,options.min_cov,options.proportion,options.alleles)




