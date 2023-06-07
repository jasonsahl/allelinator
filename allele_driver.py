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

try:
    from Bio import SeqIO
except:
    print("Biopython is not installed but needs to be...exiting")
    sys.exit()

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
        print("%s samples failed the proportion filter and will be removed" % len(diffs))
    return passing
        
def main(fasta,min_cov,proportion):
    passing_records = parse_fasta_by_coverage(fasta,min_cov)
    passing_records2 = parse_zygosity(fasta,passing_records,proportion)

if __name__ == "__main__":
    #TODO: bump version when significant changes are made
    parser = OptionParser(usage="usage: %prog [options]",version="%prog 0.0.1")
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
    options, args = parser.parse_args()
    mandatories = ["fasta"]
    for m in mandatories:
        if not getattr(options,m,None):
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)
    main(options.fasta,options.min_cov,options.proportion)




