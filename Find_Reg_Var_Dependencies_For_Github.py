#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 16:09:12 2018

@author: amg171
"""


########################################################################################################################################################################################################
#All the dependency options
from optparse import OptionParser
parser=OptionParser()

parser.add_option("--protocol",dest="protocol",default="", help="Which dependency script do you want to run? Options : fixBed, combineResults, CalculateAF, CreateAFDictionary, AnnotateVariants, CreateRegionsFile, CreateGenotypeDataframe") 

#fixBed
parser.add_option("--filefixbed", dest="filefixbed", default="", help= "Input file")
parser.add_option("--outputfixbed", dest="outputfixbed", default="./", help= "Output")


(options,args) = parser.parse_args()
########################################################################################################################################################################################################
#Sanity Checks
if options.protocol not in ['fixBed','combineResults','CalculateAF','CreateAFDictionary','AnnotateVariants','CreateRegionsFile','CreateGenotypeDataframe']:
    print('Protocol option not recognized. Valid options are: fixBed, combineResults, CalculateAF, CreateAFDictionary, AnnotateVariants, CreateRegionsFile, CreateGenotypeDataframe.')
    

########################################################################################################################################################################################################
#fix bed for homer

if options.protocol == 'fixBed':
    with open(options.filefixbed,'r') as ins:
        with open(options.outputfixbed,'w') as f:
            for line in ins:
                items=line.split()
                chrom=items[1]
                new_chrom=chrom[chrom.find('chr'):chrom.find(':')]
                f.write(new_chrom+'\t'+items[2]+'\t'+items[3]+'\t'+items[4]+'\t'+items[5]+'\t'+items[6]+'\n')
        f.close()
    ins.close()


########################################################################################################################################################################################################
