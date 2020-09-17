#!/usr/bin/env python2
# -*- coding: utf-8 -*-

############################################################
############################################################
############ Find Regulatory Variants Script ###############
############################################################
################### Anthony Gacita #########################
##################### McNally Lab ##########################


#Details:
# This script takes in VCF files and uses intersection with relevant functional data to find variants that have a high probability 
# of regulating gene expression.  It is currently set up to run on Northwestern's Quest computing cluster, but can be modified to run anywhere. 

#NOTE: Right now, the Find Regulatory Variants script can be run on any VCF or group of VCFs. 

#Set-UP Instructions
#1. Gather VCF File(s). Must end in .vcf
#2. Gather relevant functional annotations in .bed format (or use default IPSC-CM ones)
#3. Run Script to generate shell and python scripts


#Required Programs:
#GATK
#picard
#bedtools
#bedops
#homer
#vcftools


#IMPORTANT NOTE ABOUT CHR ANNOTATION
#contigs must use chr naming convention in the VCF header & variant entries 


#Running Instructions
# 1. Upload this script into the same directory as your VCF files (unzipped)
# 2. Set the run specific prameters below
# 3. Running the python script will Create an .sh Script which should be submitted to Quest.


#Version History
# Version 2 Changes:
    # - Added command line options

# Ideas for next version
    # - combine all prints into an output file
    
    
###########################
####Command Line Options###   
########################### 
from optparse import OptionParser
parser=OptionParser()

parser.add_option("-w","--Working_Directory", dest="Working_Directory", default="", help="Working_Directory- Directory where the VCF files and this script is located. No last /")
parser.add_option("-o","--Output_Directory", dest="Output_Directory", default="", help="Output_Directory- Directory where you want the output files to be placed. No last / Default is Working_Directory/Output")
parser.add_option("-s","--Scripts_Directory", dest="Scripts_Directory", default="", help="Scripts_Directory- Directory where you want the scripts to be placed. No last / Default is Working_Directory/Scripts")

(options,args) = parser.parse_args()

if options.Output_Directory == '':
    options.Output_Directory= options.Working_Directory+'/Output'

if options.Scripts_Directory =='':
    options.Scripts_Directory=options.Working_Directory+'/Scripts'

###Constant Set-Up Parameters (DONT change these unless you know why)
#Tools
GATK='/home/amg171/McNally/McNally_TOOLS/GenomeAnalysisTK.jar' #full path to the GATK .jar
picard = '/home/amg171/McNally/McNally_TOOLS/picard.jar' #full path to the picard .jar

#This script has a dependency script 
dependency_script='/home/amg171/McNally/McNally_TOOLS/Find_Reg_var_Dependencies.py' #This script requires a dependency script. Put it's full path here. 


#Find Regulatory Variants Functional Datasets
Reference='/projects/genomicsshare/rna-seq_reference/McNally_Reference/GATK_Bundle/ucsc.hg19.fasta'
ATAC_Seq='/projects/genomicsshare/Enhancer_var_prediction/IPSC_CM_data/ATAC_Data/ATAC_IPSC_CM_superset.bed'
GATA4_exp='/projects/genomicsshare/Enhancer_var_prediction/IPSC_CM_data/Exp_TF_Data/GATA4_IPSC_CM.bed'
GATA4_comp='/projects/genomicsshare/Enhancer_var_prediction/IPSC_CM_data/Comp_TF_Data/GATA4_homer.bed'
TBX5_exp='/projects/genomicsshare/Enhancer_var_prediction/IPSC_CM_data/Exp_TF_Data/TBX5_IPSC_CM.bed'
TBX5_comp='/projects/genomicsshare/Enhancer_var_prediction/IPSC_CM_data/Comp_TF_Data/TBX5_homer.bed'
Hi_C= '/projects/genomicsshare/Enhancer_var_prediction/IPSC_CM_data/Hi-C_Data/AllGenes_1000_intersected_2min.bed'
GATA4='/home/amg171/McNally/GATA4.motif'
TBX5='/home/amg171/McNally/TBX5.motif'


#End Setup
################################################################################################################
################################################################################################################
################################################################################################################
import os

listdir = os.listdir(options.Working_Directory)

if not os.path.exists(options.Scripts_Directory):   
    os.makedirs(options.Scripts_Directory)

if not os.path.exists(options.Output_Directory):
    os.makedirs(options.Output_Directory)    

if not os.path.exists(options.Output_Directory+'/logs'):
    os.makedirs(options.Output_Directory+'/logs')  

instructions=open(options.Working_Directory+'/'+'instructions.txt','w')
instructions.write('This file contains the instructions for running the Find Regulatory Variants scripts.\n')
    
########################################
###### Find Regulatory Variants ########
########################################

for file in listdir:
    if file.endswith(".vcf"):
        script=open(options.Scripts_Directory+'/'+file[0:len(file)-4]+'_Find_Regulatory_Var.sh','w+')
        script.write('''#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH --mail-type=BEGIN,END,NONE,FAIL,REQUEUE
#SBATCH --output='''+options.Output_Directory+'''/logs/'''+file[0:len(file)-4]+'''_Find_Regulatory_Var.log
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --mem 80000
module load bedtools
module load bedops
module load python/anaconda3.6
module load java
module load homer
module load vcftools

''') 
        file_output= options.Output_Directory+'/'+file[0:len(file)-4]
        script.write('mkdir '+file_output+'\n')
        
        #Find new TF sites made by variants
        ##Create Alternative Reference file
        script.write('vcftools --vcf '+options.Working_Directory+'/'+file+' --remove-indels --recode --recode-INFO-all --out '+file_output+'/'+file[0:len(file)-4]+'_NO_INDEL\n')
        script.write('java -jar '+picard+' SortVcf I='+file_output+'/'+file[0:len(file)-4]+'_NO_INDEL.recode.vcf O='+file_output+'/'+file[0:len(file)-4]+'_sorted.vcf R=='+Reference+' MAX_RECORDS_IN_RAM=200000\n')
        script.write('java -jar '+GATK+' -T FastaAlternateReferenceMaker -U ALLOW_SEQ_DICT_INCOMPATIBILITY -R '+Reference+' -V '+file_output+'/'+file[0:len(file)-4]+'_sorted.vcf -o '+file_output+'/'+file[0:len(file)-4]+'_alt_reference.fasta\n')
        
        ##Scan the Alternative File for New TF Binding Sites
        script.write('scanMotifGenomeWide.pl '+GATA4+' '+file_output+'/'+file[0:len(file)-4]+'_alt_reference.fasta -bed > '+file_output+'/'+file[0:len(file)-4]+'_alt_GATA4.bed\n')
        script.write('python '+dependency_script+' --protocol fixBed --filefixbed '+file_output+'/'+file[0:len(file)-4]+'_alt_GATA4.bed --outputfixbed '+file_output+'/'+file[0:len(file)-4]+'_alt_GATA4_fixed.bed\n')
        script.write('scanMotifGenomeWide.pl '+TBX5+' '+file_output+'/'+file[0:len(file)-4]+'_alt_reference.fasta -bed > '+file_output+'/'+file[0:len(file)-4]+'_alt_TBX5.bed\n')
        script.write('python '+dependency_script+' --protocol fixBed --filefixbed '+file_output+'/'+file[0:len(file)-4]+'_alt_TBX5.bed --outputfixbed '+file_output+'/'+file[0:len(file)-4]+'_alt_TBX5_fixed.bed\n')
        
        ##Compare Alternative TF Binding Sites to Reference binding sites
        script.write('bedtools intersect -a '+file_output+'/'+file[0:len(file)-4]+'_alt_GATA4_fixed.bed -b '+GATA4_comp+' -v > '+file_output+'/'+file[0:len(file)-4]+'_alt_new_GATA4.bed\n')
        script.write('bedtools intersect -a '+file_output+'/'+file[0:len(file)-4]+'_alt_TBX5_fixed.bed -b '+TBX5_comp+' -v > '+file_output+'/'+file[0:len(file)-4]+'_alt_new_TBX5.bed\n')
        script.write('echo found alternative TF sites\n')
        
        ##Convert VCF to Bed
        script.write('vcf2bed <'+options.Working_Directory+'/'+file+'> '+file_output+'/'+file[0:len(file)-4]+'.bed\n')
        script.write('echo VCF converted to bed\n')
     
        #ATAC Overlap
        script.write('bedtools intersect -a '+file_output+'/'+file[0:len(file)-4]+'.bed -b '+ATAC_Seq+' -u > '+file_output+'/'+file[0:len(file)-4]+'_ATAC.bed\n')
        script.write('bedtools intersect -a '+file_output+'/'+file[0:len(file)-4]+'.bed -b '+ATAC_Seq+' -v > '+file_output+'/'+file[0:len(file)-4]+'_noATAC.bed\n')
        script.write('echo ATAC Overlap Completed\n')
        
        ##Experimental TF Overlap
        script.write('bedtools intersect -a '+file_output+'/'+file[0:len(file)-4]+'_ATAC.bed -b '+TBX5_exp+' -u > '+file_output+'/'+file[0:len(file)-4]+'_ATAC_expTBX5.bed\n')
        script.write('bedtools intersect -a '+file_output+'/'+file[0:len(file)-4]+'_ATAC.bed -b '+TBX5_exp+' -v > '+file_output+'/'+file[0:len(file)-4]+'_ATAC_noexpTBX5.bed\n')
        
        script.write('bedtools intersect -a '+file_output+'/'+file[0:len(file)-4]+'_ATAC.bed -b '+GATA4_exp+' -u > '+file_output+'/'+file[0:len(file)-4]+'_ATAC_expGATA4.bed\n')
        script.write('bedtools intersect -a '+file_output+'/'+file[0:len(file)-4]+'_ATAC.bed -b '+GATA4_exp+' -v > '+file_output+'/'+file[0:len(file)-4]+'_ATAC_noexpGATA4.bed\n')
        script.write('echo Experimental TF Overlap Completed\n')
        
        ##Comp TF Overlap
        script.write('bedtools intersect -a '+file_output+'/'+file[0:len(file)-4]+'_ATAC_expTBX5.bed -b '+TBX5_comp+' '+file_output+'/'+file[0:len(file)-4]+'_alt_new_TBX5.bed -u > '+file_output+'/'+file[0:len(file)-4]+'_ATAC_expTBX5_compTBX5.bed\n')
        script.write('bedtools intersect -a '+file_output+'/'+file[0:len(file)-4]+'_ATAC_expTBX5.bed -b '+TBX5_comp+' -v > '+file_output+'/'+file[0:len(file)-4]+'_ATAC_expTBX5_nocompTBX5.bed\n')
        
        script.write('bedtools intersect -a '+file_output+'/'+file[0:len(file)-4]+'_ATAC_expGATA4.bed -b '+GATA4_comp+' '+file_output+'/'+file[0:len(file)-4]+'_alt_new_GATA4.bed -u > '+file_output+'/'+file[0:len(file)-4]+'_ATAC_expGATA4_compGATA4.bed\n')
        script.write('bedtools intersect -a '+file_output+'/'+file[0:len(file)-4]+'_ATAC_expGATA4.bed -b '+GATA4_comp+' -v > '+file_output+'/'+file[0:len(file)-4]+'_ATAC_expGATA4_nocompGATA4.bed\n')
        script.write('echo Computational TF Overlap Completed\n')

        ##Merge the TF Datasets
        script.write('cat '+file_output+'/'+file[0:len(file)-4]+'_ATAC_expGATA4_compGATA4.bed '+file_output+'/'+file[0:len(file)-4]+'_ATAC_expTBX5_compTBX5.bed > '+file_output+'/'+file[0:len(file)-4]+'_ATAC_exp_comp_anyTF.bed\n')
        script.write('bedtools sort -i '+file_output+'/'+file[0:len(file)-4]+'_ATAC_exp_comp_anyTF.bed > '+file_output+'/'+file[0:len(file)-4]+'_ATAC_exp_comp_anyTF_sorted.bed\n')
        script.write('uniq '+file_output+'/'+file[0:len(file)-4]+'_ATAC_exp_comp_anyTF_sorted.bed '+file_output+'/'+file[0:len(file)-4]+'_ATAC_exp_comp_anyTF_nodups.bed\n')
        
        ##Overlap with Hi-C and Report Overlap
        script.write('bedtools intersect -a '+file_output+'/'+file[0:len(file)-4]+'_ATAC_exp_comp_anyTF_nodups.bed -b '+Hi_C+' -wa -wb > '+file_output+'/'+file[0:len(file)-4]+'_ATAC_exp_comp_anyTF_HiC.bed\n')
        script.write('bedtools intersect -a '+file_output+'/'+file[0:len(file)-4]+'_ATAC_exp_comp_anyTF_nodups.bed -b '+Hi_C+' -v > '+file_output+'/'+file[0:len(file)-4]+'_ATAC_exp_comp_anyTF_noHiC.bed\n')
        script.write('bedtools closest -a '+file_output+'/'+file[0:len(file)-4]+'_ATAC_exp_comp_anyTF_noHiC.bed -b '+Hi_C+' -d > '+file_output+'/'+file[0:len(file)-4]+'_ATAC_exp_comp_anyTF_noHiC_closest.bed\n')
        script.write('echo Hi-C TF Overlap Completed\n')

        ##Create Report of number of lines
        script.write('for f in '+file_output+'/*; do wc -l $f >> '+file_output+'/line_counts.txt; done\n')

        
print('Find Regulatory Variants script successfully created for '+file+'\n')
instructions.write('Step 1: find -name *_Find_Regulatory_Var.sh -exec sbatch {} \;\n')

print('Check the instructions.txt file for running instructions')



   
