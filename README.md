Scripts assocaited with Gacita et al. "Integrative epigenomic analysis identifies enhancer modifying variants linked to cardiomyopathy genes" 2020. 


Files in directory:
1. Find_Regulatory_Variants_Full_Pipeline- python script for creating shell scripts to excute the Find Regulatory Variants pipeline. 
2. Find_Regulatory_Variants_Dependencies- python dependency script thats needed for the Find Regualtory Variants pipeline.
3. Annotation_Datasets- directory containing annotation data derived from IPSC-CMs for use with the Find Regulatory Variants pipeline.
	-AllGenes_1000_intersected_2min.bed- Promoter enhancer interactions from promoter capture Hi-C data in IPSC-CMs. 
	-ATAC_IPSC_CM_superset.bed- ATAC-seq peaks from IPSC-CMs.
	-GATA4_IPSC_CM.bed- GATA4 ChIP-seq peaks from IPSC-CMs.
	-TBX5_IPSC_CM.bed- TBX5 ChIP-seq peaks from IPSC-CMs.
	-GATA4_homer.bed- Genome wide locations of the GATA4 binding motif downloaded from HOMER.
	-TBX5_homer.bed- Genome wide locations of the TBX5 binding motif downloaded from HOMER.
	-GATA4.motif- Position weight matrix of the GATA4 motif downloaded from HOMER.
	-TBX5.motif- Position weight matrix of the TBX5 motif downloaded from HOMER



############################################################
############################################################
############ Find Regulatory Variants Script ###############
############################################################
################### Anthony Gacita #########################
##################### McNally Lab ##########################

Details:
This script takes in VCF files and uses intersection with relevant functional data to find variants that have a high probability of regulating gene expression.  It is currently set up to run on Northwestern's Quest computing cluster, but can be modified to run anywhere. 
 
Set-UP Instructions
1. Gather VCF File(s). Must end in .vcf and must use 'chr' naming convention in the header and variant entries.
2. Gather relevant functional annotations in .bed format (or use default IPSC-CM ones).
3. Upload Find_Regulatory_Variants_Full_Pipline.py, annotation file, and VCF files to Quest.
4. Run python script to generate shell (.sh) scripts. 
	Example Command:
	python Find_Regulatory_Variants_Full_Pipeline.py -w /Directory/to/VCF_File -o /Directory/to/Output/Location -s /Directory/to/Shell/Scripts/Location
5. Submit .sh scripts to Quest. 

Required Programs:
GATK (v.3.6.0)
picard (v.2.8.0)
bedtools
bedops
homer
vcftools


NOTE:
See the python script to:
	- Update PATHs of annotation files
	- Update PATH of depedency script
	- Update PATHs of GATK and picard jar files
