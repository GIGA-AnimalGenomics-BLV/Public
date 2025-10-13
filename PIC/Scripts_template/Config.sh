#!/bin/bash

# Define Variable
#########################################################
### Versions of software avalaible in $SingularityDir ###
### Can be replaced				     ###
#########################################################
SingularityDir=<PATH/TO/SINGULARITIES>
export NXF_SINGULARITY_CACHEDIR=$SingularityDir
export APPTAINER_BIND="PATH/TO/BIND/TO/SINGULARITIES/"

# 1. Ownmade Tools/binaries path
ToolsDir=<PATH/TO/TOOLS/DIRECTORY>
concatenateFASTQ=${ToolsDir}/concatenateFastq.pl
resynchronizePairs=${ToolsDir}/resynchronizePaired_forpy3.py

# 2. Environment Definition
Project=			 # name of the project 

echo $Project
# PATH to directory containing references
bowtie2HostVirus="path/to/bowtie2/virusHost/index/prefix"
bowtie2LTR="path/to/bowtie2/LTR/index/prefix"
geneBedFile="path/to/geneInfos.sorted.txt"

# HTLV:
LTR5len=36			# Length of the LTR5 sequence, starting from the LTR to the end of the 5'LTR primer
LTR3len=45			# Length of the LTR3 sequence, starting from the 3'LTR primer start to the LTR end
virus="HTLV_ATK"		# Name given to the indexed viral sequence
linkLen=22			# Length of the linker sequence

# Define directories :
ClonalityDir=<PATH/TO/CLONALITY/DIRECTORY>				# Directory for Outputs
SampleSheetDir=<PATH/TO/SAMPLESHEET/DIRECTORY>				# Directory containing the samplesheet

# Define files :
files=${SampleSheetDir}/HLTV1_samplesheet.txt				# SampleSheet containing information for alignment
outdir=${ClonalityDir}/results/${Project}				# 

mkdir -p ${outdir}							

echo ${outdir}
echo ${ClonalityDir}


