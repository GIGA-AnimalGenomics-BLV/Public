#!/bin/bash
#SBATCH --partition=<PARTITION_NAME>
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1  # Line of the samplesheet to process
#SBATCH --job-name=<JOBNAME>
#SBATCH --mail-user=<USERMAIL>
#SBATCH --mail-type=FAIL,END
#SBATCH --output=<PATH/TO/LOGS/DIRECTORY/OUTPUTS>/Step1_Clonality_%A_%a.log

######################
## Coloring message ##
######################
terminalColorClear='\033[0m'
terminalColorEmphasis='\033[1;32m'
terminalColorError='\033[1;31m'
terminalColorMessage='\033[1;33m'
terminalColorWarning='\033[1;34m'
 
echoDefault() {
    echo -e "${terminalColorClear}$1${terminalColorClear}"
}
 
echoMessage() {
    echo -e "${terminalColorMessage}$1${terminalColorClear}"
}
 
echoWarning() {
    echo -e "${terminalColorWarning}$1${terminalColorClear}"
}
 
echoError() {
    echo -e "${terminalColorError}$1${terminalColorClear}"
}

################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
   echoMessage
   echoWarning "This script ...."
   echoMessage
   echoWarning "Syntax: sh clonality_v2.0_step1.sh [-v|h]"
   echoWarning "options:"
   echoMessage "-h     Print this Help."
   echoMessage
}

### Variables to be added as parameters on the command line
echo $SLURM_ARRAY_TASK_ID

while getopts h option

do
 case "${option}"
 in
 h) Help    
 exit;;
 \?)  echoError "Error: Invalid option"
	echoError "For help, type sh clonality_v2.0_step1.sh -h."
     exit;;
 esac
done

if [ -v "$VIRUS" ]; then
	echo "VIRUS is not defined" 
	exit 0; 
	else echo "VIRUS is set to '$VIRUS'"; 
fi

# Load Environment
module load slurm
module load EasyBuild
module load java/8.0.111
module load R/4.2.2-foss-2022b

# Define Variable
source $PWD/Config.sh

# 3. Extract parameters from sampleSheet
r1=`sed -n "$SLURM_ARRAY_TASK_ID"p $files |  awk -F'\t' 'BEGIN {FS="\t"}; {print $12}' | sed 's/\"//g'` # R1 file
r2=`sed -n "$SLURM_ARRAY_TASK_ID"p $files |  awk -F'\t' 'BEGIN {FS="\t"}; {print $13}' | sed 's/\"//g'` # R2 file
name=`sed -n "$SLURM_ARRAY_TASK_ID"p $files |  awk -F'\t' 'BEGIN {FS="\t"}; {print $7}' | sed 's/\"//g'` # sample name
linker=`sed -n "$SLURM_ARRAY_TASK_ID"p $files |  awk -F'\t' 'BEGIN {FS="\t"}; {print $14}' | sed 's/\"//g'` #linker sequence
runDate=`sed -n ${SLURM_ARRAY_TASK_ID}p $files |  awk -F'\t' 'BEGIN {FS="\t"}; {print $15}' | sed 's/\"//g'` #runDate
species=`sed -n "$SLURM_ARRAY_TASK_ID"p $files |  awk -F'\t' 'BEGIN {FS="\t"}; {print $1}' | sed 's/\"//g'` #species (OAR,BTA or HSA)

# 4. Create final output directory name
name="$name"=`echo "$runDate" | sed 's/\;/\_/g' | sed 's/\"//g'`=`echo "$linker" | sed 's/\;/\_/g' | sed 's/\"//g'`
echo $name

# 5. Create final output directory
GROUP=`sed -n "$SLURM_ARRAY_TASK_ID"p $files |  awk -F'\t' 'BEGIN {FS="\t"}; {print $5}'`
mkdir -p ${outdir}/${GROUP}/Mapping
dir=${outdir}/${GROUP}/Mapping/${name}
echo ${dir}

# 5. Mapping parameters
minQ=30
thread=$SLURM_JOB_CPUS_PER_NODE

### |---------------------------------- CLONALITY ----------------------------------| ###

printf "date: `date +"%m-%d-%y\ntime start: %T"`\n\n\n"

printf "############ START CLONALITY ANALYSIS ############\n"
printf "sample: $name\n"
printf "R1: $r1\n"
printf "R2: $r2\n"
printf "Dir: $dir\n"
printf "linker: $linker\n"
printf "Species: $species\n"
printf "Threads: $thread\n"
printf "#####################################\n\n\n"

printf "\n# 1. ----------- CREATE DIRECTORIES #\n\n"


if [ -d $dir ]; then
  echo "Directory exists."
  rm -r $dir
fi
mkdir -p $dir
cd $dir

mkdir -p raw
mkdir -p linkerBowtie
mkdir -p 1_read_linkesrFound
mkdir -p 2_read_provirusFound
mkdir -p 3_read_hostGenome
mkdir -p 4_read_hostAlign
mkdir -p 5_randomTag
mkdir -p 6_checkUp
mkdir -p 7_IGVmapping

printf "\n# 2. ----------- TRANSFERT FILES LOCALY #\n\n"

# NB: If Nextseq500 put in R1 the fastq1;fastq2;fastq3;etc

perl $concatenateFASTQ $r1 raw/"$name"_R1.fastq.gz
perl $concatenateFASTQ $r2 raw/"$name"_R2.fastq.gz

gzip -dc raw/"$name"_R1.fastq.gz > raw/R1.fastq
gzip -dc raw/"$name"_R2.fastq.gz > raw/R2.fastq

R1=raw/R1.fastq
R2=raw/R2.fastq

printf "\n# 3. ----------- SELECT READS WITH LINKER SEQUENCE #\n\n"

printf "\n# 3.1. ----------- Index linker sequence(s) #\n"

# When several runs are done on the same sample we often end-up with different linkers
# Given a sample sheet with linkers separated by ; (linker1;linker2;linker3)
# Create an index containing all possible linker sequences

cd linkerBowtie

linkerSplit=$(echo $linker | tr ";" "\n")

touch linker.fa
for j in $linkerSplit
	do
		echo ">linker_$j" >> linker.fa
		echo "$j" >> linker.fa
		j=$((j+1))
	done

singularity exec ${SingularityDir}/bowtie_v1_3_1.sif bowtie-build -q linker.fa linker

cd ..

printf "\n# 3.2. ----------- Extract reads with linker sequence(s) #\n"

# First extract the 9nt to 9nt+linkerLength nucleotides from R2
# Align those sequences to the bowtie linker index (-l Seed length of 7, -N 0 Report alignment without)
# Find the ID of the reads with the right linkers (egrep | read_id)
# Filter the R1/R2 fastq files to get the one with the ID from read_id_linker.txt

cd 1_read_linkesrFound

singularity exec ${SingularityDir}/fastx_toolkit_v0_0_13.sif fastx_trimmer -Q33 -f 9 -l $linkLen -i ../raw/R2.fastq -o read_linkers.fastq
singularity exec ${SingularityDir}/bowtie_v1_3_1.sif bowtie -l 7 -p $thread -v 0 ../linkerBowtie/linker read_linkers.fastq | egrep `echo $linker | sed 's/\;/\|/g'` | awk '{print $1}' > readID_linkerFound.txt
singularity exec ${SingularityDir}/bbmap_v39_13.sif filterbyname.sh in=../raw/R1.fastq in2=../raw/R2.fastq out=R1_linkerFound.fastq out2=R2_linkerFound.fastq names=readID_linkerFound.txt include=t overwrite=t

cd ..

printf "\n# 4. ----------- SELECT READS WITH ALIGNING ONTO THE PROVIRUS #\n\n"

printf "\n# 4.1. ----------- Extract sequence supposed to align (R1 reads) #\n"

cd 2_read_provirusFound

## Extract part of read R1 that must align on the virus
# No mismatches allowed
singularity exec ${SingularityDir}/fastx_toolkit_v0_0_13.sif \
		fastx_trimmer -Q33 -f 8 -l $LTR3len \
		-i ../1_read_linkesrFound/R1_linkerFound.fastq \
		-o R1_virusLTR3_edge.fastq
		
singularity exec ${SingularityDir}/fastx_toolkit_v0_0_13.sif \
		fastx_trimmer -Q33 -f 8 -l $LTR5len \
		-i ../1_read_linkesrFound/R1_linkerFound.fastq \
		-o R1_virusLTR5_edge.fastq

printf "\n# 4.2. ----------- Align onto the LTR sequence #\n"

# -L length of seed substrings; must be >3, <32 (22)
# -N max # mismatches in seed alignment; can be 0 or 1 (0)
# -x bowtie2-index
# -U unpaired_data.fastq
# Find the read ID and separate LTR3 from LTR5 [grep "LTR3" | awk '{print $1}' | grep -v "^@"]
singularity exec ${SingularityDir}/bowtie2_v2_5_4.sif \
		bowtie2 -L 14 -p $thread \
		-N 1 -x $bowtie2LTR \
		-U R1_virusLTR3_edge.fastq | grep "LTR3" | awk '{print $1}' | grep -v "^@" \
		> readID_LTR3Found.txt
		
singularity exec ${SingularityDir}/bowtie2_v2_5_4.sif  \
		bowtie2 -L 14 -p $thread \
		-N 1 -x $bowtie2LTR \
		-U R1_virusLTR5_edge.fastq | grep "LTR5" | awk '{print $1}' | grep -v "^@" \
		> readID_LTR5Found.txt

printf "\n# 4.3. ----------- Extract reads with LTR & Linker sequence #\n"

# Interstect R1/R2.linker.fastq with R1.LTR3.fastq
singularity exec ${SingularityDir}/bbmap_v39_13.sif \
		filterbyname.sh \
		in=../1_read_linkesrFound/R1_linkerFound.fastq \
		in2=../1_read_linkesrFound/R2_linkerFound.fastq \
		out=R1_LTR3Found_linkerFound.fastq \
		out2=R2_LTR3Found_linkerFound.fastq \
		names=readID_LTR3Found.txt \
		include=t \
		overwrite=t
		
singularity exec ${SingularityDir}/bbmap_v39_13.sif \
		filterbyname.sh \
		in=../1_read_linkesrFound/R1_linkerFound.fastq \
		in2=../1_read_linkesrFound/R2_linkerFound.fastq \
		out=R1_LTR5Found_linkerFound.fastq \
		out2=R2_LTR5Found_linkerFound.fastq \
		names=readID_LTR5Found.txt \
		include=t overwrite=t

# To properly remove the LTR sequence, you need to keep the everything from the end of LTRlen (+1) to the end.
# To obtain LTRlen +1 use ((x++))
((LTR5len++))
((LTR3len++))
((linkLen++))

cd ..

printf "\n# 5. ----------- EXTRACT HOST SEQUENCES FROM SELECTED READS #\n\n"

printf "\n# 5.1. ----------- Trim reads to remove LTR sequences (R1 reads) #\n"

cd 3_read_hostGenome

singularity exec ${SingularityDir}/fastx_toolkit_v0_0_13.sif  \
		fastx_trimmer -Q33 -f $linkLen \
		-i ../2_read_provirusFound/R2_LTR3Found_linkerFound.fastq \
		-o R2_LTR3Found_linkerFound_Host.fastq
singularity exec ${SingularityDir}/fastx_toolkit_v0_0_13.sif  \
		fastx_trimmer -Q33 -f $linkLen \
		-i ../2_read_provirusFound/R2_LTR5Found_linkerFound.fastq \
		-o R2_LTR5Found_linkerFound_Host.fastq
#
singularity exec ${SingularityDir}/fastx_toolkit_v0_0_13.sif \
		fastx_trimmer -Q33 -f $LTR3len \
		-i ../2_read_provirusFound/R1_LTR3Found_linkerFound.fastq \
		-o R1_LTR3Found_linkerFound_Host.fastq
singularity exec ${SingularityDir}/fastx_toolkit_v0_0_13.sif \
		fastx_trimmer -Q33 -f $LTR5len \
		-i ../2_read_provirusFound/R1_LTR5Found_linkerFound.fastq \
		-o R1_LTR5Found_linkerFound_Host.fastq

printf "\n# 5.2. ----------- Take into account overlapping fragments #\n"

## In case of overlapping reads (short insert) the end of one read is often the beginning of the LTR
# In order to remove them, cutadapt take the begining of this sequence and remove everything that is after
# - m minimum length (discard reads too short after trimming)
# -g front adapter ligated at 5'end - this is removing adapters and extra sequence from the front
# -n 2 looking for at the max 2 of the given sequences/adapters to remove. But if its just one, that will be removed as well. In the command befor already expected adapter/liker/LTR is removed from 5 prime end of th read. So this is another LTR or LTR missed from the previous command. Also to remove LTR from R2?

# R1
if [ "$species" == "Oar" ] || [ "$species" == "Bta" ]
then
singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -n 2 -m 10 -g CCGGCAAACA R1_LTR3Found_linkerFound_Host.fastq | \
	singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -n 2 -m 10 -g GGCCGTTTGT - > R1_LTR3Found_linkerFound_Host.trimmed.fastq
singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -n 2 -m 10 -g GAAAGTATGT R1_LTR5Found_linkerFound_Host.fastq | \
	singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -n 2 -m 10 -g CTTTCATACA - > R1_LTR5Found_linkerFound_Host.trimmed.fastq
elif [ "$species" == "Hsa" ]
then
# TGACAATGACCATGAGCCCCAAATATCCC			LTR5
# GACAGCCCATCCTATAGCACTCTCAGGAGAGAAATTTAGTACACA	LTR3
singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -n 2 -m 10 -g TTAGTACACA R1_LTR3Found_linkerFound_Host.fastq | \
	singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -n 2 -m 10 -g AATCATGTGT - > R1_LTR3Found_linkerFound_Host.trimmed.fastq
singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -n 2 -m 10 -g TGACAATGAC R1_LTR5Found_linkerFound_Host.fastq | \
	singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -n 2 -m 10 -g ACTGTTACTG - > R1_LTR5Found_linkerFound_Host.trimmed.fastq
elif [ "$species" == "HSA_HIV" ]
then
singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -n 2 -m 10 -g GCCTGTACTG R1_LTR3Found_linkerFound_Host.fastq | \
	singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -n 2 -m 10 -g CGGACATGAC - > R1_LTR3Found_linkerFound_Host.trimmed.fastq
singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -n 2 -m 10 -g TGGAAGGGCT R1_LTR5Found_linkerFound_Host.fastq | \
	singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -n 2 -m 10 -g ACCTTCCCGA - > R1_LTR5Found_linkerFound_Host.trimmed.fastq
fi

# R2
if [ "$species" == "Oar" ] || [ "$species" == "Bta" ]
then
singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -m 10 -a CAAACGGCCAGAGAG R2_LTR3Found_linkerFound_Host.fastq | \
	singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -m 10 -a GTTTGCCGGTCTCTC - > R2_LTR3Found_linkerFound_Host.trimmed.fastq
singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -m 10 -a ACATACTTTCTAGTA R2_LTR5Found_linkerFound_Host.fastq | \
	singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -m 10 -a TGTATGAAAGATCAT - > R2_LTR5Found_linkerFound_Host.trimmed.fastq
elif [ "$species" = "Hsa" ]
then
singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -m 10 -a ACACATGATT R2_LTR3Found_linkerFound_Host.fastq | \
	singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -m 10 -a TGTGTACTAA - > R2_LTR3Found_linkerFound_Host.trimmed.fastq
singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -m 10 -a ACTGTTACTG R2_LTR5Found_linkerFound_Host.fastq | \
	singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -m 10 -a TGACAATGAC - > R2_LTR5Found_linkerFound_Host.trimmed.fastq
elif [ "$species" == "HSA_HIV" ]
then
singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -n 2 -m 10 -g GTCATGTCCG R2_LTR3Found_linkerFound_Host.fastq | \
	singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -n 2 -m 10 -g CAGTACAGGC - > R2_LTR3Found_linkerFound_Host.trimmed.fastq
singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -n 2 -m 10 -g ACCTTCCCGA R2_LTR5Found_linkerFound_Host.fastq | \
	singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -n 2 -m 10 -g TGGAAGGGCT - > R2_LTR5Found_linkerFound_Host.trimmed.fastq
fi 

printf "\n# 5.3. ----------- Trim reads to remove Linker sequences (R1 reads) #\n"

# Separate the linkers in an array (based on ";")
# Should be unique

IFS=";" read -ra linkers <<< "$linker"

for link in `echo ${linkers[@]} | tr [:space:] '\n' | awk '!a[$0]++'`
	do
 		echo 'Remove linker: '$link
 		revLinker=$(echo $link | rev)
 		revCompLinker=$(echo $revLinker | tr "[ATGCatgcNn]" "[TACGtacgNn]")
 		singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -m 10 -a $revLinker R1_LTR3Found_linkerFound_Host.trimmed.fastq | \
		singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt -m 10 --quiet -a $revCompLinker - > R1_LTR3Found_linkerFound_Host.trimmed2.fastq
 		singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt --quiet -m 10 -a $revLinker R1_LTR5Found_linkerFound_Host.trimmed.fastq | \
		singularity exec ${SingularityDir}/cutadapt_v5_0.sif cutadapt -m 10 --quiet -a $revCompLinker - > R1_LTR5Found_linkerFound_Host.trimmed2.fastq
	done

printf "\n# 5.4. ----------- Cleanup: Resynchronize pairs of reads #\n"

# Intersect R1/R2 in case of too short R1 after trimming

python $resynchronizePairs R1_LTR3Found_linkerFound_Host.trimmed2.fastq R2_LTR3Found_linkerFound_Host.trimmed.fastq
python $resynchronizePairs R1_LTR5Found_linkerFound_Host.trimmed2.fastq R2_LTR5Found_linkerFound_Host.trimmed.fastq

R1_LTR3="3_read_hostGenome/R1_LTR3Found_linkerFound_Host.trimmed2.fastq_pairs_R1.fastq"
R2_LTR3="3_read_hostGenome/R2_LTR3Found_linkerFound_Host.trimmed.fastq_pairs_R2.fastq"
R1_LTR5="3_read_hostGenome/R1_LTR5Found_linkerFound_Host.trimmed2.fastq_pairs_R1.fastq"
R2_LTR5="3_read_hostGenome/R2_LTR5Found_linkerFound_Host.trimmed.fastq_pairs_R2.fastq"

cd ..

printf "\n# 6. ----------- ALIGN ONTO THE HOST GENOME #\n\n"

printf "\n# 6.1. ----------- Bowtie2: Default parameters #\n"

cd 4_read_hostAlign

# -N 1 = One mismatch allowed
# --no-mixed = If only one mate is aligning, do not report the alignment

printf "\n# 6.2. ----------- Bowtie2: Mulimapping allowed #\n"

printf "\n### LTR3_BEST ###\n"
singularity exec ${SingularityDir}/bowtie2_v2_5_4.sif \
		bowtie2 -p $thread \
		--very-sensitive \
		-N 1 \
		--no-mixed \
		-x $bowtie2HostVirus \
		-1 ../$R1_LTR3 \
		-2 ../$R2_LTR3 \
		> LTR3_candidateIS_bowtie2_BEST.sam
printf "\n\n\n\n\n\n"
printf "\n### LTR5_BEST ###\n"
singularity exec ${SingularityDir}/bowtie2_v2_5_4.sif \
		bowtie2 -p $thread \
		--very-sensitive \
		-N 1 \
		--no-mixed \
		-x $bowtie2HostVirus \
		-1 ../$R1_LTR5 \
		-2 ../$R2_LTR5 \
		> LTR5_candidateIS_bowtie2_BEST.sam
printf "\n\n\n\n\n\n"

#module load SAMtools/1.17-GCC-12.2.0
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools view -Sb LTR3_candidateIS_bowtie2_BEST.sam -o LTR3_candidateIS_bowtie2_BEST.bam
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools sort LTR3_candidateIS_bowtie2_BEST.bam -o LTR3_candidateIS_bowtie2_BEST.sorted.bam
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools index LTR3_candidateIS_bowtie2_BEST.sorted.bam
 
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools view -Sb LTR5_candidateIS_bowtie2_BEST.sam -o LTR5_candidateIS_bowtie2_BEST.bam
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools sort LTR5_candidateIS_bowtie2_BEST.bam -o LTR5_candidateIS_bowtie2_BEST.sorted.bam
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools index LTR5_candidateIS_bowtie2_BEST.sorted.bam

# 

printf "\n### LTR3_ALTERN ###\n"
singularity exec ${SingularityDir}/bowtie2_v2_5_4.sif bowtie2 -p $thread --very-sensitive -N 1 --no-mixed -k 11 -x $bowtie2HostVirus -1 ../$R1_LTR3 -2 ../$R2_LTR3 > LTR3_candidateIS_bowtie2_ALTERN.sam
printf "\n\n\n\n\n\n"
printf "\n### LTR5_ALTERN ###\n"
singularity exec ${SingularityDir}/bowtie2_v2_5_4.sif bowtie2 -p $thread --very-sensitive -N 1 --no-mixed -k 11 -x $bowtie2HostVirus -1 ../$R1_LTR5 -2 ../$R2_LTR5 > LTR5_candidateIS_bowtie2_ALTERN.sam
printf "\n\n\n\n\n\n"

singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools view -Sb LTR3_candidateIS_bowtie2_ALTERN.sam -o LTR3_candidateIS_bowtie2_ALTERN.bam
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools sort LTR3_candidateIS_bowtie2_ALTERN.bam -o LTR3_candidateIS_bowtie2_ALTERN.sorted.bam
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools index LTR3_candidateIS_bowtie2_ALTERN.sorted.bam
 
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools view -Sb LTR5_candidateIS_bowtie2_ALTERN.sam -o LTR5_candidateIS_bowtie2_ALTERN.bam
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools sort LTR5_candidateIS_bowtie2_ALTERN.bam -o LTR5_candidateIS_bowtie2_ALTERN.sorted.bam
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools index LTR5_candidateIS_bowtie2_ALTERN.sorted.bam

cd .. 

printf "\n# 7. ----------- EXTRACT RANDOM TAG SEQUENCE (R2 READS) #\n\n"

# -Q33 for quality minimal == 33
# -l last base to keep (1 to 8) 
# -t output in a tabulated format (=fasta)
# Then retrieve from the formated output the first and third column (ID	"\t"	randomTag)

cd 5_randomTag

singularity exec ${SingularityDir}/fastx_toolkit_v0_0_13.sif fastx_trimmer -Q33 -l 8 -i ../raw/R2.fastq | \
singularity exec ${SingularityDir}/fastx_toolkit_v0_0_13.sif fastq_to_fasta -Q33 | \
singularity exec ${SingularityDir}/fastx_toolkit_v0_0_13.sif fasta_formatter -t | awk '{print $1"\t"$3}' > R2_randomTag.txt

cd ..

printf "\n# 8. ----------- CLONALITY RESULTS TABLES AND STATISTICS #\n\n"
Rscript /scratch/GIGA/USER/u215238/Clonality_AVDB/Scripts/R/run_PIC.R "4_read_hostAlign/LTR3_candidateIS_bowtie2_BEST.sorted.bam" "4_read_hostAlign/LTR5_candidateIS_bowtie2_BEST.sorted.bam" "4_read_hostAlign/LTR3_candidateIS_bowtie2_ALTERN.sorted.bam" "4_read_hostAlign/LTR5_candidateIS_bowtie2_ALTERN.sorted.bam" "5_randomTag/R2_randomTag.txt" $name $geneBedFile "raw/R1.fastq" $virus 30 25

printf "\n# 9. ----------- Prepare data for IGV visualisation #\n\n"

cd 7_IGVmapping

printf "\n# 9.1. ----------- Get BAM pre-trimming Filtering (linker and virus Ok) #\n"

singularity exec ${SingularityDir}/bowtie2_v2_5_4.sif bowtie2 -p $thread --very-sensitive -N 1 -x $bowtie2HostVirus -1 ../2_read_provirusFound/R1_LTR3Found_linkerFound.fastq -2 ../2_read_provirusFound/R2_LTR3Found_linkerFound.fastq > LTR3_candidateIS_PreTrimming_bowtie2.sam
singularity exec ${SingularityDir}/bowtie2_v2_5_4.sif bowtie2 -p $thread --very-sensitive -N 1 -x $bowtie2HostVirus -1 ../2_read_provirusFound/R1_LTR5Found_linkerFound.fastq -2 ../2_read_provirusFound/R2_LTR5Found_linkerFound.fastq > LTR5_candidateIS_PreTrimming_bowtie2.sam

# SAM TO BAM and INDEX
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools view -bS LTR3_candidateIS_PreTrimming_bowtie2.sam > LTR3_candidateIS_PreTrimming_bowtie2.bam
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools view -bS LTR5_candidateIS_PreTrimming_bowtie2.sam > LTR5_candidateIS_PreTrimming_bowtie2.bam
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools merge  "$name"_candidateIS_PreTrimming_bowtie2.bam LTR3_candidateIS_PreTrimming_bowtie2.bam LTR5_candidateIS_PreTrimming_bowtie2.bam
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools sort "$name"_candidateIS_PreTrimming_bowtie2.bam -o "$name"_candidateIS_PreTrimming_bowtie2.sorted.bam
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools index "$name"_candidateIS_PreTrimming_bowtie2.sorted.bam

printf "\n# 9.2. ----------- Get BAM post-trimming Filtering #\n"

# SAM TO BAM and INDEX
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools merge "$name"_candidateIS_PostTrimming_bowtie2.bam ../4_read_hostAlign/LTR3_candidateIS_bowtie2_BEST.bam ../4_read_hostAlign/LTR5_candidateIS_bowtie2_BEST.bam
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools sort "$name"_candidateIS_PostTrimming_bowtie2.bam -o "$name"_candidateIS_PostTrimming_bowtie2.sorted.bam
singularity exec ${SingularityDir}/SAMtools_v1_21.sif samtools index "$name"_candidateIS_PostTrimming_bowtie2.sorted.bam

rm LTR*_candidateIS_*_bowtie2.bam
rm LTR*_candidateIS_PreTrimming_bowtie2.sam

cd ..

printf "\n# 10. ----------- CHECK-UP INFORMATIONS #\n\n"

printf "\n# 10.1. ----------- Picard Tools: Insert Size #\n"

cd 6_checkUp

singularity exec ${SingularityDir}/Picard-Tools_v2_18.sif picard-tools CollectInsertSizeMetrics VERBOSITY=ERROR I=../7_IGVmapping/"$name"_candidateIS_PreTrimming_bowtie2.sorted.bam O="$name"_candidateIS_PreTrimming_insertSizeMetrics.txt H="$name"_candidateIS_PreTrimming_insertSizeHistogram.pdf M=0.5
singularity exec ${SingularityDir}/Picard-Tools_v2_18.sif picard-tools CollectInsertSizeMetrics VERBOSITY=ERROR I=../7_IGVmapping/"$name"_candidateIS_PostTrimming_bowtie2.sorted.bam O="$name"_candidateIS_PostTrimming_insertSizeMetrics.txt H="$name"_candidateIS_PostTrimming_insertSizeHistogram.pdf M=0.5

cd ..

printf "\n# 11. ----------- TIDY-UP #\n\n"

mv raw/R1.fastq raw/"$name"_R1.fastq
mv raw/R2.fastq raw/"$name"_R2.fastq

rm raw/"$name"_R1.fastq
rm raw/"$name"_R2.fastq

printf "\n# ! ----------- ! -----------  DONE !!!  ----------- ! ----------- ! #\n\n"

printf "date: `date +"%m-%d-%y\ntime end: %T"`\n\n\n"
