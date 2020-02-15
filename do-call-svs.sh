#! /bin/bash

###############################################################################

set -o errexit
set -o nounset
set -o pipefail
set -o xtrace

#export PATH=.:/home/edrabek/bin:${PATH:-}
#export PYTHONPATH=/home/edrabek/lib/python:${PYTHONPATH:-}
#export LC_ALL=C

export PATH=/usr/local/packages/python-2.7.14/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/packages/python-2.7.14/lib

###############################################################################
    #ABOUT SCRIPT

#July 2018, Kara Moser

#ABOUT: Generate a catolouge of high-confidence SVs supported by a de novo 
#       assembly of PacBio reads and read-mapping information of the reads
#       to a reference genome

#REQUIREMENTS: RHEL7 machine

#INPUT1: 1) PacBio reads (*fastq)
#	 2) Polished assembly (quiver 2x, pilon)

#OUTPUT: 1) Bam files of PacBio reads aligned to reference + assembly for 
#           visualization purposes
#	 2) VCF file of structural variants detected

###############################################################################

#Paths to necessary software

ngmlr=/usr/local/packages/ngmlr-0.2.6/ngmlr

#Scripts have been modified from the github source code to report inversions
assemblytics=/local/projects-t3/p_falciparum/kmoser/thesis/assemblies/Assemblytics-master

sniffles=/home/kara.moser/bin/Sniffles-1.0.6/bin/sniffles-core-1.0.6/sniffles

survivor=/usr/local/packages/survivor/bin/SURVIVOR

vcftools=/usr/local/packages/vcftools-0.1.15/bin

###############################################################################

#options

sample=$1
country=$2

work_d=/local/scratch/kmoser/${sample}_structural_variants

#data directories
data_d=/local/projects-t3/p_falciparum/samples/$country/$sample

ref=/local/projects-t3/p_falciparum/auxiliary_files/reference.fa

###############################################################################

#Preparation of bam files for structural variant identification

###############################################################################

mkdir -p $work_d/read_alignments
cd $work_d/read_alignments

#Prepare PacBio reads

g=combined.fastq
if [[ ! -e $g ]]; then
  cat $data_d/PACBIO_DATA/*fastq > $g
fi

#Align pacbio reads to the assembly (nice for visualizing things with IGV/Ribbon)

f=pacbio_reads_vs_assembly.sorted.bam
  if [[ ! -e $f ]]; then

  $ngmlr -t 4 -r $data_d/canu_assembly/final_assembly.fasta \
         -q $g -o pacbio_reads_vs_assembly.sam

  samtools view -S pacbio_reads_vs_assembly.sam -q 0 -b -o pacbio_reads_vs_assembly.bam
  rm pacbio_reads_vs_assembly.sam

  samtools sort pacbio_reads_vs_assembly.bam -o pacbio_reads_vs_assembly.sorted.bam
  samtools index pacbio_reads_vs_assembly.sorted.bam

fi

#Align pacbio reads to 3D7 (needed for Sniffles/structural variant detection)

h=pacbio_reads_vs_3D7.sorted.bam
if [[ ! -e $h ]]; then

  $ngmlr -t 4 -r $ref \
         -q $g -o pacbio_reads_vs_3D7.sam


  samtools view -S pacbio_reads_vs_3D7.sam -q 0 -b -o pacbio_reads_vs_3D7.bam
  rm pacbio_reads_vs_3D7.sam

  samtools sort pacbio_reads_vs_3D7.bam -o pacbio_reads_vs_3D7.sorted.bam
  samtools index pacbio_reads_vs_3D7.sorted.bam

  rm pacbio_reads_vs_3D7.bam

fi

###############################################################################

#Identifying structural variants with Sniffles

###############################################################################

ra=$work_d/read_alignments

mkdir -p $work_d/svs_from_sniffles
cd $work_d/svs_from_sniffles

#Identifying structural variants with sniffles
i=structural_variants.vcf
if [[ ! -e $i ]]; then
 $sniffles -m $ra/$h -v $i
fi

#Filtering structural variants by removoing regions with low MQ and low read support
j=structural_variants.filtered.vcf
if [[ ! -e $j ]]; then

 samtools view -H $ra/pacbio_reads_vs_3D7.sorted.bam > $ra/lowMQ.sam
 samtools view $ra/pacbio_reads_vs_3D7.sorted.bam | awk '$5<5 {print $0}' >>  $ra/lowMQ.sam
 samtools view -S -b -h $ra/lowMQ.sam > $ra/lowMQ.bam
 samtools depth $ra/lowMQ.bam >  $ra/lowMQ.cov

 $survivor bincov $ra/lowMQ.cov 10 2 > $ra/lowMQ.bed

 #Filter (ignore size restrictions and a MAF, + SV's not supported by at least 5 reads)
 $survivor filter $i $ra/lowMQ.bed -1 -1 -1 5 $j

 #Make sure vcf file is filtered (potential bug means this is critical for downstream
 # merge step
 cat $j | $vcftools/vcf-sort > structural_variants.filtered.sorted.vcf

fi

###############################################################################

#Identifying structural variants with Assemblytics

###############################################################################

mkdir -p $work_d/svs_from_assemblytics
cd $work_d/svs_from_assemblytics

k=${sample}.Assemblytics_structural_variants.bed
if [[ ! -e $k ]]; then

  #preparing delta file for assemblytics
  nucmer -maxmatch -l 100 -c 500 $ref $data_d/canu_assembly/final_assembly.fasta -prefix ${sample}_vs_3D7
  gzip ${sample}_vs_3D7.delta

  #Identifying SVs with assemblytics:
    #Name of file
    #Uniq anchor length
    #Path to R scripts
    #Min and max SV size reported
  $assemblytics/Assemblytics ${sample}_vs_3D7.delta.gz ${sample} 1000 $assemblytics 30 10000000

  $survivor convertAssemblytics $k 30 ${sample}.assemblytics.vcf

  #convertAssemblytics will output inversions and translocations, but won't label them correctly
  #this checks if a SV labeled NA is an interchromosomal event (then label as TRA)
  #if not, it's an inversion (then label as INV)
  /home/kara.moser/daily/1807/structural_variant_catologue/fix_assemblytics_file.pl \
    ${sample}.assemblytics.vcf \
    ${sample}.assemblytics.fixed.vcf

  #Make sure vcf file is filtered
  cat ${sample}.assemblytics.fixed.vcf | $vcftools/vcf-sort > structural_variants.filtered.sorted.vcf

fi

###############################################################################

#Creating pooled data set containing variants in both sniffles and assemblytics

###############################################################################

mkdir -p $work_d/merged
cd $work_d/merged

#Creating input file
echo "../svs_from_sniffles/structural_variants.filtered.sorted.vcf" > sample_files
echo "../svs_from_assemblytics/structural_variants.filtered.sorted.vcf" >> sample_files

#Merge structural variants NOTE: MAY WANT TO TWEAK THESE DEPENDING ON YOUR DATASET
  #Breakpoints of the SVs from each dataset need to be within 1300 bp of eachother,
  #Needs to be in both files (2)
  #Does not need to be the same type of SV (0) 
   #assemb. may categorize something as a DUP when sniffles calls it an INS
   #was causing true NF54 errors to be dropped
  #Doesn't  matter if they're on different strands (0)
  #Not estimating distance based on the size of SV (0)
  #Minimum size of SV is the default (30)

$survivor merge sample_files 1300 2 0 0 0 30 merged_1300_2000_30.vcf

#Generating stats from merged file
  #No minimum or maximum SV size specfied (-1, -1) 
  #No minimum number of reads needed to support an SV (-1)

$survivor stats merged_1300_2000_30.vcf -1 -1 -1 merged_1300_2000_30.stats

#get stats
for i in INS DEL DUP INV INVDUP TRA; do
  < merged_1300_2000_30.vcf grep -v "#" \
  | grep -v "PFC10_API_IRAB\|M7661" \
  | grep "<${i}>" \
  | cut -f8 \
  | cut -d';' -f3 \
  | perl -pe 's/AVGLEN=//g' \
  > dist_${i}.txt
done

###############################################################################

echo  "Read through script."

###############################################################################

