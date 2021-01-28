#!/bin/bash

# Preparation and setup of required files

FILES=$(pwd)
WKDIR=$(echo $FILES | sed 's:/required_files::g')

# Ask for raw data file format

echo 'Are the data in bam or fastq format?'
read -p 'Specify file format (bam or fastq): ' FORMAT

read -p 'Do you want to do a quality control of the raw data (yes or no): ' QCRAW

read -p 'How many threads (cores) should be used for the analysis (use 1, if you are not sure): ' THREAD


GENOME=$WKDIR/required_files/C_auris_WGS.fasta
FEATURES=$WKDIR/required_files/C_auris_annotation.gff3
ADAPT1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
mkdir $WKDIR/QC
PICARD=$WKDIR/required_files/picard.jar


# QC of raw data

if [ $QCRAW == 'yes' ]
then
	mkdir $WKDIR/QC_raw
	echo 'Quality control of raw data:'
	if [ $FORMAT == 'bam' ]
	then
		for i in $WKDIR/*.bam
		do
			fastqc -o $WKDIR/QC_raw $i
		done
	else
		for SNAME in $(ls $WKDIR | egrep '(\.f.*q$)|(q\.gz$)')
		do
			i=$WKDIR/$SNAME
			fastqc -o $WKDIR/QC_raw $i
		done
	fi
	multiqc -o $WKDIR/QC_raw $WKDIR/QC_raw
else
	echo 'No QC of raw data done.'
fi


# Adapter removal with cutadapt and mapping of all files with NGM

for SNAME in $(ls $WKDIR | egrep '(\.f.*q$)|(q\.gz$)')
do
	i=$WKDIR/$SNAME # sets sample name and file path
	
	cutadapt -q 30 -a $ADAPT1 $i > $i.trimmed.fq 2>$WKDIR/QC/Cutadapt_$SNAME.txt
  rm $i
  
  ngm -q $i.trimmed.fq -r $GENOME -o $i.trimmed.fq.bam -b -t $THREAD
  rm $i.trimmed.fq
  
  samtools sort -@ $THREAD $i.trimmed.fq.bam -o $i.trimmed.fq.bam.sort.bam   # sort .bam files using samtools
  rm $i.trimmed.fq.bam
  
  # Removal of optical duplicates
	java -jar $PICARD MarkDuplicates REMOVE_SEQUENCING_DUPLICATES=true I=$i.trimmed.fq.bam.sort.bam O=$i.trimmed.fq.bam.sort.bam.markdup.bam M=$WKDIR/QC/$SNAME.markdup.metrics.txt
	rm $i.trimmed.fq.bam.sort.bam
  
  #Quality control and statistics about mapped samples
	samtools flagstat $i.trimmed.fq.bam.sort.bam.markdup.bam >> $WKDIR/QC/$SNAME.final.flagstat_analysis.txt   # flagstat analysis

	fastqc -o $WKDIR/QC $i.trimmed.fq.bam.sort.bam.markdup.bam

done

multiqc -s -o $WKDIR/QC $WKDIR/QC


# Preparation of coverage files for visualization in IGV

mkdir $WKDIR/IGV_files

for i in $WKDIR/*.markdup.bam
do
	samtools index $i
	SNAME=$(echo $i | sed 's:/.*/::g')
	bamCoverage -b $i -o $WKDIR/IGV_files/$SNAME.bw --normalizeUsing CPM -p $THREAD
done


# Count reads

mkdir $WKDIR/count
mkdir $WKDIR/diff_expr_analysis

for i in *.markdup.bam
  do
  htseq-count -f bam -s no -t gene -i ID $i $FEATURES > $i.count.txt
  mv $i.count.txt $WKDIR/count
done

for i in $WKDIR/count/*.count.txt
do
	head -n -5 $i > $i.crop.txt  # clear count files for flags
done


cp $WKDIR/count/*.crop.txt $WKDIR/diff_expr_analysis
cp $FILES/edgeR_analysis.R $WKDIR/diff_expr_analysis
cp $FILES/Targets.txt $WKDIR/diff_expr_analysis




