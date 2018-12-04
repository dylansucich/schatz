#!/bin/bash
# Run script in the same folder as the paired end reads
# bwa index ecoli.fa


readnames=()
if [ ! -d "aln_out" ]; then
		mkdir "aln_out"
fi

for fname in t*; do
	# if [[ $fname =~ ^-?[0-9]+$ ]]; then
	# 	echo "Already in proper format"

	readnames+=("${fname:0:2}")
	echo "${fname:0:2}"
	# cp "$fname" "aln_out/${fname:0:2}" #${fname:4:}"
	
	# base=${fname%}
	# echo $base 
    # base=${fname%_R1*}  
done

 # bwa-0.7.5a/bwa mem -t 4 human_g1k_v37.fasta "${base}_R1.fastq.gz"  "${base}_R2.fastq.gz" >"$destdir/${base}_R1_R2.sam"

for fname in *.fa; do
	reference=$fname
done
# for i in "$@"; do
#     args+=("$i")
# done

readnames=("${readnames[@]}")
echo
echo "Use these Read files to order more easily in bwa"
echo "${readnames[@]}"
echo
echo "Reference for index: " "$reference"
echo


bwa index ecoli.fa

pwd

maxproc=3

#it is running each timepoint2x...
for fname in *.fq
do
	while [ $(jobs | wc -l) -ge "$maxproc" ]
    do
        sleep 1
    done
	base=${fname:0:2}
	echo "${base}"

	bwa aln -t 4 -I "$reference" "${base}.1.fq" > "aln_out/${base}.1.sai"
	echo "${base}.1.fq" " Aligned with BWA aln"
	bwa aln -t 4 -I "$reference" "${base}.2.fq" > "aln_out/${base}.2.sai"
	echo "${base}.2.fq" " Aligned with BWA aln"
	bwa sampe "$reference" "aln_out/${base}.1.sai" "aln_out/${base}.1.sai" "${base}.1.fq" "${base}.2.fq" > "aln_out/${base}.sam"
	echo "${base}.1.fq" "${base}.2.fq" " Aligned with BWA sampe to SAM file"
	samtools view -bSU SAM "aln_out/${base}.sam" | samtools sort -o "aln_out/${base}_sorted.bam"
	echo "${base}" "SAM file sorted and converted to BAM file with samtools sort piped to view"
done


bwa aln -I "$reference" "t10reads/t10.1.fq" > "t10.1.sai"
echo "t10.1.fq" " Aligned with BWA aln"
bwa aln -I "$reference" "t10reads/t10.2.fq" > "t10.2.sai"
echo "t10.2.fq" " Aligned with BWA aln"
bwa sampe "$reference" "t10.1.sai" "t10.2.sai" "t10reads/t10.1.fq" "t10reads/t10.2.fq" > "t10.sam"
echo "t10.1.fq" "t10.2.fq" " Aligned with BWA sampe to SAM file"
samtools view -bSU SAM "t10.sam" | samtools sort -o "t10_sorted.bam"
echo "t10" "SAM file sorted and converted to BAM file with samtools sort piped to view"



