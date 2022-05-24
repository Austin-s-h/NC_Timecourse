#in order to run this program, you need to have all fastq's in a directory called fastq 
# and the appropriate genome files (.gtf and HiSat2 build) in a folder called genome. 
# Run this from directory containing both folders. 
# Additionally, it would be possible to pigs-d into a directory on the NVMe drive, 
# read time shouldn't be limiting, process the files, and then mv all of the files off of the drive upon completion.

# Now with automatic file parsing!
echo "Program Started: $(date)" > timelog.txt

cd fastq
var=( $(ls) )
if [ "${var[0]##*.}" = "gz" ]; then
       find ./*.fastq.gz -type f -maxdepth 1 | cut -c3- | rev | cut -c4- | rev | paste -sd ';\n' > filePairs.txt
   fi
   
if [ "${var[0]##*.}" = "fastq" ]; then
       find ./*.fastq -type f -maxdepth 1 | cut -c3- | paste -sd ';\n' > filePairs.txt
   fi
cd ..

source /etc/profile.d/apps.sh

parallel gunzip ::: ./fastq/*.gz

wait

echo "Files Unzipped: $(date)" >> timelog.txt

mkdir stats

fastqc ./fastq/*.fastq -o ./stats -t 62

for FILE in ./fastq/*.fastq
    do
        NUM=$(wc -l $FILE | cut -f1 -d" ")
        NUMR=$(echo "$NUM"/4 | bc)
        echo -e "${FILE}\t${NUMR}" >> ./stats/readCounts.txt
    done
wait

cd fastq

for FILEPAIR in $(cat filePairs.txt)
do
    F1=$(echo $FILEPAIR | awk 'BEGIN{FS=";"}{print $1}')
    F2=$(echo $FILEPAIR | awk 'BEGIN{FS=";"}{print $2}')
    echo "Paired End mode ${FILEPAIR}"
  	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  	--minimum-length=25 -j 0 \
  	-o "trimmed_${F1}"  -p "trimmed_${F2}" \
  	${F1} ${F2} >"${F1}_cutadapt_report.txt"
done

wait
echo "Adapters Trimmed: $(date)" >> timelog.txt

mkdir trimmedFastq

for FILE in ./fastq/trimmed*
	do
	   {
		mv $FILE ./trimmedFastq
	   } &
	done
wait

for FILE in ./fastq/*_cutadapt_report.txt
	do
	{
		mv $FILE ./stats
	} &
	done
wait

cat ./stats/*_cutadapt_report.txt > ./stats/CombinedTrimmingStats.txt

for FILE in ./trimmedFastq/*.fastq
    do
        NUM=$(wc -l $FILE | cut -f1 -d" ")
        NUMR=$(echo "$NUM"/4 | bc)
        echo -e "${FILE}\t${NUMR}" >> ./stats/readCountsTrimmed.txt
    done
wait

mkdir BAM

for FILEPAIR in $(cat filePairs.txt)
do
    F1=$(echo $FILEPAIR | awk 'BEGIN{FS=";"}{print $1}')
    F2=$(echo $FILEPAIR | awk 'BEGIN{FS=";"}{print $2}')
    echo "Aligning in Paired End mode ${FILEPAIR}"
    	FN=$(basename $F1) 
    	hisat2 \
    	-p 62 \
    	--phred33 \
    	--rna-strandness RF \
    	--dta \
    	--no-unal \
    	-x /data/Austin/workdir/genome/hg38/hg38_hisat2/genome \
        --summary-file ./stats/${FN/hisatSummary} \
    	-1 $F1 -2 $F2 |
    	samtools view -bS -@ 62 - |
    	samtools sort -@ 62 - > "./BAM/${FN/fastq/bam}" &&
    	samtools index -@ 62 "./BAM/${FN/fastq/bam}"
done

wait

echo "Reads Aligned: $(date)" >> timelog.txt
mkdir counts

	featureCounts \
	./BAM/*.bam \
	-T 62 \
	-t exon \
	-s 2 \
	-g gene_id \
	-a /data/Austin/workdir/genome/hg38/annotation/hg38.ensGene.gtf \
	-o ./counts/featureCounts.txt \

echo "Reads Counted at $(date)" >> timelog.txt

multiqc ./stats --cl_config "extra_fn_clean_exts: { '_trimmed.bam'}"

echo "MultiQC Report Generated at $(date)...Complete!" >> timelog.txt

