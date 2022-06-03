echo "Program Started: $(date)" > timelog.txt
#in order to run this program, you need to have all fastq's in a directory called fastq and the appropriate genome files (.gtf and HiSat2 build) in a folder called genome. Run this from directory containing both folders. Additionally, it would be possible to pigs-d into a directory on the NVMe drive, read time shouldn't be limiting, process the files, and then mv all of the files off of the drive upon completion.
source /etc/profile.d/apps.sh

parallel gunzip ::: ./fastq/*.gz

wait

echo "Files Unzipped: $(date)" >> timelog.txt
mkdir stats

fastqc ./fastq/*.fastq -o ./stats -t 16

for FILE in ./fastq/*.fastq
    do
        NUM=$(wc -l $FILE | cut -f1 -d" ")
        NUMR=$(echo "$NUM"/4 | bc)
        echo -e "${FILE}\t${NUMR}" >> ./stats/readCounts.txt
    done
wait

for FILE in ./fastq/*.fastq
    do  
        {
        cutadapt \
        -j 0 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        --minimum-length=25 \
        -o "${FILE/.fastq/_trimmed.fastq}" \
        $FILE > "${FILE/.fastq/_trimmingStats.txt}" 
        } &
    done
wait
echo "Adapters Trimmed: $(date)" >> timelog.txt
mkdir trimmedFastq
for FILE in ./fastq/*_trimmed.fastq
	do
	   {
		mv $FILE ./trimmedFastq
	   } &
	done
wait

for FILE in ./fastq/*_trimmingStats.txt
	do
	{
		mv $FILE ./stats
	} &
	done
wait

cat ./stats/*_trimmingStats.txt > ./stats/CombinedTrimmingStats.txt

for FILE in ./trimmedFastq/*.fastq
    do
        NUM=$(wc -l $FILE | cut -f1 -d" ")
        NUMR=$(echo "$NUM"/4 | bc)
        echo -e "${FILE}\t${NUMR}" >> ./stats/readCountsTrimmed.txt
    done
wait

mkdir BAM

for FILE in ./trimmedFastq/*.fastq
    do
	{
    	FN=$(basename $FILE) 
    	hisat2 \
    	-p 16 \
    	--phred33 \
    	--rna-strandness R \
    	--dta \
    	--no-unal \
    	-x /data/Austin/workdir/genome/HiSat2_ENSEMBL_galGal6/ENSEMBL_galGal6 \
        --summary-file ./stats/${FN/hisatSummary} \
    	-U $FILE |
    	samtools view -bS -@ 16 - |
    	samtools sort -@ 16 - > "./BAM/${FN/fastq/BAM}" &&
    	samtools index -@ 16 "./BAM/${FN/fastq/BAM}"
	}
    done
wait
echo "Reads Aligned: $(date)" >> timelog.txt
mkdir counts

	featureCounts \
	./BAM/*.BAM \
	-T 16 \
	-t exon \
	-s 2 \
	-g gene_id \
	-a /data/Austin/workdir/genome/Gallus_gallus.GRCg6a.99.gtf \
	-o ./counts/featureCounts.txt \

echo "Reads Counted at $(date)" >> timelog.txt

multiqc ./stats --cl_config "extra_fn_clean_exts: { '_trimmed.BAM'}"

echo "MultiQC Report Generated at $(date)...Complete!" >> timelog.txt

