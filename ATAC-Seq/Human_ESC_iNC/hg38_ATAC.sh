# Now with automatic file parsing!

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

echo "Program Started: $(date)" > timelog.txt

parallel gunzip ::: ./fastq/*.gz

echo "Files Unzipped: $(date)" >> timelog.txt

cd fastq/

for FILEPAIR in $(cat filePairs.txt)
do
    F1=$(echo $FILEPAIR | awk 'BEGIN{FS=";"}{print $1}')
    F2=$(echo $FILEPAIR | awk 'BEGIN{FS=";"}{print $2}')
    echo "Paired End mode ${FILEPAIR}"
  	cutadapt -a CTGTCTCTTATACACATCT \
  	-A AGATGTGTATAAGAGACAG \
  	--minimum-length=25 -j 0 \
  	-o "trimmed_${F1}"  -p "trimmed_${F2}" \
  	${F1} ${F2} >"${F1}_cutadapt_report.txt"
done

wait

cd ../

echo "Adapters Trimmed: $(date)" >> timelog.txt

mkdir trimmedFastq
for FILE in ./fastq/trimmed*
	do
	   {
		mv $FILE ./trimmedFastq
	   } &
	done
wait

mkdir BAM

mv ./fastq/filePairs.txt ./trimmedFastq/filePairs.txt

cd trimmedFastq/

for FILEPAIR in $(cat filePairs.txt)
do
    F1=$(echo $FILEPAIR | awk 'BEGIN{FS=";"}{print $1}')
    F2=$(echo $FILEPAIR | awk 'BEGIN{FS=";"}{print $2}')
    echo "beginning alignment to GRCh38 of file ${FILEPAIR}"
	bowtie2 --local --very-sensitive-local \
		--no-unal --no-mixed --no-discordant \
		--threads 62 \
		-x /data/Austin/workdir/genome/hg38/hg38_bt2/GRCh38 \
		-I 10 -X 1000 \
		-1 "trimmed_${F1}" \
		-2 "trimmed_${F2}" \
		-S ${F1/_R1.fastq/_toGRCh38.SAM}
	echo "complete"

done
wait

cd ../

echo "Reads Aligned: $(date)" >> timelog.txt


for FILE in ./trimmedFastq/*.SAM; do
    sed -i '/^@PG/d' $FILE
    samtools view -F 1804 -b -@ 62 $FILE |
    samtools sort -@ 62 > "${FILE/.SAM/.bam}" &&
    samtools index -@ 62 "${FILE/.SAM/.bam}"

done

wait

for FILE in ./trimmedFastq/*.bam
	do
	{
		mv $FILE ./BAM
	} &
	done
wait

for FILE in ./trimmedFastq/*.bai
	do
	{
		mv $FILE ./BAM
	} &
	done
wait

for file in ./BAM/*.bam

do
	java -Xmx10g -jar $PICARD MarkDuplicates \
			VALIDATION_STRINGENCY=SILENT \
			I="${file}" \
			O="${file/.bam/_nodups.bam}" \
			M="${file/.bam/_dupsMarkedStats.txt}" \
			MAX_RECORDS_IN_RAM=2000000 \
			REMOVE_DUPLICATES=true 
done
wait

echo "Duplicates marked: $(date)" >> timelog.txt

for FILE in ./BAM/*_nodups.bam
do
	samtools sort -@ 62 $FILE -O BAM -o "${FILE/_nodups.bam/_sorted_nodups.bam}" &&
	samtools index "${FILE/_nodups.bam/_sorted_nodups.bam}"
done
wait

echo "Duplicates removed: $(date)" >> timelog.txt

for FILE in ./BAM/*_sorted_nodups.bam
do
	bamCoverage --bam "$FILE" \
			--outFileName "${FILE/.bam/.bw}" \
			--outFileFormat bigwig \
			--binSize 5 \
			--numberOfProcessors 62 \
			--normalizeUsing RPGC \
			--effectiveGenomeSize 2913022398 \
			--extendReads \
			--blackListFileName /data/Austin/workdir/genome/hg38/blacklist/UCSC_hg38_blacklist.bed
done
wait

echo "BigWigs made: $(date)" >> timelog.txt

mkdir BW

for FILE in ./BAM/*.bw
	do
	{
		mv $FILE ./BW
	} &
	done
wait


for FILE in ./BAM/*_nodups.bam

	do 
    macs2 callpeak -t "$FILE" \
		-n "$FILE" \
		-f BAMPE \
		-g 2913022398 \
		-q 0.05 \
		--call-summits \
		--nomodel --shift 37 --extsize 73 # ATAC-Seq Shift
done

echo "Peaks called: $(date)" >> timelog.txt

mkdir Peaks

for FILE in ./BAM/*.narrowPeak
	do
	{
		mv $FILE ./Peaks
	} &
	done
wait

for FILE in ./BAM/*.r
	do
	{
		mv $FILE ./Peaks
	} &
	done
wait

for FILE in ./BAM/*.bed
	do
	{
		mv $FILE ./Peaks
	} &
	done
wait

for FILE in ./BAM/*.xls
	do
	{
		mv $FILE ./Peaks
	} &
	done
wait

multiqc .

echo "Complete! $(date)" >> timelog.txt
