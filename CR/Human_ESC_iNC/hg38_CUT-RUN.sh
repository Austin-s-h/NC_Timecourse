# Start with fastq's in a folder named fastq with SE or PE reads with _1 or _2.fastq
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
	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	--minimum-length=25 -j 62 \
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
    echo "beginning alignment to hg38 of file ${FILEPAIR}"
	bowtie2 --local --very-sensitive-local \
		--no-unal --no-mixed --no-discordant \
		--threads 62 \
		-x /data/Austin/workdir/genome/hg38/hg38_bt2/GRCh38 \
		-I 10 -X 1000 \
		-1 "trimmed_${F1}" \
		-2 "trimmed_${F2}" \
		-S ${F1/.fastq/_toGRCh38.SAM}
	echo "complete"

done
wait

cd ../

echo "Reads Aligned: $(date)" >> timelog.txt


for FILE in ./trimmedFastq/*.SAM; do

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
	java -jar $PICARD MarkDuplicates \
			VALIDATION_STRINGENCY=SILENT \
			I="${file}" \
			O="${file/.bam/_dupsMarked.bam}" \
			M="${file/.bam/_dupsMarkedStats.txt}"
done
wait

echo "Duplicates marked: $(date)" >> timelog.txt

for FILE in ./BAM/*_dupsMarked.bam
do
	samtools view -F 1804 -f 2 -b -@ 62 $FILE |
	samtools sort -@ 62 > "${FILE/_dupsMarked.bam/_nodups.bam}" &&
	samtools index "${FILE/_dupsMarked.bam/_nodups.bam}"
done
wait

echo "Duplicates removed: $(date)" >> timelog.txt

for FILE in ./BAM/*_nodups.bam
do
	bamCoverage --bam "$FILE" \
			--outFileName "$FILE".bw \
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
		--call-summits
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
