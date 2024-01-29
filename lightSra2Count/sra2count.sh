#!/usr/sh

###########################################################################
#sra2count/lightSra2Count.sh
#
#	 Copyright (c) 2024, Noguchi Yuki (Jun Suzuki lab)
#	 This software is released under the MIT License, see LICENSE (https://opensource.org/license/mit/).
#    @citation: Noguchi, Y., Onodera, Y., Maruoka, M., Miyamoto, T., Kosako, H., Suzuki., J. 2024. In vivo CRISPR screening directly targeting testicular cells. Cell Genomics.
#    @author:  Noguchi Yuki
#    @contact: nyuhki21@gmail.com,jsuzuki@icems.kyoto-u.ac.jp
#
##REFERENCE
#1. FastQC:https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
#2.	SRA toolkit:https://hpc.nih.gov/apps/sratoolkit.html
#3. TrimGalore:https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
#4. Pertea, M., Kim, D., Pertea, G.M., Leek, J.T., and Salzberg, S.L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nat Protoc 11, 1650–1667. 10.1038/nprot.2016.095.
#5. Pertea, M., Pertea, G.M., Antonescu, C.M., Chang, T.C., Mendell, J.T., and Salzberg, S.L. (2015). StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nat Biotechnol 33, 290–295. 10.1038/nbt.3122.
#6. Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., and 1000 Genome Project Data Processing Subgroup (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078–2079. 10.1093/bioinformatics/btp352.
###########################################################################

function usage {
    cat <<EOM
Usage: $(basename "$0") [OPTION]...
    -h          Display help
    -i VALUE    input sra_file.txt
    -s VALUE    Species (Homo Sapience:hs, Mus Musculus:mm)
    -d VALUE    Direction of sequence (pair end:paired, single end:single)
    -c VALUE	CPU core
EOM
    exit 2
}

while getopts ":i:s:d:c:h" optKey; do
    case "$optKey" in
		i)
          i=$(echo "${OPTARG}")
          ;;
        s)
          s=$(echo "${OPTARG}")
          ;;
        d)
          d=$(echo "${OPTARG}")
          ;;
		c)
          c=$(echo "${OPTARG}")
          ;;
        '-h'|'--help'|* )
          usage
          ;;
    esac
done

echo "####Configuration####"
echo "input_list=$i"
echo "species=$s"
echo "direction=$d"
echo "core=$c"

echo "####Start lightSra2Count...####"

echo "#preprocessing..."
mkdir ./result
LIST=$(cat ${i})
for SampleID in `echo ${LIST}`
    do
		DIR=./result/${SampleID}
		if [ -e $DIR/${SampleID}.sra ]; then
			echo "${SampleID}.sra exists."
		else
			prefetch ${SampleID} -p
		fi
		mv ./${SampleID} ./result
		if [ -e $DIR/${SampleID}.fastq ]; then
			echo "${SampleID}.sra have been converted."
		elif [ -e $DIR/${SampleID}_1.fastq ]; then
			echo "${SampleID}.sra have been converted."
		else
			if [ "${d}" = "paired" ]; then
				fasterq-dump $DIR/${SampleID}.sra --outdir $DIR --split-files --threads ${c} --progress
				echo "${SampleID}.sra have been converted."
			else
				fasterq-dump $DIR/${SampleID}.sra --outdir $DIR --threads ${c} --progress
				echo "${SampleID}.sra have been converted."
			fi
		fi
done

echo "##Trimming..."
for SampleID in `echo ${LIST}`
	do
		DIR=./result/${SampleID}
		mkdir $DIR/${SampleID}_fastqc \
		$DIR/${SampleID}_fastqc/1st \
		$DIR/${SampleID}_fastqc/2nd
		if [ "${d}" = "paired" ]; then
			fastqc --nogroup -o $DIR/${SampleID}_fastqc/1st $DIR/${SampleID}_1.fastq
			fastqc --nogroup -o $DIR/${SampleID}_fastqc/1st $DIR/${SampleID}_2.fastq
			./TrimGalore-0.6.10/trim_galore \
				-j ${c} \
				-q 30 \
				--paired \
				-o $DIR \
				$DIR/${SampleID}_1.fastq \
				$DIR/${SampleID}_2.fastq
			fastqc --nogroup -o $DIR/${SampleID}_fastqc/2nd $DIR/${SampleID}_1_val_1.fq
			fastqc --nogroup -o $DIR/${SampleID}_fastqc/2nd $DIR/${SampleID}_2_val_2.fq
			echo "###Mapping..."
			if [ "${s}" = "hs" ]; then
				./hisat2-2.2.1/hisat2 -q -p ${c} --dta \
					-x ./index/hg38 \
					-1 $DIR/${SampleID}_1_val_1.fq \
					-2 $DIR/${SampleID}_2_val_2.fq \
						| samtools view -@ ${c} -b -  \
						| samtools sort -@ ${c} - \
							> $DIR/${SampleID}.bam
				samtools index $DIR/${SampleID}.bam
				./stringtie/stringtie -p ${c} \
					$DIR/${SampleID}.bam \
					-G ./index/hg38.gtf \
					-o $DIR/${SampleID}.gtf
			else
					./hisat2-2.2.1/hisat2 -q -p ${c} --dta \
					-x ./index/mm10 \
					-1 $DIR/${SampleID}_1_val_1.fq \
					-2 $DIR/${SampleID}_2_val_2.fq \
						| samtools view -@ ${c} -b -  \
						| samtools sort -@ ${c} - \
							> $DIR/${SampleID}.bam
				samtools index $DIR/${SampleID}.bam
				./stringtie/stringtie -p ${c} \
					$DIR/${SampleID}.bam \
					-G ./index/mm10.gtf \
					-o $DIR/${SampleID}.gtf
			fi
		else
			fastqc --nogroup -o $DIR/${SampleID}_fastqc/1st $DIR/${SampleID}.fastq
			./TrimGalore-0.6.10/trim_galore \
				-j ${c} \
				-q 30 \
				-o $DIR \
				$DIR/${SampleID}.fastq
			fastqc --nogroup -o $DIR/${SampleID}_fastqc/2nd $DIR/${SampleID}.fq
			echo "###Mapping..."
			if [ "${s}" = "hs" ]; then
				./hisat2-2.2.1/hisat2 -q -p ${c} --dta \
					-x ./index/hg38 \
					-U $DIR/${SampleID}_trimmed.fq \
						| samtools view -@ ${c} -b -  \
						| samtools sort -@ ${c} - \
							> $DIR/${SampleID}.bam
				samtools index $DIR/${SampleID}.bam
				./stringtie/stringtie -p ${c} \
					$DIR/${SampleID}.bam \
					-G ./index/hg38.gtf \
					-o $DIR/${SampleID}.gtf
			else
				./hisat2-2.2.1/hisat2 -q -p ${c} --dta \
					-x ./index/mm10 \
					-U $DIR/${SampleID}_trimmed.fq \
						| samtools view -@ ${c} -b -  \
						| samtools sort -@ ${c} - \
							> $DIR/${SampleID}.bam
				samtools index $DIR/${SampleID}.bam
				./stringtie/stringtie -p ${c} \
					$DIR/${SampleID}.bam \
					-G ./index/mm10.gtf \
					-o $DIR/${SampleID}.gtf
			fi
		fi
done

echo "####Counting..."
#preprocess for ballgown
mkdir ./summary \
	  ./result/ballgown
for SampleID in `echo ${LIST}`
	do
		DIR=./result/${SampleID}
		mkdir ./result/ballgown/${SampleID}
		echo $DIR/${SampleID}.gtf \
			>> merge_list.txt
		echo ${SampleID} ./result/ballgown/${SampleID}/${SampleID}.bg.gtf \
			>> ballgown_PATH.txt
done

if [ "${s}" = "hs" ]; then
	./stringtie/stringtie --merge -p ${c} \
		-G ./index/hg38.gtf \
		-o ./summary/merge.gtf \
		merge_list.txt
else
	./stringtie/stringtie --merge -p ${c} \
		-G ./index/mm10.gtf \
		-o ./summary/merge.gtf \
		merge_list.txt
fi

#generation of count matrix
for SampleID in `echo ${LIST}`
	do
		DIR=./result/${SampleID}
		./stringtie/stringtie -p ${c} -e -B \
			$DIR/${SampleID}.bam \
			-G ./summary/merge.gtf \
			-o ./result/ballgown/${SampleID}/${SampleID}.bg.gtf \
			-A ./summary/${SampleID}.tab
done

python ./stringtie/prepDE.py3 -i ballgown_PATH.txt
mv gene_count_matrix.csv ./result/gene_raw_count.csv
mv ./result/transcript_count_matrix.csv ./result/transcript_raw_count.csv

#estimation of annotation accuracy
if [ "${s}" = "hs" ]; then
	./gffcompare/gffcompare \
		-r ./index/hg38.gtf \
		-G \
		-o ./summary/gffcompare ./summary/merge.gtf
else
	./gffcompare/gffcompare \
		-r ./index/mm10.gtf \
		-G \
		-o ./summary/gffcompare ./summary/merge.gtf
fi

rm merge_list.txt ballgown_PATH.txt
mv ./summary ./result
mv sra_file.txt ./result

echo "#####done."
exit=0
