# lightSra2Count
## Overview
`lightSra2Count.sh` is an HISAT2-StringTie pipeline-based shell script designed for processing RNA-seq data with low machine memory. It automates various steps, including quality control, trimming, alignment, and transcript assembly. The script integrates multiple bioinformatics tools to facilitate a streamlined RNA-seq data analysis pipeline.
## Prerequisites (developer's environment)
Before using lightSra2Count.sh, ensure the following tools should be installed:\

•FastQC v0.12.1 \
•SRA Toolkit v3.0.5 \
•TrimGalore v0.6.10 \
•HISAT2 v2.2.1 \
•StringTie v2.2.1 \
•SAMtools v1.17

## Installation
1. Clone or download the script from the repository.
2. Ensure that all prerequisite software is installed and accessible in empty directory which I prepared for.
3. Prepare for list.txt describing sra accession number as described below:
   [list.txt] \
   SRRXXXXX\
   SRRXXXXX\
   ...

## Usage
Run lightSra2Count.sh from the command line, specifying necessary arguments and options. Detailed usage instructions should be obtained from the script's help message or documentation.
```bash
sh lightSra2Count.sh [OPTIONS]
OPTIONS:
    -h          Display help
    -i VALUE    input sra_file.txt
    -s VALUE    Species (Homo Sapience:hs, Mus Musculus:mm)
    -d VALUE    Direction of sequence (pair end:paired, single end:single)
    -c VALUE	 CPU core
```

## Authors
Noguchi Yuki (Jun Suzuki lab)\
Email: nyuhki21@gmail.com, jsuzuki@icems.kyoto-u.ac.jp

## Citation
If you use lightSra2Count.sh in your research, please cite:\
Noguchi, Y., Onodera, Y., Maruoka, M., Miyamoto, T., Kosako, H., Suzuki, J. 2024. "In vivo CRISPR screening directly targeting testicular cells." Cell Genomics.

## License
This software is released under the MIT License. Please refer to the LICENSE file for more details.

## References
##REFERENCE
1. FastQC:https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ 
2.	SRA toolkit:https://hpc.nih.gov/apps/sratoolkit.html 
3. TrimGalore:https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/ 
4. Pertea, M., Kim, D., Pertea, G.M., Leek, J.T., and Salzberg, S.L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nat Protoc 11, 1650–1667. 10.1038/nprot.2016.095. 
5. Pertea, M., Pertea, G.M., Antonescu, C.M., Chang, T.C., Mendell, J.T., and Salzberg, S.L. (2015). StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nat Biotechnol 33, 290–295. 10.1038/nbt.3122. 
6. Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., and 1000 Genome Project Data Processing Subgroup (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078–2079. 10.1093/bioinformatics/btp352.
