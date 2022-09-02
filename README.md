Publication Release: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7044924.svg)](https://doi.org/10.5281/zenodo.7044924)

# Neural Crest Timecourse Analysis
This repository is broken down into 3 sections
- RNA-Seq Data Processing and Analysis (./RNA-Seq)
- CUT&RUN Data Processing and Analysis (./CUT-RUN)
- ATAC-Seq Data Processing and Analysis (./ATAC-Seq)

In each of these folders, you will find invidual Rmd notebooks for various analyses including...
- Differential Expression / Binding
- Transcription Factor motif enrichment using HOMER and ChromVAR
- hint-ATAC footprinting
- Calculation of fragment length distribution
- Calculation of putative Peak to Gene Enhancer-Promoter linkages (RNA-Seq and ATAC-Seq integration)
- Genomic Overlapping and Clustering

The processing scripts for each genomic sample type are found as individual bash scripts within their respective folders. These contain...
- File Parsing
- Adapter/Quality Trimming
- Alignment
- Duplication Marking and Removal
- Calculation of BigWigs
- Calling Peaks
- Quality Control

If you reference any of these scripts or utilze methods developed from this codebase, please cite
CITATION_HERE
