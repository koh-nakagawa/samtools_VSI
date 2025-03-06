# samtools_VSI
samtools view, sort, index pipeline

This script automates the conversion of a SAM file to BAM format, sorting, and indexing using samtools. The results are stored in a subdirectory named YYMMDD_samtools_VSI_results (where YYMMDD represents the current date) within the specified output directory. In addition, the executed commands are logged in a file named CMD.txt within the results directory.This script automates the conversion of a SAM file to BAM format, sorting, and indexing using samtools. The results are stored in a subdirectory named YYMMDD_samtools_VSI_results (where YYMMDD represents the current date) within the specified output directory. In addition, the executed commands are logged in a file named CMD.txt within the results directory.

## Requirements
- Python 3.X
- samtools
If using a conda environment, you can install the required tools with the following command:
```
$ conda create -n asm_env -c bioconda minimap2 samtools
```
