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

## Usage
The script requires the following command line arguments:
- `-i`: Path to the input SAM file
- `-o`: Path to the output directory (this directory must exist)
- `-@`: Number of threads to use (default is 1, Optional)
For example, run the script as follows:
```
$ python samtools_VSI.py -i /path/to/input.sam -o /path/to/output -@ 4
```
In this example, a subdirectory named YYMMDD_samtools_VSI_results (with today's date in the format YYMMDD) will be created inside /path/to/output. The subdirectory will contain the following files:
- `alignment.bam`
- `alignment.sorted.bam`
- `alignment.sorted.bam.bai`
- `CMD.txt` (recording the executed commands)

## Notes
- The output directory must exist; otherwise, the script will exit with an error.
- If the results subdirectory already exists, the script will not overwrite it and will exit with an error.
- All error messages are displayed in English.

## Script Overview
Input Check
Verifies that the input SAM file and the output directory exist.

Result Directory Creation
Creates a subdirectory named YYMMDD_samtools_VSI_results inside the specified output directory.

SAM to BAM Conversion, BAM Sorting, and Indexing
Executes the samtools commands sequentially and saves the outputs in the results subdirectory.
Each command is logged into CMD.txt prior to execution.

Error Handling
If any step fails, the script exits with an appropriate error message in English.
