# Pipeline483 Overview

#### Pipeline project for COMP483 following **Track 1**: Differential expression in HCMV transcriptomes from patient donors 2- and 6-days post-infection (dpi).

This pipeline is built specifically for use with the following samples from Cheng et al. 2017 (<https://www.ncbi.nlm.nih.gov/pubmed/29158406>):   
Donor 1 (2dpi): <https://www.ncbi.nlm.nih.gov/sra/SRX2896360>  
Donor 1 (6dpi): <https://www.ncbi.nlm.nih.gov/sra/SRX2896363>  
Donor 3 (2dpi): <https://www.ncbi.nlm.nih.gov/sra/SRX2896374>  
Donor 3 (6dpi): <https://www.ncbi.nlm.nih.gov/sra/SRX2896375>  

# Instructions for Installation and Use

## Dependencies (Install If Needed)
To run the pipeline, users need to have the following programs [and packages] installed:

- python3 [Biopython, pandas]
- R & Rscript [sleuth, dplyr]
- kallisto
- datasets (from NCBI)
- BLAST+

## Clone the GitHub Pipeline483 repository.
From the command line, `cd` to the directory in which you want to put the cloned repository. Then continue on the command line to create the clone:

```
git clone https://github.com/TaFisc/Pipeline483.git
```

The cloned repository will include the main Python script (wrapper.py) for running the pipeline, supporting scripts, and a directory (sample_data) containing sample test data.

## Run the Code
The entire pipeline runs with a single call to a Python wrapper script that requires two arguments, `--input` and `--email`. The email address entered with `--email` will be used for retreiving records with Entrez. The `--input` argument should include the **FULL** path to the directory containing the donor transcriptome FASTQ files. 

**For general use**, run the following from the main cloned directory `Pipeline483` with insertions specific to user system:
```
time python wrapper.py --input <full path to FASTQ directory> --email <your email>
```
**For use with the sample data**, assuming the `git clone` command was just completed, perform the following steps:
```
cd Pipeline483/sample_data
pwd
```
Copy the output file path (the FULL path) for use with `--input` argument.
```
cd ..
```
Run the following command with the copied path and your email.
```
time python wrapper.py --input <copied path to sample_data directory> --email <your email>
```

# Overview of Pipeline Steps
*Step 1 is not part of the automated pipeline. Details of how it was carried out are included below, as well as the steps to create the sample data. Steps 2-5 are automated by running the main script wrapper.py.*

#### Step 1: Retrieve transcriptomes
Step 1.1: Download the transcriptomes from NCBI's SRA to a new directory.
```
#make a directory to hold the transcriptome files
mkdir HCMV_transcriptomes
cd HCMV_transcriptomes

#download the files

#for Donor 1 (2dpi) SRR5660030
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030

#for Donor 1 (6dpi) SRR5660033
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033

#for Donor 3 (2dpi) SRR5660044
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044

#for Donor 3 (6dpi) SRR5660045
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045
```
Step 1.2: Convert SRA downloads to paired end fastq files.
```
#use fastq-dump to uncompress the files
#-I appends _1 and _2 after file names
#--split-files breaks paired-end reads into 2 files

for sample in SRR*;
do
	fastq-dump -I --split-files $sample;
done	
```
Step 1.3: Subsample each FASTQ to 10000 reads to create sample test data.
```
#make a directory to hold sample data
mkdir sample_data/

#for each FASTQ file, write the first 40,000 lines to a new file in sample_data
for file in *.fastq;
do
	head -n 40000 $file > sample_data/$file;
done
```
The directory `sample_data` was then moved to the repository `Pipeline483` using `mv` for use as test data.

#### Step 2: Build a transcriptome index

#### Step 3: Quantify the TPM

#### Step 4: Determine differentially expressed genes

#### Step 5: Check for other *Betaherpesvirinae* strains with most differentially expressed gene (as protein)