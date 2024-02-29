# Pipeline483
Pipeline project for COMP483 following **Track 1**: Differential expression in HCMV transcriptomes from patient donors 2- and 6-days post-infection (dpi).

Sample details...

## Dependencies (Install If Needed)
To run the pipeline, users need to have the following programs [and packages] installed:

- python3 [Biopython, pandas]
- R & Rscript [sleuth, dplyr]
- kallisto
- datasets (from NCBI)
- BLAST+

## Installation

#### Clone the GitHub Pipeline483 repository.
From the command line, `cd` to the local directory in which you want to put the cloned repository. Then continue on the command line to create the clone:

```
git clone https://github.com/TaFisc/Pipeline483.git
```

The cloned repository will include the main Python script (wrapper.py) to run the pipeline, supporting scripts, and a directory (sample_data) containing sample test data.

## Run the Code
The entire pipeline runs with a single call to a Python wrapper script that requires two arguments, `--input` and `--email`. The email address entered with `--email` will be used for retreiving records with Entrez. The `--input` argument should include the FULL path to the directory containing the donor transcriptome FASTQ files. (To obtain the full path, `cd` to the directory containing the FASTQ files, `pwd`, and copy the path for use with `--input`.) 

From the main cloned directory, run the following command with insertions specific to your system:
```
time python wrapper.py --input <full path to FASTQ directory> --email <your email>
```

## Overview of Pipeline Steps

#### Step 1: Retrieve transcriptomes

#### Step 2: Build a transcriptome index

#### Step 3: Quantify the TPM

#### Step 4: Determine differentially expressed genes

#### Step 5: Check for other *Betaherpesvirinae* strains with most differentially expressed gene (as protein)