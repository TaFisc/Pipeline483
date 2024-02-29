# Pipeline483
Pipeline project for COMP483 following **Track 1**: Differential expression in HCMV transcriptomes from patient donors 2- and 6-days post-infection (dpi).

## Download the Code and Sample Data


## Install Dependencies (If Not Already Done)
Biopython
kallisto
sleuth

## Run the Code (Using Sample Data)
The entire pipeline runs with a single call to a Python wrapper script that requires two arguments, `--input` and `--email`. The email address entered with `--email` will be used for retreiving records with Entrez. The `--input` argument should include the FULL path to the directory containing the donor transcriptome FASTQ files. (To obtain the full path, `cd` to the directory containing the FASTQ files, `pwd`, and copy the output for use with `--input`.) 
\n\n
From the main cloned directory, run the following command:
```
time python wrapper.py --input <full path to input folder> --email <your email>
```

## Overview of Pipeline Steps

### Step 1: Retrieve Transcriptomes

### Step 2: Build a Transcriptome index

### Step 3: Quantify the TPM

### Step 4: Determine Differentially Expressed Genes

### Step 5: Check for Other *Betaherpesvirinae* Strains with Most Differentially Expressed Protein