'''Wrapper script to complete COMP483 Pipeline Project, Track 1'''

#import needed packages
import argparse
import sys
import os
from Bio import Entrez
from Bio import SeqIO
import pandas as pd

#function to parse command line arguments
def check_arg(args=None):
        parser = argparse.ArgumentParser(description='Pipeline Project: Track 1')
        parser.add_argument('-i', '--input', help='full path to folder with transcript files', required='True')
        parser.add_argument('-e', '--email', help='email for Entrez', required='True')
        return parser.parse_args(args)

#retrieve command line arguments and assign to variables
args = check_arg(sys.argv[1:])
infolder = args.input
Entrez_email = args.email

#store sample-relevant information
#create list of samples
sample_list = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']
#create corresponding list of conditions
sample_conditions = ['2dpi', '6dpi', '2dpi', '6dpi']

#make a directory for all output files and project log
out_directory = 'PipelineProject_Taylor_Fischer'
os.system('mkdir ' +out_directory)
#move into new directory
os.chdir(out_directory)

logfile = open('PipelineProject.log', 'w')

#__________________________________________________
##Step 2: Build Transcriptome Index

#get reference record in genbank format
Entrez.email = Entrez_email
CDSoutfile = 'transcriptome.faa'
transcriptome_ref = 'NC_006273.2'

#based on biostars: https://www.biostars.org/p/230441/
handle = Entrez.efetch(db='nucleotide', id=transcriptome_ref, rettype='gb')
with open(CDSoutfile, 'w') as nfh:
        for rec in SeqIO.parse(handle, "genbank"):
                if rec.features:
                        for feature in rec.features:
                                if feature.type == "CDS":
                                        f_name = feature.qualifiers['gene'][0]
                                        f_seq = feature.location.extract(rec).seq
                                        nfh.write(f'>{f_name}\n{f_seq}\n')

#Build index with Kallisto
#create command to send to os.system
kallisto_index_CMD = 'kallisto index -i index.idx ' + CDSoutfile
#send to system
os.system(kallisto_index_CMD)

#Calculate the number of CDS in the HCMV genome and output to logfile
#number of CDS will be half the number of lines in CDSoutfile
wc_CMD = 'wc -l ' +CDSoutfile
CDS_lines = os.popen(wc_CMD).read() #os.popen() runs a command and captures output as a string
#output is '<number> <filename>', but just want the number
CDS_lines = CDS_lines.split(' ')[0]
CDS_count = int(CDS_lines)/2

#write to logfile
logfile.write(f'The HCMV genome ({transcriptome_ref}) has {str(round(CDS_count))} CDS.\n')

#__________________________________________________
##Step 3: Quantify the TPM

#create a directory to hold kallisto results
os.system('mkdir kallisto_results')

#for each sample, create kallisto quant command that calculated TPM and run
for sample in sample_list:
        #build path to read1 fastq file
        fq_path1 = infolder + '/' + sample + '_1.fastq'
        #build path to read2 fastq file
        fq_path2 = infolder + '/' + sample + '_2.fastq'
        #build kallisto quant command (30 bootstraps, 2 threads)
        quant_CMD = 'time kallisto quant -i index.idx -o kallisto_results/' \
                +sample +' -b 30 -t 2 ' +fq_path1 + ' ' +fq_path2
        #run command
        os.system(quant_CMD)
        
#gather information for logfile

#write header row to logfile
logfile.write('\nsample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n')

#for each sample, access abundance.tsv file and pull into pd dataframe for access
for i in range(len(sample_list)):
        #build file path
        abund_file = 'kallisto_results/' + sample_list[i] + '/abundance.tsv'
        #read in abundance table, a tab-separated file
        abund_table = pd.read_csv(abund_file, sep='\t')
        #pull out tpm column into a list
        tpms = abund_table['tpm']
        #create list for sample output
        output = []
        #sample ID
        output.append(sample_list[i])
        #condition
        output.append(sample_conditions[i])
        #min_tpm
        output.append(str(tpms.min()))
        #med_tpm
        output.append(str(tpms.median()))
        #mean_tpm
        output.append(str(tpms.mean()))
        #max_tpm
        output.append(str(tpms.max()))
        #write to logfile
        logfile.write('\t'.join(output) + '\n')

#__________________________________________________
##Step 4: Find differentially expressed genes

#prepare sample table with sample metadata and locations for use with sleuth
with open('sample_table.txt', 'w') as st:
        st.write(' '.join(['sample', 'condition', 'path']) + '\n')
        for i in range(len(sample_list)):
                path = 'kallisto_results/' +sample_list[i]
                st.write(' '.join([sample_list[i], sample_conditions[i], path]) + '\n')

#build command for running the separate 'sleuth.R' script
#script will run sleuth to calculate differentially expressed genes and use dplyr to extract significant transcripts (FDR < 0.05)
sleuth_CMD = 'Rscript ../sleuth.R'
#run sleuth_CMD
os.system(sleuth_CMD)

#prepare output for logfile

#write header row to logfile
logfile.write('\ntarget_id\ttest_stat\tpval\tqval\n')

#read in file fdr05_results.txt containing significant transcripts (space-separated file)
sig_transcripts = pd.read_csv('fdr05_results.txt', sep=' ')
#pull desired output columns into lists for easy iteration
#subset data to only desired columns for output
sig_transcripts_out = sig_transcripts[['target_id','test_stat','pval','qval']]
#create iterator to pull out each row
for i in range(sig_transcripts_out.shape[0]):
        row= sig_transcripts_out.iloc[i, :]
        logfile.write('\t'.join([str(entry) for entry in row]) + '\n')

#__________________________________________________
##Step 5: Find other Betaherpesvirinae strains with most differentially expressed gene

#retrieve PROTEIN sequence of most differentially GENE

#pull out top gene (gene name) from earlier sleuth results
top_gene_id = sig_transcripts_out.loc[0,'target_id'] #1st row in the dataframe

#use Biopython to read in nucl sequence from transcriptome.faa and convert to protein
handle = open(CDSoutfile)

#loop through each individual sequence, and if it is top_gene, transcribe to protein & write to fasta file
for record in SeqIO.parse(handle, 'fasta'):
        if record.id == top_gene_id:
                top_prot_seq = record.seq.translate().strip('*')
                with open('top_protein.fasta', 'w') as fa:
                        fa.write(f'>{record.id}\n{str(top_prot_seq)}\n')

#taxon of interest
taxon = 'betaherpesvirinae'

#gather Betaherpesvirinae NCBI sequences to build database for blast+ searches

gather_db_seqs_CMD = 'datasets download virus genome taxon ' +taxon +' --refseq --include genome'
#send to terminal
os.system(gather_db_seqs_CMD)

#need to unzip results from NCBI, so build CMD
unzip_NCBI_CMD = 'unzip ncbi_dataset.zip'
#send to terminal
os.system(unzip_NCBI_CMD)

#prepare command to build database of NUCL sequences
make_db_CMD = 'makeblastdb -in ncbi_dataset/data/genomic.fna -out db/' +taxon +' -title ' +taxon +' -dbtype nucl'
#send to terminal
os.system(make_db_CMD)

#run blast+
#since query is PROT and subject is NUCL, use tblastn

#prepare command to run tblastn
#limit searches to 1 hsp per query-subject pair with -max_hsps flag
blast_CMD = 'tblastn -query top_protein.fasta -db db/' +taxon +' -max_hsps 1 -out myresults.csv \
        -outfmt "10 sacc pident length qstart qend sstart send bitscore evalue stitle"'
#send to terminal
os.system(blast_CMD)

#prepare output to logfile
#write header row to logfile
logfile.write('\nsacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n')

#read in 1st 10 rows (top 10 hits) of file 'myresults.csv' containing blast hits
top_blast_results = pd.read_csv('myresults.csv', nrows=10, header=None)
#create iterator to pull out each row
for i in range(top_blast_results.shape[0]):
        row= top_blast_results.iloc[i, :]
        logfile.write('\t'.join([str(entry) for entry in row]) + '\n')

#at very end, close the logfile to write contents and save file
logfile.close()