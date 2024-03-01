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

#make a directory for all output files and the project log
out_directory = 'PipelineProject_Taylor_Fischer'
os.system('mkdir ' +out_directory)
#move into new directory
os.chdir(out_directory)

#create the log file in a way that it stays open for multiple calls to write
logfile = open('PipelineProject.log', 'w')

#__________________________________________________
##Step 2: Build Transcriptome Index

#set specific variable
Entrez.email = Entrez_email #email for Entrez retrieval
CDSoutfile = 'transcriptome.faa' #file name for reference transcriptome
transcriptome_ref = 'NC_006273.2' #specific NCBI reference genome

#get reference record in genbank format
#based on biostars: https://www.biostars.org/p/230441/
#create a handle to fetch the reference genome from the nucleotide databse in gb format
handle = Entrez.efetch(db='nucleotide', id=transcriptome_ref, rettype='gb')

#create the output file for the reference transcriptome (CDS)
with open(CDSoutfile, 'w') as CDSf:
        #for each record returned by parsing the handle
        for rec in SeqIO.parse(handle, "genbank"):
                #if the record has features
                if rec.features:
                        for feature in rec.features:
                                #if feature type is CDS
                                if feature.type == "CDS":
                                        #write the gene name and sequence to CDSf
                                        f_name = feature.qualifiers['gene'][0]
                                        f_seq = feature.location.extract(rec).seq
                                        CDSf.write(f'>{f_name}\n{f_seq}\n') #manual FASTA format

#Build index with Kallisto
#create command to send to os.system; -i sets name of index
kallisto_index_CMD = 'kallisto index -i index.idx ' + CDSoutfile
#send to system
os.system(kallisto_index_CMD)

#Calculate the number of CDS in the HCMV genome and output to logfile
#since each record takes 2 lines, number of CDS will be half the number of lines in CDSoutfile

#create line count command
wc_CMD = 'wc -l ' +CDSoutfile
#send to os with ability to save output as a variable
CDS_lines = os.popen(wc_CMD).read() #os.popen() runs a command and captures output as a string
#output is '<number> <filename>', but just want the number
CDS_lines = CDS_lines.split(' ')[0] #split on space and save entry 0
#final count is CDS_lines converted to integer and divided by 2
CDS_count = int(CDS_lines)/2

#write to logfile
logfile.write(f'The HCMV genome ({transcriptome_ref}) has {str(round(CDS_count))} CDS.\n')

#__________________________________________________
##Step 3: Quantify the TPM

#create a directory to hold kallisto results
os.system('mkdir kallisto_results')

#for each sample, create kallisto quant command that calculates TPM and send to os
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

#write header row to logfile, tab-delimited
logfile.write('\nsample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n')

#for each sample, access abundance.tsv file and pull into pd dataframe
for i in range(len(sample_list)): #loop over indices in sample_list
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
        #min_tpm, convert to string
        output.append(str(tpms.min()))
        #med_tpm
        output.append(str(tpms.median()))
        #mean_tpm
        output.append(str(tpms.mean()))
        #max_tpm
        output.append(str(tpms.max()))
        #write to logfile, tab-delimited
        logfile.write('\t'.join(output) + '\n')

#__________________________________________________
##Step 4: Find differentially expressed genes

#prepare sample table with sample metadata and locations for use with sleuth
#columns needed are sample, condition, path
with open('sample_table.txt', 'w') as st:
        #write header row, space-delimited
        st.write(' '.join(['sample', 'condition', 'path']) + '\n')
        for i in range(len(sample_list)): #loop over indices of sample_list
                path = 'kallisto_results/' +sample_list[i] #build path to file
                #write sample info, space-delimited
                st.write(' '.join([sample_list[i], sample_conditions[i], path]) + '\n')

#build command for running the separate 'sleuth.R' script
#script will run sleuth to calculate differentially expressed genes and use dplyr to extract significant transcripts (FDR < 0.05)
sleuth_CMD = 'Rscript ../sleuth.R'
#run sleuth_CMD
os.system(sleuth_CMD)

#prepare output for logfile

#write header row to logfile, tab-delimited
logfile.write('\ntarget_id\ttest_stat\tpval\tqval\n')

#read in file fdr05_results.txt containing significant transcripts (space-separated file)
sig_transcripts = pd.read_csv('fdr05_results.txt', sep=' ')
#subset data to only desired columns for output
sig_transcripts_out = sig_transcripts[['target_id','test_stat','pval','qval']]
#create iterator to pull out each row
for i in range(sig_transcripts_out.shape[0]): #shape[0] is number of rows
        row= sig_transcripts_out.iloc[i, :] #pull out row i, all columns
        #write each entry (converted to str) in the row to logfile, tab-delimited
        logfile.write('\t'.join([str(entry) for entry in row]) + '\n')

#__________________________________________________
##Step 5: Find other Betaherpesvirinae strains with most differentially expressed gene

#retrieve PROTEIN sequence of most differentially expressed GENE

#pull out top gene (gene name) from earlier sleuth results
top_gene_id = sig_transcripts_out.loc[0,'target_id'] #1st entry (row 0) is most differentially expressed

#use Biopython to read in nucl sequence from transcriptome.faa created earlier and convert to protein seq

#create handle 
handle = open(CDSoutfile)
#loop through each individual sequence, and if id is top_gene_id, transcribe to protein & write to fasta file
for record in SeqIO.parse(handle, 'fasta'):
        #check if id is top_gene_id
        if record.id == top_gene_id:
                #translate nucleotide seq to protein, strip final *
                top_prot_seq = record.seq.translate().strip('*')
                #create a new .fasta file to hold top protein seq
                with open('top_protein.fasta', 'w') as fa:
                        #write id and sequence to .fasta
                        fa.write(f'>{record.id}\n{str(top_prot_seq)}\n') #manual fasta format

#taxon of interest
taxon = 'betaherpesvirinae'

#gather Betaherpesvirinae NCBI sequences to build database for blast+ searches
#look in full nr nucleotide database including genome (but no --refseq flag)
gather_db_seqs_CMD = 'datasets download virus genome taxon ' +taxon +' --include genome'
#send to terminal
os.system(gather_db_seqs_CMD)

#need to unzip results from NCBI (ncbi_dataset.zip), so build CMD
unzip_NCBI_CMD = 'unzip ncbi_dataset.zip'
#send to terminal
os.system(unzip_NCBI_CMD)
#creates folders and file ncbi_dataset/data/genomic.fna

#prepare command to build nucleotide database of NUCL sequences, save output in new db folder
make_db_CMD = 'makeblastdb -in ncbi_dataset/data/genomic.fna -out db/' +taxon +' -title ' +taxon +' -dbtype nucl'
#send to terminal
os.system(make_db_CMD)

#run blast+
#since query is PROT and subject is NUCL, use tblastn

#prepare command to run tblastn
#limit searches to 1 hsp per query-subject pair with -max_hsps flag
# -outfmt "6 ..." is tab-delimited (.tsv)
blast_CMD = 'tblastn -query top_protein.fasta -db db/' +taxon +' -max_hsps 1 -out myresults.tsv \
        -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"'
#send to terminal
os.system(blast_CMD)

#prepare output to logfile
#write header row to logfile
logfile.write('\nsacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n')

#read in file 'myresults.csv' containing blast hits
top_blast_results = pd.read_csv('myresults.tsv', header=None, sep='\t') #no header, tab-separated
#create iterator to pull out each row for 1st 10 rows
for i in range(10):
        row= top_blast_results.iloc[i, :] #pull out row i, all columns
        #write each entry (converted to str) in row to logfile, tab-delimited
        logfile.write('\t'.join([str(entry) for entry in row]) + '\n')

#at very end, close the logfile to write contents and save file
logfile.close()