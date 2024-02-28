'''Wrapper script to complete COMP483 Pipeline Project, Track 1'''

#import needed packages
import argparse
import sys
import os
from Bio import Entrez
from Bio import SeqIO
import pandas as pd #need?

#function to parse command line arguments
def check_arg(args=None):
	parser = argparse.ArgumentParser(description='Pipeline Project: Track 1')
        parser.add_argument('-f', '--format',
		help='options "sample" or "full"',
		required='True'
                )
        parser.add_argument('-i', '--input',
		help='full path to folder with transcript files',
		required='True'
                )
        parser.add_argument('-o', '--output',
		help='output file name',
		required='True'
		)
        parser.add_argument('-e', '--email',
                help='email for Entrez',
                required='True'
                )
	return parser.parse_args(args)

#retrieve command line arguments and assign to variables
args = check_arg(sys.argv[1:])
form = args.format #need?
infolder = args.input
outfile = args.output #need?
Entrez_email = args.email

#temporary args
infolder = '/home/tfischer1/Pipeline483/sample_data'
Entrez_email = ''

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
                                        nfh.write(">%s from %s\n%s\n" % (
                                        feature.qualifiers['gene'][0],
                                        rec.name,
                                        feature.location.extract(rec).seq))

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

#first, store sample-relevant information

#create list of samples
sample_list = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']
#create corresponding list of conditions
sample_conditions = ['2dpi', '6dpi', '2dpi', '6dpi']

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



logfile.close()