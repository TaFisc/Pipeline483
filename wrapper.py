'''Wrapper script to complete COMP483 Pipeline Project, Track 1'''

#import needed packages
import argparse
import sys
import os
from Bio import Entrez
from Bio import SeqIO

#function to parse command line arguments
def check_arg(args=None):
	parser = argparse.ArgumentParser(description='Pipeline Project: Track 1')
        parser.add_argument('-f', '--format',
		help='options "sample" or "full"',
		required='True'
                )
        parser.add_argument('-i', '--input',
		help='path to input file',
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
infile = args.input #need?
outfile = args.output #need?
Entrez_email = args.email

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

#to run kallisto, need to go up 1 directory from PipelineProject_Taylor_Fischer to where sample files are
os.chdir('..')
#create command to send to os.system
kallisto_quant_CMD = 'for sample in *_1.fastq; do time kallisto quant -i '



logfile.close()