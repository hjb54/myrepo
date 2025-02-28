import sys
import argparse
import os
from Bio import SeqIO

# Function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description="Download, convert to fastq, and rename SRR files.")
    parser.add_argument("-i", "--input", help="input file", required=True)
    parser.add_argument("-o", "--output", help="output file", required=True)
    return parser.parse_args(args)

# Retrieve command line arguments
arguments = check_arg(sys.argv[1:])
infile = arguments.input
outfile = arguments.output

s1 = "" #empty str vars for the SRR names
s2 = ""
with open (infile, 'r') as file:
    data = file.readlines()
    s1 = data[0].strip()
    s2 = data[1].strip()

mypath = 'PipelineProject_Hannah_Brown'
if os.path.exists(mypath) and os.path.isdir(mypath): #checks if appropriate directory exists and is a directory
    os.chdir('PipelineProject_Hannah_Brown') #moves to the directory
else:
    os.system('mkdir PipelineProject_Hannah_Brown') #makes the apropriate directory
    os.chdir('PipelineProject_Hannah_Brown') #moves to the directory

def getfile(srr):
    hold_download = ('wget https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc='+srr) #holds the download command
    os.system(hold_download) #runs the download command in the terminal
    name = ('mv fastq?acc='+(srr)+' '+ srr) # renaming file
    os.system(name)
    split = ('fasterq-dump --split-files '+ srr) #split the files, run in terminal
    os.system(split)
    return()

getfile(s1)
getfile(s2)
os.system('mv' + s1 + '_1.fastq' + 'SSR1_1.fastq') # renaming with a simpler notation, prevents silly spelling errors and allows for quick recognition
os.system('mv' + s1 + '_2.fastq' + 'SSR1_2.fastq')
os.system('mv' + s2 + '_1.fastq' + 'SSR2_1.fastq')
os.system('mv' + s2 + '_2.fastq' + 'SSR2_2.fastq')

print('SRRdown.py has run.\nFull fastq files are downloaded and directory has been created.\nContinue to wrapper.py.')