import sys, os

mypath = 'PipelineProject_Hannah_Brown'
if os.path.exists(mypath) and os.path.isdir(mypath): # checks if appropriate directory exists and is a directory
    os.chdir(mypath) # moves to the directory
else:
    print("Directory does not exist, must run SRRdown.py prior to wrapper.py")
    sys.exit()

# run shell command, not necessary but helps with debugging as it prints
def run_bash(command): 
    print(f"Running command: {command}")
    os.system(command)

# write to log file
def write_log(name, in_1, in_2): 
    file = 'PipelineProject.log'
    with open(file, 'a') as ifile: # using append instead of write to avoid overwriting content
        #series of if/else to properly write different things to log file
        if name == "spades":
            ifile.write(f'SPAdes command: {in_1}.\n') 
        elif name == "contigs":
            ifile.write(f'There are {in_1} contigs > 1000 bp in the assembly.\n')
            ifile.write(f'There are {in_2} bp in the assembly.\n')
        elif name == "blast":
            ifile.write(f'{in_1}\n')
        elif name == "hits":
            with open(f'{in_1}', 'r') as hold:
                for line in hold:
                    ifile.write(line)
        elif name == "SSR1":
            ifile.write(f'Donor 1 (2dpi) had {in_1} read pairs before Bowtie2 filtering and {in_2} read pairs after.\n')
        else: # will write SSR2 read counts to log
            ifile.write(f'Donor 1 (6dpi) had {in_1} read pairs before Bowtie2 filtering and {in_2} read pairs after.\n')

# count reads from flagstat output
def read_c(file):
    with open(file, 'r') as f:
        line = f.readline().split(' ') #space delimeted split
        return int(line[0]) if line[0].isdigit() else 0 # converts first element to int, returns 0 if not digit

# function for mapping reads, converting, filtering, calling read_c, and write_log
def map(name, fq_1, fq_2, index):
    run_bash(f'bowtie2 -x {index} -1 {fq_1} -2 {fq_2} -S align_{name}.sam') # mapping reads with bowtie2, -x index, -S save to output, -1 FWD, -2 REV
    run_bash(f'samtools view -Sb align_{name}.sam > align_{name}.bam') # convert SAM to BAM (-Sb)
    run_bash(f'samtools view -b -F 4 align_{name}.bam > filtered_align_{name}.bam') # filter out unmapped reads (-F 4), bam file (-b)
    
    run_bash(f'samtools flagstat align_{name}.bam > og_count_{name}.txt') # original read count (unfiltered)
    run_bash(f'samtools flagstat filtered_align_{name}.bam > map_count_{name}.txt') # mapped read count (filtered)
    
    og_read = read_c(f'og_count_{name}.txt') # get read counts
    mapped = read_c(f'map_count_{name}.txt')

    write_log(name, og_read, mapped) # write counts to log file
    # converting bam to fastq using samtools. https://www.metagenomics.wiki/tools/samtools/converting-bam-to-fastq 
    run_bash(f'samtools sort -n filtered_align_{name}.bam -o sorted_{name}.bam ') # sorts the bam file by name
    run_bash(f'samtools fastq -@ 8 sorted_{name}.bam -1 bow_{name}_1.fastq -2 bow_{name}_2.fastq -0 /dev/null -s /dev/null -n') # save fastq reads in paired read files 

# function for processing assembly fasta into contigs and contig-related data
def contigs(fq):
    contigs = [] # initalize list for all contigs
    with open(fq, 'r') as f:
        contig = '' #initalize str for contig seq
        for line in f:
            if line.startswith('>'): # header = new contig
                if contig: # appends contig to list if contig var isnt empty
                    contigs.append(contig)
                contig = '' # empty str for next contig
            else:
                contig += line.strip() # strips and adds to current contig string
        if contig:
            contigs.append(contig) # makes sure all contigs are appended
    contig_1000 = [contig for contig in contigs if len(contig) > 1000] # list of all contigs that are > 1000
    c_1000 = len(contig_1000) # lenth of the list OR how many contigs > 1000
    tot_bp = sum(len(contig) for contig in contig_1000) # sum of the contig lengths for all contigs > 1000
    write_log("contigs", c_1000, tot_bp) # write to log file
    longest_c = max(contigs, key= len) # stores longest contig, finds max by length
    with open("longest_contig.fasta", 'w') as file: # writing the longest contig in fasta format for blast
        file.write(">longest_contig\n")
        file.write(longest_c + '\n')

def main():
    run_bash('>PipelineProject.log') # clear log at beginning of new run
    hcmv = "NC_006273.2" # ncbi accession number so i can use f'{hcmv}
    hcmv_g = "NC_006273.2.fasta" # ncmv fasta file name
    hcmv_i = "HCMV_index"
    # sample fastqs
    ssr1_1 = 'SSR1_1.fastq' # _1 = FWD
    ssr1_2 = 'SSR1_2.fastq' # _2 = REV
    ssr2_1 = 'SSR2_1.fastq' 
    ssr2_2 = 'SSR2_2.fastq'

    # check if fastq files exist - had issues with files being read in 
    for fq_file in [ssr1_1, ssr1_2, ssr2_1, ssr2_2]:
        if not os.path.exists(fq_file): 
            print(f"Error: {fq_file} does not exist.")
            sys.exit()

    run_bash(f'wget -q -O {hcmv_g} "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={hcmv}&rettype=fasta"') # download hcmv genome
    if os.path.exists(hcmv_g): # check if the genome file exists - had issues with the genome downloading
        run_bash(f'bowtie2-build {hcmv_g} {hcmv_i}') # create bowtie2 index

        map('SSR1', ssr1_1, ssr1_2, hcmv_i) # running map function with the paired fastq files and index
        map('SSR2', ssr2_1, ssr2_2, hcmv_i)

        spds_cmd = ('spades.py --rna -1 bow_SSR1_1.fastq -2 bow_SSR1_2.fastq -1 bow_SSR2_1.fastq -2 bow_SSR2_2.fastq -k 99 -o assembly_out') # spades command using bowtie2 output in fastq format
        run_bash(spds_cmd) # running spades command
        write_log("spades", spds_cmd, "") # writing command to log file
        contigs('assembly_out/contigs.fasta') # calling contig function with contigs.fasta from the assembly

        run_bash('datasets download virus genome taxon betaherpesvirinae --refseq --include genome') # download virus genome for betaherepesvirinae from refseq, including genome 
        run_bash('unzip ncbi_dataset.zip') # unzip
        run_bash('mv ncbi_dataset/data/genomic.fna betaherpesvirinae.fna') # move genomic data to betaherpesvirinae.fna
        run_bash('makeblastdb -in betaherpesvirinae.fna -out betaherpesvirinae -title betaherpesvirinae -dbtype nucl') # make blastdb from the genomic data
        write_log("blast", "sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle", "") # write header row, \t is tab character
        # run blast with longest contig fasta against the created db, specify format, save results to a hold file to prevent append issues
        run_bash('blastn -query longest_contig.fasta -db betaherpesvirinae -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" -max_target_seqs 10 -max_hsps 1 -out hold_hits.tsv')
        write_log("hits", 'hold_hits.tsv', '') # send the hold file to write_log funct
    else:
        print(f"Failed to download the genome file {hcmv_g}.")
        sys.exit()
    
if __name__ == "__main__":
    main()

