# **Python Pipeline Project**

## *Description*<br />
This project is a pipeline analysis of transcriptomes from one patient donor 2- and 6-days post infection (dpi) from SRA. <br />

## *Useage*<br />
This project consists of 2 .py files: <br />
1) SRRdown.py- **Requires: SRA Toolkit.** This file creates the pipeline directory. Please see ***Testing Instructions*** prior to running SRRdown.py.<br />
Downloads, converts to paired-end fastq, and renames SRR files from a text file containing the SRR numbers in the following format:

   >**SRR5660030**<br />
   >**SRR5660033**<br />

    A text file named *Input_SRR* is provided, containing the SRR numbers in proper format. <br />
    To execute SRRdown.py run the following command from the command line or terminal: <br />
    >python SRRdown.py -i your_SRR_file <br />
2) wrapper.py- **Requires: Biopython, ncbi-datasets-cli, unzip, bowtie2, samtools, SPAdes.py, and blast+.** The primary code file handeling all pipeline analysis. <br />
For ***non-testing*** use: with SRRdown.py and wrapper.py in the same location, run SRRdown.py first. SRRdown.py will create the directory and move to it. wrapper.py will move to the directory automatically.<br />
After running SRRdown.py, wrapper.py can be run with the following command in the command line or terminal:
    >**python wrapper.py** <br />
<br />

**Outputs are written to \PipelineProject_Hannah_Brown\PipelineProject.log** <br />

### *Testing Instructions* and Sample Data <br />
For testing, please download the following files from the respository: <br />
>**SSR1_1.fastq** <br />
>**SSR1_2.fastq** <br />
>**SSR2_1.fastq** <br />
>**SSR2_2.fastq** <br /> 

  ***DO NOT RUN SRRdown.py!*** <br />
  Instead of running SRRdown.py to create a directory, it is easiest to create the directory in the command line or terminal. This avoids the full file download, and does not require any renaming or moving of the full files. <br />
  *If SRRdown.py is run and you would like to use the sample data, the full data and sample data will have the same file name. You can either move the full data out of the directory and move the sample data into the directory, or you can rename the sample data, move it into the directory, and edit wrapper.py lines 88-91 to the correct file name.*<br />
  To create a directory, run the following code in the command line or terminal: <br />
>**mkdir PipelineProject_Hannah_Brown** <br />

  To move the sample files into the directory, run the following command in the command line or terminal: <br />
>**mv SSR1_1.fastq PipelineProject_Hannah_Brown/** <br />
>**mv SSR1_2.fastq PipelineProject_Hannah_Brown/** <br />
>**mv SSR2_1.fastq PipelineProject_Hannah_Brown/** <br />
>**mv SSR2_2.fastq PipelineProject_Hannah_Brown/** <br />

  The sample data can now be tested in wrapper.py by executing the file with the follwoing command in the command line or terminal: <br />
>**python wrapper.py** <br />

  **Outputs are written to \PipelineProject_Hannah_Brown\PipelineProject.log** 


 
