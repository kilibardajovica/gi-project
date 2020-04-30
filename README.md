# gi-project
Illumina paired-end read sequencing simulator and BWA-MEM alignment quality evaluation

# Description

## Simulate function parameters:
  - refGenomeFile : name of the FASTA file where reference genome is located
  - avgQuality : average nucleotid quality (value between 0 and 40)
  - coverage : percentage of the reference genome covered with reads
  - readLength : single read size
  - insertSize : single fragment size ( 2 x readLength + inner distance between reads)
  - errRateSNV : error rate for single nucleotid variation
  - errRateInsert : error rate for insertion mutation
  - errRateDelete : error rate for deletion mutation
 
## EvaluateBWAMEM function parameter:
  -  fileName : name of the input file without extension


# Instructions

There are two ways to run simulator:
  - Execute python script 'script.py' in terminal with appropriate arguments
   
    e.g. python script.py -g refGenomeFile -c coverage -q avgQuality -r readLength -f insertSize -s errRateSNV -i errRateInsert -d errRateDelete
    
  - Using Jupyter Notebook:
  
    - Run simulate function. As output two FASTQ files and SAM file are created. 

      e.g. simulate('Test.fa',25,1,4,10,0,0.001,0.001)
    - Execute BWA-MEM in terminal (If it is the first time FASTA file is used it is necessary to create index)

      e.g. (optional) bwa index Test.fa

       bwa mem Test.fa Test_1.fastq Test_2.fastq > Test_bwa.sam

    - Run evaluateBWAMEM function

      e.g. evaluateBWAMEM('Test')
    
## Video presentation

  https://www.youtube.com/watch?v=TgmBEx5Wjw8
