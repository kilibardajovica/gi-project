import random as rand
import os
import numpy as np
import subprocess
import sys, getopt

nucleotids = ['G', 'C', 'A', 'T']
complementary = {
  "C": "G",
  "G": "C",
  "T": "A", 
  "A": "T"
}

#Check if simulator input parameters are valid
def checkArguments(avgQuality, coverage, readLength, insertSize, errRateSNV, errRateInsert, errRateDelete):
    if(avgQuality < 0 or avgQuality > 40):
        print("Average quaility value must be between 0 and 40!")
        return False
    if(readLength*2 > insertSize):
        print("Insert size must be at least two times read length!")
        return False
    if(errRateInsert < 0 or errRateInsert > 1):
        print("Insetion error rate is probability value and must be between 0 and 1!")
        return False
    if(errRateDelete < 0 or errRateDelete > 1):
        print("Deletion error rate is probability value and must be between 0 and 1!")
        return False
    if(errRateSNV < 0 or errRateSNV > 1):
        print("SNV error rate is probability value and must be between 0 and 1!")
        return False
    return True

#Read all sequences from fasta file
def parseFastaFile(file):
    name, seq = None, []
    for line in file:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

#Parse reference genome
def parseGenomeFile(fileName):
    refGenome = {}
    with open(fileName) as file:
        for name, seq in parseFastaFile(file):
            refGenome[name[1:]] = list(seq)
    return refGenome

DeletionMutations = {}
InsertionMutations = {}
SNVMutations = {}

def insertMutations(refGenome, errRateSNV, errRateInsert, errRateDelete):
    totalMutations = 0
    for name,sequence in refGenome.items():            
        SNVnum = round(len(sequence)*errRateSNV)
        DeleteNum = round(len(sequence)*errRateDelete)
        InsertNum = round(len(sequence)*errRateInsert)
        
        totalMutations += SNVnum
        totalMutations += DeleteNum
        totalMutations += InsertNum
        
        delList = []
        insertList = []
        snvList = []
        
            
        #Deletion mutations
        for i in range(DeleteNum):
            mutationPosition = rand.randint(0,len(sequence)-1)
            del sequence[mutationPosition]
            delList.append(mutationPosition)
        
        #Insertion mutations
        for i in range(InsertNum):
            mutationPosition = rand.randint(0,len(sequence)-1)
            insertionNucleotid = rand.randint(0,3)
            sequence.insert(mutationPosition,nucleotids[insertionNucleotid])
            insertList.append(mutationPosition)
            for deletion in delList:
                if(deletion >= mutationPosition):
                    deletion += 1
            
        #SNV mutations
        for i in range(SNVnum):
            mutationPosition = rand.randint(0,len(sequence)-1)
            mutationNucleotid = rand.randint(0,3)
            while(sequence[mutationPosition] == nucleotids[mutationNucleotid]):
                mutationNucleotid = rand.randint(0,3)
            sequence[mutationPosition] = nucleotids[mutationNucleotid]
            snvList.append(mutationPosition)
            
        DeletionMutations[name] = delList
        InsertionMutations[name] = insertList
        SNVMutations[name] = snvList
    
    return totalMutations
                                
#generate nucleotids qualities using normal distribution
def generateQuality(avgQuality, readLength):
    quality = []
    randomNormal = np.random.normal(avgQuality,5,readLength)
    for i in range(len(randomNormal)):
        if(randomNormal[i] < 0):
            randomNormal[i] = 0
        if(randomNormal[i] > 40):
            randomNormal[i] = 40
        quality.append(chr(int(round(randomNormal[i])) + 33))
    return "".join(quality)

#calculate position of the base in reference genome
def getOriginalPosition(position,name):
    insertCnt = 0
    deleteCnt = 0
    for insertion in InsertionMutations[name]:
        if(insertion < position):
            insertCnt += 1
    for deletion in DeletionMutations[name]:
        if(deletion < position):
            deleteCnt += 1
    return position - insertCnt + deleteCnt

def sequenceReads(refGenomeFile,refGenome,avgQuality,coverage,readLength,insertSize):
    #get file name
    fileNameWithExtension = os.path.basename(refGenomeFile)
    fileName = os.path.splitext(fileNameWithExtension)[0]
    #create FASTQ and SAM files
    samFile = open(fileName + "_sim.sam","w")
    firstReadFile = open(fileName + "_1.fastq","w")
    secondReadFile = open(fileName + "_2.fastq","w")
    
    totalReads = 0
    
    for name,sequence in refGenome.items():
        pairedReadsNum = round((coverage*len(sequence))/(2*readLength))
        totalReads += pairedReadsNum
        for i in range(pairedReadsNum):
            insertPosition = rand.randint(0,len(sequence)-insertSize)
            #adjust insert position so neither of reads start with mutated base
            while((insertPosition + insertSize < len(sequence))):
                if((insertPosition in SNVMutations[name]) or ((insertPosition + insertSize - readLength) in SNVMutations[name])):
                    insertPosition += 1
                elif ((insertPosition in InsertionMutations[name]) or ((insertPosition + insertSize - readLength) in InsertionMutations[name])):
                    insertPosition += 1
                else:
                    break
            
            seqId = "read_" + str(i) + ":" + str(insertPosition + 1) + ":" + str(insertPosition + insertSize + 1)
            
            #sequence first read and write it to files
            firstRead = []
            for j in range(readLength):
                firstRead.append(sequence[insertPosition + j])
            quality = generateQuality(avgQuality, readLength)
            firstReadFile.write("@" + str(seqId) + "/1\n" + "".join(firstRead) + "\n+\n" + quality + "\n")
            tab = "\t"
            pos = getOriginalPosition(insertPosition,name)
            samFile.write(seqId+tab+str(pos+1)+tab+"".join(firstRead)+tab+quality+"\n")
            
            #sequence second read and write it to files (it is inverse and complementary)
            secondRead = []
            originalSecond = []
            for j in range(readLength):
                if(sequence[insertPosition + insertSize - j - 1] in nucleotids):
                    secondRead.append(complementary[sequence[insertPosition + insertSize - j - 1]])
                else:
                    secondRead.append(sequence[insertPosition + insertSize - j - 1])
                originalSecond.append(sequence[insertPosition + insertSize - readLength + j])
            quality = generateQuality(avgQuality, readLength)
            secondReadFile.write("@" + str(seqId) + "/2\n" + "".join(secondRead) + "\n+\n" + quality + "\n")
            pos = getOriginalPosition(insertPosition+insertSize-readLength,name)
            samFile.write(seqId+tab+str(pos+1)+tab+"".join(originalSecond)+tab+quality+"\n")
            
    samFile.close()
    firstReadFile.close()
    secondReadFile.close()
    return totalReads

def simulate(refGenomeFile, avgQuality, coverage, readLength, insertSize, errRateSNV, errRateInsert, errRateDelete):
    if(checkArguments(avgQuality, coverage, readLength, insertSize, errRateSNV, errRateInsert, errRateDelete)):
        refGenome = parseGenomeFile(refGenomeFile)
        insertMutations(refGenome,errRateSNV,errRateInsert,errRateDelete)
        sequenceReads(refGenomeFile,refGenome,avgQuality,coverage,readLength,insertSize)
        return True
    else:
        return False

def evaluateBWAMEM(fileName):
    simulatedSam = open(fileName+"_sim.sam", "r")
    bwamemSam = open(fileName+"_bwa.sam", "r")
    alignments = {}
    aligned = 0
    for read in simulatedSam.readlines():
        split = read.split("\t")
        alignments[split[0]] = split[1]
    for read in bwamemSam.readlines():
        if(read[0] == "@"):
            continue
        split = read.split("\t")
        if(split[3] == alignments[split[0]]):
            aligned += 1
    simulatedSam.close()
    bwamemSam.close()
    print("Correctly aligned "+str(aligned)+" of "+str(len(alignments))+" reads, rate "+str(aligned/len(alignments)*100)+"%")


def execute(argv):
    refGenomeFile = ''
    coverage = 1
    quality = 40
    readSize = 150
    insertSize = 500
    errSNV = 0
    errInsert = 0
    errDelete = 0

    try:
        opts, args = getopt.getopt(argv,"hg:c:q:r:f:s:i:d:",["ifile=","coverage=","quality=","read=","fragment=","snv=","insert=","delete="])
    except getopt.GetoptError:
        print('script.py -g <regGenomeFile> -c <coverage> -q <quality> -r <readSize> -f <insertSize> -s <errRateSNV> -i <errRateInsertion> -d <errRateDeletion>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('script.py -g <regGenomeFile> -c <coverage> -q <quality> -r <readSize> -f <insertSize> -s <errRateSNV> -i <errRateInsertion> -d <errRateDeletion>')
            sys.exit()
        elif opt in ("-g", "--ifile"):
            refGenomeFile = arg
        elif opt in ("-c", "--coverage"):
            coverage = int(arg)
        elif opt in ("-q", "--quality"):
            quality = int(arg)
        elif opt in ("-r", "--read"):
            readSize = int(arg)
        elif opt in ("-f", "--fragment"):
            insertSize = int(arg)
        elif opt in ("-s", "--snv"):
            errSNV = int(arg)
        elif opt in ("-i", "--insert"):
            errInsert = int(arg)
        elif opt in ("-d", "--delete"):
            errDelete = int(arg)

    success = simulate(refGenomeFile,quality,coverage,readSize,insertSize,errSNV,errInsert,errDelete)
    if(success == False):
        sys.exit(2)
    fileNameWithExtension = os.path.basename(refGenomeFile)
    fileName = os.path.splitext(fileNameWithExtension)[0]
    subprocess.run(["bwa", "index", refGenomeFile])
    subprocess.run(["bwa", "mem", refGenomeFile, fileName + "_1.fastq", fileName + "_2.fastq"], stdout=open(fileName + "_bwa.sam", "w"))
    evaluateBWAMEM(fileName)

execute(sys.argv[1:])
