import concurrent.futures
import csv
from datetime import date
import multiprocessing as mp
from operator import itemgetter
import os
import pandas as pd
import random
import sys
import time
from virus_read_detection_blastn_search import (BLASTn_CHO_DB_search,
                                                BLASTn_RVDB_DB_search, 
                                                BLASTn_RVDB_DB_search_result_parser)

# Global
reads= []
COUNT = 0


# mm/dd/y
today = date.today().strftime("%b-%d-%Y")


# 4 physical cores * 2 threads per core = 8 logical cores
print()
print(f"Number of CPU logical cores (threads) in the system: {mp.cpu_count()}") 


#---------------------------------------------------------------------------------------------------
"""Input: 
   1: refGenomeFile: refernce genome file, type: fasta
   2: refGenomeDirectory: refernce genome directory path, type: string
   
   Functionality: open and read reference genome file
   
   Output: 
   1: single string, type: string"""
def readGenomeFile(refGenomeFileName):
    # refGenomeDirectory = '/home/gridsan/raeuf/Raeuf_Notebook/Code/RefSeq/'+refGenomeFileName
    # print(refGenomeFileName)
    refGenomeDirectory='/home/gridsan/alexsli/Alexander_Notebook/RefSeq/'+refGenomeFileName
    genome= ''
    # print(refGenomeDirectory)
    with open(refGenomeDirectory, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
#---------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------
"""Input: 
   1: s: genome fragment, type: string.
   
   Functionality: 
   Generate 10% error per fragment length. Then create a list of randome numbers with size errorNum 
   with those numbers are between 0 and length(s). Then loop the string characters and replace them 
   using their index.
   
   Output: 
   1: snew: same genome fragment with errors, type: string."""
def errorFunc(s, errorRate):

    errorNum = round(len(s)*(errorRate/100))
    
    alphabetNoA = ['T', 'C', 'G']
    alphabetNoT = ['A', 'C', 'G']
    alphabetNoC = ['A', 'G', 'T']
    alphabetNoG = ['A', 'C', 'T']

    s_list= list (s)
    # generate unique (without replacement) random numbers with sample size equal to errorNum, 
    # which is 10% of the length of the input fragment. The generated random numbers are within the 
    # range 0 and the length of the input fragment.
    cps = random.sample(range(0, len(s)), errorNum)
    for i  in cps:
        randomACGT= random.randint(0, 3)
        if randomACGT == 3:
            s_list[i] = ''
        else:
            if s_list[i] == 'A':
                s_list[i] = alphabetNoA[randomACGT]
            elif s_list[i] == 'T':
                s_list[i] = alphabetNoT[randomACGT]
            elif s_list[i] == 'C':
                s_list[i] = alphabetNoC[randomACGT]
            else:
                s_list[i] = alphabetNoG[randomACGT]
    snew= "".join(s_list)
    return snew
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
""" Input: genome, readLenMin, readLenMax, and errorStatus
    genome: refernce genome string, type: string
    readLenMin: smallest read length size, type: integer
    readLenMax: largest read length size, type: integer
    errorStatus: fragment with or without error, string
    Functionality: generate reads with or without errors
    Output: list of reads type: list """
def readsFragmentationGenerator(genome, readLenMin, readLenMax, errorRate, errorStatus):
    global reads
    readLen = random.randint(readLenMin, readLenMax)
    startPos = random.randint(0,len(genome)-readLen)
    read1, read2, read3 = (genome[:startPos],
                           genome[startPos:startPos+readLen],
                           genome[startPos+readLen:])
    if len(read2) !=0 and len(read2) >=readLenMin:
    # if len(read2) !=0:

        if errorStatus == 'YES':
            read2WithError = errorFunc(read2, errorRate)
            reads.append(read2WithError)
        else:
            reads.append(read2)
    for read in [read1, read3]:
        if len(read) !=0 and len(read) >=readLenMin:
        # if len(read) !=0:

            if len(read) > readLenMax:
                readsFragmentationGenerator(read, readLenMin, readLenMax, errorRate, errorStatus)

            else:
                if errorStatus == 'YES':
                    readWithError = errorFunc(read, errorRate)
                    reads.append(readWithError)
                else:
                    reads.append(read)
    return reads
#---------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------
def readsPickerForThreading(threads):
    try:
        if threads == 1:
            randomIndexes= random.randrange(0, len(reads))
            result= [reads[randomIndexes]]
            del reads[randomIndexes]
            return result

        if threads > 1:
            if len(reads) >= threads:
                randomIndexes= random.sample(range(len(reads)), threads)
                result= list(itemgetter(*randomIndexes)(reads))
                for index in sorted(randomIndexes, reverse=True):
                    del reads[index]
                return result

            if len(reads) < threads and len(reads) > 1:
                randomIndexes= random.sample(range(len(reads)), len(reads))
                result= list(itemgetter(*randomIndexes)(reads))
                for index in sorted(randomIndexes, reverse=True):
                    del reads[index]
                return result

            if len(reads) == 1:
                result = [reads[0]]
                del reads[0]
                return result

    except:
        pass
#---------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------
"""Input: 
   1: fragments: genome fragments, type: list
   
   Functionality: 
   Loop through a list of fragments of DNA, create files, and write each fragment into a different 
   file, one fragment at a time.
   
   Output: 
   1: files correspand to the number of threads, type: .fa type"""
def generateFragmentsFiles(fragments, fragmentsCountLengthFile):
    global COUNT
    counter=0
    try:
        for fragment in fragments:
            COUNT+=1
            counter+=1
            readFileName= 'fragmentsFile_'+str(counter)+'_Threading.fa'
            readFile = open(readFileName,"w") 
            readFile.write(f">read{str(COUNT)}\n") 
            readFile.write(fragment+"\n\n") 
            readFile.close()
            readCountLen = [[COUNT,len(fragment)]]
            my_df = pd.DataFrame(readCountLen)
            my_df.to_csv(fragmentsCountLengthFile, mode='a', index=False, header=False)
    except:
        pass
#---------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    
    # start_time= time.time() # timer for entire experiment
    threadsNumber= int(sys.argv[1]) # number of threads

    hGenomeFileName = sys.argv[2] # host genome fasta file name
    hGenomeCount = int(sys.argv[3]) # host genome cell count

    vGenomeFileName = sys.argv[4] # virus genome fasta file name
    vGenomeCount = int(sys.argv[5]) # virus genome count

    minReadLen= int(sys.argv[6])
    maxReadLen= int(sys.argv[7])

    readError= sys.argv[8].upper() # yes or no
    percentError= int(sys.argv[9]) # 10 and higher
    
    vGenomeSymbol= sys.argv[10] # PCV1, MVM, FeLV
    vhRatio= sys.argv[11]
    simulationNumber =sys.argv[12] # simulation number

    vTiter = '9.1x10^9 copies per mL'
    sampleVolm= '400 microL'
    ONTVirusInputSample= '0.01'

    hReadGenome = readGenomeFile(hGenomeFileName).upper() # read host refrence genome file
    hGenomeSize = len(hReadGenome)

    vReadGenome = readGenomeFile(vGenomeFileName).upper() # read virus refrence genome file
    vGenomeSize = len(vReadGenome)


    if maxReadLen > vGenomeSize:
        print(f"Can not enter desired read length (max) > {vGenomeSize} bases. Please try again!")
    else:

        hGenomes = hGenomeCount * [hReadGenome]        
        [readsFragmentationGenerator(i, minReadLen, maxReadLen, percentError, readError) for i in hGenomes]

        vGenomes = vGenomeCount * [vReadGenome]        
        [readsFragmentationGenerator(j, minReadLen, maxReadLen, percentError, readError) for j in vGenomes]



        start_time= time.time() # timer for sequencing

        resultDirectory= f"{vGenomeSymbol}_{str(threadsNumber)}Threads_{readError}{percentError}%ReadErrorRate_{vhRatio}VirusHostRatio_Simulation{simulationNumber}" # result directory
        os.mkdir(resultDirectory) # create directory
        os.chdir(resultDirectory) # set directory for result

        finalResultFileName= f"{vGenomeSymbol}_FragmentsCountLength_{str(threadsNumber)}Threads_{vhRatio}VirusHostRatio_Simulation{simulationNumber}.csv"
        finalResultFile_csv= finalResultFileName.split('.csv')[0]+'_RVDB_result.csv'

        fields_csv= ['alignment_num','sbjct_title', 'sbjct_id', 'sbjct_accession', 'sbjct_length', 'query_id', 'query_length', 'score', 'evalue', 'identities', 'positives', 'gaps', 'align_length', 'strand', 'query_start', 'query_end', 'query', 'match', 'sbjct', 'sbjct_start', 'sbjct_end', 'detect_time_min']
        with open(finalResultFile_csv, 'a') as csvfile:  
            csvwriter = csv.writer(csvfile)  
            csvwriter.writerow(fields_csv)  


        finalResultFile_txt= open(finalResultFileName.split('.csv')[0]+'_RVDB_result.txt', 'w')
        print('\n******************')
        print(f'Date: {today}')
        print('******************')
        finalResultFile_txt.write('******************')
        finalResultFile_txt.write(f'\nDate: {today}')
        finalResultFile_txt.write('\n******************')

        print(f'Virus titer: {vTiter}')
        print(f'Sample volume (microL): {sampleVolm}')
        print(f'ONT virus input sample (microgram): {ONTVirusInputSample}')

        finalResultFile_txt.write(f'\nVirus titer: {vTiter}')
        finalResultFile_txt.write(f'\nSample volume (microL): {sampleVolm}')
        finalResultFile_txt.write(f'\nONT virus input sample (microgram):{ONTVirusInputSample}')
        finalResultFile_txt.write(f'\nVirus:{ONTVirusInputSample}')


        vGenomeFileName= vGenomeFileName.split(".fasta")[0]
        print(f'Virus name: {vGenomeFileName}')
        print(f'Virus genome length (bp): {vGenomeSize}')
        finalResultFile_txt.write(f'\nVirus name: {vGenomeFileName}\n')
        finalResultFile_txt.write(f'Virus genome length (bp): {vGenomeSize}\n')
        finalResultFile_txt.write(f'Virus count: {vGenomeCount}\n')


        hGenomeFileName= hGenomeFileName.split(".fna")[0]
        print(f'Host name: {hGenomeFileName}')
        print(f'Host genome length (bp): {hGenomeSize}')
        finalResultFile_txt.write(f'Host name: {hGenomeFileName}\n')
        finalResultFile_txt.write(f'Host genome length (bp): {hGenomeSize}\n')
        finalResultFile_txt.write(f'Host count: {hGenomeCount}\n')


        print(f'\nVirus to host ratio: {vhRatio}')
        finalResultFile_txt.write(f'\nVirus to host ratio: {vhRatio}')

        print(f'\nNumber of total fragments: {len(reads)}')
        print(f'{percentError}% Error per read: {readError}\n\n')
        finalResultFile_txt.write(f'\nNumber of total fragments: {len(reads)}')
        finalResultFile_txt.write(f'\n{percentError}% Error per read: {readError}\n\n')

        for _ in range(len(reads)):
            COUNT+=0
            if reads:
                pickedReads= readsPickerForThreading(threadsNumber)

                generateFragmentsFiles(pickedReads, finalResultFileName)
                readFilesName= [fname for fname in os.listdir('.') if '_Threading.fa' in fname]
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    blastn_RVDB_results= executor.map(BLASTn_RVDB_DB_search, readFilesName)
                    for result in blastn_RVDB_results:
                        if result:
                            alignment_num= 0
                            for RVDB_result in result:
                                alignment_num+=1

                                BLASTn_RVDB_DB_search_result_parser(finalResultFile_txt,
                                                                    RVDB_result, 
                                                                    alignment_num)
                                with open(finalResultFile_csv, 'a') as csvfile:  
                                    csvwriter = csv.writer(csvfile)  
                                    detection_time = (time.time() - start_time)/60
                                    RVDB_result_alignment_num_detection_time= [alignment_num]+RVDB_result+[detection_time]

                                    csvwriter.writerow(RVDB_result_alignment_num_detection_time)  
                            if ('minute' in RVDB_result[0].lower() or 'mvm' in RVDB_result[0].lower()): ###Change this condition to match the virus you want to detect. The simulation will stop when you find an alignment to a virus with minute or mvm in the description. 
                                break
                    else:
                        time.sleep(0.2) ###Not sure why Raeuf included this wait time. Alex presumes 2 rationales: 1) realistic simulation, 2) limitations by MIT supercloud. If 1), then Alex suggests removing it. 
                        continue
                    break
        end_time= time.time() - start_time
        print(f"Processing took {end_time/60} min using {str(threadsNumber)} parallel threads.")
        finalResultFile_txt.write(f"Processing took {end_time/60} min using {str(threadsNumber)} parallel threads.")
        finalResultFile_txt.close()
#---------------------------------------------------------------------------------------------------





# Must activate the bashrch file first before running this code: source ~/.bashrch 

# Example:
#--------
# python3 -B virus_read_detection_main.py 100 'CHO_GCF_000223135.1_CriGri_1.0_genomic.fna' 5 'Minute virus of mice.fasta' 1538 50 4000 yes 10 MVM 10power2 3


# Breakdown:
#----------
# threadsNumber= int(sys.argv[1]) 100 # number of threads

# hGenomeFileName = sys.argv[2] # 'CHO_GCF_000223135.1_CriGri_1.0_genomic.fna'
# hGenomeCount = int(sys.argv[3]) # 6

# vGenomeFileName = sys.argv[4] # 'Porcine circovirus 1.fasta'
# vGenomeCount = int(sys.argv[5]) # 2253

# minReadLen= int(sys.argv[6]) 50
# maxReadLen= int(sys.argv[7]) 1700

# readError= sys.argv[8].upper() # yes
# percentError= int(sys.argv[9]) # 10

# vGenomeSymbol= sys.argv[10] # PCV1
# simulationNumber =sys.argv[11] # 1



# Blast result example:
#---------------------
"""****Alignment:2****
sbjct_title: acc|GENBANK|V01172.1|Gardner-Arnstein Feline Leukemia oncovirus B envelope gene adjacent 3' long terminal repeat.|Gardner-Arnstein feline leukemia oncovirus B|VRL|18-APR-2005
sbjct_id: gnl|BL_ORD_ID|2469051
sbjct_accession: 2469051
sbjct_length: 2581
query_id: read52

query_length: 55
This is the length of the fragment that we blasted, here: 
'TTACAGACCCCCTACCCCAAAATTTAGCCAGCTAGTGCAGTGGAGCCTTTTCCCA'. This Fragment has a length of 55.

score: 39.0
evalue: 1.41168e-11

identities: 45
This is the number of bases that fully matched, 45 based from the 55 bases exact matched.

positives: 45

gaps: 0
This mean no deletion.

align_length: 48
This is the length of the fragment that was aligined, however not fully matched, only 45 
(i.e see 'positives' from above ) from the 48 were exact match. 

strand: ('Plus', 'Plus')

query_start: 5
This mean the alignment started at position 5 (base 1 indexing) of our fragment 

query_end: 52
This mean the alignment ended at position 52 (base 1 indexing) of our fragment 

query: AGACCCCCTACCCCAAAATTTAGCCAGCTAGTGCAGTGGAGCCTTTTC
match: |||||||||||||||||||||||||||||| |||||||| ||| ||||
sbjct: AGACCCCCTACCCCAAAATTTAGCCAGCTATTGCAGTGGTGCCATTTC

sbjct_start: 2104
This mean the alignment started at position 2104 (base 1 indexing) of the fragment in RVDB database.

sbjct_end: 2151
This mean the alignment started at position 2151 (base 1 indexing) of the fragment in RVDB database.
*****************"""



