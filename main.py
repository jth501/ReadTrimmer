#!usr/bin/env python3

import argparse
import sys
import re
import gzip #detect if files are gunzip - uncomplete
import datetime
import string

#create parser
def parserfunc():
    my_parser = argparse.ArgumentParser(prog="ReadTrimmer", description="This is a read trimmer")

    #add the arguments
    my_parser.add_argument("-files","-f", nargs="+",help="the file needed")
    my_parser.add_argument("-slidingwindow", nargs="?", default="4", type = int, help="size of the sliding window")
    my_parser.add_argument("-startcut", nargs="?",  default="8", type = int, help ="number of leading nucleotides to remove")
    my_parser.add_argument("-endcut", nargs="?", default="8", type = int, help="number of trailing nucleotides to remove")
    my_parser.add_argument("-minlength", nargs="?", default="70",type = int ,help ="minimum length of the read required")
    my_parser.add_argument("-qualitythreshold", nargs="?", default="20",type = int ,help ="min quality level of read")

    #execute the parser method
    #parameters here are currenty tests but can be rewritten to otger file for testing
    args = my_parser.parse_args("-f testfile.txt testfile2.txt".split())
   
    return args
    
def decompress(textwrapper):
    for i in textwrapper:    
        checkcompress = re.search(r"\w+.g[un]z[ip]",i)
       # checkcompress = re.search(r"\w+.txt",i) ##test case
        if checkcompress is not None:
            print("is a zipped file")
            return textwrapper

def dict_creation(dict_file):
    """To create a dictionary from the infile data"""
    # Open the dictionary text file
    try:
        infile = open(dict_file, 'r')
    except IOError as e:
        print("Can't open file, Reason: " + str(e))
        sys.exit(1)

    # Create a dict by reading the imported text file
    # # Key : scores
    # # Values : coding keys (list of two elements)

    phred_dict = dict()
    headlines = None
    for line in infile:
        headline = re.search(r'^\w+', line)
        if headline is None: 
            line = line.split()
            phred_dict[line[0]] = list(line[2:])
    return phred_dict


def guess_encoding(id_line):
    """Function to detect which phred has to be used (+33 or +64)"""
    if id_line[0:3] == '@HWI':
        dictionary = 1
    else:
        dictionary = 0
    return dictionary
    
    
def translation_scores(quality_line, phred_dict, encoding_dict):
    """Function to translates Ascii characters into scores"""
    quality_scores = list()
    #To translate characters into scores
    for char in quality_line:
        for key,val in phred_dict.items():
            if val[encoding_dict] == char:                   
                quality_scores.append(int(key))     
    return quality_scores

def lefttrim(read, leading, quality_scores, quality_line, window_size, quality_threshold):
    """ remove X nucleotides from 5 prime end
        return trimmed 5 prime
    """
    # Remove X nucleotides
    ltrim = read[leading:]
    quality_scores = quality_scores[leading:]
    quality_line = quality_line[leading:]
    
    # To calculate the average of a windows from 5'       
    quality_average = 0
    for value in range(len(quality_scores)-(window_size-1)):  
        for i in range(window_size):  
            quality_average += quality_scores[value+i]      
        result = quality_average /window_size
        quality_average = 0
        
        # To trim window until finding the FIRST ONE with good quality
        if result >= quality_threshold:
            quality_line = quality_line[value:]
            quality_scores = quality_scores[value:]    
            trimmed_seq = ltrim[value:]
            break
    return trimmed_seq,quality_scores, quality_line
    
def righttrim(read, trailing, quality_scores, quality_line,window_size,quality_threshold):
    """ remove X nucleotides from 3 prime end
        return trimmed 3 prime
    """
    # Parameters definition
    
    # Remove X nucleotides
    rtrim = read[:-trailing]
    quality_scores = quality_scores[:-trailing]
    quality_line = quality_line[:-trailing]
    
    # To calculate the average of a window from 3'
    quality_average = 0
    for value in range(len(quality_scores)-1,window_size-2,-1):    
        for i in range(window_size):         
            quality_average += quality_scores[value-i] 
        result = quality_average /window_size      
        quality_average = 0

        # To trim window until finding the first one with good quality
        if result >= quality_threshold:
            quality_line = quality_line[:value+1]
            quality_scores = quality_scores[:value+1]      
            trimmed_seq = rtrim[:value+1]
            break
    return trimmed_seq,quality_scores, quality_line

def checklen(trimmed, minlen):
    """ checks the length of the trimmed read
        returns read of greater than minimum length
    """
    #print(trimmed)
    for seq in trimmed:
        print(seq)
        print(len(seq))
        if len(seq) > minlen:
            print(0)
            return True
        else: 
            return False # write index to log file
        
def meanquality(lenchecked, qualityscore, qualitythreshold):
    """check the meanquality of the trimmed read"""
    sumquality = 0
    for i in qualityscore:
        sumquality += i
    meanquality = int(sumquality/(len(qualityscore)))
    if meanquality > qualitythreshold:
        return True
    else:
        return False
    
def pairedend(read, leading, trailing, qualityscores, windowsize, qualthresh):
    forward = read[0][leading:]
    reverse = read[1][:-trailing]
    #cut from left and right at the same time until trim stops - should be equal
    
    #trim forward read from the 5 prime side
    trimmed = []
    for value in range(0,len(qualityscores[0])):
        if sum(qualityscores[0][value:value+3])/3 > qualthresh:
            trimmedf = forward
            break
        else:
            trimmedf = forward[value:]
    trimmed.append(trimmedf)
    
def run():
    """Location to run all the functions, open the file from parser
    """
    start = parserfunc()
    phred_dict = dict_creation('coding_keys.txt')
    try:
        starter = decompress(start.files)
        open_opt = gzip.open if starter is not None else open
        #checks number of files provided and makes TextWrapper list to iterate through belpow
        if len(start.files) == 2:
            fastq = [open_opt(start.files[0], "r"), open_opt(start.files[1], "r")]
            size, num = (4 ,1)
        elif len(start.files) == 1:
            fastq = [open_opt(start.files[0], "r")]
            size, num = (2 ,0)
        else:
            sys.exit("Only single file/two paired files allowed as input")
            
        #create output file to write to and initialise read   
        outfile = open("outfile.txt", "w")
        logfile = open("log_file.txt", "w")
        read = []

        # Counter to keep track of trimmed and removed reads
        trimmed_reads = 0
        removed_reads = 0
        
        #Counter to keep track of statistics
        number_entries = 0
        number_A = 0
        number_C = 0
        number_T = 0
        number_G = 0
        length_entries = list()


        #check if fastq format in first read and exit if so - unfinished
        for line in map(list,zip(fastq[0],fastq[num])):
            line = [x.strip() for x in line]
            read.append(line)
            #read length is 2 for single but 4 for paired
            if len(read) == size:
                if size == 2:
                    read = ["".join(x) for i in read for x in i]
                if size == 4:
                    pass
                number_entries += 1
                length_entries.append(len(read[1]))

                #still need to fix ability for functions to read two reads at once
                initial_read = read[1]
                number_A = initial_read.count('A')
                number_C = initial_read.count('C')
                number_T = initial_read.count('T')
                number_G = initial_read.count('G')

                dictionary = guess_encoding(read[0])
                #convert ascii values for each read/single read to decimal value
                quality_coversion = translation_scores(read[3], phred_dict, dictionary) #returns list of decimal score
                #take paired read
                if len(read[1]) == 2:
                    completetrim = pairedend(read[1], start.startcut, start.endcut,quality_coversion,start.slidingwindow,start.qualitythreshold)
                elif len(read[1]) == 1:
                    completeleft = lefttrim(read[1],start.startcut,quality_coversion, read[3], start.slidingwindow, start.qualitythreshol ) # third parameter will be from arg parse
                    completetrim = righttrim(completeleft[0], start.endcut, completeleft[1],completeleft[2], start.slidingwindow,start.qualitythreshold )
                    if len(completetrim[0]) != len(completetrim[1]):
                        sys.exit("Error - sequence length doesn't equal quality length")
                else:
                    sys.exit("there is an error in file processing, try again")
            
                if completetrim == None:
                    print("The reads could not be processed")
                #TODO: need to format the output into two files for two reads
                #To check if the read has been trimmed or removed
                if len(completetrim[0]) != len(initial_read):
                    trimmed_reads += 1
                if meanquality(completetrim[0],quality_coversion[0][:len(completetrim[0])],qualitythreshold=30) or checklen(completetrim[0], minlen=50) is False:
                    print("check why checklen is false")
                    removed_reads += 1
                else:
                    outfile.write("{0}\n{1}\n{2}\n{3}\n".format(read[0],completetrim[0], read[2],read[3][:len(completetrim[0])]))
                    #write to second output file if there are two reads
                print(read[1], 'length', len(read[1]))
                quality_coversion = translation_scores(read[3], phred_dict, dictionary)
                print(quality_coversion, 'length', len(quality_coversion))
                completeleft = lefttrim(read[1], 8,quality_coversion, read[3], start.slidingwindow, start.startcut ) # third parameter will be from arg parse
                completetrim = righttrim(completeleft[0], 8, completeleft[1],completeleft[2], start.slidingwindow, start.endcut )
                print(completetrim[0])
                if len(completetrim[0]) != len(completetrim[1]):
                    sys.exit("Error - sequence length doesn't equal quality length")
                if checklen(completetrim[0], minlen=50) is not None:
                    pass
                if meanquality(completetrim[0],completetrim[1],qualitythreshold=40) is not None:
                    print(completetrim[0])

                #To check if the read has been trimmed or removed
                if len(completetrim[0]) != len(initial_read):
                    trimmed_reads += 1
                if meanquality(completetrim[0],completetrim[1],qualitythreshold=40) or checklen(completetrim[0], minlen=50) is None:
                    removed_reads += 1
                read = []

        # To put the results in a log file
        now = datetime.datetime.now()
        print((str(now)), '\t','All reads from this FILE were CORRECTLY TRIMMED!', file = logfile)
        print('Number of trimmed reads:\t', trimmed_reads, file = logfile )
        print('Number of removed reads:\t', removed_reads, file = logfile)
        print('\nSTATISTICS (file without trimming):', file = logfile )

        print('Number of entries:\t',number_entries, file = logfile )
        print('Number of each bases:\t', file = logfile )
        print('\tNumber A:\t',number_A, file = logfile )
        print('\tNumber C:\t',number_C, file = logfile )
        print('\tNumber T:\t',number_T, file = logfile )
        print('\tNumber G:\t',number_G, file = logfile )

        print('Length of each entry:\t', file = logfile )
        for i in range(1,number_entries+1):
            print('\tEntry\t', i, ':\t', length_entries[i-1], file = logfile)



        print('Average length of entries:\t', file = logfile )
        print('Quality average length of entries:\t', file = logfile )
        print('Best 10% quality entries:\t', file = logfile )
        print('Worst 10% quality entries:\t', file = logfile )





    except FileNotFoundError as e:
        print("file not found", e) 
        # To put the results in a log file
        now = datetime.datetime.now()
        print((str(now)), '\t','ERROR: File NOT found!!', file = logfile)

         
        
  
if __name__ == "__main__":
    run()
else:
    print("There is an error")
    sys.exit(1)
