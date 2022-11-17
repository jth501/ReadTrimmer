#!usr/bin/env python3

import argparse
import sys
import re
import gzip #detect if files are gunzip - uncomplete
import string

#create parser
def parserfunc():
    my_parser = argparse.ArgumentParser(prog="ReadTrimmer", description="This is a read trimmer")

    #add the arguments
    my_parser.add_argument("-files","-f", nargs="+",help="the file needed")
    my_parser.add_argument("-slidingwindow", nargs="?", default="4", type = int, help="size of the sliding window")
    my_parser.add_argument("-startcut", nargs="?",  default="8", type = int, help ="number of leading nucleotides to remove")
    my_parser.add_argument("-endcut", nargs="?", default="8", type = int, help="number of trailing nucleotides to remove")
    my_parser.add_argument("-minlength", nargs="?", default="100",type = int ,help ="minimum length of the read required")
    my_parser.add_argument("-minlength", nargs="?", default="20",type = int ,help ="min quality level of read")

    #execute the parser method
    #parameters here are currenty tests but can be rewritten to otger file for testing
    args = my_parser.parse_args("-f testfile.txt ".split())
   
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

def translation_scores(quality_line, phred_dict):
    """Function to translates Ascii characters into scores"""
    quality_scores = list()
    #To translate characters into scores
    for char in quality_line:
        for key,val in phred_dict.items():
            if val[0] == char:                      #Previously the function need to know WHICH PHRED DICT needs to use [0] or [1]
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
    if len(trimmed) > minlen:
        return trimmed
        
def meanquality(lenchecked, qualityscore, qualitythreshold):
    """check the meanquality of the trimmed read"""
    sumquality = 0
    meanquality = sumquality/(len(qualityscore))
    for i in qualityscore:
        sumquality += i
    if meanquality > qualitythreshold:
        return lenchecked
    
    
def run():
    """Location to run all the functions, open the file from parser
    """
    start = parserfunc()
    phred_dict = dict_creation('coding_keys.txt')
    
    try:
        starter = decompress(start.files)
        print(starter)
        open_opt = gzip.open if starter is not None else open
        #checks number of files provided and makes TextWrapper list to iterate through belpow
        if len(start.files) == 2:
            fastq = [open_opt(start.files[0], "r"), open_opt(start.files[1], "r")]
            size, num = (4 ,1)
        elif len(start.files[0]) == 1:
            fastq = [open_opt(start.files[0], "r")]
            size, num = (2 ,0)
        else:
            sys.exit("Only single file/two paired files allowed as input")
            
        #create output file to write to and initialise read   
        outfile = open("outfile.txt", "w")
        read = []
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
                print(read)
                #still need to fix ability for functions to read two reads at once
                quality_coversion = translation_scores(read[3], phred_dict)
                completeleft = lefttrim(read[1], 8,quality_coversion, read[3], start.slidingwindow, start.startcut ) # third parameter will be from arg parse
                completetrim = righttrim(completeleft[0], 8, completeleft[1],completeleft[2], start.slidingwindow, start.endcut )
                print(completetrim[0])
                if len(completetrim[0]) != len(completetrim[1]):
                    sys.exit("Error - sequence length doesn't equal quality length")
                if checklen(completetrim[0], minlen=50) is not None:
                    pass
                if meanquality(completetrim[0],completetrim[1],qualitythreshold=40) is not None:
                    print(completetrim[0])
                read = []
                break
    except FileNotFoundError as e:
        sys.exit("file not found", e)  
    except IOError as x:
        print(x) #could append this to the log file
        
  
if __name__ == "__main__":
    run()
else:
    print("There is an error")
    sys.exit(1)
