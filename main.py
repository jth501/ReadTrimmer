#!usr/bin/env python3

import argparse
import sys
import re
import random

#create parser
def parserfunc():
    my_parser = argparse.ArgumentParser(prog="ReadTrimmer", description="This is a read trimmer")

    #add the arguments
    my_parser.add_argument("-file1","-f1", type = str,help="the file needed", required=True)
    my_parser.add_argument("-file2", "-f2", help ="the second file needed", required=True)
    my_parser.add_argument("-slidingwindow", nargs="?", default="4", type = int, help="size of the sliding window")
    my_parser.add_argument("-startcut", nargs="?",  default="8", type = int, help ="number of leading nucleotides to remove")
    my_parser.add_argument("-endcut", nargs="?", default="8", type = int, help="number of trailing nucleotides to remove")
    my_parser.add_argument("-minlength", nargs="?", default="100",type = int ,help ="minimum length of the read required")

    #execute the parser method
    args = my_parser.parse_args()
   
    return  args
    

def dict_creation(dict_file):
    """CODE TO CREATE A DICTIONARY FROM THE INPUT FILE"""
    # Open the dictionary text file
    import sys
    try:
        infile = open(dict_file, 'r')
    except IOError as e:
        print("Can't open file, Reason: " + str(e))
        sys.exit(1)

    # Create a dict by reading the imported text file
    # # Key : scores
    # # Values : coding keys (list of two elements)

    phred_dict = dict()
    import re
    headlines = None
    for line in infile:
        headline = re.search(r'^\w+', line)
        if headline is None: 
            line = line.split()
            phred_dict[line[0]] = list(line[2:])
    return phred_dict

def translation_scores(quality_line):
    """Function that translates Ascii characters into scores"""
    quality_scores = list()
    #To translate characters into scores
    for char in quality_line:
        for key,val in phred_dict.items():
            if val[0] == char:                      #Previously the function need to know WHICH PHRED DICT needs to use
                quality_scores.append(int(key))       
    return quality_scores


def lefttrim(read, leading, quality_scores, sequence, quality_lines ):
    """ remove X nucleotides from 3 prime end
        return trimmed 3 prime
    """
    # Parameters definition
    window_size = 4
    quality_thereshold = 20 
    ltrim = read[leading:]
    
    # To calculate the average of a windows from 5'       
    quality_average = 0
    for value in range(len(quality_scores)-1,window_size-2,-1):    
        for i in range(window_size):         
            quality_average += quality_scores[value-i] 
        result = quality_average /window_size      
        quality_average = 0
   
        # To trim window until finding the first one with good quality
        if result >= quality_thereshold:
            quality_line = quality_line[:value+1]
            quality_scores = quality_scores[:value+1]      
            trimmed_seq = trimmed_seq[:value+1]
            break
   

    return quality_scores, quality_line, trimmed_seq


def righttrim(read, trailing, quality_scores,  sequence, quality_lines ):
    # Parameters definition
    window_size = 4
    quality_thereshold = 20 
    
    ltrim = read[leading:]
    rtrim = read[trailing:]
    
      # To calculate the average of a window from 3'
    quality_average = 0
    for value in range(len(quality_scores)-(window_size-1)):  
        for i in range(window_size):  
            quality_average += quality_scores[value+i]      
        result = quality_average /window_size
        quality_average = 0
        
        # To trim window until finding the FIRST ONE with good quality
        if result >= quality_thereshold:
            quality_line = quality_line[value:]
            quality_scores = quality_scores[value:]    
            trimmed_seq = sequence[value:]
            break
    return quality_scores, quality_line, trimmed_seq


def checklen(trimmed, minlen):
    if len(trimmed) > minlen:
        return trimmed
        
def meanquality(lenchecked, qualityscore, qualitythreshold):
    """check the meanquality of the trimmed read"""
    meanquality = (sumquality(qualityscore.count()))
    for i in qualityscore:
        sumquality += i
    if meanquality > qualitythreshold:
        return lenchecked
    
    
def run():
    """Location to run all the functions, open the file from parser
    """
    phred_dict = dict_creation('coding_keys.txt')
    file = "testfile.txt"
    try:
        # with open(start.file1, "r") as read1, open(start.file2, "r") as r2, open("output.txt", "w") as outfile:
        fastq = open(file, "r")
        read = []
        for line in fastq:
            line = line.strip()
            read.append(line)
            if len(read) == 4:
                quality_coversion = translation_scores(read[3], phred_dict)
                completeleft = lefttrim(read[1], 8,quality_coversion ) # third parameter will be from arg parse
                completetrim = righttrim(completeleft, 8, quality_coversion)
                print(completeleft)
                print(completetrim)
                quality_window = slidingWindow_funct(quality_conversion, completetrim, read[3]) 
                if checklen(completetrim, minlen=50) is not None:
                    pass
                if meanquality(completetrim,quality_coversion,qualitythreshold=40):
                    print(completetrim)
                read = []
    except FileNotFoundError:
        print("file not found")  
        
  
if __name__ == "__main__":
    run()
else:
    print("There is an error")
    sys.exit(1)
