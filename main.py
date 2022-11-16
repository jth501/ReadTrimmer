#!usr/bin/env python3

import argparse
import sys
import re
import gzip #detect if files are gunzip - uncomplete

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


def lefttrim(read, leading, quality_scores, quality_line ):
    """ remove X nucleotides from 5 prime end
        return trimmed 5 prime
    """
    # Parameters definition
    window_size = 4
    quality_threshold = 20 
    
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
    


def righttrim(read, trailing, quality_scores, quality_line):
    """ remove X nucleotides from 3 prime end
        return trimmed 3 prime
    """
    # Parameters definition
    window_size = 4
    quality_threshold = 20 
    
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
    return trimmed_seq,quality_scores


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
    phred_dict = dict_creation('coding_keys.txt')
    file = "testfile.txt"
    try:
        #with open(start.file1, "r") as read1, open(start.file2, "r") as r2, open("output.txt", "w") as outfile:
        fastq = open(file, "r")
        read = []
        for line in fastq:
            line = line.strip()
            read.append(line)
            if len(read) == 4:
                quality_coversion = translation_scores(read[3], phred_dict)
                completeleft = lefttrim(read[1], 8,quality_coversion, read[3] ) # third parameter will be from arg parse
                completetrim = righttrim(completeleft[0], 8, completeleft[1],completeleft[2] )
                if len(completetrim[0]) != len(completetrim[1]):
                    sys.exit("Error - sequence length doesn't equal quality length")
              
                if checklen(completetrim[0], minlen=50) is not None:
                    pass
                if meanquality(completetrim[0],completetrim[1],qualitythreshold=40):
                    print(completetrim[0])
                read = []
                break
    except FileNotFoundError:
        print("file not found")  
        
  
if __name__ == "__main__":
    run()
else:
    print("There is an error")
    sys.exit(1)
