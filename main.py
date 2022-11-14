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
    
def lefttrim(read, leading, qualityscore):
    """ remove X nucleotides from 3 prime end
        return trimmed 3 prime
    """
    ltrim = read[leading:]
    for nucleotide in range(len(ltrim)):
        #print(nucleotide)
        for number in qualityscore:
            #print(number)
            if number > 40:
                break
            elif number < 40:
                ltrim_highqual = ltrim.replace(ltrim[nucleotide], " ")
                #print(ltrim_highqual)
        #break

    return ltrim_highqual

def laiasfunction(read):
    #returns list of scores - converted from phred scores
    randscore = []
    for i in read:
        randscore.append(random.randint(1,45))
    return randscore

def righttrim(read, trailing):
    rtrim = read[trailing:]
    for nucleotide in range(len(rtrim)):
        print(nucleotide)
        for number in qualityscore:
            #print(number)
            if number < 40:
                rtrim_highqual = rtrim.replace(rtrim[nucleotide], " ")
                #print(ltrim_highqual)
            if number > 40:
                break
        break
    
    return rtrim_highqual

def run():
    file = "testfile.txt"
    try:
        # with open(start.file1, "r") as read1, open(start.file2, "r") as r2, open("output.txt", "w") as outfile:
        fastq = open(file, "r")
        read = []
        for line in fastq:
            line = line.strip()
            read.append(line)
            if len(read) == 4:
                quality_coversion = laiasfunction(read[3])
                completeleft = lefttrim(read[1], 8,quality_coversion )
                print(completeleft)
                read = []
    except FileNotFoundError:
        print("file not found")  
        
  
if __name__ == "__main__":
    run()
else:
    print("There is an error")
    sys.exit(1)