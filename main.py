#!usr/bin/env python3

import argparse
import sys
import os

#create parser
def parserfunc():
    my_parser = argparse.ArgumentParser(prog="ReadTrimmer", description="This is a read trimmer")

    #add the arguments
    my_parser.add_argument("-file1","-f1", type = str,help="the file needed", required=True)
    my_parser.add_argument("-file2", "-f2", help ="the second file needed", required=True)
    my_parser.add_argument("SlidingWindow", nargs="?", default="4", type = int, help="size of the sliding window")
    my_parser.add_argument("startcut", nargs="?" ,type = int, help ="number of leading nucleotides to remove")
    my_parser.add_argument("endcut", nargs="?",type = int, help="number of trailing nucleotides to remove")
    my_parser.add_argument("minlength", nargs="?", default="100",type = int ,help ="minimum length of the read required")

    #execute the parser method
    args = my_parser.parse_args()
   
    return  args


def main():
    start = parserfunc()
    
    with open(start.file1, "r") as read1, open(start.file2, "r") as r2, open("output.txt", "w") as outfile:
        for line in read1:
            if line.startswith("@"):
                print(line)
            break

if __name__ == "__main__":
    main()
else:
    print("There is an error")
    print(__name__)
    sys.exit(1)