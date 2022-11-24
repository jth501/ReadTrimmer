#!user/bin/env python3
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
    #parameters here are currenty tests but can be rewritten to other file for testing
    args = my_parser.parse_args("-f testfile.txt testfile2.txt".split())
   
    return args
    
def decompress(textwrapper):
    """ check if files are fasta files - returns false and exits if not
        if true files are checked for compressiion status and decompressed if necessary
        returns 
    """
    fqcheck = (".fq", ".fastq", ".fq.gz", ".fastq.gz",".txt") #remove testcase txt
    for i in textwrapper:
        if i.endswith(fqcheck):
            print("file is in fastq format")
        else:
             return False
        checkcompress = re.search(r"\w+.g[un]z[ip]",i)
        # checkcompress = re.search(r"\w+.txt",i) ##test case
        if checkcompress is None:
            print("not a zipped file")
        else:
            return True

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
    """Function to translates Ascii characters into scores
        returns list of scores
    """
    #make score a list if single input
    if len(quality_line) != 2:
        quality_line = [quality_line]
    quality_scores = list()
    #To translate characters into scores for read/ reads
    for qualscore in quality_line:
        score = []
        for char in qualscore:
            for key,val in phred_dict.items():
                if val[encoding_dict] == char:                   
                    score.append(int(key))
        quality_scores.append(score)  
        score = []
    return quality_scores

def lefttrim(read, quality_scores, quality_line, leading, window_size, quality_threshold):
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
    
def righttrim(read,quality_scores, quality_line,trailing, window_size,quality_threshold):
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
    for seq in trimmed:
        if len(seq) > minlen:
            return True
        else: 
            return False # write index to log file
        
def meanquality(qualityscore, qualitythreshold):
    """ calculates average qulaity of all the scores from a read
        returns read if average is greater than quality threshold level
    """
    end = False
    for score in qualityscore:
        sumquality = 0
        for i in score:
            sumquality += int(i)
        meanquality = int(sumquality/(len(qualityscore)))
        if meanquality > qualitythreshold:
            end = True
    return end
    
def pairedend(read, leading, trailing, qualityscores, windowsize, qualthresh):
    #remove adapters
    forward = read[0][leading:(-leading-1)]
    reverse = read[1][trailing:(-trailing-1)]
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
    
    #trim reverse read from the 5 prime side
    for value in range(0,len(qualityscores[1])):
        if sum(qualityscores[1][-(value+3):])/3 > qualthresh:
            trimmedr = reverse
            break
        else:
            trimmedr = reverse[:-(value)]
    trimmed.append(trimmedr)

    #check if reads are the same length to discard later
    if len(trimmedf) != len(trimmedr):
        trimmed.append(False)
    return trimmed

def run():
    """Location to run all the functions, open the file from parser
    """
    start = parserfunc()
    phred_dict = dict_creation('coding_keys.txt')
    logfile = open("log_file.txt", "w")
    try:
        #check if file is fastq and if compressed change open function
        starter = decompress(start.files)
        if starter is False:
            sys.exit("The file is not in fastq format \nOnly .fq, .fastq, .fq.gz, .fastq.gz format accepted")
        open_opt = gzip.open if starter is True else open
        
        #checks number of files provided and makes list to iterate through them
        if len(start.files) == 2:
            fastq = [open_opt(start.files[0], "r"), open_opt(start.files[1], "r")]
            size, num = (4 ,1)
        elif len(start.files) == 1:
            fastq = [open_opt(start.files[0], "r")]
            size, num = (2 ,0)
        else:
            sys.exit("Only single file/two paired files allowed as input")
      
        # Counter to keep track of trimmed and removed reads
        trimmed_reads, removed_reads = (0,0)
        #Counter to keep track of statistics
        number_entries, number_A, number_C, number_T ,number_G = (0,0,0,0,0)
        length_entries = list()
        conversion_list = list()
        quality_list = list()
      
        #create output file to write to and initialise read   
        outfile = open("outfile.txt", "w")
        outfile2 = open("outfile2.txt", "w")
        read = []
        
        #map reads together if paired - else read single read
        for line in map(list,zip(fastq[0],fastq[num])):
            line = [x.strip() for x in line]
            read.append(line)
            #read length is 2 for single but 4 for paired
            if len(read) == len(start.files)*2:   
                dictionary = guess_encoding(read[0])
                print(read)
                #Average quality for each entry
                (sum_1, sum_2,sum_list)= (0,0,0)
                
                if  len((start.files))*2 == 2:
                    read = ["".join(x) for i in read for x in i] 
                    #convert reads and trim
                    quality_conversion = translation_scores(read[3], phred_dict, dictionary) #returns list of decimal score
                    completeleft = lefttrim(read[1],quality_conversion[0],read[3], 
                                            start.startcut, start.slidingwindow, start.qualitythreshold) # third parameter will be from arg parse
                
                    completetrim = righttrim(completeleft[0], completeleft[1], completeleft[2],
                                            start.endcut, start.slidingwindow, start.qualitythreshold)
                    #statistics
                    for i in quality_conversion:
                        sum_1 += i
                    average_quality_1= sum_1/len(quality_conversion)
                    quality_list.append(average_quality_1)
                    
                    number_entries += 1
                    initial_read = read[1]
                    number_A = initial_read.count('A')
                    number_C = initial_read.count('C')
                    number_T = initial_read.count('T')
                    number_G = initial_read.count('G')
                    
                    #TODO: check if id numbers oare the same in identity line
                    outfile.write("{0}\n{1}\n{2}\n{3}\n".format("".join(read[0]),"".join(completetrim[0]), 
                                                                    "".join(read[2]),"".join(read[3][:len(completetrim[0])])))
                elif len(start.files)*2 == 4:
                    
                    #convert ascii values for each read/single read to decimal value
                    quality_conversion = translation_scores(read[3], phred_dict, dictionary) #returns list of decimal score
                    completetrim = pairedend(read[1], start.startcut, start.endcut,
                                             quality_conversion,start.slidingwindow, start.qualitythreshold)
                    
                    for i in quality_conversion[0]:
                        sum_1 += i
                    average_quality_1= sum_1/len(quality_conversion[0])
                    quality_list.append(average_quality_1)
                    for i in quality_conversion[1]:
                        sum_2 += i
                    average_quality_2= sum_2/len(quality_conversion[1])
                    quality_list.append(average_quality_2)
                    
                    number_entries += 2 
                    
                    #still need to fix ability for functions to read two reads at once
                    #Number of bases in forward entries
                    initial_read_f = read[1][0]
                    number_A = initial_read_f.count('A')
                    number_C = initial_read_f.count('C')
                    number_T = initial_read_f.count('T')
                    number_G = initial_read_f.count('G')

                    #Number of bases in reverse entries
                    initial_read_r = read[1][1]
                    number_A += initial_read_r.count('A')
                    number_C += initial_read_r.count('C')
                    number_T += initial_read_r.count('T')
                    number_G += initial_read_r.count('G')
                    
                    #Length of each entry
                    length_entries.append(len(initial_read_f))
                    length_entries.append(len(initial_read_r))    
                    #take paired reads and trim or single read
                    outfile.write("{0}\n{1}\n{2}\n{3}\n".format("".join(read[0][0]),"".join(completetrim[0]), 
                                                                    "".join(read[2][1]),"".join(read[3][0][:len(completetrim[0])])))
                    outfile2.write("{0}\n{1}\n{2}\n{3}\n".format("".join(read[0][1]),"".join(completetrim[1]),
                                                                     "".join(read[2][1]),"".join(read[3][1][:len(completetrim[0])])))  
                else:
                    sys.exit("there is an error in file processing, try again")
                    
                if len(completetrim[0]) != len(completetrim[1]):
                        sys.exit("Error - sequence length doesn't equal quality length")
                    
                if completetrim == None:
                    print("The reads could not be processed")

                #TODO: need to format the output into two files for two reads
                #TODO: remove adapter from right and left of read - the user input replaces adapter removal as the reads have adapters
                #To check if the read has been trimmed or removed
                if len(completetrim[0][0]) != len(initial_read_f)  :
                    trimmed_reads += 1
                if len(completetrim[0][1]) != len(initial_read_r)  :
                    trimmed_reads += 1
            
                trimmedqual = quality_conversion[:len(completetrim[0])] 
                #TODO: check if this works
                if meanquality(trimmedqual,qualitythreshold=30) is False:
                    print("read removed mean")
                    removed_reads += 1
                elif checklen(completetrim, minlen=50) is False:
                    print("read removed too short ")
                    removed_reads += 1
                read = []
        for i in fastq:
            i.close()
        outfile.close()
        outfile2.close()

        #To calculate statistics results
        #Average length
        sum = 0    
        for i in length_entries:
            sum += i
        average_length = sum/number_entries
        #Average quality of each entry
        for i in quality_list:
            sum_list += i
        total_average = sum_list / len(quality_list)
        original_list = quality_list[:]
        quality_list.sort()
        #Best and worst 10%
        percentage = int(number_entries*0.1)
        best10 = quality_list[-percentage:]
        worst10 = quality_list[:percentage]
        
        # To print the results in the log file
        now = datetime.datetime.now()
        print((str(now)), '\t','All reads from this FILE were CORRECTLY TRIMMED!', file = logfile)
        print('\nNumber of trimmed reads:\t', trimmed_reads, file = logfile )
        print('\nNumber of removed reads:\t', removed_reads, file = logfile)
        print('\nSTATISTICS (file without trimming):', file = logfile )
        print('\nNumber of entries:\t',number_entries, file = logfile )
        print('\nNumber of each bases:\t', file = logfile )
        print('\tNumber A:\t',number_A, file = logfile )
        print('\tNumber C:\t',number_C, file = logfile )
        print('\tNumber T:\t',number_T, file = logfile )
        print('\tNumber G:\t',number_G, file = logfile )
        print('\nLength of each entry:\t', file = logfile )
        for i in range(1,number_entries+1):
            print('\tEntry\t', i, ':\t', length_entries[i-1], file = logfile)
        print('\nAverage length of entries:\t', round(average_length,2), file = logfile )
        print('\nQuality average of entries:\t',round(total_average,2), file = logfile )
        print('\nBest 10% quality entries:\t', file = logfile )
        for j in best10:
            pos = original_list.index(j)
            print('\tEntry:',pos+1, '\tAverage quality:',round(j,2), file=logfile )
        print('\nWorst 10% quality entries:\t', file = logfile )
        for j in worst10:
            pos = original_list.index(j)
            print('\tEntry:',pos+1, '\tAverage quality:',round(j,2), file=logfile)
        
        
    except FileNotFoundError as e:
        now = datetime.datetime.now()
        print((str(now)), '\t','ERROR: File NOT found!!', file = logfile)


if __name__ == "__main__":
    run()
    if run is not None:
        print("Successful run!")
else:
    print("There is an error")
    sys.exit(1)