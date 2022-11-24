 

#PSEUDOCODE OF THE MAIN FUNCTION ‘Run’ 
Run function (): 
    Variables definition #O(1)  
    Creation of the log file #O(1) 
    Try:
    Decompress function (infile/s) #O(n) where n is the number of infiles 
    Check the number of files provided #O(n) where n is the number of infiles
        Make a list to iterate through them #O(1) 
    Counter and list definitions for statistical calculations #O(m) where m is the number of reads 
    Creation of two output files #O(1)  
    Writing in log file column names #O(1)  
    Map reads together if paired – else read single end reads #O(n*m) where n is number of infiles and m is the number of reads 
        Read definition  #O(m) where m is the number of reads 
        Guess_encoding function (read[0]) #O(m) where n is the number of reads 
        If read length == 2: #SINGLE END READS   
            Translation_scores function(read[3], phred_dict, dictionary) #O(m) where m   is the number of reads 
            Lefttrim function(read[1],quality_conversion[0],read[3],start.startcut,   start.slidingwindow, start.qualitythreshold) 
            #O(m*l*k) where m is the number of reads, l is the length of the sliding window and k is the number of bases  
            Righttrim function (completeleft[0], completeleft[1], completeleft[2], start.endcut, start.slidingwindow, start.qualitythreshold) 
            #O(m*l*k) where m is the number of reads, l is the length of the sliding window and k is the number of bases 

            #Statistical calculations  
            Average quality of the read #O(m*k) where m is the number of reads and k the number of bases 
            Bases counters increases #O(k) where k is the number of bases 
            Length of the entry #O(k) where k is the number of bases 
            Average length of all entries #O((m*k)/m) where m is the number of reads and k is the number of bases 
            Average quality of all entries #O((m*k)/m) where m is the number of reads and k is the number of bases 

        If read length == 4: #PAIRED END READS 
        #(The process to trim paired end reads is the same explained above but with two entries simultaneously) 
        Meanquality function (trimmedqual,qualitythreshold=30) #O(m) where m is the number of reads 
        Checklen function (completetrim, minlen=50) #O(m) where m is the number of reads 
        Print statistical results in the log file 
    Except FileNotFOundError: 
        Print error in log file 