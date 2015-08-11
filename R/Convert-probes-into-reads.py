# Changes directory to the folder where the probes' information is stored
import os
os.chdir("C:/Users/krue/Documents/Kevin-Logs/20130320_MDM/bovine-genome-indexed/Bovine.probe_tab")
#print(os.getcwd())

# Processes the file

fr = open("./Bovine.probe_tab") # opens the probe information file stream
fo = open("./Bovine.probe_reads.fastq", "w")
line = fr.readline() # reads and skips the first (=header) line
line = fr.readline() # reads the first probe line
while (line != ""):
    # reading the probe information
    line_split = line.split("\t") # split the line by tab character
    probe_ID = "%s_%s_%s_%s" % (line_split[0], line_split[1], line_split[2], line_split[3]) # concatenates the first 4 column
    probe_seq = line_split[4] # stores the sequence of the probe 
    probe_target_strand = line_split[5].strip() # stored the strandedness of the probe target # NOTE: need to strip the EOL character
    # writing the pseudo-FASTQ read
    print >>fo, "@%s:%s" % (probe_ID, probe_target_strand) # Read header line
    print >>fo, probe_seq # Read sequence line
    print >>fo, "+" # Typical + separator line
    print >>fo, "".join("I" for i in range(len(probe_seq))) # Creates a string of "I"s of same length as probe sequence
    line = fr.readline() # reads the next probe line    
fr.close() # closes the files stream
fo.close()
