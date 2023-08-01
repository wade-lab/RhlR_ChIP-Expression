# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 11:18:46 2023

@author: Nick
"""
# This is just a simple script to normalize the .gff ChIP data I visualize in a genome browser
# Just run the script and enter the inputs that are prompted
# This script combines both plus and minus reads and normalizes both by their combined sum
# When you copy the file name as a path in windows 11 it will automatically include quotes, but those will be read
# into the character string from the input function and not be a functional path, so make sure to only input the raw path
gff = open(input("Enter Path of .gff file (No Quotes!): "), "r")
enrichment = []
for line in gff:
    line_list = line.split("\t")
    if len(line_list) < 9: #There's often a single blank line at the end of gffs which gives an out of range error when I try to slice from a "list" made from it
        continue
    enrichment.append(int(float(line_list[5]))) #need to float values first since int() doesn't like decimals
pos_enrichment = []
for value in enrichment:
    pos_enrichment.append(abs((value)))#I'm summing positive and negative numbers so need to use absolute values to account for them all
reads_sum = sum(pos_enrichment)
normalized_reads = []
for number in enrichment:
    normalized_reads.append(100000000*(number/reads_sum)) #normalizing by reads/total reads and multiplying that by 100 million
output = open(input("Enter Path/Name of Normalized File: "),"w")
line_number = 0
gff.seek(0) #Need to make sure the gff is being read from the start for the next for loop
new_label = input("Enter New Track Label for File (No Quotes!): ")
for line in gff:
    line_list = line.split("\t")
    if len(line_list) <9: #skipping lines that aren't complete again
        continue
    output.write(line_list[0]+"\t"+line_list[1]+"\t"+new_label+"\t"+line_list[3]+"\t"+line_list[4]+"\t"+str(normalized_reads[line_number])+"\t"+line_list[6]+"\t"+line_list[7]+"\t"+line_list[8])
    line_number = (line_number +1) #since I calculate values line-by line I can just pull one line at a time
    #for loop will start at 0, and go 1-by-1 along the range of lines, I start the line number as 0, then add 1 each iteration
    #so that the line and new normalized read line up
gff.close()
output.close()

#make sure to compare the old gff with your new normalized gff, with automatice scaling in signalmap the two should
#have the same relative sizes for their peaks, just with different numbers