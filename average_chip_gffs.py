# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 15:00:34 2023

@author: Nick
"""

#file1 = open(input("Enter Path of First .gff file (No Quotes!): "), "r")
#file2 = open(input("Enter Path of Second .gff file (No Quotes!): "), "r")
file1 = open("C:/RhlR_project/rhlr_gffs/New folder/561_rep1_targets.gff","r")
file2 = open("C:/RhlR_project/rhlr_gffs/New folder/561_rep2_targets.gff","r")

avg_genome = []
for line1 in file1:
    line1 = line1.strip("\n") #just a precaution against unwanted line-breaks sneaking into the function
    gfftab1 = line1.split("\t") #need to have each column as its own string so that we can splice from lines of different length
    coordinate1 = int(gfftab1[3]) #splice the coordinate column and make it an integer so it's useful
    file2.seek(0)
    for line2 in file2:
        line2 = line2.strip("\n") #just a precaution against unwanted line-breaks sneaking into the function
        gfftab2 = line2.split("\t") #need to have each column as its own string so that we can splice from lines of different length
        coordinate2 = int(gfftab2[3]) #splice the coordinate column and make it an integer so it's useful
        if coordinate1 == coordinate2:
            if int((gfftab1[5])) < 0:
                if int((gfftab2[5])) <0:
                    avg_genome.append(gfftab1[0]+"\t"+gfftab1[1]+"\t"+gfftab1[2]+"\t"+gfftab1[3]+"\t"+gfftab1[4]+"\t"+str((int((gfftab1[5]))+int((gfftab2[5])))/2)+"\t"+gfftab1[6]+"\t"+gfftab1[7]+"\t"+gfftab1[8]+"\n")
            if int((gfftab1[5])) > 0:
                if int((gfftab2[5])) > 0:
                    avg_genome.append(gfftab1[0]+"\t"+gfftab1[1]+"\t"+gfftab1[2]+"\t"+gfftab1[3]+"\t"+gfftab1[4]+"\t"+str((int((gfftab1[5]))+int((gfftab2[5])))/2)+"\t"+gfftab1[6]+"\t"+gfftab1[7]+"\t"+gfftab1[8]+"\n")


output = open("C:/RhlR_project/rhlr_gffs/New folder/561_combined.txt","w")
for read in avg_genome:
    if read != 0:
        output.write(read)

file1.close()
file2.close()
output.close()