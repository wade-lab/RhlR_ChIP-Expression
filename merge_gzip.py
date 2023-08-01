import os

sequence_dir = input("Enter Path of the File Containing Files to be Merged ONLY (include a \ at the end of the path): ") #make sure to include a \ at the end of the path
sequence_list = os.listdir(sequence_dir)

import gzip #don't feel like unzipping so many large files so how about this

def combine_rnaseq(item):
    sequence_file = gzip.open(sequence_dir + item)
    combined_run = gzip.open(sequence_dir + "Merged_file.fastq.gz", "a") 
    for line in sequence_file:
        combined_run.write(line)
    combined_run.close()

for file in sequence_list:
    if file.endswith(".fastq.gz"):
        combine_rnaseq(file)
        print(file + " added")



print("done")
