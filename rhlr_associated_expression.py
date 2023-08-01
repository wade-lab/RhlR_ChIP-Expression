#count the total expression of transcripts associated with each site
#this ended up being super simple

#expression_table = open("C:/RhlR_project/rna_transcripts_operons_allvsall.txt", "r")
expression_table = open(input("Enter Path and Name of Site Associated Gene List: "), "r")
#real_list = open("C:/RhlR_project/Real_Sites.txt", "r")
real_list = open(input("Enter Path and Name of Site List: "), "r")
sites=[]
for a in real_list: 
    sites.append(str(int(a))) #making a list form of all the sites I consider  "real", converting to int first cleans it up a bit
real_transcripts = []
for b in expression_table: 
    tab_list = b.split("\t") #create a list from each line of the input file
    if tab_list[0] == "Site Number": #want to preserve the header line
        real_transcripts.append(b)
        continue
    if tab_list[0] in sites: #then check if the line containing the site positions for each gene is one of the "real" sites
        real_transcripts.append(b) #and save those
output = open(input("Enter Name of Output File: "), "w")
for c in real_transcripts: #then just write the new table of lines to a new file
    output.write(c)

output.close()
real_list.close()
expression_table.close()
