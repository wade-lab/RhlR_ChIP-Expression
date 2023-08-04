#In order to run this code will first need to convert the .wig rockhopper output to a .gff (use gff writer from wig)
folder = input("Enter Path to Folder Containing Needed Files: ")
peaks = open(input("Enter Path and Name of File Containing Peak Coordinates: ","r") #make sure the folder contains the list of peak coordinates
peak_list = []
for a in peaks:
    peak_list.append(int(a)) #creating list form of the peak file so I have the coordinated of each peak
#first determining which genes in the file are associated with an operon
operon_table = open(input("Enter Rockhopper Output File Path: "),"r")  #"C:/RhlR_project/NC_008463_operons_allvsall.txt"
operon_binding = []
operon_starts_plus = []
operon_starts_minus = []
for c in operon_table:
    tab_list = c.split("\t")
    if len(tab_list)<=1: #python gives an error in the last couple lines (which only have a \n) of the rockhopper output since I split the line to make a list then try to slice from a nonexistent part of that list
        continue
    if tab_list[1] =="Stop": #want to preserve the first line and also make sure python doesn't try to make it into integers
        operon_binding.append("Site Number"+ "\t" + "Site Location" + "\t" + c)
        continue
    if tab_list[2] == "+":
        for f in peak_list:
            if f in range((int(tab_list[0])-500), int(tab_list[0])):
                operon_binding.append(str(peak_list.index(f)+1) + "\t" + str(f) + "\t"+ c)
                operon_starts_plus.append(f) #need a list of the location of the starting gene of each operon so I can avoid listing them twice later
    if tab_list[2] == "-":
        for f in peak_list:
            if f in range(int(tab_list[1]), (int(tab_list[1])+500)):
                operon_binding.append(str(peak_list.index(f)+1) + "\t" + str(f) + "\t"+ c)
                operon_starts_minus.append(f) #need a list of the location of the starting gene of each operon so I can avoid listing them twice later
operon_output = open(folder+ "rna_operon_table_allvsall2.txt", "w") #now writing the rhlr associated genes to file for later reference
for line in operon_binding:
    operon_output.write(line)

peaks.close()
operon_output.close()
operon_table.close()

transcript = open(folder + "NC_008463_transcripts_allvsall.txt","r")
bound_expressed = []
operons = open(folder+"rna_operon_table_allvsall2.txt", "r")
peaks = open(folder+"RhlR_peak_coordinates.txt","r")
for c in transcript:
    tab_list = c.split("\t") 
    if len(tab_list)<=1: #python gives an error in the last couple lines (which only have a \n) of the rockhopper output since I split the line to make a list then try to slice from a nonexistent part of that list
        continue
    if tab_list[1] =="Translation Start": #want to preserve the first line and also avoid errors hence the continue
        bound_expressed.append("Site Number"+ "\t" + "Site Location" + "\t" + "Internal binding" + "\t" + c)
        continue
    operons.seek(0)
    #need to make operon check strand specific to the genes it takes from the transcript list to prevent pulling of overlapping genes on opposing strand!!!!
    for r in operons: #first want to see if the aligned gene is a part of an operon, so checking each one for each gene
        operon_list = r.split("\t")
        if operon_list[2] =="Start": #will get errors if following lines try to run on strings
            continue
        operon_range = range((int(operon_list[2])),(int(operon_list[3])))#determine the sequence covered by the operon
        if operon_list[4] == "+": #need to know which strand we're using
            if tab_list[1] == '': #many lines don't have a value for transcript start/end to need to search exclusively for thigns that exist
                if int(tab_list[0]) in operon_range:
                    if tab_list[4] == "+": #need to make sure both operon and genome lines are on the same strand in case of overlapping genes on the opposing strand
                        bound_expressed.append(str(operon_list[0]) + "\t" + str(operon_list[1]) + "\t"+ "false" + "\t" + c)
                        #print(tab_list) #just a debugging line
            if tab_list[1] != '':
                if int(tab_list[1]) in operon_range:
                    if tab_list[4] == "+":
                        bound_expressed.append(str(operon_list[0]) + "\t" + str(operon_list[1]) + "\t"+ "false" + "\t" + c)
        if operon_list[4] == "-":
            if tab_list[2] == '':
                if int(tab_list[3]) in operon_range:
                    if tab_list[4] == "-":
                        bound_expressed.append(str(operon_list[0]) + "\t" + str(operon_list[1]) + "\t"+ "false" + "\t" + c)
            if tab_list[2] != '':
                if int(tab_list[2]) in operon_range:
                    if tab_list[4] == "-":
                        bound_expressed.append(str(operon_list[0]) + "\t" + str(operon_list[1]) + "\t"+ "false" + "\t" + c)
    if tab_list[4] == "+": #need range going different ways depending on the strand 
        if tab_list[1] != '': #now just need to run same logic for translated transcripts if there is no output for transcripts
            for f in peak_list:
                if f in range(int(tab_list[1]), int(tab_list[2])): #first checking for internal binding, which isn't accounted for in operon search so check for it before checking whether in operon
                    bound_expressed.append(str(peak_list.index(f)+1) + "\t" + str(f) + "\t"+ "true" + "\t" +c)
                if f in operon_starts_plus: #made this list earlier so we can skip genes that are in an operon since we've already located them
                    continue
                #the code to check for genes in operons doesn't look for internal binding, however, so I skip those genes after checking for it
                if f in range((int(tab_list[1])-500), int(tab_list[1])):
                    bound_expressed.append(str(peak_list.index(f)+1) + "\t" + str(f) + "\t"+ "false" + "\t" + c)
                
        elif tab_list[0] != '': #need to make sure there's something there
            for f in peak_list: #need to run for each peak so that I can define the site to associate with the transcript
                if tab_list[3] != '': #not all the lines with a transcript start have an end to need to check for one to use for upper bound of range
                    if f in range(int(tab_list[0]), int(tab_list[3])): #also want to see if there is binding site inside of the transcript
                        bound_expressed.append(str(peak_list.index(f)+1) + "\t"+ str(f) + "\t" + "true" + "\t" + c)
                if f in operon_starts_plus:
                    continue
                if f in range((int(tab_list[0])-500), int(tab_list[0])): #-500 of + strand gene is 500 upstream
                    bound_expressed.append(str(peak_list.index(f)+1) + "\t" + str(f) + "\t"+ "false" + "\t" +c) #since I have the line that was in the range set to f can then use index to save which number site, and also save the site value
                
    if tab_list[4] == "-": #now just repeat the same logic again but for - strand
        if tab_list[1] != '':
            for f in peak_list:
                if f in range(int(tab_list[2]),int(tab_list[1])):
                    bound_expressed.append(str(peak_list.index(f)+1) + "\t" + str(f) + "\t"+"true" + "\t" + c)
                if f in operon_starts_minus:
                    continue
                if f in range(int(tab_list[1]),(int(tab_list[1])+500)):
                    bound_expressed.append(str(peak_list.index(f)+1) + "\t" + str(f) + "\t"+ "false" + "\t" + c)                
        elif tab_list[0]!= '':
            for f in peak_list:
                if tab_list[3] != '':
                    if f in range(int(tab_list[3]), int(tab_list[0])):
                        bound_expressed.append(str(peak_list.index(f)+1) + "\t" + str(f) + "\t"+ "true" + "\t" +c)
                if f in operon_starts_minus:
                    continue
                if f in range(int(tab_list[0]),(int(tab_list[0])+500)):
                    bound_expressed.append(str(peak_list.index(f)+1) + "\t" + str(f) + "\t"+ "false" + "\t" + c)
                    
rna_transcript = open(folder + "Binding_associated_transcripts.txt", "w") #"C:/RhlR_project/rna_transcripts_operons_allvsall2.txt"
for line in bound_expressed:
    rna_transcript.write(line)

#operon overlap bug fixed!

transcript.close()
rna_transcript.close()
peaks.close()
