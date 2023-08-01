alignment = open(input("Enter Path of BLAST File: "), "r")
#alignment = open("C:\RhlR_project\WGR5MYDW013-Alignment.txt", "r")

sites = []
ranges = []
strains = []
for line in alignment:
    if "Sequence ID: " in line: #this is the sequence ID of PA01 reference genome, tax id 208964
        #print(next(alignment))
        strain = line.split(" ")
        strain_id = strain[2]
        target = next(alignment) #will save the following line, which should contain the range I care about
        if "Range" in target: #I only care about lines in the alignment that have the range, and only they have "Range" in them
            duplicate = False #need to reset value each iteration so I don't just save every line
            target = target.strip("\n")
            seprange = target.split(" ") #treat each space as a separator since there are spaces on both sites of the start/end numbers
            #site_range = range(int(seprange[2]),int(seprange[4])) #then create a range using the start/end pair
            # for item in ranges: #first for loop needs to call each list in the ranges list
            #     for number in item: #then for each of those lists checking the numbers they contain
            #         if number in site_range: #can then check if any of the numbers in each of those previously found ranges are present
            #             duplicate = True #if that is the case will change the duplicate value to true so i can avoid saving repeat values
            #             break #once found can break loop
            #     if duplicate == True:
            #         break #value will change to true if I found a duplicate coordinate, so break the next highest loop too
            # if duplicate == True:
            #     continue #don't want to break the top loop, though so just skipping the rest of the current iteration
            sites.append(seprange[2] + "\t" + seprange[4] + "\t"+strain_id) #saving new lines containing the range for the binding site
            ranges.append(range(int(seprange[2]),int(seprange[4]))) #making a nested list of all sequence ranges as well

#since multiblast is comparing each of the sequences I provided in sequence to a group of strains and I am only selecting
#from one of those strains, each of the ranges I found should corresponse with one of the artifact sites and so the 
#number of lines will be the number of artifacts that I provided that aligned to PAO1

artifact_file_name = input("Enter Path/Name of Output File: ")
#artifact_file_name = "C:\RhlR_project\pao1_artifact_new2.txt"
target_file = open(artifact_file_name, "w")
for site in sites:
    target_file.write(site +"\n")


alignment.close()
target_file.close()
#Should then be able to identify which of the identified site ranges match up with the location of an MVAT binding site

artifacts = open(artifact_file_name, "r")
#mvat = open("C:/RhlR_project/Castang_MvaT_bound_regions.txt", "r")
mvat = open(input("Enter Path/Name of File Showing MvaT Binding: "), "r")
#this file uses the pa01 artifacts so must be the right one, though seems existing output file is full of duplicates
#which I must have filtered out in excel
#now need to make sure it filters out automatically as a part of this script so it's more user friendly for publishing

bound_artifacts = []
for line in mvat:
    line = line.strip("\n")
    mvat_split = line.split("\t")
    mvat_range = range(int(mvat_split[0]),int(mvat_split[1]))   
    artifacts.seek(0)
    for unit in artifacts:
        unit = unit.strip("\n")
        arti_split = unit.split("\t")
        arti_range = range(int(arti_split[0]),int(arti_split[1]))
        for number in arti_range: #checking each base in the artifact site
            if number in mvat_range: #if one of those numbers is in the range of an mvat site
                bound_artifacts.append(unit) #save that artifact
                break #then break the bottom loop to check for the next artifact
                

#with open("C:/RhlR_project/artifacts_bound_by_mvat3.txt", "w") as save_file:
with open(input("Enter Path/Name of Final Output File: "), "w") as save_file:
    for site in bound_artifacts:
        save_file.write(site+"\n")

save_file.close()
#target.close()
artifacts.close()
mvat.close()

