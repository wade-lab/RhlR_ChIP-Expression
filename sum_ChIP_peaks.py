import os   
folder = input("Enter Path of Folder Containing Necessary Files: ")
gff_list = os.listdir(folder)
#gff files
peaks = open(folder+"RhlR_peak_coordinates.txt","r")#will throw error if you don't include a slash at the path end

peak_table_table = []
peak_table_table.append(["Peaks"])
for t in peaks:
    peak_table_table.append([t])

def find_peak_coverage(item):
    gff = open(folder + item,"r")
    coverage = [0] * 7000000
    peaks = open(folder+"RhlR_peak_coordinates.txt","r")
    for x in gff:
        x = x.strip("\n") #just a precaution against unwanted line-breaks sneaking into the function
        gfftab = x.split("\t") #need to have each column as its own string so that we can splice from lines of different length
        coordinate = int(gfftab[3]) #splice the coordinate column and make it an integer so it's useful
        coverage[coordinate] += abs(int(gfftab[5]))
    peak_table = [item]
    for y in peaks:
        y = y.strip("\n") #another line to clear any errant line breaks
        site = int(y) #can only do math with numbers!
        site_above = site + 101 #python counts starting from 0 so need to account for peak center being offset by 1 (not that it should make much difference)
        site_below = site -100
        site_coverage = range(site_below, site_above) #need the range above and below the site of interest
        coverage_area = [] 
        for z in site_coverage:
            coverage_area.append(coverage[z]) #these are essential lines, using the range make a list of all the value in that range
        coverage_size = sum(coverage_area) #can then find the sum of that list to find the sum of coverage in that range
        coverage_normalized = (1000000 * coverage_size)/ (sum(coverage)) #need to modify the value by bases/million relative to the overall coverage of the file
        peak_table.append(str(coverage_normalized)) #now can just add the values into the target file
    length = len(peak_table_table) #want to have a length for the table of tables
    for n in range(length): #then can call each of the index numbers for that range
        peak_table_table[n].append(peak_table[n]) #and add to each of those sub-tables the same line value from the table for a given strain
        #this means I'll have a table with sub tables containing the value for each row for every sample run
    peaks.close()
    gff.close()
    print(item + "done")

for file in gff_list: #now that the full function is defined can run it for each detected gff
    find_peak_coverage(file)

peak_sizes = open(folder+"Peak_size_output.txt","w")
##peak_sizes.close() #need this if there isn't already a file, will overwrite anything already present if there is a file
##peak_sizes = open("C:\\Users\\Nick\\Documents\\peak_size_table.txt","a")
for whatever in peak_table_table: #then can write to file each table
    for something in whatever: #first need to call each value in each table and add tabs
        peak_sizes.write(str(something) + "\t")
    peak_sizes.write("\n") #and for each of the tables (which are my lines) add a line break to separate them
    
        
#instead add each peak as a list so append([x]) and then append the corresponding number in the list to that!!!

peak_sizes.close()

print("super done")





