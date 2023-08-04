#/Users/jtw03/Desktop/Rockhopper_Results/genomeBrowserFiles/
#/Users/jtw03/Desktop/Rockhopper/Rockhopper/Rockhopper_Results/genomeBrowserFiles/
#/Users/jtw03/Documents/Projects/NET-seq analysis/
#/Users/jtw03/Documents/My papers/PhoB/Fastq files/s70 wiggles/

## Make pile up lists

wiggle_plus=input('Enter file path and name for plus strand .wig file:')
wiggle_minus=input('Enter file path and name for minus strand .wig file:')
normalize=input('Do you want to normalize total number of reads to 100 million? Answer yes or no:')
#read_length=int(input('What length are the mapped reads?'))
gff_name=input('Enter file path and name for .gff file (name must end with .txt or .gff):')
gff_label=input('Enter track label for .gff file:')
cutoff=int(input('Enter the lowest value to be plotted in .gff file:'))
interval=int(input('Enter the interval between plotted positions (ie. 2 for every other position):'))


f=open(wiggle_plus,'r')

f.readline()
f.readline()

plus=[]

for a in f:
    plus.append(float(a.split()[0]))

f.close()

f=open(wiggle_minus,'r')

f.readline()
f.readline()

minus=[]

for a in f:
    minus.append(float(a.split()[0]))

f.close()

##Normalize to 100 million reads

if normalize=='yes':
    reads=max(max(plus),max(minus))
    norm=100000.0/reads
    for x in range(len(plus)):
        plus[x]=int(plus[x]*norm)
        minus[x]=int(minus[x]*norm)

print ('Rep 1 pile-up complete')

## Write pile-up gff

gff=open(gff_name,'w')

for x in range (0,len(plus),interval):
    if plus[x]>cutoff:
        gff.write('NA\tAgilent\t')
        gff.write(gff_label)
        gff.write('\t')
        gff.write(str(x+1))
        gff.write('\t')
        gff.write(str(x+1))
        gff.write('\t')
        gff.write(str(plus[x]))
        gff.write('\t.\t.\t.\n')
    

for x in range (0,len(minus),interval):
    if minus[x]>cutoff:
        gff.write('NA\tAgilent\t')
        gff.write(gff_label)
        gff.write('\t')
        gff.write(str(x+1))
        gff.write('\t')
        gff.write(str(x+1))
        gff.write('\t')
        gff.write(str(minus[x]*-1))
        gff.write('\t.\t.\t.\n')

gff.close()
