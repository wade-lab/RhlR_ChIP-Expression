#import scipy #use this if including bionmial test
genome_length = int(input("Enter Length of Genome: ")) #6537648 bp according to https://pseudomonas.com/strain/show?id=109 for the PA14 reference genome
gene_list = open(input("Enter Path for Gff List of Genes: "), "r") #C:\RhlR_project\RhlR_ChIP_gffs\PA14_genes.gff
base_list = [0] * genome_length #making a list the length of the genome
genome_range = range(0, genome_length)
for gene in gene_list:
    gene_tabs = gene.split("\t")
    if len(gene_tabs) <= 1: #a lot of files seem to have a last line with a sinle line break that python tries to run this on
        continue #so will skip line that aren't usable
    gene_range = range(int(gene_tabs[3]), int(gene_tabs[4])) #creating a list that contains a number for every base in the site
    for base in gene_range: #for each base in the range associated with that gene
        if (base-1) in genome_range:  #if that gene 
            base_list[(base-1)] = 1 #set the object in the list of genome length to 1 (from 0)-python counts starting from 0 instead of 1, so need to account for that
intergenic = base_list.count(0) #count how many objects (bases) in the genome length list are not associated with a gene
intragenic = base_list.count(1) #also count how many objects (bases) in the genome length list are associated with a gene
proportion_intergenic = (intergenic / genome_length) #the quotient of objects still set to 0 divided by genome length gives the proportion of the genome which is coding
print("Genome Length = " + str(genome_length))
print("Intergenic = " + str(intergenic))
print("Intragenic = " + str(intragenic))
print("Proportion Intergenic = " + str(proportion_intergenic))
k = input("Enter Number of Intergenic Sites: ")
n = input("Enter Total Number of Sites: ")
relative_intergenic = (int(k)/int(n))/proportion_intergenic
#prob_binom = scipy.stats.binomtest(k=k, n=n, p = proportion_intergenic) #running a binomial test for the probability of 40 bases to be intergenic with 168 tries given the probability of a base to be intergenic
#k=32, n=40
#this is essentially measuring the chance that there would be 32 individual bases (so the center of the chip peaks) positioned intergenically

#print("p = " + str(prob_binom.pvalue))
print("Enrichment of Intergenic Binding from Expected = "+str(relative_intergenic))
gene_list.close()