#R version 4.2.2
#will need to use python script for summing chip peaks to generate the chip data used here
#the average of the chip coverage output for the two duplicates of the chip-seq was used here
library(ggplot2)
library(ggrepel)
library(scales)
#First, change the path in the below function to the location of this file and its companion files
setwd("C:/RhlR_project/Figures/Unlabeled/")
#Ensure the following three files are in the directory: all_transcripts_expression.txt, peaks_average.txt, rhlr_associated_expression_selective.txt
#Then run this file
# Data Set-Up ------------------------------------------------------------- 
#First need to select only sites that are associated with RhlR binding
peaks_avg = read.table("peaks_average.txt", sep = "\t", header= TRUE) #Load the file containin ChIP data
rhlr_ratio = (peaks_avg[6]/peaks_avg[2]) #Next will find the ratio of WT enrichment / rhlr deletion enrichment and store in a new file
real_peaks = subset(rhlr_ratio, subset = X222.avg > 1) #Deletion enrichment that is lower than WT enrichment will give numbers >1 which I subset for
real_rows = rownames(real_peaks) #then saving the names of this subset to a vector
peaks_avg_real = peaks_avg[real_rows,] #subset row names are also the index numbers for the rows, which I can use to extrac the rows of the original file that have enrichment reduced by rhlr deletion
#peaks_avg and peaks_avg_real will both be used to make the ChIP-seq plots
#peaks_avg_real will be important as it will provide the like of "real" sites for python to select associated genes for

#run python script operon_transcript_expression.py for relavant rockhopper rna-seq alignment .txt
#this will generate a table containing all genes associated with any binding site
#create a new text file by copying the line from the python output file that contain only the expression data
all_expression = read.table("all_transcripts_expression.txt", header= TRUE,
                            sep = "\t")
#additionally, run the file "site_associated_expression.py"
#this will make a file outlining expression for all the transcripts associated 'real' rhlr binding 
#binding sites (including genes in operon with an upstream rhlr site and genes with an internal binding site)
#the rockhopper output has a LOT of columns which R has a hard time parsing with read.table, will need to
#Using a file with smaller number of columns, loaded in the next section call "rhlr_associated_expression_selective.txt"

# Gene of Interest Frames -------------------------------------------------
#construction of a dataframe containing data for labeling geom of ChIP data
#Unfortunately a few genes are repeated in the main expression data table so they can't be
#used for line names, which I could just use as a vector to make labels from the same table
#I manually went through the dataframe to find the right indexes for these genes for both tables :(
rhlr_associated_expression = read.table("rhlr_associated_expression_selective.txt", header= TRUE,
                                        sep = "\t")
#vector of the genes I label later from first to last in genome c("phzA1", "lasB","rhlA", "lecB", "hcnA", "phzA2")
peaks_avg_real_labels = peaks_avg_real[c(6,9,12,15,25,28),]
rownames(peaks_avg_real_labels) <- c("phzA1", "lasB","rhlA", "lecB", "hcnA", "phzA2")

rhlr_associated_expression_labels = rhlr_associated_expression[c(16,26,33,38,67,81),]
rownames(rhlr_associated_expression_labels) <- c("phzA1", "lasB","rhlA", "lecB", "hcnA", "phzA2")

# RNA-seq Graphs ----------------------------------------------------------
#Showing the expression associated with sites
  #pqse plot
plotobject = ggplot(all_expression, aes(x=Expression.222, y=Expression.271 ), 
                    element_text(family = "TT Arial", size = 14))
pqse_expression = plotobject + geom_point(color= "gray") + 
  scale_x_log10(
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))  +
  xlab("WT Expression") + ylab("Δ*pqsE* Expression") + 
  theme_classic(base_size = 14) + theme( axis.title.y = ggtext::element_markdown())+
  geom_point(data = (rhlr_associated_expression), aes(x=Expression.222, y=Expression.271),
             color="purple", size = 1.5) #+
  #geom_text_repel(data = rhlr_associated_expression_labels,
   #         aes(Expression.222, Expression.271, label = rownames(rhlr_associated_expression_labels)),
    #        fontface = "italic", direction = "y"))
ggsave("pqse_expression.tiff", pqse_expression, height = 4, width = 4, dpi = 200)
        #comment out the geom_text function to remove point labels
  #rhlr plot
plotobject = ggplot(all_expression, aes(x=Expression.222, y=Expression.151 ),
                    element_text(family = "TT Arial", size = 14))
rhlr_expression = plotobject + geom_point(color= "gray") + 
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)))  +
  theme_classic(base_size = 14) + theme( axis.title.y = ggtext::element_markdown())+
  xlab("WT Expression") + ylab("Δ*rhlR* Expression") + 
  geom_point(data = (rhlr_associated_expression_only), aes(x=Expression.222, y=Expression.151), color="red", size = 1.5) #+
  #geom_text_repel(data = rhlr_associated_expression_labels, 
            #aes(Expression.222, Expression.151, label = rownames(rhlr_associated_expression_labels)),
            #fontface = "italic", direction = "y")
ggsave("rhlr_expression.tiff", rhlr_expression, height = 4, width = 4, dpi = 200)
  #rhli plot
plotobject = ggplot(all_expression, aes(x=Expression.222, y=Expression.154 ),
                    element_text(family = "TT Arial", size = 14))
rhli_expression = plotobject + geom_point(color= "gray") + 
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)))  +
  theme_classic(base_size = 14) + theme( axis.title.y = ggtext::element_markdown())+
  xlab("WT Expression") + ylab("Δ*rhlI* Expression") + 
  geom_point(data = (rhlr_associated_expression_only), aes(x=Expression.222, y=Expression.154), color="blue", size = 1.5) #+
  #geom_text_repel(data = rhlr_associated_expression_labels, 
            #aes(Expression.222, Expression.154, label = rownames(rhlr_associated_expression_labels)),
            #fontface = "italic", direction = "y"))
ggsave("rhli_expression.tiff", rhli_expression, height = 4, width = 4, dpi = 200)
  #pqse catalytically dead
plotobject = ggplot(all_expression, aes(x=Expression.222, y=Expression.568 ),
                    element_text(family = "TT Arial", size = 14))
cat_ded_pqse_expression = plotobject + geom_point(color= "gray") + 
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)))  +
  theme_classic(base_size = 14) + theme( axis.title.y = ggtext::element_markdown())+
  xlab("WT Expression") + ylab("*pqsE* (D73A) Expression") + 
  geom_point(data = (rhlr_associated_expression_only), aes(x=Expression.222, y=Expression.568), color="cyan", size = 1.5) #+
  #geom_text_repel(data = rhlr_associated_expression_labels, 
            #aes(Expression.222, Expression.568, label = rownames(rhlr_associated_expression_labels)),
            #fontface = "italic", direction = "y")
ggsave("cat_ded_pqse_expression.tiff", cat_ded_pqse_expression, height = 4, width = 4, dpi = 200)
  #pqse NI vs pqse deletion
plotobject = ggplot(all_expression, aes(x=Expression.544, y=Expression.271 ),
                    element_text(family = "TT Arial", size = 14))
pqse_NI_expression = plotobject + geom_point(color= "gray") + 
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)))  +
  theme_classic(base_size = 14) +theme( axis.title.y = ggtext::element_markdown(), axis.title.x = ggtext::element_markdown())+
  xlab(" *pqsE* Non-Interacting Expression") + ylab("Δ*pqsE* Expression") + 
  geom_point(data = (rhlr_associated_expression_only), aes(x=Expression.544, y=Expression.271), color="pink", size = 1.5) #+
  #geom_text_repel(data = rhlr_associated_expression_labels, 
            #aes(Expression.544, Expression.271, label = rownames(rhlr_associated_expression_labels)),
            #fontface = "italic", direction = "y")
#need to update table with 544 data for this graph!!!
ggsave("pqse_NI_expression.tiff", pqse_NI_expression, height = 4, width = 4, dpi = 200)

# ChIP-seq Graphs ---------------------------------------------------------
#Making plot of binding associated with 
  #pqse
plotobject = ggplot(peaks_avg_real, aes(x=X222.avg, y=X561.avg ),
                    element_text(family = "TT Arial", size = 14))
dpqse_enrichment = plotobject + geom_point(color="purple") + geom_abline() +
  theme_classic(base_size = 14)+ theme( axis.title.y = ggtext::element_markdown())+
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x, n=6),
    labels = trans_format("log10", math_format(10^.x)), limits = c(1,100000)) + 
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x, n=6),
    labels = trans_format("log10", math_format(10^.x)), limits = c(1,100000))  +
  xlab("WT Enrichment") + ylab("Δ*pqsE* Enrichment") #+
  #geom_text_repel(data = peaks_avg_real_labels, 
            #aes(X222.avg, X561.avg, label = rownames(peaks_avg_real_labels)),
            #fontface = "italic", direction = "y")
ggsave("dpqse_enrichment.tiff", dpqse_enrichment, height = 4, width = 4, dpi = 200)
  #rhlr
plotobject = ggplot(peaks_avg, aes(x=X222.avg, y=X151.avg ),
                    element_text(family = "TT Arial", size = 14))
drhlr_enrichment = plotobject + geom_point() + geom_abline() +theme_classic(base_size = 14)+
  theme( axis.title.y = ggtext::element_markdown())+
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x, n=6),
    labels = trans_format("log10", math_format(10^.x)), limits = c(1,100000)) + 
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x, n=6),
    labels = trans_format("log10", math_format(10^.x)), limits = c(1,100000))+
  xlab("WT Enrichment") + ylab("Δ*rhlR* Enrichment")+
  geom_point(data=peaks_avg_real,aes(x=X222.avg,y=X151.avg),color = "red") #+
  #geom_text_repel(data = peaks_avg_real_labels, 
            #aes(X222.avg, X151.avg, label = rownames(peaks_avg_real_labels)),
            #fontface = "italic", direction = "y")
ggsave("drhlr_enrichment.tiff", drhlr_enrichment, height = 4, width = 4, dpi = 200)
  #rhli
plotobject = ggplot(peaks_avg_real, aes(x=X222.avg, y=X154.avg ),
                    element_text(family = "TT Arial", size = 14))
drhli_enrichment = plotobject + geom_point(color="blue") + geom_abline() +theme_classic(base_size = 14)+
  theme( axis.title.y = ggtext::element_markdown())+
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x, n=6),
    labels = trans_format("log10", math_format(10^.x)), limits = c(1,100000)) + 
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x, n=6),
    labels = trans_format("log10", math_format(10^.x)), limits = c(1,100000))  +
  xlab("WT Enrichment") + ylab("Δ*rhlI* Enrichment") #+
  #geom_text_repel(data = peaks_avg_real_labels, 
            #aes(X222.avg, X154.avg, label = rownames(peaks_avg_real_labels)),
            #fontface = "italic", direction = "y")
ggsave("drhli_enrichment.tiff", drhli_enrichment, height = 4, width = 4, dpi = 200)
#rhli/pqse double deletion
plotobject = ggplot(peaks_avg_real, aes(x=X222.avg, y=X278.avg ),
                    element_text(family = "TT Arial", size = 14))
rhli_pqse_dd = plotobject + geom_point(color="orange") + geom_abline() +
  theme_classic(base_size = 14)+
  theme( axis.title.y = ggtext::element_markdown())+
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x, n=6),
    labels = trans_format("log10", math_format(10^.x)), limits = c(1,100000)) + 
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x, n=6),
    labels = trans_format("log10", math_format(10^.x)), limits = c(1,100000))  +
  xlab("WT Enrichment") + ylab("Δ*rhlI* Δ*pqsE* Enrichment") #+
  #geom_text_repel(data = peaks_avg_real_labels, 
            #aes(X222.avg, X278.avg, label = rownames(peaks_avg_real_labels)),
            #fontface = "italic", direction = "y")
ggsave("drhli_pqse.tiff", rhli_pqse_dd, height = 4, width = 4, dpi = 200)
  #pqse catalytically dead mutant
plotobject = ggplot(peaks_avg_real, aes(x=X222.avg, y=X568.avg ),
                    element_text(family = "TT Arial", size = 14))
catded_enrichment = plotobject + geom_point(color="cyan") + geom_abline() +theme_classic(base_size = 14)+
  theme( axis.title.y = ggtext::element_markdown(), axis.title.x = ggtext::element_markdown())+
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x, n=6),
    labels = trans_format("log10", math_format(10^.x)), limits = c(1,100000)) + 
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x, n=6),
    labels = trans_format("log10", math_format(10^.x)), limits = c(1,100000))  +
  xlab("WT Enrichment") + ylab("*pqsE* (D73A) Enrichment") #+
  #geom_text_repel(data = peaks_avg_real_labels, 
            #aes(X222.avg, X568.avg, label = rownames(peaks_avg_real_labels)),
            #fontface = "italic", direction = "y")
ggsave("catded_enrichment.tiff", catded_enrichment, height = 4, width = 4, dpi = 200)
  #pqse NI but against pqse deletion instead
plotobject = ggplot(peaks_avg_real, aes(x=X544.avg, y=X561.avg ),
                    element_text(family = "TT Arial", size = 14))
pqse_NI_enrichment = plotobject + geom_point(color="pink") + geom_abline() +theme_classic(base_size = 14)+
  theme( axis.title.y = ggtext::element_markdown(),axis.title.x = ggtext::element_markdown())+
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x, n=6),
    labels = trans_format("log10", math_format(10^.x)), limits = c(1,100000)) + 
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x, n=6),
    labels = trans_format("log10", math_format(10^.x)), limits = c(1,100000))  +
  xlab("*pqsE* Non-Interacting Enrichment") + ylab("Δ*pqsE* Enrichment") #+
  #geom_text_repel(data = peaks_avg_real_labels, 
            #aes(X544.avg, X561.avg, label = rownames(peaks_avg_real_labels)),
            #fontface = "italic", direction = "y")
ggsave("pqseNI_enrichment.tiff", pqse_NI_enrichment, height = 4, width = 4, dpi = 200)
 
