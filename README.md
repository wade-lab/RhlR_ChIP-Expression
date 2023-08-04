Main Data Workflow:
---------------------
-  List of genotype of strains in various datasets:
- 151 = rhlr deletion
- 154 = rhli deletion
- 158 = rhli/r double deletion (not used)
- 222 = wt
- 271 = pqse deletion (used by RNA-seq)
- 278 = rhli/pqse double deletion
- 544 = pqse non-interacting
- 561 = pqse deletion (used by ChIP)
- 568 = pqsE catalytically inactive (D73A)
Prepare the initial table of ChIP-seq data:
- Download:
   - sum_ChIP_peaks.py
   - RhlR_peak_coordinates.txt
- Extract the Zipped ChIP-seq_Gff folder into a new folder. This contains .gff form data from an alignment of ChIP-seq sequencing output
- Run the sum_ChIP-peaks.py file
- Input when prompted, this should include the path to the extracted folder, and the path/name of RhlR_peak_coordinates.txt
  - This will produce a new file titled Peak_size_output.txt in the file containing the ChIP .gffs. This will be a table containing peak values
    for each peak and each strain file produced.
- Open the file in Excel, create new columns from the average of the two replicates for each strain number, and save the columns (including the peak coordinate column to a table in a new file, name it peaks_average.txt
  -This will be a key table for future steps
Prepare the initial table of RNA-seq data:
- Download:
   - Rockhopper (https://cs.wellesley.edu/~btjaden/Rockhopper/)
   - merge_gzip.py
   - binding_transcript_matching.py
   - RNA_seq_Gff.zip
- Extract the folder in RNA_seq_Gff.zip
- For each folder in run 1 and run 3 plus or minus containing RNA-seq duplicates, run merge_gzip and input that folder when prompted
- Once all duplicates are merged (renaming merged files appropriately), move all from each run with all other merged files from the same run (mainly for convenience)
- Run Rockhopper, with one experiment for each of the strains, simultaneously. Add both pair ended reads for each replicate, but uncheck the strand specific option in the parameter options. Use the Pseudomonas aeruginosa UCBPP-PA14 replicon for a reference.
- Rockhopper will generate two relevant files, NC_008463_operons.txt and NC_008463_transcripts.txt. Ensure they are named NC_008463_operons_allvsall.txt and NC_008463_operons_allvsall.txt, move them to a folder with RhlR_peak_coordinates.txt
- Run binding_transcript_matching.py, when prompted input the path to the folder containing the files from the previous step
- this will provide two output files, one listing all operons with an internal or upstream binding site, one with the final table containing all genes with an internal or upstream binding site or in an operon with an upstream binding site
Processing of Data
- Download:
   - rhlr_project_r_code.R
   - rhlr_associated_expression.py
- Open rhlr_project_r_code.R, set the working directory to the folder containing peaks_average.txt, and run lines 13-18 to generate peaks_average_real.txt in the working directory
- Run rhlr_associated_expression.py, input prompted paths
- From the output file, create a new version containing columns containing site locations, expression values for each site, and qvalue 222 vs 151
- Move the new file to the R working directory
- Run the R file, it will generate many of the scatterplots used in this study

Additional Code:
----------------------
Determine the presence of MvaT binding sites overlapping with ChIP-sites not affected by RhlR deletion
- Download:
   - Castang_MvaT_bound_regions.txt (this data is modified from Castang et al. 2008 doi: 10.1073/pnas.0808215105 supplementary table 1)
   - Galaxy105-[222_comb_peaks.fasta]_artifact.txt
   - find_sequence_from_multiblast.py
- 
