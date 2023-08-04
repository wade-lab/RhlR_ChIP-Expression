This file outlines the process needed to replicate the data and several graphs in the paper.

Prepare the initial table of ChIP-seq data:
- Download:
   - sum_ChIP_peaks.py
   - RhlR_peak_coordinates.txt
- Extract the Zipped ChIP-seq_Gff folder into a new folder. This contains .gff form data from an alignment of ChIP-seq sequencing output
- Run the sum_ChIP-peaks.py file
- Input when prompted, this should include the path to the extracted folder, and the path/name of RhlR_peak_coordinates.txt
  - This will produce a new file titled Peak_size_output.txt in the file containing the ChIP .gffs. This will be a table containing peak values
    for each peak and each strain file produced.
- Open the file in Excel, create new columns from the average of the two replicates for each strain number, and save the columns (including the peak coordinate column to a table in a new file
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
