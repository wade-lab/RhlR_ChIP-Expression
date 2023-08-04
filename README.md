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
- 
