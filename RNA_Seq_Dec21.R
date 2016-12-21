##RNA Seq differential expression analysis
require(readr)
require(plyr)
require(data.table)
require(dplyr)
rna_seq = read_delim('~/Documents/RNA_Seq/RNA_seq_Analysis_Dec21/tpm.diff.cell.diff.stim.txt', delim = "\t")
=======