<<<<<<< HEAD
##RNA Seq differential expression analysis
require(readr)
require(plyr)
require(data.table)
require(dplyr)
rna_seq = read_delim('~/Documents/RNA_Seq/RNA_seq_Analysis_Dec21/tpm.diff.cell.diff.stim.txt', delim = "\t")
names(rna_seq)
##make delta expression by subtracting the expression from pandemrix stimulation the control expression 

=======
##RNA Seq differential expression analysis
>>>>>>> 63055afdfdadd9d6e28467ffb0eb39417ee088f3
