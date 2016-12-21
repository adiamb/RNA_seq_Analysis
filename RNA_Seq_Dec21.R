##RNA Seq differential expression analysis
require(readr)
require(plyr)
require(data.table)
require(dplyr)
rna_seq = read_delim('~/Documents/RNA_Seq/RNA_seq_Analysis_Dec21/tpm.diff.cell.diff.stim.txt', delim = "\t")
names(rna_seq)
##create delta expression from subtracting the values of controls from px stimulations
rna_seq$DeltaCD14 = rna_seq$cd14_pandemrix.tpm - rna_seq$cd14_control.tpm
rna_seq$DeltaCD4 = rna_seq$cd4_pandemrix.tpm - rna_seq$cd4_control.tpm
rna_seq$DeltaCD8 = rna_seq$cd8_pandemrix.tpm - rna_seq$cd8_control.tpm
rna_seq$DeltaCD19 = rna_seq$cd19_pandemrix.tpm - rna_seq$cd19_control.tpm
rna_seq$DeltaNK = rna_seq$nk_pandemrix.tpm - rna_seq$nk_control.tpm
##subset the dataframe to include only the deltas
rna_seq2 = select(rna_seq, gene_id, DeltaCD4, DeltaCD8, DeltaCD19, DeltaCD14, DeltaNK)
#check for missing values
apply(rna_seq2, 2, FUN = function(x)(sum(is.na(x))))
##remove row1 _ missing value
rna_seq2 = rna_seq2[-1,]     
rna_seq2 = arrange(rna_seq2, desc(DeltaCD4, DeltaCD8))
