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
rna_seq2 = arrange(rna_seq2, desc(DeltaCD8))
rna_seq2
library(limma)
tot_seq = select(rna_seq, gene_id, cd4_pandemrix.tpm, cd4_control.tpm, cd8_pandemrix.tpm, cd8_control.tpm, cd14_pandemrix.tpm, cd14_control.tpm, cd19_pandemrix.tpm, cd19_control.tpm, nk_pandemrix.tpm, nk_control.tpm)
tot_seq = tot_seq[-1,]
tot_seq = as.data.frame(tot_seq)
rownames(tot_seq) = tot_seq$gene_id
#tot_seq$gene_id = NULL
##log convert the columns to account for the skewness
totlog=apply(tot_seq, 2, FUN = function(x)(log(x)))
class(totlog)
totlog = as.data.frame(totlog)
totlog$gene_id = tot_seq$gene_id #assign gene names back to this frame
##inf values should be assigned the lowest in each sample
for (i in colnames(totlog[1:10])){
  totlog[i][totlog[i]== -Inf] = totlog[i][totlog[i]!= -Inf] %>% min()
}
require(ggplot2)
plot1=ggplot(totlog, aes(cd4_pandemrix.tpm, cd4_control.tpm, label = totlog$gene_id))+geom_point(alpha = 1/5)+geom_text(angle=60, check_overlap = T, hjust = -0.2, nudge_x = 0.1)+theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size =15), axis.title.x = element_text(size = 18, face = "bold"), axis.title.y = element_text(size = 18, face = "bold"))
ggsave('~/Documents/RNA_Seq/RNA_seq_Analysis_Dec21/plot5.pdf', plot = plot5, width = 13.6, height = 8, dpi = 800, units = "in")
plot2=ggplot(totlog, aes(cd8_pandemrix.tpm, cd8_control.tpm, label = totlog$gene_id))+geom_point(alpha = 1/5)+geom_text(angle=60, check_overlap = T, hjust = -0.2, nudge_x = 0.1)+theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size =15), axis.title.x = element_text(size = 18, face = "bold"), axis.title.y = element_text(size = 18, face = "bold"))
plot3=ggplot(totlog, aes(cd19_pandemrix.tpm, cd19_control.tpm, label = totlog$gene_id))+geom_point(alpha = 1/5)+geom_text(angle=60, check_overlap = T, hjust = -0.2, nudge_x = 0.1)+theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size =15), axis.title.x = element_text(size = 18, face = "bold"), axis.title.y = element_text(size = 18, face = "bold"))
plot4=ggplot(totlog, aes(cd14_pandemrix.tpm, cd14_control.tpm, label = totlog$gene_id))+geom_point(alpha = 1/5)+geom_text(angle=60, check_overlap = T, hjust = -0.2, nudge_x = 0.1)+theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size =15), axis.title.x = element_text(size = 18, face = "bold"), axis.title.y = element_text(size = 18, face = "bold"))
plot5=ggplot(totlog, aes(nk_pandemrix.tpm, nk_control.tpm, label = totlog$gene_id))+geom_point(alpha = 1/5)+geom_text(angle=60, check_overlap = T, hjust = -0.2, nudge_x = 0.1)+theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size =15), axis.title.x = element_text(size = 18, face = "bold"), axis.title.y = element_text(size = 18, face = "bold"))
##heatmaps for top 50 upreg and 50 downreg genes
setDT(totlog)
rnaseq_melt=melt(totlog, id.vars = "gene_id")
write_csv(rna_seq1, file ='~/Documents/RNA_Seq/')
