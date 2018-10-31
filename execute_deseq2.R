rm(list=ls())
library("DESeq2")

setwd(file.path("//salt","wheat_rnaseq", "PeterBuchner_RNA", "timeseries_exp", "Analysis"))
#counts <- read.table("TGAC_counts.tab", row.names = 1, comment="", as.is=T)
counts <- read.table("IWGSC_counts.tab", row.names = 1, comment="", as.is=T)
sampleInfo <- data.frame(row.names=c("Node100T1_A","Node100T1_B","Node100T1_C",
                                     "Node100T2_A","Node100T2_B","Node100T2_C",
                                     "Node100T3_A","Node100T3_B","Node100T3_C",
                                     "Node100T4_A","Node100T4_B","Node100T4_C",
                                     "Node100T5_A","Node100T5_B","Node100T5_C",
                                     "Node100T6_A","Node100T6_B","Node100T6_C",
                                     "Node100T7_A","Node100T7_B","Node100T7_C",
                                     "Node100T8_A","Node100T8_B","Node100T8_C",
                                     "Node200T1_A","Node200T1_B","Node200T1_C",
                                     "Node200T2_A","Node200T2_B","Node200T2_C",
                                     "Node200T3_A","Node200T3_B","Node200T3_C",
                                     "Node200T4_A","Node200T4_B","Node200T4_C",
                                     "Node200T5_A","Node200T5_B","Node200T5_C",
                                     "Node200T6_A","Node200T6_B","Node200T6_C",
                                     "Node200T7_A","Node200T7_B","Node200T7_C",
                                     "Node200T8_A","Node200T8_B","Node200T8_C"),
                         time=factor(c(rep(c(rep("T1",3),rep("T2",3),rep("T3",3),rep("T4",3),
                                             rep("T5",3),rep("T6",3),rep("T7",3),rep("T8",3)),2))),
                         N=factor(c(rep("100N",24), rep("200N", 24))))

colnames(counts) <- rownames(sampleInfo)

ddsMat_T5 <- DESeqDataSetFromMatrix(countData = counts, colData=sampleInfo, design = ~ N + N:time) #Model to find genes activating/deactivating at T5
ddsMat_N_effect <- DESeqDataSetFromMatrix(countData = counts, colData=sampleInfo, design = ~ N + time) #Model to find genes having N and time dependencies regardless of the other
ddsMat_interaction <- DESeqDataSetFromMatrix(countData = counts, colData=sampleInfo, design = ~ N * time) #Model to find genes with N and time interacting (e.g. time shifts) 

dds_T5 <- ddsMat_T5[rowSums(counts(ddsMat_T5)) > 1,]
dds_N_effect <- ddsMat_N_effect[rowSums(counts(ddsMat_N_effect)) > 1,]
dds_interaction <- ddsMat_interaction[rowSums(counts(ddsMat_interaction)) > 1,]

dds_interaction <- DESeq(dds_interaction)

#Plot the expression of the three highly expressed Nitrate Transporters identified by Peter
source(file.path("//salt","wheat_rnaseq","NRgeneV1","parse_iwgsc_homeologs_file.R"))
gene_homeologs <- find_homeologs("TraesCS7A01G301700") %>% gsub("01G", "02G", .)

sample_factors_and_counts <- sampleInfo
sample_factors_and_counts$homeologA <- counts(dds_interaction, normalized = T)[gene_homeologs[1],]
sample_factors_and_counts$homeologB <- counts(dds_interaction, normalized = T)[gene_homeologs[2],]
sample_factors_and_counts$homeologD <- counts(dds_interaction, normalized = T)[gene_homeologs[3],]
sample_factors_and_counts$time_daa <- c(rep(0,3),rep(7,3),rep(14,3),rep(21,3),rep(25,3),rep(28,3),
                                        rep(31,3),rep(34,3), rep(0,3),rep(7,3),rep(14,3),rep(21,3),
                                        rep(28,3),rep(32,3),rep(35,3),rep(38,3))

#Plot with time points (T1, T2, ..., T8) in the x axis
expression_plot_A <- ggplot(sample_factors_and_counts, aes(x = time, y = homeologA, color = N, group = N)) +
                      geom_point() +
                      expand_limits(y=0) +
                      ggtitle(gene_homeologs[1]) +
                      theme(plot.title = element_text(hjust = 0.5))

expression_plot_B <- ggplot(sample_factors_and_counts, aes(x = time, y = homeologB, color = N, group = N)) +
                      geom_point() +
                      expand_limits(y=0) +
                      ggtitle(gene_homeologs[2]) +
                      theme(plot.title = element_text(hjust = 0.5))

expression_plot_D <- ggplot(sample_factors_and_counts, aes(x = time, y = homeologD, color = N, group = N)) +
                      geom_point() +
                      expand_limits(y=0) +
                      ggtitle(gene_homeologs[3]) +
                      theme(plot.title = element_text(hjust = 0.5))

x11();grid.arrange(expression_plot_A, expression_plot_B, expression_plot_D, nrow=1, ncol=3)

#Plot with days after anthesis in the x axis
expression_plot_A <- ggplot(sample_factors_and_counts, aes(x = time_daa, y = homeologA, color = N, group = N)) +
                      geom_point() +
                      expand_limits(y=0) +
                      ggtitle(gene_homeologs[1]) +
                      theme(plot.title = element_text(hjust = 0.5)) +
                      scale_x_continuous("time (DAA)", breaks=c(0,7,14,21,25,28,31,32,34,35,38), labels=c("0","7","14","21","25","28","31","32","34","35","38"), minor_breaks = NULL)

expression_plot_B <- ggplot(sample_factors_and_counts, aes(x = time_daa, y = homeologB, color = N, group = N)) +
                      geom_point() +
                      expand_limits(y=0) +
                      ggtitle(gene_homeologs[2]) +
                      theme(plot.title = element_text(hjust = 0.5)) +
                      scale_x_continuous("time (DAA)", breaks=c(0,7,14,21,25,28,31,32,34,35,38), labels=c("0","7","14","21","25","28","31","32","34","35","38"), minor_breaks = NULL)

expression_plot_D <- ggplot(sample_factors_and_counts, aes(x = time_daa, y = homeologD, color = N, group = N)) +
                      geom_point() +
                      expand_limits(y=0) +
                      ggtitle(gene_homeologs[3]) +
                      theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous("time (DAA)", breaks=c(0,7,14,21,25,28,31,32,34,35,38), labels=c("0","7","14","21","25","28","31","32","34","35","38"), minor_breaks = NULL)

x11();grid.arrange(expression_plot_A, expression_plot_B, expression_plot_D, nrow=1, ncol=3)

#Find genes with time effect
dds_time_effect <- DESeq(dds_N_effect, test="LRT", reduced = ~ N)
deg_time_effect <- subset(results(dds_time_effect, alpha=0.05), padj<0.05)
head(deg_time_effect[order(deg_time_effect$padj),])

gene_id <- "TRIAE_CS42_1AL_TGACv1_000064_AA0001950"
deg <- plotCounts(dds_time_effect, gene_id, intgroup = c("N","time"), returnData = TRUE)
ggplot(deg, aes(x = time, y = count, color = N, group = N)) +
  geom_point() + geom_smooth(se = FALSE, method = "loess") + scale_y_log10() +
  ggtitle(paste(gene_id, " (padj = ",format(deg_time_effect[rownames(deg_time_effect)==gene_id,]$padj, nsmall=2) ,")")) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("time_effect_deg.pdf", device = "pdf", path = "results")
ggsave("time_effect_deg.png", device = "png", path = "results", dpi=300)

#Find genes with nitrogen effect
dds_N_effect <- DESeq(dds_N_effect, test="LRT", reduced = ~ time)
deg_N_effect <- subset(results(dds_N_effect, alpha=0.05), padj<0.05)
head(deg_N_effect[order(deg_N_effect$padj),])

gene_id <- "TRIAE_CS42_3DL_TGACv1_249674_AA0854070"
deg <- plotCounts(dds_time_effect, gene_id, intgroup = c("N","time"), returnData = TRUE)
ggplot(deg, aes(x = time, y = count, color = N, group = N)) +
  geom_point() + geom_smooth(se = FALSE, method = "loess") + scale_y_log10() +
  ggtitle(paste(gene_id, " (padj = ",format(deg_N_effect[rownames(deg_N_effect)==gene_id,]$padj, nsmall=2) ,")")) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("N_effect_deg.pdf", device = "pdf", path = "results")
ggsave("N_effect_deg.png", device = "png", path = "results", dpi=300)

#Find genes with nitrogen and time effect
dds_NplusTime_effect <- DESeq(dds_T5, test="LRT", reduced = ~ time)
deg_NplusTime_effect <- subset(results(dds_NplusTime_effect, alpha=0.05), padj<0.05)
head(deg_NplusTime_effect[order(deg_NplusTime_effect$padj),])
gene_id <- "TRIAE_CS42_7DL_TGACv1_603261_AA1979450"
deg <- plotCounts(dds_NplusTime_effect, gene_id, intgroup = c("N","time"), returnData = TRUE)
ggplot(deg, aes(x = time, y = count, color = N, group = N)) +
  geom_point() + geom_smooth(se = FALSE, method = "loess") + scale_y_log10() +
  ggtitle(paste(gene_id, " (padj = ",format(deg_NplusTime_effect[rownames(deg_NplusTime_effect)==gene_id,]$padj, nsmall=2) ,")")) +
  theme(plot.title = element_text(hjust = 0.5))

gene_id <- "TRIAE_CS42_1AL_TGACv1_000064_AA0001950"
deg <- plotCounts(dds_NplusTime_effect, gene_id, intgroup = c("N","time"), returnData = TRUE)
ggplot(deg, aes(x = time, y = count, color = N, group = N)) +
  geom_point() + geom_smooth(se = FALSE, method = "loess") + scale_y_log10() +
  ggtitle(paste(gene_id, " (padj = ",format(deg_NplusTime_effect[rownames(deg_NplusTime_effect)==gene_id,]$padj, nsmall=2) ,")")) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("NplusTime_effect_deg.pdf", device = "pdf", path = "results")
ggsave("NplusTime_effect_deg.png", device = "png", path = "results", dpi=300)

#Find genes differentially expressed at senescence (T5) but that may also be differentially expressed in another time point
dds_T5_bothN <- DESeq(dds_N_effect)
deg_T5 <- subset(results(dds_T5_bothN, name="timeT5", alpha=0.05), padj<0.05)
head(deg_T5[order(deg_T5$padj),])

gene_id <- "TRIAE_CS42_2DS_TGACv1_179582_AA0608070"
deg <- plotCounts(dds_T5_bothN, gene_id, intgroup = c("N","time"), returnData = TRUE)
ggplot(deg, aes(x = time, y = count, color = N, group = N)) +
  geom_point() + geom_smooth(se = FALSE, method = "loess") + scale_y_log10() +
  ggtitle(paste(gene_id, " (padj = ",format(deg_T5[rownames(deg_T5)==gene_id,]$padj, nsmall=2) ,")")) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("deg_T5_overexpressed.pdf", device = "pdf", path = "results")
ggsave("deg_T5_overexpressed.png", device = "png", path = "results", dpi=300)


#Find genes with shifted expression pattern (N and time interaction)
dds_interaction <- DESeq(dds_interaction, test="LRT", reduced = ~ N + time)
deg_interaction <- subset(results(dds_interaction, alpha=0.05), padj<0.05)
head(deg_interaction[order(deg_interaction$padj),])

gene_id <- "TRIAE_CS42_7AL_TGACv1_556210_AA1758400"
deg <- plotCounts(dds_interaction, gene_id, intgroup = c("N","time"), returnData = TRUE)
ggplot(deg, aes(x = time, y = count, color = N, group = N)) +
  geom_point() + geom_smooth(se = FALSE, method = "loess") + scale_y_log10() +
  ggtitle(paste(gene_id, " (padj = ",format(deg_interaction[rownames(deg_interaction)==gene_id,]$padj, nsmall=2) ,")")) +
  theme(plot.title = element_text(hjust = 0.5))

#Find genes activating/deactivating at T5
dds_T5 <- DESeq(dds_T5)

deg100T2 <- subset(results(dds_T5, name="N100N.timeT2", alpha=0.05), padj<0.05)
deg100T3 <- subset(results(dds_T5, name="N100N.timeT3", alpha=0.05), padj<0.05)
deg100T4 <- subset(results(dds_T5, name="N100N.timeT4", alpha=0.05), padj<0.05)
deg100T5 <- subset(results(dds_T5, name="N100N.timeT5", alpha=0.05), padj<0.05) #Senescence
deg200T2 <- subset(results(dds_T5, name="N200N.timeT2", alpha=0.05), padj<0.05)
deg200T3 <- subset(results(dds_T5, name="N200N.timeT3", alpha=0.05), padj<0.05)
deg200T4 <- subset(results(dds_T5, name="N200N.timeT4", alpha=0.05), padj<0.05)
deg200T5 <- subset(results(dds_T5, name="N200N.timeT5", alpha=0.05), padj<0.05) #Senescence

deg100NT5only <- setdiff(rownames(deg100T5), union(rownames(deg100T4), union(rownames(deg100T2), rownames(deg100T3))))
deg200NT5only <- setdiff(rownames(deg200T5), union(rownames(deg200T4), union(rownames(deg200T2), rownames(deg200T3))))
deg_bothN_T5only <- intersect(deg100NT5only, deg200NT5only)

deg100T5vsT2 <- subset(results(dds_T5, contrast=c(0,0,-1,0,0,0,0,0,1,0,0,0,0,0,0,0), alpha=0.05), padj<0.05)
deg100T5vsT3 <- subset(results(dds_T5, contrast=c(0,0,0,0,-1,0,0,0,1,0,0,0,0,0,0,0), alpha=0.05), padj<0.05)
deg100T5vsT4 <- subset(results(dds_T5, contrast=c(0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0), alpha=0.05), padj<0.05)
deg200T5vsT2 <- subset(results(dds_T5, contrast=c(0,0,0,-1,0,0,0,0,0,1,0,0,0,0,0,0), alpha=0.05), padj<0.05)
deg200T5vsT3 <- subset(results(dds_T5, contrast=c(0,0,0,0,0,-1,0,0,0,1,0,0,0,0,0,0), alpha=0.05), padj<0.05)
deg200T5vsT4 <- subset(results(dds_T5, contrast=c(0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0), alpha=0.05), padj<0.05)

deg100NT5only_ <- intersect(rownames(deg100T5), intersect(rownames(deg100T5vsT4), intersect(rownames(deg100T5vsT2), rownames(deg100T5vsT3))))
deg200NT5only_ <- intersect(rownames(deg200T5), intersect(rownames(deg200T5vsT4), intersect(rownames(deg200T5vsT2), rownames(deg200T5vsT3))))
deg_bothN_T5only_ <- intersect(deg100NT5only_, deg200NT5only_)
changeT5 <- intersect(deg_bothN_T5only, deg_bothN_T5only_)

gene_id <- "TRIAE_CS42_7DS_TGACv1_622658_AA2043420"
deg <- plotCounts(dds, gene_id, intgroup = c("N","time"), returnData = TRUE)
ggplot(deg, aes(x = time, y = count, color = N, group = N)) +
geom_point() + geom_smooth(se = FALSE, method = "loess") + scale_y_log10() +
  ggtitle(paste(gene_id, " (padj = ",degTCT5[rownames(degTCT5)==gene_id,]$padj ,")")) +
  theme(plot.title = element_text(hjust = 0.5))

rld <- rlog(dds_T5)
#plotPCA(rld, intgroup = "condition")
#Heatmap of the genes activated/deactivated at T5
changeT5_mat <- assay(rld)[which(rownames(rld)%in% changeT5),]
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
colnames(changeT5_mat) <- rownames(sampleDistMatrix)
#Heatmap of the raw values to get distances of the genes based on absolute count values
pheatmap(changeT5_mat, scale="none", main="Clustering of genes activating/deactivating at T5 (no scaling)")
#Heatmap of the scaled values to get distances of the genes based on their expression profile/deviation from mean (similar values should be obtained using correlation as distance measure)
pheatmap(changeT5_mat, scale="row", main="Clustering of genes activating/deactivating at T5 (with scaling)")



degs <- subset(res, padj < 0.05)
dim(degs)
degs <- res[order(res$padj),]
write.table(degs, file="paragon_degs.tab", quote = F, sep = "\t", row.names = T)

