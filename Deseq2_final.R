#DeSeq2 Installation
#BiocManager::install("DESeq2")
#BiocManager::install("vsn")
#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")
#BiocManager::install("AnnotationHub")
#devtools::install_github("teunbrand/ggh4x")
#BiocManager::install("ReactomePA")

args = commandArgs(trailingOnly=TRUE)

WD <- getwd()
if (!is.null(WD)) setwd(WD)

library("DESeq2")
library("vsn")
library("ggplot2")
library("dplyr")
library("tidyverse")
library("clusterProfiler")
library("biomaRt")
library("pathview")
library("pheatmap")
library("ggh4x")
library("AnnotationHub")
library("ggnewscale")
library("ReactomePA")
library(DOSE)
library(biomaRt)
library(openxlsx)



####### Create Enrich output Directory ##########################
dir.create("3_Analysis/2_Deseq2/Enrich_results", showWarnings = T)
dir.create("3_Analysis/2_Deseq2/Plots/Kegg_Pathway", showWarnings = T)

##################################################################
####################LOADING COUNT DATA#############################
counts <- read.table(args[1], sep = "\t", header = T)
#counts <- read.table("3_Analysis/2_Deseq2/counts.txt", sep = "\t", header = T)

counts <- counts[,-c(2:6)]
row.names(counts) <- counts$Geneid
counts <- counts[-1]
samplenames <- colnames(counts)
###################################################################

############## Extract sample names  ############################
samplenames <- as.data.frame(strsplit(samplenames, split = "[.]"))
colnames(samplenames) <- NULL
samplenames <- as.character(samplenames[5,])
colnames(counts) <- samplenames
##############################################################

########Combination Extraction##################
metatables <- list.files(path = "2_Combinations/")




######################################################## NOW IT's TIME FOR A LOOP ####################################################################
for (i in metatables) {
  path <- paste("2_Combinations",print(i), sep = "/")
  metatable <- read.delim(path, header=T, sep="\t")
 #metatable <- read.table("2_Combinations/control_Sample1", header = T)
  group_name <- strsplit(i, split = "_")
  counts_com <- counts[,match(metatable[,1], colnames(counts))]
  rownames(metatable) <- metatable[,1]
  nrow(counts_com)
  #Verification-1
  verification_1 <- as.character(all(rownames(metatable) %in% colnames(counts_com)))
  verify_1 <- sprintf("%s is %s for case 1", i,verification_1)
  #Verification-2
  verification_2 <- as.character(all(rownames(metatable) == colnames(counts_com)))
  verify_2 <- sprintf("%s is %s for case 2", i,verification_2)
  print(verify_1)
  print(verify_2)
  
  #construction of a DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = counts_com,
                                colData = metatable,
                                design = ~ condition)
  dds
  
  #Filtering the readcounts
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds
  
  #set the factor level for the ctrl vs trtd comparision
  #dds$condition <- relevel(dds$condition, ref = "Control")
  dds$condition <- relevel(dds$condition, ref = metatable[1,2])
  
  
  #Run DESeq
  dds <- DESeq(dds)
  
  #Data transformations and visualization
  #rlog transformation
  rld <- rlog(dds, blind=FALSE)
  #jpeg(paste(i,"3_Analysis/2_Deseq2/rplot_PCA.jpg",sep = "_"), width = 350, height = 350)
  #pca <- plotPCA(rld, intgroup=c("condition", "sample_name"))
  #dev.off()
  
  #plot dispersion estimate
  #par(mfrow=c(1,1))
  #plotDispEsts(dds)
  
  #get results of Deseq2
  res <- results(dds)
  res
  summary(res)
  
  ##to obtain Basemean value for two groups
  a <- sapply( levels(dds$condition), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$condition == lvl] ) )
  #colnames(a) <- c("BaseMean_Control", "BaseMean_Treated")
  
  res_2 <- as.data.frame(res)
  b <- merge(res_2,a, by=0, all=TRUE)
  #b <- na.omit(as.data.frame(b))
  
  ############################################ DGE TABLE GENERATION ##################################################
  # add a column of NAs
  b$Regulation <- "NO"
  # if log2Foldchange > 0.0 and pvalue < 0.05, set as "UP" 
  b$Regulation[b$log2FoldChange > 0.0] <- "UP"
  # if log2Foldchange < -0.0 and pvalue < 0.05, set as "DOWN"
  b$Regulation[b$log2FoldChange < -0.0] <- "DOWN"
  b$Significance <- "No" 
  b$Significance[b$pvalue < 0.05] <- "Yes"
  colnames(b)[1]  <- "Gene_IDs" 
  b <- b[,-c(2,4,5)]
  
  ######################################## CREATE MULTIPLE EXCEL SHEETS ###############################################
  output_path_DGE <- file.path("3_Analysis", "2_Deseq2", "DGE_files", paste(i,"DGE.csv",sep = "_"))
  
  write.csv(b, file = output_path_DGE , row.names = F)
  ###################################################################################################################
  
  ##################################### PLOT TITLE and AXIS NAMES ####################################################
  ctrl <- as.data.frame(group_name)[1,]
  trtd <- as.data.frame(group_name)[2,]
  title <- paste(ctrl, "vs", trtd, sep = " ")
  # add a column of NAs
  b$Expression <- "NO"
  # if log2Foldchange > 0.0 and pvalue < 0.05, set as "UP" 
  b$Expression[b$log2FoldChange > 0.0 & b$pvalue < 0.05] <- "UP"
  # if log2Foldchange < -0.0 and pvalue < 0.05, set as "DOWN"
  b$Expression[b$log2FoldChange < -0.0 & b$pvalue < 0.05] <- "DOWN"
  #####################################################################################################################
  mycolors <- c("red", "limegreen", "black")
  names(mycolors) <- c("UP", "DOWN", "NO")
  
  ############################################### VOLCANO PLOT ########################################################
  ggplot(data=b, aes(x=log2FoldChange, y=-log10(pvalue), col=Expression)) + geom_point(size = 4) + theme_minimal(base_size = 16) + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red") + scale_colour_manual(values = mycolors) + ggtitle(title) +theme(plot.title = element_text(hjust = 0.5))
  
  output_path_PLot <- file.path("3_Analysis", "2_Deseq2", "Plots")
  
  ggsave(path = output_path_PLot, filename = paste(ctrl, trtd, "Volcanoplot.png", sep = "_"), device='png', dpi=1000)
  ######################################################################################################################
  
  ############################################### SCATTERED PLOT ##########################################################
  ggplot(b, aes(x= log10(b[,5]), y=log10(b[,6]), col=Expression)) + theme_minimal(base_size = 16) + geom_point(size = 4) + 
    scale_colour_manual(values = mycolors) + ggtitle(title) + xlab(trtd) + ylab(ctrl) + theme(plot.title = element_text(hjust = 0.5))  
  ggsave(path = output_path_PLot, filename = paste(ctrl, trtd, "Scatteredplot.png", sep = "_"), device='png', dpi=1000)
  #######################################################################################################################
  
  
  ###################################################### HEATMAP ########################################################
  library("pheatmap")
  #HeatMap of normalized count
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:50]
 
  heatmap_ann <- as.data.frame(colData(dds)[,c(1,2)])
  output_path_heatmapPLot <- file.path("3_Analysis", "2_Deseq2", "Plots/", ctrl)
  ##Sample Based heatmap
  pheatmap(assay(rld)[select,], cluster_rows=F, color = colorRampPalette(c("green1", "white", "red"))(20), show_rownames=T,
           cluster_cols=F, annotation_col=heatmap_ann, filename = paste(output_path_heatmapPLot, trtd, "Heatmap.png", sep = "_"))
  
  ##Extract significant ones
  sig <- as.data.frame(na.omit(b))
  sig
  sig <- sig[sig$pvalue < 0.05,]
  summary(sig)
  sum(sig$pvalue < 0.05, na.rm=TRUE) 
  
  # #### get differential expressed gene matrix 
  # #logFC cutoff
  # lfc.cutoff <- 0.0
  # df <- as.data.frame(sig)
  # #df.top <- df[ (df$baseMean > 50) & (abs(df$log2FoldChange) > lfc.cutoff),]
  # df.top <- df[(abs(df$log2FoldChange) > lfc.cutoff),]
  # df.top
  # df.top <- df.top[order(df.top$log2FoldChange, decreasing = T),]
  # head(df.top)
  # dim(df.top)
  # #Filtered rlog based on significant genes HeatMap
  # fil_genes <- a[rownames(df.top),]
  # select_2 <- order(rowMeans(fil_genes),
  #                 decreasing=TRUE)[1:50]  
  
  
  ############################################ Values for Heatmap2 ###############################################################
  sig_up <- sig[sig$Expression == "UP",]
  sig_up <- sig_up[order(sig_up$log2FoldChange, decreasing = T),]
  
  sig_down <- sig[sig$Expression == "DOWN",]
  sig_down <- sig_down[order(sig_down$log2FoldChange, decreasing = F),]
  
  sig_up_down_top50 <-rbind(sig_up[1:25,], sig_down[1:25,])
  rownames(sig_up_down_top50) <- sig_up_down_top50$Gene_IDs
  
  #################### Add increment to the Zero values ##################################################################### 
  sig_up_down_top50[5][sig_up_down_top50[5] == "0"] <- 0+1
  
  sig_up_down_top50[6][sig_up_down_top50[6] == "0"] <- 0+1
  
  heatmap_ann <- as.data.frame(colData(dds)[,c(1,2)])
  #HeatMap of normalized count
  output_path_heatmapPLot <- file.path("3_Analysis", "2_Deseq2", "Plots/", ctrl)
  ##Sample Based heatmap
  pheatmap(log10(sig_up_down_top50[5:6]), cluster_rows=T, color = colorRampPalette(c("green1", "black", "red"))(10), show_rownames=T,
           cluster_cols=T, border_color = "black", cellwidth = 100, cellheight = 10, treeheight_row = 200, filename = paste(output_path_heatmapPLot, trtd, "groupbased_Heatmap.png", sep = "_"))
  
  ##Basemean Heatmap
 
  #pheatmap(log2(fil_genes[1:50,]), cluster_rows=F, color = colorRampPalette(c("green1", "white", "red"))(10), show_rownames=T,
           #cluster_cols=F)
  #######################################################################################################################
  
  if ( args[3] == "NA") {
    print("Completed with STAR pipeline")
  } else {
    
  ################################################# OrgDb Assign ########################################################
  ah <- AnnotationHub()
  orgs <- subset(ah, ah$rdataclass == "OrgDb")
  orgdb <- query(orgs, args[2])[[1]]
  #orgdb <- query(orgs, "Homo sapiens")[[1]]
  
  
  
  ########################### Fetch the genes for the GO and KEGG analysis ##############################################
  Gene_List_ORA <- b$Gene_IDs
  ########################################################################################################################
  
  
  
  ############################################## GO Enrichment Analysis #################################################
  ######################### GO Analysis ORA ########################
  
  ORA_GO_rich <-  enrichGO(gene         = Gene_List_ORA,
                           OrgDb        = orgdb,
                           keyType      = 'ENSEMBL',
                           ont          = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           readable      = TRUE)
  
  output_path_GO_ORA <- file.path("3_Analysis", "2_Deseq2", "Enrich_results", paste(i,"GO_ORA_enrich.tsv",sep = "_"))
  
  head(ORA_GO_rich)
  write.table(ORA_GO_rich, file = output_path_GO_ORA , row.names = F, sep = "\t")
  
  mycolors_GO <- c("red", "green", "blue")
  names(mycolors_GO) <- c("BP", "CC", "MF")
  
  ORA_GO_rich_result <- ORA_GO_rich@result
  ORA_GO_rich_result$Description <- factor(ORA_GO_rich_result$Description, levels = ORA_GO_rich_result$Description)
  
  ###################################### GO Filtration ##################################################
  ORA_GO_rich_result <-  filter(ORA_GO_rich_result, pvalue < 0.05, Count > 10)
  
  ############################# GO Visualization ##############################################
  ORA_GO_BP <-  ORA_GO_rich_result[ORA_GO_rich_result$ONTOLOGY == "BP",]
  ORA_GO_CC <- ORA_GO_rich_result[ORA_GO_rich_result$ONTOLOGY == "CC",]
  ORA_GO_MF <- ORA_GO_rich_result[ORA_GO_rich_result$ONTOLOGY == "MF",]
  
  ORA_GO_rich_result <- rbind(ORA_GO_BP[1:25,], ORA_GO_CC[1:25,], ORA_GO_MF[1:25,])
  ORA_GO_rich_result <- na.omit(ORA_GO_rich_result)
  ##############################################################################################
  
  ####################################### GO ORA Plot ##################################################################
  ggplot(ORA_GO_rich_result, aes(x=-log10(p.adjust), y=Description, fill=ONTOLOGY)) + 
    geom_bar(stat = "identity", col="black") + scale_fill_manual(values = mycolors_GO) + ggtitle(title) +
    facet_grid2(
      ONTOLOGY ~ ., scales="free_y", space="free_y",
      strip = strip_themed(
        background_y = list(element_rect(fill = "red"),
                            element_rect(fill = "green"),
                            element_rect(fill = "blue"))
      )
    )
  ggsave(path = output_path_PLot, filename = paste(ctrl, trtd, "GO_ORA_plot.png", sep = "_"), device='png', dpi=1000, width = 10, height = 10)
  #########################################################################################################################
  #dotplot(ORA_GO_rich, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
  
  ####################################################################
  
  ######################### GO Analysis GSEA ########################
  
  ## To Fetch FC and Gene ID to prepare for GSEA
  ## feature 1: Fetch Fold changes
  geneList_gsea <- b[,2]
  
  ## feature 2: Assign gene names
  names(geneList_gsea) = as.character(b[,1])
  
  ## feature 3: decreasing order
  geneList_gsea = sort(geneList_gsea, decreasing = TRUE)
  gseGO <- gseGO(geneList     = geneList_gsea,
                 OrgDb        = orgdb,
                 ont          = "ALL",
                 keyType = "ENSEMBL",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 eps = 0)
  
  gseGO <- setReadable(gseGO, OrgDb = orgdb, keyType="ENSEMBL")
  
  output_path_GO_gsea <- file.path("3_Analysis", "2_Deseq2", "Enrich_results", paste(i,"GO_GSEA_enrich.tsv",sep = "_"))
  
  write.table(gseGO, output_path_GO_gsea, row.names = F, sep = "\t")
  ####################################################################
  ##############################################################################################################################
  if (nrow(gseGO@result) != 0){
  ################################################# GO Gene Set Enrichment Analysis ######################################
  require(DOSE)
  Go_Gsea_path <- file.path("3_Analysis", "2_Deseq2", "Plots", paste(ctrl, trtd, "GO_GSEA.png", sep = "_"))
  GSEA_Dot_GO <- dotplot(gseGO, showCategory=10, split=".sign") + facet_grid(.~.sign)
  ggplot2::ggsave(file.path(Go_Gsea_path), plot = GSEA_Dot_GO, device='png', dpi=1000, width = 10, height = 10)
  ########################################################################################################################
  }
  else { 
    print("There was no Enrichment in GO ORA")
  }
  # ###################### PLOT-2 ##########################
  # ## count the gene number
  # gseGO_count <- gseGO@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
  # 
  # ## merge with the original dataframe
  # gseGO_dot_df<- left_join(gseGO@result, gseGO_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
  # gseGO_dot_df = gseGO_dot_df[1:50,] ## small dataset
  # gseGO_dot_df$type = "Activated"
  # gseGO_dot_df$type[gseGO_dot_df$NES < 0] = "Repressed"
  # 
  # ### GGPLOT ###########
  # ggplot(gseGO_dot_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
  #   geom_point(aes(size = GeneRatio, color = p.adjust)) +
  #   theme_bw(base_size = 14) +
  #   scale_colour_gradient(limits=c(0, 0.05), low="red") +
  #   ylab(NULL) +
  #   ggtitle(title) + facet_grid(.~type)
  # ggsave(path = output_path_PLot, filename = paste(ctrl, trtd, "GO_GSEA_plot", sep = "_"), device='png', dpi=1000, width = 10, height = 10)
  
  ############################################## KEGG ENRICHMENT ANALYSIS ######################################################
  ######################### KEGG Analysis ORA ########################
  #mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  # Kegg_genes <- getBM(filters="ensembl_gene_id",
  #                     attributes=c("ensembl_gene_id", "entrezgene_id"),
  #                     values=b$Gene_IDs,
  #                     mart=mart)
  
  ################### ENSEMBL id to ENTREZ id CONVERSION #####################
  Kegg_genes <- select(orgdb, keys=b$Gene_IDs, columns="ENTREZID", keytype="ENSEMBL")
  nrow(Kegg_genes)
  Kegg_genes <- na.omit(Kegg_genes)
  
  
  kegg_ORA <- enrichKEGG(gene   = Kegg_genes$ENTREZID,
                         organism = args[3],
                         keyType = "kegg",
                         pvalueCutoff = 0.05)
  kegg_ORA <- setReadable(kegg_ORA, OrgDb = orgdb, keyType="ENTREZID")
  
  output_path_KEGG_ORA <- file.path("3_Analysis", "2_Deseq2", "Enrich_results", paste(i,"KEGG_ORA_enrich.tsv",sep = "_"))
  
  write.table(kegg_ORA, output_path_KEGG_ORA, row.names = F, sep = "\t")
  #############################################################################
  
  
  ####################### KEGG Analysis GSEA ##########################
  ##### Fetch gene and FC for Kegg GSEA
  kegg_gene_list <- b[1:2]
  ##Convert to entrez ID
  kegg_gene_list$EntrezID <- Kegg_genes$ENTREZID[match(kegg_gene_list$Gene_IDs, Kegg_genes$ENSEMBL)]
  nrow(kegg_gene_list)
  #### Remove NA IDs
  kegg_gene_list <- na.omit(kegg_gene_list)
  nrow(kegg_gene_list)
  ## To Fetch FC and Gene ID to prepare for GSEA
  ## feature 1: Fetch Fold changes
  kegg_gene_gsea = kegg_gene_list[,2]
  
  ## feature 2: Assign gene names
  names(kegg_gene_gsea) = as.character(kegg_gene_list[,3])
  
  ## feature 3: decreasing order
  kegg_gene_gsea = sort(kegg_gene_gsea, decreasing = TRUE)
  
  kegg_GSEA <- gseKEGG(geneList = kegg_gene_gsea,
                      organism = args[3],
                      keyType = "kegg",
                      pvalueCutoff = 0.05,
                      eps = 0,
                      verbose = T)
  kegg_GSEA <- setReadable(kegg_GSEA, OrgDb = orgdb, keyType="ENTREZID")
  
  output_path_KEGG_GSEA <- file.path("3_Analysis", "2_Deseq2", "Enrich_results", paste(i,"KEGG_GSEA_enrich.tsv",sep = "_"))
  
  write.table(kegg_GSEA, output_path_KEGG_GSEA, row.names = F, sep = "\t")
  
  ####################################################################################################################################
  if (nrow(kegg_ORA) != 0){
  ################################################### KEGG ORA BAR PLOT #####################################
  Kegg_ORA_Bar_path <- file.path("3_Analysis", "2_Deseq2", "Plots", paste(ctrl, trtd, "Kegg_ORA_Barplot.png", sep = "_"))
  ORA_bar_Kegg <- barplot(kegg_ORA, 
          drop = TRUE, 
          showCategory =30, 
          title = paste(ctrl, trtd, sep = " vs "),
          font.size = 10)
  ggplot2::ggsave(file.path(Kegg_ORA_Bar_path), plot = ORA_bar_Kegg, device='png', dpi=1000, width = 10, height = 10)
  
  ##########################################################################################################
  ################################################# KEGG ORA DOT PLOT ######################################
  Kegg_ORA_dot_path <- file.path("3_Analysis", "2_Deseq2", "Plots", paste(ctrl, trtd, "Kegg_ORA_Dotplot.png", sep = "_"))
  ORA_dot_kegg <- dotplot(kegg_ORA,
          showCategory =30,
          title = paste(ctrl, trtd, sep = " vs "))
  ggplot2::ggsave(file.path(Kegg_ORA_dot_path), plot = ORA_dot_kegg, device='png', dpi=1000, width = 10, height = 10)
  ############################################################################################################
  }
  else { 
    print("There was no Enrichment in KEGG ORA")
  }
  
  if (nrow(kegg_GSEA) != 0){
  ################################################### KEGG GSEA VISUALIZATION #########################################################
  
  Kegg_GSEA_dot_path <- file.path("3_Analysis", "2_Deseq2", "Plots", paste(ctrl, trtd, "Kegg_GSEA_Dotplot.png", sep = "_"))
  GSEA_dot_Kegg <- dotplot(kegg_GSEA, showCategory = 10, title = paste(ctrl, trtd, sep = " vs ") , split=".sign") + facet_grid(.~.sign)
  ggplot2::ggsave(file.path(Kegg_GSEA_dot_path), plot = GSEA_dot_Kegg, device='png', dpi=1000, width = 10, height = 10)
  #####################################################################################################################################
  }
  else { 
    print("There was no Enrichment in KEGG GSEA")
  }
  ##################################################### KEGG VISUALIZATION #################################################################
  #output_path_pathview <- file.path("3_Analysis", "2_Deseq2", "Plots", "Kegg_Pathway")
  
  #pathways <- pathview(gene.data  = kegg_gene_gsea,
                       #pathway.id = kegg_GSEA_results$ID[1:10],
                       #species    = "hsa",
                       #kegg.dir = output_path_pathview,
                       #limit      = list(gene=round(max(abs(geneList)), digits = 1), cpd=1))
  #######################################################################
  #########################################################################################################################################
  
  if ( args[3] == "hsa") {
  ##################################################### REACTOME PATHWAY ANALYSIS ########################################################
  
  ################# Convert Ensemble Id to Entrez ID #############################
  Reac_genes <- select(orgdb, keys=b$Gene_IDs, columns="ENTREZID", keytype="ENSEMBL")
  nrow(Reac_genes)
  Reac_genes <- na.omit(Reac_genes)
  
  #########################################  REACTOME ORA Analysis ######################################################
  Reac_ORA <-  enrichPathway(gene=Reac_genes$ENTREZID, organism = "human", pvalueCutoff = 0.05, readable=TRUE)

  output_path_Reac_ORA <- file.path("3_Analysis", "2_Deseq2", "Enrich_results", paste(i,"Reactome_ORA_enrich.tsv",sep = "_"))
  
  write.table(Reac_ORA, output_path_Reac_ORA, row.names = F, sep = "\t")
  #####################################################################################################################
  
  #########################################  REACTOME GSEA Analysis ######################################################
  Reac_gene_list <- b[1:2]
  ##Convert to entrez ID
  Reac_gene_list$EntrezID <- Reac_genes$ENTREZID[match(Reac_gene_list$Gene_IDs, Reac_genes$ENSEMBL)]
  nrow(Reac_gene_list)
  #### Remove NA IDs
  Reac_gene_list <- na.omit(Reac_gene_list)
  nrow(Reac_gene_list)
  ## To Fetch FC and Gene ID to prepare for GSEA
  ## feature 1: Fetch Fold changes
  Reac_gene_gsea = Reac_gene_list[,2]
  
  ## feature 2: Assign gene names
  names(Reac_gene_gsea) = as.character(Reac_gene_list[,3])
  
  ## feature 3: decreasing order
  Reac_gene_gsea = sort(Reac_gene_gsea, decreasing = TRUE)
  
  Reac_GSEA <- gsePathway(geneList = Reac_gene_gsea,
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH", 
                          verbose = T,
                          eps = 0)
  Reac_GSEA <- setReadable(Reac_GSEA, OrgDb = orgdb, keyType="ENTREZID")
  
  output_path_Reac_GSEA <- file.path("3_Analysis", "2_Deseq2", "Enrich_results", paste(i,"Reactome_GSEA_enrich.tsv",sep = "_"))
  
  write.table(Reac_GSEA, output_path_Reac_GSEA, row.names = F, sep = "\t")
  #####################################################################################################################
  
  if (nrow(Reac_ORA) != 0){
  ######################################### REACTOME ORA BAR PLOT #####################################################
  Reac_ORA_Bar_path <- file.path("3_Analysis", "2_Deseq2", "Plots", paste(ctrl, trtd, "Reactome_ORA_barplot.png", sep = "_"))
  ORA_bar_Reac <- barplot(Reac_ORA, 
          drop = TRUE, 
          showCategory =30, 
          title = paste(ctrl, trtd, sep = " vs "),
          font.size = 10)
  ggplot2::ggsave(file.path(Reac_ORA_Bar_path), plot = ORA_bar_Reac, device='png', dpi=1000, width = 10, height = 10)
  
  ##########################################################################################################################
  ################################################# Reactome ORA DOT PLOT ##################################################
  Reac_ORA_dot_path <- file.path("3_Analysis", "2_Deseq2", "Plots", paste(ctrl, trtd, "Reactome_ORA_Dotplot.png", sep = "_"))
  ORA_dot_Reac <- dotplot(Reac_ORA,
                          showCategory = 30,
                          title = paste(ctrl, trtd, sep = " vs "),
                          font.size = 10)
  ggplot2::ggsave(file.path(Reac_ORA_dot_path), plot = ORA_dot_Reac, device='png', dpi=1000, width = 10, height = 10)
  ###########################################################################################################################
  }
  else { 
    print("There was no Enrichment in REACTOME ORA")
  }
  
  if (nrow(Reac_GSEA) != 0){
  ################################################### Reactome GSEA VISUALIZATION ####################################################
  
  Reac_GSEA_dot_path <- file.path("3_Analysis", "2_Deseq2", "Plots", paste(ctrl, trtd, "Reactome_GSEA_Dotplot.png", sep = "_"))
  GSEA_dot_Reac <- dotplot(Reac_GSEA, showCategory = 10, title = paste(ctrl, trtd, sep = " vs ") , split=".sign") + facet_grid(.~.sign)
  ggplot2::ggsave(file.path(Reac_GSEA_dot_path), plot = GSEA_dot_Reac, device='png', dpi=1000, width = 10, height = 10)
  #####################################################################################################################################
  }
  else { 
    print("There was no Enrichment in REACTOME GSEA")
  }
  
  ################################################## DO  analysis #################################################################
  ################# Convert Ensemble Id to Entrez ID #############################
  DO_genes <- select(orgdb, keys=b$Gene_IDs, columns="ENTREZID", keytype="ENSEMBL")
  nrow(DO_genes)
  DO_genes <- na.omit(DO_genes)
  
  ############################################# DO ORA ANALYSIS ###############################################################
  DO_ORA <- enrichDO(gene = DO_genes$ENTREZID,
           ont           = "DO",
           pvalueCutoff  = 0.05,
           pAdjustMethod = "BH",
           qvalueCutoff  = 0.05,
           readable      = T)
  output_path_DO_ORA <- file.path("3_Analysis", "2_Deseq2", "Enrich_results", paste(i,"Disease_Ontology_ORA_enrich.tsv",sep = "_"))
  
  write.table(DO_ORA, output_path_DO_ORA, row.names = F, sep = "\t")
  #####################################################################################################################
  
  #########################################  Disease Ontology GSEA Analysis ######################################################
  DO_gene_list <- b[1:2]
  ##Convert to entrez ID
  DO_gene_list$EntrezID <- DO_genes$ENTREZID[match(DO_gene_list$Gene_IDs, DO_genes$ENSEMBL)]
  nrow(DO_gene_list)
  #### Remove NA IDs
  DO_gene_list <- na.omit(DO_gene_list)
  nrow(DO_gene_list)
  ## To Fetch FC and Gene ID to prepare for GSEA
  ## feature 1: Fetch Fold changes
  DO_gene_gsea = DO_gene_list[,2]
  
  ## feature 2: Assign gene names
  names(DO_gene_gsea) = as.character(DO_gene_list[,3])
  
  ## feature 3: decreasing order
  DO_gene_gsea = sort(DO_gene_gsea, decreasing = TRUE)
  
  DO_GSEA <- gseDO(DO_gene_gsea,
        pvalueCutoff  = 0.05,
        pAdjustMethod = "BH",
        verbose = T,
        eps = 0)
  
  DO_GSEA <- setReadable(DO_GSEA, OrgDb = orgdb, keyType="ENTREZID")
  
  output_path_DO_GSEA <- file.path("3_Analysis", "2_Deseq2", "Enrich_results", paste(i,"Disease_Ontology_GSEA_enrich.tsv",sep = "_"))
  
  write.table(DO_GSEA, output_path_DO_GSEA, row.names = F, sep = "\t")
  ######################################################################################################################################
  
  if (nrow(DO_ORA) != 0){
  ######################################### Disease Ontology ORA BAR PLOT #####################################################
  DO_ORA_Bar_path <- file.path("3_Analysis", "2_Deseq2", "Plots", paste(ctrl, trtd, "Disease_Ontology_ORA_barplot.png", sep = "_"))
  ORA_bar_DO <- barplot(DO_ORA, 
                          drop = TRUE, 
                          showCategory =30, 
                          title = paste(ctrl, trtd, sep = " vs "),
                          font.size = 10)
  ggplot2::ggsave(file.path(DO_ORA_Bar_path), plot = ORA_bar_DO, device='png', dpi=1000, width = 10, height = 10)
  
  ##########################################################################################################################
  ################################################# Disease Ontology ORA DOT PLOT ##################################################
  DO_ORA_dot_path <- file.path("3_Analysis", "2_Deseq2", "Plots", paste(ctrl, trtd, "Disease_Ontology_ORA_Dotplot.png", sep = "_"))
  ORA_dot_DO <- dotplot(DO_ORA,
                          showCategory = 30,
                          title = paste(ctrl, trtd, sep = " vs "),
                          font.size = 10)
  ggplot2::ggsave(file.path(DO_ORA_dot_path), plot = ORA_dot_DO, device='png', dpi=1000, width = 10, height = 10)
  ###########################################################################################################################
  }
  else { 
    print("There was no Enrichment in Disease Ontology ORA")
  }
  if (nrow(DO_GSEA) != 0){
  ################################################### Disease Ontology GSEA VISUALIZATION ####################################################
  
  DO_GSEA_dot_path <- file.path("3_Analysis", "2_Deseq2", "Plots", paste(ctrl, trtd, "Disease_Ontology_GSEA_Dotplot.png", sep = "_"))
  GSEA_dot_DO <- dotplot(DO_GSEA, showCategory = 10, title = paste(ctrl, trtd, sep = " vs ") , split=".sign") + facet_grid(.~.sign)
  ggplot2::ggsave(file.path(DO_GSEA_dot_path), plot = GSEA_dot_DO, device='png', dpi=1000, width = 10, height = 10)
  #####################################################################################################################################
  }
  else { 
    print("There was no Enrichment in Disease Ontology GSEA")
  }
  
  
  ################################################## Network Cancer Gene  analysis #################################################################
  ################# Convert Ensemble Id to Entrez ID #############################
  NCG_genes <- select(orgdb, keys=b$Gene_IDs, columns="ENTREZID", keytype="ENSEMBL")
  nrow(NCG_genes)
  NCG_genes <- na.omit(NCG_genes)
  
  ############################################# Network of Cancer Gene ORA ANALYSIS ###############################################################
  NCG_ORA <- enrichNCG(gene = NCG_genes$ENTREZID,
                     pvalueCutoff  = 0.05,
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.05,
                     readable      = T)
  output_path_NCG_ORA <- file.path("3_Analysis", "2_Deseq2", "Enrich_results", paste(i,"Network_Cancer_Gene_ORA_enrich.tsv",sep = "_"))
  
  write.table(NCG_ORA, output_path_NCG_ORA, row.names = F, sep = "\t")
  #####################################################################################################################
  
  #########################################  Network Cancer Gene GSEA Analysis ######################################################
  NCG_gene_list <- b[1:2]
  ##Convert to entrez ID
  NCG_gene_list$EntrezID <- NCG_genes$ENTREZID[match(NCG_gene_list$Gene_IDs, NCG_genes$ENSEMBL)]
  nrow(NCG_gene_list)
  #### Remove NA IDs
  NCG_gene_list <- na.omit(NCG_gene_list)
  nrow(NCG_gene_list)
  ## To Fetch FC and Gene ID to prepare for GSEA
  ## feature 1: Fetch Fold changes
  NCG_gene_gsea = NCG_gene_list[,2]
  
  ## feature 2: Assign gene names
  names(NCG_gene_gsea) = as.character(NCG_gene_list[,3])
  
  ## feature 3: decreasing order
  NCG_gene_gsea = sort(NCG_gene_gsea, decreasing = TRUE)
  
  NCG_GSEA <- gseNCG(NCG_gene_gsea,
                   pvalueCutoff  = 0.05,
                   pAdjustMethod = "BH",
                   verbose = T,
                   eps = 0)
  
  NCG_GSEA <- setReadable(NCG_GSEA, OrgDb = orgdb, keyType="ENTREZID")
  
  output_path_NCG_GSEA <- file.path("3_Analysis", "2_Deseq2", "Enrich_results", paste(i,"Network_Cancer_Gene_GSEA_enrich.tsv",sep = "_"))
  
  write.table(NCG_GSEA, output_path_NCG_GSEA, row.names = F, sep = "\t")
  #############################################################################################################################
  
  if (nrow(NCG_ORA) != 0){
  ######################################### Network Cancer Gene ORA BAR PLOT #####################################################
  NCG_ORA_Bar_path <- file.path("3_Analysis", "2_Deseq2", "Plots", paste(ctrl, trtd, "Network_Cancer_Gene_ORA_barplot.png", sep = "_"))
  ORA_bar_NCG <- barplot(NCG_ORA, 
                        drop = TRUE, 
                        showCategory =30, 
                        title = paste(ctrl, trtd, sep = " vs "),
                        font.size = 10)
  ggplot2::ggsave(file.path(NCG_ORA_Bar_path), plot = ORA_bar_NCG, device='png', dpi=1000, width = 10, height = 10)
  
  ##########################################################################################################################
  ################################################# Network Cancer Gene ORA DOT PLOT ##################################################
  NCG_ORA_dot_path <- file.path("3_Analysis", "2_Deseq2", "Plots", paste(ctrl, trtd, "Network_Cancer_Gene_ORA_Dotplot.png", sep = "_"))
  ORA_dot_NCG <- dotplot(NCG_ORA,
                        showCategory = 30,
                        title = paste(ctrl, trtd, sep = " vs "),
                        font.size = 10)
  ggplot2::ggsave(file.path(NCG_ORA_dot_path), plot = ORA_dot_NCG, device='png', dpi=1000, width = 10, height = 10)
  ###########################################################################################################################
  }
  else { 
    print("There was no Enrichment in Network Cancer Gene ORA")
  }
  
  if (nrow(NCG_GSEA) != 0){
  ################################################### Network Cancer Gene GSEA VISUALIZATION ####################################################
  
  NCG_GSEA_dot_path <- file.path("3_Analysis", "2_Deseq2", "Plots", paste(ctrl, trtd, "Network_Cancer_Gene_GSEA_Dotplot.png", sep = "_"))
  GSEA_dot_NCG <- dotplot(NCG_GSEA, showCategory = 10, title = paste(ctrl, trtd, sep = " vs ") , split=".sign") + facet_grid(.~.sign)
  ggplot2::ggsave(file.path(NCG_GSEA_dot_path), plot = GSEA_dot_NCG, device='png', dpi=1000, width = 10, height = 10)
  ###############################################################################################################################################
  }
  else { 
    print("There was no Enrichment in Network Cancer Gene GSEA")
  }
    } else {
    print("Completed with Enrichment Analysis")
    }
  }
}



