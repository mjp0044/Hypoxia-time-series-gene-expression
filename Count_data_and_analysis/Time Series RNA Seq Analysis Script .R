#Hypoxia Time Series RNA Seq Analysis
#Matthew J. Powers
#Last edited 11-06-24
#Built on R version 4.3.3 "Angel Food Cake"
#Also tested on R version 4.4.0 "Puppy Cup"


#Gene expression and GO analysis packages installed using BiocManager::install("packagename")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")


library(remotes) #To install packages from github using 'install_github' function as below
#install_github("HuntsmanCancerInstitute/hciR")

library(hciR) #To utilize the featureCounts read in function, and format to count matrix

#Use to install gene expression and GO term testing packages
#BiocManager::install("ggkegg")

#Gene expression and cluster analysis packages
library(DESeq2) #Load DESeq2
library(edgeR) #Load edgeR
library(maSigPro) #Load maSigPro to analyze expression data over time course
library(mclust) #To use Mclust method with maSigPro

#Batch effect correction
library(sva) #Also loads ComBat-seq for batch correction 

#Go term testing
library(topGO) #GO analysis
library(rrvgo) #Vizualize results from GO enrichment analysis
library(Rgraphviz) #Visualize results from topGO

#Install from local database
#install.packages("./org.Tcalifornicus.eg.db", type = "source", repos=NULL)
library(org.Tcalifornicus.eg.db) #Load custom OrgDb package for T. californicus
Tcalif_orgdb_object <- org.Tcalifornicus.eg.db #Assign package to an object for use in rrvgo

#All other packages with required/preferred functions
library(ggpubr) #Also loads basic ggplot2
library(ggrepel) #For added repel geom in ggplot
library(cowplot) #Pretty ggplots!
library(reshape2) #Data wrangling
library(dplyr) #Data wrangling
library(tidyverse) #Data wrangling
library(stringr) #data wrangling for splitting column characters
library(MetBrewer) #Pretty colors!
library(MoMAColors) #Pretty colors!
library(palettetown) #Pretty ggplots!
library(scales) #Be able to preview pretty colors
library(patchwork) #stitching graphs together
library(ggiraph) #Connected interactive graphs
library(MESS) #data wrangling and modeling
library(stats) #data wrangling and modeling
library(lme4) #data wrangling and modeling
library(sparrow) #Data wrangling (specifically scale_row function) NOTE: the results function in sparrow blocks the
# use of the results function from deseq2. Specify DESeq2::results for pairwise comparisons
library(lmerTest) #Get null hypothesis test from lme4 lmer models
library(emmeans) #data wrangling and modeling
library(metafor) # To use the escalc function to calculate log odds ratio
library(epitools) # To use the oddsratio.fisher function to calculate the odds ratio
library(multcomp) #automatically assign pairwise significance groupings
library(multcompView) #See groupings from multcomp
library(betareg) #To model proportions
library(reactablefmtr) #load reactable and expansion for table export and creation
library(ggVennDiagram) #Make pretty venn diagrams
library(venn) #Make venn diagrams with more than 4 groups
library(pheatmap) #Generate heatmaps 
library(ggkegg) #Generate kegg pathways under ggplot framework
library(tidygraph) #Used with ggkegg
#Create custom '%notin%' function for querying datasets
`%notin%` <- Negate(`%in%`)
#Set theme globally
theme_set(theme_cowplot())



# Summarize read totals at each stage of the pipeline
     pipe_tots <- read.csv(file = "Pipeline Totals.csv", header = TRUE)
     
   # Remove percent symbols to convert to numeric for summary function  
     pipe_tots$Percent.retained.after.trim <-as.numeric(gsub("\\%","", pipe_tots$Percent.retained.after.trim))
     pipe_tots$Percent.mapped <-as.numeric(gsub("\\%","", pipe_tots$Percent.mapped))
     pipe_tots$Percent.counted.mapped <-as.numeric(gsub("\\%","", pipe_tots$Percent.counted.mapped))
   
    #Generate summary stats of each column   
     summary(pipe_tots)

#Read in main count data .txt files from featureCounts using the hciR package and save it as a count matrix
#Read in as a data frame instead of a tibble. 
#All feature counts txt files must be in a folder called "feature_counts" within your R working directory (with no other txt files in that folder!)
    ft_counts <- as.data.frame(read_featureCounts(path = "./feature_counts", pattern = '.txt', stats = FALSE))

#Rename rownames using the geneid in the first column and then delete first column
    rownames(ft_counts) <- ft_counts[,1]
    ft_counts[,1] <- NULL

#Convert counts to numeric for DESeq2
    ft_counts %>% 
      mutate_all(as.integer) -> ft_counts 
    str(ft_counts) #Check coding
    
#Re-order count matrix to be ordered left to right in ascending timepoint on oxygen curve for plotting (in impulse and others)
    ft_counts = ft_counts[,c(19,20,21,22,23,24, 7,8,9,10,11,12, 1,2,3,4,5,6, 13,14,15,16,17,18, 25,26,27,28,29,30)]


#Read in the metadata for the DGE programs, specifying sample ID as rownames
#Code info as factors for DESeq2. Can change time point to numeric later for maSigPro
    coldata <- read.csv(file="timeseries_coldata.csv", header = TRUE, row.names = 1)
    
    coldata %>% 
      mutate_all(as.factor) -> coldata
    str(coldata)#Check coding


##Filtering
    ft_counts.nomito <- ft_counts[-c(1),] #Remove the mitochondrial DNA
    
    #Remove genes where we didn't have at least 6 samples with 10 counts
    #Chose six based on DESeq2 recommendations of matching the smallest number of replicates for a given grouping in the data set
    #Ours is 6 replicates per time point
      ft_counts_filtered <- ft_counts.nomito[rowSums(ft_counts.nomito >= 10) >= 6, ] 
      
      write.csv(ft_counts_filtered, file="raw_counts_no_mito_filtered.csv") #Export raw counts after filtering
    
    #Batch effect correction among replicates using Combat-seq
    
    # Basic usage (users need to input at least two parameters - 
    # a raw count matrix from RNA-Seq studies, without any normalization or transformation, and a vector for batch separation):
      combat_batch <- rep(c(0,1,2,4,5,6), 5) #Specify the replicates as batches
      combat_group <- c(rep(1, 6), rep(2, 6), rep(3, 6), rep(4, 6), rep(5, 6))#Specify the time points as biological variables
      ft_filtered_combat <- ComBat_seq(as.matrix(ft_counts_filtered), batch=combat_batch, group=combat_group)
    #Export raw counts after filtering 
      write.csv(ft_filtered_combat, file="raw_counts_filtered_batch_corrected.csv") 
    
    # Remove intermediate data frames
      rm(ft_counts.nomito)
      rm(ft_counts_filtered)
    
#Sanity check: make sure rownames of the info key match the colnames of the count matrix
    all(rownames(coldata) == colnames(ft_filtered_combat)) #needs to be TRUE


##Normalization using DEseq2
#Create DESeq object
    dds <- DESeqDataSetFromMatrix(countData = ft_filtered_combat,
                                  colData = coldata,
                                  design = ~ Time)
    dds #View object
    
    dds2 <-DESeq(dds) #run the normalization and save as new object
    
    rm(dds) #Get rid of old object
    
    ddsNorm <- counts(dds2, normalized=TRUE) #Get the normalized counts out for maSigPro. 
    
    write.csv(ddsNorm, file = "normalized_counts_from_DESeq2.csv")

    
  # Create scaled by row version of normalized counts for visualizing expressions patterns and 
  # plotting without swamping from high gene counts
    ddsNormScaled <- scale_rows(ddsNorm, center =  TRUE, scale = TRUE)
    
############################ maSigPro ################################

#Get edesign object from coldata above and reassign time as.numeric for maSigPro
    hypoxia.edesign <- coldata
    hypoxia.edesign$Time <- as.numeric(hypoxia.edesign$Time)
    hypoxia.edesign$Replicate <- as.character(hypoxia.edesign$Replicate)
    hypoxia.edesign$Group <- as.character(hypoxia.edesign$Group)
    str(hypoxia.edesign)
    
#Sanity check: make sure rownames of the info key match the colnames of the count matrix
    all(rownames(hypoxia.edesign) == colnames(ddsNorm)) #needs to be TRUE
    
#Create a regression matrix for the full regression model for the single series:
    hypoxia.design <- make.design.matrix(hypoxia.edesign, degree = 3)
#View comparisons
    hypoxia.design$groups.vector
    
# run maSigPro
#compute a regression fit for each gene.
#By default maSigPro corrects this p-value for multiple comparisons by applying the linear step-up (B-H) false discovery rate (FDR) 
#procedure (Benjamini and Hochberg, 1995)
#The level of FDR control is given by the function parameter Q.
  ss.hypoxia <- p.vector(ddsNorm, hypoxia.design, Q = 0.05, MT.adjust = "BH", counts = TRUE)
  ss.hypoxia$i # returns the number of significant genes
  
  #Send genes and pvalues to dataframe. 
  #Below, in synthesis section, will combine with SD.2.2 annotationall.masigpro <- as.data.frame(hypoxia.all$sig.genes$Group$sig.pvalues)
  all.masigpro <- as.data.frame(ss.hypoxia$p.vector)
  all.masigpro$padj_maSigPro <- ss.hypoxia$p.adjusted
  all.masigpro$Gene <- rownames(all.masigpro) #create column with gene names and reset row names
  rownames(all.masigpro) <- NULL
  
  
#Once significant genes have been found, maSigPro applies a variable selection procedure to find significant variables for each gene
  hypoxia.tstep <- T.fit(ss.hypoxia, alfa = 0.05)
  hypoxia.sigs <- get.siggenes(hypoxia.tstep, rsq = 0.1, vars = "groups")
  
  #Send genes and pvalues to dataframe. Below, in synthesis section, will combine with SD.2.2 annotationall.masigpro <- as.data.frame(hypoxia.all$sig.genes$Group$sig.pvalues)
  all.masigpro$Gene <- rownames(all.masigpro) #create column with gene names and reset row names
  rownames(all.masigpro) <- NULL
  
# Cluster analysis in maSigPro
# Use see.genes() to visualize the result of a group of genes
#Cluster.data set to 9 by default
#Can start with k = 1 to group all together to visualize profiles as groups, then try splitting into more. 
  
  #Recode edesign as integers for the see.genes function
    hypoxia.sigs$sig.genes$Group$edesign$Replicate <- as.integer(hypoxia.sigs$sig.genes$Group$edesign$Replicate)
    hypoxia.sigs$sig.genes$Group$edesign$Group <- as.integer(hypoxia.sigs$sig.genes$Group$edesign$Group)
    
    str(hypoxia.sigs$sig.genes$Group$edesign) #check coding
    
  #Cluster genes in k clusters and save to object sigs.clust
    sigs.clust <- see.genes(hypoxia.sigs$sig.genes$Group, show.fit =T, dis = hypoxia.design$dis, 
              cluster.method="hclust" ,cluster.data = 1, distance = "cor")
  
  
  #Output to pdf
    pdf(file='maSigPro hypoxia single series data gene clusters.pdf', width = 9, height = 8)
    see.genes(hypoxia.sigs$sig.genes$Group, show.fit =T, dis = hypoxia.design$dis, 
              cluster.method="hclust" ,cluster.data = 1, distance = "cor")
    dev.off()
  
  
  #Isolate sig genes clusters in a data frame
   sig.gene.clusters <- as.data.frame(sigs.clust$cut)
   sig.gene.clusters$genes <- row.names(sig.gene.clusters)
   names(sig.gene.clusters)[1] <- "clusters"
  #Export clusters
    write.csv(sigs.clust$cut, file = "Clustered_genes_maSigPro.csv")
  
  
  # Make matrix with the significant genes and their expression values
    sigs.masigpro <- as.data.frame(hypoxia.sigs$sig.genes$Group$sig.profiles) 
    
    
# Assign cluster groupings from clusters data frame  
    sigs.masigpro$clusters <- sig.gene.clusters$clusters
    
  #Save gene names to vector
    sigs.masigpro.names <- row.names(sigs.masigpro)
    sigs.masigpro.names.cluster1 <- row.names(sigs.masigpro[sigs.masigpro$clusters == "1",])
    sigs.masigpro.names.cluster2 <- row.names(sigs.masigpro[sigs.masigpro$clusters == "2",])
    sigs.masigpro.names.cluster3 <- row.names(sigs.masigpro[sigs.masigpro$clusters == "3",])
    sigs.masigpro.names.cluster4 <- row.names(sigs.masigpro[sigs.masigpro$clusters == "4",])
    sigs.masigpro.names.cluster5 <- row.names(sigs.masigpro[sigs.masigpro$clusters == "5",])
    sigs.masigpro.names.cluster6 <- row.names(sigs.masigpro[sigs.masigpro$clusters == "6",])
    sigs.masigpro.names.cluster7 <- row.names(sigs.masigpro[sigs.masigpro$clusters == "7",])
    sigs.masigpro.names.cluster8 <- row.names(sigs.masigpro[sigs.masigpro$clusters == "8",])
    sigs.masigpro.names.cluster9 <- row.names(sigs.masigpro[sigs.masigpro$clusters == "9",])
    
  # Add gene name to its own column for indexing later. Will ignore it in melting step below. 
    sigs.masigpro$Gene <- row.names(sigs.masigpro)
    
    
          # # Combine clusters 6 and 7
          #   #Recode clusters as factor
          #     sigs.masigpro$clusters <- as.factor(sigs.masigpro$clusters)
          #   # Rename them to be the same
          #     sigs.masigpro %>% 
          #      mutate(clusters = fct_recode(clusters,"67" = "6",
          #                                  "67" = "7")) -> sigs.masigpro
          # #Save gene names to vector for combined clusters 6 and 7  
          #     sigs.masigpro.names.cluster67 <- row.names(sigs.masigpro[sigs.masigpro$clusters == "67",])
  
  
  #Send coefficients and p-values to a single data frame in case we need it for later
    sigs.masigpro.pvalues.and.betas <- merge(hypoxia.sigs$sig.genes$Group$sig.pvalues, hypoxia.sigs$sig.genes$Group$coefficients, by = 'row.names')
    names(sigs.masigpro.pvalues.and.betas)[1] <- "Gene" #Rename geneId columnn
    
  # Assign cluster groupings from clusters data frame  
    sigs.masigpro.pvalues.and.betas$clusters <- sig.gene.clusters$clusters
    
              # # Combine clusters 6 and 7
              #   #Recode clusters as factor
              #   sigs.masigpro.pvalues.and.betas$clusters <- as.factor(sigs.masigpro.pvalues.and.betas$clusters)
              #   # Rename them to be the same
              #   sigs.masigpro.pvalues.and.betas %>% 
              #     mutate(clusters = fct_recode(clusters,"67" = "6",
              #                                  "67" = "7")) -> sigs.masigpro.pvalues.and.betas
  
  #Write the sig genes expression profiles to a csv
    write.csv(sigs.masigpro, file = "Sig_genes_maSigPro_expression_profiles.csv", row.names = FALSE)
  #Write the sig genes pvalues and coefficients to a csv
    write.csv(sigs.masigpro.pvalues.and.betas, file = "Sig_genes_maSigPro_pvalues_and_coefficients.csv", row.names = FALSE)
  
  #Make melted data frame for plotting
    sigs.masigpro.melt <- sigs.masigpro[, -c(32)] #Leave off last column
    sigs.masigpro.melt$genes <- row.names(sigs.masigpro.melt) #Send gene names to a column to use as id for melting
    rownames(sigs.masigpro.melt) <- NULL #Make row.names default numbers
    sigs.masigpro.melt <- melt(sigs.masigpro.melt, id = c("genes", "clusters")) #Melt the data frame lengthwise
    names(sigs.masigpro.melt)[3:4] <- c("sampleID", "expression") #Rename the melted columns
    
    
    sigs.masigpro.melt[c('time', 'replicate')] <- str_split_fixed(sigs.masigpro.melt$sampleID, '_', 2)
    #Relevel group variable to put in the desired order for plotting
    sigs.masigpro.melt$time<- factor(sigs.masigpro.melt$time,
                                     levels = c("C", "35","05", "A", "R"))
  
    sigs.masigpro.melt$gene.rep <- paste(sigs.masigpro.melt$genes,sigs.masigpro.melt$replicate)
  

  
#Replicate summary
  sigs.averaged <- sigs.masigpro.melt |>
    group_by(time, clusters) |> 
    summarise(mean_exp = mean(expression) |> round(4),
              sd_exp = sd(expression) |> round(4),
              n = n(),
              range = paste(range(expression) |> round(4), collapse = ' - ')
    ) |> 
    mutate(se_exp = sd_exp/sqrt(n)) 
  sigs.averaged$se_exp <- round(sigs.averaged$se_exp, 4) #Round se column to 2 digits
  
# Recode time points for plotting more clearly
  #Rename 35 and 05 to H and P for hypoxia and pcrit respectively
    sigs.averaged %>% 
      mutate(time = fct_recode(time,"Mild-hypox" = "35",
                              "Pcrit" = "05",
                              "Normoxia" = "C",
                              "Anoxia" = "A",
                              "Recovery" = "R")) -> sigs.averaged
    
  #Write to csv
    write.csv(sigs.averaged, file = "Data average example maSigPro.csv", row.names = FALSE)


#Plot expression levels of significant genes across time points for each of five clusters 
      expression.curve.avgs.clust.1 <- sigs.averaged |>
        filter(clusters == 1) |>
        ggplot(aes(x=time, y=mean_exp)) +
        geom_line(aes(group=clusters), linewidth=2, alpha=1)+
        geom_point(size = 5)+
        geom_errorbar(aes(ymax = mean_exp + se_exp, ymin = mean_exp - se_exp), width=0.3) +
        #geom_label_repel(data = subset(sigs.averaged, time == "A"), aes(label=genes), nudge_x = 2, max.overlaps = 20, size = 3)+
        scale_y_continuous(name="Normalized Expression", breaks = seq(300, 550, by = 50), limits = c(300, 550))+
        scale_x_discrete(name="")+
        theme(axis.text.x = element_text(size = 9),
              axis.title = element_text(size = 13),
              axis.text.y = element_text(size=10),
              legend.position = "bottom", legend.direction = "horizontal")+
        annotate("text", 2, 550, label = "n = 182 genes", size = 4)+
        ggtitle("Cluster 1")
      
      expression.curve.avgs.clust.1
      
      expression.curve.avgs.clust.2 <- sigs.averaged |>
        filter(clusters == 2) |>
        ggplot(aes(x=time, y=mean_exp)) +
        geom_line(aes(group=clusters), linewidth=2, alpha=1)+
        geom_point(size = 5)+
        geom_errorbar(aes(ymax = mean_exp + se_exp, ymin = mean_exp - se_exp), width=0.3) +
        #geom_label_repel(data = subset(sigs.averaged, time == "A"), aes(label=genes), nudge_x = 2, max.overlaps = 20, size = 3)+
        scale_y_continuous(name="Normalized Expression", breaks = seq(600, 1400, by = 100), limits = c(600, 1450))+
        scale_x_discrete(name="")+
        theme(axis.text.x = element_text(size = 9),
              axis.title = element_text(size = 13),
              axis.text.y = element_text(size=10),
              legend.position = "bottom", legend.direction = "horizontal")+
        annotate("text", 2, 1400, label = "n = 264 genes", size = 4)+
        ggtitle("Cluster 2")
      
      expression.curve.avgs.clust.2
      
      expression.curve.avgs.clust.3 <- sigs.averaged |>
        filter(clusters == 3) |>
        ggplot(aes(x=time, y=mean_exp)) +
        geom_line(aes(group=clusters), linewidth=2, alpha=1)+
        geom_point(size = 5)+
        geom_errorbar(aes(ymax = mean_exp + se_exp, ymin = mean_exp - se_exp), width=0.3) +
        #geom_label_repel(data = subset(sigs.averaged, time == "A"), aes(label=genes), nudge_x = 2, max.overlaps = 20, size = 3)+
        scale_y_continuous(name="Normalized Expression", breaks = seq(500, 1100, by = 100), limits = c(500, 1150))+
        scale_x_discrete(name="")+
        theme(axis.text.x = element_text(size = 9),
              axis.title = element_text(size = 13),
              axis.text.y = element_text(size=10),
              legend.position = "bottom", legend.direction = "horizontal")+
        annotate("text", 2, 1100, label = "n = 65 genes", size = 4)+
        ggtitle("Cluster 3")
      
      expression.curve.avgs.clust.3
      
      expression.curve.avgs.clust.4 <- sigs.averaged |>
        filter(clusters == 4) |>
        ggplot(aes(x=time, y=mean_exp)) +
        geom_line(aes(group=clusters), linewidth=2, alpha=1)+
        geom_point(size = 5)+
        geom_errorbar(aes(ymax = mean_exp + se_exp, ymin = mean_exp - se_exp), width=0.3) +
        #geom_label_repel(data = subset(sigs.averaged, time == "A"), aes(label=genes), nudge_x = 2, max.overlaps = 20, size = 3)+
        scale_y_continuous(name="Normalized Expression", breaks = seq(225, 350, by = 25), limits = c(225, 350))+
        scale_x_discrete(name="Time Points")+
        theme(axis.text.x = element_text(size = 9),
              axis.title = element_text(size = 13),
              axis.text.y = element_text(size=10),
              legend.position = "bottom", legend.direction = "horizontal")+
        annotate("text", 2, 350, label = "n = 156 genes", size = 4)+
        ggtitle("Cluster 4")
      
      expression.curve.avgs.clust.4
      
      expression.curve.avgs.clust.5 <- sigs.averaged |>
        filter(clusters == 5) |>
        ggplot(aes(x=time, y=mean_exp)) +
        geom_line(aes(group=clusters), linewidth=2, alpha=1)+
        geom_point(size = 5)+
        geom_errorbar(aes(ymax = mean_exp + se_exp, ymin = mean_exp - se_exp), width=0.3) +
        #geom_label_repel(data = subset(sigs.averaged, time == "A"), aes(label=genes), nudge_x = 2, max.overlaps = 20, size = 3)+
        scale_y_continuous(name="Normalized Expression", breaks = seq(200, 350, by = 50), limits = c(180, 350))+
        scale_x_discrete(name="Time Points")+
        theme(axis.text.x = element_text(size = 9),
              axis.title = element_text(size = 13),
              axis.text.y = element_text(size=10),
              legend.position = "bottom", legend.direction = "horizontal")+
        annotate("text", 2, 350, label = "n = 125 genes", size = 4)+
        ggtitle("Cluster 5")
      
      expression.curve.avgs.clust.5
      
      expression.curve.avgs.clust.6 <- sigs.averaged |>
        filter(clusters == 6) |>
        ggplot(aes(x=time, y=mean_exp)) +
        geom_line(aes(group=clusters), linewidth=2, alpha=1)+
        geom_point(size = 5)+
        geom_errorbar(aes(ymax = mean_exp + se_exp, ymin = mean_exp - se_exp), width=0.3) +
        #geom_label_repel(data = subset(sigs.averaged, time == "A"), aes(label=genes), nudge_x = 2, max.overlaps = 20, size = 3)+
        scale_y_continuous(name="Normalized Expression", breaks = seq(1000, 1900, by = 100), limits = c(1000, 1900))+
        scale_x_discrete(name="Time Points")+
        theme(axis.text.x = element_text(size = 9),
              axis.title = element_text(size = 13),
              axis.text.y = element_text(size=10),
              legend.position = "bottom", legend.direction = "horizontal")+
        annotate("text", 2, 1900, label = "n = 164 genes", size = 4)+
        ggtitle("Cluster 6")

      expression.curve.avgs.clust.6

      expression.curve.avgs.clust.7 <- sigs.averaged |>
        filter(clusters == 7) |>
        ggplot(aes(x=time, y=mean_exp)) +
        geom_line(aes(group=clusters), linewidth=2, alpha=1)+
        geom_point(size = 5)+
        geom_errorbar(aes(ymax = mean_exp + se_exp, ymin = mean_exp - se_exp), width=0.3) +
        #geom_label_repel(data = subset(sigs.averaged, time == "A"), aes(label=genes), nudge_x = 2, max.overlaps = 20, size = 3)+
        scale_y_continuous(name="Normalized Expression", breaks = seq(700, 1500, by = 100), limits = c(700, 1500))+
        scale_x_discrete(name="Time Points")+
        theme(axis.text.x = element_text(size = 9),
              axis.title = element_text(size = 13),
              axis.text.y = element_text(size=10),
              legend.position = "bottom", legend.direction = "horizontal")+
        annotate("text", 2, 1500, label = "n = 99 genes", size = 4)+
        ggtitle("Cluster 7")

      expression.curve.avgs.clust.7
      
      expression.curve.avgs.clust.8 <- sigs.averaged |>
        filter(clusters == 8) |>
        ggplot(aes(x=time, y=mean_exp)) +
        geom_line(aes(group=clusters), linewidth=2, alpha=1)+
        geom_point(size = 5)+
        geom_errorbar(aes(ymax = mean_exp + se_exp, ymin = mean_exp - se_exp), width=0.3) +
        #geom_label_repel(data = subset(sigs.averaged, time == "A"), aes(label=genes), nudge_x = 2, max.overlaps = 20, size = 3)+
        scale_y_continuous(name="Normalized Expression", breaks = seq(1000, 1600, by = 100), limits = c(990, 1600))+
        scale_x_discrete(name="Time Points")+
        theme(axis.text.x = element_text(size = 9),
              axis.title = element_text(size = 13),
              axis.text.y = element_text(size=10),
              legend.position = "bottom", legend.direction = "horizontal")+
        annotate("text", 2, 1600, label = "n = 236 genes", size = 4)+
        ggtitle("Cluster 8")
      
      expression.curve.avgs.clust.8
      
      expression.curve.avgs.clust.9 <- sigs.averaged |>
        filter(clusters == 9) |>
        ggplot(aes(x=time, y=mean_exp)) +
        geom_line(aes(group=clusters), linewidth=2, alpha=1)+
        geom_point(size = 5)+
        geom_errorbar(aes(ymax = mean_exp + se_exp, ymin = mean_exp - se_exp), width=0.3) +
        #geom_label_repel(data = subset(sigs.averaged, time == "A"), aes(label=genes), nudge_x = 2, max.overlaps = 20, size = 3)+
        scale_y_continuous(name="Normalized Expression", breaks = seq(400, 1400, by = 200), limits = c(400, 1400))+
        scale_x_discrete(name="Time Points")+
        theme(axis.text.x = element_text(size = 9),
              axis.title = element_text(size = 13),
              axis.text.y = element_text(size=10),
              legend.position = "bottom", legend.direction = "horizontal")+
        annotate("text", 2, 1400, label = "n = 56 genes", size = 4)+
        ggtitle("Cluster 9")
      
      expression.curve.avgs.clust.9
      # 
      # expression.curve.avgs.clust.67 <- sigs.averaged |>
      #   filter(clusters == "67") |>
      #   ggplot(aes(x=time, y=mean_exp)) +
      #   geom_line(aes(group=clusters), linewidth=2, alpha=1)+
      #   geom_point(size = 5)+
      #   geom_errorbar(aes(ymax = mean_exp + se_exp, ymin = mean_exp - se_exp), width=0.3) +
      #   #geom_label_repel(data = subset(sigs.averaged, time == "A"), aes(label=genes), nudge_x = 2, max.overlaps = 20, size = 3)+
      #   scale_y_continuous(name="Normalized Expression", breaks = seq(1100, 2100, by = 125), limits = c(1100, 2100))+
      #   scale_x_discrete(name="Time Points")+
      #   theme(axis.text.x = element_text(size = 9),
      #         axis.title = element_text(size = 13),
      #         axis.text.y = element_text(size=10),
      #         legend.position = "bottom", legend.direction = "horizontal")+
      #   annotate("text", 2, 1200, label = "n = 243 genes", size = 4)+
      #   ggtitle("Cluster 6")
      # 
      # expression.curve.avgs.clust.67
  
      
    #Export cluster plots as combined figure
      jpeg(filename = "maSigPro cluster patterns 9 clusters.jpg", width = 13, height = 12, units = "in", res = 300)
      expression.curve.avgs.clust.1 + expression.curve.avgs.clust.2 + expression.curve.avgs.clust.3 + 
        expression.curve.avgs.clust.4 + expression.curve.avgs.clust.5 + expression.curve.avgs.clust.6 +
        expression.curve.avgs.clust.7 + expression.curve.avgs.clust.8 + expression.curve.avgs.clust.9 +
        plot_layout(ncol = 3, nrow = 3) + plot_annotation(tag_levels = 'A')
      dev.off()
      
    # #Export just 6 and 7 cluster plots as combined figure for reference before combining 
    #   jpeg(filename = "maSigPro cluster patterns 6 and 7.jpg", width = 9, height = 5, units = "in", res = 300)
    #   expression.curve.avgs.clust.6 + expression.curve.avgs.clust.7 +
    #     plot_layout(ncol = 2, nrow = 1) + plot_annotation(tag_levels = 'A')
    #   dev.off()
      
    # Generate heat map of genes
      #Make color palette for plotting
        br_pal_heat <- met.brewer("OKeeffe1")
        show_col(br_pal_heat)
        my_pal_heat <- br_pal_heat[c(11,9,8,7,6,5,4,3,1)]
        show_col(my_pal_heat)
        
      # Create cluster annotation object for pheatmap
        my_gene_col <- sigs.masigpro[, c(31,32)]
        my_gene_col$clusters <- sub("^", "Cluster_", my_gene_col$clusters)
        my_gene_col$clusters <- as.factor(my_gene_col$clusters)
        my_gene_col$Gene <- NULL
        my_gene_col.ordered <- my_gene_col[order(my_gene_col$clusters), ]
        
      #Create ordered version of the sigs.masigpro data frame based on cluster
        sigs.masigpro.ordered <- sigs.masigpro[order(sigs.masigpro$clusters), ]
        
       
      # Sample colors from met.brewer to create color list for pheatmap annotation
        br_ann_colors <- met.brewer("Redon")
        show_col(br_ann_colors)
       
      # Generate color list 
        ann_colors = list(
          clusters = c(Cluster_1="#5b859e", Cluster_2="#1e395f", Cluster_3= "#75884b", Cluster_4= "#1e5a46",
                       Cluster_5= "#d48f90", Cluster_6= "#732f30", Cluster_7= "#d8b847", Cluster_8= "#b38711", 
                       Cluster_9= "#59385c"))
        
      ### Run pheatmap
        clusters.heatmap<- pheatmap(sigs.masigpro.ordered[1:30], color = my_pal_heat, cluster_rows = F, cluster_cols = F,
                 annotation_colors = ann_colors, cutree_rows = 9,
                 show_rownames=F, border_color=NA, fontsize = 9, scale="row", annotation_row = my_gene_col, 
                 fontsize_row = 10, height=20, 
                 labels_col = c(paste0("Normoxia_", 1:6), paste0("Mild-Hypoxia_", 1:6), paste0("Pcrit_", 1:6), 
                                paste0("Anoxia_", 1:6), paste0("Recovery_", 1:6))
        )
        
      #Export heatmap
        jpeg(filename = "maSigPro clusters heatmap.jpg", width = 5, height = 7, units = "in", res = 600)
        clusters.heatmap
        dev.off()
      
############################ Deseq2 ########################
    
  ## Make contrasts 
        
      # For anoxia vs control time point
          res_anox_vs_cont <- DESeq2::results(dds2, name="anoxia_vs_control", contrast =c("Time", "3", "0"))
          res_anox_vs_cont
        
        # Output tummary table of up and down genes  
          summary(res_anox_vs_cont)
          
        # How many adjusted p-values were less than 0.1?
          sum(res_anox_vs_cont$padj < 0.1, na.rm=TRUE)
        
        # Export results as ordered dataframe by padjusted  
          write.csv(res_anox_vs_cont[order(res_anox_vs_cont$padj),], file = "Sig_genes_DESeq2_anox_vs_cont.csv")
          
        # Get significant gene names for those under 0.1 padj
          sig_names_deseq2_anox_vs_control <- na.omit(rownames(res_anox_vs_cont)[(res_anox_vs_cont$padj < 0.1)])
           
      
    # For anoxia vs recovery time point
          res_anox_vs_recov <- DESeq2::results(dds2, name="anoxia_vs_recovery", contrast =c("Time", "4", "3"))
          res_anox_vs_recov
          
        # Output tummary table of up and down genes  
          summary(res_anox_vs_recov)
          
        # How many adjusted p-values were less than 0.1?
          sum(res_anox_vs_recov$padj < 0.1, na.rm=TRUE)

        # Export results as ordered dataframe by padjusted     
          write.csv(res_anox_vs_recov[order(res_anox_vs_recov$padj),], file = "Sig_genes_DESeq2_anox_vs_recov.csv")
          
        # Get significant gene names for those under 0.1 padj
          sig_names_deseq2_anox_vs_recov <- na.omit(rownames(res_anox_vs_recov)[(res_anox_vs_recov$padj < 0.1)])
          
    # For pcrit vs control time point
          res_pcrit_vs_cont <- DESeq2::results(dds2, name="pcrit_vs_control", contrast =c("Time", "2", "0"))
          res_pcrit_vs_cont
          
        # Output tummary table of up and down genes  
          summary(res_pcrit_vs_cont)
          
        # How many adjusted p-values were less than 0.1?
          sum(res_pcrit_vs_cont$padj < 0.1, na.rm=TRUE)
          
        # Export results as ordered dataframe by padjusted     
          write.csv(res_pcrit_vs_cont[order(res_pcrit_vs_cont$padj),], file = "Sig_genes_DESeq2_pcrit_vs_control.csv")
          
        # Get significant gene names for those under 0.1 padj
          sig_names_deseq2_pcrit_vs_cont <- na.omit(rownames(res_pcrit_vs_cont)[(res_pcrit_vs_cont$padj < 0.1)])
          
    # For recovery vs control time point
          res_recov_vs_cont <- DESeq2::results(dds2, name="recovery_vs_control", contrast =c("Time", "4", "0"))
          res_recov_vs_cont
          
        # Output tummary table of up and down genes  
          summary(res_recov_vs_cont)
          
        # How many adjusted p-values were less than 0.1?
          sum(res_recov_vs_cont$padj < 0.1, na.rm=TRUE)
          
        # Export results as ordered dataframe by padjusted     
          write.csv(res_recov_vs_cont[order(res_recov_vs_cont$padj),], file = "Sig_genes_DESeq2_recovery_vs_control.csv")
          
        # Get significant gene names for those under 0.1 padj
          sig_names_deseq2_recov_vs_cont <- na.omit(rownames(res_recov_vs_cont)[(res_recov_vs_cont$padj < 0.1)])
         
           
     # For 3.5 vs control time point
          res_hypoxia_vs_cont <- DESeq2::results(dds2, name="hypoxia_vs_control", contrast =c("Time", "1", "0"))
          res_hypoxia_vs_cont
          
          # Output tummary table of up and down genes  
          summary(res_hypoxia_vs_cont)
          
          # How many adjusted p-values were less than 0.1?
          sum(res_hypoxia_vs_cont$padj < 0.1, na.rm=TRUE)
          
          # Export results as ordered dataframe by padjusted     
          write.csv(res_hypoxia_vs_cont[order(res_hypoxia_vs_cont$padj),], file = "Sig_genes_DESeq2_hypoxia_vs_control.csv")
          
          # Get significant gene names for those under 0.1 padj
          sig_names_deseq2_hypox_vs_cont <- na.omit(rownames(res_hypoxia_vs_cont)[(res_hypoxia_vs_cont$padj < 0.1)])
        
      # For 3.5 vs anoxia time point
          res_hypox_vs_anox <- DESeq2::results(dds2, name="hypoxia_vs_anoxia", contrast =c("Time", "1", "3"))
          res_hypox_vs_anox
          
          # Output tummary table of up and down genes  
          summary(res_hypox_vs_anox)
          
          # How many adjusted p-values were less than 0.1?
          sum(res_hypox_vs_anox$padj < 0.1, na.rm=TRUE)
          
          # Export results as ordered dataframe by padjusted     
          write.csv(res_hypox_vs_anox[order(res_hypox_vs_anox$padj),], file = "Sig_genes_DESeq2_hypoxia_vs_anoxia.csv")
          
          # Get significant gene names for those under 0.1 padj
          sig_names_deseq2_hypoxia_vs_anoxia <- na.omit(rownames(res_hypox_vs_anox)[(res_hypox_vs_anox$padj < 0.1)])
          
        
##################### Top GO ##############################
  
    #Reorganize gene names to get rid of the PA annotation from gene to go file
      #Read in data
        genes2go <- read.table(file = "Genes2GOterms_SDv2.2.txt", sep = "\t")  
      
    # Adjust sequence name to match my naming scheme (remove the -PA from TCALIF)   
        genes2go$V1 <- str_split_i(genes2go$V1, "-", 1)
     
      # Write it out as a new text file to be read in with readMappings function   
        write.table(genes2go, file = "Genes2GOterms_SDv2.2_noPA.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
        
        
    #Read in gene to GO mappings
      geneID2GO <- readMappings(file = "Genes2GOterms_SDv2.2_noPA.txt")
      
      geneUniverse <- names(geneID2GO) #assign gene names to a list. This will be our global list of genes. 
     
    
      
      
      
####  Top GO maSigPro analysis  ####
      
    #Get first column from maSigPro results w/ gene names and pvalues. Send to list. 
      GeneInt.maSigPro <- as.character(sigs.masigpro.pvalues.and.betas$Gene) 
      
    #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
      geneList.maSigPro <- factor(as.integer(geneUniverse %in% GeneInt.maSigPro)) 
      names(geneList.maSigPro) <- geneUniverse
      
    # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
    # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
    # needed to perform the desired enrichment analysis.
      maSigPro_GO <- new("topGOdata", description = "GO analysis of maSigPro genes", ontology = "BP", nodeSize = 10,
                   allGenes = geneList.maSigPro, annot = annFUN.gene2GO, gene2GO = geneID2GO)
      
    # Once we have an object of class topGOdata we can start with the enrichment analysis.
    # Run the classic fisher exact test to find enriched go terms
        resultFisher.masigpro <- runTest(maSigPro_GO, algorithm = "weight01", statistic = "fisher")
        
        resultFisher.masigpro #View results summary
        
    # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
    # values. In the following example, we list the top 10 significant GO terms identified by the elim method. At
    # the same time we also compare the ranks and the p-values of these GO terms with the ones obtained by the
    # classic method   
    #Adjusted to look just at the classic fisher results since we have just counts
        allRes.masigpro <- GenTable(maSigPro_GO, fisher = resultFisher.masigpro, ranksOf = "fisher", topNodes = 150, numChar = 500)
        colnames(allRes.masigpro)[6] <- "p-value"
        allRes.masigpro$`p-value` <- as.numeric(allRes.masigpro$`p-value`)
       
      # Filter nodes with pvalue greater than 0.05 
        allRes.masigpro <- allRes.masigpro[allRes.masigpro$`p-value` < 0.05,]
        
        
      #Extract names of significant genes in GenTable result. Add to column in allRes.masigpro
        allRes.masigpro$genes <- sapply(allRes.masigpro$GO.ID, function(x)
        {
          genes<-genesInTerm(maSigPro_GO, x)
          genes[[1]][genes[[1]] %in% sigs.masigpro.names] # myGenes is the queried gene list
        })
    
      # Create and export enrichment result as a formatted and interactable table
        maSigPro.allsigs.table <- reactable(allRes.masigpro[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                  wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = TRUE,
                  columns = list(Term = colDef(minWidth = 280),
                                 GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                 )
                  ) %>%
          add_title("maSigPro - all significant genes", font_size = 20) %>% 
          add_subtitle("555 of 849 significant genes", font_size = 16, font_style = "italic") 
        
        
        maSigPro.allsigs.table #View table
        
      #Save table as an html file
        save_reactable_test(maSigPro.allsigs.table, "maSigPro all sig genes table.html")
    
    # investigate how the significant GO terms are distributed over the GO graph
        jpeg(filename = "TopGO subgraph maSigPro.jpg", width = 10, height = 10, units = "in", res = 900)
        par(cex = 0.9)
        showSigOfNodes(maSigPro_GO, score(resultFisher.masigpro), firstSigNodes = 25, useInfo = 'all')
        dev.off()   
       
        
        
    ### Cluster 1 
        
        #Get first column from maSigPro results w/ gene names and pvalues. Send to list. 
          GeneInt.maSigPro.cluster1 <- as.character(sig.gene.clusters[sig.gene.clusters$clusters == "1", "genes"]) 
        
        #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
          geneList.maSigPro.cluster1 <- factor(as.integer(geneUniverse %in% GeneInt.maSigPro.cluster1)) 
          names(geneList.maSigPro.cluster1) <- geneUniverse
        
        # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
        # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
        # needed to perform the desired enrichment analysis.
          maSigPro_GO_cluster1 <- new("topGOdata", description = "Cluster 1 maSigPro genes", ontology = "BP", nodeSize = 10,
                             allGenes = geneList.maSigPro.cluster1, annot = annFUN.gene2GO, gene2GO = geneID2GO)
        
        # Once we have an object of class topGOdata we can start with the enrichment analysis.
        # Run the classic fisher exact test to find enriched go terms
          resultFisher.masigpro.cluster1 <- runTest(maSigPro_GO_cluster1, algorithm = "weight01", statistic = "fisher")
        
          resultFisher.masigpro.cluster1 #View results summary
        
        # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
        # values.  
        #Adjusted to look just at the classic fisher results since we have just counts
          allRes.masigpro.cluster1 <- GenTable(maSigPro_GO_cluster1, fisher = resultFisher.masigpro.cluster1, 
                                               ranksOf = "fisher", topNodes = 150, numChar = 500)
          colnames(allRes.masigpro.cluster1)[6] <- "p-value"
          allRes.masigpro.cluster1$`p-value` <- as.numeric(allRes.masigpro.cluster1$`p-value`)
          
        # Filter nodes with pvalue greater than 0.05 
          allRes.masigpro.cluster1 <- allRes.masigpro.cluster1[allRes.masigpro.cluster1$`p-value` < 0.05,]
        
        
        #Extract names of significant genes in GenTable result. Add to column in allRes.masigpro.cluster1
          allRes.masigpro.cluster1$genes <- sapply(allRes.masigpro.cluster1$GO.ID, function(x)
          {
            genes<-genesInTerm(maSigPro_GO_cluster1, x)
            genes[[1]][genes[[1]] %in% sigs.masigpro.names.cluster1] # myGenes is the queried gene list
          })
          
        # Create and export enrichment result as a formatted and interactable table
          maSigPro.cluster1.table <- reactable(allRes.masigpro.cluster1[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                              wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = FALSE,
                                              columns = list(Term = colDef(minWidth = 400),
                                                             GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                              )
          ) %>%
            add_title("maSigPro - cluster 1", font_size = 20) %>% 
            add_subtitle("112 of 182 significant genes", font_size = 16, font_style = "italic") 
        
        
          maSigPro.cluster1.table #View table
        
        #Save table as an html file
          save_reactable_test(maSigPro.cluster1.table, "maSigPro cluster 1 table.html")
        
        # investigate how the significant GO terms are distributed over the GO graph
          jpeg(filename = "TopGO subgraph maSigPro cluster 1f.jpg", width = 10, height = 10, units = "in", res = 900)
          par(cex = 0.9)
          showSigOfNodes(maSigPro_GO_cluster1, score(resultFisher.masigpro.cluster1), firstSigNodes = 10, useInfo = 'all')
          dev.off() 
          
          
      ### Cluster 2 
          
          #Get first column from maSigPro results w/ gene names and pvalues. Send to list. 
            GeneInt.maSigPro.cluster2 <- as.character(sig.gene.clusters[sig.gene.clusters$clusters == "2", "genes"]) 
          
          #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
            geneList.maSigPro.cluster2 <- factor(as.integer(geneUniverse %in% GeneInt.maSigPro.cluster2)) 
            names(geneList.maSigPro.cluster2) <- geneUniverse
          
          # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
          # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
          # needed to perform the desired enrichment analysis.
            maSigPro_GO_cluster2 <- new("topGOdata", description = "Cluster 2 maSigPro genes", ontology = "BP", nodeSize=10,
                                        allGenes = geneList.maSigPro.cluster2, annot = annFUN.gene2GO, gene2GO = geneID2GO)
          
          # Once we have an object of class topGOdata we can start with the enrichment analysis.
          # Run the classic fisher exact test to find enriched go terms
            resultFisher.masigpro.cluster2 <- runTest(maSigPro_GO_cluster2, algorithm = "weight01", statistic = "fisher")
          
            resultFisher.masigpro.cluster2 #View results summary
          
          # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
          # values.  
          #Adjusted to look just at the classic fisher results since we have just counts
            allRes.masigpro.cluster2 <- GenTable(maSigPro_GO_cluster2, fisher = resultFisher.masigpro.cluster2, 
                                                 ranksOf = "fisher", topNodes = 150, numChar = 500)
            colnames(allRes.masigpro.cluster2)[6] <- "p-value"
            allRes.masigpro.cluster2$`p-value` <- as.numeric(allRes.masigpro.cluster2$`p-value`)
            
          # Filter nodes with pvalue greater than 0.05 
            allRes.masigpro.cluster2 <- allRes.masigpro.cluster2[allRes.masigpro.cluster2$`p-value` < 0.05,]
          
          
          #Extract names of significant genes in GenTable result. Add to column in allRes.masigpro.cluster1
            allRes.masigpro.cluster2$genes <- sapply(allRes.masigpro.cluster2$GO.ID, function(x)
            {
              genes<-genesInTerm(maSigPro_GO_cluster2, x)
              genes[[1]][genes[[1]] %in% sigs.masigpro.names.cluster2] # myGenes is the queried gene list
            })
          
          # Create and export enrichment result as a formatted and interactable table
            maSigPro.cluster2.table <- reactable(allRes.masigpro.cluster2[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                                 wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = FALSE,
                                                 columns = list(Term = colDef(minWidth = 400),
                                                                GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                                 )
            ) %>%
              add_title("maSigPro - cluster 2", font_size = 20) %>% 
              add_subtitle("179 of 264 significant genes", font_size = 16, font_style = "italic") 
          
          
            maSigPro.cluster2.table #View table
          
          #Save table as an html file
            save_reactable_test(maSigPro.cluster2.table, "maSigPro cluster 2 table.html")
          
          # investigate how the significant GO terms are distributed over the GO graph
            jpeg(filename = "TopGO subgraph maSigPro cluster 2.jpg", width = 10, height = 10, units = "in", res = 900)
            par(cex = 0.9)
            showSigOfNodes(maSigPro_GO_cluster2, score(resultFisher.masigpro.cluster2), firstSigNodes = 25, useInfo = 'all')
            dev.off() 
          
            
            
      ### Cluster 3 
            
          #Get first column from maSigPro results w/ gene names and pvalues. Send to list. 
            GeneInt.maSigPro.cluster3 <- as.character(sig.gene.clusters[sig.gene.clusters$clusters == "3", "genes"]) 
          
          #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
            geneList.maSigPro.cluster3 <- factor(as.integer(geneUniverse %in% GeneInt.maSigPro.cluster3)) 
            names(geneList.maSigPro.cluster3) <- geneUniverse
          
          # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
          # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
          # needed to perform the desired enrichment analysis.
            maSigPro_GO_cluster3 <- new("topGOdata", description = "Cluster 3 maSigPro genes", ontology = "BP", nodeSize = 10,
                                        allGenes = geneList.maSigPro.cluster3, annot = annFUN.gene2GO, gene2GO = geneID2GO)
          
          # Once we have an object of class topGOdata we can start with the enrichment analysis.
          # Run the classic fisher exact test to find enriched go terms
            resultFisher.masigpro.cluster3 <- runTest(maSigPro_GO_cluster3, algorithm = "weight01", statistic = "fisher")
            
            resultFisher.masigpro.cluster3 #View results summary
          
          # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
          # values.   
          #Adjusted to look just at the classic fisher results since we have just counts
            allRes.masigpro.cluster3 <- GenTable(maSigPro_GO_cluster3, fisher = resultFisher.masigpro.cluster3, 
                                                 ranksOf = "fisher", topNodes = 150, numChar = 500)
            colnames(allRes.masigpro.cluster3)[6] <- "p-value"
            allRes.masigpro.cluster3$`p-value` <- as.numeric(allRes.masigpro.cluster3$`p-value`)
          
          # Filter nodes with pvalue greater than 0.05 
            allRes.masigpro.cluster3 <- allRes.masigpro.cluster3[allRes.masigpro.cluster3$`p-value` < 0.05,]
            
          #Extract names of significant genes in GenTable result. Add to column in allRes.masigpro.cluster1
            allRes.masigpro.cluster3$genes <- sapply(allRes.masigpro.cluster3$GO.ID, function(x)
            {
              genes<-genesInTerm(maSigPro_GO_cluster3, x)
              genes[[1]][genes[[1]] %in% sigs.masigpro.names.cluster3] # myGenes is the queried gene list
            })
          
          # Create and export enrichment result as a formatted and interactable table
            maSigPro.cluster3.table <- reactable(allRes.masigpro.cluster3[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                                 wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = FALSE,
                                                 columns = list(Term = colDef(minWidth = 400),
                                                                GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                                 )
            ) %>%
              add_title("maSigPro - cluster 3", font_size = 20) %>% 
              add_subtitle("44 of 65 significant genes", font_size = 16, font_style = "italic") 
          
          
            maSigPro.cluster3.table #View table
          
          #Save table as an html file
            save_reactable_test(maSigPro.cluster3.table, "maSigPro cluster 3 table.html")
          
          # investigate how the significant GO terms are distributed over the GO graph
            jpeg(filename = "TopGO subgraph maSigPro cluster 3.jpg", width = 10, height = 10, units = "in", res = 1100)
            par(cex = 0.9)
            showSigOfNodes(maSigPro_GO_cluster3, score(resultFisher.masigpro.cluster3), firstSigNodes = 25, useInfo = 'all')
            dev.off()
            
            
            
       ### Cluster 4
            
          #Get first column from maSigPro results w/ gene names and pvalues. Send to list. 
            GeneInt.maSigPro.cluster4 <- as.character(sig.gene.clusters[sig.gene.clusters$clusters == "4", "genes"]) 
          
          #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
            geneList.maSigPro.cluster4 <- factor(as.integer(geneUniverse %in% GeneInt.maSigPro.cluster4)) 
            names(geneList.maSigPro.cluster4) <- geneUniverse
          
          # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
          # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
          # needed to perform the desired enrichment analysis.
            maSigPro_GO_cluster4 <- new("topGOdata", description = "Cluster 4 maSigPro genes", ontology = "BP", nodeSize=10,
                                        allGenes = geneList.maSigPro.cluster4, annot = annFUN.gene2GO, gene2GO = geneID2GO)
          
          # Once we have an object of class topGOdata we can start with the enrichment analysis.
          # Run the classic fisher exact test to find enriched go terms
            resultFisher.masigpro.cluster4 <- runTest(maSigPro_GO_cluster4, algorithm = "weight01", statistic = "fisher")
            
            resultFisher.masigpro.cluster4 #View results summary
          
          # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
          # values.   
          #Adjusted to look just at the classic fisher results since we have just counts
            allRes.masigpro.cluster4 <- GenTable(maSigPro_GO_cluster4, fisher = resultFisher.masigpro.cluster4, 
                                                 ranksOf = "fisher", topNodes = 150, numChar = 500)
            colnames(allRes.masigpro.cluster4)[6] <- "p-value"
            allRes.masigpro.cluster4$`p-value` <- as.numeric(allRes.masigpro.cluster4$`p-value`)
          
          # Filter nodes with pvalue greater than 0.05 
            allRes.masigpro.cluster4 <- allRes.masigpro.cluster4[allRes.masigpro.cluster4$`p-value` < 0.05,]
            
          #Extract names of significant genes in GenTable result. Add to column in allRes.masigpro.cluster1
            allRes.masigpro.cluster4$genes <- sapply(allRes.masigpro.cluster4$GO.ID, function(x)
            {
              genes<-genesInTerm(maSigPro_GO_cluster4, x)
              genes[[1]][genes[[1]] %in% sigs.masigpro.names.cluster4] # myGenes is the queried gene list
            })
          
          # Create and export enrichment result as a formatted and interactable table
            maSigPro.cluster4.table <- reactable(allRes.masigpro.cluster4[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                                 wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = FALSE,
                                                 columns = list(Term = colDef(minWidth = 400),
                                                                GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                                 )
            ) %>%
              add_title("maSigPro - cluster 4", font_size = 20) %>% 
              add_subtitle("105 of 156 significant genes", font_size = 16, font_style = "italic") 
            
            
            maSigPro.cluster4.table #View table
            
          
          #Save table as an html file
            save_reactable_test(maSigPro.cluster4.table, "maSigPro cluster 4 table.html")
          
          # investigate how the significant GO terms are distributed over the GO graph
            jpeg(filename = "TopGO subgraph maSigPro cluster 4.jpg", width = 10, height = 10, units = "in", res = 900)
            par(cex = 0.9)
            showSigOfNodes(maSigPro_GO_cluster4, score(resultFisher.masigpro.cluster4), firstSigNodes = 10, useInfo = 'all')
            dev.off()
            
            
            
            
        ### Cluster 5 
            
          #Get first column from maSigPro results w/ gene names and pvalues. Send to list. 
            GeneInt.maSigPro.cluster5 <- as.character(sig.gene.clusters[sig.gene.clusters$clusters == "5", "genes"]) 
            
          #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
            geneList.maSigPro.cluster5 <- factor(as.integer(geneUniverse %in% GeneInt.maSigPro.cluster5)) 
            names(geneList.maSigPro.cluster5) <- geneUniverse
            
          # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
          # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
          # needed to perform the desired enrichment analysis.
            maSigPro_GO_cluster5 <- new("topGOdata", description = "Cluster 5 maSigPro genes", ontology = "BP", nodeSize = 10,
                                        allGenes = geneList.maSigPro.cluster5, annot = annFUN.gene2GO, gene2GO = geneID2GO)
            
          # Once we have an object of class topGOdata we can start with the enrichment analysis.
          # Run the classic fisher exact test to find enriched go terms
            resultFisher.masigpro.cluster5 <- runTest(maSigPro_GO_cluster5, algorithm = "weight01", statistic = "fisher")
            
            resultFisher.masigpro.cluster5 #View results summary
            
          # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
          # values.   
          #Adjusted to look just at the classic fisher results since we have just counts
            allRes.masigpro.cluster5 <- GenTable(maSigPro_GO_cluster5, fisher = resultFisher.masigpro.cluster5, 
                                                 ranksOf = "fisher", topNodes = 150, numChar = 500)
            colnames(allRes.masigpro.cluster5)[6] <- "p-value"
            allRes.masigpro.cluster5$`p-value` <- as.numeric(allRes.masigpro.cluster5$`p-value`)
            
          # Filter nodes with pvalue greater than 0.05 
            allRes.masigpro.cluster5 <- allRes.masigpro.cluster5[allRes.masigpro.cluster5$`p-value` < 0.05,]  
            
          #Extract names of significant genes in GenTable result. Add to column in allRes.masigpro.cluster1
            allRes.masigpro.cluster5$genes <- sapply(allRes.masigpro.cluster5$GO.ID, function(x)
            {
              genes<-genesInTerm(maSigPro_GO_cluster5, x)
              genes[[1]][genes[[1]] %in% sigs.masigpro.names.cluster5] # myGenes is the queried gene list
            })
            
          # Create and export enrichment result as a formatted and interactable table
            maSigPro.cluster5.table <- reactable(allRes.masigpro.cluster5[,1:7], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                                 wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = FALSE,
                                                 columns = list(Term = colDef(minWidth = 400),
                                                                GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                                 )
            ) %>%
              add_title("maSigPro - cluster 5", font_size = 20) %>% 
              add_subtitle("82 of 125 significant genes", font_size = 16, font_style = "italic") 
            
            
            maSigPro.cluster5.table #View table
            
          #Save table as an html file
            save_reactable_test(maSigPro.cluster5.table, "maSigPro cluster 5 table.html")
            
          # investigate how the significant GO terms are distributed over the GO graph
            jpeg(filename = "TopGO subgraph maSigPro cluster 5.jpg", width = 10, height = 10, units = "in", res = 900)
            par(cex = 0.9)
            showSigOfNodes(maSigPro_GO_cluster5, score(resultFisher.masigpro.cluster5), firstSigNodes = 10, useInfo = 'all')
            dev.off()
            
            
            
            
        ### Cluster 6
            
          #Get first column from maSigPro results w/ gene names and pvalues. Send to list. 
            GeneInt.maSigPro.cluster6 <- as.character(sig.gene.clusters[sigs.masigpro$clusters == "6", "genes"]) 
            
          #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
            geneList.maSigPro.cluster6 <- factor(as.integer(geneUniverse %in% GeneInt.maSigPro.cluster6)) 
            names(geneList.maSigPro.cluster6) <- geneUniverse
            
          # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
          # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
          # needed to perform the desired enrichment analysis.
            maSigPro_GO_cluster6 <- new("topGOdata", description = "Cluster 6 maSigPro genes", ontology = "BP", nodeSize=10,
                                        allGenes = geneList.maSigPro.cluster6, annot = annFUN.gene2GO, gene2GO = geneID2GO)
            
          # Once we have an object of class topGOdata we can start with the enrichment analysis.
          # Run the classic fisher exact test to find enriched go terms
            resultFisher.masigpro.cluster6 <- runTest(maSigPro_GO_cluster6, algorithm = "weight01", statistic = "fisher")
            
            resultFisher.masigpro.cluster6 #View results summary
            
          # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
          # values. 
          #Adjusted to look just at the classic fisher results since we have just counts
            allRes.masigpro.cluster6 <- GenTable(maSigPro_GO_cluster6, fisher = resultFisher.masigpro.cluster6, 
                                                 ranksOf = "fisher", topNodes = 150, numChar = 500)
            colnames(allRes.masigpro.cluster6)[6] <- "p-value"
            allRes.masigpro.cluster6$`p-value` <- as.numeric(allRes.masigpro.cluster6$`p-value`)
            
          # Filter nodes with pvalue greater than 0.05 
            allRes.masigpro.cluster6 <- allRes.masigpro.cluster6[allRes.masigpro.cluster6$`p-value` < 0.05,]
            
          #Extract names of significant genes in GenTable result. Add to column in allRes.masigpro.cluster1
            allRes.masigpro.cluster6$genes <- sapply(allRes.masigpro.cluster6$GO.ID, function(x)
            {
              genes<-genesInTerm(maSigPro_GO_cluster6, x)
              genes[[1]][genes[[1]] %in% sigs.masigpro.names.cluster6] # myGenes is the queried gene list
            })
            
          # Create and export enrichment result as a formatted and interactable table
            maSigPro.cluster6.table <- reactable(allRes.masigpro.cluster6[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                                 wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = FALSE,
                                                 columns = list(Term = colDef(minWidth = 400),
                                                                GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                                 )
            ) %>%
              add_title("maSigPro - cluster 6", font_size = 20) %>% 
              add_subtitle("124 of 164 significant genes", font_size = 16, font_style = "italic") 
            
            
            maSigPro.cluster6.table #View table
            
          #Save table as an html file
            save_reactable_test(maSigPro.cluster6.table, "maSigPro cluster 6 table.html")
            
          # investigate how the significant GO terms are distributed over the GO graph
            jpeg(filename = "TopGO subgraph maSigPro cluster 6.jpg", width = 10, height = 10, units = "in", res = 900)
            par(cex = 0.9)
            showSigOfNodes(maSigPro_GO_cluster6, score(resultFisher.masigpro.cluster6), firstSigNodes = 10, useInfo = 'all')
            dev.off()
            
            
        ### Cluster 7
            
          #Get first column from maSigPro results w/ gene names and pvalues. Send to list. 
            GeneInt.maSigPro.cluster7 <- as.character(sig.gene.clusters[sigs.masigpro$clusters == "7", "genes"]) 
            
          #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
            geneList.maSigPro.cluster7 <- factor(as.integer(geneUniverse %in% GeneInt.maSigPro.cluster7)) 
            names(geneList.maSigPro.cluster7) <- geneUniverse
            
          # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
          # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
          # needed to perform the desired enrichment analysis.
            maSigPro_GO_cluster7 <- new("topGOdata", description = "Cluster 7 maSigPro genes", ontology = "BP", nodeSize=10,
                                        allGenes = geneList.maSigPro.cluster7, annot = annFUN.gene2GO, gene2GO = geneID2GO)
            
          # Once we have an object of class topGOdata we can start with the enrichment analysis.
          # Run the classic fisher exact test to find enriched go terms
            resultFisher.masigpro.cluster7 <- runTest(maSigPro_GO_cluster7, algorithm = "weight01", statistic = "fisher")
            
            resultFisher.masigpro.cluster7 #View results summary
            
          # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
          # values. 
          #Adjusted to look just at the classic fisher results since we have just counts
            allRes.masigpro.cluster7 <- GenTable(maSigPro_GO_cluster7, fisher = resultFisher.masigpro.cluster7, 
                                                 ranksOf = "fisher", topNodes = 150, numChar = 500)
            colnames(allRes.masigpro.cluster7)[6] <- "p-value"
            allRes.masigpro.cluster7$`p-value` <- as.numeric(allRes.masigpro.cluster7$`p-value`)
            
          # Filter nodes with pvalue greater than 0.05 
            allRes.masigpro.cluster7 <- allRes.masigpro.cluster7[allRes.masigpro.cluster7$`p-value` < 0.05,]
            
          #Extract names of significant genes in GenTable result. Add to column in allRes.masigpro.cluster1
            allRes.masigpro.cluster7$genes <- sapply(allRes.masigpro.cluster7$GO.ID, function(x)
            {
              genes<-genesInTerm(maSigPro_GO_cluster7, x)
              genes[[1]][genes[[1]] %in% sigs.masigpro.names.cluster7] # myGenes is the queried gene list
            })
            
          # Create and export enrichment result as a formatted and interactable table
            maSigPro.cluster7.table <- reactable(allRes.masigpro.cluster7[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                                 wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = FALSE,
                                                 columns = list(Term = colDef(minWidth = 400),
                                                                GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                                 )
            ) %>%
              add_title("maSigPro - cluster 7", font_size = 20) %>% 
              add_subtitle("61 of 99 significant genes", font_size = 16, font_style = "italic") 
            
            
            maSigPro.cluster7.table #View table
            
          #Save table as an html file
            save_reactable_test(maSigPro.cluster7.table, "maSigPro cluster 7 table.html")
            
          # investigate how the significant GO terms are distributed over the GO graph
            jpeg(filename = "TopGO subgraph maSigPro cluster 7.jpg", width = 10, height = 10, units = "in", res = 900)
            par(cex = 0.9)
            showSigOfNodes(maSigPro_GO_cluster7, score(resultFisher.masigpro.cluster7), firstSigNodes = 10, useInfo = 'all')
            dev.off()
            
            
        ### Cluster 8
            
          #Get first column from maSigPro results w/ gene names and pvalues. Send to list. 
            GeneInt.maSigPro.cluster8 <- as.character(sig.gene.clusters[sigs.masigpro$clusters == "8", "genes"]) 
            
          #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
            geneList.maSigPro.cluster8 <- factor(as.integer(geneUniverse %in% GeneInt.maSigPro.cluster8)) 
            names(geneList.maSigPro.cluster8) <- geneUniverse
            
          # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
          # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
          # needed to perform the desired enrichment analysis.
            maSigPro_GO_cluster8 <- new("topGOdata", description = "Cluster 8 maSigPro genes", ontology = "BP", nodeSize=10,
                                        allGenes = geneList.maSigPro.cluster8, annot = annFUN.gene2GO, gene2GO = geneID2GO)
            
          # Once we have an object of class topGOdata we can start with the enrichment analysis.
          # Run the classic fisher exact test to find enriched go terms
            resultFisher.masigpro.cluster8 <- runTest(maSigPro_GO_cluster8, algorithm = "weight01", statistic = "fisher")
            
            resultFisher.masigpro.cluster8 #View results summary
            
          # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
          # values. 
          #Adjusted to look just at the classic fisher results since we have just counts
            allRes.masigpro.cluster8 <- GenTable(maSigPro_GO_cluster8, fisher = resultFisher.masigpro.cluster8, 
                                                 ranksOf = "fisher", topNodes = 150, numChar = 500)
            colnames(allRes.masigpro.cluster8)[6] <- "p-value"
            allRes.masigpro.cluster8$`p-value` <- as.numeric(allRes.masigpro.cluster8$`p-value`)
            
          # Filter nodes with pvalue greater than 0.05 
            allRes.masigpro.cluster8 <- allRes.masigpro.cluster8[allRes.masigpro.cluster8$`p-value` < 0.05,]
            
          #Extract names of significant genes in GenTable result. Add to column in allRes.masigpro.cluster1
            allRes.masigpro.cluster8$genes <- sapply(allRes.masigpro.cluster8$GO.ID, function(x)
            {
              genes<-genesInTerm(maSigPro_GO_cluster8, x)
              genes[[1]][genes[[1]] %in% sigs.masigpro.names.cluster8] # myGenes is the queried gene list
            })
            
          # Create and export enrichment result as a formatted and interactable table
            maSigPro.cluster8.table <- reactable(allRes.masigpro.cluster8[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                                 wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = FALSE,
                                                 columns = list(Term = colDef(minWidth = 400),
                                                                GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                                 )
            ) %>%
              add_title("maSigPro - cluster 8", font_size = 20) %>% 
              add_subtitle("166 of 236 significant genes", font_size = 16, font_style = "italic") 
            
            
            maSigPro.cluster8.table #View table
            
          #Save table as an html file
            save_reactable_test(maSigPro.cluster8.table, "maSigPro cluster 8 table.html")
            
          # investigate how the significant GO terms are distributed over the GO graph
            jpeg(filename = "TopGO subgraph maSigPro cluster 8.jpg", width = 10, height = 10, units = "in", res = 900)
            par(cex = 0.9)
            showSigOfNodes(maSigPro_GO_cluster8, score(resultFisher.masigpro.cluster8), firstSigNodes = 10, useInfo = 'all')
            dev.off()
            
            
            
        ### Cluster 9
            
          #Get first column from maSigPro results w/ gene names and pvalues. Send to list. 
            GeneInt.maSigPro.cluster9 <- as.character(sig.gene.clusters[sigs.masigpro$clusters == "9", "genes"]) 
            
          #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
            geneList.maSigPro.cluster9 <- factor(as.integer(geneUniverse %in% GeneInt.maSigPro.cluster9)) 
            names(geneList.maSigPro.cluster9) <- geneUniverse
            
          # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
          # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
          # needed to perform the desired enrichment analysis.
            maSigPro_GO_cluster9 <- new("topGOdata", description = "Cluster 9 maSigPro genes", ontology = "BP", nodeSize=10,
                                        allGenes = geneList.maSigPro.cluster9, annot = annFUN.gene2GO, gene2GO = geneID2GO)
            
          # Once we have an object of class topGOdata we can start with the enrichment analysis.
          # Run the classic fisher exact test to find enriched go terms
            resultFisher.masigpro.cluster9 <- runTest(maSigPro_GO_cluster9, algorithm = "weight01", statistic = "fisher")
            
            resultFisher.masigpro.cluster9 #View results summary
            
          # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
          # values. 
          #Adjusted to look just at the classic fisher results since we have just counts
            allRes.masigpro.cluster9 <- GenTable(maSigPro_GO_cluster9, fisher = resultFisher.masigpro.cluster9, 
                                                 ranksOf = "fisher", topNodes = 150, numChar = 500)
            colnames(allRes.masigpro.cluster9)[6] <- "p-value"
            allRes.masigpro.cluster9$`p-value` <- as.numeric(allRes.masigpro.cluster9$`p-value`)
            
          # Filter nodes with pvalue greater than 0.05 
            allRes.masigpro.cluster9 <- allRes.masigpro.cluster9[allRes.masigpro.cluster9$`p-value` < 0.05,]
            
          #Extract names of significant genes in GenTable result. Add to column in allRes.masigpro.cluster1
            allRes.masigpro.cluster9$genes <- sapply(allRes.masigpro.cluster9$GO.ID, function(x)
            {
              genes<-genesInTerm(maSigPro_GO_cluster9, x)
              genes[[1]][genes[[1]] %in% sigs.masigpro.names.cluster9] # myGenes is the queried gene list
            })
            
          # Create and export enrichment result as a formatted and interactable table
            maSigPro.cluster9.table <- reactable(allRes.masigpro.cluster9[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                                 wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = FALSE,
                                                 columns = list(Term = colDef(minWidth = 400),
                                                                GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                                 )
            ) %>%
              add_title("maSigPro - cluster 9", font_size = 20) %>% 
              add_subtitle("27 of 56 significant genes", font_size = 16, font_style = "italic") 
            
            
            maSigPro.cluster9.table #View table
            
          #Save table as an html file
            save_reactable_test(maSigPro.cluster9.table, "maSigPro cluster 9 table.html")
            
          # investigate how the significant GO terms are distributed over the GO graph
            jpeg(filename = "TopGO subgraph maSigPro cluster 9.jpg", width = 10, height = 10, units = "in", res = 900)
            par(cex = 0.9)
            showSigOfNodes(maSigPro_GO_cluster9, score(resultFisher.masigpro.cluster9), firstSigNodes = 10, useInfo = 'all')
            dev.off()
        
            

#### Top GO DESeq2 analysis ####
        
    # Anoxia vs control 
        
      #Using sig_names_deseq list that we made before
      #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
        geneList.deseq.anox.control <- factor(as.integer(geneUniverse %in% sig_names_deseq2_anox_vs_control)) 
        names(geneList.deseq.anox.control) <- geneUniverse
        
      # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
      # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
      # needed to perform the desired enrichment analysis.
        DESeq_GO_anox_cont <- new("topGOdata", description = "GO analysis of DESeq Anox vs control genes", ontology = "BP", nodeSize=10,
                        allGenes = geneList.deseq.anox.control, annot = annFUN.gene2GO, gene2GO = geneID2GO)
        
      # Once we have an object of class topGOdata we can start with the enrichment analysis.
      # Run the classic fisher exact test to find enriched go terms
        resultFisher.deseq.anox.cont <- runTest(DESeq_GO_anox_cont, algorithm = "weight01", statistic = "fisher")
        
        resultFisher.deseq.anox.cont #View results summary
        
      # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
      # values.   
      #Adjusted to look just at the classic fisher results since we have just counts
        allRes.deseq.anox.cont <- GenTable(DESeq_GO_anox_cont, fisher = resultFisher.deseq.anox.cont, ranksOf = "fisher", 
                                           topNodes = 150, numChar = 500)
        colnames(allRes.deseq.anox.cont)[6] <- "p-value"
        allRes.deseq.anox.cont$`p-value` <- as.numeric(allRes.deseq.anox.cont$`p-value`)
        
      # Filter nodes with pvalue greater than 0.05 
        allRes.deseq.anox.cont <- allRes.deseq.anox.cont[allRes.deseq.anox.cont$`p-value` < 0.05,]
        
      #Extract names of significant genes in GenTable result. Add to column in allRes.deseq.anox.cont
        allRes.deseq.anox.cont$genes <- sapply(allRes.deseq.anox.cont$GO.ID, function(x)
        {
          genes<-genesInTerm(DESeq_GO_anox_cont, x)
          genes[[1]][genes[[1]] %in% sig_names_deseq2_anox_vs_control] 
        })
        
      # Create and export enrichment result as a formatted and interactable table
        deseq.anox.cont.table <- reactable(allRes.deseq.anox.cont[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                         wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = FALSE,
                                         columns = list(Term = colDef(minWidth = 400),
                                                        GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                         )
        ) %>%
          add_title("DESeq2 - Anoxia vs control significant genes", font_size = 20) %>% 
          add_subtitle("569 of 865 significant genes", font_size = 16, font_style = "italic") 
        
        
        deseq.anox.cont.table #View table
        
      #Save table as an html file
        save_reactable_test(deseq.anox.cont.table, "DESeq2 anoxia vs control table.html")
        
      # investigate how the significant GO terms are distributed over the GO graph
        jpeg(filename = "TopGO subgraph DESeq2 anoxia vs control.jpg", width = 10, height = 10, units = "in", res = 900)
        par(cex = 0.9)
        showSigOfNodes(DESeq_GO_anox_cont, score(resultFisher.deseq.anox.cont), firstSigNodes = 25, useInfo = 'all')
        dev.off() 
        
        
    # Anoxia vs recovery 
        
      #Using sig_names_deseq list that we made before
      #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
        geneList.deseq.anox.recov <- factor(as.integer(geneUniverse %in% sig_names_deseq2_anox_vs_recov)) 
        names(geneList.deseq.anox.recov) <- geneUniverse
        
      # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
      # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
      # needed to perform the desired enrichment analysis.
        DESeq_GO_anox_recov <- new("topGOdata", description = "GO analysis of DESeq Anox vs recovery genes", ontology = "BP", nodeSize=10,
                                  allGenes = geneList.deseq.anox.recov, annot = annFUN.gene2GO, gene2GO = geneID2GO)
        
      # Once we have an object of class topGOdata we can start with the enrichment analysis.
      # Run the classic fisher exact test to find enriched go terms
        resultFisher.deseq.anox.recov <- runTest(DESeq_GO_anox_recov, algorithm = "weight01", statistic = "fisher")
        
        resultFisher.deseq.anox.recov #View results summary
        
      # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
      # values.  
      #Adjusted to look just at the classic fisher results since we have just counts
        allRes.deseq.anox.recov <- GenTable(DESeq_GO_anox_recov, fisher = resultFisher.deseq.anox.recov, ranksOf = "fisher", 
                                           topNodes = 150, numChar = 500)
        colnames(allRes.deseq.anox.recov)[6] <- "p-value"
        allRes.deseq.anox.recov$`p-value` <- as.numeric(allRes.deseq.anox.recov$`p-value`)
        
      # Filter nodes with pvalue greater than 0.05 
        allRes.deseq.anox.recov <- allRes.deseq.anox.recov[allRes.deseq.anox.recov$`p-value` < 0.05,]
        
        
      #Extract names of significant genes in GenTable result. Add to column in allRes.deseq.anox.recov
        allRes.deseq.anox.recov$genes <- sapply(allRes.deseq.anox.recov$GO.ID, function(x)
        {
          genes<-genesInTerm(DESeq_GO_anox_recov, x)
          genes[[1]][genes[[1]] %in% sig_names_deseq2_anox_vs_recov] 
        })
        
      # Create and export enrichment result as a formatted and interactable table
        deseq.anox.recov.table <- reactable(allRes.deseq.anox.recov[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                           wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = FALSE,
                                           columns = list(Term = colDef(minWidth = 400),
                                                          GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                           )
        ) %>%
          add_title("DESeq2 - Anoxia vs recovery significant genes", font_size = 20) %>% 
          add_subtitle("411 of 634 significant genes", font_size = 16, font_style = "italic") 
        
        
        deseq.anox.recov.table #View table
        
      #Save table as an html file
        save_reactable_test(deseq.anox.recov.table, "DESeq2 anoxia vs recovery table.html")
        
      # investigate how the significant GO terms are distributed over the GO graph
        jpeg(filename = "TopGO subgraph DESeq2 anoxia vs recovery.jpg", width = 10, height = 10, units = "in", res = 900)
        par(cex = 0.9)
        showSigOfNodes(DESeq_GO_anox_recov, score(resultFisher.deseq.anox.recov), firstSigNodes = 25, useInfo = 'all')
        dev.off()
        
        
        
    # Pcrit vs control 
        
      #Using sig_names_deseq list that we made before
      #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
        geneList.deseq.pcrit.cont <- factor(as.integer(geneUniverse %in% sig_names_deseq2_pcrit_vs_cont)) 
        names(geneList.deseq.pcrit.cont) <- geneUniverse
        
      # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
      # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
      # needed to perform the desired enrichment analysis.
        DESeq_GO_pcrit_cont <- new("topGOdata", description = "GO analysis of DESeq pcrit vs control genes", ontology = "BP", nodeSize=10,
                                   allGenes = geneList.deseq.pcrit.cont, annot = annFUN.gene2GO, gene2GO = geneID2GO)
        
      # Once we have an object of class topGOdata we can start with the enrichment analysis.
      # Run the classic fisher exact test to find enriched go terms
        resultFisher.deseq.pcrit.cont <- runTest(DESeq_GO_pcrit_cont, algorithm = "weight01", statistic = "fisher")
        
        resultFisher.deseq.pcrit.cont #View results summary
        
      # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
      # values. 
      #Adjusted to look just at the classic fisher results since we have just counts
        allRes.deseq.pcrit.cont <- GenTable(DESeq_GO_pcrit_cont, fisher = resultFisher.deseq.pcrit.cont, ranksOf = "fisher", 
                                            topNodes = 150, numChar = 500)
        colnames(allRes.deseq.pcrit.cont)[6] <- "p-value"
        allRes.deseq.pcrit.cont$`p-value` <- as.numeric(allRes.deseq.pcrit.cont$`p-value`)
        
      # Filter nodes with pvalue greater than 0.05 
        allRes.deseq.pcrit.cont <- allRes.deseq.pcrit.cont[allRes.deseq.pcrit.cont$`p-value` < 0.05,]
        
        
      #Extract names of significant genes in GenTable result. Add to column in allRes.deseq.pcrit.cont
        allRes.deseq.pcrit.cont$genes <- sapply(allRes.deseq.pcrit.cont$GO.ID, function(x)
        {
          genes<-genesInTerm(DESeq_GO_pcrit_cont, x)
          genes[[1]][genes[[1]] %in% sig_names_deseq2_pcrit_vs_cont] 
        })
        
      # Create and export enrichment result as a formatted and interactable table
        deseq.pcrit.cont.table <- reactable(allRes.deseq.pcrit.cont[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                            wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = FALSE,
                                            columns = list(Term = colDef(minWidth = 400),
                                                           GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                            )
        ) %>%
          add_title("DESeq2 - Pcrit vs control significant genes", font_size = 20) %>% 
          add_subtitle("310 of 489 significant genes", font_size = 16, font_style = "italic") 
        
        
        deseq.pcrit.cont.table #View table
        
      #Save table as an html file
        save_reactable_test(deseq.pcrit.cont.table, "DESeq2 pcrit vs control table.html")
        
      # investigate how the significant GO terms are distributed over the GO graph
        jpeg(filename = "TopGO subgraph DESeq2 pcrit vs control.jpg", width = 10, height = 10, units = "in", res = 900)
        par(cex = 0.9)
        showSigOfNodes(DESeq_GO_pcrit_cont, score(resultFisher.deseq.pcrit.cont), firstSigNodes = 20, useInfo = 'all')
        dev.off()
        
        
     # Hypoxia vs control 
        
      #Using sig_names_deseq list that we made before
      #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
        geneList.deseq.hypox.cont <- factor(as.integer(geneUniverse %in% sig_names_deseq2_hypox_vs_cont)) 
        names(geneList.deseq.hypox.cont) <- geneUniverse
        
      # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
      # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
      # needed to perform the desired enrichment analysis.
        DESeq_GO_hypox_cont <- new("topGOdata", description = "GO analysis of DESeq hypoxia vs control genes", ontology = "BP", nodeSize=10,
                                   allGenes = geneList.deseq.hypox.cont, annot = annFUN.gene2GO, gene2GO = geneID2GO)
        
      # Once we have an object of class topGOdata we can start with the enrichment analysis.
      # Run the classic fisher exact test to find enriched go terms
        resultFisher.deseq.hypox.cont <- runTest(DESeq_GO_hypox_cont, algorithm = "weight01", statistic = "fisher")
        
        resultFisher.deseq.hypox.cont #View results summary
        
      # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
      # values. 
      #Adjusted to look just at the classic fisher results since we have just counts
        allRes.deseq.hypox.cont <- GenTable(DESeq_GO_hypox_cont, fisher = resultFisher.deseq.hypox.cont, ranksOf = "fisher", 
                                            topNodes = 150, numChar = 500)
        colnames(allRes.deseq.hypox.cont)[6] <- "p-value"
        allRes.deseq.hypox.cont$`p-value` <- as.numeric(allRes.deseq.hypox.cont$`p-value`)
        
      # Filter nodes with pvalue greater than 0.05 
        allRes.deseq.hypox.cont <- allRes.deseq.hypox.cont[allRes.deseq.hypox.cont$`p-value` < 0.05,]
        
        
      #Extract names of significant genes in GenTable result. Add to column in allRes.deseq.hypox.cont
        allRes.deseq.hypox.cont$genes <- sapply(allRes.deseq.hypox.cont$GO.ID, function(x)
        {
          genes<-genesInTerm(DESeq_GO_hypox_cont, x)
          genes[[1]][genes[[1]] %in% sig_names_deseq2_hypox_vs_cont] 
        })
        
      # Create and export enrichment result as a formatted and interactable table
        deseq.hypox.cont.table <- reactable(allRes.deseq.hypox.cont[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                            wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = FALSE,
                                            columns = list(Term = colDef(minWidth = 400),
                                                           GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                            )
        ) %>%
          add_title("DESeq2 - Hypoxia vs control significant genes", font_size = 20) %>% 
          add_subtitle("126 of 203 significant genes", font_size = 16, font_style = "italic") 
        
        
        deseq.hypox.cont.table #View table
        
      # Save table as an html file
        save_reactable_test(deseq.hypox.cont.table, "DESeq2 hypoxia vs control table.html")
        
      # investigate how the significant GO terms are distributed over the GO graph
        jpeg(filename = "TopGO subgraph DESeq2 hypoxia vs control.jpg", width = 10, height = 10, units = "in", res = 900)
        par(cex = 0.9)
        showSigOfNodes(DESeq_GO_hypox_cont, score(resultFisher.deseq.hypox.cont), firstSigNodes = 25, useInfo = 'all')
        dev.off()
        
        
        
    # Recovery vs control 
        
      #Using sig_names_deseq list that we made before
      #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
        geneList.deseq.recov.cont <- factor(as.integer(geneUniverse %in% sig_names_deseq2_recov_vs_cont)) 
        names(geneList.deseq.recov.cont) <- geneUniverse
        
      # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
      # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
      # needed to perform the desired enrichment analysis.
        DESeq_GO_recov_cont <- new("topGOdata", description = "GO analysis of DESeq recovery vs control genes", ontology = "BP", nodeSize=10,
                                   allGenes = geneList.deseq.recov.cont, annot = annFUN.gene2GO, gene2GO = geneID2GO)
        
      # Once we have an object of class topGOdata we can start with the enrichment analysis.
      # Run the classic fisher exact test to find enriched go terms
        resultFisher.deseq.recov.cont <- runTest(DESeq_GO_recov_cont, algorithm = "weight01", statistic = "fisher")
        
        resultFisher.deseq.recov.cont #View results summary
        
      # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
      # values. 
      #Adjusted to look just at the classic fisher results since we have just counts
        allRes.deseq.recov.cont <- GenTable(DESeq_GO_recov_cont, fisher = resultFisher.deseq.recov.cont, ranksOf = "fisher", 
                                            topNodes = 150, numChar = 500)
        colnames(allRes.deseq.recov.cont)[6] <- "p-value"
        allRes.deseq.recov.cont$`p-value` <- as.numeric(allRes.deseq.recov.cont$`p-value`)
        
      # Filter nodes with pvalue greater than 0.05 
        allRes.deseq.recov.cont <- allRes.deseq.recov.cont[allRes.deseq.recov.cont$`p-value` < 0.05,]
        
        
      #Extract names of significant genes in GenTable result. Add to column in allRes.deseq.recov.cont
        allRes.deseq.recov.cont$genes <- sapply(allRes.deseq.recov.cont$GO.ID, function(x)
        {
          genes<-genesInTerm(DESeq_GO_recov_cont, x)
          genes[[1]][genes[[1]] %in% sig_names_deseq2_recov_vs_cont] 
        })
        
      # Create and export enrichment result as a formatted and interactable table
        deseq.recov.cont.table <- reactable(allRes.deseq.recov.cont[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                            wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = FALSE,
                                            columns = list(Term = colDef(minWidth = 400),
                                                           GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                            )
        ) %>%
          add_title("DESeq2 - Recovery vs control significant genes", font_size = 20) %>% 
          add_subtitle("513 of 787 significant genes", font_size = 16, font_style = "italic") 
        
        
        deseq.recov.cont.table #View table
        
      # Save table as an html file
        save_reactable_test(deseq.recov.cont.table, "DESeq2 recovery vs control table.html")
        
      # investigate how the significant GO terms are distributed over the GO graph
        jpeg(filename = "TopGO subgraph DESeq2 recovery vs control.jpg", width = 10, height = 10, units = "in", res = 900)
        par(cex = 0.9)
        showSigOfNodes(DESeq_GO_recov_cont, score(resultFisher.deseq.recov.cont), firstSigNodes = 25, useInfo = 'all')
        dev.off()
        
        
        
    # Hypoxia vs anoxia 
        
        #Using sig_names_deseq list that we made before
        #Make named vector with the gene universe where genes of interest are coded with a 1 so the 'new' function knows to focus on them. 
        geneList.deseq.hypox.anox <- factor(as.integer(geneUniverse %in% sig_names_deseq2_hypoxia_vs_anoxia)) 
        names(geneList.deseq.hypox.anox) <- geneUniverse
        
        # We now have all data necessary to build an object of type topGOdata. This object will contain all gene
        # identifiers and their scores, the GO annotations, the GO hierarchical structure and all other information
        # needed to perform the desired enrichment analysis.
        DESeq_GO_hypox_anox <- new("topGOdata", description = "GO analysis of DESeq hypoxia vs anoxia genes", 
                                   ontology = "BP", nodeSize=10,
                                   allGenes = geneList.deseq.hypox.anox, annot = annFUN.gene2GO, gene2GO = geneID2GO)
        
        # Once we have an object of class topGOdata we can start with the enrichment analysis.
        # Run the classic fisher exact test to find enriched go terms
        resultFisher.deseq.hypox.anox <- runTest(DESeq_GO_hypox_anox, algorithm = "weight01", statistic = "fisher")
        
        resultFisher.deseq.hypox.anox #View results summary
        
        # GenTable is an easy to use function for analysing the most significant GO terms and the corresponding p
        # values. 
        #Adjusted to look just at the classic fisher results since we have just counts
        allRes.deseq.hypox.anox <- GenTable(DESeq_GO_hypox_anox, fisher = resultFisher.deseq.hypox.anox, ranksOf = "fisher", 
                                            topNodes = 150, numChar = 500)
        colnames(allRes.deseq.hypox.anox)[6] <- "p-value"
        allRes.deseq.hypox.anox$`p-value` <- as.numeric(allRes.deseq.hypox.anox$`p-value`)
        
        # Filter nodes with pvalue greater than 0.05 
        allRes.deseq.hypox.anox <- allRes.deseq.hypox.anox[allRes.deseq.hypox.anox$`p-value` < 0.05,]
        
        
        #Extract names of significant genes in GenTable result. Add to column in allRes.deseq.hypox.anox
        allRes.deseq.hypox.anox$genes <- sapply(allRes.deseq.hypox.anox$GO.ID, function(x)
        {
          genes<-genesInTerm(DESeq_GO_hypox_anox, x)
          genes[[1]][genes[[1]] %in% sig_names_deseq2_hypoxia_vs_anoxia] 
        })
        
        # Create and export enrichment result as a formatted and interactable table
        deseq.hypox.anox.table <- reactable(allRes.deseq.hypox.anox[,1:6], defaultPageSize = 150, theme = pff(centered = FALSE, font_color = "black"), 
                                            wrap = FALSE, bordered = TRUE, compact = TRUE, striped = TRUE, highlight = TRUE, fullWidth = FALSE,
                                            columns = list(Term = colDef(minWidth = 400),
                                                           GO.ID = colDef(cell = pill_buttons(colors = "darkgreen"), minWidth = 120)
                                            )
        ) %>%
          add_title("DESeq2 - Hypoxa vs anoxia significant genes", font_size = 20) %>% 
          add_subtitle("79 of 103 significant genes", font_size = 16, font_style = "italic") 
        
        
        deseq.hypox.anox.table #View table
        
        # Save table as an html file
        save_reactable_test(deseq.hypox.anox.table, "DESeq2 hypoxia vs anoxia table.html")
        
        # investigate how the significant GO terms are distributed over the GO graph
        jpeg(filename = "TopGO subgraph DESeq2 hypoxia vs anoxia.jpg", width = 10, height = 10, units = "in", res = 900)
        par(cex = 0.9)
        showSigOfNodes(DESeq_GO_hypox_anox, score(resultFisher.deseq.hypox.anox), firstSigNodes = 25, useInfo = 'all')
        dev.off()
        
        
        
 
        
#### rrvgo visualization ####
        
    #Make sure you load the T.californicus OrgDb package before doing this. 
      
        
    # Cluster 1 from maSigPro
      #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
        simMatrix.cluster1 <- calculateSimMatrix(allRes.masigpro.cluster1$GO.ID,
                                                 orgdb = Tcalif_orgdb_object,
                                                 ont="BP",
                                                 method="Wang", 
                                                 keytype = "GID")
        
      #Create groupings of reduced terms for easier visualization
        reducedTerms.cluster1 <- reduceSimMatrix(simMatrix.cluster1,
                                                 threshold=0.8,
                                                 orgdb=Tcalif_orgdb_object,
                                                 keytype = "GID")
        
      # Make scatter plot depicting groups and distance between terms  
        jpeg(filename = "rrvgo scatter maSigPro clust 1.jpg", width = 16, height = 11, units = "in", res = 500)
        scatterPlot(simMatrix.cluster1, reducedTerms.cluster1)
        dev.off()
        
      # Make treemap of GO terms clustered under their parent terms
        jpeg(filename = "rrvgo treemap maSigPro clust 1.jpg", width = 16, height = 11, units = "in", res = 500)
        treemapPlot(reducedTerms.cluster1)
        dev.off()
        
    # Cluster 2 from maSigPro
      #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
        simMatrix.cluster2 <- calculateSimMatrix(allRes.masigpro.cluster2$GO.ID,
                                                 orgdb = Tcalif_orgdb_object,
                                                 ont="BP",
                                                 method="Wang", 
                                                 keytype = "GID")
        
      #Create groupings of reduced terms for easier visualization
        reducedTerms.cluster2 <- reduceSimMatrix(simMatrix.cluster2,
                                                 threshold=0.8,
                                                 orgdb=Tcalif_orgdb_object,
                                                 keytype = "GID")
        
      # Make scatter plot depicting groups and distance between terms  
        jpeg(filename = "rrvgo scatter maSigPro clust 2.jpg", width = 16, height = 11, units = "in", res = 500)
        scatterPlot(simMatrix.cluster2, reducedTerms.cluster2)
        dev.off()
        
      # Make treemap of GO terms clustered under their parent terms
        jpeg(filename = "rrvgo treemap maSigPro clust 2.jpg", width = 16, height = 11, units = "in", res = 500)
        treemapPlot(reducedTerms.cluster2)
        dev.off()
        
    # Cluster 3 from maSigPro
      #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
        simMatrix.cluster3 <- calculateSimMatrix(allRes.masigpro.cluster3$GO.ID,
                                                 orgdb = Tcalif_orgdb_object,
                                                 ont="BP",
                                                 method="Wang", 
                                                 keytype = "GID")
        
      #Create groupings of reduced terms for easier visualization
        reducedTerms.cluster3 <- reduceSimMatrix(simMatrix.cluster3,
                                                 threshold=0.8,
                                                 orgdb=Tcalif_orgdb_object,
                                                 keytype = "GID")
        
      # Make scatter plot depicting groups and distance between terms  
        jpeg(filename = "rrvgo scatter maSigPro clust 3.jpg", width = 16, height = 11, units = "in", res = 500)
        scatterPlot(simMatrix.cluster3, reducedTerms.cluster3)
        dev.off()
        
      # Make treemap of GO terms clustered under their parent terms
        jpeg(filename = "rrvgo treemap maSigPro clust 3.jpg", width = 16, height = 11, units = "in", res = 500)
        treemapPlot(reducedTerms.cluster3)
        dev.off()
        
    # Cluster 4 from maSigPro
      #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
        simMatrix.cluster4 <- calculateSimMatrix(allRes.masigpro.cluster4$GO.ID,
                                                 orgdb = Tcalif_orgdb_object,
                                                 ont="BP",
                                                 method="Wang", 
                                                 keytype = "GID")
        
      #Create groupings of reduced terms for easier visualization
        reducedTerms.cluster4 <- reduceSimMatrix(simMatrix.cluster4,
                                                 threshold=0.8,
                                                 orgdb=Tcalif_orgdb_object,
                                                 keytype = "GID")
        
      # Make scatter plot depicting groups and distance between terms  
        jpeg(filename = "rrvgo scatter maSigPro clust 4.jpg", width = 16, height = 11, units = "in", res = 500)
        scatterPlot(simMatrix.cluster4, reducedTerms.cluster4)
        dev.off()
        
      # Make treemap of GO terms clustered under their parent terms
        jpeg(filename = "rrvgo treemap maSigPro clust 4.jpg", width = 16, height = 11, units = "in", res = 500)
        treemapPlot(reducedTerms.cluster4)
        dev.off()
        
        
    # Cluster 5 from maSigPro
      #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
        simMatrix.cluster5 <- calculateSimMatrix(allRes.masigpro.cluster5$GO.ID,
                                                 orgdb = Tcalif_orgdb_object,
                                                 ont="BP",
                                                 method="Wang", 
                                                 keytype = "GID")
        
      #Create groupings of reduced terms for easier visualization
        reducedTerms.cluster5 <- reduceSimMatrix(simMatrix.cluster5,
                                                 threshold=0.8,
                                                 orgdb=Tcalif_orgdb_object,
                                                 keytype = "GID")
        
      # Make scatter plot depicting groups and distance between terms  
        jpeg(filename = "rrvgo scatter maSigPro clust 5.jpg", width = 16, height = 11, units = "in", res = 500)
        scatterPlot(simMatrix.cluster5, reducedTerms.cluster5)
        dev.off()
        
      # Make treemap of GO terms clustered under their parent terms
        jpeg(filename = "rrvgo treemap maSigPro clust 5.jpg", width = 16, height = 11, units = "in", res = 500)
        treemapPlot(reducedTerms.cluster5)
        dev.off()
        
        
    # Cluster 6 from maSigPro
      #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
        simMatrix.cluster6 <- calculateSimMatrix(allRes.masigpro.cluster6$GO.ID,
                                                 orgdb = Tcalif_orgdb_object,
                                                 ont="BP",
                                                 method="Wang", 
                                                 keytype = "GID")
        
      #Create groupings of reduced terms for easier visualization
        reducedTerms.cluster6 <- reduceSimMatrix(simMatrix.cluster6,
                                                 threshold=0.8,
                                                 orgdb=Tcalif_orgdb_object,
                                                 keytype = "GID")
        
      # Make scatter plot depicting groups and distance between terms  
        jpeg(filename = "rrvgo scatter maSigPro clust 6.jpg", width = 16, height = 11, units = "in", res = 500)
        scatterPlot(simMatrix.cluster6, reducedTerms.cluster6)
        dev.off()
        
      # Make treemap of GO terms clustered under their parent terms
        jpeg(filename = "rrvgo treemap maSigPro clust 6.jpg", width = 16, height = 11, units = "in", res = 500)
        treemapPlot(reducedTerms.cluster6)
        dev.off()
        
    # Cluster 7 from maSigPro
      #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
        simMatrix.cluster7 <- calculateSimMatrix(allRes.masigpro.cluster7$GO.ID,
                                                 orgdb = Tcalif_orgdb_object,
                                                 ont="BP",
                                                 method="Wang", 
                                                 keytype = "GID")
        
      #Create groupings of reduced terms for easier visualization
        reducedTerms.cluster7 <- reduceSimMatrix(simMatrix.cluster7,
                                                 threshold=0.8,
                                                 orgdb=Tcalif_orgdb_object,
                                                 keytype = "GID")
        
      # Make scatter plot depicting groups and distance between terms  
        jpeg(filename = "rrvgo scatter maSigPro clust 7.jpg", width = 16, height = 11, units = "in", res = 500)
        scatterPlot(simMatrix.cluster7, reducedTerms.cluster7)
        dev.off()
        
      # Make treemap of GO terms clustered under their parent terms
        jpeg(filename = "rrvgo treemap maSigPro clust 7.jpg", width = 16, height = 11, units = "in", res = 500)
        treemapPlot(reducedTerms.cluster7)
        dev.off()
        
    # Cluster 8 from maSigPro
      #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
        simMatrix.cluster8 <- calculateSimMatrix(allRes.masigpro.cluster8$GO.ID,
                                                 orgdb = Tcalif_orgdb_object,
                                                 ont="BP",
                                                 method="Wang", 
                                                 keytype = "GID")
        
      #Create groupings of reduced terms for easier visualization
        reducedTerms.cluster8 <- reduceSimMatrix(simMatrix.cluster8,
                                                 threshold=0.8,
                                                 orgdb=Tcalif_orgdb_object,
                                                 keytype = "GID")
        
      # Make scatter plot depicting groups and distance between terms  
        jpeg(filename = "rrvgo scatter maSigPro clust 8.jpg", width = 16, height = 11, units = "in", res = 500)
        scatterPlot(simMatrix.cluster8, reducedTerms.cluster8)
        dev.off()
        
      # Make treemap of GO terms clustered under their parent terms
        jpeg(filename = "rrvgo treemap maSigPro clust 8.jpg", width = 16, height = 11, units = "in", res = 500)
        treemapPlot(reducedTerms.cluster8)
        dev.off()
        
    # Cluster 9 from maSigPro
      #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
        simMatrix.cluster9 <- calculateSimMatrix(allRes.masigpro.cluster9$GO.ID,
                                                 orgdb = Tcalif_orgdb_object,
                                                 ont="BP",
                                                 method="Wang", 
                                                 keytype = "GID")
        
      #Create groupings of reduced terms for easier visualization
        reducedTerms.cluster9 <- reduceSimMatrix(simMatrix.cluster9,
                                                 threshold=0.8,
                                                 orgdb=Tcalif_orgdb_object,
                                                 keytype = "GID")
        
      # Make scatter plot depicting groups and distance between terms  
        jpeg(filename = "rrvgo scatter maSigPro clust 9.jpg", width = 16, height = 11, units = "in", res = 500)
        scatterPlot(simMatrix.cluster9, reducedTerms.cluster9)
        dev.off()
        
      # Make treemap of GO terms clustered under their parent terms
        jpeg(filename = "rrvgo treemap maSigPro clust 9.jpg", width = 16, height = 11, units = "in", res = 500)
        treemapPlot(reducedTerms.cluster9)
        dev.off()
        
        
        
    # Pcrit vs control from deseq2
      #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
        simMatrix.pcrit <- calculateSimMatrix(allRes.deseq.pcrit.cont$GO.ID,
                                                 orgdb = Tcalif_orgdb_object,
                                                 ont="BP",
                                                 method="Wang", 
                                                 keytype = "GID")
        
      #Create groupings of reduced terms for easier visualization
        reducedTerms.pcrit <- reduceSimMatrix(simMatrix.pcrit,
                                                 threshold=0.8,
                                                 orgdb=Tcalif_orgdb_object,
                                                 keytype = "GID")
      
      # Make scatter plot depicting groups and distance between terms  
        jpeg(filename = "rrvgo deseq pcrit vs control scatter.jpg", width = 16, height = 11, units = "in", res = 500)
        scatterPlot(simMatrix.pcrit, reducedTerms.pcrit)
        dev.off()
        
      # Make treemap of GO terms clustered under their parent terms
        jpeg(filename = "rrvgo deseq pcrit vs control treemap.jpg", width = 16, height = 11, units = "in", res = 500)
        treemapPlot(reducedTerms.pcrit)
        dev.off()
        
        
     # hypoxia vs control from deseq2
      #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
        simMatrix.hypox <- calculateSimMatrix(allRes.deseq.hypox.cont$GO.ID,
                                              orgdb = Tcalif_orgdb_object,
                                              ont="BP",
                                              method="Wang", 
                                              keytype = "GID")
        
      #Create groupings of reduced terms for easier visualization
        reducedTerms.hypox <- reduceSimMatrix(simMatrix.hypox,
                                              threshold=0.8,
                                              orgdb=Tcalif_orgdb_object,
                                              keytype = "GID")
        
      # Make scatter plot depicting groups and distance between terms  
        jpeg(filename = "rrvgo deseq hypox vs control scatter.jpg", width = 16, height = 11, units = "in", res = 500)
        scatterPlot(simMatrix.hypox, reducedTerms.hypox)
        dev.off()
        
      # Make treemap of GO terms clustered under their parent terms
        jpeg(filename = "rrvgo deseq hypox vs control treemap.jpg", width = 16, height = 11, units = "in", res = 500)
        treemapPlot(reducedTerms.hypox)
        dev.off()
        
    # Anoxia vs control from deseq2
      #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
        simMatrix.anox <- calculateSimMatrix(allRes.deseq.anox.cont$GO.ID,
                                              orgdb = Tcalif_orgdb_object,
                                              ont="BP",
                                              method="Wang", 
                                              keytype = "GID")
        
      #Create groupings of reduced terms for easier visualization
        reducedTerms.anox <- reduceSimMatrix(simMatrix.anox,
                                              threshold=0.8,
                                              orgdb=Tcalif_orgdb_object,
                                              keytype = "GID")
        
      # Make scatter plot depicting groups and distance between terms  
        jpeg(filename = "rrvgo deseq anox vs control scatter.jpg", width = 16, height = 11, units = "in", res = 500)
        scatterPlot(simMatrix.anox, reducedTerms.anox)
        dev.off()
        
      # Make treemap of GO terms clustered under their parent terms
        jpeg(filename = "rrvgo deseq anox vs control treemap.jpg", width = 16, height = 11, units = "in", res = 500)
        treemapPlot(reducedTerms.anox)
        dev.off()
        
    # Recovery vs control from deseq2
      #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
        simMatrix.recov <- calculateSimMatrix(allRes.deseq.recov.cont$GO.ID,
                                             orgdb = Tcalif_orgdb_object,
                                             ont="BP",
                                             method="Wang", 
                                             keytype = "GID")
        
      #Create groupings of reduced terms for easier visualization
        reducedTerms.recov <- reduceSimMatrix(simMatrix.recov,
                                             threshold=0.8,
                                             orgdb=Tcalif_orgdb_object,
                                             keytype = "GID")
        
       # Make scatter plot depicting groups and distance between terms  
        jpeg(filename = "rrvgo deseq recov vs control scatter.jpg", width = 16, height = 11, units = "in", res = 500)
        scatterPlot(simMatrix.recov, reducedTerms.recov)
        dev.off()
        
      # Make treemap of GO terms clustered under their parent terms
        jpeg(filename = "rrvgo deseq recov vs control treemap.jpg", width = 16, height = 11, units = "in", res = 500)
        treemapPlot(reducedTerms.recov)
        dev.off()
        
    # Recovery vs anoxia from deseq2
      #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
        simMatrix.anox.recov <- calculateSimMatrix(allRes.deseq.anox.recov$GO.ID,
                                              orgdb = Tcalif_orgdb_object,
                                              ont="BP",
                                              method="Wang", 
                                              keytype = "GID")
        
      #Create groupings of reduced terms for easier visualization
        reducedTerms.anox.recov <- reduceSimMatrix(simMatrix.anox.recov,
                                              threshold=0.8,
                                              orgdb=Tcalif_orgdb_object,
                                              keytype = "GID")
        
      # Make scatter plot depicting groups and distance between terms  
        jpeg(filename = "rrvgo deseq recov vs anoxia scatter.jpg", width = 16, height = 11, units = "in", res = 500)
        scatterPlot(simMatrix.anox.recov, reducedTerms.anox.recov)
        dev.off()
        
      # Make treemap of GO terms clustered under their parent terms
        jpeg(filename = "rrvgo deseq recov vs anoxia treemap.jpg", width = 16, height = 11, units = "in", res = 500)
        treemapPlot(reducedTerms.anox.recov)
        dev.off()
        
    # Hypoxia vs anoxia from deseq2
      #Create sim matrix using method Wang since it is most recent and keytype GID since we dont have ENTREZIDs
        simMatrix.hypox.anox <- calculateSimMatrix(allRes.deseq.hypox.anox$GO.ID,
                                                   orgdb = Tcalif_orgdb_object,
                                                   ont="BP",
                                                   method="Wang", 
                                                   keytype = "GID")
        
      #Create groupings of reduced terms for easier visualization
        reducedTerms.hypox.anox <- reduceSimMatrix(simMatrix.hypox.anox,
                                                   threshold=0.8,
                                                   orgdb=Tcalif_orgdb_object,
                                                   keytype = "GID")
        
      # Make scatter plot depicting groups and distance between terms  
        jpeg(filename = "rrvgo deseq hypox vs anoxia scatter.jpg", width = 16, height = 11, units = "in", res = 500)
        scatterPlot(simMatrix.hypox.anox, reducedTerms.hypox.anox)
        dev.off()
        
      # Make treemap of GO terms clustered under their parent terms
        jpeg(filename = "rrvgo deseq hypox vs anoxia treemap.jpg", width = 16, height = 11, units = "in", res = 500)
        treemapPlot(reducedTerms.hypox.anox)
        dev.off()
        
    
        
##################### synthesis ###################
        
      # Read in output from DESeq pairwise comparisons
        Sigs.anox.vs.control <- read.csv(file = "Sig_genes_DESeq2_anox_vs_cont.csv", header = TRUE)
        
        Sigs.anox.vs.recov <- read.csv(file = "Sig_genes_DESeq2_anox_vs_recov.csv", header = TRUE)
        
        Sigs.pcrit.vs.cont <- read.csv(file = "Sig_genes_DESeq2_pcrit_vs_control.csv", header = TRUE)
        
        Sigs.hypox.vs.cont <- read.csv(file = "Sig_genes_DESeq2_hypoxia_vs_control.csv", header = TRUE)
        
        Sigs.recov.vs.cont <- read.csv(file = "Sig_genes_DESeq2_recovery_vs_control.csv", header = TRUE)
        
        Sigs.hypox.vs.anox <- read.csv(file = "Sig_genes_DESeq2_hypoxia_vs_anoxia.csv", header = TRUE)
        
      # Order data frames by gene ID to combine easily with cbind
        Sigs.hypox.vs.cont <- Sigs.hypox.vs.cont[order(Sigs.hypox.vs.cont$X), ]
        Sigs.pcrit.vs.cont <- Sigs.pcrit.vs.cont[order(Sigs.pcrit.vs.cont$X), ]
        Sigs.anox.vs.control <- Sigs.anox.vs.control[order(Sigs.anox.vs.control$X), ]
        Sigs.recov.vs.cont <- Sigs.recov.vs.cont[order(Sigs.recov.vs.cont$X), ]
        Sigs.anox.vs.recov <- Sigs.anox.vs.recov[order(Sigs.anox.vs.recov$X), ]
        Sigs.hypox.vs.anox <- Sigs.hypox.vs.anox[order(Sigs.hypox.vs.anox$X), ]
        
      ## RUN THIS BEFORE OVERWRITING FRAMES WITH JUST SIGNIFICANT GENES BELOW 
      #Combine gene names and select columns from above sets to make a list of genes with the fold changes relative to controls for plotting
      #Leave off anoxia vs recov for now. Just want to plot fold change relative to control group at this point 
        Sigs.folds <- cbind(Sigs.hypox.vs.cont[c(1,3)], Sigs.pcrit.vs.cont[3], Sigs.anox.vs.control[3], Sigs.recov.vs.cont[3]) 
      #Rename columns 
        names(Sigs.folds) <- c("Gene", "fold_hypoxia", "fold_pcrit", "fold_anoxia", "fold_recovery")
      #Reset row names
        row.names(Sigs.folds) <- NULL
        
        
        
      # Get just significant genes
        Sigs.anox.vs.control <- Sigs.anox.vs.control[Sigs.anox.vs.control$padj < 0.1 & !is.na(Sigs.anox.vs.control$padj), ]
        names(Sigs.anox.vs.control)[1] <- "Gene" #Rename column
        
        Sigs.anox.vs.recov <- Sigs.anox.vs.recov[Sigs.anox.vs.recov$padj < 0.1 & !is.na(Sigs.anox.vs.recov$padj), ]
        names(Sigs.anox.vs.recov)[1] <- "Gene" # Rename column
        
        Sigs.pcrit.vs.cont <- Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$padj < 0.1 & !is.na(Sigs.pcrit.vs.cont$padj), ]
        names(Sigs.pcrit.vs.cont)[1] <- "Gene" # Rename column
        
        Sigs.hypox.vs.cont <- Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$padj < 0.1 & !is.na(Sigs.hypox.vs.cont$padj), ]
        names(Sigs.hypox.vs.cont)[1] <- "Gene" # Rename column
        
        Sigs.recov.vs.cont <- Sigs.recov.vs.cont[Sigs.recov.vs.cont$padj < 0.1 & !is.na(Sigs.recov.vs.cont$padj), ]
        names(Sigs.recov.vs.cont)[1] <- "Gene" # Rename column
        
        Sigs.hypox.vs.anox <- Sigs.hypox.vs.anox[Sigs.hypox.vs.anox$padj < 0.1 & !is.na(Sigs.hypox.vs.anox$padj), ]
        names(Sigs.hypox.vs.anox)[1] <- "Gene" # Rename column
        
        
      # Create overlapping data frame with all the significant genes from maSigPro, ImpulseDE2, and EBSeq-HMM
        combined_sigs <- merge(sigs.masigpro.pvalues.and.betas[c(1,2,12)], Sigs.hypox.vs.cont[c(1,3,7)],  by = "Gene", all = TRUE)
        names(combined_sigs)[2:5] <- c("padj_maSigPro" , "clusters_maSigPro", "Fold_change_h_vs_c", "padj_h_vs_c")
        
        combined_sigs <- merge(combined_sigs, Sigs.pcrit.vs.cont[c(1,3,7)], by = "Gene", all = TRUE)
        names(combined_sigs)[6:7] <- c("Fold_change_p_vs_c" , "padj_p_vs_c")
        
        combined_sigs <- merge(combined_sigs, Sigs.anox.vs.control[c(1,3,7)], by = "Gene", all = TRUE)
        names(combined_sigs)[8:9] <- c("Fold_change_a_vs_c" , "padj_a_vs_c")
        
        combined_sigs <- merge(combined_sigs, Sigs.recov.vs.cont[c(1,3,7)], by = "Gene", all = TRUE)
        names(combined_sigs)[10:11] <- c("Fold_change_r_vs_c" , "padj_r_vs_c")
        
        combined_sigs <- merge(combined_sigs, Sigs.anox.vs.recov[c(1,3,7)], by = "Gene", all = TRUE)
        names(combined_sigs)[12:13] <- c("Fold_change_a_vs_r" , "padj_a_vs_r")
        
        combined_sigs <- merge(combined_sigs, Sigs.hypox.vs.anox[c(1,3,7)], by = "Gene", all = TRUE)
        names(combined_sigs)[14:15] <- c("Fold_change_h_vs_a" , "padj_h_vs_a")
        
        
      #If you want to get unique gene names for DESeq comparisons
        # Genes unique to hypoxia vs control
          combined_sigs[is.na(combined_sigs$padj_p_vs_c) & is.na(combined_sigs$padj_a_vs_c)
                        & is.na(combined_sigs$padj_r_vs_c) & is.na(combined_sigs$padj_a_vs_r)
                        & !is.na(combined_sigs$padj_h_vs_c), 1] 
          
        # Genes unique to pcrit vs control
          combined_sigs[is.na(combined_sigs$padj_h_vs_c) & is.na(combined_sigs$padj_a_vs_c)
                        & is.na(combined_sigs$padj_r_vs_c) & is.na(combined_sigs$padj_a_vs_r)
                        & !is.na(combined_sigs$padj_p_vs_c), 1] 
          
        # Genes unique to anoxia vs control
          combined_sigs[is.na(combined_sigs$padj_p_vs_c) & is.na(combined_sigs$padj_h_vs_c)
                        & is.na(combined_sigs$padj_r_vs_c) & is.na(combined_sigs$padj_a_vs_r)
                        & !is.na(combined_sigs$padj_a_vs_c), 1] 
          
        # Genes unique to recovery vs control
          combined_sigs[is.na(combined_sigs$padj_p_vs_c) & is.na(combined_sigs$padj_a_vs_c)
                        & is.na(combined_sigs$padj_h_vs_c) & is.na(combined_sigs$padj_a_vs_r)
                        & !is.na(combined_sigs$padj_r_vs_c), 1] 
          
        # Genes unique to anoxia vs recovery
          combined_sigs[is.na(combined_sigs$padj_p_vs_c) & is.na(combined_sigs$padj_a_vs_c)
                        & is.na(combined_sigs$padj_r_vs_c) & is.na(combined_sigs$padj_h_vs_c)
                        & !is.na(combined_sigs$padj_a_vs_r), 1]
          
        # Genes shared by anoxia vs recovery and hypoxia vs control
          combined_sigs[is.na(combined_sigs$padj_p_vs_c) & is.na(combined_sigs$padj_a_vs_c)
                        & !is.na(combined_sigs$padj_r_vs_c) & !is.na(combined_sigs$padj_h_vs_c)
                        & !is.na(combined_sigs$padj_a_vs_r), 1] 
        
        
        
        
      # Make list of genes in each group for venn diagram of time course programs
        # Deseq comparison
          gene.venn.list.deseq <- list(
            Mild_vs_normoxia = combined_sigs[!is.na(combined_sigs$padj_h_vs_c), 1],
            Anoxia_vs_normoxia = combined_sigs[!is.na(combined_sigs$padj_a_vs_c), 1],
            Pcrit_vs_normoxia = combined_sigs[!is.na(combined_sigs$padj_p_vs_c), 1],
            Recovery_vs_normoxia = combined_sigs[!is.na(combined_sigs$padj_r_vs_c), 1],
            Anoxia_vs_recovery = combined_sigs[!is.na(combined_sigs$padj_a_vs_r), 1],
            Hypoxia_vs_anoxia = combined_sigs[!is.na(combined_sigs$padj_h_vs_a), 1]
            
          )
          
  
      ## Get the total number of unique significant genes in the six comparisons above
         length(unique(c(gene.venn.list.deseq$Mild_vs_normoxia, 
                         gene.venn.list.deseq$Anoxia_vs_normoxia, 
                         gene.venn.list.deseq$Pcrit_vs_normoxia,
                         gene.venn.list.deseq$Recovery_vs_normoxia, 
                         gene.venn.list.deseq$Anoxia_vs_recovery,
                         gene.venn.list.deseq$Hypoxia_vs_anoxia)))
         
        
        
       # Draw venn diagram with ggvenn package for time course packages
        
        jpeg(filename = "Venn Diagram of Sig Genes from deseq comparisons.jpg", width = 10, height = 8, units = "in", res = 300)
        ggVennDiagram(gene.venn.list.deseq[c(1,3,2,5)], label_percent_digit = 1, label_size = 5) +
          scale_x_continuous(expand = expansion(mult = .2)) +
          scale_fill_distiller(name = "Counts:", palette = "Blues", direction = 1) + 
          theme(legend.position = "bottom", legend.direction = "horizontal")
        dev.off()
        
        jpeg(filename = "Venn Diagram of Sig Genes from deseq CONTROL comparisons.jpg", width = 10, height = 8, units = "in", res = 300)
        ggVennDiagram(gene.venn.list.deseq[c(1,3,2,4)], label_percent_digit = 1, label_size = 5) +
          scale_x_continuous(expand = expansion(mult = .2)) +
          scale_fill_distiller(name = "Counts:", palette = "Greens", direction = 1) + 
          theme(legend.position = "bottom", legend.direction = "horizontal")
        dev.off()
        
        
      # Draw venn diagram of all six groups using venn package 
        #Make color palette for plotting
        br_pal <- moma.colors("Klein")
        show_col(br_pal)
        
        #Draw diagram
        jpeg(filename = "Venn Diagram All DESeq2 comparisons.jpg", width = 10, height = 10, units = "in", res = 300)
        venn(gene.venn.list.deseq, ilabels = "counts", zcolor = br_pal[c(2,8,5,4,10,1)], 
             box = FALSE, opacity = 0.4, ilcs = 1.1, sncs = 1)
        dev.off()
        
    # Add functional annotation to combined sigs table
          SD_annot_2017 <- read.csv(file = "SDv2.2_BlastTable.csv", header = TRUE)
        
        # Adjust sequence name to match my naming scheme (remove the -PA and the IF from TCALIF)
          SD_annot_2017$SeqName <- str_split_i(SD_annot_2017$SeqName, "-", 1)
          
        # Get rid of extra columns because of my OCD
          SD_annot_2017 <- SD_annot_2017[c(1:2)]
          
        # Rename seqname to have the same column header as combined_sigs frame
          names(SD_annot_2017)[1] <- "Gene"
          
        # Merge the annotation with our combined_sigs data frame, keeping only the annotations that match to our entries
          combined_sigs <- merge(combined_sigs, SD_annot_2017, by = "Gene", all.x = TRUE)
        
        #Move annotation columns to the front with the gene names  
          combined_sigs %>% relocate(Description, .after=Gene) -> combined_sigs
          
        # rename rows to have gene name for indexing and reference later
          row.names(combined_sigs) <- combined_sigs$Gene
          
        # Export master significant genes list as excel sheet
          write.csv(combined_sigs, file = "Master Significant Genes list with Annotation.csv", row.names = FALSE)
          
          
      # Merge the annotation with our all.masigpro data frame, keeping only the annotations that match to our entries
          #Follows most of the same steps as above with some modifications
          combined_all <- merge(all.masigpro, SD_annot_2017, by = "Gene", all.x = TRUE) # Merge all genes with annotation info
          combined_all2 <- merge(combined_sigs, combined_all[,c(1,3,4)], by = "Gene", all = T) #Merge data frames (padj and description columns with retain with .x and .y notation)
          combined_all2$Description.x <- combined_all2$Description.y #Copy annotation info from all genes into first description column
          combined_all2 <- combined_all2[,-c(18)] #Remove second description column
          colnames(combined_all2)[2] <- "Description" #Get rid of .x at the end of column name
          # mutate .x pvalue column to fill in NA's with info from .y pvalue column using coalesce function
          combined_all2 %>% 
            mutate(padj_maSigPro.x = coalesce(padj_maSigPro.x, padj_maSigPro.y)) -> combined_all2 
          combined_all2 <- combined_all2[,-c(17)]
          rm(combined_all) # Remove intermediate data frame
          #Export as csv
          write.csv(combined_all2, file = "All_genes_pvalues_maSigPro.csv", row.names = FALSE)
          
    #Exploring the expression pattern of select genes    
          
            #Make melted data frame for plotting
            ddsNorm.melt <- as.data.frame(ddsNorm)
            ddsNorm.melt$Gene <- row.names(ddsNorm.melt) #Send gene names to a column to use as id for melting
            rownames(ddsNorm.melt) <- NULL #Make row.names default numbers
            ddsNorm.melt <- melt(ddsNorm.melt, id = c("Gene")) #Melt the data frame lengthwise
            names(ddsNorm.melt)[2:3] <- c("sampleID", "expression") #Rename the melted columns
            
            
            ddsNorm.melt[c('time', 'replicate')] <- str_split_fixed(ddsNorm.melt$sampleID, '_', 2)
            #Rename and relevel group variable to put in the desired order for plotting
            ddsNorm.melt$time<- factor(ddsNorm.melt$time,
                                             levels = c("C", "35","05", "A", "R"))
            ddsNorm.melt %>%
              mutate(time = fct_recode(time, "H" = "35",
                                            "P" = "05")) -> ddsNorm.melt
  
            ddsNorm.melt$gene.rep <- paste(ddsNorm.melt$Gene, ddsNorm.melt$replicate)
            
          # Make averaged data frame and export  
            expression.averaged <- ddsNorm.melt |>
              group_by(time, Gene) |> 
              summarise(mean_exp = mean(expression) |> round(4),
                        sd_exp = sd(expression) |> round(4),
                        n = n(),
                        range = paste(range(expression) |> round(4), collapse = ' - ')
              ) |> 
              mutate(se_exp = sd_exp/sqrt(n)) 
            expression.averaged$se_exp <- round(expression.averaged$se_exp, 4) #Round se column to 2 digits
            
            write.csv(expression.averaged, file = "Counts_averaged_by_time_point_ALL_GENES.csv", row.names = FALSE)
         
          #Generate wider version wwhere we average the 6 replicates for individual gene heatmaps later
            expression.averaged.wide <- pivot_wider(expression.averaged[,c(1:3)], names_from=time, values_from=c(mean_exp))
            
      # Probing individual genes and GO terms from enrichment analysis   
        
        #Plot individual genes
        # Generate expression plot   
          jpeg(filename = "TCAL_17457 Gfpt1.jpg", width = 5, height = 4, units = "in", res = 300) 
            expression.averaged |>
              filter(Gene == "TCAL_09981") |>
              ggplot(aes(x=time, y=mean_exp)) +
              geom_line(aes(group=Gene), linewidth=2, alpha=1)+
              geom_point(size = 5)+
              geom_errorbar(aes(ymax = mean_exp + se_exp, ymin = mean_exp - se_exp), width=0.3) +
              scale_y_continuous(name="Normalized Expression")+
              scale_x_discrete(name="Time Points")+
              theme(axis.text.x = element_text(size = 11),
                    axis.title = element_text(size = 15),
                    axis.text.y = element_text(size=10),
                    legend.position = "bottom", legend.direction = "horizontal",
                    title = element_text(size = 13))+
              ggtitle("TCAL_17457 Gfpt1")
            dev.off()
          
      # Pull gene names from GO terms in enrichment test
        # Inspect html tables from GO Enrichment output for terms of interest and enter in quotes at the end of lines below
          
          #All sig genes
          allRes.masigpro$genes$`GO:0016137`
          
          #Clusters
          int.1 <-allRes.masigpro.cluster1$genes$`GO:`
          int.2 <-allRes.masigpro.cluster2$genes$`GO:`
          int.3 <-allRes.masigpro.cluster3$genes$`GO:0045721`
          int.4 <-allRes.masigpro.cluster4$genes$`GO:0042546`
          int.5 <-allRes.masigpro.cluster5$genes$`GO:0043482`
          int.6 <- allRes.masigpro.cluster6$genes$`GO:0006030`
          int.7 <- allRes.masigpro.cluster7$genes$`GO:0006030`
          int.8 <- allRes.masigpro.cluster8$genes$`GO:0006030`
          int.9 <- allRes.masigpro.cluster9$genes$`GO:0006000`
          
          #Pairwise DESeq2 comparisons
          int.hypox.cont <- allRes.deseq.hypox.cont$genes$`GO:1901700`
          int.pcrit.cont <- allRes.deseq.pcrit.cont$genes$`GO:0015767`
          int.anox.cont <- allRes.deseq.anox.cont$genes$`GO:0001678`
          int.anox.recov <- allRes.deseq.anox.recov$genes$`GO:0043981`
          int.recov.cont <- allRes.deseq.recov.cont$genes$`GO:0001523`
          int.hypox.anox <- allRes.deseq.hypox.anox$genes$`GO:0001523`
          
          
      # Pull out parts of the combined_sigs master results data frame that correspond to go terms of interest saved from above
          combined_sigs[int.1, c(1:2)]
          combined_sigs[int.2, c(1:2)]
          combined_sigs[int.3, c(1:2)]
          combined_sigs[int.4, c(1:2)]
          combined_sigs[int.5, c(1:2)]
          combined_sigs[int.6, c(1:2)]
          combined_sigs[int.7, c(1:2)]
          combined_sigs[int.8, c(1:2)]
          combined_sigs[int.9, c(1:2)]
        
        # Pairwise DESeq2 comparisons  
          combined_sigs[int.hypox.cont, c(1:2)]
          combined_sigs[int.pcrit.cont, c(1:2)]
          combined_sigs[int.anox.cont, c(1:2)]
          combined_sigs[int.anox.recov, c(1:2)]
          combined_sigs[int.recov.cont, c(1:2)]
          combined_sigs[int.hypox.anox, c(1:2)]
          
          
####   Functional analysis ####
          
  ### Generate interactive line graphs to explore expression patterns across similar gene groups
          
        #Make melted data frame for plotting
          ddsNormScaled.melt <- as.data.frame(ddsNormScaled)
          ddsNormScaled.melt$Gene <- row.names(ddsNormScaled.melt) #Send gene names to a column to use as id for melting
          rownames(ddsNormScaled.melt) <- NULL #Make row.names default numbers
          ddsNormScaled.melt <- melt(ddsNormScaled.melt, id = c("Gene")) #Melt the data frame lengthwise
          names(ddsNormScaled.melt)[2:3] <- c("sampleID", "zscore") #Rename the melted columns
          ddsNormScaled.melt[c('time', 'replicate')] <- str_split_fixed(ddsNormScaled.melt$sampleID, '_', 2) #Split info into new columns
      
          #Rename and relevel group variable to put in the desired order for plotting
          ddsNormScaled.melt$time<- factor(ddsNormScaled.melt$time,
                                     levels = c("C", "35","05", "A", "R"))
          ddsNormScaled.melt %>% 
            mutate(time = fct_recode(time,"H" = "35",
                                         "P" = "05")) -> ddsNormScaled.melt
          
          ddsNormScaled.melt$gene.rep <- paste(ddsNormScaled.melt$Gene, ddsNormScaled.melt$replicate)
          
        # Make averaged data frame of gene expression z-score scalings and export  
          zscores.averaged <- ddsNormScaled.melt |>
            group_by(time, Gene) |> 
            summarise(mean_zscore = mean(zscore) |> round(4),
                      sd_zscore = sd(zscore) |> round(4),
                      n = n(),
                      range = paste(range(zscore) |> round(4), collapse = ' - ')
            ) |> 
            mutate(se_zscore = sd_zscore/sqrt(n)) 
          zscores.averaged$se_zscore <- round(zscores.averaged$se_zscore, 4) #Round se column to 2 digits
          
          write.csv(zscores.averaged, file = "zscores_averaged_by_time_point_ALL_GENES.csv", row.names = FALSE)
          
          
      
          
          
    ### Pull out select genes based on desired groupings and plot zscores
          
      # Cuticle and chitin related genes
        
        # Read in data and reformat gene names to remove PA notation if present
          cuticles <- read.csv(file = "Chitin_and_cuticle_processes_gene_list_SDv2.2_noPA.csv", header = TRUE)
          names(cuticles)[1] <- "Gene"
          cuticles$tip <- paste0(cuticles$Gene,": ", cuticles$Description)
          
        # Add coloring instructions for lines based on whether gene is in significant list
          cuticles$sig <- NA
          for(i in 1:nrow(cuticles)) {
          if (cuticles[i,1] %in% combined_sigs$Gene) {
            cuticles[i,7] <- "sig"
          } else {
            cuticles[i,7] <- "nonsig"
          }
          }
          
        # Make line plot   
          cuticle.plot <- merge(cuticles, zscores.averaged,  by = "Gene", all = FALSE) |>
            ggplot(aes(x=time, y=mean_zscore, group = Gene)) +
            facet_wrap(~Process,  nrow = 2, strip.position = "bottom")+
            geom_line(aes(colour = sig),linewidth=1)+
            geom_point_interactive(aes(data_id = Gene, tooltip = tip), color = "grey20", size = 3)+
            #geom_errorbar(aes(ymax = mean_zscore + se_zscore, ymin = mean_zscore - se_zscore), width=0.3) +
            scale_y_continuous(name="Z-scores", breaks = seq(-2,2, by =0.4))+
            scale_x_discrete(name="Time Points")+
            theme(axis.text.x = element_text(size = 9),
                  axis.title = element_text(size = 11),
                  axis.text.y = element_text(size=9),
                  legend.position = "bottom", legend.direction = "horizontal",
                  title = element_text(size= 10))+
            scale_color_manual(name = "Significance:", values = c("grey40", "dodgerblue"))+
            ggtitle("Cuticle genes")
          cuticle.plot
          
        # Export as jpeg
          jpeg(filename = "cuticle genes time series line plots.jpg", width = 12, height = 10, units = "in", res = 300)
          cuticle.plot
          dev.off()
       
        # Generate interactive figure   
          girafe(
            ggobj = cuticle.plot,
            options = list(opts_hover_inv(css = "opacity:0.1;stroke:yellow"), 
                           opts_tooltip(offx =0, offy = -50, css = 'font-size:larger;background-color:#000000;color:white')
            ),
            height_svg = 5,
            width_svg = 10
          )
          
          
      # Glycolysis related
          
        # Read in data and reformat gene names to remove PA notation if present
          glycolysis <- read.csv(file = "Glycolysis_and_related_processes_gene_list_SDv2.2_noPA.csv", header = TRUE)
          names(glycolysis)[1] <- "Gene"
          glycolysis$tip <- paste0(glycolysis$Gene,": ", glycolysis$Description)
          
        
          
        # Add coloring instructions for lines based on whether gene is in significant list
          glycolysis$sig <- NA
          for(i in 1:nrow(glycolysis)) {
            if (glycolysis[i,1] %in% combined_sigs$Gene) {
              glycolysis[i,8] <- "sig"
            } else {
              glycolysis[i,8] <- "nonsig"
            }
          }
          
        # Make line plot   
          glycolysis.plot <- merge(glycolysis, zscores.averaged,  by = "Gene", all = FALSE) |>
            ggplot(aes(x=time, y=mean_zscore, group = Gene)) +
            facet_wrap(~Process,  nrow = 2, strip.position = "bottom")+
            geom_line(aes(colour = sig),linewidth=1)+
            geom_point_interactive(aes(data_id = Gene, tooltip = tip), color="grey20", size = 3)+
            #geom_errorbar(aes(ymax = mean_zscore + se_zscore, ymin = mean_zscore - se_zscore), width=0.3) +
            scale_y_continuous(name="Z-scores", breaks = seq(-2,2, by =0.4))+
            scale_x_discrete(name="Time Points")+
            theme(axis.text.x = element_text(size = 9),
                  axis.title = element_text(size = 11),
                  axis.text.y = element_text(size=9),
                  legend.position = "bottom", legend.direction = "horizontal",
                  title = element_text(size= 10))+
            scale_color_manual(name = "Significance:", values = c("grey40", "dodgerblue"))+
            ggtitle("Glycolysis genes")
          glycolysis.plot
          
        # Export as jpeg
          jpeg(filename = "Glycolysis genes time series line plots.jpg", width = 12, height = 10, units = "in", res = 300)
          glycolysis.plot
          dev.off()
          
        # Generate interactive figure   
          girafe(
            ggobj = glycolysis.plot,
            options = list(opts_hover_inv(css = "opacity:0.1;stroke:yellow"), 
                           opts_tooltip(offx = 0, offy = -50, css = 'font-size:larger;background-color:#000000;color:white')
            ),
            height_svg = 5,
            width_svg = 10
          )
          
          
      # Antioxidant related genes
          
        # Read in data and reformat gene names to remove PA notation if present
          antioxidants <- read.csv(file = "Antioxidant_genes_noPA.csv", header = TRUE)
          names(antioxidants)[1] <- "Gene"
          antioxidants$tip <- paste0(antioxidants$Gene,": ", antioxidants$Description)
          antioxidants$Gene <-str_split_i(antioxidants$Gene, "-", 1)
          
        # Add coloring instructions for lines based on whether gene is in significant list
          antioxidants$sig <- NA
          for(i in 1:nrow(antioxidants)) {
            if (antioxidants[i,1] %in% combined_sigs$Gene) {
              antioxidants[i,7] <- "sig"
            } else {
              antioxidants[i,7] <- "nonsig"
            }
          }
          
        # Make line plot   
          antioxidants.plot <- merge(antioxidants, zscores.averaged,  by = "Gene", all = FALSE) |>
            ggplot(aes(x=time, y=mean_zscore, group = Gene)) +
            facet_wrap(~Process,  nrow = 2, strip.position = "bottom")+
            geom_line(aes(colour = sig),linewidth=1)+
            geom_point_interactive(aes(data_id = Gene, tooltip = tip), color = "grey20", size = 3)+
            scale_y_continuous(name="Z-scores", breaks = seq(-2,2, by =0.4))+
            scale_x_discrete(name="Time Points")+
            theme(axis.text.x = element_text(size = 9),
                  axis.title = element_text(size = 11),
                   axis.text.y = element_text(size=9),
                   legend.position = "bottom", legend.direction = "horizontal",
                   title = element_text(size= 10))+
            scale_color_manual(name = "Significance:", values = c("grey40", "dodgerblue"))+
             ggtitle("Antioxidant genes")
             antioxidants.plot
                                   
        # Export as jpeg
         jpeg(filename = "antioxidants genes time series line plots.jpg", width = 12, height = 10, units = "in", res = 300)
         antioxidants.plot
         dev.off()
         
        # Generate interactive figure   
         girafe(
           ggobj = antioxidants.plot,
           options = list(opts_hover_inv(css = "opacity:0.1;stroke:yellow"), 
                          opts_tooltip(offx = 0, offy = -50, css = 'font-size:larger;background-color:#000000;color:white')
           ),
           height_svg = 5,
           width_svg = 10
         )
         
         
      # Carotenoid processing genes
         
        # Read in data and reformat gene names to remove PA notation if present
         carotenoids <- read.csv(file = "Carotenoid_processes_gene_list_SDv2.2_noPA_final.csv", header = TRUE)
         names(carotenoids)[1] <- "Gene"
         carotenoids$tip <- paste0(carotenoids$Gene,": ", carotenoids$Description)
         carotenoids$Gene <-str_split_i(carotenoids$Gene, "-", 1)
         
        # Add coloring instructions for lines based on whether gene is in significant list
         carotenoids$sig <- NA
         for(i in 1:nrow(carotenoids)) {
           if (carotenoids[i,1] %in% combined_sigs$Gene) {
             carotenoids[i,8] <- "sig"
           } else {
             carotenoids[i,8] <- "nonsig"
           }
         }
         
        # Make line plot   
         carotenoids.plot <- merge(carotenoids, zscores.averaged,  by = "Gene", all = FALSE) |>
           ggplot(aes(x=time, y=mean_zscore, group = Gene)) +
           facet_wrap(~Process,  nrow = 2, strip.position = "bottom")+
           geom_line(aes(colour = sig),linewidth=1)+
           geom_point_interactive(aes(data_id = Gene, tooltip = tip), color = "grey20", size = 3)+
           scale_y_continuous(name="Z-scores", breaks = seq(-2,2, by =0.4))+
           scale_x_discrete(name="Time Points")+
           theme(axis.text.x = element_text(size = 9),
                 axis.title = element_text(size = 11),
                 axis.text.y = element_text(size=9),
                 legend.position = "bottom", legend.direction = "horizontal",
                 title = element_text(size= 10))+
           scale_color_manual(name = "Significance:", values = c("grey40", "dodgerblue"))+
           ggtitle("Carotenoid genes")
         carotenoids.plot
         
        # Export as jpeg
         jpeg(filename = "carotenoids genes time series line plots.jpg", width = 12, height = 10, units = "in", res = 300)
         carotenoids.plot
         dev.off()
         
        # Generate interactive figure   
         girafe(
           ggobj = carotenoids.plot,
           options = list(opts_hover_inv(css = "opacity:0.1;stroke:yellow"), 
                          opts_tooltip(offx = 0, offy = -50, css = 'font-size:larger;background-color:#000000;color:white')
           ),
           height_svg = 5,
           width_svg = 10
         )
         
         
     # Mito targeted genes
         
        # Read in data and reformat gene names to remove PA notation if present
         mitos <- read.csv(file = "Mito-target_proteins_SDv2.2_noPA.csv", header = TRUE)
         names(mitos)[1] <- "Gene"
         mitos$tip <- paste0(mitos$Gene,": ", mitos$Description)
         mitos$Gene <-str_split_i(mitos$Gene, "-", 1)
         
        # Add coloring instructions for lines based on whether gene is in significant list
         mitos$sig <- NA
         for(i in 1:nrow(mitos)) {
           if (mitos[i,1] %in% combined_sigs$Gene) {
             mitos[i,5] <- "sig"
           } else {
             mitos[i,5] <- "nonsig"
           }
         }
         
        # Make line plot   
         mitos.plot <- merge(mitos, zscores.averaged,  by = "Gene", all = FALSE) |>
           ggplot(aes(x=time, y=mean_zscore, group = Gene)) +
           facet_wrap(~Process,  nrow = 2, strip.position = "bottom")+
           geom_line(aes(colour = sig),linewidth=1)+
           geom_point_interactive(aes(data_id = Gene, tooltip = tip), color = "grey20", size = 2)+
           scale_y_continuous(name="Z-scores", breaks = seq(-2,2, by =0.4))+
           scale_x_discrete(name="Time Points")+
           theme(axis.text.x = element_text(size = 9),
                 axis.title = element_text(size = 11),
                 axis.text.y = element_text(size=9),
                 legend.position = "bottom", legend.direction = "horizontal",
                 title = element_text(size= 10))+
           scale_color_manual(name = "Significance:", values = c("grey40", "dodgerblue"))+
           ggtitle("Mitochondrial genes")
         mitos.plot
         
        # Export as jpeg
         jpeg(filename = "Mitochondrial genes time series line plots.jpg", width = 12, height = 10, units = "in", res = 300)
         mitos.plot
         dev.off()
         
        # Generate interactive figure   
         girafe(
           ggobj = mitos.plot,
           options = list(opts_hover_inv(css = "opacity:0.1;stroke:yellow"), 
                          opts_tooltip(offx = 0, offy = -50, css = 'font-size:larger;background-color:#000000;color:white')
           ),
           height_svg = 5,
           width_svg = 10
         )
         
         
   #### Exploring functional groups ####
      
    ## Response to oxygen GO term 
         
        # Pull out all gene IDs associated with GO term "Response to Oxygen containing compound" GO:1901700
         genes2go %>%
           filter(str_detect(V2, "GO:1901700")) -> all.rto.genes
         
        # How many are in just the masigpro list and which clusters?
         sum(sigs.masigpro$Gene %in% all.rto.genes$V1)  
         sigs.masigpro[sigs.masigpro$Gene %in% all.rto.genes$V1, c(32,31)]
         
          # Create small table of significant gene totals for each comparison to the control
           response.to.oxgyen <- data.frame(times = c("Hypoxia vs control","Pcrit vs control","Anoxia vs control",
                                                      "Recovery vs control"),
                                            sigs = c(18,42,72,79),
                                            total = rep(1022,4),
                                            all_annotated_genes = rep(11461, 4), #all annotated genes WITH a GO term
                                            all_sigs = c(203,489,865,787)
           )
           response.to.oxgyen$all_other_sigs = response.to.oxgyen$all_sigs - response.to.oxgyen$sigs
           response.to.oxgyen$ns = response.to.oxgyen$total - response.to.oxgyen$sigs
           response.to.oxgyen$all_other_ns = response.to.oxgyen$all_annotated_genes - (response.to.oxgyen$all_other_sigs + response.to.oxgyen$ns + response.to.oxgyen$sigs)
          
          # Reorder factor to plot how we want
           response.to.oxgyen$times <- factor(response.to.oxgyen$times, 
                                              levels = c("Recovery vs control","Anoxia vs control","Pcrit vs control",
                                                         "Hypoxia vs control"))
          
        # Calculate odds ratio with confidence intervals using epitools instead of metafor
          # Hypoxia vs control
           # Create contingency table matrix
           tapw <- c("Oxygen", "All_others")
           outc <- c("Significant", "Nonsignificant")	
           dat <- matrix(c(response.to.oxgyen[1,2], response.to.oxgyen[1,7], response.to.oxgyen[1,6], response.to.oxgyen[1,8])
                         ,2,2,byrow=TRUE)
           dimnames(dat) <- list("Groups" = tapw, "Comparisons" = outc)
           oddsratio.fisher(dat)
           
        # Make simple bar plot of the number of genes in response to oxygen compound GO term
          # Create and export figure
           jpeg(filename = "Response to oxygen genes DESeq2 totals.jpg", width = 6.5, height = 4, units = "in", res = 300)
           response.to.oxgyen |>
             ggplot() +
             geom_segment(aes(x = 0, xend = sigs, y = times, yend = times), linewidth = 2)+
             geom_point(aes(x = sigs, y = times), size = 9, col = "darkblue")+
             geom_text(aes(x = sigs, y = times, label = sigs), color = "white") +
             scale_x_continuous(expand = expansion(mult = c(0, 0.09)), breaks = seq(0, 90, by = 10))+
             xlab("Number of genes") + ylab("") +
             theme(axis.text.x = element_text(size = 12),
                   axis.title.y = element_blank(),
                   title = element_text(size = 9),
                   panel.grid.major.x = element_line(color = "lightgrey"))+
             ggtitle("Significant DEG's out of 1022 in 'Response to Oxygen' GO TERM")
           dev.off()
         
        # How many significant genes are shared between the each time point's RTO GO term?
           # Create venn diagram first
           # Make list of genes in each group for venn diagram of time course programs
           # Deseq comparison
           RTO.venn.list <- list(
             hypoxia_vs_control = allRes.deseq.hypox.cont$genes$`GO:1901700`,
             pcrit_vs_control = allRes.deseq.pcrit.cont$genes$`GO:1901700`,
             anoxia_vs_control = allRes.deseq.anox.cont$genes$`GO:1901700`,
             recovery_vs_control = allRes.deseq.recov.cont$genes$`GO:1901700`,
             anoxia_vs_recovery = allRes.deseq.anox.recov$genes$`GO:1901700`,
             hypoxia_vs_anoxia = allRes.deseq.hypox.anox$genes$`GO:1901700`
           )
           
           # Convert list to data frame to pull out shared and unique genes
           RTO.df <- data.frame(unlist(RTO.venn.list))
           RTO.df$groups <- rownames(RTO.df)
           rownames(RTO.df) <- NULL
           names(RTO.df)[1] <- "Gene"
           RTO.df$groups <- gsub("\\d+$", "", RTO.df$groups) #Remove numbers at end of group names
           
        # Draw venn diagram of all five groups using venn package  
           jpeg(filename = "Venn Diagram Response to Oxygen GO Term Genes.jpg", width = 8.5, height = 5.5, units = "in", res = 300)
           ggVennDiagram(RTO.venn.list[c(1,2,3,4)], label_percent_digit = 1, label_size = 4) +
             scale_x_continuous(expand = expansion(mult = .2)) +
             scale_fill_distiller(name = "Counts:", palette = "GnBu", direction = 1) + 
             theme(legend.position = "bottom", legend.direction = "horizontal", 
                   plot.margin = unit(c(0,0,0,0), "cm"))
           dev.off()
           
          
       # Pull out genes unique to each group
           # NOTE: %notin% function is a custom function assigned at top of script
        # Genes unique to hypoxia vs control
         RTO.venn.list$hypoxia_vs_control[RTO.venn.list$hypoxia_vs_control %notin% 
                                                                     c(RTO.venn.list$pcrit_vs_control,
                                                                      RTO.venn.list$anoxia_vs_control,
                                                                      RTO.venn.list$recovery_vs_control,
                                                                      RTO.venn.list$anoxia_vs_recovery,
                                                                      RTO.venn.list$hypoxia_vs_anoxia)
           ]
         
        # Genes unique to pcrit vs control
         RTO.venn.list$pcrit_vs_control[RTO.venn.list$pcrit_vs_control %notin% c(RTO.venn.list$hypoxia_vs_control,
                                                                                        RTO.venn.list$anoxia_vs_control,
                                                                                        RTO.venn.list$recovery_vs_control,
                                                                                 RTO.venn.list$anoxia_vs_recovery,
                                                                                 RTO.venn.list$hypoxia_vs_anoxia)
         ]
         
        # Genes unique to anoxia vs control
         RTO.venn.list$anoxia_vs_control[RTO.venn.list$anoxia_vs_control %notin% c(RTO.venn.list$pcrit_vs_control,
                                                                                        RTO.venn.list$hypoxia_vs_control,
                                                                                        RTO.venn.list$recovery_vs_control,
                                                                                   RTO.venn.list$anoxia_vs_recovery,
                                                                                   RTO.venn.list$hypoxia_vs_anoxia)
         ]
         
        # Genes unique to recovery vs control
         RTO.venn.list$recovery_vs_control[RTO.venn.list$recovery_vs_control %notin% c(RTO.venn.list$pcrit_vs_control,
                                                                                          RTO.venn.list$anoxia_vs_control,
                                                                                          RTO.venn.list$hypoxia_vs_control,
                                                                                       RTO.venn.list$anoxia_vs_recovery,
                                                                                       RTO.venn.list$hypoxia_vs_anoxia)
         ]

         
        # Pull out genes shared among all groups 
         RTO.venn.list$hypoxia_vs_control[RTO.venn.list$hypoxia_vs_control %in% RTO.venn.list$pcrit_vs_control &
                                                RTO.venn.list$hypoxia_vs_control %in% RTO.venn.list$anoxia_vs_control &
                                              RTO.venn.list$hypoxia_vs_control %in% RTO.venn.list$recovery_vs_control &
                                            RTO.venn.list$hypoxia_vs_control %in% RTO.venn.list$anoxia_vs_recovery &
                                            RTO.venn.list$hypoxia_vs_control %in% RTO.venn.list$hypoxia_vs_anoxia
                                              ]
         
        # Pull out all genes that are in hypoxia, pcrit, and anoxia but NOT recovery
         unique(c(RTO.venn.list$hypoxia_vs_control, 
           RTO.venn.list$pcrit_vs_control, 
           RTO.venn.list$anoxia_vs_control)[c(RTO.venn.list$hypoxia_vs_control, 
           RTO.venn.list$pcrit_vs_control, 
           RTO.venn.list$anoxia_vs_control) %notin% RTO.venn.list$recovery_vs_control
           ])
         
       # Query masigpro clusters to see if any of the RTO genes from the DESeq2 contrasts were found significant and clustered by masigpro 
        deseq.RTO.in.masigpro <- data.frame(sigs.masigpro[sigs.masigpro$Gene %in% 
                                                            unique(c(RTO.venn.list$hypoxia_vs_control, 
                                                                     RTO.venn.list$pcrit_vs_control, 
                                                                      RTO.venn.list$anoxia_vs_control, 
                                                                      RTO.venn.list$recovery_vs_control)), c(32,31)
                                                          ])
        
      # Grab gene annotations from combined_sigs frame
        deseq.RTO.in.masigpro <- merge(deseq.RTO.in.masigpro, combined_sigs[,1:2], by = "Gene", all = FALSE)
        
        table(deseq.RTO.in.masigpro$clusters) #Tally genes in each cluster

      # How many different genes grouped under response to oxygen were detected by maSigPro
        # Combine list of genes from DESeq comparisons found in significant genes by maSigPro with those just found by maSigPro
          length(unique(c(deseq.RTO.in.masigpro$Gene, sigs.masigpro[sigs.masigpro$Gene %in% all.rto.genes$V1, c(32)])))
          
        # Pull out the functions and clusters from combined sigs data frame of the above list of genes
          combined_sigs[combined_sigs$Gene %in% c(deseq.RTO.in.masigpro$Gene, sigs.masigpro[sigs.masigpro$Gene %in% all.rto.genes$V1, c(32)]), c(1,4)]
          
        
    ## Chitin and cuticle genes
        
      # Sum the number of chitin/cuticle genes in our significant gene lists
        sum(combined_sigs$Gene %in% cuticles$Gene)
        combined_sigs[combined_sigs$Gene %in% cuticles$Gene, c(1,2)] #Which genes and description
        
        # How many are in the masigpro list and which clusters?
          sum(sigs.masigpro$Gene %in% cuticles$Gene)  
          sigs.masigpro[sigs.masigpro$Gene %in% cuticles$Gene, c(32,31)]
        
        # How many are in the hypoxia vs control DEG list?
          sum(Sigs.hypox.vs.cont$Gene %in% cuticles$Gene)  
          Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% cuticles$Gene, c(1,3,7)]
        
        # How many are in the pcrit vs control DEG list?
          sum(Sigs.pcrit.vs.cont$Gene %in% cuticles$Gene)  
          Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% cuticles$Gene, c(1,3,7)]
        
        # How many are in the anoxia vs control DEG list?
          sum(Sigs.anox.vs.control$Gene %in% cuticles$Gene)  
          Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% cuticles$Gene, c(1,3,7)]
        
        # How many are in the recovery vs control DEG list?
          sum(Sigs.recov.vs.cont$Gene %in% cuticles$Gene)  
          Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% cuticles$Gene, c(1,3,7)]
            
          
        # Make simple bar plot of the number of cuticle genes
          #Make data frame for plotting
          cuticles.totals <- data.frame(times = c("Hypoxia vs control","Pcrit vs control","Anoxia vs control",
                                                     "Recovery vs control"),
                                           sigs = c(9,21,21,29),
                                           total = rep(224,4)
          )
          # Reorder factor to plot how we want
          cuticles.totals$times <- factor(cuticles.totals$times, 
                                             levels = c("Recovery vs control","Anoxia vs control","Pcrit vs control",
                                                        "Hypoxia vs control"))
          # Create and export figure
          jpeg(filename = "Chitin and cuticle genes DESeq2 totals.jpg", width = 6.5, height = 4, units = "in", res = 300)
          cuticles.totals |>
            ggplot() +
            geom_segment(aes(x = 0, xend = sigs, y = times, yend = times), linewidth = 2)+
            geom_point(aes(x = sigs, y = times), size = 9, col = "darkblue")+
            geom_text(aes(x = sigs, y = times, label = sigs), color = "white") +
            scale_x_continuous(expand = expansion(mult = c(0, 0.09)), breaks = seq(0, 30, by = 5))+
            xlab("Number of genes") + ylab("") +
            theme(axis.text.x = element_text(size = 12),
                  axis.title.y = element_blank(),
                  title = element_text(size = 9),
                  panel.grid.major.x = element_line(color = "lightgrey"))+
            ggtitle("Significant DEG's out of 223 chitin and cuticle genes")
          dev.off()
          
        # How many significant genes are shared between the each time point's RTO GO term?
        # Create venn diagram first
        # Make list of genes in each group for venn diagram of time course programs
        # Deseq comparison
          cuticles.venn.list <- list(
            hypoxia_vs_control = Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% cuticles$Gene, 1],
            pcrit_vs_control = Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% cuticles$Gene, 1],
            anoxia_vs_control = Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% cuticles$Gene, 1],
            recovery_vs_control =  Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% cuticles$Gene, 1],
            Anoxia_vs_recovery = Sigs.anox.vs.recov[Sigs.anox.vs.recov$Gene %in% cuticles$Gene, 1],
            Hypoxia_vs_anoxia = Sigs.hypox.vs.anox[Sigs.hypox.vs.anox$Gene %in% cuticles$Gene, 1]
          )
        
          
        # Draw venn diagram of first four groups using venn package  
          jpeg(filename = "Venn Diagram Chitin and Cuticle Genes.jpg", width = 8.5, height = 5.5, units = "in", res = 300)
          ggVennDiagram(cuticles.venn.list[c(1,2,3,4)], label_percent_digit = 1, label_size = 4) +
            scale_x_continuous(expand = expansion(mult = .2)) +
            scale_fill_distiller(name = "Counts:", palette = "GnBu", direction = 1) + 
            theme(legend.position = "bottom", legend.direction = "horizontal", 
                  plot.margin = unit(c(0,0,0,0), "cm"))
          dev.off()
          
  # Pull out genes unique to each group
    # NOTE: %notin% function is a custom function assigned at top of script
    # Genes unique to hypoxia vs control
      cuticles.venn.list$hypoxia_vs_control[cuticles.venn.list$hypoxia_vs_control %notin% 
                                                                  c(cuticles.venn.list$pcrit_vs_control,
                                                                    cuticles.venn.list$anoxia_vs_control,
                                                                    cuticles.venn.list$recovery_vs_control,
                                                                    cuticles.venn.list$Anoxia_vs_recovery,
                                                                    cuticles.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to pcrit vs control
      cuticles.venn.list$pcrit_vs_control[cuticles.venn.list$pcrit_vs_control %notin% c(cuticles.venn.list$hypoxia_vs_control,
                                                                                                cuticles.venn.list$anoxia_vs_control,
                                                                                                cuticles.venn.list$recovery_vs_control,
                                                                                        cuticles.venn.list$Anoxia_vs_recovery,
                                                                                        cuticles.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to anoxia vs control
      cuticles.venn.list$anoxia_vs_control[cuticles.venn.list$anoxia_vs_control %notin% c(cuticles.venn.list$pcrit_vs_control,
                                                                                                  cuticles.venn.list$hypoxia_vs_control,
                                                                                                  cuticles.venn.list$recovery_vs_control,
                                                                                          cuticles.venn.list$Anoxia_vs_recovery,
                                                                                          cuticles.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to recovery vs control
      cuticles.venn.list$recovery_vs_control[cuticles.venn.list$recovery_vs_control %notin% c(cuticles.venn.list$pcrit_vs_control,
                                                                                                      cuticles.venn.list$anoxia_vs_control,
                                                                                                      cuticles.venn.list$hypoxia_vs_control,
                                                                                              cuticles.venn.list$Anoxia_vs_recovery,
                                                                                              cuticles.venn.list$Hypoxia_vs_anoxia)
      ]
      
      
    # Pull out genes shared among all groups 
      cuticles.venn.list$hypoxia_vs_control[cuticles.venn.list$hypoxia_vs_control %in% cuticles.venn.list$pcrit_vs_control &
                                                  cuticles.venn.list$hypoxia_vs_control %in% cuticles.venn.list$anoxia_vs_control &
                                                  cuticles.venn.list$hypoxia_vs_control %in% cuticles.venn.list$recovery_vs_control &
                                              cuticles.venn.list$hypoxia_vs_control %in% cuticles.venn.list$Anoxia_vs_recovery &
                                            cuticles.venn.list$hypoxia_vs_control %in% cuticles.venn.list$Hypoxia_vs_anoxia 
      ]
      
    # Pull out all genes that are in hypoxia, pcrit, and anoxia but NOT recovery
      unique(c(cuticles.venn.list$hypoxia_vs_control, 
               cuticles.venn.list$pcrit_vs_control, 
               cuticles.venn.list$anoxia_vs_control)[c(cuticles.venn.list$hypoxia_vs_control, 
                                                           cuticles.venn.list$pcrit_vs_control, 
                                                           cuticles.venn.list$anoxia_vs_control) %notin% cuticles.venn.list$recovery_vs_control
               ])
      
    # Query masigpro clusters to see if any of the genes from the DESeq2 contrasts were found significant and clustered by masigpro 
     cuticles.in.masigpro <- data.frame(sigs.masigpro[sigs.masigpro$Gene %in% 
                                                          unique(c(cuticles.venn.list$hypoxia_vs_control, 
                                                                   cuticles.venn.list$pcrit_vs_control, 
                                                                   cuticles.venn.list$anoxia_vs_control, 
                                                                   cuticles.venn.list$recovery_vs_control,
                                                                   cuticles.venn.list$Anoxia_vs_recovery,
                                                                   cuticles.venn.list$Hypoxia_vs_anoxia)), c(32,31)
      ])
    # Grab gene annotations from combined_sigs frame
     cuticles.in.masigpro <- merge(cuticles.in.masigpro, combined_sigs[,1:2], by = "Gene", all = FALSE)
     
      table(cuticles.in.masigpro$clusters) #Tally genes in each cluster
      
    #Generate pathway
      cuticles.pathway <- highlight_entities("tcf00520", c("131877547", "131886326", "131885557", "131884437", "131879841",
                                                           "131887274", "131880429", "131890320"), 
                                             fill_color = "paleturquoise", legend_name="Significant")
      cuticles.pathway
      
    #Export as png
      ggkeggsave(filename="cuticles pathway.png", cuticles.pathway, dpi=300)
      
    ### Generate heatmap of just the sig genes from the expression averaged wide dataframe
      cuticles.masigpro<- as.data.frame(expression.averaged.wide[expression.averaged.wide$Gene %in% combined_sigs[combined_sigs$Gene %in% cuticles$Gene, 1],])
      row.names(cuticles.masigpro) <- cuticles.masigpro$Gene
      cuticles.masigpro <- cuticles.masigpro[,-c(1)]
      
      
      jpeg(file = "cuticles Sigs Heatmap.jpg", width = 5, height = 25, units= "in", res = 300)
      pheatmap(cuticles.masigpro, cluster_rows = F, cluster_cols = F,
               show_rownames=T, border_color=NA,
               scale = "row", 
               labels_col = c("Normoxia", "Mild Hypox", "Pcrit", "Anoxia", "Recovery"),
               fontsize = 9,fontsize_row = 10, height=20, main = "cuticles")
      dev.off()
      
      
  ## Carotenoid genes
      
    # Sum the number of carotenoid genes in our significant gene lists
      sum(combined_sigs$Gene %in% carotenoids$Gene)
      combined_sigs[combined_sigs$Gene %in% carotenoids$Gene, c(1,2)] #Which genes and description
      
    # How many are in the masigpro list and which clusters?
      sum(sigs.masigpro$Gene %in% carotenoids$Gene)  
      sigs.masigpro[sigs.masigpro$Gene %in% carotenoids$Gene, c(32,31)]
      
     # How many are in the hypoxia vs control DEG list?
      sum(Sigs.hypox.vs.cont$Gene %in% carotenoids$Gene)  
      Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% carotenoids$Gene, c(1,3,7)]
      
    # How many are in the pcrit vs control DEG list?
      sum(Sigs.pcrit.vs.cont$Gene %in% carotenoids$Gene)  
      Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% carotenoids$Gene, c(1,3,7)]
      
    # How many are in the anoxia vs control DEG list?
      sum(Sigs.anox.vs.control$Gene %in% carotenoids$Gene)  
      Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% carotenoids$Gene, c(1,3,7)]
      
    # How many are in the recovery vs control DEG list?
      sum(Sigs.recov.vs.cont$Gene %in% carotenoids$Gene)  
      Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% carotenoids$Gene, c(1,3,7)]
      
      
    # Make simple bar plot of the number of genes in response to oxygen compound GO term
    #Make data frame for plotting
      carotenoids.totals <- data.frame(times = c("Hypoxia vs control","Pcrit vs control","Anoxia vs control",
                                              "Recovery vs control"),
                                    sigs = c(1,4,6,7),
                                    total = rep(34,4)
      )
    # Reorder factor to plot how we want
      carotenoids.totals$times <- factor(carotenoids.totals$times, 
                                      levels = c("Recovery vs control","Anoxia vs control","Pcrit vs control",
                                                 "Hypoxia vs control"))
    # Create and export figure
      jpeg(filename = "Carotenoids genes DESeq2 totals.jpg", width = 6.5, height = 4, units = "in", res = 300)
      carotenoids.totals |>
        ggplot() +
        geom_segment(aes(x = 0, xend = sigs, y = times, yend = times), linewidth = 2)+
        geom_point(aes(x = sigs, y = times), size = 9, col = "darkblue")+
        geom_text(aes(x = sigs, y = times, label = sigs), color = "white") +
        scale_x_continuous(expand = expansion(mult = c(0, 0.09)), breaks = seq(0, 10, by = 2))+
        xlab("Number of genes") + ylab("") +
        theme(axis.text.x = element_text(size = 12),
              axis.title.y = element_blank(),
              title = element_text(size = 9),
              panel.grid.major.x = element_line(color = "lightgrey"))+
        ggtitle("Significant DEG's out of 22 carotenoid genes")
      dev.off()
      
    # How many significant genes are shared between the each time point's RTO GO term?
    # Create venn diagram first
    # Make list of genes in each group for venn diagram of time course programs
    # Deseq comparison
      carotenoids.venn.list <- list(
        hypoxia_vs_control = Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% carotenoids$Gene, 1],
        pcrit_vs_control = Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% carotenoids$Gene, 1],
        anoxia_vs_control = Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% carotenoids$Gene, 1],
        recovery_vs_control =  Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% carotenoids$Gene, 1],
        Anoxia_vs_recovery = Sigs.anox.vs.recov[Sigs.anox.vs.recov$Gene %in% carotenoids$Gene, 1],
        Hypoxia_vs_anoxia = Sigs.hypox.vs.anox[Sigs.hypox.vs.anox$Gene %in% carotenoids$Gene, 1]
      )
      
      
    # Draw venn diagram of the four groups using venn package  
      jpeg(filename = "Venn Diagram Carotenoid Genes.jpg", width = 8.5, height = 5.5, units = "in", res = 300)
      ggVennDiagram(carotenoids.venn.list[c(1,2,3,4)], label_percent_digit = 1, label_size = 4) +
        scale_x_continuous(expand = expansion(mult = .2)) +
        scale_fill_distiller(name = "Counts:", palette = "GnBu", direction = 1) + 
        theme(legend.position = "bottom", legend.direction = "horizontal", 
              plot.margin = unit(c(0,0,0,0), "cm"))
      dev.off()
      
    # Pull out genes unique to each group
    # NOTE: %notin% function is a custom function assigned at top of script
    # Genes unique to hypoxia vs control
      carotenoids.venn.list$hypoxia_vs_control[carotenoids.venn.list$hypoxia_vs_control %notin% 
                                              c(carotenoids.venn.list$pcrit_vs_control,
                                                carotenoids.venn.list$anoxia_vs_control,
                                                carotenoids.venn.list$recovery_vs_control,
                                                carotenoids.venn.list$Anoxia_vs_recovery,
                                                carotenoids.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to pcrit vs control
      carotenoids.venn.list$pcrit_vs_control[carotenoids.venn.list$pcrit_vs_control %notin% c(carotenoids.venn.list$hypoxia_vs_control,
                                                                                        carotenoids.venn.list$anoxia_vs_control,
                                                                                        carotenoids.venn.list$recovery_vs_control,
                                                                                        carotenoids.venn.list$Anoxia_vs_recovery,
                                                                                        carotenoids.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to anoxia vs control
      carotenoids.venn.list$anoxia_vs_control[carotenoids.venn.list$anoxia_vs_control %notin% c(carotenoids.venn.list$pcrit_vs_control,
                                                                                          carotenoids.venn.list$hypoxia_vs_control,
                                                                                          carotenoids.venn.list$recovery_vs_control,
                                                                                          carotenoids.venn.list$Anoxia_vs_recovery,
                                                                                          carotenoids.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to recovery vs control
      carotenoids.venn.list$recovery_vs_control[carotenoids.venn.list$recovery_vs_control %notin% c(carotenoids.venn.list$pcrit_vs_control,
                                                                                              carotenoids.venn.list$anoxia_vs_control,
                                                                                              carotenoids.venn.list$hypoxia_vs_control,
                                                                                              carotenoids.venn.list$Anoxia_vs_recovery,
                                                                                              carotenoids.venn.list$Hypoxia_vs_anoxia)
      ]
      
      
    # Pull out genes shared among all groups 
      carotenoids.venn.list$hypoxia_vs_control[carotenoids.venn.list$hypoxia_vs_control %in% carotenoids.venn.list$pcrit_vs_control &
                                              carotenoids.venn.list$hypoxia_vs_control %in% carotenoids.venn.list$anoxia_vs_control &
                                              carotenoids.venn.list$hypoxia_vs_control %in% carotenoids.venn.list$recovery_vs_control &
                                              carotenoids.venn.list$hypoxia_vs_control %in% carotenoids.venn.list$Anoxia_vs_recovery &
                                              carotenoids.venn.list$hypoxia_vs_control %in% carotenoids.venn.list$Hypoxia_vs_anoxia 
      ]
      
    # Pull out all genes that are in hypoxia, pcrit, and anoxia but NOT recovery
      unique(c(carotenoids.venn.list$hypoxia_vs_control, 
               carotenoids.venn.list$pcrit_vs_control, 
               carotenoids.venn.list$anoxia_vs_control)[c(carotenoids.venn.list$hypoxia_vs_control, 
                                                       carotenoids.venn.list$pcrit_vs_control, 
                                                       carotenoids.venn.list$anoxia_vs_control) %notin% carotenoids.venn.list$recovery_vs_control
               ])
      
    # Query masigpro clusters to see if any of the genes from the DESeq2 contrasts were found significant and clustered by masigpro 
      carotenoids.in.masigpro <- data.frame(sigs.masigpro[sigs.masigpro$Gene %in% 
                                                         unique(c(carotenoids.venn.list$hypoxia_vs_control, 
                                                                  carotenoids.venn.list$pcrit_vs_control, 
                                                                  carotenoids.venn.list$anoxia_vs_control, 
                                                                  carotenoids.venn.list$recovery_vs_control,
                                                                  carotenoids.venn.list$Anoxia_vs_recovery,
                                                                  carotenoids.venn.list$Hypoxia_vs_anoxia)), c(32,31)
      ])
    # Grab gene annotations from combined_sigs frame
      carotenoids.in.masigpro <- merge(carotenoids.in.masigpro, combined_sigs[,1:2], by = "Gene", all = FALSE)
      
      table(carotenoids.in.masigpro$clusters) #Tally genes in each cluster
      
    ### Generate heatmap of just the sig genes from the expression averaged wide dataframe
      carotenoids.masigpro<- as.data.frame(expression.averaged.wide[expression.averaged.wide$Gene %in% combined_sigs[combined_sigs$Gene %in% carotenoids$Gene, 1],])
      row.names(carotenoids.masigpro) <- carotenoids.masigpro$Gene
      carotenoids.masigpro <- carotenoids.masigpro[,-c(1)]
      
      
      jpeg(file = "carotenoids Sigs Heatmap.jpg", width = 5, height = 7, units= "in", res = 300)
      pheatmap(carotenoids.masigpro, cluster_rows = F, cluster_cols = F,
               show_rownames=T, border_color=NA,
               scale = "row", 
               labels_col = c("Normoxia", "Mild Hypox", "Pcrit", "Anoxia", "Recovery"),
               fontsize = 9,fontsize_row = 10, height=20, main = "carotenoids")
      dev.off()
      
      
  ## Antioxidant genes
      
    # Sum the number of antioxidant genes in our significant gene lists
      sum(combined_sigs$Gene %in% antioxidants$Gene)
      combined_sigs[combined_sigs$Gene %in% antioxidants$Gene, c(1,2)] #Which genes and description
      
    # How many are in the masigpro list and which clusters?
      sum(sigs.masigpro$Gene %in% antioxidants$Gene)  
      sigs.masigpro[sigs.masigpro$Gene %in% antioxidants$Gene, c(32,31)]
      
    # How many are in the hypoxia vs control DEG list?
      sum(Sigs.hypox.vs.cont$Gene %in% antioxidants$Gene)  
      Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% antioxidants$Gene, c(1,3,7)]
      
    # How many are in the pcrit vs control DEG list?
      sum(Sigs.pcrit.vs.cont$Gene %in% antioxidants$Gene)  
      Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% antioxidants$Gene, c(1,3,7)]
      
    # How many are in the anoxia vs control DEG list?
      sum(Sigs.anox.vs.control$Gene %in% antioxidants$Gene)  
      Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% antioxidants$Gene, c(1,3,7)]
      
    # How many are in the recovery vs control DEG list?
      sum(Sigs.recov.vs.cont$Gene %in% antioxidants$Gene)  
      Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% antioxidants$Gene, c(1,3,7)]
      
      
    # Make simple bar plot of the number of genes in response to oxygen compound GO term
    #Make data frame for plotting
      antioxidants.totals <- data.frame(times = c("Hypoxia vs control","Pcrit vs control","Anoxia vs control",
                                                 "Recovery vs control"),
                                       sigs = c(3,6,8,8),
                                       total = rep(34,4)
      )
    # Reorder factor to plot how we want
      antioxidants.totals$times <- factor(antioxidants.totals$times, 
                                         levels = c("Recovery vs control","Anoxia vs control","Pcrit vs control",
                                                    "Hypoxia vs control"))
    # Create and export figure
      jpeg(filename = "Antioxidant genes DESeq2 totals.jpg", width = 6.5, height = 4, units = "in", res = 300)
      antioxidants.totals |>
        ggplot() +
        geom_segment(aes(x = 0, xend = sigs, y = times, yend = times), linewidth = 2)+
        geom_point(aes(x = sigs, y = times), size = 9, col = "darkblue")+
        geom_text(aes(x = sigs, y = times, label = sigs), color = "white") +
        scale_x_continuous(expand = expansion(mult = c(0, 0.09)), breaks = seq(0, 10, by = 2))+
        xlab("Number of genes") + ylab("") +
        theme(axis.text.x = element_text(size = 12),
              axis.title.y = element_blank(),
              title = element_text(size = 9),
              panel.grid.major.x = element_line(color = "lightgrey"))+
        ggtitle("Significant DEG's out of 48 antioxidant genes")
      dev.off()
      
    # How many significant genes are shared between the each time point's RTO GO term?
    # Create venn diagram first
    # Make list of genes in each group for venn diagram of time course programs
    # Deseq comparison
      antioxidants.venn.list <- list(
        hypoxia_vs_control = Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% antioxidants$Gene, 1],
        pcrit_vs_control = Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% antioxidants$Gene, 1],
        anoxia_vs_control = Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% antioxidants$Gene, 1],
        recovery_vs_control =  Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% antioxidants$Gene, 1],
        Anoxia_vs_recovery = Sigs.anox.vs.recov[Sigs.anox.vs.recov$Gene %in% antioxidants$Gene, 1],
        Hypoxia_vs_anoxia = Sigs.hypox.vs.anox[Sigs.hypox.vs.anox$Gene %in% antioxidants$Gene, 1]
      )
      
      
    # Draw venn diagram of all five groups using venn package  
      jpeg(filename = "Venn Diagram Antioxidant Genes.jpg", width = 8.5, height = 5.5, units = "in", res = 300)
      ggVennDiagram(antioxidants.venn.list[c(1,2,3,4)], label_percent_digit = 1, label_size = 4) +
        scale_x_continuous(expand = expansion(mult = .2)) +
        scale_fill_distiller(name = "Counts:", palette = "GnBu", direction = 1) + 
        theme(legend.position = "bottom", legend.direction = "horizontal", 
              plot.margin = unit(c(0,0,0,0), "cm"))
      dev.off()
      
    # Pull out genes unique to each group
    # NOTE: %notin% function is a custom function assigned at top of script
    # Genes unique to hypoxia vs control
      antioxidants.venn.list$hypoxia_vs_control[antioxidants.venn.list$hypoxia_vs_control %notin% 
                                                 c(antioxidants.venn.list$pcrit_vs_control,
                                                   antioxidants.venn.list$anoxia_vs_control,
                                                   antioxidants.venn.list$recovery_vs_control,
                                                   antioxidants.venn.list$Anoxia_vs_recovery,
                                                   antioxidants.venn.list$Hypoxia_vs_anoxia)
      ]
      
     # Genes unique to pcrit vs control
      antioxidants.venn.list$pcrit_vs_control[antioxidants.venn.list$pcrit_vs_control %notin% c(antioxidants.venn.list$hypoxia_vs_control,
                                                                                              antioxidants.venn.list$anoxia_vs_control,
                                                                                              antioxidants.venn.list$recovery_vs_control,
                                                                                              antioxidants.venn.list$Anoxia_vs_recovery,
                                                                                              antioxidants.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to anoxia vs control
      antioxidants.venn.list$anoxia_vs_control[antioxidants.venn.list$anoxia_vs_control %notin% c(antioxidants.venn.list$pcrit_vs_control,
                                                                                                antioxidants.venn.list$hypoxia_vs_control,
                                                                                                antioxidants.venn.list$recovery_vs_control,
                                                                                                antioxidants.venn.list$Anoxia_vs_recovery,
                                                                                                antioxidants.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to recovery vs control
      antioxidants.venn.list$recovery_vs_control[antioxidants.venn.list$recovery_vs_control %notin% c(antioxidants.venn.list$pcrit_vs_control,
                                                                                                    antioxidants.venn.list$anoxia_vs_control,
                                                                                                    antioxidants.venn.list$hypoxia_vs_control,
                                                                                                    antioxidants.venn.list$Anoxia_vs_recovery,
                                                                                                    antioxidants.venn.list$Hypoxia_vs_anoxia)
      ]
      
      
    # Pull out genes shared among all groups 
      antioxidants.venn.list$hypoxia_vs_control[antioxidants.venn.list$hypoxia_vs_control %in% antioxidants.venn.list$pcrit_vs_control &
                                                 antioxidants.venn.list$hypoxia_vs_control %in% antioxidants.venn.list$anoxia_vs_control &
                                                 antioxidants.venn.list$hypoxia_vs_control %in% antioxidants.venn.list$recovery_vs_control &
                                                 antioxidants.venn.list$hypoxia_vs_control %in% antioxidants.venn.list$Anoxia_vs_recovery &
                                                 antioxidants.venn.list$hypoxia_vs_control %in% antioxidants.venn.list$Hypoxia_vs_anoxia 
      ]
      
    # Pull out all genes that are in hypoxia, pcrit, and anoxia but NOT recovery
      unique(c(antioxidants.venn.list$hypoxia_vs_control, 
               antioxidants.venn.list$pcrit_vs_control, 
               antioxidants.venn.list$anoxia_vs_control)[c(antioxidants.venn.list$hypoxia_vs_control, 
                                                          antioxidants.venn.list$pcrit_vs_control, 
                                                          antioxidants.venn.list$anoxia_vs_control) %notin% antioxidants.venn.list$recovery_vs_control
               ])
      
    # Query masigpro clusters to see if any of the genes from the DESeq2 contrasts were found significant and clustered by masigpro 
      antioxidants.in.masigpro <- data.frame(sigs.masigpro[sigs.masigpro$Gene %in% 
                                                            unique(c(antioxidants.venn.list$hypoxia_vs_control, 
                                                                     antioxidants.venn.list$pcrit_vs_control, 
                                                                     antioxidants.venn.list$anoxia_vs_control, 
                                                                     antioxidants.venn.list$recovery_vs_control,
                                                                     antioxidants.venn.list$Anoxia_vs_recovery,
                                                                     antioxidants.venn.list$Hypoxia_vs_anoxia)), c(32,31)
      ])
    # Grab gene annotations from combined_sigs frame
      antioxidants.in.masigpro <- merge(antioxidants.in.masigpro, combined_sigs[,1:2], by = "Gene", all = FALSE)
      
      table(antioxidants.in.masigpro$clusters) #Tally genes in each cluster
      
      
      
### Glycolysis and related pathways genes
      
    # split glycolysis frame into individual processes
      glycolysis.only <- unique(glycolysis[glycolysis$Process == "Glycolysis", ])
      tca.only <- unique(glycolysis[glycolysis$Process == "TCA", ])
      pentose.only <- unique(glycolysis[glycolysis$Process == "Pentose_Phosphate", ])
      pyruvate.only <- unique(glycolysis[glycolysis$Process == "Pyruvate", ])
      starch.only <- unique(glycolysis[glycolysis$Process == "Starch_Sucrose", ])
      fructose.only <- unique(glycolysis[glycolysis$Process == "Fructose_Mannose", ])
      
     
   #Glycolysis only 
    # Sum the number of glycolysis genes in our significant gene lists
      sum(combined_sigs$Gene %in% glycolysis.only$Gene)
      combined_sigs[combined_sigs$Gene %in% glycolysis.only$Gene, c(1,2)] #Which genes and description
      
    # How many are in the masigpro list and which clusters?
      sum(sigs.masigpro$Gene %in% glycolysis.only$Gene)  
      sigs.masigpro[sigs.masigpro$Gene %in% glycolysis.only$Gene, c(32,31)]
      glycolysis.only[glycolysis.only$Gene %in% sigs.masigpro[sigs.masigpro$Gene %in% glycolysis.only$Gene, 32], c(7,4)]
      
    # How many are in the hypoxia vs control DEG list?
      sum(Sigs.hypox.vs.cont$Gene %in% glycolysis.only$Gene)  
      Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% glycolysis.only$Gene, c(1,3,7)]
      
    # How many are in the pcrit vs control DEG list?
      sum(Sigs.pcrit.vs.cont$Gene %in% glycolysis.only$Gene)  
      Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% glycolysis.only$Gene, c(1,3,7)]
      
    # How many are in the anoxia vs control DEG list?
      sum(Sigs.anox.vs.control$Gene %in% glycolysis.only$Gene)  
      Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% glycolysis.only$Gene, c(1,3,7)]
      
    # How many are in the recovery vs control DEG list?
      sum(Sigs.recov.vs.cont$Gene %in% glycolysis.only$Gene)  
      Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% glycolysis.only$Gene, c(1,3,7)]
      
      
    # Make simple bar plot of the number of genes in response to oxygen compound GO term
    #Make data frame for plotting
      glycolysis.totals <- data.frame(times = c("Hypoxia vs control","Pcrit vs control","Anoxia vs control",
                                                  "Recovery vs control"),
                                        sigs = c(1,3,12,7),
                                        total = rep(150,4)
      )
    # Reorder factor to plot how we want
      glycolysis.totals$times <- factor(glycolysis.totals$times, 
                                          levels = c("Recovery vs control","Anoxia vs control","Pcrit vs control",
                                                     "Hypoxia vs control"))
    # Create and export figure
      jpeg(filename = "Glycolysis genes DESeq2 totals.jpg", width = 6.5, height = 4, units = "in", res = 300)
      glycolysis.totals |>
        ggplot() +
        geom_segment(aes(x = 0, xend = sigs, y = times, yend = times), linewidth = 2)+
        geom_point(aes(x = sigs, y = times), size = 9, col = "darkblue")+
        geom_text(aes(x = sigs, y = times, label = sigs), color = "white") +
        scale_x_continuous(expand = expansion(mult = c(0, 0.3)), breaks = seq(0, 15, by = 3))+
        xlab("Number of genes") + ylab("") +
        theme(axis.text.x = element_text(size = 12),
              axis.title.y = element_blank(),
              title = element_text(size = 9),
              panel.grid.major.x = element_line(color = "lightgrey"))+
        ggtitle("Significant DEG's out of 61 glycolysis genes")
      dev.off()
      
    # How many significant genes are shared between the each time point's RTO GO term?
    # Create venn diagram first
    # Make list of genes in each group for venn diagram of time course programs
    # Deseq comparison
      glycolysis.venn.list <- list(
        hypoxia_vs_control = Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% glycolysis.only$Gene, 1],
        pcrit_vs_control = Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% glycolysis.only$Gene, 1],
        anoxia_vs_control = Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% glycolysis.only$Gene, 1],
        recovery_vs_control =  Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% glycolysis.only$Gene, 1],
        Anoxia_vs_recovery = Sigs.anox.vs.recov[Sigs.anox.vs.recov$Gene %in% glycolysis.only$Gene, 1],
        Hypoxia_vs_anoxia = Sigs.hypox.vs.anox[Sigs.hypox.vs.anox$Gene %in% glycolysis.only$Gene, 1]
      )
      
      
    # Draw venn diagram of all five groups using venn package  
      jpeg(filename = "Venn Diagram Glycolysis Genes.jpg", width = 8.5, height = 5.5, units = "in", res = 300)
      ggVennDiagram(glycolysis.venn.list[c(1,2,3,4)], label_percent_digit = 1, label_size = 4) +
        scale_x_continuous(expand = expansion(mult = .2)) +
        scale_fill_distiller(name = "Counts:", palette = "GnBu", direction = 1) + 
        theme(legend.position = "bottom", legend.direction = "horizontal", 
              plot.margin = unit(c(0,0,0,0), "cm"))
      dev.off()
      
    # Pull out genes unique to each group
    # NOTE: %notin% function is a custom function assigned at top of script
    # Genes unique to hypoxia vs control
      glycolysis.venn.list$hypoxia_vs_control[glycolysis.venn.list$hypoxia_vs_control %notin% 
                                                  c(glycolysis.venn.list$pcrit_vs_control,
                                                    glycolysis.venn.list$anoxia_vs_control,
                                                    glycolysis.venn.list$recovery_vs_control,
                                                    glycolysis.venn.list$Anoxia_vs_recovery,
                                                    glycolysis.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to pcrit vs control
      glycolysis.venn.list$pcrit_vs_control[glycolysis.venn.list$pcrit_vs_control %notin% c(glycolysis.venn.list$hypoxia_vs_control,
                                                                                                glycolysis.venn.list$anoxia_vs_control,
                                                                                                glycolysis.venn.list$recovery_vs_control,
                                                                                                glycolysis.venn.list$Anoxia_vs_recovery,
                                                                                                glycolysis.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to anoxia vs control
      glycolysis.venn.list$anoxia_vs_control[glycolysis.venn.list$anoxia_vs_control %notin% c(glycolysis.venn.list$pcrit_vs_control,
                                                                                                  glycolysis.venn.list$hypoxia_vs_control,
                                                                                                  glycolysis.venn.list$recovery_vs_control,
                                                                                                  glycolysis.venn.list$Anoxia_vs_recovery,
                                                                                                  glycolysis.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to recovery vs control
      glycolysis.venn.list$recovery_vs_control[glycolysis.venn.list$recovery_vs_control %notin% c(glycolysis.venn.list$pcrit_vs_control,
                                                                                                      glycolysis.venn.list$anoxia_vs_control,
                                                                                                      glycolysis.venn.list$hypoxia_vs_control,
                                                                                                      glycolysis.venn.list$Anoxia_vs_recovery,
                                                                                                      glycolysis.venn.list$Hypoxia_vs_anoxia)
      ]
      
      
    # Pull out genes shared among all groups 
      glycolysis.venn.list$hypoxia_vs_control[glycolysis.venn.list$hypoxia_vs_control %in% glycolysis.venn.list$pcrit_vs_control &
                                                  glycolysis.venn.list$hypoxia_vs_control %in% glycolysis.venn.list$anoxia_vs_control &
                                                  glycolysis.venn.list$hypoxia_vs_control %in% glycolysis.venn.list$recovery_vs_control &
                                                  glycolysis.venn.list$hypoxia_vs_control %in% glycolysis.venn.list$Anoxia_vs_recovery &
                                                  glycolysis.venn.list$hypoxia_vs_control %in% glycolysis.venn.list$Hypoxia_vs_anoxia 
      ]
      
    # Pull out all genes that are in hypoxia, pcrit, and anoxia but NOT recovery
      unique(c(glycolysis.venn.list$hypoxia_vs_control, 
               glycolysis.venn.list$pcrit_vs_control, 
               glycolysis.venn.list$anoxia_vs_control)[c(glycolysis.venn.list$hypoxia_vs_control, 
                                                           glycolysis.venn.list$pcrit_vs_control, 
                                                           glycolysis.venn.list$anoxia_vs_control) %notin% glycolysis.venn.list$recovery_vs_control
               ])
      
    # Query masigpro clusters to see if any of the genes from the DESeq2 contrasts were found significant and clustered by masigpro 
      glycolysis.in.masigpro <- data.frame(sigs.masigpro[sigs.masigpro$Gene %in% 
                                                             unique(c(glycolysis.venn.list$hypoxia_vs_control, 
                                                                      glycolysis.venn.list$pcrit_vs_control, 
                                                                      glycolysis.venn.list$anoxia_vs_control, 
                                                                      glycolysis.venn.list$recovery_vs_control,
                                                                      glycolysis.venn.list$Anoxia_vs_recovery,
                                                                      glycolysis.venn.list$Hypoxia_vs_anoxia)), c(32,31)
      ])
    # Grab gene annotations from combined_sigs frame
      glycolysis.in.masigpro <- merge(glycolysis.in.masigpro, combined_sigs[,1:2], by = "Gene", all = FALSE)
      
      table(glycolysis.in.masigpro$clusters) #Tally genes in each cluster
      
    
      
    #Generate pathway
      glyco.pathway <- highlight_entities("tcf00010", c("131878644", "131880287", "131879670", "131877236",
                                                        "131879841", "131888003", "131883197", "131893264", 	
                                                        "131877162", "131878133", "131887499", "131889496"), 
                                          fill_color = "paleturquoise", legend_name="Significant")
      glyco.pathway
      
    #Export as png
      ggkeggsave(filename="glyco pathway.png", glyco.pathway, dpi=300)
      
  ### Generate heatmap of just the sig genes from the expression averaged wide dataframe
      glyco.masigpro<- as.data.frame(expression.averaged.wide[expression.averaged.wide$Gene %in% combined_sigs[combined_sigs$Gene %in% glycolysis.only$Gene, 1],])
      row.names(glyco.masigpro) <- glyco.masigpro$Gene
      glyco.masigpro <- glyco.masigpro[,-c(1)]
      
      
      jpeg(file = "Glyco Sigs Heatmap.jpg", width = 5, height = 9, units= "in", res = 300)
      pheatmap(glyco.masigpro, cluster_rows = F, cluster_cols = F,
                                  show_rownames=T, border_color=NA,
                                scale = "row", 
                                labels_col = c("Normoxia", "Mild Hypox", "Pcrit", "Anoxia", "Recovery"),
                                fontsize = 9,fontsize_row = 10, height=20, main = "Glycolysis")
      dev.off()
      
      
      
   # TCA only
      
      # Sum the number of antioxidant genes in our significant gene lists
      sum(combined_sigs$Gene %in% tca.only$Gene)
      combined_sigs[combined_sigs$Gene %in% tca.only$Gene, c(1,2)] #Which genes and description
      
      # How many are in the masigpro list and which clusters?
      sum(sigs.masigpro$Gene %in% tca.only$Gene)  
      sigs.masigpro[sigs.masigpro$Gene %in% tca.only$Gene, c(32,31)]
      
      # How many are in the hypoxia vs control DEG list?
      sum(Sigs.hypox.vs.cont$Gene %in% tca.only$Gene)  
      Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% tca.only$Gene, c(1,3,7)]
      
      # How many are in the pcrit vs control DEG list?
      sum(Sigs.pcrit.vs.cont$Gene %in% tca.only$Gene)  
      Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% tca.only$Gene, c(1,3,7)]
      
      # How many are in the anoxia vs control DEG list?
      sum(Sigs.anox.vs.control$Gene %in% tca.only$Gene)  
      Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% tca.only$Gene, c(1,3,7)]
      
      # How many are in the recovery vs control DEG list?
      sum(Sigs.recov.vs.cont$Gene %in% tca.only$Gene)  
      Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% tca.only$Gene, c(1,3,7)]
      
      
      # Make simple bar plot of the number of genes in response to oxygen compound GO term
      #Make data frame for plotting
      tca.totals <- data.frame(times = c("Hypoxia vs control","Pcrit vs control","Anoxia vs control",
                                                "Recovery vs control"),
                                      sigs = c(0,1,2,1),
                                      total = rep(28,4)
      )
      # Reorder factor to plot how we want
      tca.totals$times <- factor(tca.totals$times, 
                                        levels = c("Recovery vs control","Anoxia vs control","Pcrit vs control",
                                                   "Hypoxia vs control"))
      # Create and export figure
      jpeg(filename = "TCA genes DESeq2 totals.jpg", width = 6.5, height = 4, units = "in", res = 300)
      tca.totals |>
        ggplot() +
        geom_segment(aes(x = 0, xend = sigs, y = times, yend = times), linewidth = 2)+
        geom_point(aes(x = sigs, y = times), size = 9, col = "darkblue")+
        geom_text(aes(x = sigs, y = times, label = sigs), color = "white") +
        scale_x_continuous(expand = expansion(mult = c(0.1, 0.5)), breaks = seq(0, 3, by = 1))+
        xlab("Number of genes") + ylab("") +
        theme(axis.text.x = element_text(size = 12),
              axis.title.y = element_blank(),
              title = element_text(size = 9),
              panel.grid.major.x = element_line(color = "lightgrey"))+
        ggtitle("Significant DEG's out of 28 TCA genes")
      dev.off()
      
      # How many significant genes are shared between the each time point's RTO GO term?
      # Create venn diagram first
      # Make list of genes in each group for venn diagram of time course programs
      # Deseq comparison
      # Deseq comparison
        tca.venn.list <- list(
          hypoxia_vs_control = Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% tca.only$Gene, 1],
          pcrit_vs_control = Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% tca.only$Gene, 1],
          anoxia_vs_control = Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% tca.only$Gene, 1],
          recovery_vs_control =  Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% tca.only$Gene, 1],
          Anoxia_vs_recovery = Sigs.anox.vs.recov[Sigs.anox.vs.recov$Gene %in% tca.only$Gene, 1],
          Hypoxia_vs_anoxia = Sigs.hypox.vs.anox[Sigs.hypox.vs.anox$Gene %in% tca.only$Gene, 1]
        )
      
      
      # Draw venn diagram of all five groups using venn package  
      jpeg(filename = "Venn Diagram TCA Genes.jpg", width = 8.5, height = 5.5, units = "in", res = 300)
      ggVennDiagram(tca.venn.list[c(1,2,3,4)], label_percent_digit = 1, label_size = 4) +
        scale_x_continuous(expand = expansion(mult = .2)) +
        scale_fill_distiller(name = "Counts:", palette = "GnBu", direction = 1) + 
        theme(legend.position = "bottom", legend.direction = "horizontal", 
              plot.margin = unit(c(0,0,0,0), "cm"))
      dev.off()
      
    # Pull out genes unique to each group
    # NOTE: %notin% function is a custom function assigned at top of script
    # Genes unique to hypoxia vs control
      tca.venn.list$hypoxia_vs_control[tca.venn.list$hypoxia_vs_control %notin% 
                                                c(tca.venn.list$pcrit_vs_control,
                                                  tca.venn.list$anoxia_vs_control,
                                                  tca.venn.list$recovery_vs_control,
                                                  tca.venn.list$Anoxia_vs_recovery,
                                                  tca.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to pcrit vs control
      tca.venn.list$pcrit_vs_control[tca.venn.list$pcrit_vs_control %notin% c(tca.venn.list$hypoxia_vs_control,
                                                                                            tca.venn.list$anoxia_vs_control,
                                                                                            tca.venn.list$recovery_vs_control,
                                                                                            tca.venn.list$Anoxia_vs_recovery,
                                                                                            tca.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to anoxia vs control
      tca.venn.list$anoxia_vs_control[tca.venn.list$anoxia_vs_control %notin% c(tca.venn.list$pcrit_vs_control,
                                                                                              tca.venn.list$hypoxia_vs_control,
                                                                                              tca.venn.list$recovery_vs_control,
                                                                                              tca.venn.list$Anoxia_vs_recovery,
                                                                                              tca.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to recovery vs control
      tca.venn.list$recovery_vs_control[tca.venn.list$recovery_vs_control %notin% c(tca.venn.list$pcrit_vs_control,
                                                                                                  tca.venn.list$anoxia_vs_control,
                                                                                                  tca.venn.list$hypoxia_vs_control,
                                                                                                  tca.venn.list$Anoxia_vs_recovery,
                                                                                                  tca.venn.list$Hypoxia_vs_anoxia)
      ]
      
      
    # Pull out genes shared among all groups 
      tca.venn.list$hypoxia_vs_control[tca.venn.list$hypoxia_vs_control %in% tca.venn.list$pcrit_vs_control &
                                                tca.venn.list$hypoxia_vs_control %in% tca.venn.list$anoxia_vs_control &
                                                tca.venn.list$hypoxia_vs_control %in% tca.venn.list$recovery_vs_control &
                                                tca.venn.list$hypoxia_vs_control %in% tca.venn.list$Anoxia_vs_recovery &
                                                tca.venn.list$hypoxia_vs_control %in% tca.venn.list$Hypoxia_vs_anoxia 
      ]
      
    # Pull out all genes that are in hypoxia, pcrit, and anoxia but NOT recovery
      unique(c(tca.venn.list$hypoxia_vs_control, 
               tca.venn.list$pcrit_vs_control, 
               tca.venn.list$anoxia_vs_control)[c(tca.venn.list$hypoxia_vs_control, 
                                                         tca.venn.list$pcrit_vs_control, 
                                                         tca.venn.list$anoxia_vs_control) %notin% tca.venn.list$recovery_vs_control
               ])
      
    # Query masigpro clusters to see if any of the genes from the DESeq2 contrasts were found significant and clustered by masigpro 
      tca.in.masigpro <- data.frame(sigs.masigpro[sigs.masigpro$Gene %in% 
                                                           unique(c(tca.venn.list$hypoxia_vs_control, 
                                                                    tca.venn.list$pcrit_vs_control, 
                                                                    tca.venn.list$anoxia_vs_control, 
                                                                    tca.venn.list$recovery_vs_control,
                                                                    tca.venn.list$Anoxia_vs_recovery,
                                                                    tca.venn.list$Hypoxia_vs_anoxia)), c(32,31)
      ])
    # Grab gene annotations from combined_sigs frame
      tca.in.masigpro <- merge(tca.in.masigpro, combined_sigs[,1:2], by = "Gene", all = FALSE)
      
      table(tca.in.masigpro$clusters) #Tally genes in each cluster
        
      
    #Generate pathway
      tca.pathway <- highlight_entities("tcf00020", c("131884131", "131888003"), 
                                          fill_color = "paleturquoise", legend_name="Significant")
      tca.pathway
      
    #Export as png
      ggkeggsave(filename="tca pathway.png", tca.pathway, dpi=900)
      
    ### Generate heatmap of just the sig genes from the expression averaged wide dataframe
      tca.masigpro<- as.data.frame(expression.averaged.wide[expression.averaged.wide$Gene %in% combined_sigs[combined_sigs$Gene %in% tca.only$Gene, 1],])
      row.names(tca.masigpro) <- tca.masigpro$Gene
      tca.masigpro <- tca.masigpro[,-c(1)]
      
      
      jpeg(file = "TCA Sigs Heatmap.jpg", width = 5, height = 2, units= "in", res = 300)
      pheatmap(tca.masigpro, cluster_rows = F, cluster_cols = F,
               show_rownames=T, border_color=NA,
               scale = "row", 
               labels_col = c("Normoxia", "Mild Hypox", "Pcrit", "Anoxia", "Recovery"),
               fontsize = 9,fontsize_row = 10, height=20, main = "TCA Cycle")
      dev.off()
      
  # Pyruvate pathway only
      
    # Sum the number of antioxidant genes in our significant gene lists
      sum(combined_sigs$Gene %in% pyruvate.only$Gene)
      combined_sigs[combined_sigs$Gene %in% pyruvate.only$Gene, c(1,2)] #Which genes and description
      
    # How many are in the masigpro list and which clusters?
      sum(sigs.masigpro$Gene %in% pyruvate.only$Gene)  
      sigs.masigpro[sigs.masigpro$Gene %in% pyruvate.only$Gene, c(32,31)]
      
    # How many are in the hypoxia vs control DEG list?
      sum(Sigs.hypox.vs.cont$Gene %in% pyruvate.only$Gene)  
      Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% pyruvate.only$Gene, c(1,3,7)]
      
    # How many are in the pcrit vs control DEG list?
      sum(Sigs.pcrit.vs.cont$Gene %in% pyruvate.only$Gene)  
      Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% pyruvate.only$Gene, c(1,3,7)]
      
    # How many are in the anoxia vs control DEG list?
      sum(Sigs.anox.vs.control$Gene %in% pyruvate.only$Gene)  
      Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% pyruvate.only$Gene, c(1,3,7)]
      
    # How many are in the recovery vs control DEG list?
      sum(Sigs.recov.vs.cont$Gene %in% pyruvate.only$Gene)  
      Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% pyruvate.only$Gene, c(1,3,7)]
      
      
    # Make simple bar plot of the number of genes in response to oxygen compound GO term
    #Make data frame for plotting
      pyruvate.totals <- data.frame(times = c("Hypoxia vs control","Pcrit vs control","Anoxia vs control",
                                         "Recovery vs control"),
                               sigs = c(1,3,8,7),
                               total = rep(28,4)
      )
    # Reorder factor to plot how we want
      pyruvate.totals$times <- factor(pyruvate.totals$times, 
                                 levels = c("Recovery vs control","Anoxia vs control","Pcrit vs control",
                                            "Hypoxia vs control"))
    # Create and export figure
      jpeg(filename = "Pyruvate pathway genes DESeq2 totals.jpg", width = 6.5, height = 4, units = "in", res = 300)
      pyruvate.totals |>
        ggplot() +
        geom_segment(aes(x = 0, xend = sigs, y = times, yend = times), linewidth = 2)+
        geom_point(aes(x = sigs, y = times), size = 9, col = "darkblue")+
        geom_text(aes(x = sigs, y = times, label = sigs), color = "white") +
        scale_x_continuous(expand = expansion(mult = c(0.1, 0.3)), breaks = seq(0, 10, by = 2))+
        xlab("Number of genes") + ylab("") +
        theme(axis.text.x = element_text(size = 12),
              axis.title.y = element_blank(),
              title = element_text(size = 9),
              panel.grid.major.x = element_line(color = "lightgrey"))+
        ggtitle("Significant DEG's out of 54 pyruvate genes")
      dev.off()
      
    # How many significant genes are shared between the each time point's RTO GO term?
    # Create venn diagram first
    # Make list of genes in each group for venn diagram of time course programs
    # Deseq comparison
      pyruvate.venn.list <- list(
        hypoxia_vs_control = Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% pyruvate.only$Gene, 1],
        pcrit_vs_control = Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% pyruvate.only$Gene, 1],
        anoxia_vs_control = Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% pyruvate.only$Gene, 1],
        recovery_vs_control =  Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% pyruvate.only$Gene, 1],
        Anoxia_vs_recovery = Sigs.anox.vs.recov[Sigs.anox.vs.recov$Gene %in% pyruvate.only$Gene, 1],
        Hypoxia_vs_anoxia = Sigs.hypox.vs.anox[Sigs.hypox.vs.anox$Gene %in% pyruvate.only$Gene, 1]
      )
      
      
    # Draw venn diagram of all five groups using venn package  
      jpeg(filename = "Venn Diagram Pyruvate Genes.jpg", width = 8.5, height = 5.5, units = "in", res = 300)
      ggVennDiagram(pyruvate.venn.list[c(1,2,3,4)], label_percent_digit = 1, label_size = 4) +
        scale_x_continuous(expand = expansion(mult = .2)) +
        scale_fill_distiller(name = "Counts:", palette = "GnBu", direction = 1) + 
        theme(legend.position = "bottom", legend.direction = "horizontal", 
              plot.margin = unit(c(0,0,0,0), "cm"))
      dev.off()
      
    # Pull out genes unique to each group
    # NOTE: %notin% function is a custom function assigned at top of script
    # Genes unique to hypoxia vs control
      pyruvate.venn.list$hypoxia_vs_control[pyruvate.venn.list$hypoxia_vs_control %notin% 
                                         c(pyruvate.venn.list$pcrit_vs_control,
                                           pyruvate.venn.list$anoxia_vs_control,
                                           pyruvate.venn.list$recovery_vs_control,
                                           pyruvate.venn.list$Anoxia_vs_recovery,
                                           pyruvate.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to pcrit vs control
      pyruvate.venn.list$pcrit_vs_control[pyruvate.venn.list$pcrit_vs_control %notin% c(pyruvate.venn.list$hypoxia_vs_control,
                                                                              pyruvate.venn.list$anoxia_vs_control,
                                                                              pyruvate.venn.list$recovery_vs_control,
                                                                              pyruvate.venn.list$Anoxia_vs_recovery,
                                                                              pyruvate.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to anoxia vs control
      pyruvate.venn.list$anoxia_vs_control[pyruvate.venn.list$anoxia_vs_control %notin% c(pyruvate.venn.list$pcrit_vs_control,
                                                                                pyruvate.venn.list$hypoxia_vs_control,
                                                                                pyruvate.venn.list$recovery_vs_control,
                                                                                pyruvate.venn.list$Anoxia_vs_recovery,
                                                                                pyruvate.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to recovery vs control
      pyruvate.venn.list$recovery_vs_control[pyruvate.venn.list$recovery_vs_control %notin% c(pyruvate.venn.list$pcrit_vs_control,
                                                                                    pyruvate.venn.list$anoxia_vs_control,
                                                                                    pyruvate.venn.list$hypoxia_vs_control,
                                                                                    pyruvate.venn.list$Anoxia_vs_recovery,
                                                                                    pyruvate.venn.list$Hypoxia_vs_anoxia)
      ]
      
      
    # Pull out genes shared among all groups 
      pyruvate.venn.list$hypoxia_vs_control[pyruvate.venn.list$hypoxia_vs_control %in% pyruvate.venn.list$pcrit_vs_control &
                                         pyruvate.venn.list$hypoxia_vs_control %in% pyruvate.venn.list$anoxia_vs_control &
                                         pyruvate.venn.list$hypoxia_vs_control %in% pyruvate.venn.list$recovery_vs_control &
                                         pyruvate.venn.list$hypoxia_vs_control %in% pyruvate.venn.list$Anoxia_vs_recovery &
                                         pyruvate.venn.list$hypoxia_vs_control %in% pyruvate.venn.list$Hypoxia_vs_anoxia 
      ]
      
    # Pull out all genes that are in hypoxia, pcrit, and anoxia but NOT recovery
      unique(c(pyruvate.venn.list$hypoxia_vs_control, 
               pyruvate.venn.list$pcrit_vs_control, 
               pyruvate.venn.list$anoxia_vs_control)[c(pyruvate.venn.list$hypoxia_vs_control, 
                                                  pyruvate.venn.list$pcrit_vs_control, 
                                                  pyruvate.venn.list$anoxia_vs_control) %notin% pyruvate.venn.list$recovery_vs_control
               ])
      
    # Query masigpro clusters to see if any of the genes from the DESeq2 contrasts were found significant and clustered by masigpro 
      pyruvate.in.masigpro <- data.frame(sigs.masigpro[sigs.masigpro$Gene %in% 
                                                    unique(c(pyruvate.venn.list$hypoxia_vs_control, 
                                                             pyruvate.venn.list$pcrit_vs_control, 
                                                             pyruvate.venn.list$anoxia_vs_control, 
                                                             pyruvate.venn.list$recovery_vs_control,
                                                             pyruvate.venn.list$Anoxia_vs_recovery,
                                                             pyruvate.venn.list$Hypoxia_vs_anoxia)), c(32,31)
      ])
    # Grab gene annotations from combined_sigs frame
      pyruvate.in.masigpro <- merge(pyruvate.in.masigpro, combined_sigs[,1:2], by = "Gene", all = FALSE)
      
      table(pyruvate.in.masigpro$clusters) #Tally genes in each cluster
        
    #Generate pathway
      pyruvate.pathway <- highlight_entities("tcf00620", c("131878644", "131880287", "131879670",
                                                           "131881313", "131884131", "131888003"), 
                                        fill_color = "paleturquoise", legend_name="Significant")
      pyruvate.pathway
      
    #Export as png
      ggkeggsave(filename="pyruvate pathway.png", pyruvate.pathway, dpi=300)
      
    ### Generate heatmap of just the sig genes from the expression averaged wide dataframe
      pyruvate.masigpro<- as.data.frame(expression.averaged.wide[expression.averaged.wide$Gene %in% combined_sigs[combined_sigs$Gene %in% pyruvate.only$Gene, 1],])
      row.names(pyruvate.masigpro) <- pyruvate.masigpro$Gene
      pyruvate.masigpro <- pyruvate.masigpro[,-c(1)]
      
      
      jpeg(file = "pyruvate Sigs Heatmap.jpg", width = 5, height = 7, units= "in", res = 300)
      pheatmap(pyruvate.masigpro, cluster_rows = F, cluster_cols = F,
               show_rownames=T, border_color=NA,
               scale = "row", 
               labels_col = c("Normoxia", "Mild Hypox", "Pcrit", "Anoxia", "Recovery"),
               fontsize = 9,fontsize_row = 10, height=20, main = "Pyruvate Cycle")
      dev.off()
      
      
      
  # Pentose phosphate pathway only
      
    # Sum the number of antioxidant genes in our significant gene lists
      sum(combined_sigs$Gene %in% pentose.only$Gene)
      combined_sigs[combined_sigs$Gene %in% pentose.only$Gene, c(1,2)] #Which genes and description
      
    # How many are in the masigpro list and which clusters?
      sum(sigs.masigpro$Gene %in% pentose.only$Gene)  
      sigs.masigpro[sigs.masigpro$Gene %in% pentose.only$Gene, c(32,31)]
      
    # How many are in the hypoxia vs control DEG list?
      sum(Sigs.hypox.vs.cont$Gene %in% pentose.only$Gene)  
      Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% pentose.only$Gene, c(1,3,7)]
      
    # How many are in the pcrit vs control DEG list?
      sum(Sigs.pcrit.vs.cont$Gene %in% pentose.only$Gene)  
      Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% pentose.only$Gene, c(1,3,7)]
      
    # How many are in the anoxia vs control DEG list?
      sum(Sigs.anox.vs.control$Gene %in% pentose.only$Gene)  
      Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% pentose.only$Gene, c(1,3,7)]
      
    # How many are in the recovery vs control DEG list?
      sum(Sigs.recov.vs.cont$Gene %in% pentose.only$Gene)  
      Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% pentose.only$Gene, c(1,3,7)]
      
      
    # Make simple bar plot of the number of genes in response to oxygen compound GO term
    #Make data frame for plotting
      pentose.totals <- data.frame(times = c("Hypoxia vs control","Pcrit vs control","Anoxia vs control",
                                              "Recovery vs control"),
                                    sigs = c(0,1,4,5),
                                    total = rep(28,4)
      )
    # Reorder factor to plot how we want
      pentose.totals$times <- factor(pentose.totals$times, 
                                      levels = c("Recovery vs control","Anoxia vs control","Pcrit vs control",
                                                 "Hypoxia vs control"))
    # Create and export figure
      jpeg(filename = "Pentose phosphate pathway genes DESeq2 totals.jpg", width = 6.5, height = 4, units = "in", res = 300)
      pentose.totals |>
        ggplot() +
        geom_segment(aes(x = 0, xend = sigs, y = times, yend = times), linewidth = 2)+
        geom_point(aes(x = sigs, y = times), size = 9, col = "darkblue")+
        geom_text(aes(x = sigs, y = times, label = sigs), color = "white") +
        scale_x_continuous(expand = expansion(mult = c(0.1, 0.3)), breaks = seq(0, 10, by = 2))+
        xlab("Number of genes") + ylab("") +
        theme(axis.text.x = element_text(size = 12),
              axis.title.y = element_blank(),
              title = element_text(size = 9),
              panel.grid.major.x = element_line(color = "lightgrey"))+
        ggtitle("Significant DEG's out of 32 pentose-phosphate genes")
      dev.off()
      
    # How many significant genes are shared between the each time point's RTO GO term?
    # Create venn diagram first
    # Make list of genes in each group for venn diagram of time course programs
    # Deseq comparison
      pentose.venn.list <- list(
        hypoxia_vs_control = Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% pentose.only$Gene, 1],
        pcrit_vs_control = Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% pentose.only$Gene, 1],
        anoxia_vs_control = Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% pentose.only$Gene, 1],
        recovery_vs_control =  Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% pentose.only$Gene, 1],
        Anoxia_vs_recovery = Sigs.anox.vs.recov[Sigs.anox.vs.recov$Gene %in% pentose.only$Gene, 1],
        Hypoxia_vs_anoxia = Sigs.hypox.vs.anox[Sigs.hypox.vs.anox$Gene %in% pentose.only$Gene, 1]
      )
      
      
    # Draw venn diagram of all five groups using venn package  
      jpeg(filename = "Venn Diagram Pentose-phosphate Genes.jpg", width = 8.5, height = 5.5, units = "in", res = 300)
      ggVennDiagram(pentose.venn.list[c(1,2,3,4)], label_percent_digit = 1, label_size = 4) +
        scale_x_continuous(expand = expansion(mult = .2)) +
        scale_fill_distiller(name = "Counts:", palette = "GnBu", direction = 1) + 
        theme(legend.position = "bottom", legend.direction = "horizontal", 
              plot.margin = unit(c(0,0,0,0), "cm"))
      dev.off()
      
    # Pull out genes unique to each group
    # NOTE: %notin% function is a custom function assigned at top of script
    # Genes unique to hypoxia vs control
      pentose.venn.list$hypoxia_vs_control[pentose.venn.list$hypoxia_vs_control %notin% 
                                              c(pentose.venn.list$pcrit_vs_control,
                                                pentose.venn.list$anoxia_vs_control,
                                                pentose.venn.list$recovery_vs_control,
                                                pentose.venn.list$Anoxia_vs_recovery,
                                                pentose.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to pcrit vs control
      pentose.venn.list$pcrit_vs_control[pentose.venn.list$pcrit_vs_control %notin% c(pentose.venn.list$hypoxia_vs_control,
                                                                                        pentose.venn.list$anoxia_vs_control,
                                                                                        pentose.venn.list$recovery_vs_control,
                                                                                        pentose.venn.list$Anoxia_vs_recovery,
                                                                                        pentose.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to anoxia vs control
      pentose.venn.list$anoxia_vs_control[pentose.venn.list$anoxia_vs_control %notin% c(pentose.venn.list$pcrit_vs_control,
                                                                                          pentose.venn.list$hypoxia_vs_control,
                                                                                          pentose.venn.list$recovery_vs_control,
                                                                                          pentose.venn.list$Anoxia_vs_recovery,
                                                                                          pentose.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to recovery vs control
      pentose.venn.list$recovery_vs_control[pentose.venn.list$recovery_vs_control %notin% c(pentose.venn.list$pcrit_vs_control,
                                                                                              pentose.venn.list$anoxia_vs_control,
                                                                                              pentose.venn.list$hypoxia_vs_control,
                                                                                              pentose.venn.list$Anoxia_vs_recovery,
                                                                                              pentose.venn.list$Hypoxia_vs_anoxia)
      ]
      
      
    # Pull out genes shared among all groups 
      pentose.venn.list$hypoxia_vs_control[pentose.venn.list$hypoxia_vs_control %in% pentose.venn.list$pcrit_vs_control &
                                              pentose.venn.list$hypoxia_vs_control %in% pentose.venn.list$anoxia_vs_control &
                                              pentose.venn.list$hypoxia_vs_control %in% pentose.venn.list$recovery_vs_control &
                                              pentose.venn.list$hypoxia_vs_control %in% pentose.venn.list$Anoxia_vs_recovery &
                                              pentose.venn.list$hypoxia_vs_control %in% pentose.venn.list$Hypoxia_vs_anoxia 
      ]
      
    # Pull out all genes that are in hypoxia, pcrit, and anoxia but NOT recovery
      unique(c(pentose.venn.list$hypoxia_vs_control, 
               pentose.venn.list$pcrit_vs_control, 
               pentose.venn.list$anoxia_vs_control)[c(pentose.venn.list$hypoxia_vs_control, 
                                                       pentose.venn.list$pcrit_vs_control, 
                                                       pentose.venn.list$anoxia_vs_control) %notin% pentose.venn.list$recovery_vs_control
               ])
      
    # Query masigpro clusters to see if any of the genes from the DESeq2 contrasts were found significant and clustered by masigpro 
      pentose.in.masigpro <- data.frame(sigs.masigpro[sigs.masigpro$Gene %in% 
                                                         unique(c(pentose.venn.list$hypoxia_vs_control, 
                                                                  pentose.venn.list$pcrit_vs_control, 
                                                                  pentose.venn.list$anoxia_vs_control, 
                                                                  pentose.venn.list$recovery_vs_control,
                                                                  pentose.venn.list$Anoxia_vs_recovery,
                                                                  pentose.venn.list$Hypoxia_vs_anoxia)), c(32,31)
      ])
    # Grab gene annotations from combined_sigs frame
      pentose.in.masigpro <- merge(pentose.in.masigpro, combined_sigs[,1:2], by = "Gene", all = FALSE)
      
      table(pentose.in.masigpro$clusters) #Tally genes in each cluster
      
    #Generate pathway
      pentose.pathway <- highlight_entities("tcf00030", c("131877162", "131887117", "131878396", "131883400", "131892879",
                                                          "131878020", "131887499", "100576547", "131888855"), 
                                             fill_color = "paleturquoise", legend_name="Significant")
      pentose.pathway
      
    #Export as png
      ggkeggsave(filename="pentose pathway.png", pentose.pathway, dpi=900)
      
    ### Generate heatmap of just the sig genes from the expression averaged wide dataframe
      pentose.masigpro<- as.data.frame(expression.averaged.wide[expression.averaged.wide$Gene %in% combined_sigs[combined_sigs$Gene %in% pentose.only$Gene, 1],])
      row.names(pentose.masigpro) <- pentose.masigpro$Gene
      pentose.masigpro <- pentose.masigpro[,-c(1)]
      
      
      jpeg(file = "pentose Sigs Heatmap.jpg", width = 5, height = 7, units= "in", res = 300)
      pheatmap(pentose.masigpro, cluster_rows = F, cluster_cols = F,
               show_rownames=T, border_color=NA,
               scale = "row", 
               labels_col = c("Normoxia", "Mild Hypox", "Pcrit", "Anoxia", "Recovery"),
               fontsize = 9,fontsize_row = 10, height=20, main = "pentose Cycle")
      dev.off()
      
  # Starch sucrose pathway only
      
    # Sum the number of antioxidant genes in our significant gene lists
      sum(combined_sigs$Gene %in% starch.only$Gene)
      combined_sigs[combined_sigs$Gene %in% starch.only$Gene, c(1,2)] #Which genes and description
      
    # How many are in the masigpro list and which clusters?
      sum(sigs.masigpro$Gene %in% starch.only$Gene)  
      sigs.masigpro[sigs.masigpro$Gene %in% starch.only$Gene, c(32,31)]
      
    # How many are in the hypoxia vs control DEG list?
      sum(Sigs.hypox.vs.cont$Gene %in% starch.only$Gene)  
      Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% starch.only$Gene, c(1,3,7)]
      
    # How many are in the pcrit vs control DEG list?
      sum(Sigs.pcrit.vs.cont$Gene %in% starch.only$Gene)  
      Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% starch.only$Gene, c(1,3,7)]
      
    # How many are in the anoxia vs control DEG list?
      sum(Sigs.anox.vs.control$Gene %in% starch.only$Gene)  
      Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% starch.only$Gene, c(1,3,7)]
      
    # How many are in the recovery vs control DEG list?
      sum(Sigs.recov.vs.cont$Gene %in% starch.only$Gene)  
      Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% starch.only$Gene, c(1,3,7)]
      
      
    # Make simple bar plot of the number of genes in response to oxygen compound GO term
    #Make data frame for plotting
      starch.totals <- data.frame(times = c("Hypoxia vs control","Pcrit vs control","Anoxia vs control",
                                             "Recovery vs control"),
                                   sigs = c(4,9,13,10),
                                   total = rep(76,4)
      )
    # Reorder factor to plot how we want
      starch.totals$times <- factor(starch.totals$times, 
                                     levels = c("Recovery vs control","Anoxia vs control","Pcrit vs control",
                                                "Hypoxia vs control"))
    # Create and export figure
      jpeg(filename = "Starch-sucrose pathway genes DESeq2 totals.jpg", width = 6.5, height = 4, units = "in", res = 300)
      starch.totals |>
        ggplot() +
        geom_segment(aes(x = 0, xend = sigs, y = times, yend = times), linewidth = 2)+
        geom_point(aes(x = sigs, y = times), size = 9, col = "darkblue")+
        geom_text(aes(x = sigs, y = times, label = sigs), color = "white") +
        scale_x_continuous(expand = expansion(mult = c(0, 0.3)), breaks = seq(0, 10, by = 2))+
        xlab("Number of genes") + ylab("") +
        theme(axis.text.x = element_text(size = 12),
              axis.title.y = element_blank(),
              title = element_text(size = 9),
              panel.grid.major.x = element_line(color = "lightgrey"))+
        ggtitle("Significant DEG's out of 38 starch-sucrose genes")
      dev.off()
      
    # How many significant genes are shared between the each time point's RTO GO term?
    # Create venn diagram first
    # Make list of genes in each group for venn diagram of time course programs
    # Deseq comparison
      starch.venn.list <- list(
        hypoxia_vs_control = Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% starch.only$Gene, 1],
        pcrit_vs_control = Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% starch.only$Gene, 1],
        anoxia_vs_control = Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% starch.only$Gene, 1],
        recovery_vs_control =  Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% starch.only$Gene, 1],
        Anoxia_vs_recovery = Sigs.anox.vs.recov[Sigs.anox.vs.recov$Gene %in% starch.only$Gene, 1],
        Hypoxia_vs_anoxia = Sigs.hypox.vs.anox[Sigs.hypox.vs.anox$Gene %in% starch.only$Gene, 1]
      )
      
      
    # Draw venn diagram of all five groups using venn package  
      jpeg(filename = "Venn Diagram Pentose-phosphate Genes.jpg", width = 8.5, height = 5.5, units = "in", res = 300)
      ggVennDiagram(starch.venn.list[c(1,2,3,4)], label_percent_digit = 1, label_size = 4) +
        scale_x_continuous(expand = expansion(mult = .2)) +
        scale_fill_distiller(name = "Counts:", palette = "GnBu", direction = 1) + 
        theme(legend.position = "bottom", legend.direction = "horizontal", 
              plot.margin = unit(c(0,0,0,0), "cm"))
      dev.off()
      
    # Pull out genes unique to each group
    # NOTE: %notin% function is a custom function assigned at top of script
    # Genes unique to hypoxia vs control
      starch.venn.list$hypoxia_vs_control[starch.venn.list$hypoxia_vs_control %notin% 
                                                c(starch.venn.list$pcrit_vs_control,
                                                  starch.venn.list$anoxia_vs_control,
                                                  starch.venn.list$recovery_vs_control,
                                                  starch.venn.list$Anoxia_vs_recovery,
                                                  starch.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to pcrit vs control
      starch.venn.list$pcrit_vs_control[starch.venn.list$pcrit_vs_control %notin% c(starch.venn.list$hypoxia_vs_control,
                                                                                            starch.venn.list$anoxia_vs_control,
                                                                                            starch.venn.list$recovery_vs_control,
                                                                                            starch.venn.list$Anoxia_vs_recovery,
                                                                                            starch.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to anoxia vs control
      starch.venn.list$anoxia_vs_control[starch.venn.list$anoxia_vs_control %notin% c(starch.venn.list$pcrit_vs_control,
                                                                                              starch.venn.list$hypoxia_vs_control,
                                                                                              starch.venn.list$recovery_vs_control,
                                                                                              starch.venn.list$Anoxia_vs_recovery,
                                                                                              starch.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to recovery vs control
      starch.venn.list$recovery_vs_control[starch.venn.list$recovery_vs_control %notin% c(starch.venn.list$pcrit_vs_control,
                                                                                                  starch.venn.list$anoxia_vs_control,
                                                                                                  starch.venn.list$hypoxia_vs_control,
                                                                                                  starch.venn.list$Anoxia_vs_recovery,
                                                                                                  starch.venn.list$Hypoxia_vs_anoxia)
      ]
      
      
    # Pull out genes shared among all groups 
      starch.venn.list$hypoxia_vs_control[starch.venn.list$hypoxia_vs_control %in% starch.venn.list$pcrit_vs_control &
                                                starch.venn.list$hypoxia_vs_control %in% starch.venn.list$anoxia_vs_control &
                                                starch.venn.list$hypoxia_vs_control %in% starch.venn.list$recovery_vs_control &
                                                starch.venn.list$hypoxia_vs_control %in% starch.venn.list$Anoxia_vs_recovery &
                                                starch.venn.list$hypoxia_vs_control %in% starch.venn.list$Hypoxia_vs_anoxia 
      ]
      
    # Pull out all genes that are in hypoxia, pcrit, and anoxia but NOT recovery
      unique(c(starch.venn.list$hypoxia_vs_control, 
               starch.venn.list$pcrit_vs_control, 
               starch.venn.list$anoxia_vs_control)[c(starch.venn.list$hypoxia_vs_control, 
                                                         starch.venn.list$pcrit_vs_control, 
                                                         starch.venn.list$anoxia_vs_control) %notin% starch.venn.list$recovery_vs_control
               ])
      
    # Query masigpro clusters to see if any of the genes from the DESeq2 contrasts were found significant and clustered by masigpro 
      starch.in.masigpro <- data.frame(sigs.masigpro[sigs.masigpro$Gene %in% 
                                                           unique(c(starch.venn.list$hypoxia_vs_control, 
                                                                    starch.venn.list$pcrit_vs_control, 
                                                                    starch.venn.list$anoxia_vs_control, 
                                                                    starch.venn.list$recovery_vs_control,
                                                                    starch.venn.list$Anoxia_vs_recovery,
                                                                    starch.venn.list$Hypoxia_vs_anoxia)), c(32,31)
      ])
    # Grab gene annotations from combined_sigs frame
      starch.in.masigpro <- merge(starch.in.masigpro, combined_sigs[,1:2], by = "Gene", all = FALSE)
      
      table(starch.in.masigpro$clusters) #Tally genes in each cluster
      
    #Generate pathway
      starch.pathway <- highlight_entities("tcf00500", c("131877115", "131882533", "131879841", "131890320", 
                                                         "131887366"), 
                                            fill_color = "paleturquoise", legend_name="Significant")
      starch.pathway
      
    #Export as png
      ggkeggsave(filename="starch pathway.png", starch.pathway, dpi=300)
      
    ### Generate heatmap of just the sig genes from the expression averaged wide dataframe
      starch.masigpro<- as.data.frame(expression.averaged.wide[expression.averaged.wide$Gene %in% combined_sigs[combined_sigs$Gene %in% starch.only$Gene, 1],])
      row.names(starch.masigpro) <- starch.masigpro$Gene
      starch.masigpro <- starch.masigpro[,-c(1)]
      
      
      jpeg(file = "starch Sigs Heatmap.jpg", width = 5, height = 10, units= "in", res = 300)
      pheatmap(starch.masigpro, cluster_rows = F, cluster_cols = F,
               show_rownames=T, border_color=NA,
               scale = "row", 
               labels_col = c("Normoxia", "Mild Hypox", "Pcrit", "Anoxia", "Recovery"),
               fontsize = 9,fontsize_row = 10, height=20, main = "starch Cycle")
      dev.off()
      
      
      
  #Fructose and mannose only 
    # Sum the number of antioxidant genes in our significant gene lists
      sum(combined_sigs$Gene %in% fructose.only$Gene)
      combined_sigs[combined_sigs$Gene %in% fructose.only$Gene, c(1,2)] #Which genes and description
      
    # How many are in the masigpro list and which clusters?
      sum(sigs.masigpro$Gene %in% fructose.only$Gene)  
      sigs.masigpro[sigs.masigpro$Gene %in% fructose.only$Gene, c(32,31)]
      
    # How many are in the hypoxia vs control DEG list?
      sum(Sigs.hypox.vs.cont$Gene %in% fructose.only$Gene)  
      Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% fructose.only$Gene, c(1,3,7)]
      
    # How many are in the pcrit vs control DEG list?
      sum(Sigs.pcrit.vs.cont$Gene %in% fructose.only$Gene)  
      Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% fructose.only$Gene, c(1,3,7)]
      
    # How many are in the anoxia vs control DEG list?
      sum(Sigs.anox.vs.control$Gene %in% fructose.only$Gene)  
      Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% fructose.only$Gene, c(1,3,7)]
      
    # How many are in the recovery vs control DEG list?
      sum(Sigs.recov.vs.cont$Gene %in% fructose.only$Gene)  
      Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% fructose.only$Gene, c(1,3,7)]
      
      
    # Make simple bar plot of the number of genes in response to oxygen compound GO term
    #Make data frame for plotting
      fructose.totals <- data.frame(times = c("Hypoxia vs control","Pcrit vs control","Anoxia vs control",
                                                "Recovery vs control"),
                                      sigs = c(0,0,5,5),
                                      total = rep(32,4)
      )
    # Reorder factor to plot how we want
      fructose.totals$times <- factor(fructose.totals$times, 
                                        levels = c("Recovery vs control","Anoxia vs control","Pcrit vs control",
                                                   "Hypoxia vs control"))
    # Create and export figure
      jpeg(filename = "Fructose genes DESeq2 totals.jpg", width = 6.5, height = 4, units = "in", res = 300)
      fructose.totals |>
        ggplot() +
        geom_segment(aes(x = 0, xend = sigs, y = times, yend = times), linewidth = 2)+
        geom_point(aes(x = sigs, y = times), size = 9, col = "darkblue")+
        geom_text(aes(x = sigs, y = times, label = sigs), color = "white") +
        scale_x_continuous(expand = expansion(mult = c(0.07, 0.3)), breaks = seq(0, 10, by = 2))+
        xlab("Number of genes") + ylab("") +
        theme(axis.text.x = element_text(size = 12),
              axis.title.y = element_blank(),
              title = element_text(size = 9),
              panel.grid.major.x = element_line(color = "lightgrey"))+
        ggtitle("Significant DEG's out of 32 fructose and mannose genes")
      dev.off()
      
    # How many significant genes are shared between the each time point's RTO GO term?
    # Create venn diagram first
    # Make list of genes in each group for venn diagram of time course programs
    # Deseq comparison
      fructose.venn.list <- list(
        hypoxia_vs_control = Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% fructose.only$Gene, 1],
        pcrit_vs_control = Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% fructose.only$Gene, 1],
        anoxia_vs_control = Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% fructose.only$Gene, 1],
        recovery_vs_control =  Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% fructose.only$Gene, 1],
        Anoxia_vs_recovery = Sigs.anox.vs.recov[Sigs.anox.vs.recov$Gene %in% fructose.only$Gene, 1],
        Hypoxia_vs_anoxia = Sigs.hypox.vs.anox[Sigs.hypox.vs.anox$Gene %in% fructose.only$Gene, 1]
      )
      
      
    # Draw venn diagram of all five groups using venn package  
      jpeg(filename = "Venn Diagram Fructose Genes.jpg", width = 8.5, height = 5.5, units = "in", res = 300)
      ggVennDiagram(fructose.venn.list[c(1,2,3,4)], label_percent_digit = 1, label_size = 4) +
        scale_x_continuous(expand = expansion(mult = .2)) +
        scale_fill_distiller(name = "Counts:", palette = "GnBu", direction = 1) + 
        theme(legend.position = "bottom", legend.direction = "horizontal", 
              plot.margin = unit(c(0,0,0,0), "cm"))
      dev.off()
      
    # Pull out genes unique to each group
    # NOTE: %notin% function is a custom function assigned at top of script
    # Genes unique to hypoxia vs control
      fructose.venn.list$hypoxia_vs_control[fructose.venn.list$hypoxia_vs_control %notin% 
                                                c(fructose.venn.list$pcrit_vs_control,
                                                  fructose.venn.list$anoxia_vs_control,
                                                  fructose.venn.list$recovery_vs_control,
                                                  fructose.venn.list$Anoxia_vs_recovery,
                                                  fructose.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to pcrit vs control
      fructose.venn.list$pcrit_vs_control[fructose.venn.list$pcrit_vs_control %notin% c(fructose.venn.list$hypoxia_vs_control,
                                                                                            fructose.venn.list$anoxia_vs_control,
                                                                                            fructose.venn.list$recovery_vs_control,
                                                                                            fructose.venn.list$Anoxia_vs_recovery,
                                                                                            fructose.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to anoxia vs control
      fructose.venn.list$anoxia_vs_control[fructose.venn.list$anoxia_vs_control %notin% c(fructose.venn.list$pcrit_vs_control,
                                                                                              fructose.venn.list$hypoxia_vs_control,
                                                                                              fructose.venn.list$recovery_vs_control,
                                                                                              fructose.venn.list$Anoxia_vs_recovery,
                                                                                              fructose.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to recovery vs control
      fructose.venn.list$recovery_vs_control[fructose.venn.list$recovery_vs_control %notin% c(fructose.venn.list$pcrit_vs_control,
                                                                                                  fructose.venn.list$anoxia_vs_control,
                                                                                                  fructose.venn.list$hypoxia_vs_control,
                                                                                                  fructose.venn.list$Anoxia_vs_recovery,
                                                                                                  fructose.venn.list$Hypoxia_vs_anoxia)
      ]
      
      
    # Pull out genes shared among all groups 
      fructose.venn.list$hypoxia_vs_control[fructose.venn.list$hypoxia_vs_control %in% fructose.venn.list$pcrit_vs_control &
                                                fructose.venn.list$hypoxia_vs_control %in% fructose.venn.list$anoxia_vs_control &
                                                fructose.venn.list$hypoxia_vs_control %in% fructose.venn.list$recovery_vs_control &
                                                fructose.venn.list$hypoxia_vs_control %in% fructose.venn.list$Anoxia_vs_recovery &
                                                fructose.venn.list$hypoxia_vs_control %in% fructose.venn.list$Hypoxia_vs_anoxia 
      ]
      
    # Pull out all genes that are in hypoxia, pcrit, and anoxia but NOT recovery
      unique(c(fructose.venn.list$hypoxia_vs_control, 
               fructose.venn.list$pcrit_vs_control, 
               fructose.venn.list$anoxia_vs_control)[c(fructose.venn.list$hypoxia_vs_control, 
                                                         fructose.venn.list$pcrit_vs_control, 
                                                         fructose.venn.list$anoxia_vs_control) %notin% fructose.venn.list$recovery_vs_control
               ])
      
    # Query masigpro clusters to see if any of the genes from the DESeq2 contrasts were found significant and clustered by masigpro 
      fructose.in.masigpro <- data.frame(sigs.masigpro[sigs.masigpro$Gene %in% 
                                                           unique(c(fructose.venn.list$hypoxia_vs_control, 
                                                                    fructose.venn.list$pcrit_vs_control, 
                                                                    fructose.venn.list$anoxia_vs_control, 
                                                                    fructose.venn.list$recovery_vs_control,
                                                                    fructose.venn.list$Anoxia_vs_recovery,
                                                                    fructose.venn.list$Hypoxia_vs_anoxia)), c(32,31)
      ])
    # Grab gene annotations from combined_sigs frame
      fructose.in.masigpro <- merge(fructose.in.masigpro, combined_sigs[,1:2], by = "Gene", all = FALSE)
      
      table(fructose.in.masigpro$clusters) #Tally genes in each cluster
      
    # Grab gene annotations from combined_sigs frame
      starch.in.masigpro <- merge(starch.in.masigpro, combined_sigs[,1:2], by = "Gene", all = FALSE)
      
      table(starch.in.masigpro$clusters) #Tally genes in each cluster
      
    #Generate pathway
      fructose.pathway <- highlight_entities("tcf00051", c("131879671", "131889496", "131884437", "131879841", "131877224"), 
                                           fill_color = "paleturquoise", legend_name="Significant")
      fructose.pathway
      
    #Export as png
      ggkeggsave(filename="fructose pathway.png", fructose.pathway, dpi=900)
      
    ### Generate heatmap of just the sig genes from the expression averaged wide dataframe
      fructose.masigpro<- as.data.frame(expression.averaged.wide[expression.averaged.wide$Gene %in% combined_sigs[combined_sigs$Gene %in% fructose.only$Gene, 1],])
      row.names(fructose.masigpro) <- fructose.masigpro$Gene
      fructose.masigpro <- fructose.masigpro[,-c(1)]
      
      
      jpeg(file = "fructose Sigs Heatmap.jpg", width = 5, height = 7, units= "in", res = 300)
      pheatmap(fructose.masigpro, cluster_rows = F, cluster_cols = F,
               show_rownames=T, border_color=NA,
               scale = "row", 
               labels_col = c("Normoxia", "Mild Hypox", "Pcrit", "Anoxia", "Recovery"),
               fontsize = 9,fontsize_row = 10, height=20, main = "fructose Cycle")
      dev.off()
      
      
      
  # Mitochondria targeted genes
      
    # Sum the number of antioxidant genes in our significant gene lists
      sum(combined_sigs$Gene %in% mitos$Gene)
      combined_sigs[combined_sigs$Gene %in% mitos$Gene, c(1,2)] #Which genes and description
      
    # How many are in the masigpro list and which clusters?
      sum(sigs.masigpro$Gene %in% mitos$Gene)  
      sigs.masigpro[sigs.masigpro$Gene %in% mitos$Gene, c(32,31)]
      
    # How many are in the hypoxia vs control DEG list?
      sum(Sigs.hypox.vs.cont$Gene %in% mitos$Gene)  
      Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% mitos$Gene, c(1,3,7)]
      
    # How many are in the pcrit vs control DEG list?
      sum(Sigs.pcrit.vs.cont$Gene %in% mitos$Gene)  
      Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% mitos$Gene, c(1,3,7)]
      
    # How many are in the anoxia vs control DEG list?
      sum(Sigs.anox.vs.control$Gene %in% mitos$Gene)  
      Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% mitos$Gene, c(1,3,7)]
      
    # How many are in the recovery vs control DEG list?
      sum(Sigs.recov.vs.cont$Gene %in% mitos$Gene)  
      Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% mitos$Gene, c(1,3,7)]
      
      
    # Make simple bar plot of the number of genes in response to oxygen compound GO term
    #Make data frame for plotting
      mitos.totals <- data.frame(times = c("Hypoxia vs control","Pcrit vs control","Anoxia vs control",
                                            "Recovery vs control"),
                                  sigs = c(6,18,32,24),
                                  total = rep(28,4)
      )
    # Reorder factor to plot how we want
      mitos.totals$times <- factor(mitos.totals$times, 
                                    levels = c("Recovery vs control","Anoxia vs control","Pcrit vs control",
                                               "Hypoxia vs control"))
    # Create and export figure
      jpeg(filename = "Mito targeted genes DESeq2 totals.jpg", width = 6.5, height = 4, units = "in", res = 300)
      mitos.totals |>
        ggplot() +
        geom_segment(aes(x = 0, xend = sigs, y = times, yend = times), linewidth = 2)+
        geom_point(aes(x = sigs, y = times), size = 9, col = "darkblue")+
        geom_text(aes(x = sigs, y = times, label = sigs), color = "white") +
        scale_x_continuous(expand = expansion(mult = c(0, 0.2)), breaks = seq(0, 35, by = 5))+
        xlab("Number of genes") + ylab("") +
        theme(axis.text.x = element_text(size = 12),
              axis.title.y = element_blank(),
              title = element_text(size = 9),
              panel.grid.major.x = element_line(color = "lightgrey"))+
        ggtitle("Significant DEG's out of 600 mitochondria-targeted genes")
      dev.off()
      
    # How many significant genes are shared between the each time point's RTO GO term?
    # Create venn diagram first
    # Make list of genes in each group for venn diagram of time course programs
    # Deseq comparison
      mitos.venn.list <- list(
        hypoxia_vs_control = Sigs.hypox.vs.cont[Sigs.hypox.vs.cont$Gene %in% mitos.only$Gene, 1],
        pcrit_vs_control = Sigs.pcrit.vs.cont[Sigs.pcrit.vs.cont$Gene %in% mitos.only$Gene, 1],
        anoxia_vs_control = Sigs.anox.vs.control[Sigs.anox.vs.control$Gene %in% mitos.only$Gene, 1],
        recovery_vs_control =  Sigs.recov.vs.cont[Sigs.recov.vs.cont$Gene %in% mitos.only$Gene, 1],
        Anoxia_vs_recovery = Sigs.anox.vs.recov[Sigs.anox.vs.recov$Gene %in% mitos.only$Gene, 1],
        Hypoxia_vs_anoxia = Sigs.hypox.vs.anox[Sigs.hypox.vs.anox$Gene %in% mitos.only$Gene, 1]
      )
      
      
    # Draw venn diagram of all five groups using venn package  
      jpeg(filename = "Venn Diagram mito targeted Genes.jpg", width = 8.5, height = 5.5, units = "in", res = 300)
      ggVennDiagram(mitos.venn.list[c(1,2,3,4)], label_percent_digit = 1, label_size = 4) +
        scale_x_continuous(expand = expansion(mult = .2)) +
        scale_fill_distiller(name = "Counts:", palette = "GnBu", direction = 1) + 
        theme(legend.position = "bottom", legend.direction = "horizontal", 
              plot.margin = unit(c(0,0,0,0), "cm"))
      dev.off()
      
    # Pull out genes unique to each group
    # NOTE: %notin% function is a custom function assigned at top of script
    # Genes unique to hypoxia vs control
      mitos.venn.list$hypoxia_vs_control[mitos.venn.list$hypoxia_vs_control %notin% 
                                                c(mitos.venn.list$pcrit_vs_control,
                                                  mitos.venn.list$anoxia_vs_control,
                                                  mitos.venn.list$recovery_vs_control,
                                                  mitos.venn.list$Anoxia_vs_recovery,
                                                  mitos.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to pcrit vs control
      mitos.venn.list$pcrit_vs_control[mitos.venn.list$pcrit_vs_control %notin% c(mitos.venn.list$hypoxia_vs_control,
                                                                                            mitos.venn.list$anoxia_vs_control,
                                                                                            mitos.venn.list$recovery_vs_control,
                                                                                            mitos.venn.list$Anoxia_vs_recovery,
                                                                                            mitos.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to anoxia vs control
      mitos.venn.list$anoxia_vs_control[mitos.venn.list$anoxia_vs_control %notin% c(mitos.venn.list$pcrit_vs_control,
                                                                                              mitos.venn.list$hypoxia_vs_control,
                                                                                              mitos.venn.list$recovery_vs_control,
                                                                                              mitos.venn.list$Anoxia_vs_recovery,
                                                                                              mitos.venn.list$Hypoxia_vs_anoxia)
      ]
      
    # Genes unique to recovery vs control
      mitos.venn.list$recovery_vs_control[mitos.venn.list$recovery_vs_control %notin% c(mitos.venn.list$pcrit_vs_control,
                                                                                                  mitos.venn.list$anoxia_vs_control,
                                                                                                  mitos.venn.list$hypoxia_vs_control,
                                                                                                  mitos.venn.list$Anoxia_vs_recovery,
                                                                                                  mitos.venn.list$Hypoxia_vs_anoxia)
      ]
      
      
    # Pull out genes shared among all groups 
      mitos.venn.list$hypoxia_vs_control[mitos.venn.list$hypoxia_vs_control %in% mitos.venn.list$pcrit_vs_control &
                                                mitos.venn.list$hypoxia_vs_control %in% mitos.venn.list$anoxia_vs_control &
                                                mitos.venn.list$hypoxia_vs_control %in% mitos.venn.list$recovery_vs_control &
                                                mitos.venn.list$hypoxia_vs_control %in% mitos.venn.list$Anoxia_vs_recovery &
                                                mitos.venn.list$hypoxia_vs_control %in% mitos.venn.list$Hypoxia_vs_anoxia 
      ]
      
    # Pull out all genes that are in hypoxia, pcrit, and anoxia but NOT recovery
      unique(c(mitos.venn.list$hypoxia_vs_control, 
               mitos.venn.list$pcrit_vs_control, 
               mitos.venn.list$anoxia_vs_control)[c(mitos.venn.list$hypoxia_vs_control, 
                                                         mitos.venn.list$pcrit_vs_control, 
                                                         mitos.venn.list$anoxia_vs_control) %notin% mitos.venn.list$recovery_vs_control
               ])
      
    # Query masigpro clusters to see if any of the genes from the DESeq2 contrasts were found significant and clustered by masigpro 
      mitos.in.masigpro <- data.frame(sigs.masigpro[sigs.masigpro$Gene %in% 
                                                           unique(c(mitos.venn.list$hypoxia_vs_control, 
                                                                    mitos.venn.list$pcrit_vs_control, 
                                                                    mitos.venn.list$anoxia_vs_control, 
                                                                    mitos.venn.list$recovery_vs_control,
                                                                    mitos.venn.list$Anoxia_vs_recovery,
                                                                    mitos.venn.list$Hypoxia_vs_anoxia)), c(32,31)
      ])
    # Grab gene annotations from combined_sigs frame
      mitos.in.masigpro <- merge(mitos.in.masigpro, combined_sigs[,1:2], by = "Gene", all = FALSE)
      
      table(mitos.in.masigpro$clusters) #Tally genes in each cluster
      