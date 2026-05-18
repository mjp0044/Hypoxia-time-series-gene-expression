# Hypoxia-time-series-gene-expression

## Interactive results at: https://mjp0044.github.io/Hypoxia-time-series-gene-expression/



# Analysis highlights 
 
This experiment consisted of measuring gene expression over the course of a typical hypoxia exposure event in the marine crustacean *Tigriopus californicus*. 

Below is an illustration showing the longitudinal sampling over time. There were 5 points at which we sampled RNA. 

- Normoxia (100% oxygen)
- Mild hypoxia (1/2 way through oxygen depletion)
- *P*<sub>crit</sub> (the point at which this species can no longer control its respiratory rate)
- 1 hour of anoxia (0% oxygen)
- 2 hours into recovery (back at 100% oxygen)

<img src="Figures and Tables/Fig 1.png" width="500">

*Fig 1: Experimental setup of the hypoxia course.*


---------------------------------------------------------------------------------------------------------------------------

## 1. Identifying genes with longitudinal differential expression

The first thing I'd like to highlight is the statistical analysis, because I think the package and approach is pretty cool. 

The heaviest load of the analysis in R is carried by the gene expression packages `DESeq2` and `maSigPro`. 

Following batch correction with `Combat-seq` and filtering genes how with low counts, I used `DESeq2` to normalize my count matrix (this was in the form of genes as rownames and mRNA counts for each sample in the columns. 

```r
#Create DESeq object
    dds <- DESeqDataSetFromMatrix(countData = ft_filtered_combat,
                                  colData = coldata,
                                  design = ~ Time)
    dds #View object
    
    dds2 <-DESeq(dds) #run the normalization and save as new object
    
    rm(dds) #Get rid of old object
    
    ddsNorm <- counts(dds2, normalized=TRUE) #Get the normalized counts out for maSigPro.

> head(ddsNorm[,c(1:10)])
                   C_1         C_2         C_3         C_4          C_5         C_6        35_1        35_2        35_3        35_4
TCAL_00001 50985.71117 47423.38269 38153.12984 47922.59540 58051.463491 60699.60141 48459.15642 31447.46699 44186.09381 49988.85621
TCAL_00002   260.49746   220.62677   236.55538   251.41254   352.902264   346.15194   251.60103   301.42577   244.76859   265.58548
TCAL_00003    21.53719    18.19602    13.19082    14.01929    19.179471    20.97891    21.96517    15.45773    16.57795     8.39796
TCAL_00004   229.73004   251.33256   246.22865   234.58940   267.553619   203.11395   239.62003   198.37422   258.42102   202.60078
TCAL_00006    16.40929    20.47052    25.50225    26.16934     8.630762    23.83966    19.96834    21.89845    26.32969    36.74107
TCAL_00008   303.57184   325.25390   327.13235   346.74369   353.861238   234.58230   354.43795   365.83298   359.83908   271.88395
```


Most people use `DESeq2` to identify differentially expressed genes between two groups (e.g., hot vs cold treatments or between two time points). 

However, my goal was to find genes with significantly changing trajectories over multiple time points. `DESeq2` can do this with consecutive pairwise contrasts, but the package `maSigPro` is built specifically for time series analyses. 

`maSigPro` isn't unique in this regard, but it does have a lot of cool features that made it ideal for my needs. For example, it can identify genes that respond differently between treatments (i.e., have different expression patterns over time) but can also handle single series without groups. 

My data takes the form of the latter. It's just 6 replicates of gene expression data over the hypoxia exposure course. 



![](https://github.com/mjp0044/Hypoxia-time-series-gene-expression/blob/85d4ad4b4d69f3bccf42f32cd593a32d4166318a/Figures/Fig%202%20maSigPro%20cluster%20patterns%209%20clusters.jpg)
*Fig 2: The expression patterns of genes identified to significantly respond over the hypoxia course by maSigPro. Genes were clustered into groups using the default “hclust” (hierarchical cluster analysis) method using Ward.D aggregation and the default k = 9 clusters by maSigPro. Black dots are the mean expression at each time point and error bars represent the standard error of the mean.*


---------------------------------------------------------------------------------------------------------------------------


![](https://github.com/mjp0044/Hypoxia-time-series-gene-expression/blob/11d05ee6aabb724b562a904c7acbb84016ecf49a/Figures/Fig%203%20Pcrit%20only%20genes%20time%20series%20line%20plots.jpg)
*Fig 3: Genes not marked as significant by maSigPro but identified as significantly different in expression between Pcrit and Normoxia levels by DESeq2. Gene-standardized z-scores plotted on the y-axis to show relative change per gene for plotting. Time points are N = normoxia, MH = mild hypoxia, P = Pcrit, A = anoxia, and R = recovery. 


---------------------------------------------------------------------------------------------------------------------------


