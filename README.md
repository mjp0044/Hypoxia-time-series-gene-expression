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

Following batch correction with `Combat-seq` and filtering genes how with low counts, I used `DESeq2` to normalize my count matrix.

**ft_filtered_combat** = count matrix
**coldata** = metadata file for DESeq2 to understand treatments and replicates

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

Above you can see a snippet of the matrix with the gene IDs on the row names and the sample replicates on the column names (For context, C = normoxia and 35 = mild hypoxia; numbers indicate replicates)

Most people use `DESeq2` to identify differentially expressed genes between two groups (e.g., hot vs cold treatments or between two time points). 

However, my goal was to find genes with significantly changing trajectories over multiple time points. `DESeq2` can do this with consecutive pairwise contrasts, but the package `maSigPro` is built specifically for time series analyses. 

`maSigPro` isn't unique in this regard, but it does have a lot of cool features that made it ideal for my needs. For example, it can identify genes that respond differently between treatments (i.e., have different expression patterns over time) but can also handle single series without groups. 

My data takes the form of the latter. There was no separate treatments, all the organisms got the same hypoxia exposure course. We just wanted to know how this species would respond in general to hypoxia. 

The key thing with using `maSigPro` (and most differential expression programs) is you need to provide an accurate meta data file for the program to understand your data structure. 

Now that you've seen a snippet of the count matrix, you can see the meta data file for maSigPro: 

```r
> head(hypoxia.edesign, 12)
     Time Replicate Group
C_1     1         1     1
C_2     1         2     1
C_3     1         3     1
C_4     1         4     1
C_5     1         5     1
C_6     1         6     1
35_1    2         1     1
35_2    2         2     1
35_3    2         3     1
35_4    2         4     1
35_5    2         5     1
35_6    2         6     1
```

You can see that samples are in the rows and match the columns of the count matrix, time points are indicates by the Time column, replicate structure in the Replicate column, and finally, the treatments (or lack thereof) indicated in the Group column. For the Group column, all the values just say "1" to indicate we have no separate treatments. 


First thing to do when using `maSigPro` is to convert your meta data file into a design object. 

```r
#Sanity check: make sure rownames of the info key match the colnames of the count matrix
    all(rownames(hypoxia.edesign) == colnames(ddsNorm)) #needs to be TRUE
    
#Create a regression matrix for the full regression model for the single series:
    hypoxia.design <- make.design.matrix(hypoxia.edesign, degree = 3)
```

Next, we can fit our model. `maSigPro` generates a regression matrix to model gene expression over time. This model can be modified to include polynomial terms that let you capture changes in expression that are not just linear, but also changes that look quadratic or even cubic. 

<img src="Figures and Tables/masigpro method.png" width="600">


Because we had 5 time points, I included both a 2<sup>nd</sup> order (quadratic) and 3<sup>rd</sup> order (cubic) polynomial term. That was specified above when making the design matrix. The design matrix looks like this with those terms included (denoted Time2 and Time3):

```r
> head(hypoxia.design)
$dis
     Time Time2 Time3
C_1     1     1     1
C_2     1     1     1
C_3     1     1     1
C_4     1     1     1
C_5     1     1     1
C_6     1     1     1
35_1    2     4     8
35_2    2     4     8
35_3    2     4     8
35_4    2     4     8
35_5    2     4     8
35_6    2     4     8
```

With our model design structure specified, now I can fit the model using the `p.vector` function. 

```r
# By default maSigPro corrects this p-value for multiple comparisons by applying the linear step-up (B-H) false discovery rate (FDR) 
# The level of FDR control is given by the function parameter Q.
# counts = FALSE because I'm using normalized count data
  ss.hypoxia <- p.vector(ddsNorm, hypoxia.design, Q = 0.05, MT.adjust = "BH", counts = FALSE)
```

Let's see how many significant genes we found! 

```r
ss.hypoxia$i # returns the number of significant genes
[1] 1347
```

# Next we can view how the epxression of those 1,347 genes changed over time using a cluster analysis. 

I used the built-in clustering algorithm included in the `maSigPro` package. This algorithm lets me set a hierarchical clustering method and specify any number of centroids. 

I tested values of k = 6 through k = 11 and plotted the gap statistic, silhouette score, and elbow (wss). To do this, I generated a version of my data that is scaled so that I could work with z-scores and a kmeans approach for calculating these statistics. If we don't use scaled data, then our clustering doesn't match the type of pattern recognition that the hierarchical approach uses. 

```r
# Prepare data
    # expr_mat is your replicate-level matrix
     expr_mat <- hypoxia.sigs$sig.genes$Group$sig.profiles  
    
    # Make a condition label for each column (strip "_1", "_2", etc.)
     sample_groups <- gsub("_[0-9]+$", "", colnames(expr_mat))
    
    # Collapse replicates by taking row means
      expr_avg <- sapply(unique(sample_groups), function(grp) {
        rowMeans(expr_mat[, sample_groups == grp, drop = FALSE])
      })
    
    # Scale for clustering
      expr_scaled <- scale(expr_avg)
```

Now that we have scaled expression data, we can compute wss, silhouette, and gap statistics. 

```r
  ks <- 6:11 # Cluster numbers to test
    
  # WSS (Elbow method)
    wss <- sapply(ks, function(k) {
      kmeans(expr_scaled, centers = k, nstart = 25)$tot.withinss
    })
    
  # Heuristic: choose k where % reduction in WSS drops the most ("elbow")
  # Compute relative reduction in WSS
    rel_diff <- diff(wss) / wss[-length(wss)]
    best_k_wss <- ks[which.min(rel_diff) + 1]  # +1 because diff reduces length by 1
    
    df_wss <- data.frame(k = ks, WSS = wss)
    
  # Silhouette
    avg_sil <- sapply(ks, function(k) {
      km <- kmeans(expr_scaled, centers = k, nstart = 25)
      sil <- silhouette(km$cluster, dist(expr_scaled))
      mean(sil[, 3])
    })
    best_k_sil <- ks[which.max(avg_sil)]
    df_sil <- data.frame(k = ks, avg_sil = avg_sil)
    
  #  Gap statistic
    set.seed(1974) # Go sounders!
    gap_stat <- clusGap(expr_scaled, FUN = kmeans, K.max = max(ks),
                        B = 50, nstart = 25)
    gap_df <- data.frame(
      k = ks,
      gap = gap_stat$Tab[ks, "gap"],
      SE.sim = gap_stat$Tab[ks, "SE.sim"]
    )
    best_k_gap <- ks[which.max(gap_df$gap)]
```

And then we can generate our plots: 

```r
 # Plotting
  # WSS
    wss.plot <- ggplot(df_wss, aes(x = k, y = WSS)) +
      geom_line() + geom_point(size = 3) +
      geom_vline(xintercept = best_k_wss, linetype = "dashed", color = "red") +
      scale_x_continuous(breaks = ks) +
      scale_y_continuous(limits = c(100, 500), breaks = seq(100, 500, by =100))+
      ylab("Within-cluster SS") + xlab("Number of clusters (k)") +
      ggtitle(paste0("Elbow plot (best k = ", best_k_wss, ")"))
    wss.plot
    
  # Silhouette
    silhouette.plot <- ggplot(df_sil, aes(x = k, y = avg_sil)) +
      geom_line() + geom_point(size = 3) +
      geom_vline(xintercept = best_k_sil, linetype = "dashed", color = "red") +
      scale_x_continuous(breaks = ks) +
      scale_y_continuous(limits = c(0.55, 0.8), breaks = seq(0.55, 0.8, by =0.05))+
      ylab("Average silhouette width") + xlab("Number of clusters (k)") +
      ggtitle(paste0("Silhouette method (best k = ", best_k_sil, ")"))
    silhouette.plot
    
  # Gap statistic
    gap.plot <- ggplot(gap_df, aes(x = k, y = gap)) +
      geom_line() + geom_point(size = 3) +
      geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), width = 0.2) +
      geom_vline(xintercept = best_k_gap, linetype = "dashed", color = "red") +
      scale_x_continuous(breaks = ks) +
      scale_y_continuous(limits = c(2.8, 3.3), breaks = seq(2.8, 3.3, by =0.1))+
      ylab("Gap statistic") + xlab("Number of clusters (k)") +
      ggtitle(paste0("Gap statistic (best k = ", best_k_gap, ")"))
    gap.plot
    
  # Combine plots vertically
    combined_plot <- wss.plot / silhouette.plot / gap.plot
    combined_plot
```
<img src="Figures and Tables/Clustering statistics.jpg" width="400">

As you can see above, our statistics don't agree on a single best value for k. They often don't because each can be biased toward a lower number of clusters for parsimony or a greater number for increased structure. 

But that's ok. We can make decisions based on the interprebility of the patterns we see in the clusters. By that I mean we can visibly choose a value for k within a reasoable range based on the statistics and also based on whether adding or subtracting values of k creates redundant patterns or lumps interesting clusters that we don't want to lose into a larger group. 

In the end, I chose k = 9 because it produced enough interesting clusters on which I could further explore based on their patterns of expression over time and it has a reasonable value for wss (>400), silhouette (>0.6), and gap (3rd highest value). 


![](https://github.com/mjp0044/Hypoxia-time-series-gene-expression/blob/85d4ad4b4d69f3bccf42f32cd593a32d4166318a/Figures/Fig%202%20maSigPro%20cluster%20patterns%209%20clusters.jpg)
*Fig 2: The expression patterns of genes identified to significantly respond over the hypoxia course by maSigPro. Genes were clustered into groups using the default “hclust” (hierarchical cluster analysis) method using Ward.D aggregation and the default k = 9 clusters by maSigPro. Black dots are the mean expression at each time point and error bars represent the standard error of the mean.*


---------------------------------------------------------------------------------------------------------------------------


![](https://github.com/mjp0044/Hypoxia-time-series-gene-expression/blob/11d05ee6aabb724b562a904c7acbb84016ecf49a/Figures/Fig%203%20Pcrit%20only%20genes%20time%20series%20line%20plots.jpg)
*Fig 3: Genes not marked as significant by maSigPro but identified as significantly different in expression between Pcrit and Normoxia levels by DESeq2. Gene-standardized z-scores plotted on the y-axis to show relative change per gene for plotting. Time points are N = normoxia, MH = mild hypoxia, P = Pcrit, A = anoxia, and R = recovery. 


---------------------------------------------------------------------------------------------------------------------------


