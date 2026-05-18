# Hypoxia-time-series-gene-expression

## Interactive results at: https://mjp0044.github.io/Hypoxia-time-series-gene-expression/



# Analysis highlights 

## 1. Identifying genes with longitudinal differential expression
 
This experiment consisted of measuring gene expression over the course of a typical hypoxia exposure event in the marine crustacean *Tigriopus californicus*. 

Below is an illustration showing the longitudinal sampling over time. There were 5 points at which we sampled RNA. 

- Normoxia (100% oxygen)
- Mild hypoxia (1/2 way through oxygen depletion)
- *P*<sub>crit</sub> (the point at which this species can no longer control its respiratory rate)
- 1 hour of anoxia (0% oxygen)
- 2 hours into recovery (back at 100% oxygen)

<img src="figures/DO_timeseries.png" width="600">

*Fig 1: Experimental setup of the hypoxia course.*


---------------------------------------------------------------------------------------------------------------------------



![](https://github.com/mjp0044/Hypoxia-time-series-gene-expression/blob/85d4ad4b4d69f3bccf42f32cd593a32d4166318a/Figures/Fig%202%20maSigPro%20cluster%20patterns%209%20clusters.jpg)
*Fig 2: The expression patterns of genes identified to significantly respond over the hypoxia course by maSigPro. Genes were clustered into groups using the default “hclust” (hierarchical cluster analysis) method using Ward.D aggregation and the default k = 9 clusters by maSigPro. Black dots are the mean expression at each time point and error bars represent the standard error of the mean.*


---------------------------------------------------------------------------------------------------------------------------


![](https://github.com/mjp0044/Hypoxia-time-series-gene-expression/blob/11d05ee6aabb724b562a904c7acbb84016ecf49a/Figures/Fig%203%20Pcrit%20only%20genes%20time%20series%20line%20plots.jpg)
*Fig 3: Genes not marked as significant by maSigPro but identified as significantly different in expression between Pcrit and Normoxia levels by DESeq2. Gene-standardized z-scores plotted on the y-axis to show relative change per gene for plotting. Time points are N = normoxia, MH = mild hypoxia, P = Pcrit, A = anoxia, and R = recovery. 


---------------------------------------------------------------------------------------------------------------------------


