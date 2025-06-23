# Hypoxia-time-series-gene-expression

## Gene expression patterns normalized and plotted as Z-scores across the five time points found below. 
Interactive results at: https://mjp0044.github.io/Hypoxia-time-series-gene-expression/




## Results of the cluster analysis following identification of significant genes using maSigPro
 

![](https://github.com/mjp0044/Hypoxia-time-series-gene-expression/blob/8d8f83fc64ae410017d69020e6ca1bc90e9d8c5c/Figures/Fig%201%20Experiment%20setup.png)
*Fig 1: Experimental setup of the hypoxia course.*


---------------------------------------------------------------------------------------------------------------------------



![](https://github.com/mjp0044/Hypoxia-time-series-gene-expression/blob/85d4ad4b4d69f3bccf42f32cd593a32d4166318a/Figures/Fig%202%20maSigPro%20cluster%20patterns%209%20clusters.jpg)
*Fig 2: The expression patterns of genes identified to significantly respond over the hypoxia course by maSigPro. Genes were clustered into groups using the default “hclust” (hierarchical cluster analysis) method using Ward.D aggregation and the default k = 9 clusters by maSigPro. Black dots are the mean expression at each time point and error bars represent the standard error of the mean.*


---------------------------------------------------------------------------------------------------------------------------


![](https://github.com/mjp0044/Hypoxia-time-series-gene-expression/blob/11d05ee6aabb724b562a904c7acbb84016ecf49a/Figures/Fig%203%20Pcrit%20only%20genes%20time%20series%20line%20plots.jpg)
*Fig 3: Genes not marked as significant by maSigPro but identified as significantly different in expression between Pcrit and Normoxia levels by DESeq2. Gene-standardized z-scores plotted on the y-axis to show relative change per gene for plotting. Time points are N = normoxia, MH = mild hypoxia, P = Pcrit, A = anoxia, and R = recovery. 


---------------------------------------------------------------------------------------------------------------------------


![](https://github.com/mjp0044/Hypoxia-time-series-gene-expression/blob/bc70dd11516628d7f612e2149e5941d8f3cbfa93/Figures/Fig%204%20Pathway%20summary%20figure.png)
*Fig 4: A simplified diagram uniting the starch-sucrose pathway, glycolysis, chitin synthesis pathway, and chitin metabolism and showing genes that significantly responded to the hypoxia course. This diagram is a reduced version of Figure S11, which shows the entirety of the KEGG pathways for glycolysis, starch-sucrose metabolism, pentose-phosphate metabolism, fructose-mannose metabolism, TCA cycle, and amino sugar metabolism. The heatmaps for significant genes are colored based on gene expression with red indicating an increase and blue a decrease. The five squares correspond to the five hypoxia levels over the exposure course: normoxia, mild hypoxia, Pcrit, anoxia, and recovery. Inset on each heat map indicates the cluster or pairwise comparison in which that gene was identified.
