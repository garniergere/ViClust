## ViClust
**Vi**sualizing and **Clust**ering Sequenom data in batch 

(R scripts for SNP data from Sequenom: standalone R and Galaxy version using Python wrapper, with example files)

## What is ViClust doing?

This tool imports the raw xml file from the Typer-4_0_20 Sequenom equipment software and allows: 

**1)** to visualize the original clusters proposed by the Sequenom method based on a mixture of Gaussian distributions (Johansen *et al.* 2013) -> (**Sequenom clusters**)

**2)** to visualize the same clusters after proposed additional and more of less stringent filters, based on the different allele signal magnitudes -> (**Filter 1 clusters**),

**3)** to perform an alternative hierarchical clustering method (Hclust) using the Ward algorithm (Ward 1963) which relaxes the constraining assumption of normal distribution for allele signals, then to visualize the new clusters -> (**alternative Hclust clusters**), and to export the proposed genotypes assignment file.

The xml file is loaded then parsed to extract raw data and original clustering information, and angle and magnitude values are recomputed based on raw peak height data from the different allele signals. The program gives SNPs plots across the 3 different alternative options in batch, allows comparing and validating each SNPs more easily, and exports automatically the proposed genotype assignments, making the whole process much less time consuming.

**Sequenom clusters** are plotted with x-axis data representing the converted angle of the signal and y-axis data the magnitude of the signal, both derived from the xml file.

**Filter 1 clusters** are similar to those obtained from the Sequenom Software except that datapoints have been excluded after a filter for magnitude chosen by the user (default value is 6, can be changed, filtering stage is skipped if the value is set to "0")
Warning: At this stage, if the user has originally excluded datapoints (i.e. before exporting the xml datafile), those should be identified globally with grey symbols, but we recommend using the raw xml data file to avoid errors.

**Alternative Hclust clusters** are those obtained after the hierarchical clustering method. 
The Hclust method will be applied if the hclust.condition is "Yes"in the list of argument (see examples). If the Hclust method is chosen and filter 1 has been previously skipped (i.e. set to 0), the filter default value of “6” is nevertheless applied on magnitude, since it is required before performing the hierarchical clustering. 
Warning: At this stage, if the user has originally excluded datapoints, those are identified with pale color symbols (grey, pale green, pale pink) corresponding to the genotypes assignment to different clusters. However, we recommend using the raw xml data file to allow an optimal application of the alternative clustering method.

Genotype data used for each clustering method and plots are provided by the program:

**Rawdat.txt** contains the original data in a user-friendly format (see the User's guide from the Typer-4_0_20 Sequenom software for variable description)
**res.filter1.txt** contains the data which passed the filter on magnitude, the Sequenom original call (ori.genotyteId), the transformed data (angle and magnitude) used for plotting clusters with a unique identifier for individuals (sampleId.new) to deal with replicates, and an index variable indicating whether data points passed the filter on magnitude (call.filter1).
**res.with.hclust.txt** also contains the transformed data with the proposed new call in up to 3 groups (HClust.groups variable).
These files can easily be further imported and processed in any software (such as R).
The different methods can be applied on either all uploaded data or a smaller group of SNPs and/or individuals (which are loaded as text files by the user, see argument list and examples).

What could change with the results from the Hclust method is the way in which individuals with intermediate angle values between clusters are being attributed to a particular cluster (and thus a genotype call).
The difference between the Hclust and the Sequenom methods is that Hclust is not assuming any underlying distribution for the data, favoring a criterion of minimum variance within groups. Based on our experience, this may improve objectively the clustering in a number of cases.

Both methods are however subject to caution, especially if you observe sub-structuring in any of the optimal angle value ranges for clustering, and/or intermediate values outside the optimal ranges  (see the Typer-4_0_20 User Guide for more information on optimal range values).

**Examples**
See examples of plots and result files provided in the Example.files folder where the results from the Hclust method can be compared to the Sequenom proposed clustering. 
See also the ViClust.2.Matrix.R script for reformatting ViClust output into a classical genotype format.
See Rscript.and.folder.description.txt for further details.

**Help tips:** If anything goes wrong, please 1) check the format/content of the loaded xml file, 2) check the format/content of the SNP list or of the individuals' list if any is used (these should be simple text files, one line per SNP or per individual).

**Literature cited** 
Ward JH. 1963. Hierarchical Grouping to Optimize an Objective Function. Journal of the American Statistical Association. 58:236-244. </br>
Johansen PP, Andersen JD, Børsting C, Morling N (2013) Evaluation of the iPLEX Sample ID Plus Panel designed for the Sequenom MassARRAY system. A SNP typing assay developed for human identification and sample tracking based on the SNP for ID panel. Forensic Science International: Genetics 7: 482–487. https://doi.org/10.1016/j.fsigen.2013.04.009


## How to cite the ViClust tool?

If you are using ViClust, please cite:

Bouteiller XP, Barraquand F, Garnier-Géré P, Harmand N, Laizet Y, Raimbault A, Segura R, Lassois L, Monty A, Verdu C, Mariette S & Porté AJ (2018) No evidence for genetic differentiation in juvenile traits between Belgian and French populations of the invasive tree <i>Robinia pseudoacacia </i>. Plant Ecology and Evolution 151 (1): 5–17, https://doi.org/10.5091/plecevo.2018.1403

## Contact 
Pauline Garnier-Géré ( pauline.garniergere@inra.fr )
