This workshop is a follow-up of the Microbiome analysis using QIIME2 workshop. The result from the previous workshop will be used to demonstrate basic analyses of microbiota data to determine if and how communities differ by variables of interest using R. This all-day workshop will consist of lectures and hands on training to analyze from raw dataset through publication-quality statistics and visualizations


## Get Started
### Download and Install
* Base R: http://cran.mtu.edu/

NOTE: If you need to update to the most recent version of R on Windows you can do so using the installr package. Instructions here. For OSX and Ubuntu, download from CRAN using the link above.

* RStudio: https://www.rstudio.com/products/rstudio/download3/


## Project folder: Stay organized
All of our analyses will be organized into a “Project”.

Make a new project by selecting File->New project. Select “New Directory” and “Empty Project”. Name the project “Microbiome_Analysis” and save the project to your Desktop. Place all of your files for this analysis in the folder created on the Desktop

Create a new R script (File->New file->R script) to save your code. This file will automatically be saved in the project folder.

Now your screen should look like this

* Upper left: Where you type and save the code you want to run.
* Upper right: Files you load into and create in R. To view one, click on it and it will open in the upper left pane.
* Lower left: The console. Where commands and outputs run (similar to the one mothur window).
* Lower right: Variable. Explore the different tabs.

## Data description

In this workshop, We will work with data processed through Qiime2 and can be downloaded from the following links
* [Metadata](https://data.qiime2.org/2018.4/tutorials/moving-pictures/sample_metadata.tsv) 
* [Dada2 feature table](https://docs.qiime2.org/2018.4/data/tutorials/moving-pictures/table.qza)
* [taxonomy assignments](https://docs.qiime2.org/2018.4/data/tutorials/moving-pictures/taxonomy.qza)
* [midpointed-tree](https://docs.qiime2.org/2018.4/data/tutorials/moving-pictures/rooted-tree.qza)

## Install packages
Open RStudio on your computer. If you have not already downloaded these packages, go to the lower right quadrant of your screen and open the Package tab. Click “download” and search for the package you want to download.

DECIPHER
adespatial
adespatial
ape
DESeq2
devtools
ggdendro
ggplot2
grid
gridExtra
knitr
MicrobeR
microbiomeSeq
pander
philr
phyloseq
plotly
png
qiime2R
ranacapa
tidyverse
vegan

To install most packages using the following command eg. **knitr**

```
install.packages("knitr")
```
For other packages, 
Copy and paste the following into your console.

```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
```
for example to install "phyloseq"

```
BiocManager::install("phyloseq")
```

For MicrobeR and qiime2R, use the following to install the packages

```
install_github("jbisanz/MicrobeR")
```

### Load Packages
The library command tells R to open the package you want to use. You need to do this every time you open R.
```
library(DECIPHER) # This package will help in importing, maintaining, analyzing, manipulating, and exporting a massive amount of sequences.
library(ape) # Analyses of Phylogenetics and Evolution package. Required for tree calculations to be used with phyloseq
library(DESeq2) # This package will help analyze "differential expression" in the microbiota alongside phyloseq
library(ggplot2) # Graphing package used in phyloseq. To edit the default setting of a plot, you need to use functions in this package.
library(phyloseq) # The phyloseq package seeks to address issues with multiple microbiome analysis packages by providing a set of functions that internally manage the organizing, linking, storing, and analyzing of phylogenetic sequencing data. In general, this package is used for UniFrac analyses.
library(plotly) # A package to create interactive web graphics of use in 3D plots
library(vegan) # The vegan package provides tools for descriptive community ecology. It has most basic functions of diversity analysis, community ordination and dissimilarity analysis. In general, this package is used for Bray-Curtis and Jaccard analyses.
library(philr) # This package provides functions for the analysis of compositional data 
library(tidyverse) # This package is designed to make it easy to install and load multiple 'tidyverse' packages in a single step
library(adespatial) # Tools for the multiscale spatial analysis of multivariate data
library(devtools) # Make package development easier by providing R functions that simplify and expedite common tasks
library(qiime2R) # A package for importing qiime artifacts into an R session
library(MicrobeR) # Data visualization
library(microbiome) # Data analysis and visualization
library(microbiomeSeq) # Data analysis and visualization
library("pander") # provide a minimal and easy tool for rendering R objects into Pandoc's markdown
library(ranacapa) # Data analysis 
library(grid) # support data visualization
library(gridExtra)  # support data visualization
library(knitr) # Provides a general-purpose tool for dynamic report generation in R using Literate Programming techniques.
library(png) # Figure download
library("ggdendro") #set of tools for dendrograms and tree plots using 'ggplot2'
library(ggpubr) # publication quality figures, based on ggplot2
library(RColorBrewer) # nice color options
```

### Load Data
In the code, the text before = is what the file will be called in R. Make this short but unique as this is how you will tell R to use this file in later commands.

Convert qiime artifacts directly to phyloseq
```
phy<-qza_to_phyloseq("table.qza", "rooted-tree.qza", "taxonomy.qza","Moving Pictures sample-metadata (QIIME 2 2018.4) - sample-metadata(1).tsv")
```
Or if you want to have more control over the object adding more or less data, you can make it yourself as below:
```
# Importing ASVs abundance file
ASVs<-read_qza("table.qza")
# Importing metadata
metadata<-read.table("Moving Pictures sample-metadata (QIIME 2 2018.4) - sample-metadata(1).tsv", , sep='\t', header=T, row.names=1, comment="")
metadata<-metadata[-1,]#remove the second line that specifies the data type
# Importing tree
tree<-read_qza("rooted-tree.qza")
# Importing taxonomy
taxonomy<-read_qza("taxonomy.qza")
tax_table<-do.call(rbind, strsplit(as.character(taxonomy$data$Taxon), "; "))
colnames(tax_table)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rownames(tax_table)<-taxonomy$data$Feature.ID

# Creating phyloseq object
physeq<-phyloseq(
  otu_table(ASVs$data, taxa_are_rows = T),
  phy_tree(tree$data),
  tax_table(tax_table),
  sample_data(metadata)
)
```

Summarizing the phyloseq object to check for feature of data 

```
# check for features of data  
summarize_phyloseq(physeq)
summary(sample_sums(physeq))
```

Rarefy the phyloseq object to even depth prior various analysis
```
physeq_rarefy <- rarefy_even_depth(physeq, rngseed=1, sample.size=0.9*min(sample_sums(physeq)), replace=F)
```

Check the taxa prevalence at Phylum level 

```
plot_taxa_prevalence(physeq, "Phylum")
```

Note:
Print version information about R, the OS and attached or loaded packages
```
sessionInfo()
```
## Composition plots

```
physeq_fam <- microbiome::aggregate_top_taxa(physeq, "Family", top = 10)
physeq.fam.rel <- microbiome::transform(physeq_fam, "compositional")
plot_composition(physeq.fam.rel,sample.sort = "BodySite", x.label = "SampleID") + theme(legend.position = "bottom") + scale_fill_brewer("Family", palette = "Paired") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Relative abundance") + theme(legend.title = element_text(size = 18))
```
## Alpha diversities
Alpha diversity measures are used to identify within individual taxa richness and evenness. The commonly used metrics/indices are Shannon, Inverse Simpson, Simpson, Gini, Observed and Chao1. These indices do not take into account the phylogeny of the taxa identified in sequencing. Phylogenetic diversity (Faith’s PD) uses phylogenetic distance to calculate the diversity of a given sample.

One has to consider the sequencing depth (how much of the taxa have been sampled) for each sample. If there is a large difference, then it is important to normalize the samples to equal sampling depth. First, we look at the sampling depth (no. of reads per sample).

We can plot the rarefaction curve for the observed ASVs in the entire data set. This is a way to check how has the richness captured in the sequencing effort.


```
ggrare(physeq, step = 50, color="BodySite", label = "Sample", se = TRUE)
```
phyloseq also allows you to easily plot alpha diversity, both by sample and by group

```
plot_richness(physeq_rarefy, measures="Shannon")
```
```
plot_richness(physeq_rarefy, x="BodySite", measures=c("Shannon", "simpson", "Observed"), color = "BodySite") + geom_boxplot() + theme_bw()
```
generate a csv file of the richness estimates using
```
richness <- estimate_richness(physeq_rarefy,measures=c("Shannon", "simpson", "Observed"))
write.csv(richness, file = "alpha_div.csv")
```

In addition to plotting you can also run anova test using the following option (it will also write the richness ina seperate csv file)
```
plot_anova_diversity(physeq, method = c("richness", "simpson", "shannon"), grouping_column = "BodySite", pValueCutoff = 0.05)
```

### Alpha statistics

Overall, for alpha-diversity:

* **ANOVA, t-test, or general linear models** with the normal distribution are used when the data is roughly normal
* **Kruskal-Wallis, Wilcoxon rank sum test, or general linear models** with another distribution are used when the data is not normal

To test for normalcy statistically, we can run the Shapiro-Wilk test of normality.
```
shannon <- estimate_richness(physeq,measures=c("Shannon"))
shapiro.test(shannon$Shannon)
hist(shannon$Shannon, main="Shannon diversity", xlab="", breaks=10)
```

#### Normally distributed metrics
we will use Shannon’s diversity as an example. we will test BodySite, which is a categorical variable with more than 2 levels. Thus, we run ANOVA. If age were only two levels, we could run a t-test

```
#Run the ANOVA and save it as an object
aov.shannon.bodysite <- aov(shannon$Shannon ~ BodySite, data=metadata)
#Call for the summary of that ANOVA, which will include P-values
summary(aov.shannon.bodysite)
```
To do all the pairwise comparisons between groups and correct for multiple comparisons, we run Tukey’s honest significance test of our ANOVA.

```
TukeyHSD(aov.shannon.bodysite)
```

#### Non-normally distributed metrics
we use Kruskal-Wallis (non-parametric equivalent of ANOVA). If we have only two levels, we would run Wilcoxon rank sum test (non-parametric equivalent of t-test)

```
kruskal.test(shannon$Shannon ~ BodySite, data=metadata)
```
We can test pairwise within the age groups with Wilcoxon Rank Sum Tests. This test has a slightly different syntax than our other tests
```
pairwise.wilcox.test(shannon$Shannon, metadata$BodySite, p.adjust.method="fdr")
```



