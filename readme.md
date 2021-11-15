This workshop is a follow-up of the Microbiome analysis using QIIME2 workshop. The result from the previous workshop will be used to demonstrate basic analyses of microbiota data to determine if and how communities differ by variables of interest using R. This all-day workshop will consist of lectures and hands on training to analyze from raw dataset through publication-quality statistics and visualizations


## 1. Get Started
### 1.1 Download and Install
* Base R: <http://cran.mtu.edu/>

NOTE: If you need to update to the most recent version of R on Windows you can do so using the installr package. Instructions here. For OSX and Ubuntu, download from CRAN using the link above.

* RStudio: <https://www.rstudio.com/products/rstudio/download3/>


## 2. Project folder: Stay organized
All of our analyses will be organized into a “Project”.

Make a new project by selecting File->New project. Select “New Directory” and “Empty Project”. Name the project “Microbiome_Analysis” and save the project to your Desktop. Place all of your files for this analysis in the folder created on the Desktop

Create a new R script (File->New file->R script) to save your code. This file will automatically be saved in the project folder.

Now your screen should look like this

* Upper left: Where you type and save the code you want to run.
* Upper right: Files you load into and create in R. To view one, click on it and it will open in the upper left pane.
* Lower left: The console. Where commands and outputs run (similar to the one mothur window).
* Lower right: Variable. Explore the different tabs.

## 3. Data description

In this workshop, We will work with data processed through Qiime2 and can be downloaded from the following links
* [Metadata](https://data.qiime2.org/2018.4/tutorials/moving-pictures/sample_metadata.tsv) 
* [Dada2 feature table](https://docs.qiime2.org/2018.4/data/tutorials/moving-pictures/table.qza)
* [taxonomy assignments](https://docs.qiime2.org/2018.4/data/tutorials/moving-pictures/taxonomy.qza)
* [midpointed-tree](https://docs.qiime2.org/2018.4/data/tutorials/moving-pictures/rooted-tree.qza)

## 4. Install packages
Open RStudio on your computer. If you have not already downloaded these packages, go to the lower right quadrant of your screen and open the Package tab. Click “download” and search for the package you want to download.

* DECIPHER  
* adespatial  
* ape  
* DESeq2  
* devtools  
* ggdendro
* ggplot2
* grid
* gridExtra
* knitr
* MicrobeR
* microbiomeSeq
* microbiome
* pander
* philr
* phyloseq
* plotly
* png
* qiime2R
* ranacapa
* tidyverse
* vegan
* ggpubr
* RColorBrewer
* microbiomeutilities

```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install(c("DECIPHER","DESeq2", "philr", "phyloseq"))

BiocManager::install("preprocessCore")
BiocManager::install("GO.db")
BiocManager::install("impute")

install.packages(c("adespatial","ape","devtools","ggdendro","gridExtra","knitr","MicrobeR","pander","plotly","png","tidyverse","vegan"))


library(devtools)

devtools::install_github("gauravsk/ranacapa")
devtools::install_github("umerijaz/microbiomeSeq") 

install.packages("remotes")

remotes::install_github("jbisanz/qiime2R")
remotes::install_github("jbisanz/MicrobeR")

New Packages:

remotes::install_github("microsud/microbiomeutilities")
BiocManager::install("microbiome")

install.packages(c("ggpubr", "RColorBrewer"))
```

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

### 4.1 Load Packages
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

library(microbiomeutilities) # some utility tools 
```

### 4.2 Load Data
In the code, the text before = is what the file will be called in R. Make this short but unique as this is how you will tell R to use this file in later commands.

Convert qiime artifacts directly to phyloseq
```
phy <- qza_to_phyloseq("table.qza", "rooted-tree.qza", "taxonomy.qza","sample_metadata.tsv")
```
Or if you want to have more control over the object adding more or less data, you can make it yourself as below:
```
# Importing ASVs abundance file
ASVs <- read_qza("table.qza")
# Importing metadata
metadata <- read.table("sample_metadata.tsv", , sep='\t', header=T, row.names=1, comment="")
metadata <- metadata[-1,] # remove the second line that specifies the data type
# Importing tree
tree <- read_qza("rooted-tree.qza")
# Importing taxonomy
taxonomy <- read_qza("taxonomy.qza")
tax_table <- do.call(rbind, strsplit(as.character(taxonomy$data$Taxon), "; "))
colnames(tax_table) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rownames(tax_table) <- taxonomy$data$Feature.ID

# Creating phyloseq object
physeq <- phyloseq(
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
print_ps(physeq)
summary(sample_sums(physeq))
```
Accessing the phyloseq object
```
ntaxa(physeq)

nsamples(physeq)

sample_names(physeq)[1:5]  

rank_names(physeq)  

sample_variables(physeq)  

otu_table(physeq)[1:5, 1:5]  

tax_table(physeq)[1:5, 1:4]
```

Distribution of reads

```
plot_read_distribution(physeq, groups = "BodySite", plot.type = "density") + theme_biome_utils()
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
## 5. Composition plots
Barplots are a one way of visualising the composition of your samples.
At Family level and relative abundance
```
physeq_fam <- microbiome::aggregate_rare(physeq, level = "Family", detection = 50/100, prevalence = 70/100)

physeq.fam.rel <- microbiome::transform(physeq_fam, "compositional")

physeq.fam.rel <- physeq %>%
  aggregate_rare(level = "Family", detection = 50/100, prevalence = 70/100) %>%
  microbiome::transform(transform = "compositional")

plot_composition(physeq.fam.rel,sample.sort = "BodySite", x.label = "SampleID") + theme(legend.position = "bottom") + scale_fill_brewer("Family", palette = "Paired") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Relative abundance") + theme(legend.title = element_text(size = 18))
```
Another option fot barplots 
```
taxa_barplot(Summarize.Taxa(ASVs$data, as.data.frame(tax_table))$Family, metadata, "BodySite")

# To make it interactive
ggplotly(taxa_barplot(Summarize.Taxa(ASVs$data, as.data.frame(tax_table))$Family, metadata, "BodySite"))

# save the plot
b.plot <- taxa_barplot(Summarize.Taxa(ASVs$data, as.data.frame(tax_table))$Family, metadata, "BodySite")

ggsave("barplot_family.png", b.plot,  width = 14, height = 10, dpi = 300)
```

## 6. Heatmap
These are a good alternative to barplots

```
taxa_heatmap(Summarize.Taxa(ASVs$data, as.data.frame(tax_table))$Family, metadata, "BodySite")
```

Another option
```
library(pheatmap)

p <- plot_taxa_heatmap(physeq,
                  subset.top = 25,
                  VariableA = c("BodySite","Subject"),
                  transformation = "log10",
                  cluster_rows = T,
                  cluster_cols = F,
                  show_colnames = F,
                  heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                )
#the plot is stored here
p$plot

# table used for plot is here
p$tax_tab[1:3,1:3]
```

Alternative option with more control and options
```
h.map <- plot_heatmap(physeq.fam.rel, method="PCoA", distance="bray", taxa.label = "Family", sample.order = unique(sample_names(physeq))) + facet_grid(~BodySite, scales = "free_x", drop = TRUE) + theme_bw() + theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1)) + theme(legend.key = element_blank(),strip.background = element_rect(colour="black", fill="white"))

# Make bacterial names italics
h.map <- h.map + theme(axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'))

# Change the color palette
h.map <- h.map + scale_fill_distiller("Abundance", palette = "RdYlBu")

# clean the x-axis
h.map <- h.map + rremove("x.text")

print(h.map)

# Saving the plot
ggsave("heatmap_family.png", h.map,  width = 14, height = 10, dpi = 300)
```

## 7. Boxplot

```
physeq_df <- microbiomeutilities::phy_to_ldf(physeq_fam, 
                                         transform.counts = "compositional")

# An additonal column Sam_rep with sample names is created for reference purpose
colnames(physeq_df)

# Box plot at Family level

ggstripchart(physeq_df, "BodySite", "Abundance", 
             facet.by = "Family", color = "BodySite",
             palette = "jco") + rremove("x.text")
```

Plot relative abundance of top taxa

```
mycols <- c("coral", "steelblue2", "slategray2", "olivedrab")

plot_taxa_boxplot(physeq,
                  taxonomic.level = "Phylum",
                  top.otu = 6, 
                  group = "BodySite",
                  add.violin= FALSE,
                  title = "Top three family", 
                  keep.other = FALSE,
                  group.order = c("gut","tongue","right palm", "left palm"),
                  group.colors = mycols,
                  dot.size = 2) + theme_biome_utils()
```

plotting taxa specified by the user

```
physeq.f <- format_to_besthit(physeq)

top_taxa(physeq.f, 5)

select.taxa <- c("d29fe3c70564fc0f69f2c03e0d1e5561:k__Bacteria", "154709e160e8cada6bfb21115acc80f5:ovatus")

p <- plot_listed_taxa(physeq.f, select.taxa, 
                 group= "BodySite",
                 group.order = c("gut","tongue","right palm", "left palm"),
                 group.colors = mycols,
                 add.violin = F,
                 dot.opacity = 0.25,
                 box.opacity = 0.25,
                 panel.arrange= "grid") + ylab("Relative abundance") + scale_y_continuous(labels = scales::percent)


# Adding statistical test with ggpubr::stat_compare_means()

# If more than two variables
comps <- make_pairs(sample_data(physeq.f)$BodySite)
print(comps)

p <- p + stat_compare_means(
  comparisons = comps,
  label = "p.format",
  tip.length = 0.05,
  method = "wilcox.test") 

p + scale_y_continuous(labels = scales::percent)
```

Plot top four genera

```
physeq.genus <- aggregate_taxa(physeq, "Genus")
top_four <- top_taxa(physeq.genus, 4)
top_four

top_genera <- plot_listed_taxa(physeq.genus, top_four, 
                 group= "BodySite",
                 group.order = c("gut","tongue","right palm", "left palm"),
                 group.colors = mycols,
                 add.violin = F,
                 dot.opacity = 0.25,
                 box.opacity = 0.25,
                 panel.arrange= "wrap")
top_genera

top_genera + stat_compare_means(
                   comparisons = comps,
                   label = "p.format",
                   tip.length = 0.05,
                   method = "wilcox.test")
```

#### Taxa overview
To check which taxa is present and how ti is distributed. 

#### Dominant taxa
Sometimes, we are interested in identifying the most dominant taxa in each sample. We may also wish to check what percent of samples within a given group are these taxa dominating.

```
physeq.gen <- aggregate_taxa(physeq,"Genus")

dom.tax <- dominant_taxa(physeq,level = "Genus", group="BodySite")
head(dom.tax$dominant_overview)
```
#### Taxa summary
On the enitre dataset

```
taxa_summary(physeq, "Phylum")
```

#### Taxa summary by groups
For group specific abundances of taxa

```
grp_abund <- get_group_abundances(physeq, 
                                  level = "Phylum", 
                                  group="BodySite",
                                  transform = "compositional")

# clean names 
grp_abund$OTUID <- gsub("p__", "",grp_abund$OTUID)
grp_abund$OTUID <- ifelse(grp_abund$OTUID == "", 
                          "Unclassified", grp_abund$OTUID)

mean.plot <- grp_abund %>% # input data
  ggplot(aes(x= reorder(OTUID, mean_abundance), # reroder based on mean abundance
             y= mean_abundance,
             fill=BodySite)) + # x and y axis 
  geom_bar(     stat = "identity", 
                position=position_dodge()) + 
  scale_fill_manual("BodySite", values=mycols) + # manually specify colors
  theme_bw() + # add a widely used ggplot2 theme
  ylab("Mean Relative Abundance") + # label y axis
  xlab("Phylum") + # label x axis
  coord_flip() # rotate plot 

mean.plot
```
#### Find samples dominated by specific taxa

```
bact_dom <- find_samples_taxa(physeq.gen, taxa = "g__Bacteroides")

#get sample dominated by "g__Bacteroides
ps.sub <- prune_samples(sample_names(physeq.gen) %in% bact_dom, physeq.gen)
```

## 8 Alpha diversities
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
a.div <- plot_richness(physeq_rarefy, x="BodySite", measures=c("Shannon", "simpson", "Observed"), color = "BodySite") + geom_boxplot() + theme_bw()

# adding statistical support
a.div + stat_compare_means(
  comparisons = comps,
  label = "p.signif",
  tip.length = 0.05,
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
    symbols = c("xxxx", "***", "**", "*", "ns")
  ),
  method = "wilcox.test")
```
generate a csv file of the richness estimates using
```
richness <- estimate_richness(physeq_rarefy,measures=c("Shannon", "simpson", "Observed"))
write.csv(richness, file = "alpha_div.csv")
```
Creating a plot with one index and stats
```
plot_diversity_stats(physeq, group = "BodySite", 
                     index = "diversity_shannon", 
                     group.order = c("gut","tongue","right palm", "left palm"),                      
                     group.colors = mycols,
                     label.format="p.format",
                     stats = TRUE) 
                     + ylab("Shannon Diversity") + xlab("")
```


### 8.1 Alpha statistics

Overall, for alpha-diversity:

* **ANOVA, t-test, or general linear models** with the normal distribution are used when the data is roughly normal
* **Kruskal-Wallis, Wilcoxon rank sum test, or general linear models** with another distribution are used when the data is not normal

To test for normalcy statistically, we can run the Shapiro-Wilk test of normality.
```
shannon <- estimate_richness(physeq,measures=c("Shannon"))
shapiro.test(shannon$Shannon)
hist(shannon$Shannon, main="Shannon diversity", xlab="", breaks=10)
```

#### 8.1.1 Normally distributed metrics
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

#### 8.1.2 Non-normally distributed metrics
we use Kruskal-Wallis (non-parametric equivalent of ANOVA). If we have only two levels, we would run Wilcoxon rank sum test (non-parametric equivalent of t-test)

```
kruskal.test(shannon$Shannon ~ BodySite, data=metadata)
```
We can test pairwise within the age groups with Wilcoxon Rank Sum Tests. This test has a slightly different syntax than our other tests
```
pairwise.wilcox.test(shannon$Shannon, metadata$BodySite, p.adjust.method="fdr")
```

## 9. Beta diversity metrices
Beta-diversity: Measures for differences between samples from different groups to identify if there are differences in the overall community composition and structure.

#### 9.1 Non - Phylogenetic beta diversity metrics
```
physeq.ord <- ordinate(physeq_rarefy, "PCoA", "bray")
b.div.bray <- plot_ordination(physeq_rarefy, physeq.ord, type= "samples", color= "BodySite") + geom_point(size=3)
b.div.bray <- b.div.bray + stat_ellipse() + ggtitle("Bray Curtis")  + theme_classic() + scale_color_brewer("Location", palette = "Set2")
print(b.div.bray)
```
#### 9.2 Phylogenetic beta diversity metrics
Weighted Unifrac will consider the abundances of different taxa.
```
# convert ot relative abundance
physeq_rel <- microbiome::transform(physeq, "compositional")
physeq.ord.wuni <- ordinate(physeq_rel, "PCoA", "unifrac", weighted=T)
b.div.wuni <- plot_ordination(physeq_rel, physeq.ord.wuni, type= "samples", color= "BodySite") + geom_point(size=3)
b.div.wuni <- b.div.wuni + stat_ellipse() + ggtitle("Weighted Unifrac")  + theme_classic() + scale_color_brewer("Location", palette = "Set2")
print(b.div.wuni)
```
### 9.3 Permanova
Permutational multivariate analysis of variance [further reading](https://onlinelibrary.wiley.com/doi/10.1002/9781118445112.stat07841)

```
otu <- abundances(physeq_rel)
meta <- meta(physeq_rel)


#Statistics - Bdiv
permanova <- adonis(t(otu) ~ BodySite, data = meta, permutations=99, method = "bray")

#P-value
print(as.data.frame(permanova$aov.tab)["BodySite", "Pr(>F)"])
```
#### 9.4 Checking the homogeneity condition
More infromation can be found by typing ```?betadisper```
```
#Pair - wise stats
dist <- vegdist(t(otu))
anova(betadisper(dist, meta$BodySite))

permutest(betadisper(dist, meta$BodySite), pairwise = TRUE)
```
## 10. Core microbiota

Subset the data to keep only gut samples.

```
physeq.gut <- subset_samples(physeq, BodySite == "gut")

# convert to relative abundance  
physeq.gut.rel <- microbiome::transform(physeq.gut, "compositional")

physeq.gut.rel2 <- prune_taxa(taxa_sums(physeq.gut.rel) > 0, physeq.gut.rel)
```

Check for the core ASVs

```
core.taxa.standard <- core_members(physeq.gut.rel2, detection = 0.001, prevalence = 50/100)
print(core.taxa.standard)
```

we only see IDs, not very informative. We can get the classification of these as below.

```
# Extract the taxonomy table
taxonomy_core <- as.data.frame(tax_table(physeq.gut.rel2))

# Subset this taxonomy table to include only core OTUs
core_taxa_id <- subset(taxonomy_core, rownames(taxonomy_core) %in% core.taxa.standard)

DT::datatable(core_taxa_id)
```
### 10.1 Core abundance and diversity
Total core abundance in each sample (sum of abundances of the core members):

```
core.abundance <- sample_sums(core(physeq.gut.rel2, detection = 0.001, prevalence = 50/100))

DT::datatable(as.data.frame(core.abundance))
```
#### 10.2 Core visualization - Heatmap

```
# Core with compositionals:
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

# Also define color palette
gray <- gray(seq(0,1,length=5))
p.core <- plot_core(physeq.gut.rel2, 
                    plot.type = "heatmap", 
                    colours = rev(brewer.pal(5, "Spectral")),
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = .5) +
  xlab("Detection Threshold (Relative Abundance (%))")
print(p.core)
```
Use the format_to_besthit function from microbiomeutilities to get the best classification of the ASVs.
```
physeq.gut.rel2.f <- microbiomeutilities::format_to_besthit(physeq.gut.rel2)

p.core <- plot_core(physeq.gut.rel2.f, 
                    plot.type = "heatmap", 
                    colours = rev(brewer.pal(5, "Spectral")),
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = .5) + 
  xlab("Detection Threshold (Relative Abundance (%))")

p.core + theme(axis.text.y = element_text(face="italic"))
```

### Abundance-Prevalence relationship

```
library(ggrepel)

physeq_comp <- microbiome::transform(physeq, "compositional")
# select gut
physeq_comp_gut <- subset_samples(physeq_comp, BodySite=="gut")

physeq_comp_gut <- core(physeq_comp_gut,detection = 0.0001, prevalence = 0.50) # reduce size for example

physeq_comp_gut <- format_to_besthit(physeq_comp_gut)

set.seed(163897)
prevalance <- plot_abund_prev(physeq_comp_gut, 
                       label.core = TRUE,
                       color = "Phylum", # NA or "blue"
                       mean.abund.thres = 0.01,
                       mean.prev.thres = 0.99,
                       dot.opacity = 0.7,
                       label.size = 3,
                       label.opacity = 1.0,
                       nudge.label=-0.15,
                       bs.iter=999, # increase for actual analysis e.g. 999
                       size = 20, # increase to match your nsamples(asv_ps)
                       replace = TRUE,
                       label.color="#5f0f40") 
prevelance <- prevalance + 
  geom_vline(xintercept = 0.95, lty="dashed", alpha=0.7) + 
  geom_hline(yintercept = 0.01,lty="dashed", alpha=0.7) +
  scale_color_brewer(palette = "Dark2")

prevelance
```


## 11. Tree plot

Subset the data to top 50 taxa for beteer visualization
```
physeq_top_50 <- subset_taxa(physeq, Kingdom=="k__Bacteria")
physeq_top_50 <- prune_taxa(names(sort(taxa_sums(physeq_top_50),TRUE)[1:50]), physeq_top_50)

# Default plot tree
plot_tree(physeq_top_50)

# Add genus labels to the tree and bootstrap values
plot_tree(physeq_top_50, label.tips="Genus", ladderize="left")

# Remove bootstrap values
plot_tree(physeq_top_50, nodelabf=nodeplotblank, label.tips="Genus", ladderize="left")

# Color the nodes by category
plot_tree(physeq_top_50, nodelabf=nodeplotblank, label.tips="Genus", ladderize="left", color="BodySite")

# Add size by abundance
plot_tree(physeq_top_50, nodelabf=nodeplotblank, label.tips="Genus", ladderize="left", color="BodySite", size="abundance")

# Convert to radial tree
plot_tree(physeq_top_50, nodelabf=nodeplotblank, label.tips="Genus", ladderize="left", color="BodySite") + coord_polar(theta="y")
```

### 11.1 Hierarchical Clustering

Heirachial clustering to visualize the distance between the samples using Weighted Unifrac and UPGMA method.

```
phy.hclust <- hclust(UniFrac(physeq_rarefy, weighted = TRUE), method="average")

ggdendrogram(phy.hclust, rotate = TRUE, theme_dendro = TRUE)
```

## 12. Microbiome network

You can plot the distances between ASVs as a network.

```
plot_net(physeq_rel, maxdist = 0.8, color = "BodySite")
#change distance to Jaccard
plot_net(physeq_rel, maxdist = 0.8, color = "BodySite", distance="jaccard")
```

### 12.1 igraph-based network

```
ig <- make_network(physeq_rel, max.dist=0.8)
plot_network(ig, physeq_rel)

# Add color label 
plot_network(ig, physeq_rel, color="BodySite", line_weight=0.4, label=NULL)

#replace the Jaccard (default) distance method with Bray-Curtis
ig <- make_network(physeq_rel, dist.fun="bray", max.dist=0.8)
plot_network(ig, physeq_rel, color="BodySite", line_weight=0.4, label=NULL)
```
Note: For co-occurrence networks of OTUs, I suggest trying Gephi or Cytoscape

## 13. Differential abundance testing

DeSeq2 to test for differential abundance between categories

```
#Convert phyloseq object ot DeSeq
bsdds <- phyloseq_to_deseq2(physeq_rarefy, ~ BodySite)
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <- apply(counts(bsdds), 1, gm_mean)
bsdds <- estimateSizeFactors(bsdds, geoMeans = geoMeans)
# DeSeq function tests for differential abundance 
bsdds <- DESeq(bsdds, test="Wald", fitType="parametric")

# Results function call creates a table of the results of the tests
res <- results(bsdds, cooksCutoff = FALSE)
alpha <- 0.01
sigtab <- res[which(res$padj < alpha), ]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(physeq_rarefy)[rownames(sigtab), ], "matrix"))
head(sigtab)

# Cleaning up the table a little for legibility
posigtab <- sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab <- posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

# Bar plot showing the log2-fold-change, showing Genus and Phylum. Uses some ggplot2 commands
sigtabgen <- subset(sigtab, !is.na(Genus))
# Phylum order
x <- tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x <- sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x <- tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
```

