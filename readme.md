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
library(microbiomeSeq) # Data analysis and visualization
library("pander") # provide a minimal and easy tool for rendering R objects into Pandoc's markdown
library(ranacapa) # Data analysis 
library(grid) # support data visualization
library(gridExtra)  # support data visualization
library(knitr) # Provides a general-purpose tool for dynamic report generation in R using Literate Programming techniques.
library(png) # Figure download
library("ggdendro") #set of tools for dendrograms and tree plots using 'ggplot2'
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

# 1. First steps

### 1.1 Format metadata file

You can format the metadata file in a spreadsheet and validate the format using a google spreadsheet plugin 'Keemei' as described on the QIIME2 website. The only required column is the sample id column, which should be first. All other columns should correspond to sample metadata. Below, we use the default filename of metadata.txt.

```
wget -O "sample-metadata.tsv" "https://data.qiime2.org/2019.1/tutorials/atacama-soils/sample_metadata.tsv"
```
  
### 1.2 Activate conda environment

You should run this workflow in a conda environment, which makes sure the correct version of the Python packages required by QIIME2 are being used. You can activate this conda environment with this command:

```
source activate qiime2-2019.1
```

## 2 Paired-end read analysis commands
To analyze these data, the sequences that you just downloaded must first be imported into an qiime artifact of type EMPPairedEndSequences.

```
qiime tools import --type EMPPairedEndSequences --input-path raw_data --output-path raw_data/paired-end-sequences.qza
```

You next can demultiplex the sequence reads. This requires the sample metadata file, and you must indicate which column in that file contains the per-sample barcodes. In this case, that column name is BarcodeSequence. In this data set, the barcode reads are the reverse complement of those included in the sample metadata file, so we additionally include the --p-rev-comp-mapping-barcodes parameter. After demultiplexing, we can generate and view a summary of how many sequences were obtained per sample.

```
qiime demux emp-paired --m-barcodes-file sample-metadata.tsv --m-barcodes-column BarcodeSequence --i-seqs raw_data/paired-end-sequences.qza --o-per-sample-sequences demux.qza --p-rev-comp-mapping-barcodes

qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
```

### 2.1 Importing already demultiplexed data
   
It's assumed that the input is raw paired-end MiSeq data in demultiplexed FASTQ format located within a folder called raw_data. The filenames are expected to look like this: 105CHE6WT_S325_L001_R1_001.fastq.gz, where each field separated by an _ character corresponds to: (1) the sample name, (2) the sample# on the run, (3) the lane number, (4) the read number (R1 = forward, R2 = reverse), and (5) the set number. You may need to re-name your files to match this format; however, QIIME2 accepts many different input formats so the format and naming scheme you're using may also be supported.

```
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path demux/data --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux/demux-paired-end.qza    

qiime demux summarize --i-data demux/demux-paired-end.qza --o-visualization demux/demux-paired-end.qzv
```
### 2.2 transferring qzv files between server and your system

* open a new terminal in your laptop and write the following scp command to transfer files between server and laptop only for linux and Mac

change the xxx to match your id's

```
scp -r studentxx@sumo.biotech.wisc.edu:/home/BIOTECH/xxxxx/qiime2_tutorial/demux/demux-paired-end.qzv /xxx/xxxx/
```
Those on Windows using MobaXterm (listed above) will have a built-in GUI for transfers.


## 3 Running DADA2 workflow :

DADA2 is a pipeline for detecting and correcting (where possible) Illumina amplicon sequence data. This quality control process will additionally filter any phiX reads (commonly present in marker gene Illumina sequence data) that are identified in the sequencing data, and will filter chimeric sequences.
To determine what values to pass to remove low quality regions of the sequences, you should review the Interactive Quality Plot tab in the demux.qzv file that was generated by qiime demux summarize above.

if you interested in getting to know more about the DADA2, [here](https://benjjneb.github.io/dada2/tutorial.html) is a tutorial from the developer

```
mkdir Analysis

qiime dada2 denoise-paired --i-demultiplexed-seqs demux/demux-paired-end.qza --p-trim-left-f 13 --p-trim-left-r 13 --p-trunc-len-f 150 --p-trunc-len-r 150 --p-n-threads 2 --o-table Analysis/table.qza --o-representative-sequences Analysis/rep-seqs.qza --o-denoising-stats Analysis/denoising-stats.qza
```

```
cp /mnt/software/workshop/data/qiime2/qiime2_tutorial/Analysis/*.qza Analysis/.
```

### 3.1 Visualizing Denoise stats

```
qiime metadata tabulate --m-input-file Analysis/denoising-stats.qza --o-visualization Analysis/denoising-stats.qzv
```

### Optional
Generate qzv files for table and representative sequences

```
qiime feature-table summarize --i-table Analysis/table.qza --o-visualization Analysis/table.qzv --m-sample-metadata-file sample-metadata.tsv            

qiime feature-table tabulate-seqs --i-data Analysis/rep-seqs.qza --o-visualization Analysis/rep-seqs.qzv
```

### 3.2 Convert stats QZA to TSV

You should take a look at the read count table to see how many reads were retained at each step of the DADA2 pipeline:

```
qiime tools export --input-path Analysis/denoising-stats.qza --output-path Analysis/dada2_output
```

## 4. Filtering resultant table

### 4.1 Summarize output table

Once a denoising pipeline has been run you can summarize the output table with the below command, which will create a visualization artifact for you to view. You should use this visualization to decide how to filter your table and also if you want to rarefy your data (which simplifies your analyses, but throws out data!) you can choose a cut-off based on this file.

```
qiime feature-table summarize --i-table Analysis/table.qza --o-visualization Analysis/table_summary.qzv
```

### 4.2 Total-frequency-based filtering

```
mkdir Analysis/filtered
```

Total-frequency-based filtering is used to filter samples or features based on how frequently they are represented in the feature table. (X is a placeholder for your choice)

```
qiime feature-table filter-samples --i-table Analysis/table.qza --p-min-frequency X --o-filtered-table Analysis/filtered/sample-frequency-filtered-table.qza
```
you can make a summary of the filtered abundance table:    

```
qiime feature-table summarize --i-table Analysis/filtered/sample-frequency-filtered-table.qza --o-visualization Analysis/filtered/sample-frequency-filtered-table_summary.qzv
```

### 4.3 Filter out rare ASVs

Based on the above summary visualization you can choose a cut-off for how frequent a variant needs to be (and optionally how many samples need to have the variant) for it to be retained. One possible choice would be to remove all ASVs that have a frequency of less than 0.1% of the mean sample depth. This cut-off excludes ASVs that are likely due to MiSeq bleed-through between runs (reported by Illumina to be 0.1% of reads). To calculate this cut-off you would identify the mean sample depth in the above visualization, multiply it by 0.001, and round to the nearest integer.

Once you've determined how you would like to filter your table you can do so with this command (X is a placeholder for your choice):

```
qiime feature-table filter-features --i-table Analysis/filtered/sample-frequency-filtered-table.qza --p-min-frequency 4 --p-min-samples 1 --o-filtered-table Analysis/filtered/table_filt.qza     
```

Since the ASVs will be in a separate FASTA file you can exclude the low frequency ASVs from the sequence file with this command:

```
qiime feature-table filter-seqs --i-data Analysis/rep-seqs.qza --i-table Analysis/filtered/table_filt.qza --o-filtered-data Analysis/filtered/rep_seqs_filt.qza
```

Finally, you can make a new summary of the filtered abundance table:

```
qiime feature-table summarize --i-table Analysis/filtered/table_filt.qza --o-visualization Analysis/filtered/table_filt_summary.qzv
```

### 5. Build quick phylogeny with FastTree

### 5.1 Making multiple-sequence alignment

We'll need to make a multiple-sequence alignment of the ASVs before running FastTree. First, we'll make a folder for the output files.

```
mkdir Analysis/phylogeny
```

We'll use MAFFT to make a de novo multiple-sequence alignment of the ASVs.

```
qiime alignment mafft --i-sequences Analysis/filtered/rep_seqs_filt.qza --o-alignment Analysis/phylogeny/rep_seqs_filt_aligned.qza
```

### 5.2 Filtering multiple-sequence alignment

Variable positions in the alignment need to be masked before FastTree is run, which can be done with this command:

```
qiime alignment mask --i-alignment Analysis/phylogeny/rep_seqs_filt_aligned.qza --o-masked-alignment Analysis/phylogeny/rep_seqs_filt_aligned_masked.qza
```

### 5.3 Running FastTree

Finally FastTree can be run on this masked multiple-sequence alignment:

```
qiime phylogeny fasttree --i-alignment Analysis/phylogeny/rep_seqs_filt_aligned_masked.qza --o-tree Analysis/phylogeny/rep_seqs_filt_aligned_masked_tree
```

### 5.4 Add root to tree

FastTree returns an unrooted tree. One basic way to add a root to a tree is to add it add it at the midpoint of the largest tip-to-tip distance in the tree, which is done with this command:

```
qiime phylogeny midpoint-root --i-tree Analysis/phylogeny/rep_seqs_filt_aligned_masked_tree.qza --o-rooted-tree Analysis/phylogeny/rep_seqs_filt_aligned_masked_tree_rooted.qza
```

The ASV/OTU phylogeny tree can be found in the phylogeny folder in the file rep_seqs_filt_aligned_masked_tree_rooted.qza and can be directly uploaded into [iTOL website](https://itol.embl.de) to visualize the tree

## 6. Generate rarefaction curves

A key quality control step is to plot rarefaction curves for all of your samples to determine if you performed sufficient sequencing. The below command will generate these plots (X is a placeholder for the maximum depth in your dataset).

```
qiime diversity alpha-rarefaction --i-table Analysis/filtered/table_filt.qza --p-max-depth X --p-steps 20 --i-phylogeny Analysis/phylogeny/rep_seqs_filt_aligned_masked_tree_rooted.qza --m-metadata-file sample-metadata.tsv --o-visualization Analysis/rarefaction/rarefaction_curves.qzv
```

For some reason, the QIIME2 default in the above curves with the metadata file (which you can see in the visualization) is to not give you the option of seeing each sample's rarefaction curve individually (even though this is the default later on in stacked barplots!), only the "grouped" curves by each metadata type. As it can be quite important in data QC to see if you have inconsistent samples, we need to rerun the above command, but this time omitting the metadata file (use the same X for the maximum depth, as above).

```
qiime diversity alpha-rarefaction --i-table Analysis/filtered/table_filt.qza --p-max-depth X --p-steps 20 --i-phylogeny Analysis/phylogeny/rep_seqs_filt_aligned_masked_tree_rooted.qza --o-visualization Analysis/rarefaction/rarefaction_curves_each_curve.qzv
```

## 7. Assign taxonomy

You can assign taxonomy to your ASVs using a Naive-Bayes approach implemented in the scikit learn Python library and the SILVA database. Note that we have trained classifiers, but you will need to generate your own if your region of interest isn't there. 

**Do not run:** Download a pre trained reference database

```
wget https://data.qiime2.org/2019.4/common/silva-132-99-515-806-nb-classifier.qza
```
```
qiime feature-classifier classify-sklearn --i-reads Analysis/filtered/rep_seqs_filt.qza --i-classifier /mnt/software/workshop/data/qiime2/qiime2_tutorial/silva-132-99-515-806-nb-classifier.qza --output-dir Analysis/taxa 
```

As with all QZA files, you can export the output file to take a look at the classifications and confidence scores:

```
qiime tools export --input-path Analysis/taxa/classification.qza --output-path Analysis/taxa
```

A more useful output is the interactive stacked bar-charts of the taxonomic abundances across samples, which can be output with this command:

```
qiime taxa barplot --i-table Analysis/filtered/table_filt.qza --i-taxonomy Analysis/taxa/classification.qza --m-metadata-file sample-metadata.tsv --o-visualization Analysis/taxa/taxa_barplot.qzv
```

Now the QIIME2 default in the above barplots is to plot each sample individually (opposite of what it was for the rarefaction curves!). However, in the visualizer you can label the samples according to their metadata types, but the plots are not summed together according to their metadata (as we could do in QIIME1 using the summarize_taxa_through_plots.py script), so we have to rerun the above command after running a new command first to group the samples by metadata - this creates a new feature table (you have to do this for each metadata category) which becomes the new input for the above same taxa barplot command (again, run each time for each metadata category). Note that CATEGORY here below is a placeholder for the text label of your category of interest from the metadata file.

```
qiime feature-table group --i-table Analysis/filtered/table_filt.qza --p-axis sample --p-mode sum --m-metadata-file sample-metadata.tsv --m-metadata-column TransectName --o-grouped-table Analysis/filtered/table_filt_transectname.qza
```

```
qiime taxa barplot --i-tabnalysis/filtered/table_filt_transectname.qza --i-taxonomy Analysis/taxa/classification.qza --m-metadata-file sample-metadata.tsv --o-visualization Analysis/taxa/taxa_barplot_transectname.qzv
```

### 7.1 Taxonomy-based filtering of tables and sequences
Taxonomy-based filtering is a very common type of feature-metadata-based filtering

```
qiime taxa filter-table --i-table Analysis/filtered/table_filt.qza --i-taxonomy Analysis/taxa/classification.qza --p-exclude mitochondria,chloroplast --o-filtered-table Analysis/filtered/table-no-contam.qza

qiime taxa filter-seqs --i-sequences Analysis/filtered/rep_seqs_filt.qza --i-taxonomy Analysis/taxa/classification.qza --p-exclude mitochondria,chloroplast --o-filtered-sequences Analysis/filtered/seq-no-contam.qza
```

## 8. Calculating diversity metrics and generating ordination plots

Common alpha and beta-diversity metrics can be calculated with a single command in QIIME2. In addition, ordination plots (such as PCoA plots for weighted UniFrac distances) will be generated automatically as well. This command will also rarefy all samples to the sample sequencing depth before calculating these metrics (X is a placeholder for the lowest reasonable sample depth; samples with depth below this cut-off will be excluded).

```
qiime diversity core-metrics-phylogenetic --i-table Analysis/filtered/table_filt.qza --i-phylogeny Analysis/phylogeny/rep_seqs_filt_aligned_masked_tree_rooted.qza --p-sampling-depth 971 --m-metadata-file sample-metadata.tsv --output-dir Analysis/diversity
```

You can then produce an alpha-diversity Shannon metric boxplot comparing the different categories in your metadata file using the following command (note that all the individual Shannon values for each sample are inside the shannon_vector.qza input file, if you wish to see them, and QIIME2 has many different metrics if you wish to use them).

```
qiime diversity alpha-group-significance --i-alpha-diversity Analysis/diversity/shannon_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization Analysis/diversity/shannon_compare_groups.qzv
```

Next we’ll analyze sample composition in the context of categorical metadata using PERMANOVA using the `beta-group-significance` command. The following commands will test whether distances between samples within a group, are more similar to each other than they are to samples from the other groups. If you call this command with the -`-p-pairwise parameter`, as we’ll do here, it will also perform pairwise tests that will allow you to determine which specific pairs of groups (e.g., tongue and gut) differ from one another. 

```
qiime diversity beta-group-significance --i-distance-matrix Analysis/diversity/unweighted_unifrac_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column TransectName --o-visualization Analysis/diversity/unweighted-unifrac-Transect-Name-significance.qzv --p-pairwise
```

## 9 Exporting the final abundance, and sequence files

Lastly, to get the BIOM file (with associated taxonomy), and FASTA file (one per ASV) for your dataset to plug into other programs you can use the commands below. Note that there are a couple of file manipulations included here since the conversion scripts expect certain slightly different headers than those produced by QIIME2 by default:

### 9.1 For FASTA:

```
qiime tools export --input-path Analysis/filtered/rep_seqs_filt.qza --output-path Analysis/filtered/rep_seq_exported
```

### 9.2 For BIOM w/taxonomy:

```
sed -i -e '1 s/Feature/#Feature/' -e '1 s/Taxon/taxonomy/' Analysis/taxa/taxonomy.tsv                

qiime tools export --input-path Analysis/filtered/table_filt.qza --output-path Analysis/filtered/table_filt_exported       

biom add-metadata -i Analysis/filtered/table_filt_exported/feature-table.biom -o Analysis/filtered/feature-table_w_tax.biom --observation-metadata-fp Analysis/taxa/taxonomy.tsv --sc-separated taxonomy         

biom convert -i Analysis/filtered/feature-table_w_tax.biom -o Analysis/filtered/feature-table_w_tax.txt --to-tsv --header-key taxonomy           
```

### 9.3 Relative abundance at different levels

```
mkdir Analysis/rel_table
```
Instead of $levels use either 1,2,3,4,5,6,7 for different taxonomic levels

```
qiime taxa collapse --i-table Analysis/filtered/table_filt.qza --i-taxonomy Analysis/taxa/classification.qza --p-level $levels --o-collapsed-table Analysis/rel_table/phyla-table_L$levels.qza     

qiime feature-table relative-frequency --i-table Analysis/rel_table/phyla-table_L$levels.qza --o-relative-frequency-table Analysis/rel_table/rel-phyla-table_L$levels.qza
      
qiime tools export --input-path Analysis/rel_table/rel-phyla-table_L$levels.qza --output-path Analysis/rel_table/rel-table_L$levels      
    
biom convert -i Analysis/rel_table/rel-table_L$levels/feature-table.biom -o Analysis/rel_table/rel-table_L$levels/rel-phyla-table_L$levels.tsv --to-tsv        
```

## Differential abundance testing
ANCOM can be applied to identify features that are differentially abundant (i.e. present in different abundances) across sample groups. As with any bioinformatics method, you should be aware of the assumptions and limitations of ANCOM before using it. 

```
mkdir ancom
```

We’ll start by creating a feature table that contains only the Site name - BAQ2420.

```
qiime feature-table filter-samples --i-table Analysis/filtered/table_filt.qza --m-metadata-file sample-metadata.tsv --p-where "SiteName='BAQ2420'" --o-filtered-table Analysis/ancom/BAQ2420-table.qza
```

adding pusedocount
```
qiime composition add-pseudocount --i-table Analysis/ancom/BAQ2420-table.qza --o-composition-table Analysis/ancom/comp-BAQ2420-table.qza
```

We can then run ANCOM on the `Description` column to determine what features differ in abundance across the Site name - BAQ2420   

```
qiime composition ancom --i-table Analysis/ancom/comp-BAQ2420-table.qza --m-metadata-file sample-metadata.tsv --m-metadata-column Description --o-visualization Analysis/ancom/ancom-Description.qzv
```

### classifiers

We can use a random forest classifier directly within QIIME via the qiime sample-classifier tool.

```
qiime sample-classifier classify-samples \
  --i-table Analysis_filtered/table_filt.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column TransectName \
  --p-optimize-feature-selection \
  --p-parameter-tuning \
  --p-estimator RandomForestClassifier \
  --p-n-estimators 100 \
  --o-visualization attacamasoil-transectname.qzv
```

### regressors

For continuous data, you can use a random forest regressor.

Important: continuous data could be age, weight, etc, for each of your samples, but it could also be quantitative levels of some metabolite, host gene expression, etc.

```
qiime sample-classifier regress-samples \
  --i-table Analysis_filtered/table_filt.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column Elevation \
  --p-optimize-feature-selection \
  --p-parameter-tuning \
  --p-estimator RandomForestRegressor \
  --p-n-estimators 100 \
  --o-visualization Analysis/atacamsoil_elevation.qzv
```
