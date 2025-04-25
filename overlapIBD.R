##librairies
library("readr")
library("igraph")
library("dplyr")
library("ggplot2")
library("ggpattern")
library("GenomicRanges")
library("RColorBrewer")
library("intervals")
library("tidyr")

##Data (see examples in example/rawdata/)
hmmIBD_file <- "test/rawdata/IBD-WGS.hmm.txt.zip"
IBD_nodes_file <- "test/rawdata/WGS-nodes.tsv"
IBD_edges_file <- "test/rawdata/WGS-edges.tsv"
chromosomes_file <- "test/rawdata/chrom-sizes.tsv"
highlight_regions_file <- "test/rawdata/Nwakanma-2013-ConservedRegions.tsv" #set to NA if you want to ignore this
highlight_genes_file <- "test/rawdata/drug-resistance-genes.tsv" #set to NA if you want to ignore this

##Options
parentalmin <- 0.35 #Minimal IBD for parental-offspring relationship
parentalmax <- 0.65 #Maximal IBD for parental-offspring relationship
identical <- 0.9 #Minimal IBD to consider two samples identical
different <- 0.2 #Maximal IBD to consider two samples different
similar <- 0.5 #Minimal IBD to consider two samples related
lowercov <- 0.5 #Minimal percentage of IBD coverage for hotspots
uppercov <- 10 #Maximal percentage of IBD coverage for hotspots

##Output
outDir <- "test/out" #Output directory
density_table <- paste0(outDir, "/IBD-density.csv") #set to NA of not run already

##Steps
#Setup
if (!dir.exists(outDir)){
  dir.create(outDir, recursive = TRUE)
}
source("src/utils.R")
#

#Parental relationships inferred from IBD
source("src/find-triades.R")
source("src/view-triades.R")
#

#Overlapping IBD fragments
source("src/IBD-hotspots.R")
#







