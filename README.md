# overlapIBD
Find overlapping IBD segments to retrieve parental offspring relationships or IBD hotspots in the population.

## Run

The pipeline can be configured and run from `overlapIBD.R`.

## Raw data

Examples of these files are available in example/rawdata.

### *hmmIBD_file*

ZIP archive containing the coordinates of the identical and different IBD segments between pairs of samples.

### *IBD_nodes_file*

Tabulated file containing sample metadata.

### *IBD_edges_file*

Tabulated file containing the value of IBD between each pair of samples.

### *chromosomes_file*

Tabulated file containing the size of each chromosome. The chromosome ID must be a number.

### *highlight_regions_file*

Tabulated file containing regions to highlight in the IBD hotspot plot. The chromosome ID must be a number.

### *highlight_genes_file*

Tabulated file containing genes to highlight in the IBD hotspot plot. The chromosome ID must be a number.
