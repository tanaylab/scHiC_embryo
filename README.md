# scHiCe9

This repository contains code for the analysis of single-cell Hi-C data from E9.5 mouse embryos.

To fully reproduce the results, you should:
1. Download the data
2. Run it through the scHi-C pipeline
3. Set up several additional misha repositories
4. Run the code in this repository.

Downloading the data
The raw Hi-C data for the project can be found in https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148793.
Please download the data using run ids, rather than sample ids, e.g:
fasterq-dump -O /path/to/download SRR24224552

scRNA-seq from E9 embryos is also available on GEO, but it is better if you use the UMI matrices directly, which will be referenced later. 
