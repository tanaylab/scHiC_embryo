# single-cell Hi-C in E9.5 mouse embryos

This repository contains code for the analysis of single-cell Hi-C data from E9.5 mouse embryos.

To fully reproduce the results, you should:
1. Download the data
2. Set up misha databases and run the lab's previous scHi-C pipeline
3. Run the code in this repository.

Note that for most of the cases, you can skip directly to step 3, as for most analyses we provide processed files that can be analysed directly, without requiring all the raw data.

## Downloading the data
The raw Hi-C data for the project can be found in https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148793.
Please download the data using run ids, rather than sample ids, e.g:
fasterq-dump -O /path/to/download SRR24224552

scRNA-seq from E9 embryos is also available on GEO, but it is better if you use the UMI matrices directly, which will be referenced later. 

## Misha databases and scHi-C pipeline
You can install and read on creating Misha databases here:
https://tanaylab.github.io/misha/

You can install and read on the scHi-C pipeline here:
https://github.com/tanaylab/schic2

Our analysis was performed on the mm9 reference, except for the ATAC analysis which requires an mm10 database.

## Running the code in this repository
Most users will prefer to use processed files we provide, instead of downloading the raw data, setting up misha databases, and running the scHi-C pipeline.
Note that some analyses still require the raw data, and have no processed files that aid in their execution.
The processed files can be found here:
TODO
Please download and extract the processed files. You should then edit the params.r file
