# single-cell Hi-C in E9.5 mouse embryos

This repository contains code for the analysis of single-cell Hi-C data from E9.5 mouse embryos.

To fully reproduce the results, you should:
1. Download the data
2. Set up misha databases and run the lab's previous scHi-C pipeline
3. Run the code in this repository.

Note that for most of the cases, you can skip directly to step 3, as for most analyses we provide preprocessed files that can be analysed directly, without requiring all the raw data.

## Download the data
The raw Hi-C data for the project can be found in https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148793.

Please download the data using run ids, rather than sample ids, e.g:  
fasterq-dump -O /path/to/download SRR24224552

scRNA-seq from E9 embryos is also available on GEO, but it is better if you use the UMI matrices directly, which will be referenced later. 

## Misha databases and scHi-C pipeline
You can install and read on creating Misha databases here:  
https://tanaylab.github.io/misha/

You can install and read on the scHi-C pipeline here:  
https://github.com/tanaylab/schic2

Our analysis was performed on the mm9 reference, except for the ATAC analysis which was performed on an mm10 database.

## Running the code in this repository
Most users will prefer to use preprocessed files we provide, instead of downloading the raw data, setting up misha databases, and running the scHi-C pipeline.  
Note that some analyses still require the raw data, and have no preprocessed files that aid in their execution.

The preprocessed files can be found here:
TODO  
They contain two root directories - *data* and *metacell*, which contain data files and the *metacell* database directory. These directories will be referenced in the *params.r* file. 

### Configuration file
You should then edit the following variables in *params.r*:
*SCHIC.MISHA.PATH* - the path to the misha mm9 database. This database should include a track for each single Hi-C cell, as created using the scHi-C pipeline.
*SHAMAN.MISHA.PATH* - a path to another misha database. This database will include tracks for pools of single cells (across cells in the same cluster), and the shaman score tracks. This can be the same as the previous mm9 database, and we allow different databases just for performance reasons.
*ATAC.MISHA.PATH* - an mm10 misha database, used just for the ATAC analysis.
*DATA.DIR* - path to the directory containing the preprocessed data files (the directory "data" in the preprocessed files).
*FIGURE.DIR* - path to the directory where the figures will be saved.
*EXPRESSION.MC.DIR* - path to the directory containing the preprocessed metacell files (the directory "metacell" in the preprocessed files)..

### Additional misha tracks
While not required for most analyses, there are several tracks that should exist in the *SCHIC.MISHA.PATH* database:

Hi-C tracks for the single-cells: for the E9.5 data, and the E10.5 pEry confirmation data. These tracks are created by the scHi-C pipeline.  
ChIP tracks for tal1 and gata1, called *tal1_bigwig* and *gata1_bigwig*.  
Hi-C tracks for the E14.5 HSC and HPC data (called *hic.SC.bonev_npc.NPC_Bonev* and *hic.SC.chen_hsc.HSC_HSC*).  
ChIP tracks for h3k4me1 and h3k27me3 from ENCDOE of 5 different tissues. The track names are *ENCODE.bing_ren.e10_5.X.h3k4me1_rep1* and *ENCODE.bing_ren.e10_5.X.h3k27me3_sum*, where X is one of forebrain, heart, hindbrain, limb and midbrain.  

### Brief explanation of the code files
*initial_clustering* - this directory contains several files with functions used for the initial analysis of the data - the clustering into C1-3.  
*initial_hic_clustering.r* - this is the script which uses the functions in *initial_clustering*.  
*create_figs.r* - contains function to load data (*load.decay.metrics* and *load.ab.and.cov*) which should be called before other functions are called, and *create_figs* which creates most of figures 1-4 and their corresponding EDFs.  
*esc_analysis.r* - the analysis comparing the embryo-proper cells to the ESCs (mostly figure 1).    
*ery_analysis.r* - the analysis comparing the embryo-proper cells to the pErys (mostly figure 2).  

Files relating to the replication model and its application on the embryo-proper cells:
*repl_model_fit.r* - the core code for optimizing the replication model.  
*repl_model_validation.r* - code for validation of the replication model - simulations and cross validation to choose meta parameters.  
*method_comparison.r* - code for calling other methods for the analysis of scHi-C data for comparison with the replication model. 
*repl_model.r* - code that is used to actually apply the replication model in our data - all kinds of wrappers before calling functions from *repl_model_fit.r*. Includes *repl.model.analysis*, which runs the analysis involving the replication until the classification of additional cells into C2.1 and C2.3 (figure 3).   


*emb_proper_analysis.r* - code for analyzing the embryo-proper clusters, starting with the classification of additional cells into C2.1 and C2.3 (more of figure 3).   

*epigenetic_analysis.r* - code for the analysis of the C2.1 and C2.3 clusters together with histone modifications data (figure 4).

Next are some utility files:  
*pooled_hic_plots.r*, *hic_plots.r* - functions to plot Hi-C matrices.  
*utils.r* - different util functions.  
*track_creation.r* - code to create misha tracks for pools of cells used in the analysis. These functions are not called by the various scripts, so in order to recreate all results you should call these functions manually to create these misha tracks if needed.  
*params.r* - configuration file. See instructions above about the variables that should be edited in this file.

*atac* - this directory contains the code for the ATAC analysis (figure 5). *pipe.r* contains the main script.  


