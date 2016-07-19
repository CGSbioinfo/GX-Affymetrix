# GX-Affymetrix
GX analysis for Affymetrix plate array

Requirements:
* R/Rstudio with the following packages: oligo, limma, knitr, genefilter, ggplot2, reshape, dendextend, colorspace, scatterplot3d
* Annotation package for you array (example: pd.hugene.2.1.st, pd.mogene.2.1.st for human and mouse)
* For report generation CGS_logo_final.jpg and UNI_logo.jpg
* Library qcc file if not present for your array (currently human and mouse are provide in this hub), available in the library files from Affymetrix. Required to plot the controls.
* Sample info file, sample_info.txt, containing two columns, `Sample.Name` (sample name used for the cels files) and `Sample.Group`  (group the sample belong to)

To run the pipeline, copy all the file in the analysis folder, the cel file should be copied in to a folder called "cels". 
Then open R/RStudio and set the working directory to the directory where all the files are. Type `source("run_pipeline.R")` and the analysis will start.
To change the parameters you need to edit the script called `affy_pipeline_v1.R`, such as the array type to switch from human to mouse for example.
Adding other types of arrays is possible but requires more in depth modification of the scripts.
