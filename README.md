# GX-Affymetrix
GX analysis for Affymetrix plate array

R library required: `oligo`, `limma`, `knitr`, `genefilter`, `ggplot2`, `reshape`, `dendextend`, `colorspace`, `scatterplot3d`
Additionally: Annotation package for you array (example: `pd.hugene.2.1.st`, `pd.mogene.2.1.st` for human and mouse)

The pipeline is run from a wrapper python script, preferably from the command line on a UNIX system.
The wrapper scripts takes one arguments, --arguments_file, which is the arguments file and it should contains the following arguments: (with the format `argument = `)
* **sample file** : containing two columns, `Sample.Name` (sample name used for the cels files) and `Sample.Group`  (group the sample belong to)
* **qcc file** :Library qcc file if not present for your array (currently human and mouse are provide in this hub, in the QCC folder), available in the library files from Affymetrix. Required to plot the controls.
* **comparison file** : specifying which groups to compare (pairwise groups), two columns, one labelled `Group1` and the other `Group2`, filled in with groups to compare, name must match the name in the `Sample.Group` column of the sample info file
* **raw data folder** : folder containing the raw data file (cels)

To run the pipeline:
```
git clone https://github.com/CGSbioinfo/GX-Affymentrix.git
```
Then to start the analysis once all the arguments are set:
```
python run_affymetrix_analysis.py --arguments_file arguments.txt
```
**Note:** that currently you will have to convert the tex file generated to a pdf file manually once the pipeline has finished running.
