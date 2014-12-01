# Analysis of *P. ramorum* epidemic isolates

Citation
========
Please use the following DOI to cite this data and associated functions/analyses.

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.13007.svg)](http://dx.doi.org/10.5281/zenodo.13007)


R markdown files
========

These files contain the steps used to perform the analyses with comments.
Reports are written in `\*.Rmd` files and rendered as HTML documents. Folders
contain cached figures and data from respective analyses. Below is a table
describing each analysis.

|File|Description|
|----|----|
|**DAPC**|Discriminant analysis of principal components over forest data.|
|**Diversity_stats**|file with genotypic, genotype, and allelic diversity. This includes AMOVA analysis|
|**ForestNursery**|Analyzing the combined forest and nursery data.|
|**Gst**|Calculating $G_{st}$ values for the data set|
|**Gst_clonecorrect**|Calculating $G_{st}$ values for clone corrected data.|
|**Mantel**|Analyzing correlations between space and genotype with mantel tests and linear models.|
|**maps**| plotting multilocus genotypes onto maps. Note: the output is not kept in the online repository as the file size with all the images becomes too large to host freely on github.|
|**mlg_distribution**| Displaying the distribution of multilocus genotypes in the data. |
|**msn**| Displays minimum spanning networks based off of Bruvo's distance. |
|**Pop_structure**|dendrograms of populations based on Nei's distance.|

R package
=========

### PramCurry

This package contains data and functions utilized in this analysis.

installation can be done from within the main directory using devtools:

```r
devtools::install("./PramCurry")
```

Data Sets
=========

  - **ramdat** - a data set of all the unique forest isolates genotyped over
    five loci from Curry County, OR from 2001 to 2014. The csv file containing
    this data is in the main directory under
    `Pramorum_Curry_County_2001_2014.csv` in GenAlEx format with xy coordinates.
    This can be accessed from within the **PramCurry** package with the command
    `data(ramdat)`
  - **for2nur** - combination of forest and nursery isolates genotyped over five
    loci. The csv file containing this data is in the main directory under
    `Pramorum_Forest_Nursery.csv` in GenAlEx format. This can be accessed from
    within the **PramCurry** package with the command `data(for2nur)`
  
Misc
=========

### Shapefiles

A folder called `shapefiles/` contains all of the shapefiles necessary to
produce figure 1.

### Packrat

This project is managed via [packrat](http://rstudio.github.io/packrat/). There
is a folder called `packrat/` that contains all the information about the R
packages utilized to perform this analysis. It also contains a folder called
`src/`, which holds all of the tarballs for each package.

#### Exceptions

The R package [rgdal](http://cran.r-project.org/web/packages/rgdal/) is a
library that will read and translate shapefiles into a format that R can
understand. It requires that you download and install
[GDAL](http://www.gdal.org/) onto your system. I have left it out because OSX
\>10.9 does not have a readily available binary for it. For users on any other
system, you should install GDAL and then install rgdal normally:
`install.packages("rgdal")`. For users of OSX \>10.9, please use
`install.packages("rgdal_0.8-16.tgz", repos = NULL)` with the file in this
repository.

### Various Files

Many of the analyses will produce either image files (\*.gif, \*.pdf, \*.png,
\*.svg), data files (\*.nex), or tables (\*.csv). These were saved in the main
folder of the repository.
