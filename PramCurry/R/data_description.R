#==============================================================================#
#' \emph{Phytophthora ramorum} forest genotype data.
#' 
#' @name ramdat
#' @docType data
#' @usage data(ramdat)
#' @description This data set contains 513 diploid samples from the forest in 
#'   Curry County, OR, USA between 2001 and 2014. A population hierarchy
#'   contains the names Pop, ZONE1, and ZONE2. Pop represents years, ZONE2
#'   represents broad-scale watershed regions, and ZONE1 represents finer-scale
#'   watershed regions. The \code{other} slot contains xy coordinates for each
#'   sample ($xy), the population hierarchy in a data frame
#'   ($population_hierarchy), and the repeat motif lengths for each of the 5
#'   loci ($REPLEN).
#' @format a \code{\link{genclone}} object.
#' @references ZN Kamvar, MM Larsen, AM Kanaskie, EM Hansen  and NJ Grünwald. 
#'   20XX. Spatial and temporal population dynamics of the sudden oak death 
#'   epidemic in Oregon Forests. Phytopathology XX:XXX-XXX.
#==============================================================================#
NULL

#==============================================================================#
#' \emph{Phytophthora ramorum} forest and nursery genotype data.
#' 
#' @name for2nur
#' @docType data
#' @usage data(for2nur)
#' @description This data set contains 729 diploid samples from the forest in 
#'   Curry County, OR, USA between 2001 and 2014 and genotypes found in
#'   Nurseries from CA and OR between 2000 and 2012. A population hierarchy
#'   contains the names SOURCE, YEAR, and STATE. SOURCE represents either
#'   watershed regions (forest data) or nursery. YEAR represents the year
#'   sampled. STATE represents the state of origin (OR or CA). The \code{other}
#'   slot contains xy coordinates for each sample ($xy), the population
#'   hierarchy in a data frame ($population_hierarchy), and the repeat motif
#'   lengths for each of the 5 loci ($REPLEN).
#' @format a \code{\link{genclone}} object.
#' @references ZN Kamvar, MM Larsen, AM Kanaskie, EM Hansen  and NJ Grünwald. 
#'   20XX. Spatial and temporal population dynamics of the sudden oak death 
#'   epidemic in Oregon Forests. Phytopathology XX:XXX-XXX.
#==============================================================================#
NULL

#==============================================================================#
#' Metadata for the forest populations of \emph{Phytophthora ramorum}.
#' 
#' @name pop_data
#' @docType data
#' @usage data(pop_data)
#' @description This data frame contains 9 columns containing metadata that 
#'   pretains to the data contained in \code{\link{ramdat}}: \itemize{ 
#'   \item{ISOLATE - }{Original names given to each isolate. Decimal notations 
#'   represent isolates collected from same host.} \item{SAMPLE - }{Indicator of
#'   host sample each isolate originated from} \item{LAT - }{Lattitude 
#'   measurements for each isolate} \item{LON - }{Longitude measurements for 
#'   each isolate} \item{NOTES - }{Notes for specific isolates} \item{ZONE1 - 
#'   }{Watershed classifications at a fine scale} \item{ZONE2 - }{Watershed 
#'   classifications at a broad scale} \item{YEAR - }{Year isolated} \item{MLG -
#'   }{Multilocus genotype assignment as defined in poppr}}
#'   
#' @format a \code{\link{data.frame}} with 9 columns.
#' @references ZN Kamvar, MM Larsen, AM Kanaskie, EM Hansen  and NJ Grünwald. 
#'   20XX. Spatial and temporal population dynamics of the sudden oak death 
#'   epidemic in Oregon Forests. Phytopathology XX:XXX-XXX.
#==============================================================================#
NULL

#==============================================================================#
#' Color palette for 70 Multilocus Genotypes in \emph{Phytophthora ramorum} 
#' forest data.
#' 
#' @name myPal
#' @docType data
#' @usage data(myPal)
#' @description This palette defines a separate color for each MLG in the data 
#'   set \code{\link{ramdat}}. This is especially useful for the function 
#'   \code{\link[ggplot2]{scale_color_manual}} in the \pkg{ggplot2} package.
#' @format a named character vector of hexadecimal colors derived from the
#'   function \code{\link{char2pal}} with the \pkg{adegenet} palette
#'   \code{\link[adegenet]{funky}}.
#' @references ZN Kamvar, MM Larsen, AM Kanaskie, EM Hansen  and NJ Grünwald. 
#'   20XX. Spatial and temporal population dynamics of the sudden oak death 
#'   epidemic in Oregon Forests. Phytopathology XX:XXX-XXX.
#==============================================================================#
NULL

#==============================================================================#
#' Color palette for 9 populations defined in forest and nursery isolates of 
#' \emph{Phytophthora ramorum}.
#' 
#' @name comparePal
#' @docType data
#' @usage data(comparePal)
#' @description This palette defines a separate color for each popualtion in the
#'   data set \code{\link{for2nur}}. This is based on the hierarchy of 
#'   \code{~SOURCE/STATE}. Forest populations are defined by colors, whereas 
#'   Nursery populations are defined by shades of grey/black. This is especially
#'   useful for the function \code{\link[ggplot2]{scale_color_manual}} in the 
#'   \pkg{ggplot2} package.
#' @format a named character vector of hexadecimal colors derived from the
#'   function \code{\link{char2pal}} with the \pkg{adegenet} palette
#'   \code{\link[adegenet]{funky}}.
#' @references ZN Kamvar, MM Larsen, AM Kanaskie, EM Hansen  and NJ Grünwald. 
#'   20XX. Spatial and temporal population dynamics of the sudden oak death 
#'   epidemic in Oregon Forests. Phytopathology XX:XXX-XXX.
#==============================================================================#
NULL
