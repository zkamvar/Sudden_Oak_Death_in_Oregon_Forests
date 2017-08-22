#' @import ggplot2 ggmap poppr ape reshape2 dplyr knitr adegenet
NULL
#==============================================================================#
# Stuff up here.
#==============================================================================#
#' Function to return numeric MLGs from mlg strings.
#' 
#' Yup
#'
#' @param x a string of MLG definitions formatted as "MLG.n"
#' where n is a number.
#' @return a numeric vector defining the multilocus genotypes.
#' @author Zhian N. Kamvar
#' @export
#' @examples
#' mlgFromString("MLG.100")
#==============================================================================#
mlgFromString <- function(x) as.numeric(substr(x, 5, nchar(x)))

#' Get bounding box from Latitute and Lognitude vectors.
#' 
#' @param LAT a numeric vector of latitudes
#' @param LON a numeric vector of longitudes
#' @return a numeric vector of four elements containing:
#' \enumerate{
#' \item lower left LON
#' \item lower left LAT
#' \item upper right LON
#' \item upper right LAT}
#' @export
get_ranges <- function(LAT, LON){
  c(min(LON, na.rm = TRUE), 
    min(LAT, na.rm = TRUE), 
    max(LON, na.rm = TRUE), 
    max(LAT, na.rm = TRUE))
}


#' Rearrange mlg.table output into a data frame for plotting.
#' @param mat a matrix derived from \code{\link[poppr]{mlg.table}}
#' @param total defaults to \code{TRUE} if there is a row for total MLG counts
#' @param type "heatmap" (default) will set all zero values to NA. "identity"
#' will keep zeroes as zero.
#' @return a molten data frame with the values of "MLG", "count", and "unique".
#' @keywords internal
prepare_mlg_df <- function(mat, total = TRUE, type = c("heatmap", "identity")){
  mat.df <- melt(mat, value.name = "count", varnames = c("Year", "MLG"))
  ARGS <- c("heatmap", "identity")
  type <- match.arg(type, ARGS)
  # Set the uncounted MLGs to NA so that it will be reflected in the heatmap.
  if (type == "heatmap"){
    mat.df$count[mat.df$count == 0] <- NA
  }
  if (total){
    unique.mlg <- colSums(mat[-nrow(mat), ] > 0) == 1
  } else {
    unique.mlg <- colSums(mat > 0) == 1
    mat.df$Year <- factor(mat.df$Year, levels = unique(mat.df$Year))
  }
  
  unique.mlg.df     <- data.frame(unique = unique.mlg)
  unique.mlg.df$MLG <- rownames(unique.mlg.df)
  
  # unique.pops <- rowSums(mat[, unique.mlg]) > 0
  # unique.mlg.pops <- mat[unique.pops, unique.mlg]
  
  mat.df <- merge(mat.df, unique.mlg.df)
  return(mat.df)
}


#' Plot MLG occurance by population (year)
#' @param mat a matrix derived from \code{\link[poppr]{mlg.table}}
#' @param total defaults to \code{TRUE} if there is a row for total MLG counts
#' @param log_transform should the counts be log transformed for easier viewing? Default to TRUE.
#' @return a ggplot object that contains a heatmap of multilocus genotypes
#' per year.
#' @export 
mlg.heatmap <- function(mat, total = TRUE, log_transform = TRUE){
  mat.df <- prepare_mlg_df(mat, total = total, type = "heatmap")
  plotlabs <- round(exp(0:round(log(max(mat.df$count, na.rm = TRUE)))))
  if (log_transform){
    mat.df$count <- log(mat.df$count)
  }
  myplot <- ggplot(mat.df, aes_string(x = "Year", y = "MLG", fill = "count")) + 
    geom_tile() + 
    geom_line(aes_string(color = "unique", group = "MLG", alpha = "unique"),
                         size = rel(5)) +
    scale_fill_continuous(low = "blue", high = "yellow", na.value = "grey95",
                          labels = plotlabs) +
    scale_color_manual(values = c("white", "black")) +
    scale_alpha_manual(values = c(0, 0.25)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  return(myplot)
}

#' Plot MLG distribution by year.
#' @param mat a matrix derived from \code{\link[poppr]{mlg.table}}
#' @param total defaults to \code{TRUE} if there is a row for total MLG counts
#' @param pal a color palette or vector of named colors to be used for MLGs.
#' @return a ggplot object that contains a barplot of MLG abundances per year.
#' @export 
mlg.barplot <- function(mat, total = TRUE, pal = funky){
  if (length(pal) == 1){
    PAL <- match.fun(pal)
    myPal <- PAL(ncol(mat))
    names(myPal) <- colnames(mat)
  } else {
    myPal <- pal
  }
  mat.df <- prepare_mlg_df(mat, total = total, type = "identity")
  myPlot <- ggplot(mat.df, aes_string(x = "Year", fill = "MLG", y = "count")) + 
    geom_bar(position = "fill", stat = "identity") +
    #   geom_text(aes(label = rowSums(mat)), mapping = aes_string(y = 1)) + 
    scale_fill_manual(values = myPal) + 
    guides(fill = guide_legend(ncol = 3)) +
    annotate("text", x = 1:nrow(mat), y = 1, hjust = 1,
             label = paste("n", rowSums(mat), sep = "="), angle = 90) +
    theme_classic() + 
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  return(myPlot)
}

#' Function to query ranges of MLG observations.
#'
#' In terms of what I use it for: To find the first and last date a MLG was found.
#' @param x a logical vector of a single MLG observations per year.
#' @return two integers indicating the indices for the first and last observed
#' instances of the MLG.
#' @keywords internal
begin_end <- function(x) range(which(x))

#' Transform MLG matrix into date range of observations per MLG
#' @param mat a TRUE/FALSE matrix of MLGs (columns) per year (rows)
#' @return a character matrix of two columns and rows equal to the
#' number of observed MLGs giving the first and last observed instances.
#' @export
MLG_range <- function(mat){
  res <- matrix(ncol = 2, nrow = ncol(mat), 
                dimnames = list(MLG  = colnames(mat),
                                seen = c("first", "last"))
  )
  first_last     <- t(apply(mat, MARGIN = 2, FUN = begin_end))
  res[, "first"] <- rownames(mat)[first_last[, 1]]
  res[, "last"]  <- rownames(mat)[first_last[, 2]]
  return(res)
}

#' Add leading zeroes onto a number range.
#' @param x a numeric vector 
#' @return a character vector where all numbers have the same
#' number of characters. 
#' @export
pad_zeroes <- function(x){
  nzeroes <- paste0("%0", max(nchar(x)), "d")
  return(sprintf(nzeroes, x))
}

#' Generate layout function
#' 
#' @param graph an igraph object
#' @param LAYOUT a layout function defined by igraph (default) OR a matrix
#' defining the layout. 
#' @return a function that will return the a layout matrix based on the layout 
#' function.
#' 
#' @details This function is used for generating a function that will return a
#' pre-defined layout. This is useful for subsetting networks but wanting all of
#' the points to be equivalent. Used with the layfun argument of plot_poppr_msn
#' @export
#' @importFrom igraph layout.auto V V<- incident
get_layout <- function(graph, LAYOUT = layout.auto){
  if (is.function(LAYOUT)){
    LAYFUN <- match.fun(LAYOUT)
    lay <- LAYFUN(graph)
  } else {
    lay <- LAYOUT[rownames(LAYOUT) %in% V(graph)$label, ]
    lay <- lay[V(graph)$label, ]
  }
  rownames(lay) <- V(graph)$label
  return(function(x) lay)
}

#' Creates a named color palette.
#' 
#' This is useful for defining a color palette that can be used by population factors.
#' 
#' @param x a vector of identifiers to be used for colors.
#' @param pal a color palette. Default is adegenet's "funky"
#' 
#' @return a named character vector of hexadecimal colors.
#' 
#' @export
#' @examples
#' char2pal(LETTERS)
char2pal <- function(x, pal = adegenet::funky){
  PAL <- match.fun(pal)
  outPal <- PAL(length(unique(x)))
  names(outPal) <- unique(x)
  return(outPal)
}

#' Removes edges not shared in both graphs
#' 
#' This function will take two graphs, \code{a} and \code{b}, where \code{b} is
#' a subgraph of \code{a} and remove all edges in \code{a} that are not present
#' in \code{b}. Note that there can be edges in \code{b} that are not present in
#' \code{a}. This function will only return edges present in \code{a}.
#' 
#' @param a an igraph graph object
#' @param b an igraph graph object with a subset of nodes from a OR a character
#'   labeling the nodes to subset.
#' @return an modified version of graph \code{a} only containing edges from
#'   graph \code{b}
#' @keywords internal
#' @importFrom igraph is.igraph delete.edges
remove_the_edges <- function(a, b){
  V(a)$name <- V(a)$label
  if (is.igraph(b)){
    V(b)$name <- V(b)$label
    vnion <- V(a)$name %in% V(b)$name
  } else {
    vnion <- V(a)$name %in% b
  }
  temp <- a
  for (i in which(!vnion)){
    temp <- delete.edges(temp, incident(temp, i))
  }
  #   zero_degree <- which(degree(temp)[vnion] == 0)
  #   for (i in zero_degree){
  #     theEdges <- incident(b, i)
  #     v_index  <- get.edges(b, theEdges)
  #     theNames <- V(b)$name[v_index]
  #     temp     <- add.edges(temp, theNames, 
  #                       attr = list(weight = E(b)[theEdges]$weight))
  #   }
  return(temp)
}

#' Highlight a subpopulation in a minimum spaning network
#' 
#' @param x a genind or genclone object.
#' @param poppr_msn a minimum spanning network derived from \code{x}.
#' @param ... arguments to be passed on to \code{plot_poppr_msn}
#' @param sublist a character string indicating the subpopulation to plot.
#' @param marker should the samples within the subgraph be marked with a polygon
#'   underneath the nodes?
#' 
#' @return NULL
#' 
#' @export
#' @importFrom igraph ecount 
population_subgraph <- function(x, poppr_msn, ..., sublist = "ALL", marker = TRUE){
  dots <- list(...)
#   ppmforms    <- formals(plot_poppr_msn)
#   poppr_args  <- names(dots) %in% names(ppmforms)
#   igraph_args <- !poppr_args
  if (sublist == "ALL"){
    plot_poppr_msn(x, poppr_msn, ...)
    invisible(return(NULL))
  }
  all_colors  <- poppr_msn$colors
  names(all_colors) <- poppr_msn$populations
  grey_inds <- !names(all_colors) %in% sublist
  if (any(names(dots) == "palette")){
    palmarker <- TRUE
    PAL <- match.fun(dots$palette)
    the_colors <- PAL(length(all_colors))
  } else {
    palmarker <- FALSE
    the_colors <- all_colors
  }
  the_colors[grey_inds] <- grey.colors(length(all_colors), 0.85, 1)[grey_inds]

  newPal <- function(x) the_colors
  if (palmarker){
    dots$palette <- newPal
  } else {
    poppr_msn <- update_my_graph(poppr_msn, newPal)
  }
  mlgs <- popsub(x, sublist)@mlg
  mlgsub <- paste("MLG", unique(mlgs), sep = ".")


  theGraph <- remove_the_edges(poppr_msn$graph, mlgsub)
  if (ecount(theGraph) > 0){
    poppr_msn$graph <- theGraph
  }

  arglist <- list()
  if (!any(names(dots) == "inds")){
    arglist$inds <- unique(mlgs)  
  }
  if (marker){
    group_color <- the_colors[!grey_inds]
    group_color <- adjustcolor(group_color, 0.3)
    arglist$mark.col <- group_color
    arglist$mark.border <- group_color
    markV <- list(which(V(poppr_msn$graph)$label %in% mlgsub))
    arglist$mark.groups <- markV
#     print(dput(markV))
  }
  do.call(plot_poppr_msn, c(list(x, poppr_msn), 
                          dots, 
                          arglist
                          )
        )
  invisible(return(NULL))
}


#==============================================================================#
# Function used to update colors in poppr msn 
#
# Public functions utilizing this function:
# plot_poppr_msn
#
# Private functions utilizing this function:
# # none
#==============================================================================#
#' Function used to update colors in poppr msn 
#' 
#' @param graphlist a graph produced from one of poppr's *.msn functions
#' @param PALETTE a function or character defining a palette
#' 
#' @return a shiny new graph
#' 
#' @export
update_my_graph <- function(graphlist, PALETTE){
  PALETTE <- match.fun(PALETTE)
  lookup  <- data.frame(old    = graphlist$populations, 
                        update = PALETTE(length(graphlist$populations)), 
                        stringsAsFactors = FALSE)
  if (nrow(lookup) > 1){
    colorlist                    <- V(graphlist$graph)$pie.color
    V(graphlist$graph)$pie.color <- lapply(colorlist, update_colors, lookup)
  } else {
    colorlist <- V(graphlist$graph)$color
    V(graphlist$graph)$color <- rep(PALETTE(1), length(colorlist))
  }
  graphlist$colors <- lookup[[2]]
  names(graphlist$colors) <- graphlist$populations
  return(graphlist)
}

#==============================================================================#
# Function used to update colors in poppr msn 
#
# Public functions utilizing this function:
# none
#
# Private functions utilizing this function:
# # update_poppr_graph
#==============================================================================#
update_colors <- function(colorvec, lookup){
#   x <- vapply(1:length(colorvec), update_single_color, colorvec, lookup, colorvec)
#   return(x)
  pops     <- lookup[[1]]
  update   <- lookup[[2]]
  names(update) <- pops
  matching <- match(pops, names(colorvec))
  return(update[matching[!is.na(matching)]])
}
#==============================================================================#
# Function used to update colors in poppr msn 
#
# Public functions utilizing this function:
# none
#
# Private functions utilizing this function:
# # update_colors
#==============================================================================#
update_single_color <- function(x, lookup, colorvec){
  update   <- lookup[[2]]
  original <- lookup[[1]]
  return(update[original %in% colorvec[x]])
}


#' match and return all duplicated elements.
#' 
#' @param x a character or numeric vector
#' @return a logical vector of elements that are duplicated.
#' 
#' @export
match_duplicates <- function(x){
  if (!anyDuplicated(x)){
    return(duplicated(x))
  } else {
    dupes <- duplicated(x)
    return(x %in% x[dupes])
  }
}


#' Merge or remove duplicate entries from a genclone object
#' 
#' This will compare duplicate names in genclone objects and merge the ones with the
#' same genotype and remove the ones with differing genotypes
#' 
#' @param x a genind object
#' @param stats if set to \code{TRUE}, the status of the merges will be printed
#' to the console.
#' @return a genind object
#' 
#' @export
merge_duplicates <- function(x, stats = TRUE){
  stopifnot(suppressWarnings(is.genind(x)))
  dupes        <- match_duplicates(indNames(x))
  tab          <- table(indNames(x)[dupes], x@mlg[dupes]) > 0
  toRemove     <- rownames(tab)[rowSums(tab) > 1]
  noUnmatched  <- !indNames(x) %in% toRemove
  toCombine    <- !duplicated(indNames(x))
  merged       <- noUnmatched & toCombine
  if (stats){
    msg <- paste("Found", nrow(tab), "duplicated genotypes:\t",
                 paste(rownames(tab), collapse = ", "), "\n",
                 length(toRemove), "did not have matching MLGs and",
                 "were removed:\t", paste(toRemove, collapse = ", "), "\n")
    message(msg)
  }
  return(x[merged, ])
}

#' Plot posterior values from DAPC analysis in adegenet
#' 
#' @param da.object an object of class "dapc"
#' @param gid an object of class "genind"
#' @param pal a color palette
#' @param cols the number of columns to display
#' @return a ggplot object with each population stacked on top of each other.
#' @export
plot_posterior <- function(da.object, gid, pal = funky, cols = 1){
  posterior <- da.object$posterior
  names(dimnames(posterior)) <- c("sample", "population")
  to_merge <- data.frame(list(sample = dimnames(posterior)$sample, 
                              oldPopulation = pop(gid)))
  post <- melt(posterior, value.name = "probability")
  post <- merge(post, to_merge)
  if (is.numeric(post$sample)){
    post$sample <- factor(post$sample, levels = unique(post$sample))
  }
  if (is.numeric(post$population)){
    post$population <- factor(post$population)
  }
  if (length(pal) == 1){
    PAL <- match.fun(pal)
    pal <- char2pal(post$population, PAL)
  }
  outPlot <- ggplot(post, aes_string(x = "sample", fill = "population", y = "probability")) + 
    geom_bar(stat = "identity", position = "fill", width = 1) + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    facet_wrap(~oldPopulation, scales = "free_x", drop = TRUE, ncol = cols) +
    scale_fill_manual(values = pal)
  return(outPlot)
}


#' Spatial test on multiple populations.
#' 
#' Performs a mantel test or linear model on multiple populations in a genclone
#' object.
#' 
#' @param gen a genclone object
#' @param xy a matrix of latitude and longitude data or a distance matrix.
#' @param stat a character selectine either "mantel" or "lm"
#' @param hierarchy a formula specifying the hierarchical levels to be used to
#'   define the population.
#' @param distance a function to caclulate the genetic distance or a genetic 
#'   distance matrix of class dist.
#' @param sample the number of replications for each mantel test (performed in
#'   C)
#' @param plot if TRUE, all possible
#' @param ncol if plot = TRUE, how many columns should be plotted?
#' @param seed numeric This is he seed to set before randomizations for the
#'   mantel test. By default, this is set to NULL so that the results are
#'   random. When given a numeric value, results are reproducible.
#' @param ... arguments passed onto the distance function.
#'   
#' @return a list containing multiple test results.
#'   
#' @export
spatial_stats <- function(gen, xy, stat = "mantel", hierarchy = NULL,
                          distance = "nei.dist", sample = 999, plot = TRUE, 
                          ncol = 3, seed = NULL, ...){

  if (length(dim(xy)) != 2 & !class(xy) %in% "dist" & !is.null(dim(xy))){
    stop("xy must be a two column matrix or data frame of geographic coordinates")
  } else if (!class(xy) %in% "dist" && dim(xy)[2] != 2){
    stop("xy must be a two column matrix or data frame of geographic coordinates")
  }
  if (!class(xy) %in% "dist"){
    dist_xy <- dist(xy)  
  } else {
    dist_xy <- xy
  }
  if (!class(distance) %in% "dist"){
    DISTFUN  <- match.fun(distance)
    dist_gen <- DISTFUN(gen, ...)
  } else {
    dist_gen <- distance
  }
  
  
  
  if (is.null(hierarchy)){
    if (stat == "mantel"){
      set.seed(seed)
      res <- mantel.randtest(dist_xy, dist_gen, nrepet = sample)
      if (plot){
        plot(res, main = paste("n", "=", nInd(gen)))
      }
    } else if (stat == "lm"){
      dist_xy_vec  <- as.vector(dist_xy)
      dist_gen_vec <- as.vector(dist_gen)
      res <- lm(dist_xy_vec ~ dist_gen_vec)
      if (plot){
        p <- ggplotRegression(res, "Total", "Genetic", "Geographic")
        print(p)
      }
    }
    
    return(res)
  }
  
  mat_gen     <- as.matrix(dist_gen)
  mat_xy      <- as.matrix(dist_xy)
  setpop(gen) <- hierarchy
  full_hier   <- gethierarchy(gen, hierarchy, combine = TRUE)
  full_hier   <- full_hier[[length(full_hier)]]
  pops        <- levels(full_hier)
  res         <- lapply(pops, do_space_test, mat_gen, mat_xy, full_hier, sample, 
                        stat, seed)
  names(res)  <- pops
  if (plot){
    if (stat == "mantel"){
      if (length(res) < ncol){
        par(mfrow = c(1, length(res)))
      } else {
        cols <- ncol
        rows <- ceiling(length(res)/ncol)
        par(mfrow = c(rows, cols))
      }
      lapply(names(res), plot_mantel, res, full_hier)
      par(mfrow = c(1, 1))
    } else if (stat == "lm"){
      plots <- lapply(names(res), function(x)
        ggplotRegression(res[[x]], x, "Genetic Distance", "Geographic Distance"))
      multiplot(plotlist = plots, cols = ncol)
    }
  }
  return(res)
}

#' @importFrom ade4 mantel.randtest
do_space_test <- function(population, mat_gen, mat_xy, pops, sample, 
                          stat = "mantel", seed){
  population   <- pops %in% population
  dist_gen_pop <- as.dist(mat_gen[population, population])
  dist_xy_pop  <- as.dist(mat_xy[population, population])
  if (stat == "mantel"){
    set.seed(seed)
    res <- mantel.randtest(dist_xy_pop, dist_gen_pop, nrepet = sample) 
  } else if (stat == "lm"){
    dist_xy_pop_vec  <- as.vector(dist_xy_pop)
    dist_gen_pop_vec <- as.vector(dist_gen_pop) 
    if (length(dist_xy_pop_vec) < 1){
      res <- lm(1 ~ 1)
    } else {
      res <- lm(dist_xy_pop_vec ~ dist_gen_pop_vec)
    }
  }
  return(res)
}

plot_mantel <- function(xname, xlist, n = ""){
  x <- xlist[[xname]]
  if (length(n) > 1){
    n <- paste0("\nn = ", sum(n %in% xname))
  }
  if (all(is.na(x$sim))){
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste0(xname, "\n", "insufficient data", n))
  } else {
    plot(x, main = paste0(xname, n))
  }
}

ggplotRegression <- function (fit, title, x, y) {
  if (is.na(fit$coef[2])){
    p <- qplot(1, 1, geom = "text", label = "insufficient data") + 
      theme(text = element_blank(), line = element_blank(), title = element_blank())
  } else {
    p <- ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
      geom_point(alpha = 0.5) +
      stat_smooth(method = "lm", col = "red") +
      ggtitle(paste0(title, "\n",
                    "Adj R2 = ", signif(summary(fit)$adj.r.squared, 5),
                    "; Intercept = ",signif(fit$coef[[1]],5 ),
                    "\n Slope = ",signif(fit$coef[[2]], 5),
                    "; P =",signif(summary(fit)$coef[2,4], 5))) +
      xlab(x) + ylab(y)
  }
  return(p)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
#' @importFrom grid pushViewport viewport grid.layout grid.newpage
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



#' Performs lazy load on a directory
#'
#' @param path a filepath containing the necessary files for lazy loading
#' @return NULL
#'
#' @details This function will go into a directory, search for all the files
#' that seem like they can be lazily loaded and attempt to load them.
#'
#' @export
#' @examples
#' \dontrun{
#' lazierLoad("Diversity_stats_cache/html")
#' }
lazierLoad <- function(path){
  files <- dir(path)
  cache_files <- sub(".rdb$", "", files[grepl(".rdb$", files)])
  for (i in cache_files){
    try(lazyLoad(paste(path, i, sep = "/"), envir = .GlobalEnv))
  } 
}

#' Test repeat length consistency.
#'
#' This function will test for consistency in the sense that all alleles are
#' able to be represented as discrete units after division and rounding.
#' @param gid a genind object
#' @param replen a numeric vector of repeat motif lengths.
#' @return a logical vector indicating whether or not the repeat motif length
#' is consistent.
#' @export
#' @examples
#' library(poppr)
#' data(nancycats)
#' test_replen(nancycats, rep(2, 9))
test_replen <- function(gid, replen){
  alleles <- lapply(gid@all.names, as.numeric)
  are_consistent <- vapply(1:nLoc(gid), consistent_replen, logical(1), 
                           alleles, replen)
  names(are_consistent) <- locNames(gid)
  return(are_consistent)
}

consistent_replen <- function(index, alleles, replen){
  !any(duplicated(round(alleles[[index]]/replen[index])))
}

#' Find and fix inconsistent repeat lengths
#'
#' Attempts to fix inconsistent repeat lengths found by \code{test_replen}
#'
#' @param gid a genind object
#' @param replen a numeric vector of repeat motif lengths.
#' @param e a number to be added to inconsistent repeat lengths to allow for 
#' proper rounding.
#' @return a numeric vector of corrected reapeat motif lengths.
#' @details Before being fed into the algorithm to calculate Bruvo's distance, 
#'   the amplicon length is divided by the repeat unit length. Because of the 
#'   amplified primer sequence attached to sequence repeat, this division does
#'   not always result in an integer and so the resulting numbers are rounded.
#'   The rounding also protects against slight mis-calls of alleles. Because we
#'   know that \deqn{\frac{(A + e) - (B + e)}{r}}{((A + e) - (B + e))/r} is
#'   equivalent to \deqn{\frac{A - B}{r}}{(A - B)/r}, we know that the primer 
#'   sequence will not alter the relationships between the alleles.
#'   Unfortunately for nucleotide repeats that have powers of 2, rounding in R
#'   is based off of the IEC 60559 standard (see \code{\link{round}}), that
#'   means that any number ending in 5 is rounded to the nearest \emph{even}
#'   digit. This function will attempt to alleviate this problem by adding a
#'   very small amount to the repeat length so that division will not result in
#'   a 0.5. If this fails, the same amount will be subtracted. If neither of
#'   these work, a warning will be issued and it is up to the user to determine
#'   if the fault is in the allele calls or the repeat lengths.
#' @export 
#' @examples
#' 
#' library(poppr)
#' data(nancycats)
#' fix_replen(nancycats, rep(2, 9))
#' # Let's start with an example of a tetranucleotide repeat motif and imagine
#' # that there are twenty alleles all 1 step apart:
#' (x <- 1:20L * 4L)
#' # These are the true lengths of the different alleles. Now, let's add the
#' # primer sequence to them. 
#' (PxP <- x + 21 + 21)
#' # Now we make sure that x / 4 is equal to 1:20, which we know each have
#' # 1 difference.
#' x/4
#' # Now, we divide the sequenc with the primers by 4 and see what happens.
#' (PxPc <- PxP/4)
#' (PxPcr <- round(PxPc))
#' diff(PxPcr) # we expect all 1s
#' 
#' # Let's try that again by adding a tiny amount to 4
#' (PxPc <- PxP/4.00001)
#' (PxPcr <- round(PxPc))
#' diff(PxPcr)
fix_replen <- function(gid, replen, e = 1e-5){
  if (length(replen) != nLoc(gid)) {
    stop(paste0("length of repeats (", length(replen), ") does not equal",
                " the number of loci (", nLoc(gid), ")."))
  }
  consistent_reps <- test_replen(gid, replen)
  names(replen)   <- locNames(gid)
  ADD <- FALSE
  SUB <- FALSE
  newReps <- replen
  while (any(!consistent_reps)){
    if (!ADD){
      newReps[!consistent_reps] <- newReps[!consistent_reps] + e
      ADD <- TRUE
    } else {
      newReps[!consistent_reps] <- newReps[!consistent_reps] - (2*e)
      SUB <- TRUE
    }
    consistent_reps <- test_replen(gid, newReps)
    if (any(!consistent_reps) & ADD & SUB){
      inconsistent <- paste(names(replen[!consistent_reps]), collapse = ", ")
      msg <- paste("The repeat lengths for", inconsistent, 
                   "are not consistent.\n\n",
                   "This might be due to inconsistent allele calls or repeat",
                   "lengths that are too large.\n",
                   "Check the alleles to make sure there are no duplicated",
                   "or similar alleles that might end up being the same after",
                   "division.\n\nOriginal repeat lengths are being returned.")
      warning(msg, immediate. = TRUE)
      consistent_reps <- TRUE
      newReps <- replen
    }
  }
  return(newReps)
}

#==============================================================================#
#' Function to reset MLG counts.
#' 
#' This is a workaround for attempting to subset populations before running
#' AMOVA on genclone objects.
#'
#' @param x genclone object.
#' @return a genclone object with reset MLG counts. 
#' @author Zhian N. Kamvar
#' @export
#' @examples
#' data(Pinf)
#' resetMLG(Pinf)
#==============================================================================#
resetMLG <- function(x){
  xhier  <- gethierarchy(x)
  xother <- other(x)
  gid    <- genind(truenames(x)$tab, pop = pop(x), ploidy = x@ploidy, type = x@type)
  res    <- as.genclone(gid, hierarchy = xhier)
  other(res) <- xother
  return(res)
}

#' Produce a table of diversity statistics
#' 
#' @param z a table of integers representing counts of MLGs (columns) per 
#' population (rows)
#' 
#' @return a numeric matrix with 4 columns giving the following statistics for 
#'   each population: \itemize{\item H - Shannon Diversity \item G Stoddart and
#'   Taylor's Diveristy (inverse Simpson's) \item unbiased Simpson Diversity \eqn{(N/(N-1))*(1 - D)} \item E.5 evenness}
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' library(poppr)
#' data(Pinf)
#' tab <- mlg.table(Pinf, bar = FALSE)
#' get_stats(tab)
get_stats <- function(z){
  mat <- matrix(nrow = nrow(z), ncol = 4, 
                dimnames = list(Pop = rownames(z), 
                                Index = c("H", "G", "Hexp", "E.5")
                )
  )
  N     <- rowSums(z)
  H     <- vegan::diversity(z)
  G     <- vegan::diversity(z, "inv")
  Simp  <- vegan::diversity(z, "simp")
  nei   <- (N/(N-1)) * Simp
  E.5   <- (G - 1)/(exp(H) -1)
  mat[] <- c(H, G, nei, E.5)
  return(mat)
}

#' @importFrom vegan diversity
boot_stats <- function(x, i){
  res        <- numeric(4)
  names(res) <- c("H", "G", "Hexp", "E.5")
  z     <- table(x[i])
  N     <- sum(z)
  H     <- vegan::diversity(z)
  G     <- vegan::diversity(z, "inv")
  Simp  <- vegan::diversity(z, "simp")
  nei   <- (N/(N-1)) * Simp
  E.5   <- (G - 1)/(exp(H) -1)
  res[] <- c(H, G, nei, E.5)
  return(res)
}

extract_samples <- function(x) rep(1:length(x), x)

#' Perform a bootstrap analysis on diversity statistics
#' 
#' @param tab a table produced from the \pkg{poppr} function \code{\link[poppr]{mlg.table}}. MLGs in columns and populations in rows
#' @param n an integer > 0 specifying the number of bootstrap replicates to perform (corresponds to \code{R} in the function \code{\link[boot]{boot}}.
#' @param ... other parameters passed on to \code{\link[boot]{boot}}.
#' 
#' @return a list of objects of class "boot". 
#' @seealso \code{\link{boot_ci}}
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' library(poppr)
#' data(Pinf)
#' tab <- mlg.table(Pinf, bar = FALSE)
#' do_boot(tab, 10L)
#' \dontrun{
#' # This can be done in a parallel fasion (OSX uses "multicore", Windows uses "snow")
#' system.time(do_boot(tab, 10000L, parallel = "multicore", ncpus = 4L))
#' system.time(do_boot(tab, 10000L))
#' }
#' @importFrom boot boot
do_boot <- function(tab, n, ...){
  res <- apply(tab, 1, function(x) boot::boot(extract_samples(x), boot_stats, n, ...))
  return(res)
}

get_ci <- function(x, lb, ub){
  res <- apply(x$t, 2, quantile, c(lb, ub), na.rm = TRUE)
  return(res)
}

get_all_ci <- function(res, ci = 95){
  lower_bound  <- (100 - ci)/200
  upper_bound  <- 1 - lower_bound
  funval       <- matrix(numeric(8), nrow = 2)
  CI           <- vapply(res, FUN = get_ci, FUN.VALUE = funval, 
                         lower_bound, upper_bound)
  dCI          <- dimnames(CI)
  dimnames(CI) <- list(CI    = dCI[[1]], 
                       Index = c("H", "G", "Hexp", "E.5"),
                       Pop   = dCI[[3]])
  
  return(CI)
}

#' Perform bootstrap statistics, calculate and plot confidence intervals.
#' 
#' @param tab a genind object OR a matrix produced from \code{\link[poppr]{mlg.table}}.
#' @param n an integer defining the number of bootstrap replicates (defaults to 1000).
#' @param ci the percent for confidence interval.
#' @param total argument to be passed on to \code{\link[poppr]{mlg.table}} if \code{tab} is a genind object.
#' @param ... parameters to be passed on to \code{\link[boot]{boot}}
#' 
#' @return an array of 3 dimensions giving the lower and upper bound, the index
#' measured, and the population. This also prints barplots for each population
#' faceted by each index using \pkg{ggplot2}. This plot can be retrieved by using
#' \code{p <- last_plot()}
#' 
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' library(poppr)
#' data(Pinf)
#' boot_ci(Pinf, n = 100)
#' \dontrun{
#' # This can be done in a parallel fasion (OSX uses "multicore", Windows uses "snow")
#' system.time(boot_ci(tab, 10000L, parallel = "multicore", ncpus = 4L))
#' system.time(boot_ci(tab, 10000L))
#' }
#' 
boot_ci <- function(tab, n = 1000, ci = 95, total = TRUE, ...){
  if (!is.matrix(tab) & is.genind(tab)){
    tab <- mlg.table(tab, total = total, bar = FALSE)
  }
  res  <- do_boot(tab, n, ...)
  orig <- get_stats(tab)
  orig <- melt(orig)
  orig$Pop <- factor(orig$Pop)
  CI   <- get_all_ci(res, ci = ci)
  samp <- vapply(res, "[[", FUN.VALUE = res[[1]]$t, "t")
  dimnames(samp) <- list(NULL, 
                         Index = c("H", "G", "Hexp", "E.5"),
                         Pop = rownames(tab))
  sampmelt <- melt(samp)
  sampmelt$Pop <- factor(sampmelt$Pop)
  pl <- ggplot(sampmelt, aes_string(x = "Pop", y = "value", group = "Pop")) + 
    geom_boxplot() + 
    geom_point(aes_string(color = "Pop", x = "Pop", y = "value"), 
               size = 5, pch = 16, data = orig) +
    xlab("Population") + labs(color = "Observed") +
    facet_wrap(~Index, scales = "free_y") + myTheme
  print(pl)
  return(CI)
}

myTheme <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


#' Calculate allelic or genotypic diversity by population. 
#' 
#' @param dat a genind or genclone object
#' @param hier a hierarchy to set for clone correction
#' @param combine logical. Should the hierarchy be combined to create the population after clone correction (default, FALSE, no).
#' @param lev see \code{\link[poppr]{locus_table}}.
#' @param cc logical. If TRUE (default), the data is clone corrected by the hierarchy
#' before stats are calculated.
#' 
#' @return a list containing \itemize{
#' \item{pops - }{a list of locus tables by population}
#' \item{total - }{a locus table for the pooled population}
#' }
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' library(poppr)
#' data(Pinf)
#' calc_loc_table(Pinf, hier = ~Continent/Country)
calc_loc_table <- function(dat, hier, combine = FALSE, lev = "allele", cc = TRUE){
  if (cc){
    dat <- clonecorrect(dat, hier = hier, combine = combine)
  }
  theTable <- locus_table(dat, information = FALSE, lev = lev)
  theHier <- gethierarchy(dat, hier, combine = combine)
  if (combine){
    setpop(dat) <- hier
    theHier <- levels(theHier[[length(theHier)]])
  } else {
    setpop(dat) <- as.formula(paste0("~", names(theHier)[1]))
    theHier <- levels(theHier[[1]])
  }
  loc_tab <- lapply(theHier, 
                    function(x){
                      locus_table(dat, 
                                  population = x, 
                                  information = FALSE,
                                  lev = lev)
                    })
  names(loc_tab) <- theHier
  return(list(pops = loc_tab, total = theTable))
}

#' Visualize locus tables
#' 
#' @param tab a list of locus tables produced from \code{\link{calc_loc_table}}.
#' @return a ggplot
#' @export
#' @author Zhian N. Kamvar
#' @examples
#' library(poppr)
#' data(Pinf)
#' plot_loc_table(calc_loc_table(Pinf, hier = ~Continent/Country))
plot_loc_table <- function(tab){
  #   tab <- theList[["pops"]]
  #   tot <- melt(theList["total"])
  tab <- melt(tab)
  ggplot(tab, aes_string(y = "value", x = "L2", fill = "L2", linetype = NA)) +
    geom_bar(stat = "identity") +
    geom_hline(aes_string(yintercept = "value", linetype = "L1"), 
               data = tab[tab$L1 == "total", ]) + 
    facet_grid(summary ~ locus, scales = "free_y") +
    scale_linetype_discrete(labels = "", name = "Pooled") +
    myTheme
}



amova_pair <- function(pops, dat, hier, ...){
  dat <- popsub(dat, pops)
  if (any(table(pop(dat)) < 3)){
    return(NULL)
  }
  poppr.amova(dat, hier, ...)
}

pairwise_amova <- function(x, hier, ...){
  pops         <- x@pop.names
  pop_combs    <- combn(pops, 2)
  xlist        <- apply(pop_combs, 2, amova_pair, x, hier, ...)
  names(xlist) <- apply(pop_combs, 2, paste, collapse = " : ")
  return(xlist)
}

pairwise_amova_test <- function(pairlist, nrepet = 99){
  res <- lapply(names(pairlist), print_amova_pairs, pairlist, nrepet)
  names(res) <- names(pairlist)
  return(res)
}

#' @importFrom ade4 randtest
print_amova_pairs <- function(pairname, pairlist, nrepet){
  thePair <- pairlist[[pairname]]
  if (!is.null(thePair)){
    cat(pairname, "\n")
    theTest <- randtest(thePair, nrepet = nrepet)
    try(plot(theTest))
    return(theTest)
  } else {
    return(NULL)
  }
}


getpval <- function(x){
  if (is.null(x$pvalue)){
    return(as.numeric(rep(NA, 4)))
  } else {
    pvals <- x$pvalue
  }
  pvals[x$rep == 0] <- NA
  return(pvals)
}

getphi <- function(x){
  if (is.null(x)) return(numeric(4))
  return(x$statphi$Phi)
}

refactor <- function(strings, factors){
  u_factors <- unique(factors)
  u_factors[match(strings, u_factors)]
}

make_amova_table <- function(am, amt){
  tot <- nrow(am$results)
  res <- data.frame(list(am$results[-tot, c("Df", "Sum Sq")], 
                         Percent = am$componentsofcovariance[-tot, 2],
                         Pval    = rev(amt$pvalue), 
                         Phi     = rev(am$statphi$Phi[-tot])))
  res <- as.matrix(res)
  colnames(res) <- c("d.f.", "Sum of Squares", "Percent variation", "P", 
                     "Phi statistic")
  names(dimnames(res)) <- c("levels", "statistic")
  return(res)
}

make_amova_printable <- function(amtab, amtabcc){
  am_array <- array(dim      = c(dim(amtab), 2),
                    dimnames = c(dimnames(amtab), 
                                 list(c("full", "clone-corrected"))))
  am_array[, , 1] <- amtab
  am_array[, , 2] <- amtabcc
  tabfun <- function(x){
    x <- paste0(paste0(signif(x, 3), collapse = " ("), ")")
    return(x)
  }
  res <- apply(am_array, c(1, 2), tabfun)
  return(res)
}

ci2string <- function(x, sig = 3, colString = "-"){
  if (all(is.na(x))) return("(-)")
  paste0("(", paste(signif(x, sig), collapse = colString),")")
}

add_ci <- function(dtable, ciarray, sig = 3, colString = "-"){
  
  dtable <- data.frame(lapply(dtable, function(x, sig){
    if (is.numeric(x)) x <- signif(x, sig)
    return(x)
  }, sig))
  for (i in colnames(ciarray)){
    to_add <- apply(ciarray[, i, ], 2, ci2string, sig, colString)
    dtable[[i]] <- paste(signif(dtable[[i]], sig), to_add)
  }
  return(dtable)
}
