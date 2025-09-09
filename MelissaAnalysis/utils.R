# All Functions Script 

# Section 1: 

# Sections 2: 

# Sections 3:

# Sections 4:

# Sections 5:

### Needed for all the diffbox's below ###
extract_args <- function(obj, arg){
  if (!"someparams" %in% methods::slotNames(obj)){
    stop("The object doesn't have a 'someparams' slot!")
  } else {
    res <- obj@someparams[[arg]]
    return(res)
  }
}

set_newlevels0 <- function(data, newlevels, factorNames){
  oldlevels <- unique(as.vector(data[[factorNames]]))
  newlevels <- oldlevels[match(newlevels, oldlevels)]
  newlevels <- newlevels[!is.na(newlevels)]
  data[[factorNames]] <- factor(data[[factorNames]], levels=newlevels)
  return(data)
}

set_newlevels <- function(data, newlevels, factorNames){
  oldlevels <- unique(as.vector(data[[factorNames]]))
  newlevels <- unique(newlevels)  # <-- ensure no duplicates
  newlevels <- oldlevels[match(newlevels, oldlevels)]
  newlevels <- newlevels[!is.na(newlevels)]
  data[[factorNames]] <- factor(data[[factorNames]], levels = newlevels)
  return(data)
}


#' @keywords internal
get_comparelist <- function(data, classgroup, controlgroup){
  groups <- get_classlevels(sampleda=data, classgroup=classgroup)
  if (!is.null(controlgroup)){
    groups <- setdiff(groups, controlgroup)
    tmplen <- length(groups)
    comparelist <- matrix(data=c(rep(controlgroup, tmplen), groups), nrow=tmplen)
  }else{
    comparelist <- get_compareclass(classlevels=groups)
  }
  comparelist <- split(comparelist, slice.index(comparelist, 1))
  names(comparelist) <- NULL
  return(comparelist)
}

keep_know_taxa <- function(featurelist, removeUnknown=TRUE){
  if (isTRUE(removeUnknown)){
    tmpflag <- grep("__un_",featurelist)
    if (length(tmpflag)>0){
      featurelist <- featurelist[-tmpflag]
    }
  }
  return(featurelist)
}

### Microbiome Process LDA ggdiff box left panel function ###

plotdiffbox <- function(obj, sampleda, factorNames, 
                        dodge_width = 0.6, box_width = 0.05, 
                        factorLevels = NULL, featurelist = NULL, 
                        box_notch = TRUE) {
  
  # Check inputs
  if (missing(sampleda) || is.null(sampleda)) {
    stop("The 'sampleda' should be provided!")
  }
  if (missing(factorNames) || is.null(factorNames)) {
    stop("The 'factorNames' should be provided!")
  }
  
  # Select relevant columns from sampleda
  sampleda <- sampleda[, match(factorNames, colnames(sampleda)), drop = FALSE]
  
  # Merge data and reshape
  featureda <- merge(obj, sampleda, by = 0) %>%
    column_to_rownames(var = "Row.names") %>%
    tidyr::pivot_longer(!factorNames, names_to = "feature", values_to = "value")
  
  # Set feature list and factor levels if provided
  if (!is.null(featurelist)) {
    featureda$feature <- factor(featureda$feature, levels = featurelist)
  }
  if (!is.null(factorLevels)) {
    featureda <- setfactorlevels(featureda, factorLevels)
  }
  
  # Aesthetic mapping for boxplot
  boxmapping <- aes_string(x = "feature", y = "value", color = factorNames)
  
  # Plot
  p <- ggplot(featureda, mapping = boxmapping) +
    geom_boxplot(notch = box_notch,
                 width = 10 * box_width,
                 outlier.color = NA,
                 position = position_dodge(width = dodge_width)) +
    theme_bw() +
    xlab(NULL) +
    theme(axis.text.x = element_text(size = 12),
          plot.margin = unit(c(2, 0, 2, 2), "mm"),
          panel.grid = element_blank())
  
  return(p)
}


### Microbiome Process full ggdiffbox funxtion ###

ggdiffbox <- function(obj, geom = "boxplot", 
                      box_notch = TRUE, box_width = 0.05, dodge_width = 0.6,
                      addLDA = TRUE, factorLevels = NULL, featurelist = NULL,
                      removeUnknown = TRUE, colorlist = NULL, l_xlabtext = NULL, ...) {
  
  featureda <- obj@originalD
  classname <- extract_args(obj, "classgroup")
  normalization <- extract_args(obj, "normalization")
  
  if (!is.null(normalization)) {
    featureda <- featureda / normalization
  }
  
  sampleda <- obj@sampleda
  tmpgroup <- unique(as.vector(sampleda[[classname]]))
  
  nodedfres <- obj@result
  nodedfres <- set_newlevels(data = nodedfres, newlevels = tmpgroup, factorNames = classname)
  
  if (is.null(featurelist)) {
    featurelist <- unique(as.vector(nodedfres$f))
  }
  
  params <- list(...)
  if (!is.null(params$removeUnkown) && inherits(params$removeUnkown, "logical")) {
    message("The `removeUnkown` has been deprecated, Please use `removeUnknown` instead!")
    removeUnknown <- params$removeUnkown
  }
  
  featurelist <- keep_know_taxa(featurelist, removeUnknown = removeUnknown)
  nodedfres <- nodedfres[nodedfres$f %in% featurelist, , drop = FALSE]
  nodedfres <- nodedfres[order(nodedfres[, 2], decreasing = TRUE), , drop = FALSE]
  nodedfres$f <- factor(nodedfres$f, levels = as.vector(nodedfres$f))
  featurelist <- as.vector(nodedfres$f)
  
  featureda <- featureda[, match(featurelist, colnames(featureda)), drop = FALSE]
  
  if (is.null(colorlist)) {
    colorlist <- get_cols(length(tmpgroup))
  }
  if (is.null(names(colorlist))) {
    names(colorlist) <- tmpgroup
  }
  
  if (is.null(extract_args(obj, "standard_method"))) {
    if (is.null(l_xlabtext)) {
      xlabtext <- "abundance"
    } else {
      xlabtext <- l_xlabtext
    }
  } else {
    xlabtext <- "abundance"
  }
  
  p <- plotdiffbox(obj = featureda, sampleda = sampleda, factorNames = classname, 
                   factorLevels = factorLevels, featurelist = featurelist, 
                   geom = geom, box_notch = box_notch, dodge_width = dodge_width, 
                   box_width = box_width) +
    coord_flip() +
    scale_fill_manual(values = colorlist) +
    ylab(xlabtext)
  
  if (addLDA) {
    if ("LDAmean" %in% colnames(obj@mlres)) {
      effectsizename <- "LDA"
    } else {
      effectsizename <- "MeanDecreaseAccuracy"
    }
    colorlist <- colorlist[unique(as.vector(nodedfres[[classname]]))]
    p2 <- ggeffectsize.data.frame(obj = nodedfres, factorName = classname,
                                  effectsizename = effectsizename,
                                  setFacet = FALSE, ...) +
      scale_color_manual(values = colorlist) +
      theme(legend.position = "none",
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            plot.margin = unit(c(2, 2, 2, 0), "mm"))
    
    p <- p + p2 + plot_layout(guides = 'collect', widths = c(3, 2))
  }
  
  return(p)
}

### Microbiome Process ggeffect size dataframe ###

ggeffectsize.data.frame <- function(obj, 
                                    factorName, 
                                    effectsizename,
                                    factorLevels=NULL,
                                    linecolor="grey50",
                                    linewidth=0.4,
                                    lineheight=0.2,
                                    pointsize=1.5,
                                    setFacet=TRUE,
                                    ...){
  if (effectsizename %in% "LDA"){
    xlabtext <- bquote(paste(Log[10],"(",.("LDA"), ")"))
    xtext <- "LDAmean"
    xmintext <- "LDAlower"
    xmaxtext <- "LDAupper"
  }else{
    xlabtext <- effectsizename
    xtext <- "MDAmean"
    xmintext <- "MDAlower"
    xmaxtext <- "MDAupper"
  }
  if (!is.null(factorLevels)){
    obj <- setfactorlevels(obj,factorLevels)
  }
  
  p <- ggplot(data=obj, aes_(y=~f)) +
    geom_errorbarh(aes_(xmin=as.formula(paste0("~", xmintext)), 
                        xmax=as.formula(paste0("~", xmaxtext))),
                   color=linecolor, size=linewidth, height=lineheight) +
    geom_point(aes_(x=as.formula(paste0("~",xtext)),
                    color=as.formula(paste0("~",factorName))),
               size=pointsize) 
  if (setFacet) {
    p <- p + facet_grid(as.formula(paste0(factorName," ~.")),
                        scales = "free_y", space = "free_y")
  }
  p <- p + ylab(NULL) + xlab(xlabtext) 
  tmpn <- length(unique(as.vector(obj[[factorName]])))
  message("The color has been set automatically, you can reset it manually by adding scale_color_manual(values=yourcolors)")
  p <- p + scale_color_manual(values=get_cols(tmpn))
  #}
  p <- p + theme_bw() +
    theme(
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      panel.grid.major.y = element_line(color = "grey", linewidth=0.2),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      strip.background = element_rect(colour = NA, fill = "grey")
    )
  return(p)
}



### Breaking the function up into a few separate parts 

# p1 boxplot RA
make_p1_boxplot <- function(obj, featurelist, factorLevels = NULL,
                            box_notch = TRUE, box_width = 0.05, dodge_width = 0.6,
                            colorlist = NULL, xlabtext = "abundance") {
  
  featureda <- obj@originalD
  classname <- extract_args(obj, "classgroup")
  normalization <- extract_args(obj, "normalization")
  
  if (!is.null(normalization)) {
    featureda <- featureda / normalization
  }
  
  sampleda <- obj@sampleda
  tmpgroup <- unique(as.vector(sampleda[[classname]]))
  featureda <- featureda[, match(featurelist, colnames(featureda)), drop = FALSE]
  
  p1 <- plotdiffbox(obj = featureda, sampleda = sampleda, factorNames = classname,
                    factorLevels = factorLevels, featurelist = featurelist,
                    box_notch = box_notch, dodge_width = dodge_width,
                    box_width = box_width) +
    coord_flip() +
    scale_color_manual(values = colorlist) +
    ylab(xlabtext) +
    theme(legend.position = "none",
          panel.grid.major.y = element_line(color = "grey", linewidth=0.2),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())
  
  print(p1$mapping)
  return(p1)
}

# Descending RA order 
make_p1_boxplot_ordered <- function(obj, featurelist, factorLevels = NULL,
                            box_notch = TRUE, box_width = 0.05, dodge_width = 0.6,
                            colorlist = NULL, xlabtext = "abundance") {
  
  featureda <- obj@originalD
  classname <- extract_args(obj, "classgroup")
  normalization <- extract_args(obj, "normalization")
  
  if (!is.null(normalization)) {
    featureda <- featureda / normalization
  }
  
  sampleda <- obj@sampleda
  tmpgroup <- unique(as.vector(sampleda[[classname]]))
  
  # Subset to the features you're interested in
  featureda <- featureda[, match(featurelist, colnames(featureda)), drop = FALSE]
  
  ## Sort features by mean relative abundance ---
  feature_means <- colMeans(featureda, na.rm = TRUE)
  ordered_features <- names(sort(feature_means, decreasing = TRUE))
  featurelist <- ordered_features  # overwrite with sorted order
  ## --------------------------------------------------------
  
  # Continue with plotting
  p1 <- plotdiffbox(obj = featureda, sampleda = sampleda, factorNames = classname,
                    factorLevels = factorLevels, featurelist = featurelist,
                    box_notch = box_notch, dodge_width = dodge_width,
                    box_width = box_width) +
    coord_flip() +
    scale_color_manual(values = colorlist) +
    ylab(xlabtext) +
    theme(legend.position = "none",
          panel.grid.major.y = element_line(color = "grey", linewidth = 0.2),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())
  
  return(p1)
}



# P2 Effect size plot function 
make_p2_effectsize <- function(obj, nodedfres, colorlist) {
  classname <- extract_args(obj, "classgroup")
  effectsizename <- if ("LDAmean" %in% colnames(obj@mlres)) "LDA" else "MeanDecreaseAccuracy"
  
  p2 <- ggeffectsize.data.frame(obj = nodedfres, factorName = classname,
                                effectsizename = effectsizename,
                                setFacet = FALSE) +
    scale_color_manual(values = colorlist) +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          plot.margin = unit(c(2, 2, 2, 0), "mm"))
  
  return(p2)
}

make_p2_effectsize <- function(obj, nodedfres, colorlist, xlim_range = NULL) {
  classname <- extract_args(obj, "classgroup")
  effectsizename <- if ("LDAmean" %in% colnames(obj@mlres)) "LDA" else "MeanDecreaseAccuracy"
  
  p2 <- ggeffectsize.data.frame(
    obj = nodedfres,
    factorName = classname,
    effectsizename = effectsizename,
    setFacet = FALSE
  ) +
    scale_color_manual(values = colorlist) +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      plot.margin = unit(c(2, 2, 2, 0), "mm")
    )
  # Optional x-axis limit control
  if (!is.null(xlim_range)) {
    p2 <- p2 + xlim(xlim_range)
  }
  return(p2)
}

# P3 ANCOMBS plot function 
make_p3_ancombc <- function(ancombc_df, dataset_name, tax_level, featurelist, colorlist) {
  ancombc_df$taxon <- factor(ancombc_df$taxon, levels = featurelist)
  
  p3 <- ggplot(ancombc_df, aes(x = taxon, y = lfc, fill = group)) +
    geom_col(position = position_dodge(width = 0.9)) +
    coord_flip() +
    labs(title = paste("ANCOMBC2:", dataset_name, "-", tax_level),
         x = NULL, y = "LFC (ANCOMBC2)", fill = "Group") +
    scale_fill_manual(values = colorlist, na.translate = FALSE) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      panel.grid.major.y = element_line(color = "grey", linewidth=0.2),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      strip.background = element_rect(colour = NA, fill = "grey")
    )
  
  return(p3)
}


# Wrapper function 

ggdiffbox_combined <- function(obj, ancombc_df = NULL,
                               dataset_name = NULL, tax_level = NULL,
                               box_notch = TRUE, box_width = 0.05, dodge_width = 0.6,
                               factorLevels = NULL, featurelist = NULL,
                               removeUnknown = TRUE, colorlist = NULL,
                               l_xlabtext = NULL, addLDA = TRUE) {
  
  featureda <- obj@originalD
  classname <- extract_args(obj, "classgroup")
  sampleda <- obj@sampleda
  tmpgroup <- unique(as.vector(sampleda[[classname]]))
  
  nodedfres <- obj@result
  nodedfres <- set_newlevels(data = nodedfres, newlevels = tmpgroup, factorNames = classname)
  
  if (is.null(featurelist)) {
    featurelist <- unique(as.vector(nodedfres$f))
  }
  
  featurelist <- keep_know_taxa(featurelist, removeUnknown = removeUnknown)
  nodedfres <- nodedfres[nodedfres$f %in% featurelist, , drop = FALSE]
  nodedfres <- nodedfres[order(nodedfres[, 2], decreasing = TRUE), , drop = FALSE]
  nodedfres$f <- factor(nodedfres$f, levels = as.vector(nodedfres$f))
  featurelist <- as.vector(nodedfres$f)
  
  if (is.null(colorlist)) {
    colorlist <- get_cols(length(tmpgroup))
  }
  if (is.null(names(colorlist))) {
    names(colorlist) <- tmpgroup
  }
  
  xlabtext <- ifelse(is.null(l_xlabtext), "abundance", l_xlabtext)
  
  ## Make plots
  p1 <- make_p1_boxplot(obj, featurelist, factorLevels, box_notch, box_width, dodge_width, colorlist, xlabtext)
  
  if (addLDA) {
    colorlist_sub <- colorlist[unique(as.vector(nodedfres[[classname]]))]
    p2 <- make_p2_effectsize(obj, nodedfres, colorlist_sub)
  } else {
    p2 <- NULL
  }
  
  if (!is.null(ancombc_df)) {
    p3 <- make_p3_ancombc(ancombc_df, dataset_name, tax_level, featurelist, colorlist)
  } else {
    p3 <- NULL
  }
  
  ## Combine plots
  if (!is.null(p2) && !is.null(p3)) {
    return(p1 + p2 + p3 + plot_layout(guides = 'collect', widths = c(3, 2, 2)))
  } else if (!is.null(p2)) {
    return(p1 + p2 + plot_layout(guides = 'collect', widths = c(3, 2)))
  } else if (!is.null(p3)) {
    return(p1 + p3 + plot_layout(guides = 'collect', widths = c(3, 2)))
  } else {
    return(p1)
  }
}



# Sections 6:




### From Microbiom Process Utils 

#' @author GhuangChuangYu
#' @importFrom grDevices colorRampPalette
#' @keywords internal
# this is from `ggtree`
get_cols <- function (n){
  col <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
           "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
           "#ccebc5", "#ffed6f")
  col2 <- c("#1f78b4", "#ffff33", "#c2a5cf", "#ff7f00", "#810f7c",
            "#a6cee3", "#006d2c", "#4d4d4d", "#8c510a", "#d73027",
            "#78c679", "#7f0000", "#41b6c4", "#e7298a", "#54278f")
  col3 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
            "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
            "#ffff99", "#b15928")
  colorRampPalette(col2)(n)
}

#' @keywords internal 
setfactorlevels <- function(data, factorlist){
  factornames <- intersect(colnames(data), names(factorlist))
  if (length(factornames)>0){
    for(i in factornames){
      data[[match(i,colnames(data))]] <- factor(data[[match(i, colnames(data))]], 
                                                levels=as.vector(factorlist[[match(i,names(factorlist))]]))
    }
  }
  return(data)
}

#' @keywords internal
get_otudata <- function(obj){
  otudata <- obj@otu_table
  otudata <- data.frame(otudata, check.names=FALSE)
  if(obj@otu_table@taxa_are_rows){
    otudata <- data.frame(t(otudata), check.names=FALSE)
  }
  return (otudata)
}

#' @keywords internal
checkotu <- function(obj){
  if (is.null(obj@otu_table)){
    stop("The otu table is empty!")
  }else{
    otuda <- get_otudata(obj)
    return(otuda)
  }
}

#' @keywords internal
checksample <- function(obj){
  if (is.null(obj@sam_data)){
    stop("The sample_data is empty")
  }else{
    sampleda <- get_sample(obj)
    return(sampleda)
  }
}

#' @keywords internal.
get_sample <- function(obj){
  if (is.null(obj@sam_data)){
    sampleda <- NULL
  }else{
    sampleda <- data.frame(obj@sam_data, check.names=FALSE)
  }
  return(sampleda)
}

# #' @keywords internal
# #taxlevelchar <- c("k", "p", "c", "o", "f", "g", "s", "st")

newtaxname <- function(x, y){
  y <- as.vector(y)
  x[y] <- paste(taxlevelchar[y], x[y], sep="__un_")
  x
}

#' @importFrom zoo na.locf
#' @keywords internal
filltaxname <- function(taxdf, type="species"){
  tmprownames <- rownames(taxdf)
  indexmark <- apply(taxdf, 2, function(x){nchar(x, keepNA = TRUE)})<=4
  taxdf[indexmark] <- NA
  if (any(is.na(taxdf[,1]))){
    if (type == "species"){
      prefix <- "k__"
    }else{
      prefix <- "d1__"
    }
    taxdf[is.na(taxdf[,1]), 1] <- paste0(prefix, "Unknown")
  }    
  indextmp <- apply(is.na(taxdf), 1, which)
  if(length(indextmp)==0){
    taxdf <- data.frame(taxdf, check.names=FALSE)
    return(taxdf)
  }
  taxdf <- apply(taxdf, 1, na.locf)
  taxdf <- lapply(seq_len(ncol(taxdf)), function(i) taxdf[,i])
  taxdf <- data.frame(t(mapply(newtaxname, taxdf, indextmp)), 
                      stringsAsFactors=FALSE)
  rownames(taxdf) <- tmprownames
  return(taxdf)
}

#' @keywords internal
addtaxlevel <- function(taxdf){#, type="species"){
  #if (type != "species"){
  #    taxlevelchar <- paste0("d", seq_len(ncol(taxdf)))
  #}else{
  #    taxlevelchar <- taxlevelchar[seq_len(ncol(taxdf))]
  #}
  taxlevelchar <- taxlevelchar[seq_len(length(taxdf))]
  paste(taxlevelchar, taxdf, sep="__")
}

#' @importFrom tibble column_to_rownames
#' @keywords internal
fillNAtax <- function(taxdf, type="species"){
  #taxdf <- remove_unclassfied(taxdf)
  taxdf <- remove_na_taxonomy_rank(taxdf)
  if (type!="species"){
    assign("taxlevelchar", paste0("d", seq_len(ncol(taxdf))), envir = .GlobalEnv)
  }else{
    assign("taxlevelchar", c("k", "p", "c", "o", "f", "g", "s", "st"), envir = .GlobalEnv)
  }
  #    if (any(is.na(taxdf[,1]))){
  #        if (type == "species"){
  #            prefix <- "k__"
  #        }else{
  #            prefix <- "d1__"
  #        }
  #        taxdf[is.na(taxdf[,1]), 1] <- paste0(prefix, "Unknown")
  #    }
  if (!(grepl("^k__", taxdf[1,1]) || grepl("^d1__", taxdf[1,1]))){
    tmprownames <- rownames(taxdf)
    tmpcolnames <- colnames(taxdf)
    taxdf <- t(apply(taxdf, 1, as.character))
    taxdf[is.na(taxdf)] <- ""
    taxdf <- data.frame(t(apply(taxdf, 1, addtaxlevel)),
                        stringsAsFactors=FALSE)
    rownames(taxdf) <- tmprownames
    colnames(taxdf) <- tmpcolnames
  }
  taxdf <- filltaxname(taxdf, type = type)
  taxdf <- repduplicatedtaxcheck(taxdf) #%>% column_to_rownames(var="rowname")
  attr(taxdf, "fillNAtax") <- TRUE 
  return(taxdf)
}

#' @keywords internal
remove_unclassfied <- function(taxdf){
  taxdf[grepl.data.frame("Unclassified|uncultured|Ambiguous|Unknown|unknown|metagenome|Unassig", taxdf, ignore.case=TRUE)] <- NA
  return(taxdf)
}

grepl.data.frame <- function(pattern, x, ...){
  y <- if (length(x)) {
    do.call("cbind", lapply(x, "grepl", pattern=pattern, ...))
  }else{
    matrix(FALSE, length(row.names(x)), 0)
  }
  if (.row_names_info(x) > 0L)
    rownames(y) <- row.names(x)
  y
}


#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @keywords internal
duplicatedtaxcheck <- function(taxdf){
  if (ncol(taxdf)==1){return(taxdf)}
  taxdf <- taxdf %>% rownames_to_column()
  for (i in ncol(taxdf):3){
    tmp <- split(taxdf,taxdf[,i])
    for (j in seq_len(length(tmp))){
      flag <- length(unique(as.vector(tmp[[j]][,i-1])))
      if (flag > 1){
        tmp[[j]][,i] <- paste(tmp[[j]][,i],tmp[[j]][,i-1],sep="_")
      }
    }
    taxdf <- do.call("rbind",c(tmp, make.row.names=FALSE)) 
  }
  return(taxdf)
}

#' @keywords internal
repduplicatedtaxcheck <- function(taxdf){
  for (i in seq_len(7)){
    taxdf <- duplicatedtaxcheck(taxdf) %>% 
      column_to_rownames(var="rowname")
  }
  return(taxdf)
}

## #' @keywords internal
## ## reference https://rdrr.io/cran/stackoverflow/man/match.call.defaults.html
## match.call.defaults <- function(fun) {
##     if (!is.na(fun)){
##         print(args(diff_analysis.data.frame))
##         args(diff_analysis.data.frame)
##     }else{
##         call <- evalq(match.call(expand.dots=TRUE), parent.frame(1))
##         formals <- evalq(formals(), parent.frame(1))
##         for(i in setdiff(names(formals), c(names(call)))){
##             call[i] <- list(formals[[i]])
##         }
##         match.call(sys.function(sys.parent()), call)
##     }
## }

#' @keywords internal
extract_args <- function(obj, arg){
  if (!"someparams" %in% methods::slotNames(obj)){
    stop("The object don't have someparams slot!")
  }else{
    args <- obj@someparams
    argres <- args[[arg]]
    return(argres)
  }
}

#' @importFrom utils globalVariables
utils::globalVariables('taxlevelchar')

#' @importFrom stats sd 
#' @importFrom plyr ddply
# Adapted from Rmisc
summarySE <- function (data = NULL, measurevar, groupvars = NULL, na.rm = FALSE,
                       conf.interval = 0.95, .drop = TRUE){
  length2 <- function(x, na.rm = FALSE) {
    if (na.rm)
      sum(!is.na(x))
    else length(x)
  }
  datac <- ddply(data, 
                 groupvars, 
                 .drop = .drop, 
                 .fun = function(xx, col, na.rm) {
                   c(N = length2(xx[, col], na.rm = na.rm), mean = mean(xx[, col], na.rm = na.rm), 
                     sd = sd(xx[, col], na.rm = na.rm))
                 }, measurevar, na.rm)
  datac %<>% dplyr::rename(!!measurevar:="mean")
  datac$se <- datac$sd/sqrt(datac$N)
  ciMult <- qt(conf.interval/2 + 0.5, datac$N - 1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

#' @importFrom stats qt
# reference to Rmisc
CI <- function (x, ci = 0.95, na.rm=FALSE){
  a <- mean(x, na.rm = na.rm)
  s <- sd(x, na.rm = na.rm)
  if (na.rm){
    x <- x[!is.na(x)]
  }
  n <- length(x)
  error <- qt(ci + (1 - ci)/2, df = n - 1) * s/sqrt(n)
  return(c(upper = a + error, mean = a, lower = a - error))
}

.return_wrap <- function(...){
  msg <- paste(..., collapse = "", sep = "")
  wrapped <- strwrap(msg, width = getOption("width") - 2) %>%
    glue::glue_collapse(., "\n", last = "\n")
  wrapped
}

message_wrap <- function(...){
  msg <- .return_wrap(...)
  message(msg)
}

stop_wrap <- function(...){
  msg <- .return_wrap(...)
  stop(msg, call. = FALSE)
}

warning_wrap <- function(...){
  msg <- .return_wrap(...)
  warning(msg, call. = FALSE)
}

remove_MP_internal_res <- function(x){
  x <- x[,!vapply(x, function(i)is.list(i), logical(1))]
  index <- lapply(MP_internal_res, function(i)which(grepl(paste0("^",i), colnames(x)))) %>%
    unlist() %>%
    unique()
  if (length(index) > 0) x <- x[-index]
  return(x)
}

MP_internal_res <- c("Observe", "Chao1", "ACE", "Shannon", 
                     "Simpson", "Pielou", 'PC', 'PCo', 
                     'CA', 'NMDS', 'RDA', 'CCA', "PAE", "NRI", 
                     "NTI", "PD", 'IAC', "HAED", "EAED")