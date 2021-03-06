### Plot Illumina picard metrics
### librarys, plot functions and helper functions.

# Required packages
suppressMessages(library(ggplot2))
suppressMessages(library(knitr))
suppressMessages(library(markdown))
suppressMessages(library(reshape))
suppressMessages(library(xtable))
suppressMessages(library(tools))
suppressMessages(library(brew))

# Plotting functions

#Custom colorscale used for plotting
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
GetRandomColorSet <- function(nsamples) {
  SetTextContrastColor <- function(color) {
    ifelse( mean(col2rgb(color)) > 130, "black", "white")
  }
  TextContrastColor <- unlist( lapply(colors(), SetTextContrastColor) )
  color_set <- colors()[TextContrastColor=="white"][append(-grep("grey",colors()[TextContrastColor == "white"]),-grep("gray",colors()[TextContrastColor=="white"]))]
  random_color_set <- color_set[sample(1:length(color_set),nsamples)]
  return(random_color_set)
}

plot_insert_size_metrics <- function(){
  ggplot(insert_size_metrics.table, aes(x=insert_size, y=All_Reads.fr_count)) + 
    geom_bar(stat="identity", width=1, fill="#0072B2") +
    xlab("Insert Size") + ylab("Count") +
    scale_fill_manual(name="", values=colorSet)+
    ggtitle(paste("Insert size for all reads in", samples[i], sep=" ")) +
    theme(axis.title = element_text(face="bold", size=15),
          axis.text = element_text(size=15),
          plot.title = element_text(size=15, face ="bold"))
}

plot_quality_by_cycle_metrics <- function(){
  ggplot(quality_by_cycle_metrics.table, aes(x=CYCLE, y=MEAN_QUALITY)) + 
    geom_bar(stat="identity", width=1, fill="#0072B2") +
    xlab("Cycle") + ylab("Mean Quality") +
    scale_fill_manual(name="", values=colorSet)+
    ggtitle(paste("Quality by cycle in", samples[i], sep=" ")) +
    theme(axis.title = element_text(face="bold", size=15),
          axis.text = element_text(size=15),
          plot.title = element_text(size=15, face ="bold"))
}

plot_quality_distribution_metrics <- function(){
  ggplot(quality_distribution_metrics.table, aes(x=QUALITY, y=COUNT_OF_Q)) + 
    geom_bar(stat="identity", fill="#0072B2") +
    xlab("Quality Score") + ylab("Observations") +
    scale_fill_manual(name="", values=colorSet)+
    ggtitle(paste("Quality score distribution in", samples[i], sep=" ")) +
    theme(axis.title = element_text(face="bold", size=15),
          axis.text = element_text(size=15),
          plot.title = element_text(size=15, face ="bold"))
}

plot_gcMetricsBaseQuality <- function(){
  ggplot(gc_bias_metrics.table, aes(x=GC, y=MEAN_BASE_QUALITY))+
    geom_path(colour="#0072B2")+
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 40))+
    xlab("GC%") + ylab("Mean Base Quality")+
    ggtitle(paste("Mean base quality per GC% in", samples[i], sep=" ")) +
    theme(axis.title = element_text(face="bold", size=15),
          axis.text = element_text(size=15),
          plot.title = element_text(size=15, face ="bold"))
}

limits = aes(ymax=NORMALIZED_COVERAGE + (ERROR_BAR_WIDTH/2), ymin=NORMALIZED_COVERAGE-(ERROR_BAR_WIDTH/2))

plot_gcMetricsNormalizedCoverage <- function(){
  ggplot(gc_bias_metrics.table, aes(x=GC, y=NORMALIZED_COVERAGE))+
    geom_point(colour="#0072B2")+
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 3))+
    geom_errorbar(limits, width=0.2, color="#E69F00")+
    xlab("GC%") + ylab("Normalized Coverage")+
    ggtitle(paste("Normalized coverage per GC% in", samples[i], sep=" "))+
    theme(axis.title = element_text(face="bold", size=15),
          axis.text = element_text(size=15),
          plot.title = element_text(size=15, face ="bold"))
}

plot_wgs_metrics <- function(){
  ggplot(wgsMetrics_sample.table, aes(x=coverage, y=count)) + 
    geom_bar(stat="identity", width=1, fill="#0072B2") +
    xlab("Coverage") + ylab("Count") +
    scale_fill_manual(name="", values=colorSet)+
    ggtitle(paste("Coverage histogram", samples[i], sep=" ")) +
    theme(axis.title = element_text(face="bold", size=15),
          axis.text = element_text(size=15),
          plot.title = element_text(size=15, face ="bold"))
}

plot_pctOffBait <- function(){
  ggplot(summaryTable, aes(x=sample, y=PCT_OFF_BAIT, fill=sample)) + 
    geom_bar(stat="identity") +
    xlab("Sample") + ylab("Percentage off bait") +
    scale_fill_manual(name="", values=colorSet)+
    ggtitle("Percentage off bait") +
    theme(axis.title = element_text(face="bold", size=15),
          axis.text.y = element_text(size=15),
          axis.text.x = element_blank(),
          legend.text = element_text(size=15),
          plot.title = element_text(size=15, face ="bold"))
}

plot_meanTargetCov <- function(){
  ggplot(summaryTable, aes(x=sample, y=MEAN_TARGET_COVERAGE, fill=sample)) + 
    geom_bar(stat="identity") +
    xlab("Sample") + ylab("Mean target coverage") +
    scale_fill_manual(name="", values=colorSet) +
    ggtitle("Mean target coverage") +
    theme(axis.title = element_text(face="bold", size=15),
          axis.text.y = element_text(size=15),
          axis.text.x = element_blank(),
          legend.text = element_text(size=15),
          plot.title = element_text(size=15, face ="bold"))
}

plot_pctTargetBases <- function(){
  ggplot(summaryTableMelted,aes(x = sample, y = value)) + 
    geom_bar(aes(fill=variable), stat="identity",position = "dodge") +
    xlab("Sample") + ylab("Percentage") +
    scale_fill_manual(name="", values=cbPalette) +
    ggtitle("Percentage target bases") +
    theme(axis.title = element_text(face="bold", size=15),
          axis.text.x = element_text(size=15, angle=90),
          axis.text.y = element_text(size=15),
          legend.text = element_text(size=15),
          plot.title = element_text(size=15, face ="bold"))
}

# Brew/Tex helper functions
include_graph <- function(width = 1, filename) { 
  paste("\\includegraphics[width=", width, "\\linewidth]{", filename, "}", sep = "") 
}
include_tbl <- function(tableName) {
  colNumber = ncol(tableName)
  colNameMean = mean(nchar(colnames(tableName)))
  splitTableN = floor(75/colNameMean)  
  splitTable = split(1:colNumber, rep(1:colNumber,each=splitTableN,length=colNumber))
  for (i in 1:length(splitTable)){
    tempTable = as.matrix(tableName[, unlist(splitTable[i],use.names=F), drop=FALSE])
    colnames(tempTable) = colnames(tableName)[unlist(splitTable[i],use.names=F)]
    print(xtable(tempTable), type="latex",table.placement = "H")
  }
}
subfloat_graph <- function(width, filename) { 
  paste("\\subfloat{", "\\begin{minipage}[h]{", width, "\\linewidth}\\centering", include_graph(width = 1, filename), "\\end{minipage}}", sep = "")
}
