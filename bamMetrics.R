### Plot picard stats to report
### Robert Ernst

### Load packages
library(GetoptLong)

### Command line arguments
GetoptLong(c(
  "samples=s@", "Sample names",
  "output_dir=s", "Output directory",
  "root_dir=s", "bamMetrics root directory",
  "run_name=s", "Run name to use as title for bamMetrics output"
))

#debug
#samples = c("CONTROLP25","CONTROLP26") #debug
#root_dir = "/hpc/cog_bioinf/data/robert/scripts/bamMetrics/" #debug
#output_dir = "/hpc/cog_bioinf/data/robert/testIAP/testData/bams/bamMetrics" #debug
#run_name = "bamMetrics"

### Load functions 
source(paste(root_dir,"bamMetrics_include.R",sep="/"))

### Set wd and create temp output folders
setwd(output_dir)
dir.create("pdfFigures", showWarnings=F) #create output Dir

### Set colorscale
if (length(samples) <= length(cbPalette)){
  colorSet = cbPalette
} else {
  colorSet = GetRandomColorSet(length(samples))
}

### Parse hsmetrics table
fileName = paste(output_dir,"HSMetrics_summary.txt",sep="/")
hsMetrics = FALSE
if (file.exists(fileName)){
  hsMetrics = TRUE
  summaryTable = read.table(file=fileName, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  summaryTable$sample = as.character(summaryTable$sample)
  summaryTableMelted = melt(summaryTable[,c('sample','PCT_TARGET_BASES_2X','PCT_TARGET_BASES_10X','PCT_TARGET_BASES_20X','PCT_TARGET_BASES_30X','PCT_TARGET_BASES_40X','PCT_TARGET_BASES_50X','PCT_TARGET_BASES_100X')],id.vars = 1)
  
  #Transpose and write summaryTable
  summaryTableT = t(summaryTable)
  colnames(summaryTableT) = summaryTableT[1,]
  write.table(summaryTableT, file="HSMetrics_summary.transposed.txt", col.names=FALSE, na="", quote=FALSE, sep="\t")
  
  pdfOut =  paste(output_dir,"pdfFigures", sep="/")
  
  pctOffBait <- plot_pctOffBait()
  ggsave(paste(pdfOut,"pctOffBait.pdf", sep="/"), pctOffBait, dpi = 300, width=15, height=10)
  
  meanTargetCov <- plot_meanTargetCov()
  ggsave(paste(pdfOut,"meanTargetCoverage.pdf", sep="/"), meanTargetCov, dpi = 300, width=15, height=10)
  
  pctTargetBases <- plot_pctTargetBases()
  ggsave(paste(pdfOut,"pctTargetBases.pdf", sep="/"), pctTargetBases, dpi = 300, width=20, height=10)
}

### Parse WGSMetrics table
fileName = paste(output_dir,"WGSMetrics_summary.txt",sep="/")
wgsMetrics = FALSE
if (file.exists(fileName)){
  wgsMetrics = TRUE
  summaryTable = read.table(file=fileName, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  summaryTable$sample = as.character(summaryTable$sample)
  
  #Transpose and write summaryTable
  summaryTableT = t(summaryTable)
  colnames(summaryTableT) = summaryTableT[1,]
  write.table(summaryTableT, file="WGSMetrics_summary.transposed.txt", col.names=FALSE, na="", quote=FALSE, sep="\t")
}

### Parse RNAMetrics table
fileName = paste(output_dir,"RNAMetrics_summary.txt",sep="/")
rnaMetrics = FALSE
if (file.exists(fileName)){
  rnaMetrics = TRUE
  summaryTable = read.table(file=fileName, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  summaryTable$sample = as.character(summaryTable$sample)
  
  #Transpose and write summaryTable
  summaryTableT = t(summaryTable)
  colnames(summaryTableT) = summaryTableT[1,]
  write.table(summaryTableT, file="RNAMetrics_summary.transposed.txt", col.names=FALSE, na="", quote=FALSE, sep="\t")
}

### Plot sample metrics to pdf files.
for(i in 1:length(samples)) {
  #samplePath = paste(samples[i],"QCStats",samples[i],sep="/")
  sample = samples[i]
  
  pdfOut =  paste("pdfFigures", sample, sep="/")
  quality_by_cycle_metrics = paste(sample,"/",sample,"_MultipleMetrics.txt.quality_by_cycle_metrics", sep="")
  quality_distribution_metrics = paste(sample,"/",sample,"_MultipleMetrics.txt.quality_distribution_metrics", sep="")
  
  quality_by_cycle_metrics.table = read.table(file=quality_by_cycle_metrics, head=TRUE)
  quality_distribution_metrics.table = read.table(file=quality_distribution_metrics, head=TRUE)
  
  cycleQuality <- plot_quality_by_cycle_metrics()
  ggsave(paste(pdfOut,"_cycleQuality.pdf", sep=""), cycleQuality, dpi = 300)
  
  qualityDistribution <- plot_quality_distribution_metrics()
  ggsave(paste(pdfOut,"_qualityDistribution.pdf", sep=""), qualityDistribution, dpi = 300)
  
  paired_end = FALSE
  insert_size_metrics = paste(sample,"/",sample,"_MultipleMetrics.txt.insert_size_metrics", sep="")
  if (file.exists(insert_size_metrics)){
    paired_end = TRUE
    insert_size_metrics.table = read.table(file=insert_size_metrics, skip=10, head=TRUE)
    insertSize <- plot_insert_size_metrics()
    ggsave(paste(pdfOut,"_insertSize.pdf", sep=""), insertSize, dpi = 300)
  }
}

### Generate .html based on .Rmd file
options(knitr.unnamed.chunk.label = "bamMetrics")
knit(paste(root_dir,"bamMetrics.Rmd", sep="/"),quiet=TRUE)
markdownToHTML("bamMetrics.md", "bamMetrics.html", options=c("use_xhml"), stylesheet=paste(root_dir,"bamMetrics_html.css",sep="/"))

### Generate .pdf based on .brew file
brew(paste(root_dir,"bamMetrics.brew", sep="/"), "bamMetrics.tex")
texi2dvi("bamMetrics.tex", pdf = TRUE, clean=TRUE)

### Clean R files and rename output
unlink(c("bamMetrics.tex","bamMetrics.md","pdfFigures","Rplots.pdf"), recursive=T)
file.rename("bamMetrics.html",paste(run_name,"bamMetrics.html",sep="."))
file.rename("bamMetrics.pdf",paste(run_name,"bamMetrics.pdf",sep="."))