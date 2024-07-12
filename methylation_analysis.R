# Load necessary libraries
library(minfi)
library(limma)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(ggplot2)

# Function to perform DNA methylation analysis
perform_methylation_analysis <- function(idatDir, outputDir) {
  
  # Set output directory
  if (!dir.exists(outputDir)) dir.create(outputDir)
  
  # Read IDAT files into an RGChannelSet object
  rgSet <- read.metharray.exp(base = idatDir, recursive = TRUE)
  
  # Perform preprocessing
  rgSetFiltered <- preprocessQuantile(rgSet)
  
  # Define groups (modify according to your groups)
  group <- factor(c(rep("Disease", nDiseaseSamples), rep("Normal", nNormalSamples)))
  
  # Create design matrix
  design <- model.matrix(~ group)
  
  # Perform differential methylation analysis
  fit <- lmFit(rgSetFiltered, design)
  contrastMatrix <- makeContrasts(Disease - Normal, levels = design)
  fitContrast <- contrasts.fit(fit, contrastMatrix)
  fitContrast <- eBayes(fitContrast)
  topDMRs <- topTable(fitContrast, coef = 1, number = Inf)
  
  # Save results
  write.table(topDMRs, file.path(outputDir, "topDMRs.csv"), sep = ",", quote = FALSE, row.names = TRUE)
  
  # Generate plots
  pdf(file.path(outputDir, "methylation_plots.pdf"))
  volcanoPlot(fitContrast, coef = 1, highlight = Inf)
  dev.off()
  
  # MA plot
  pdf(file.path(outputDir, "ma_plot.pdf"))
  plotMA(fitContrast, coef = 1)
  dev.off()
  
  # Heatmap of methylation patterns
  methylation_values <- assay(rgSetFiltered)
  pdf(file.path(outputDir, "methylation_heatmap.pdf"))
  heatmap(methylation_values)
  dev.off()
  
  # Print completion message
  cat("Analysis complete. Results saved in", outputDir, "\n")
}

# Command-line argument handling
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  cat("Usage: Rscript methylation_analysis.R <input_directory> <output_directory>\n")
  quit(status = 1)
}

# Input arguments
inputDir <- args[1]
outputDir <- args[2]

# Perform analysis
perform_methylation_analysis(inputDir, outputDir)

