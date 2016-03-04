# Gather all of the GigaMUGA data from all sources.
library(DOQTL)
library(mclust)
setwd("/home/rstudio/geneseek")

# These are the directories containing the data.
dirs = dir()
# Remove any *.gz files.
dirs = dirs[-grep("\\.tar\\.gz$", dirs)]
# Remove the first run because it has too many SNPs.
dirs = dirs[-grep("14Nov2012$", dirs)]

prefix = paste0("set", 1:length(dirs))

extract.raw.data(in.path = dirs, prefix = prefix, out.path = "/home/rstudio/haplo_input", array = "gigamuga")

# Read in the data and save to an Rdata file.
setwd("../haplo_input")
x = read.delim("x.txt", row.names = NULL)
y = read.delim("y.txt")
g = read.delim("g.txt")
save(x, y, g, file = "x_y_geno.Rdata")

