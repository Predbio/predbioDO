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

extract.raw.data(in.path = dirs, prefix = prefix, out.path = "../haplo_input", array = "gigamuga")

# Read in the data and save to an Rdata file.
setwd("../haplo_input")
x = read.delim("x.txt", row.names = NULL)
y = read.delim("y.txt", row.names = NULL)
g = read.delim("geno.txt", row.names = NULL)
rownames(x) = make.unique(make.names(x[,1]))
x = as.matrix(x[,-1])
rownames(y) = make.unique(make.names(y[,1]))
y = as.matrix(y[,-1])
rownames(g) = make.unique(make.names(g[,1]))
g = as.matrix(g[,-1])
save(x, y, g, file = "x_y_geno.Rdata")
rm(x,y,g)
gc()

# Filter samples with low call rates.
bad.samples = filter.samples(path = ".", thr = 0.9)
write.table(bad.samples, "bad.samples.txt", sep = "\t")

# Read in the data and save to an Rdata file.
setwd("../haplo_input")
x = read.delim("x.filt.txt", row.names = NULL)
y = read.delim("y.filt.txt", row.names = NULL)
g = read.delim("geno.filt.txt", row.names = NULL)
rownames(x) = make.unique(make.names(x[,1]))
x = as.matrix(x[,-1])
rownames(y) = make.unique(make.names(y[,1]))
y = as.matrix(y[,-1])
rownames(g) = make.unique(make.names(g[,1]))
g = as.matrix(g[,-1])
save(x, y, g, file = "x_y_geno_filtered.Rdata")

# Predict sex.
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
png("predbio_sex.png", width = 800, height = 800, res = 128)
sex = sex.predict(x = x, y = y, snps = GM_snps, plot = TRUE)
dev.off()

# Add generations.
#######################
# Get true generations.
#######################
gen = rep(16, nrow(x))
names(gen) = names(sex)

# Write out the data that will be used for haplotype reconstruction.
setwd("..")
data.filename = "haplo_input/GM_input_data.Rdata"
g[g == TRUE] = "T"
g[g == "-"] = "N"
save(x, y, g, sex, gen, file = data.filename)

###
# Haplotype Reconstruction.

# Load in the data.
setwd("/home/rstudio")
load(file = data.filename)

# Craete the data object for DOQTL.
data = list(geno = g, sex = sex, gen = gen)
rm(x, y, g, sex, gen)
gc()

# Run calc.genoprobs.
setwd("HMM")
calc.genoprob(data = data, output.dir = ".",
     plot = FALSE, array = "gigamuga", sampletype = "DO",
     method = "allele")

