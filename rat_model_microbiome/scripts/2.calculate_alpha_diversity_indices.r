
library("tools")
library("phyloseq")
library("metagenomeSeq")

##### PARAMETERS ##############################################
basefolder = "/home/filippo/Documents/cremonesi/tunisian_rats/"
metadataf = "microbiomes/rat_model_microbiome/Config/mapping_file_colitis.xlsx"
sample_prefix = "sample"
outdir = "colitis/Analysis/"
###############################################################

writeLines(" - loading data")
## loading data previously imported in phyloseq
fname = file.path(outdir, "phyloseq.RData")
load(fname)

writeLines(" - calculate alpha diversity indices")
## alpha div estimated on not normalised OTU counts: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531
alpha = estimate_richness(otu_tax_sample, split = TRUE)
alpha$"sample-id" = row.names(alpha)
alpha = relocate(alpha, `sample-id`)

writeLines(" - writing out the alpha diversity indices")
fname = file.path(basefolder, outdir, "alpha.csv")
fwrite(x = alpha, file = fname)


