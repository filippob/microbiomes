
library("tools")
library("phyloseq")
library("data.table")
library("metagenomeSeq")

##### PARAMETERS ##############################################
basefolder = "/home/filippo/Documents/cremonesi/tunisian_rats/"
otuf = "colitis/otu_table_filtered/otu_table_filtered-138.biom"
metadataf = "microbiomes/rat_model_microbiome/Config/mapping_file_colitis.xlsx"
sample_prefix = "sample"
outdir = "colitis/Analysis/"
###############################################################

writeLines(" - reading the filtered (OTU-wise) biom file into phyloseq")
fname = file.path(basefolder, otuf)
biom_otu_tax <- phyloseq::import_biom(BIOMfilename = fname)

otu = otu_table(biom_otu_tax, taxa_are_rows = TRUE)
nrow(otu) ## n. rows 
ncol(otu) ## n. col (samples) BEFORE filtering

writeLines(" - removing samples with too few total counts")
biom_otu_tax = prune_samples(sample_sums(biom_otu_tax) >= 10, biom_otu_tax)

otu = otu_table(biom_otu_tax, taxa_are_rows = TRUE)
ncol(otu) ## n. col (samples) AFTER filtering

colnames(otu) <- paste(sample_prefix,colnames(otu),sep="-")
print(head(otu))

taxa = tax_table(biom_otu_tax)

writeLines(" - reading the metadata")
fname = file.path(basefolder, metadataf)
if (file_ext(fname) == "xlsx") {
  
  metadata = readxl::read_xlsx(fname, sheet = 1, col_names = TRUE)
} else metadata = fread(fname, header = TRUE)

names(metadata)[1] <- "sample-id"
metadata <- as.data.frame(metadata)
metadata$`sample-id` = paste(sample_prefix, metadata$`sample-id`, sep="-") 
row.names(metadata) <- metadata$`sample-id`
metadata$`sample-id` <- NULL
metadata$treatment = as.factor(metadata$treatment)


## read into phyloseq
writeLines(" - add metadata to the phyloseq object")
samples = sample_data(metadata)

otu_tax_sample = phyloseq(otu,taxa,samples)
sample_data(otu_tax_sample) |> head()
sample_data(otu_tax_sample) |> nrow()

## save phyloseq object
fname = file.path(basefolder, outdir, "phyloseq.RData")
save(otu_tax_sample, file = fname)

print(paste("Annotated OTU table and metadata saved to", fname))

