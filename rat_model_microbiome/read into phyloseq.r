
library("phyloseq")
library("metagenomeSeq")

writeLines(" - reading the filtered (OTU-wise) biom file into phyloseq")
fname = "otu_table_filtered-138.biom"
biom_otu_tax <- phyloseq::import_biom(BIOMfilename = fname)

otu = otu_table(biom_otu_tax, taxa_are_rows = TRUE)
nrow(otu) ## n. rows 
ncol(otu) ## n. col (samples) BEFORE filtering

writeLines(" - removing samples with too few total counts")
biom_otu_tax = prune_samples(sample_sums(biom_otu_tax) >= 10, biom_otu_tax)

otu = otu_table(biom_otu_tax, taxa_are_rows = TRUE)
ncol(otu) ## n. col (samples) AFTER filtering

taxa = tax_table(biom_otu_tax)

writeLines(" - reading the metadata")
fname = "mapping_file_colitis.csv"
metadata = fread(fname, header = TRUE)

names(metadata)[1] <- "sample-id"
metadata <- as.data.frame(metadata)
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
fname = "phyloseq.RData"
save(otu_tax_sample, file = fname)

