
library("tools")
library("phyloseq")
library("metagenomeSeq")

##### PARAMETERS ##############################################
basefolder = "/home/filippo/Documents/cremonesi/tunisian_rats/"
metadataf = "microbiomes/rat_model_microbiome/Config/mapping_file_colitis.xlsx"
sample_prefix = "sample"
outdir = "colitis/Analysis/"
###############################################################

## repo with supporting functions
HOME <- Sys.getenv("HOME")
repo = file.path(HOME, "Documents/cremonesi/metabarcoding")
source(file.path(repo, "support_functions/dist2list.R")) ## from: https://github.com/vmikk/metagMisc/
source(file.path(repo, "support_functions/phyloseq_transform.R")) ## from: https://github.com/vmikk/metagMisc/


writeLines(" - loading data")
## loading data previously imported in phyloseq
fname = file.path(outdir, "phyloseq.RData")
load(fname)

## Preprocessing: e.g. filtering
## making normalization folder
if(!file.exists(file.path(basefolder, outdir))) dir.create(file.path(outdir), showWarnings = FALSE)

writeLines(" - CSS normalization")
otu_tax_sample_norm = phyloseq_transform_css(otu_tax_sample, norm = TRUE, log = FALSE)

## add taxonomy to normalised counts
otu_css_norm = base::as.data.frame(otu_table(otu_tax_sample_norm))
otu_css_norm$tax_id = row.names(otu_css_norm)
otu_css_norm <- relocate(otu_css_norm, tax_id)
taxonomy = as.data.frame(tax_table(otu_tax_sample_norm))
taxonomy$tax_id = row.names(taxonomy)
taxonomy <- relocate(taxonomy, tax_id)
otu_css_norm = otu_css_norm %>% inner_join(taxonomy, by = "tax_id")

writeLines(" - writing out the CSS normalized OTU table")
fwrite(x = otu_css_norm, file = file.path(basefolder, outdir, "otu_norm_CSS.csv"))

## BETA DIVERSITY ##
writeLines(" - beta diversity: distance matrices")
dist_methods <- unlist(distanceMethodList)
# print(dist_methods)

## bray-curtis
writeLines(" - calculate Bray-Curtis distances")
distances = distance(otu_tax_sample_norm, method="bray", type = "samples")
iMDS  <- ordinate(otu_tax_sample_norm, "MDS", distance=distances)
plot_ordination(otu_tax_sample_norm, iMDS, color="treatment")

writeLines(" - write out distance matrix")
dd = dist2list(distances, tri = FALSE)
sample_order = gtools::mixedsort(as.character(unique(dd$row)), decreasing = TRUE)
dd$row = factor(dd$row, levels = sample_order)
dd$col = factor(dd$col, levels = sample_order)
dx <- dd |> pivot_wider(names_from = "col", values_from = "value")
fwrite(x = dx, file = file.path(basefolder, outdir, "bray_curtis_distances.csv"))

print("DONE!")

