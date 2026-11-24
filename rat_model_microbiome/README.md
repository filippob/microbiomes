### Repo

Repository for code, workflow and metadata related to two projects on the gut microbiome of rats used as models for the following diseases:
- colitis
- gastric ulcer

### Workflow

1. `read_into_phyloseq.r` ([here](scripts/1.read_into_phyloseq.r)): metadata + OTU/ASV table $\rightarrow$ `phyloseq.RData`
2. `calculate_alpha_diversity_indices.r` ([here](scripts/2.calculate_alpha_diversity_indices.r)): `phyloseq.RData` $\rightarrow$ `alpha.csv`
3. 
