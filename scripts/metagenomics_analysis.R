### ===============================
## BINF110 Assignment 3
## Taxonomic classification of shotgun metagenomics data from omnivore and vegan gut microbiomes

## Kenneth Gamueda
## 2026-03-23

### === PACKAGES USED ===========
library(phyloseq)
library(biomformat)
library(vegan)
library(ggplot2)

### === 1 | DATA INPUT AND CLEANING ========
# Import BIOM tables (threshold = 0 and 10) from Kraken2/Bracken classification
biom_data_10 <- read_biom("bracken_results/table.biom")
physeq_10 <- import_biom(biom_data_10)
physeq_10

biom_data_0 <- read_biom("bracken_unfilt_results/table_unfilt.biom")
physeq_0 <- import_biom(biom_data_0)
physeq_0

# Attach metadata to phyloseq object
metadata <- data.frame(
  row.names = sample_names(physeq_0),
  diet = c("omn", "omn", "omn", "veg", "veg", "veg")
)
sample_data(physeq_10) <- metadata
sample_data(physeq_0) <- metadata

sample_names(physeq_0) <- gsub("_bracken_species", "", sample_names(physeq_0))

# Transpose table to make rows = samples
otu_table_10 <- as.data.frame(t(otu_table(physeq_10)))
otu_table_10
rare_curve_10 <- rarecurve(otu_table_10, step = 1)

otu_table_0 <- as.data.frame(t(otu_table(physeq_0)))
otu_table_0
rare_curve_0 <- rarecurve(otu_table_0, step = 1)

### === 2 | TAXONOMIC ABUNDANCE ========
# Convert to relative abundance
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# Phylum level
physeq_phy <- tax_glom(physeq_rel, taxrank = "Rank2")
df <- psmelt(physeq_phy)
ggplot(df, aes(x = Sample, y = Abundance, fill = Rank2)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Relative Abundance at Phylum Level",
       x = "Sample",
       y = "Relative Abundance") +
  theme_minimal()

# Genus level

# Species level (top 10 or 15)

# TO-DO make graph prettier and also remove p_, s_ etc.

### === 3 | ALPHA DIVERSITY MEASURES ========
plot_richness(physeq)
