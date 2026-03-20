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
# Import BIOM tables (threshold = 0) from Kraken2/Bracken classification
biom_data <- read_biom("../bracken_unfilt_results/table_unfilt.biom")
physeq <- import_biom(biom_data)
physeq

# Attach metadata to phyloseq object
metadata <- data.frame(
  row.names = sample_names(physeq),
  diet = c("omn", "omn", "omn", "veg", "veg", "veg")
)
sample_data(physeq) <- metadata
sample_names(physeq) <- gsub("_bracken_species", "", sample_names(physeq))
tax_table(physeq) <- gsub("^[a-z]__", "", tax_table(physeq)) # clean taxa names
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
head(tax_table(physeq))

# Transpose table to make rows = samples
otu_table <- as.data.frame(t(otu_table(physeq)))
head(otu_table)
#rare_curve_0 <- rarecurve(otu_table, step = 1)

### === 2 | TAXONOMIC ABUNDANCE ========
# Convert to relative abundance
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# Filter low abundance taxa
# physeq_rel_filtered <- filter_taxa(physeq_rel, function(x) mean(x) > 0.01, prune = TRUE)

# Phylum level
physeq_phy <- tax_glom(physeq_rel, taxrank = "Phylum")
physeq_phy <- psmelt(physeq_phy)
ggplot(physeq_phy, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Relative Abundance at Phylum Level",
       x = "Sample",
       y = "Relative Abundance") +
  theme_minimal()

# Genus level
physeq_gen <- tax_glom(physeq_rel, taxrank = "Genus")

# Species level (top 10 or 15)
physeq_spe <- tax_glom(physeq_rel, taxrank = "Species")
top_spe <- names(sort(taxa_sums(physeq_spe), decreasing = T)[1:15])
physeq_top_spe <- prune_taxa(top_spe, physeq_spe)
physeq_top_spe <- psmelt(physeq_top_spe)
physeq_top_spe$Species <- paste(physeq_top_spe$Genus, physeq_top_spe$Species, sep = " ")

ggplot(physeq_top_spe, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ diet, scales = "free_x") +
  labs(title = "Top 15 Most Abundant Species",
       x = "Sample",
       y = "Relative Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


# TO-DO make graph prettier (change font, color palette)

### === 3 | ALPHA DIVERSITY MEASURES ========
plot_richness(physeq)
