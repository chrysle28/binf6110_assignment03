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
library(dplyr)
library(sunburstR)
library(d3r)
library(ANCOMBC)
library(microbiome)
library(extrafont)
library(RColorBrewer)

### === 1 | DATA INPUT AND CLEANING ========
# Import BIOM table (threshold = 0) from Kraken2/Bracken classification
biom_data <- read_biom("../table.biom")
physeq_0 <- import_biom(biom_data)
physeq_0

# Import BIOM table (threshold = 10) from Kraken2/Bracken classification
biom_data <- read_biom("../table_10.biom")
physeq_10 <- import_biom(biom_data)
physeq_10

# Make metadata table
metadata <- data.frame(
  row.names = sample_names(physeq_0),
  Diet = c("Omnivore", "Omnivore", "Omnivore", "Omnivore", "Omnivore", "Vegan", "Vegan", "Vegan", "Vegan", "Vegan"),
  SRR = c("SRR8146935", "SRR8146936", "SRR8146938", "SRR8146969", "SRR8146970", "SRR8146963", "SRR8146968", "SRR8146973", "SRR8146977", "SRR8146978"))

# Function to format physeq objects
format_physeq <- function(physeq) {

  sample_data(physeq) <- metadata # Add metadata
  sample_names(physeq) <- gsub("_bracken_species", "", sample_names(physeq)) # Clean sample names
  tax_table(physeq) <- gsub("^[a-z]__", "", tax_table(physeq)) # Clean taxa names
  sample_names(physeq) <- paste0(sample_names(physeq), " (", metadata$SRR, ")") # Add SRR accessions
  colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") # Add taxa classifications
  tax_table(physeq)[, "Species"] <- paste(tax_table(physeq)[, "Genus"],
                                          tax_table(physeq)[, "Species"]) # Add genus to species names
  otu <- as.data.frame(t(otu_table(physeq))) # Make transposed OTU table
  return(list(physeq = physeq, otu = otu))
}

# Unfiltered abundances for alpha diversity measures
format_0 <- format_physeq(physeq_0)
physeq_0 <- format_0$physeq
otu_0 <- format_0$otu

# Filtered abundances for beta diversity and differential abundance
format_10 <- format_physeq(physeq_10)
physeq_10 <- format_10$physeq
otu_10 <- format_10$otu

# Rarefaction curve to check for sufficient sequencing depth
rare_curve <- rarecurve(otu_0, step = 1000, sample = min(rowSums(otu_0)), label = F)

rare_df <- lapply(seq_along(rare_curve), function(i){
  data.frame(Sample = rownames(otu_0)[i],
             Depth = attr(rare_curve[[i]], "Subsample"),
             Richness = rare_curve[[i]])
})
rare_df <- bind_rows(rare_df)
rare_df$Diet <- data.frame(sample_data(physeq_0))[rare_df$Sample, 1]

ggplot(rare_df, aes(x = Depth, y = Richness, group = Sample, color = Diet)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  labs(title = "Rarefaction Curves of Omnivore and Vegan Samples",
       x = "Sequencing Depth",
       y = "Observed Richness") +
  theme_minimal() +
  theme(text = element_text(family = "Open Sans"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 15))) +
  scale_color_manual(values = c("Omnivore" = "#dd1a24", "Vegan" = "#3e8f27"))

### === 2 | TAXONOMIC ABUNDANCE ========
# Convert to relative abundance
physeq_rel <- transform_sample_counts(physeq_10, function(x) x / sum(x))

# Genus level (top 10)
physeq_gen <- tax_glom(physeq_rel, taxrank = "Genus")
top_gen <- names(sort(taxa_sums(physeq_gen), decreasing = T)[1:10])
physeq_top_gen <- prune_taxa(top_gen, physeq_gen)
physeq_top_gen <- psmelt(physeq_top_gen)

ggplot(physeq_top_gen, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Diet, scales = "free_x") +
  labs(title = "Relative Abundance at Genus Level (Top 10)",
       x = "Sample",
       y = "Relative Abundance") +
  theme_minimal() +
  theme(text = element_text(family = "Open Sans"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 11))

# Species level (top 10)
physeq_spe <- tax_glom(physeq_rel, taxrank = "Species")
top_spe <- names(sort(taxa_sums(physeq_spe), decreasing = T)[1:10])
physeq_top_spe <- prune_taxa(top_spe, physeq_spe) %>%
  psmelt()

ggplot(physeq_top_spe, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Diet, scales = "free_x") +
  labs(title = "Relative Abundance at Species Level (Top 10)",
       x = "Sample",
       y = "Relative Abundance") +
  theme_minimal() +
  theme(text = element_text(family = "Open Sans"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 11)) +
  scale_fill_brewer(palette = "Paired")

# Creating sunburst plots for relative abundance in all taxa
physeq_rel <- psmelt(physeq_rel) %>%
  filter(Diet %in% c("Omnivore", "Vegan")) %>%
  mutate(across(c(Phylum, Class, Order, Family, Genus, Species),
                ~ifelse(is.na(.) | . == "", "Unclassified", .)))

# Omnivore sunburst
omni_tree <- d3_nest(
  physeq_rel %>%
    filter(Diet == "Omnivore") %>%
    group_by(Phylum, Class, Order, Family, Genus) %>%
    summarise(size = mean(Abundance), .groups = "drop"),
  value_cols = "size"
)
omni_tree <- gsub('"name":"root"', '"name":"Omnivore"', omni_tree)

omni_sb  <- sund2b(omni_tree, showLabels = TRUE)
omni_sb
htmlwidgets::saveWidget(omni_sb, "omni_sb.html")

# Vegan sunburst
vegan_tree <- d3_nest(
  physeq_rel %>%
    filter(Diet == "Vegan") %>%
    group_by(Phylum, Class, Order, Family, Genus) %>%
    summarise(size = mean(Abundance), .groups = "drop"),
  value_cols = "size"
)
vegan_tree <- gsub('"name":"root"', '"name":"Vegan"', vegan_tree)

vegan_sb  <- sund2b(vegan_tree, showLabels = TRUE)
vegan_sb
htmlwidgets::saveWidget(vegan_sb, "vegan_sb.html")

### === 3 | ALPHA DIVERSITY MEASURES ========
# Unfiltered abundances used to provide more robust estimates of richness
plot_richness(physeq_0, x = "Diet", measures = c("Observed", "Chao1", "Simpson", "Shannon", "Fisher")) +
  geom_point(aes(color = Diet), size = 2) +
  geom_boxplot(aes(fill = Diet), alpha = 0.6) +
  labs(title = "Alpha Diversity Measures by Diet") +
  theme_minimal() +
  theme(text = element_text(family = "Open Sans"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 15)),
        strip.text.x = element_text(size = 10),
        panel.spacing.x = unit(1.5, "cm")) +
  scale_fill_manual(values = c("Omnivore" = "#dd1a24", "Vegan" = "#3e8f27")) + 
  scale_color_manual(values = c("Omnivore" = "#dd1a24", "Vegan" = "#3e8f27"))

# Wilcoxon tests for alpha diversity measures
alpha_d <- estimate_richness(physeq_0, measures = c("Observed", "Chao1", "Simpson", "Shannon", "Fisher"))
alpha_d$Diet <- sample_data(physeq_0)$Diet

alpha_measures <- c("Observed", "Chao1", "Simpson", "Shannon", "Fisher")

alpha_results <- lapply(alpha_measures, function(i) {
  test <- wilcox.test(get(i) ~ Diet, data = alpha_d)
  data.frame(Measure = i,
             W = test$statistic,
             p_val = test$p.value)
})

alpha_df <- do.call(rbind, alpha_results)
alpha_df
# none of the measures differ significantly between diet

### === 4 | BETA DIVERSITY MEASURES ========
# Relative abundance of physeq_10
physeq_rel_10 <- transform_sample_counts(physeq_10, function(x) x / sum(x))

# Bray-Curtis PCoA
pcoa_bray <- ordinate(physeq_rel_10, method = "PCoA", distance = "bray")
plot_ordination(physeq_rel_10, pcoa_bray, color = "Diet", shape = "Diet") +
  labs(title = "Bray-Curtis PCoA of Omnivores vs Vegans") +
  geom_point(size = 5) +
  theme_minimal() +
  theme(text = element_text(family = "Open Sans"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 15)),
        strip.text.x = element_text(size = 10)) +
  scale_color_manual(values = c("Omnivore" = "#dd1a24", "Vegan" = "#3e8f27"))

# Bray-Curtis NMDS
nmds_bray <- ordinate(physeq_rel_10, method = "NMDS", distance = "bray")
plot_ordination(physeq_rel_10, nmds_bray, color = "Diet", shape = "Diet") +
  labs(title = "Bray-Curtis NMDS of Omnivores vs Vegans") +
  geom_point(size = 5) +
  theme_minimal() +
  theme(text = element_text(family = "Open Sans"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 15)),
        strip.text.x = element_text(size = 10)) +
  scale_color_manual(values = c("Omnivore" = "#dd1a24", "Vegan" = "#3e8f27"))

# PERMANOVA test
set.seed(14)
adonis2(phyloseq::distance(physeq_rel_10, method = "bray") ~ Diet, data = metadata)
# non-significant result (p > 0.05)

### === 5 | DIFFERENTIAL ABUNDANCE ========
ancombc_spe <- tax_glom(physeq_10, taxrank = "Species")
taxa_names(ancombc_spe) <- tax_table(physeq_10)[, "Species"]

ancombc <- ancombc2(data = ancombc_spe, fix_formula = "Diet", rand_formula = NULL, p_adj_method = "BH", pseudo_sens = T, prv_cut = 0, lib_cut = 1000, s0_perc = 0.05, group = "Diet", struc_zero = T, neg_lb = T)

ancombc_sig <- subset(ancombc$res, q_DietVegan < 0.05)
ancombc_sig # ANCOM-BC2 returns 0 significant taxa

ancombc_top <- ancombc$res %>%
  head(20) %>%
  mutate(Change = ifelse(lfc_DietVegan > 0, "Higher in Vegan", "Higher in Omnivore"))

ggplot(ancombc_top, aes(x = lfc_DietVegan, y = reorder(taxon, lfc_DietVegan))) +
  geom_point(aes(color = Change), size = 3) +
  geom_errorbar(aes(xmin = lfc_DietVegan - se_DietVegan,
                    xmax = lfc_DietVegan + se_DietVegan), color = "darkgray") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  labs(title = "Differential Abundance of Taxa in Omnivores vs Vegans",
       x = "Log Fold Change",
       y = "Species") +
  theme_minimal() +
  theme(text = element_text(family = "Open Sans"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 15)),
        strip.text.x = element_text(size = 10)) +
  scale_color_manual(values = c("Higher in Omnivore" = "#dd1a24", "Higher in Vegan" = "#3e8f27"))

