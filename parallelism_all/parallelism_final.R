# # Identify genes mutated in parallel in evolved populations of S. pneumoniae

# This analysis is intended to identify the number of independent mutations that occurred in each gene
# Many mutations occur at the same nucleotide within multiple samples in this study but such mutations are counted only once here

library(tidyverse)
library(ggsci)
library(ggrepel)
library(scales)

# ---- I/O ----
snps <- readr::read_csv(
  "~/Documents/pitt/streppneumo/post_breseq_filtering/output/snps_after_filtering.csv")

# Distinct mutation (unique position+mutation etc.)
dmut <- snps %>%
  select(position, mutation, gene, annotation, description) %>%
  distinct()
message("Distinct mutations (total): ", nrow(dmut))

# Keep only genic mutations
dmut_genes <- dmut %>%
  filter(!str_detect(annotation, "intergenic"))
message("Distinct genic mutations: ", nrow(dmut_genes))

# ---- Count distinct mutations per gene (across dataset) ----
# Group by gene only (description grabbed as first non-NA to avoid splitting a gene by minor annotation differences)
mutations_per_gene_distinct <- dmut_genes %>%
  group_by(gene) %>%
  summarize(
    n_distinct_mutations = n(),
    description = coalesce(first(na.omit(description)), first(description)),
    .groups = "drop"
  ) %>%
  arrange(desc(n_distinct_mutations))

# ---- Read gene list and merge ----
genes <- readr::read_csv(
  "~/Documents/pitt/streppneumo/variantcalling/breseqv35/TIGR4vApr2021_gff3_split.csv")

# Substitute RefSeq locus tags with the old locus tags where possible,
# this was performed upstream for the snp df
genes <- genes %>% 
  mutate(old_locus_tag = coalesce(old_locus_tag,locus_tag))
# not the best approach but it works
for (i in 1:nrow(genes)) {
  if (genes$gene[i] == genes$locus_tag[i]) {
    genes$gene <- gsub(genes$gene[i], genes$old_locus_tag[i], genes$gene)
  }
}

# Join counts (keep all genes, fill missing counts with 0)
gene_hist <- genes %>%
  left_join(mutations_per_gene_distinct, by = "gene") %>%
  mutate(n_distinct_mutations = replace_na(n_distinct_mutations, 0))

# Create significance flag (>4 distinct mutations)
gene_hist <- gene_hist %>%
  mutate(sig = n_distinct_mutations > 4) %>%
  arrange(desc(n_distinct_mutations))

# Save merged table
readr::write_csv(gene_hist, "~/Documents/pitt/streppneumo/parallelism_all/parallel.csv")

# ---- Plot ----
genes_to_label <- gene_hist %>% filter(n_distinct_mutations > 4)

p <- gene_hist %>%
  filter(n_distinct_mutations > 0) %>%
  ggplot(aes(x = start, y = n_distinct_mutations)) +
  geom_point(aes(color = sig),size = 1.5) +
  geom_text_repel(
    data = genes_to_label,
    aes(x = start, y = n_distinct_mutations, label = old_locus_tag),
    size = 3,
    max.overlaps = 30
  ) +
  theme_bw() +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "#472a7a")) +
  scale_x_continuous(labels = scales::comma) +
  theme(legend.position = "none", axis.text = element_text(size = 10)) +
  xlab("Gene Start Position") + ylab("Number of Distinct Mutations")

print(p)

ggsave(
  filename = "parallelism.png",
  plot = p,
  path = "~/Documents/pitt/streppneumo/parallelism_all",
  device = "png",
  dpi = 300,
  height = 3,
  width = 5,
  units = "in"
)

