# # # Identify genes with mutations associated with antibiotic and immune treatments in S. pneumoniae evolved populations
# This script also includes the cowplot to plot Figure 1

library(tidyverse)
library(cowplot)

setwd("~/spn-evo/treat_association/")

# Begin with snp df with snps fixed in the ancestral strain, <10%, and <100% or <200% cumulative frequency per lineage removed
snps <- read.csv("~/spn-evo/post_breseq_filtering/output/snps_after_filtering.csv",header=TRUE)

# subset only mutations that are in genes
snps$type <- "gene"
snps$type[grep("intergenic", snps$annotation)] <- "intergenic"
snps <- subset(snps, snps$type=="gene")

# Pivot the data so that passages where a mutation was not detected have frequency values of zero
w <- snps %>%
  pivot_wider(id_cols = c(description,gene), 
              names_from = lineage, values_from = freq, values_fn = list(freq = sum), 
              values_fill = list(freq = 0))
l <- w %>%
  pivot_longer(!c(description,gene), names_to = "lineage", values_to = "sum") %>%
  separate(lineage, c("drug", "immune", "rep"), sep = c(2,4))

# Create list of genes that were mutated in more than one lineage
im <- l %>%
  filter(sum > 0) %>%
  dplyr::count(gene) %>%
  filter(n > 1)
genelist <- as.character(im$gene) 


# # # Identify genes that mutated more frequently in certain treatments

# # Immune treatment

# if sum of mutation frequency for a gene is > 0, set sum to 1
l$sum[l$sum > 0] <- 1

result_fet <- length(genelist) 

for(i in 1:length(genelist)) {
  n <- l %>%
    filter(gene == genelist[i]) %>% 
    group_by(immune) %>% 
    summarize(mut = sum(sum))
  ct <- matrix(c(n$mut, 30-n$mut), nrow=3, 
               dimnames = list(c("Macrophage depleted", "Neutrophil depleted", "Not depleted"), 
                               c("Mutated", "Not Mutated")))
  # print(ct)
  r <- fisher.test(ct)
  result_fet[i] <- r$p.value
}

output <- as.data.frame(result_fet)

output <- output %>%
  rename(pvalue = result_fet) %>%
  mutate(pvalue = as.numeric(pvalue)) %>%
  mutate(transp = -log10(pvalue)) %>%
  mutate(gene = genelist) %>%
  arrange(desc(transp)) %>%
  mutate(num = 1:n()) %>%
  mutate(fdr = 0.05*num/nrow(output)) %>%
  mutate(sig = pvalue < fdr) 

# print table 
table <- l %>%
  pivot_wider(id_cols = c(description,gene), names_from = immune, values_from = sum, values_fn = list(sum = sum), values_fill = list(freq = 0))
table <- full_join(table, output) %>%
  arrange(desc(transp)) 
head(table)
# write.csv(table,file="pvalues_fet_withcontingency.csv")

df_long <- table %>%
  filter(sig == TRUE) %>%
  select(gene, M0, N0, Nd) %>%
  pivot_longer(cols = !gene, names_to = "condition")

italic_if_not_SP <- function(x) {
  labs <- vapply(x, FUN.VALUE = character(1), FUN = function(g) {
    if (is.na(g)) return("plain('NA')")
    # escape single quotes
    g_esc <- gsub("'", "\\\\'", g)
    if (grepl("SP_", g, fixed = TRUE)) {
      paste0("plain('", g_esc, "')")
    } else {
      paste0("italic('", g_esc, "')")
    }
  })
  parse(text = labs)
}

p_heatmap <- df_long %>%
  ggplot(aes(x = condition, y = gene, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#297b8e") +
  scale_y_discrete(labels = italic_if_not_SP) +
  geom_text(aes(label = value), color = "black", size = 3) +
  labs(x = "Immune Condition", y = "Gene", fill = "Count", title = "Number of Mutated Lineages") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 10,hjust = 0.5))

print(p_heatmap)

# Add total lineages per condition
# M0 (Macrophage depleted), Nd (Not depleted), and N0 (Neutrophil depleted) are condition labels
total_lineages <- c(M0 = 30, Nd = 30, N0 = 30)  

# Calculate proportions of mutated lineages per gene per condition
df_long_prop <- df_long %>%
  mutate(proportion = value / total_lineages[condition])  # Calculate proportion

# Create heatmap with proportions of lineages mutated
p_heatmap_prop <- df_long_prop %>%
  ggplot(aes(x = condition, y = gene, fill = proportion)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#297b8e", limits = c(0,1)) +  
  scale_y_discrete(labels = italic_if_not_SP) +  
  geom_text(aes(label = sprintf("%.2f", proportion)), color = "black", size = 3) +  # Show proportion as label
  labs(x = "Immune Condition", y = "Gene", fill = "Proportion", title = "Proportion of Lineages\nwith Mutations") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 10, hjust = 0.5))

# Print the heatmap
print(p_heatmap_prop)


# # # Perform Fisher's Exact Test for drug association

# Assign drug mechanisms of action
l$moa <- "CWSI"
l$moa[l$drug == "Az"] <- "PSI"
l$moa[l$drug == "Li"] <- "PSI"
l$moa[l$drug == "Ci"] <- "DSI"
l$moa[l$drug == "Le"] <- "DSI"
l$moa[l$drug == "Ri"] <- "RSI"

result_fetab <- length(genelist)

for(i in 1:length(genelist)) {
  n <- l %>%
    filter(gene == genelist[i]) %>% 
    group_by(moa) %>% 
    summarize(mut = sum(sum))
  ct <- matrix(c(n$mut, c(45-n$mut[1], 18-n$mut[2], 18-n$mut[3], 9-n$mut[4])), 
               nrow=4, dimnames = list(c("CWSI", "DSI", "PSI", "RSI"), c("Mutated", "Not Mutated")))
  print(genelist[i])
  print(ct)
  r <- fisher.test(ct)
  result_fetab[i] <- r$p.value
}

output <- as.data.frame(result_fetab)
output <- output %>%
  rename(pvalue = result_fetab) %>%
  mutate(pvalue = as.numeric(pvalue)) %>%
  mutate(transp = -log10(pvalue)) %>%
  mutate(gene = genelist) %>%
  arrange(desc(transp)) %>%
  mutate(num = 1:nrow(output)) %>%
  mutate(fdr = 0.05*num/nrow(output)) %>%
  mutate(sig = pvalue < fdr)

# Print table 
table <- l %>%
  pivot_wider(id_cols = c(description,gene), names_from = moa, values_from = sum, values_fn = list(sum = sum), values_fill = list(freq = 0))
table <- full_join(table, output) %>%
  arrange(desc(transp))
# write.csv(table,file="pvalues_fetab_withcontingency.csv")

df_long <- table %>%
  filter(sig == TRUE) %>%
  select(gene, CWSI, DSI, PSI, RSI) %>%
  pivot_longer(cols = !gene, names_to = "condition")

original_genes <- table$gene
df_long <- df_long %>% mutate(gene = factor(gene, levels = rev(original_genes)))

p_heatmap_a <- df_long %>%
  ggplot(aes(x = condition, y = gene, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#297b8e") +
  scale_y_discrete(labels = italic_if_not_SP) +
  geom_text(aes(label = value), color = "black", size = 3) +
  labs(x = "Antibiotic\nMechanism of Action", y = "Gene", fill = "Count") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10))

print(p_heatmap_a)

# Define total lineages per condition
# These represent the total number of lineages in each Mechanism of Action (MoA) group
total_lineages_moa <- c(CWSI = 45, DSI = 18, PSI = 18, RSI = 9)  

# Calculate proportions for each gene across MoAs
df_long_prop_moa <- df_long %>%
  mutate(proportion = value / total_lineages_moa[condition]) 

# Create heatmap with proportions as the metric
p_heatmap_prop_moa <- df_long_prop_moa %>%
  ggplot(aes(x = condition, y = gene, fill = proportion)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#297b8e") +
  scale_y_discrete(labels = italic_if_not_SP) +  
  geom_text(aes(label = sprintf("%.2f", proportion)), color = "black", size = 3) +  # Format proportions to 2 decimal places
  labs(x = "Antibiotic\nMechanism of Action", y = "Gene", fill = "Proportion") +
  theme_bw() +
  theme(legend.position = "none",  
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10))

# Print the heatmap
print(p_heatmap_prop_moa)

#------------ Plot Figure 1 -------------

# # Can only be used if other plots are stored as objects 

right <- plot_grid(p_heatmap_prop, p_heatmap_prop_moa, ncol = 1, align = "v", rel_heights = c(1, 1.6), labels = c("C", ""))
top <- plot_grid(p_jaccard_withlegend, right, ncol = 2, rel_widths = c(2.4,1), labels = c("A", ""))
bottom <- plot_grid(p, bottom_right, ncol = 2, rel_widths = c(1, 1.1), labels = c("B", ""))
combined <- plot_grid(top, bottom, ncol = 1, align = "v", rel_heights = c(2, 1.3))
combined

# save 
ggsave("Figure1.png", combined, width = 10, height = 11, dpi = 300, units = "in")



