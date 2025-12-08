# # # Identify genes with mutations associated with antibiotic and immune treatments in S. pneumoniae evolved populations
# This script also includes the cowplot to plot Figure 1

library(tidyverse)
library(cowplot)

setwd("~/Documents/pitt/streppneumo/treat_association/")
  
# Begin with snp df with snps fixed in the ancestral strain, <10%, and <100% or <200% cumulative frequency per lineage removed
snps <- read.csv("~/Documents/pitt/streppneumo/post_breseq_filtering/output/snps_after_filtering.csv",header=TRUE)

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

# genes_to_label <- output %>% filter(sig == TRUE) %>% select(gene)
# genes_to_label$gene <- with(genes_to_label, ifelse(gene == "SP_1583", "SP_1583", 
#                 paste0("*",gene,"*")))

# output %>%
#   ggplot(aes(x=num, y=transp, color = transp)) + 
#   geom_point(size=3) +
#   theme_bw() + 
#   xlab("Gene") + ylab('-log10P') #+
#   # geom_text_repel(aes(x=num, y=transp, label = gene), 
#   #                 data = output[output$gene %in% genes_to_label$gene, ]) +
#   # theme(legend.text = element_markdown())
# ggsave("pvalues_fet.pdf", device = "pdf", height =4, width =5)

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

# ggsave(
#   filename = "immune_association.png",
#   plot = p_heatmap,
#   device = "png",
#   dpi = 300,
#   height = 1.5,
#   width = 4
# )

# Print contingency tables for significant genes
# n <- l %>%
#   filter(gene == "SP_1583") %>% 
#   group_by(immune) %>% 
#   summarize(mut = sum(sum)) 
# ct <- t(matrix(c(n$mut, 30-n$mut), nrow=3, dimnames = list(c("M0", "N0", "Nd"), c("Mutated", "Not Mutated"))))
# print(ct)


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

output %>%
 ggplot(aes(x=num, y=transp, color = transp)) + geom_point(size=3) +
  theme_bw() + xlab("Gene") + ylab('-log10P') #+
  geom_text_repel(aes(x=num, y=transp, label = gene))
# ggsave("pvalues_fetab.pdf", device = "pdf", height =4, width =5)

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

# ggsave(
#   filename = "antibiotic_association.png",
#   plot = p_heatmap_a,
#   device = "png",
#   dpi = 300,
#   height = 2,
#   width = 4
# )

#------------ Plot Figure 1 -------------

# # Can only be used if other plots are stored as objects 

right <- plot_grid(p_heatmap, p_heatmap_a, ncol = 1, align = "v", rel_heights = c(1, 1.6), labels = c("C", ""))
top <- plot_grid(p_jaccard_withlegend, right, ncol = 2, rel_widths = c(2.4,1), labels = c("A", ""))
bottom <- plot_grid(p, bottom_right, ncol = 2, rel_widths = c(1, 1.1), labels = c("B", ""))
combined <- plot_grid(top, bottom, ncol = 1, align = "v", rel_heights = c(2, 1.3))
combined

# save 
ggsave("combined.png", combined, width = 10, height = 11, dpi = 300, units = "in")

