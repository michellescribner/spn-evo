### Allele Frequency Plots for S. pneumoniae populations

library(tidyverse)
library(viridis)
library(randomcoloR)

# Begin with filtered data frame of SNPs 
# In prior step, SNPs fixed in the ancestral strain, <10% in a sample, and <100% (15 passages) or <200% (29 passages) cumulative frequency per lineage were removed
snps <- read.csv(file="~/Documents/pitt/streppneumo/post_breseq_filtering/output/snps_after_filtering.csv", header=TRUE)
nrow(snps) 

# Paste gene and annotation data
snps$gene_annot <- paste(snps$gene, snps$annotation, sep = " ")
snps$gene_annot <- gsub("\\s*\\([^\\)]+\\)","",as.character(snps$gene_annot))
# Replace locus tag with gene name for a key gene
snps$gene <- gsub("SP_1530", "murE", snps$gene)

# How many unique mutations?
length(unique(snps$gene_annot)) 

# Assign colors only to mutations in genes that were mutated in many lineages
parallel <- snps %>%
  group_by(gene) %>%
  dplyr::count(lineage) %>%
  dplyr::count(gene)
par <- parallel %>% 
  filter(n > 6)
par_genes <- par$gene
length(par_genes) 

# Add a column that indicates if a mutation is in one of these genes or "other"
snps <- snps %>%
  mutate(grp = if_else(gene %in% par_genes, gene, "other"))

# Build named color vector matching factor levels for plotting
grp_levels <- c(sort(unique(par_genes)), "other")
n_colors <- length(grp_levels)
pal <- distinctColorPalette(k = n_colors, altCol = FALSE, runTsne = FALSE)

# name the palette vector with the factor levels
names(pal) <- grp_levels

# Pivot data so that passages where mutations weren't detected are 0 freq
snps_cast <- snps %>%
  pivot_wider(id_cols = c(lineage,drug,immune,pop,gene,gene_annot,grp), 
              names_from = gen, values_from = freq, values_fn=mean, values_fill=0)

# Filter mutations that do not reach at least 100% frequency per lineage,
# this should now be completely redundant as this was introduced into upstream processing
lineages <- unique(snps$lineage)
output <- data.frame()
for (i in 1:length(lineages)) {
  n <- snps_cast %>%
    filter(lineage == lineages[i]) %>%
    mutate(sums = rowSums(select(., -c(lineage,drug,immune,pop,gene,gene_annot,grp)))) %>%
    filter(sums > 100) %>%
    select(-c(sums)) %>%
    pivot_longer(cols = -c(lineage,drug,immune,pop,gene,gene_annot,grp), names_to = "passage", values_to = "freq")
  output <- rbind(output,n)
}
output$passage <- as.numeric(output$passage)
output$sample <- paste(output$lineage, output$passage, sep = "G")

#----- Mask frequency of low coverage samples ----------
# This must be performed after pivot_longer because pivot_longer will falsely
# set frequency to zero

coverage <- read.csv("~/Documents/pitt/streppneumo/post_breseq_filtering/raw_breseq_output/coverage_results.csv",header=TRUE)

# clean sample name data
clean_samplenames <- function(snps) {
  snps$sample <- gsub("CM", "CeM", snps$sample)
  snps$sample <- gsub("CN", "CeN", snps$sample)
  snps$sample <- gsub("MM", "MeM", snps$sample)
  snps$sample <- gsub("MN", "MeN", snps$sample)
  snps$sample <- gsub("ND", "Nd", snps$sample)
  snps$sample <- gsub("A", "Az", snps$sample)
  snps$sample <- gsub("V", "Va", snps$sample)
  snps$sample <- gsub("P", "Pe", snps$sample)
  snps$sample <- gsub("CI", "Ci", snps$sample)
  snps$sample <- gsub("I", "Im", snps$sample)
  snps$sample <- gsub("R", "Ri", snps$sample)
  # separate sample names into columns 
  snps <- separate(snps, sample, into = c("sample", "garbage"), sep = c(9), remove = FALSE)
  snps <- separate(snps, sample, into = c("drug", "immune", "pop", "gen"), sep = c(2, 4, 6), remove = FALSE)
  # remove G from generation name
  snps$gen <- gsub("G", "", snps$gen) 
  return(snps)
}
coverage <- clean_samplenames(coverage)

high_coverage <- coverage %>%
  filter(average_cov >= 30) %>%
  mutate(across(everything(), gsub, pattern = ".gd", replacement = "")) %>%
  mutate(across(everything(), gsub, pattern = "\\.", replacement = ""))

output <- output[(output$sample %in% high_coverage$sample), ]
output$sample <- NULL

#----------- Plot frequencies ----------------------

# Function to parse labels so non-SP_ names are italicized
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

# Prepare data for plotting 

setwd("~/Documents/pitt/streppneumo/allelefreq_all/facetbydrug_output_100cutoff")

plot_data <- output %>%
  filter(passage <= 29) %>%
  mutate(condition = paste(drug, lineage, sep = "_"))

# Ensure grp exists
if (!"grp" %in% colnames(plot_data)) stop("column 'grp' not found in output")

# Create a scale that parses labels (expressions) so italic/plain are rendered
colScale_parsed <- scale_colour_manual(name = "grp",
                                       values = pal,
                                       labels = italic_if_not_SP)

# Define drug column order (10 columns)
desired_order <- c("Ci","Az","Im","Va","Ri","Pe","Me","Ce","Le","Li")
drugs <- intersect(desired_order, sort(unique(plot_data$drug)))
if (length(drugs) == 0) stop("No matching drugs found in plot_data$drug")

n_drugs <- length(drugs)

# Build list of sorted lineages for each drug
lineages_by_drug <- lapply(drugs, function(d) {
  plot_data %>% filter(drug == d) %>% distinct(lineage) %>% arrange(lineage) %>% pull(lineage)
})
names(lineages_by_drug) <- drugs

# Column-major ordering (so each drug's lineages go down the same column)
max_rows <- max(sapply(lineages_by_drug, length))
cond_levels <- character()
for (r in seq_len(max_rows)) {
  for (d in drugs) {
    lin_vec <- lineages_by_drug[[d]]
    if (length(lin_vec) >= r) {
      cond_levels <- c(cond_levels, paste(d, lin_vec[r], sep = "_"))
    }
  }
}
cond_levels <- cond_levels[cond_levels %in% unique(plot_data$condition)]
plot_data$condition <- factor(plot_data$condition, levels = cond_levels)

# Main plot ---------------------------------------------------------------

lineage_labeller <- labeller(condition = function(x) sapply(strsplit(x, "_"), `[`, 2))

main_plot <- ggplot(plot_data,
                    aes(x = passage, y = freq, group = gene_annot, color = grp)) +
  geom_point(size = 0.45, alpha = 0.9) +
  geom_line(size = 0.5, alpha = 0.9) +
  colScale_parsed +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  scale_y_continuous(limits = c(0, 100)) +
  xlab("Passage") + ylab("Allele Frequency") +
  guides(col = guide_legend(nrow = 5, byrow = TRUE)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.text = element_text(size = rel(0.6), hjust = 0),
    legend.key.width = unit(0.25, "cm"),
    legend.title = element_blank(),
    strip.text = element_text(size = 6),     
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 8),
    plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
    strip.background = element_blank()
  ) +
  facet_wrap(~ condition, ncol = n_drugs, scales = "free_x", labeller = lineage_labeller)

# Save ---------------------------------------------------------------------
n_conditions <- length(cond_levels)
n_rows <- ceiling(n_conditions / n_drugs)
panel_width_in  <- 0.8   
panel_height_in <- 1.3   

out_width  <- max(8, n_drugs * panel_width_in)
out_height <- max(6, n_rows * panel_height_in)

ggsave(filename = "FigureS2.pdf",
       plot = main_plot,
       device = "pdf",
       width = out_width,
       height = out_height,
       units = "in",
       limitsize = FALSE)

#----------- Plot by Drug-----------------------------

# Plot allele frequency faceted by drug, with only key genes shown

# setwd("~/Documents/pitt/streppneumo/allelefreq_all/facetbydrug_output_100cutoff")
# 
# drugs1 <- c("Ci","Az", "Im", "Va", "Ri", "Pe", "Me")
# for (i in 1:length(drugs1)) {
#   plot <- output %>% 
#     filter(drug == drugs1[i]) %>%
#     filter(passage < 15) %>%
#     ggplot(mapping = aes(x=passage, y=freq, group=gene_annot, color=grp)) + 
#     geom_point(size =0.5) +
#     geom_line(size = 0.75) +
#     guides(col=guide_legend(nrow=4)) +
#     scale_x_continuous(breaks = c(0,5,10,15), limits = c(0,15))+
#     scale_y_continuous(limits = c(0,100)) + 
#     ylab("Frequency") +
#     xlab("Passage") + 
#     colScale +
#     theme_bw() +
#     theme(legend.position = "bottom",
#           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#           panel.background = element_rect(fill = "white"),
#           legend.text=element_text(size=rel(0.5)),legend.key.width = unit(0.25, 'cm'),
#           legend.title = element_blank()) +
#     facet_wrap(vars(lineage)) 
#   ggsave(plot, file = paste(drugs1[i], "_facet.pdf", sep =""), device = "pdf", width = 5, height = 5)
# }

# # Plot longer time
# drugs2 <- c("Ce", "Le", "Li")
# for (i in 1:length(drugs2)) {
#   plot <- output %>% 
#     filter(drug == drugs2[i]) %>%
#     filter(passage < 30) %>%
#     ggplot(mapping = aes(x=passage, y=freq, group=gene_annot, color=grp)) + 
#     geom_point(size =0.5) +
#     geom_line(size = 0.75) +
#     scale_x_continuous(breaks = c(0,10,20,30), limits = c(0,29))+
#     scale_y_continuous(limits = c(0,100)) + 
#     colScale+
#     ylab("Frequency") +
#     xlab("Passage") + 
#     theme_bw() +
#     theme(legend.position = "bottom",
#           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#           panel.background = element_rect(fill = "white"),
#           legend.text=element_text(size=rel(0.5)),legend.key.width = unit(0.25, 'cm'),
#           legend.title = element_blank()) +
#     facet_wrap(vars(lineage)) 
#   ggsave(plot, file = paste(drugs2[i], "_facet.pdf", sep=""), device = "pdf", width = 5, height = 5)
# }

# Plot individually
# 
# setwd("~/Documents/pitt/streppneumo/allelefreq_all/singlelineage_output")
# 
# # Ce, Le, Li
# muller <- function(snps, ab, im, rep) {
#   snps <- snps %>%
#     filter(drug == ab) %>%
#     filter(immune == im) %>%
#     filter(pop == rep) %>%
#     select(!c(lineage,drug,immune,pop,gene,grp)) %>%
#     mutate(sums = rowSums(select(., -gene_annot))) %>%
#     filter(sums > 100) %>%
#     select(!c(sums))
#   write.csv(snps, file = paste(ab, im, rep, "mullerinput.csv", sep="")) 
#   
#   snps <- snps %>%
#     pivot_longer(cols = !gene_annot, names_to = "passage", values_to = "freq")
#   snps$passage <- as.numeric(snps$passage)
#   
#   remove_low <- high_coverage %>%
#     filter(drug == ab) %>%
#     filter(immune == im) %>%
#     filter(pop == rep) 
#   snps <- snps[(snps$passage %in% remove_low$gen),]
#   
#   snps_plot <- ggplot(snps, aes(x=passage, y=freq, color=gene_annot))+
#     geom_point(size=0.5) + 
#     geom_line(size=1) +
#     scale_x_continuous(breaks = c(0,5,10,15,20,25,29), limits = c(0,29))+
#     scale_y_continuous(limits = c(0,100)) + 
#     scale_color_viridis_d()+
#     ylab("Allele Frequency") +
#     xlab("Passage") + 
#     theme(axis.title.x = element_text(size=10),axis.text.x = element_text(angle=0, colour = "black", vjust=1, hjust = 1, size=10), 
#           axis.text.y = element_text(colour = "black",size=10),axis.title.y = element_text(size=10), 
#           plot.title = element_text(face="bold",size = 10, hjust = 0.5), 
#           legend.position="none",strip.text.x = element_text(size=10), 
#           strip.text.y = element_text(size=10),strip.background = element_rect(colour="black"),
#           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"))
#   pdf(paste(ab,im,rep, ".pdf", sep = ""), width= 6, height = 6, useDingbats=F) 
#   print(snps_plot) 
#   dev.off()
# }
# 
# muller(snps_cast, "Ce", "M0", "L1")
# muller(snps_cast, "Ce", "M0", "L2")
# muller(snps_cast, "Ce", "M0", "L3")
# muller(snps_cast, "Ce", "N0", "L1")
# muller(snps_cast, "Ce", "N0", "L2")
# muller(snps_cast, "Ce", "N0", "L3")
# muller(snps_cast, "Ce", "Nd", "L1")
# muller(snps_cast, "Ce", "Nd", "L2")
# muller(snps_cast, "Ce", "Nd", "L3")

# muller(snps_cast, "Le", "M0", "L1")
# muller(snps_cast, "Le", "M0", "L2")
# muller(snps_cast, "Le", "M0", "L3")
# muller(snps_cast, "Le", "N0", "L1")
# muller(snps_cast, "Le", "N0", "L2")
# muller(snps_cast, "Le", "N0", "L3")
# muller(snps_cast, "Le", "Nd", "L1")
# muller(snps_cast, "Le", "Nd", "L2")
# muller(snps_cast, "Le", "Nd", "L3")

# muller(snps_cast, "Li", "M0", "L1")
# muller(snps_cast, "Li", "M0", "L2")
# muller(snps_cast, "Li", "M0", "L3")
# muller(snps_cast, "Li", "N0", "L1")
# muller(snps_cast, "Li", "N0", "L2")
# muller(snps_cast, "Li", "N0", "L3")
# muller(snps_cast, "Li", "Nd", "L1")
# muller(snps_cast, "Li", "Nd", "L2")
# muller(snps_cast, "Li", "Nd", "L3")
# 
# muller <- function(snps, ab, im, rep) {
#   snps <- snps %>%
#     filter(drug == ab) %>%
#     filter(immune == im) %>%
#     filter(pop == rep) %>%
#     select(c(gene_annot,'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14')) %>%
#     mutate(sums = rowSums(select(., -gene_annot))) %>%
#     filter(sums > 100) %>%
#     select(!c(sums)) 
#   write.csv(snps, file = paste(ab, im, rep, "mullerinput.csv", sep="")) 
#   
#   snps <- snps %>%
#     pivot_longer(cols = !gene_annot, names_to = "passage", values_to = "freq")
#   snps$passage <- as.numeric(snps$passage)
#   
#   remove_low <- high_coverage %>%
#     filter(drug == ab) %>%
#     filter(immune == im) %>%
#     filter(pop == rep) 
#   snps <- snps[(snps$passage %in% remove_low$gen),]
#   
#   snps_plot <- ggplot(snps, aes(x=passage, y=freq, color=gene_annot))+
#     geom_point(size=0.5) + 
#     geom_line(size=1) +
#     guides(col=guide_legend(ncol=2)) +
#     ylab("Allele Frequency") +
#     xlab("Passage") + 
#     scale_x_continuous(breaks = c(0,5,10,15), limits = c(0,15))+
#     scale_y_continuous(limits = c(0,100)) + 
#     scale_color_viridis_d()+
#     theme(axis.title.x = element_text(size=10),axis.text.x = element_text(angle=0, colour = "black", vjust=1, hjust = 1, size=10), 
#           axis.text.y = element_text(colour = "black",size=10),axis.title.y = element_text(size=10), 
#           plot.title = element_text(face="bold",size = 10, hjust = 0.5), 
#           legend.position="bottom",legend.key=element_rect(fill="white"),strip.text.x = element_text(size=10), 
#           strip.text.y = element_text(size=10),strip.background = element_rect(colour="black"),
#           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"))
#   pdf(paste(ab,im,rep, ".pdf", sep = ""), width= 6, height = 6, useDingbats=F) 
#   print(snps_plot) 
#   dev.off()
# }
# 
# muller(snps_cast, "Va", "M0", "L1")
# muller(snps_cast, "Va", "M0", "L2")
# muller(snps_cast, "Va", "M0", "L3")
# muller(snps_cast, "Va", "N0", "L1")
# muller(snps_cast, "Va", "N0", "L2")
# muller(snps_cast, "Va", "N0", "L3")
# muller(snps_cast, "Va", "Nd", "L1")
# muller(snps_cast, "Va", "Nd", "L2")
# muller(snps_cast, "Va", "Nd", "L3")
# 
# muller(snps_cast, "Pe", "M0", "L1")
# muller(snps_cast, "Pe", "M0", "L2")
# muller(snps_cast, "Pe", "M0", "L3")
# muller(snps_cast, "Pe", "N0", "L1")
# muller(snps_cast, "Pe", "N0", "L2")
# muller(snps_cast, "Pe", "N0", "L3")
# muller(snps_cast, "Pe", "Nd", "L1")
# muller(snps_cast, "Pe", "Nd", "L2")
# muller(snps_cast, "Pe", "Nd", "L3")
# 
# muller(snps_cast, "Me", "M0", "L1")
# muller(snps_cast, "Me", "M0", "L2")
# muller(snps_cast, "Me", "M0", "L3")
# muller(snps_cast, "Me", "N0", "L1")
# muller(snps_cast, "Me", "N0", "L2")
# muller(snps_cast, "Me", "N0", "L3")
# muller(snps_cast, "Me", "Nd", "L1")
# muller(snps_cast, "Me", "Nd", "L2")
# muller(snps_cast, "Me", "Nd", "L3")
# 
# muller(snps_cast, "Ri", "M0", "L1")
# muller(snps_cast, "Ri", "M0", "L2")
# muller(snps_cast, "Ri", "M0", "L3")
# muller(snps_cast, "Ri", "N0", "L1")
# muller(snps_cast, "Ri", "N0", "L2")
# muller(snps_cast, "Ri", "N0", "L3")
# muller(snps_cast, "Ri", "Nd", "L1")
# muller(snps_cast, "Ri", "Nd", "L2")
# muller(snps_cast, "Ri", "Nd", "L3")
# 
# muller(snps_cast, "Az", "M0", "L1")
# muller(snps_cast, "Az", "M0", "L2")
# muller(snps_cast, "Az", "M0", "L3")
# muller(snps_cast, "Az", "N0", "L1")
# muller(snps_cast, "Az", "N0", "L2")
# muller(snps_cast, "Az", "N0", "L3")
# muller(snps_cast, "Az", "Nd", "L1")
# muller(snps_cast, "Az", "Nd", "L2")
# muller(snps_cast, "Az", "Nd", "L3")
# 
# muller(snps_cast, "Ci", "M0", "L1")
# muller(snps_cast, "Ci", "M0", "L2")
# muller(snps_cast, "Ci", "M0", "L3")
# muller(snps_cast, "Ci", "N0", "L1")
# muller(snps_cast, "Ci", "N0", "L2")
# muller(snps_cast, "Ci", "N0", "L3")
# muller(snps_cast, "Ci", "Nd", "L1")
# muller(snps_cast, "Ci", "Nd", "L2")
# muller(snps_cast, "Ci", "Nd", "L3")
# 
# muller(snps_cast, "Im", "M0", "L1")
# muller(snps_cast, "Im", "M0", "L2")
# muller(snps_cast, "Im", "M0", "L3")
# muller(snps_cast, "Im", "N0", "L1")
# muller(snps_cast, "Im", "N0", "L2")
# muller(snps_cast, "Im", "N0", "L3")
# muller(snps_cast, "Im", "Nd", "L1")
# muller(snps_cast, "Im", "Nd", "L2")
# muller(snps_cast, "Im", "Nd", "L3")

