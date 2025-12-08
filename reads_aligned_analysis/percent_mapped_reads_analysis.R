## Analysis of proportion of reads that align to S. pneumoniae genome

library(tidyverse)

setwd("~/Documents/pitt/streppneumo/reads_aligned_analysis")

coverage_results <- read.csv('~/Documents/pitt/streppneumo/post_breseq_filtering/raw_breseq_output/coverage_results.csv')
coverage_results$sample <- coverage_results$sample_name

clean_samplenames <- function(snps) {
  snps$sample <- gsub(".gd", "", snps$sample)
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
  # snps$sample <- gsub("G([1-9]{1})$", "G0\\1", as.character(snps$sample))
  # snps$sample <- gsub("G([1-9]{1})R", "G0\\1", as.character(snps$sample))
  
  # separate sample names into columns 
  snps <- separate(snps, sample, into = c("drug", "immune", "pop", "gen", "garbage"), sep = c(2, 4, 6, 9), remove = FALSE)
  
  # remove samples that aren't in the drug and immune treatment list
  drugs <- c("Ce", "Le", "Li", "Ci", "Az", "Im", "Va", "Pe", "Ri", "Me")
  immune <- c("M0", "N0", "Nd")
  snps <- snps[ (snps$drug %in% drugs), ]
  snps <- snps[ (snps$immune %in% immune), ]
  
  # create new columns by pasting other columns
  snps$treat <- paste(snps$drug, snps$immune, sep = "")
  snps$lineage <- paste(snps$drug, snps$immune, snps$pop, sep = "")
  snps$sample <- paste(snps$drug, snps$immune, snps$pop, snps$gen, sep = "")
  snps$treat <- paste(snps$drug, snps$immune, sep = "")
  
  # remove G from generation name
  snps$gen <- gsub("G", "", snps$gen) 
  
  return(snps)
}

df <- clean_samplenames(coverage_results)

df$drug <- recode(df$drug, "Az" = "Azithromycin", "Ci" = "Ciprofloxacin", "Im"="Imipenem",
                    "Ce"="Cefepime","Le"="Levofloxacin", "Li"="Linezolid",
                    "Pe"="Penicillin","Me"="Meropenem", "Ri"="Rifampicin","Va"="Vancomycin")
df$immune <- recode(df$immune, "M0" = "Macrophage Depleted", "N0" = "Neutrophil Depleted", "Nd" = "Not Depleted")
df$pop <- recode(df$pop, "L1" = "Lineage 1", "L2" = "Lineage 2", "L3" = "Lineage 3")

df$gen <- as.integer(df$gen)

percent_mapped_bases_block1 <- df %>% 
  filter(drug %in% c("Cefepime", "Levofloxacin", "Linezolid")) %>%
  ggplot(mapping = aes(x=gen, y=percent_mapped_bases, color=pop)) + 
  geom_point(size =1) +
  ylab("Proportion of Bases Aligned") +
  xlab("Passage") + 
  labs(color = "Lineage") +
  scale_y_continuous(limits=c(0,1))+
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(vars(drug,immune)) 
ggsave("percent_mapped_bases_block1.pdf", plot = percent_mapped_bases_block1, device = "pdf", height =5, width =6)

percent_mapped_bases_block2 <- df %>% 
  filter(drug %in% c("Imipenem", "Ciprofloxacin", "Azithromycin")) %>%
  ggplot(mapping = aes(x=gen, y=percent_mapped_bases, color=pop)) + 
  geom_point(size =1) +
  ylab("Proportion of Bases Aligned") +
  xlab("Passage") + 
  labs(color = "Lineage") +
  scale_y_continuous(limits=c(0,1))+
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(vars(drug,immune))
ggsave("percent_mapped_bases_block2.pdf", plot = percent_mapped_bases_block2, device = "pdf", height =5, width =6)

percent_mapped_bases_block3 <- df %>% 
  filter(drug %in% c("Penicillin","Rifampicin","Vancomycin","Mereopenem")) %>%
  ggplot(mapping = aes(x=gen, y=percent_mapped_bases, color=pop)) + 
  geom_point(size =1) +
  ylab("Proportion of Bases Aligned") +
  labs(color = "Lineage") +
  xlab("Passage") + 
  scale_y_continuous(limits=c(0,1))+
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(vars(drug,immune),nrow=4)

ggsave("percent_mapped_bases_block3.pdf", plot = percent_mapped_bases_block3, device = "pdf", height =6, width =6)

###########################################3

# # Coverage 

cov <- df %>% filter(average_cov < 1000)

read_depth_block1 <- cov %>% 
  filter(drug %in% c("Cefepime", "Levofloxacin", "Linezolid")) %>%
  ggplot(mapping = aes(x=gen, y=average_cov, color=pop)) + 
  geom_point(size =1) +
  ylab("Read Depth") +
  xlab("Passage") + 
  labs(color = "Lineage") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(vars(drug,immune))
ggsave("read_depth_block1.pdf", plot = read_depth_block1, device = "pdf", height =5, width =6)


read_depth_block2 <- cov %>% 
  filter(drug %in% c("Imipenem", "Ciprofloxacin", "Azithromycin")) %>%
  ggplot(mapping = aes(x=gen, y=average_cov, color=pop)) + 
  geom_point(size =1) +
  ylab("Read Depth") +
  xlab("Passage") + 
  labs(color = "Lineage") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(vars(drug,immune))
ggsave("read_depth_block2.pdf", plot = read_depth_block2, device = "pdf", height =5, width =6)

read_depth_block3 <- cov %>% 
  filter(drug %in% c("Penicillin","Rifampicin","Vancomycin","Mereopenem")) %>%
  ggplot(mapping = aes(x=gen, y=average_cov, color=pop)) + 
  geom_point(size =1) +
  ylab("Read Depth") +
  xlab("Passage") + 
  labs(color = "Lineage") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(vars(drug,immune),nrow=4)
ggsave("read_depth_block3.pdf", plot = read_depth_block3, device = "pdf", height =6, width =6)

