# # Initial data cleaning and filtering for S. pneumoniae populations  

library(tidyverse)
library(matrixStats)

setwd("/Users/mrs/Documents/pitt/streppneumo/post_breseq_filtering")

#---------- Summary of Upstream Processing ------------------------

# Prior to this step, sequencing reads were quality trimmed and filtered using Trimmomtic v0.36, (LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36) 
# and variants were called using Breseq v0.35.0 with respect to the S. pneumoniae reference genome NC_003028.3

#------------- Post Breseq Filtering ---------------------------

#--- Read data---
snps <- read.csv("raw_breseq_output/snps_combined.csv",header=TRUE)
nrow(snps) # 973838
mc <- read.csv("raw_breseq_output/mc_combined.csv",header=TRUE)
nje <- read.csv("raw_breseq_output/nje_combined.csv",header=TRUE)

#--- Remove samples with average read depth <30X ---
# coverage was calculated using the number of mapped bases from breseq output.gd file, divided by the reference genome length (2160842)
coverage <- read.csv("raw_breseq_output/coverage_results.csv",header=TRUE)
# filter coverage file for samples with >=30X coverage
high_coverage <- coverage %>%
  filter(average_cov >= 30) %>%
  mutate(sample_name = str_remove(sample_name, "\\.gd"))
# filter variant data for only samples with coverage >= 30X
snps <- snps %>%
  filter(sample %in% high_coverage$sample)
nrow(snps) # 966959

#--- Clean sample names ---
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
  snps$sample <- gsub("G([1-9]{1})$", "G0\\1", as.character(snps$sample))
  snps$sample <- gsub("G([1-9]{1})R", "G0\\1", as.character(snps$sample))
  
  # separate sample names into columns 
  snps <- separate(snps, sample, into = c("drug", "immune", "pop", "gen", "garbage"), sep = c(2, 4, 6, 9), remove = FALSE)
  
  # remove samples that aren't in the drug and immune treatment list
  drugs <- c("Ce", "Le", "Li", "Ci", "Az", "Im", "Va", "Pe", "Ri", "Me")
  immune <- c("M0", "N0", "Nd")
  snps <- snps[ (snps$drug %in% drugs), ]
  snps <- snps[ (snps$immune %in% immune), ]
  
  # remove sample with suspect numbers of mutations (>100,000)
  snps <- subset(snps, snps$sample != "CeM0L1G29")

  # create new columns by pasting other columns
  snps$treat <- paste(snps$drug, snps$immune, sep = "")
  snps$lineage <- paste(snps$drug, snps$immune, snps$pop, sep = "")
  snps$sample <- paste(snps$drug, snps$immune, snps$pop, snps$gen, sep = "")
  snps$treat <- paste(snps$drug, snps$immune, sep = "")

  # remove G from generation name
  snps$gen <- gsub("G", "", snps$gen) 
  
  return(snps)
}
snps <- clean_samplenames(snps)

#--- Clean position column ---
snps$position <- gsub(":1", "", as.character(snps$position))

#--- Clean gene names  ---
snps$gene <- gsub(" <-", "", snps$gene)
snps$gene <- gsub(" ->", "", snps$gene)
snps$gene <- gsub(" ", "", snps$gene)

#--- Replace RefSeq locus tags with "SP_" locus tags if possible ---
# Locus tag key below was made by parsing the gff3 file
genes <- read.csv("~/Documents/pitt/streppneumo/variantcalling/TIGR4vApr2021_gff3_split.csv", header = TRUE, row.names = 1)

# coalesce finds the first non-missing element, so old locus tag column 
# will be populated with the RefSeq locus tag if no old locus tag exists
genes <- genes %>% 
  mutate(old_locus_tag = coalesce(old_locus_tag,locus_tag))
for (i in 1:nrow(genes)) {
  if (genes$locus_tag[i] %in% snps$gene) {
    print(paste("Replacing ",genes$locus_tag[i]," with ",genes$old_locus_tag[i]))
    snps$gene <- gsub(genes$locus_tag[i], genes$old_locus_tag[i], snps$gene)
} else {
    print(paste(genes$locus_tag[i]," is not found in snps dataframe"))
}
}

# remove question marks from annotations
snps$mutation <- gsub("\\?","",as.character(snps$mutation))
snps$annotation <- gsub("\\?","",as.character(snps$annotation))

# Create new columns combining gene and annotation information
snps$gene_annot_desc <- paste(snps$gene, snps$annotation, snps$description, sep = "::")
snps$desc_gene <- paste(snps$description, snps$gene, snps$locustag, sep = "::")
snps$desc_annot <- paste(snps$description, snps$annotation, sep = " ")
snps$desc_annot <- gsub("\\s*\\([^\\)]+\\)","",as.character(snps$desc_annot))

# remove % sign from frequency column
snps$freq <- gsub( "%", "", as.character(snps$freq))
# Convert frequency values to numeric
snps$freq <- as.numeric(as.character(snps$freq))

# Are all samples accounted for?
snps_check <- snps %>% 
  pivot_wider(id_cols = gen, names_from = lineage, values_from = freq, values_fn = list(freq = sum))
write.csv(snps_check, file="output/sample_check.csv")

#--- Manually filter mutations based on dynamics --- 
# These mutations occur at intermediate frequencies with trajectories that are 
# not possible given the trajectories of other mutations
# first, save table of mutations
manually_removed <- snps %>% filter(desc_annot == "accessory Sec-dependent serine-rich glycoprotein adhesin coding" | 
                                      desc_annot == "accessory Sec-dependent serine-rich glycoprotein adhesin coding " | 
                                      desc_annot == "accessory Sec-dependent serine-rich glycoprotein adhesin T4094T " | 
                                      desc_annot == "accessory Sec-dependent serine-rich glycoprotein adhesin S3029S " | 
                                      desc_annot == "accessory Sec-dependent serine-rich glycoprotein adhesin A3626A " | 
                                      desc_annot == "accessory Sec-dependent serine-rich glycoprotein adhesin G245G " | 
                                      desc_annot == "accessory Sec-dependent serine-rich glycoprotein adhesin/IS630 transposase-related protein intergenic " | 
                                      desc_annot == "accessory Sec-dependent serine-rich glycoprotein adhesin A456A " | 
                                      desc_annot == "accessory Sec-dependent serine-rich glycoprotein adhesin S453S " | 
                                      desc_annot == "hypothetical protein G13G " | 
                                      desc_annot == "hypothetical protein G13C " | 
                                      desc_annot == "hypothetical protein T11A " | 
                                      desc_annot == "hypothetical protein R8K " |  
                                      desc_annot == "hypothetical protein V60A " | 
                                      desc_annot == "UDP-galactopyranose mutase/lysyl-tRNA synthetase intergenic" |
                                      desc_annot == "UDP-galactopyranose mutase pseudogene" |
                                      desc_annot == "A/G-specific adenine glycosylase/formate--tetrahydrofolate ligase intergenic" |
                                      desc_annot == "helicase-exonuclease AddAB subunit AddA/hypothetical protein intergenic") %>%
  select(position,mutation,gene,annotation,description,desc_annot) %>%
  distinct()

write.csv(manually_removed, file="output/snps_manually_removed.csv",row.names=FALSE)

# remove  
nrow(snps) # 296862 	
snps <- subset(snps, snps$desc_annot != "accessory Sec-dependent serine-rich glycoprotein adhesin coding")
snps <- subset(snps, snps$desc_annot != "accessory Sec-dependent serine-rich glycoprotein adhesin T4094T ")
snps <- subset(snps, snps$desc_annot != "accessory Sec-dependent serine-rich glycoprotein adhesin S3029S ")
snps <- subset(snps, snps$desc_annot != "accessory Sec-dependent serine-rich glycoprotein adhesin A3626A ")
snps <- subset(snps, snps$desc_annot != "accessory Sec-dependent serine-rich glycoprotein adhesin G245G ")
snps <- subset(snps, snps$desc_annot != "accessory Sec-dependent serine-rich glycoprotein adhesin/IS630 transposase-related protein intergenic ")
snps <- subset(snps, snps$desc_annot != "accessory Sec-dependent serine-rich glycoprotein adhesin A456A ")
snps <- subset(snps, snps$desc_annot != "accessory Sec-dependent serine-rich glycoprotein adhesin S453S ")
snps <- subset(snps, snps$desc_annot != "hypothetical protein G13G ")
snps <- subset(snps, snps$desc_annot != "hypothetical protein G13C ")
snps <- subset(snps, snps$desc_annot != "hypothetical protein T11A ")
snps <- subset(snps, snps$desc_annot != "hypothetical protein R8K ")
snps <- subset(snps, snps$desc_annot != "hypothetical protein V60A ")
snps <- subset(snps, snps$desc_annot != "UDP-galactopyranose mutase/lysyl-tRNA synthetase intergenic")
snps <- subset(snps, snps$desc_annot != "UDP-galactopyranose mutase pseudogene")
snps <- subset(snps, snps$desc_annot != "A/G-specific adenine glycosylase/formate--tetrahydrofolate ligase intergenic")
snps <- subset(snps, snps$desc_annot != "A/G-specific adenine glycosylase/formate--tetrahydrofolate ligase intergenic")
snps <- subset(snps, snps$desc_annot != "helicase-exonuclease AddAB subunit AddA/hypothetical protein intergenic")
nrow(snps) # 268848

# # Remove mutations that were fixed in the ancestor
# identify mutations in passage 0 with 100% frequency
snps_g0 <- subset(snps, snps$gen == 0 & snps$freq == 100)
snps_g0_table <- snps_g0 %>%
  group_by(gene_annot_desc) %>%
  summarize(n())
write.csv(snps_g0_table, file="output/snps_fixed_in_ancestor.csv",row.names=FALSE)
snps_g0_n <- snps_g0 %>% 
  group_by(sample) %>%
  summarize(n())
# Remove any mutation that is at 100% frequency in a passage 0 sample
snps_noref <- snps[ !(snps$position %in% snps_g0$position), ]
nrow(snps_noref) # 215164

# # Remove mutations at less than 10% frequency
# some mutations at low frequencies may be real, but many may be mapping errors
snps_10 <- subset(snps_noref, snps_noref$freq > 9.9)
nrow(snps_10) # 76807
length(unique(snps_10$gene_annot_desc)) # 12360
write.csv(snps_10, file="output/snps_10.csv",row.names=FALSE)

# Remove mutations that do not sum to at least 100% over all timepoints within a lineage
# except for lineages propagated for 30 passages, in which case the threshold is 200%
snps_lin_100 <- snps_10 %>%
  pivot_wider(id_cols = c(drug,immune,pop,evidence,lineage,treat,position,mutation,gene,annotation,description,gene_annot_desc), 
              names_from = gen, values_from = freq, values_fn = sum,values_fill=0) %>%
  mutate(sums = rowSums(select(., -c(drug,immune,pop,evidence,lineage,treat,position,mutation,gene,annotation,description,gene_annot_desc)))) %>%
  filter(case_when(drug=="Ce" ~ sums > 200,
                   drug=="Le" ~ sums > 200,
                   drug=="Li" ~ sums > 200, 
                   T ~ sums > 100)) %>%
  select(-sums)

# Filter the snps_10 df so that it only includes these mutations
snps_lin_100_only <- snps_lin_100 %>%
  pivot_longer(cols = !c(drug,immune,pop,evidence,lineage,treat,position,mutation,gene,annotation,description,gene_annot_desc), names_to = "gen", values_to = "freq") %>%
  filter(freq > 0)
nrow(snps_lin_100_only) # 13185

write.csv(snps_lin_100_only, file="output/snps_after_filtering.csv",row.names=FALSE)

# order lineages
factors <- read.csv("/Users/mrs/Documents/pitt/streppneumo/jaccard_index/treat_factors.csv",header=TRUE)
snps_lin_100_only <- left_join(snps_lin_100_only,factors) %>%
  arrange(order)

# Make a table showing mutation frequencies by treatment
snps_lin_100_cast <- snps_lin_100_only %>%
  pivot_wider(id_cols = c(position,mutation,gene,annotation,description), 
              names_from = treat, values_from = freq, values_fn = sum,values_fill=0) 
write.csv(snps_lin_100_cast, file="output/snps_mutations_by_treatment.csv",row.names=FALSE)

# Make a table showing mutation frequencies by lineage
snps_lin_100_lineage_cast <- snps_lin_100_only %>%
  pivot_wider(id_cols = c(position,mutation,gene,annotation,description), 
              names_from = lineage, values_from = freq, values_fn = sum,values_fill=0) 
write.csv(snps_lin_100_lineage_cast, file="output/snps_mutations_by_lineage.csv",row.names=FALSE)


#------------ Missing Coverage--------------------

nrow(mc) # 55225

mc <- clean_samplenames(mc)

#  remove mutations that were fixed in the ancestor
mc_g0 <- subset(mc, mc$gen == 0)
mc_g0_n <- mc_g0 %>% 
  group_by(sample) %>%
  summarize(n())
mc_noref <- mc[ !(mc$gene %in% mc_g0$gene), ]
nrow(mc_noref)# 50211

write.csv(mc_noref, file = "output/mc_noref.csv",row.names=FALSE)

# How many variants were called at each timepoint? 
mut_per_day <- mc_noref %>% 
  group_by(sample) %>% 
  summarize(n())
hist(mut_per_day$'n()')

mc_noref$size <- as.numeric(mc_noref$size)
mc_noref %>%
  arrange(desc(size))

#------------- New Junction Evidence ----------------

nrow(nje) # 180221

nje <- subset(nje, nje$sample != "sample")

# append even rows to the end of odd rows
evens <- seq(from = 2, to = nrow(nje), by = 2)
nje_even <- nje[evens, ]
nje_even <- nje_even[, c(3,4,9,10,11)]
colnames(nje_even) <- c("position2", "reads..cov2", "annotation2", "gene2", "product2")
nje_odd <- nje[-evens, ]
nje_done <- cbind(nje_odd, nje_even)

nje_done <- clean_samplenames(nje_done)

nje_done$gene_annot_desc <- paste(nje_done$gene, nje_done$annotation, nje_done$product, sep = "::")
# remove % sign from frequency column
nje_done$freq <- gsub( "%", "", as.character(nje_done$freq))
#convert frequency values to numeric
nje_done$freq <- as.numeric(as.character(nje_done$freq))
nrow(nje_done) # 86835

# remove mutations in the ancestor (not removed in final table aka nje_done)
nje_g0 <- subset(nje_done, nje_done$gen == 0)
nrow(nje_g0) #4360
nje_g0_table <- nje_g0 %>%
  group_by(gene_annot_desc) %>%
  summarize(n())
nje_g0_n <- nje_g0 %>% 
  group_by(sample) %>%
  summarize(n())
nje_noref <- nje_done[ !(nje_done$gene %in% nje_g0$gene), ]

# Remove mutations with less than 10% frequency
nje_done <- subset(nje_done, nje_done$freq > 9.9)

write.csv(nje_done, file = "output/nje_done.csv",row.names=FALSE)

# pivot data frame - organizing with each mutation as the rows and the frequency of that mutation on a given day as the columns
nje_cast <- nje_done %>%
  pivot_wider(id_cols = c(gene_annot_desc, lineage), names_from = gen, values_from = freq, values_fn = sum, values_fill = 0) %>%
  mutate(sums = rowSums(select(., -c(lineage,gene_annot_desc)))) 

write.csv(nje_cast, file = "output/nje_cast.csv",row.names=FALSE)


