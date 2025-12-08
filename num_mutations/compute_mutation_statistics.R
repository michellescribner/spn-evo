# # # Compute number of mutations and mutation type table after filtering 
# for S. pneumoniae populations propagated in the presence of ABX

library(tidyverse)

setwd("~/Documents/pitt/streppneumo/num_mutations")

snps <- read.csv("~/Documents/pitt/streppneumo/post_breseq_filtering/output/snps_after_filtering.csv", header = TRUE)

# count the number of unique mutations that occurred per lineage
# mutations that occurred in multiple timepoints are counted only ONCE
nmut <- snps %>%
  select(lineage,drug,immune,pop,position,mutation,gene,annotation,description) %>%
  distinct() %>%
  group_by(lineage,drug,immune,pop) %>%
  summarize(count= n()) 
write.csv(nmut, file="number_of_mutuations.csv",row.names=FALSE)

# Tidy the metadata
nmut$drug <- recode(nmut$drug, "Az" = "Azithromycin", "Ci" = "Ciprofloxacin", "Im"="Imipenem",
                    "Ce"="Cefepime","Le"="Levofloxacin", "Li"="Linezolid",
                    "Pe"="Penicillin","Me"="Meropenem", "Ri"="Rifampicin","Va"="Vancomycin")
nmut$immune <- recode(nmut$immune, "M0" = "Macrophage Depleted", "N0" = "Neutrophil Depleted", "Nd" = "Not Depleted")
nmut$pop <- recode(nmut$pop, "L1" = "Lineage 1", "L2" = "Lineage 2", "L3" = "Lineage 3")

# Plot mutation counts by lineage
pal <- c("#287D8EFF", "#10A53DFF", "#541352FF")
nmut%>%
  ungroup %>%
  mutate(across(drug, factor, levels=c("Imipenem", "Ciprofloxacin", "Azithromycin","Cefepime", "Levofloxacin", "Linezolid","Penicillin","Rifampicin","Vancomycin","Meropenem"))) %>%
  ggplot(aes(x=drug, y=count, color = immune, group = immune, shape = pop)) + 
  geom_point(size=3, position=position_dodge(width = 1)) +
  scale_color_manual(values=pal)+
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Antibiotic Treatment") + 
  ylab('Number of Mutations') +
  labs(color = "Immune Treatment") +
  labs(shape = "Lineage")
ggsave("number_of_mutations_plot.pdf", plot = last_plot(), device = "pdf", height =4, width =7)


# # Create a mutation type table

# create table listing distinct mutations
dmut <- snps %>%
  select(position,mutation,gene,annotation,description) %>%
  distinct() 
nrow(dmut) 

dmut$annotation <- gsub("\\(","",as.character(dmut$annotation))
dmut$annotation <- gsub("\\)","",as.character(dmut$annotation))
dmut$annotation <- gsub("nt$","",as.character(dmut$annotation))

# Determine mutation types
dmut$MutationType <- "missense"
dmut$MutationType[grep("intergenic", dmut$annotation)] <- "intergenic"
dmut$MutationType[grep("\\*", dmut$annotation)] <- "nonsense"
dmut$MutationType[grep("coding", dmut$annotation)] <- "indel"
dmut$MutationType[grep("pseudogene", dmut$annotation)] <- "pseudogene"

# separate out the missense mutations, determine if synonymous
dmut1 <- dmut %>%
  filter(MutationType == "missense") %>%
  separate(annotation, into = c("one", "two","three"), sep = " ",remove=FALSE) %>%
  separate(one, into = c("start", "aa_position", "end"), sep = c(1, -1),remove=FALSE) %>%
  transform(aa_position = as.numeric(aa_position)) 

for (i in 1:nrow(dmut1)) {
  if(dmut1$start[i] == dmut1$end[i]) {
    dmut1$MutationType[i] <- "synonymous SNP"
    print(paste(dmut1$annotation[i],"is a synonymous SNP"))
  } else {
    dmut1$MutationType[i] <- "nonsynonymous SNP"
    print(paste(dmut1$annotation[i],"is a nonsynonymous SNP"))
  }
}

# recombine missense mutations with other mutations
dmut1 <- dmut1 %>%
  select(position, mutation, gene, annotation, description, MutationType)
dmut <- dmut %>% filter(MutationType != "missense")
dmut <- rbind(dmut, dmut1)

# summarize in table
table <- dmut %>%
  group_by(MutationType) %>%
  summarize(count= n())
write.csv(table, "mutation_type_table.csv",row.names=FALSE)

