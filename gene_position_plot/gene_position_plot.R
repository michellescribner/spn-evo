# Plotting the locations of mutations within genes that acquired mutations in parallel

library(tidyverse)
library(viridisLite)
library(ggrepel)
library(cowplot)

pal <- c("#287D8EFF", "#10A53DFF", "#541352FF")

setwd("~/Documents/pitt/streppneumo/gene_position_plot")

snps <- read.csv("~/Documents/pitt/streppneumo/post_breseq_filtering/output/snps_after_filtering.csv")

# Create a table of only ribonuclease Y mutations
df <- snps
# df$drug <- recode(df$drug, "Az" = "Azithromycin", "Ci" = "Ciprofloxacin", "Im"="Imipenem",
#                     "Ce"="Cefepime","Le"="Levofloxacin", "Li"="Linezolid",
#                     "Pe"="Penicillin","Me"="Meropenem", "Ri"="Rifampicin","Va"="Vancomycin")
df <- df %>% 
  filter(gene == "SP_1739") %>% 
  group_by(position,annotation,description,drug) %>%
  summarise(n = n()) %>%
  pivot_wider(id_cols = c(position,annotation,description),names_from = drug,values_from=n,values_fill = 0)
write.csv(df, file="ribonucleaseY_table.csv",row.names=FALSE)


# # Subset and clean SP_1583 and SP_1739 data
# missense and nonsense mutations
snps_a <- snps %>%
  filter(gene == "SP_1739" | gene == "SP_1583") %>%
  separate(annotation, into = c("one", "two", "three"), sep = " ") %>%
  filter(one != "coding") %>%
  separate(one, into = c("start", "aa_position", "end"), sep = c(1, -1),remove=FALSE) %>%
  select(drug,immune,pop,mutation,gene,one,start,aa_position,end,description,freq,gen) %>%
  transform(aa_position = as.numeric(aa_position)) %>%
  group_by(drug,immune,mutation,gene,one,start,aa_position,end,description) %>%
  summarize(number_of_samples = n())
snps_a$type <- "SNP"
snps_a$type[snps_a$end == "*"] <- "Stop"

# indels
snps_b <- snps %>%
  filter(gene == "SP_1739" | gene == "SP_1583") %>%
  separate(annotation, into = c("one", "two", "three"), sep = " ") %>%
  filter(one == "coding") %>%
  separate(two, into = c("start", "aa_position"), sep = 1) %>%
  separate(aa_position, into = c("aa_position", "end"), sep = "/") %>%
  separate(aa_position, into = c("aa_position", "end"), sep = "-", ) %>%
  select(drug,immune,pop,mutation,gene,start,aa_position,end,description,freq,gen) %>%
  transform(aa_position = as.numeric(aa_position)) %>%
  mutate(aa_position = aa_position/3) %>%
  mutate(one = NA) %>%
  group_by(drug,immune,mutation,gene,start,one,aa_position,end,description) %>%
  summarize(number_of_samples = n())
snps_b$type <- "Indel"

snps <- rbind(snps_a, snps_b)

# snps$drug <- recode(snps$drug, "Az" = "Azithromycin", "Ci" = "Ciprofloxacin", 
#                     "Im"="Imipenem","Ce"="Cefepime","Le"="Levofloxacin", "Li"="Linezolid", 
#                     "Pe"="Penicillin","Me"="Meropenem", "Ri"="Rifampicin","Va"="Vancomycin")
# snps$immune <- recode(snps$immune, "M0" = "Macrophage\nDepleted", "N0" = "Neutrophil \n Depleted", "Nd" = "Not Depleted")

# set gene lengths
iso_aa <- 576/3
ribo_aa <- 550

# plot SP_1583
iso_p <- snps %>% 
  filter(gene == "SP_1583") %>%
  ggplot(aes(y=drug, x=aa_position, color=immune,shape=type,size=number_of_samples,label=one)) + 
  geom_point(stroke=2) +
  scale_size_continuous(range = c(1,4), limits=c(0,40),breaks=c(0,20,40), guide = guide_legend(ncol = 3, byrow = TRUE, title.position = "top")) +
  scale_shape_manual(values=c(2,1,8),guide = guide_legend(ncol = 3, byrow = TRUE, title.position = "top")) +
  scale_color_manual(values=pal,guide = guide_legend(ncol = 3, byrow = TRUE, title.position = "top"))+
  scale_x_continuous(limits = c(0, iso_aa), breaks = c(0,100))+
  theme_minimal() + 
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.position = "right",
        legend.box = "vertical",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) +
  geom_segment(
    aes(x = 0, xend = iso_aa, y = 10, yend = 10),
    inherit.aes = FALSE,
    size = 5,
    lineend = "butt",
    linejoin = "mitre",
    colour = "light gray",
    arrow = arrow(length = unit(0.01, "cm"), type = "closed")) +
  coord_cartesian(clip = "off") +
  xlab("Amino Acid") + ylab(element_blank()) + 
  labs(color="Immune Treatment",shape="Mutation Type",size="Number of Samples") +
  geom_text(
    x = iso_aa / 2,                     
    y = 10,                       
    label = "SP_1583",    
    inherit.aes = FALSE,
    hjust = 0.5, vjust = 0.5,                 
    size = 3) 
iso_p
# ggsave(iso_p,filename= "SP_1583_cysteinehydrolase.png", dpi=300, device = "png", width = 6, height = 2)

# plot ribonuclease Y
ribo_p <- snps %>% 
  filter(gene == "SP_1739") %>%
  ggplot(aes(y=drug, x=aa_position, color=immune,shape=type,size=number_of_samples)) + 
  geom_point(stroke=2) +
  scale_size(range = c(1,4), limits=c(0,40), breaks = c(0,20,40), guide = guide_legend(ncol = 3, byrow = TRUE, title.position = "top")) +
  scale_shape_manual(values=c(1),guide = guide_legend(ncol = 3, byrow = TRUE, title.position = "top")) +
  scale_x_continuous(limits = c(0, ribo_aa))+
  scale_color_manual(values=pal,guide = guide_legend(ncol = 3, byrow = TRUE, title.position = "top"))+ 
  labs(color="Immune Treatment",shape="Mutation Type",size="Number of Samples") +
  xlab(element_blank()) + ylab(element_blank()) +
  theme_minimal() + 
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        legend.position = "none",
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) +
  geom_segment(
    aes(x = 0, xend = ribo_aa, y = 7, yend = 7),
    inherit.aes = FALSE,
    size = 5,                         
    lineend = "butt",   
    linejoin = "mitre",
    colour = "light gray",
    arrow = arrow(length = unit(0.01, "cm"), type = "closed")) +
  coord_cartesian(clip = "off") +
  geom_text(
    x = ribo_aa / 2,                     
    y = 7,                       
    label = "SP_1739",    
    inherit.aes = FALSE,
    hjust = 0.5, vjust = 0.5,                 
    size = 3) 

ribo_p
# ggsave(ribo_p,filename= "SP_1739_ribonucleaseY.png", dpi=300, device = "png", width = 6, height = 1.7)

#------------- Cowplot  --------------------
bottom_right <- plot_grid(ribo_p, iso_p, ncol = 1, rel_heights = c(1, 1.6), labels = c("D", ""))
bottom_right

