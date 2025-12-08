### Lung Titer Association with Mutations 

library(tidyverse)
library(readr)
library(matrixStats)

setwd("~/Documents/pitt/streppneumo/lungtiters/rawdata/")

pal <- c("#287D8EFF", "#10A53DFF", "#541352FF")

# Plot antibiotic concentration regime

conc <- read.csv("abconc.csv")
l <- pivot_longer(conc, cols = !gen, names_to = "lineage") %>% 
  separate(lineage, into = c("drug", "immune", "rep"), sep = c(2,4)) %>%
  unique()

l %>% 
  filter(drug %in% c("Im", "Ci", "Az")) %>%
  mutate(across(drug, factor, levels=c("Im", "Ci", "Az"))) %>%
  filter(gen < 15) %>%
  ggplot(mapping = aes(x=gen, y=value, shape=drug)) + geom_line() + geom_point(size =1) +
  scale_y_continuous(trans = "log2") + ylab("Drug Concentration (mg/kg)") + xlab("Passage") + theme_bw() +
  facet_wrap(vars(drug), nrow=4, ncol = 3) + 
  theme(strip.background = element_blank(), strip.text.x = element_blank())
ggsave("block1_drug.pdf", plot = last_plot(), device = "pdf", path = "~/Documents/pitt/streppneumo/lungtiters/", height = 1, width = 9)

l %>% 
  filter(drug %in% c("Ce", "Le", "Li")) %>%
  mutate(across(drug, factor, levels=c("Ce", "Le", "Li"))) %>%
  filter(gen < 30) %>%
  ggplot(mapping = aes(x=gen, y=value, shape=drug)) + geom_line() + geom_point(size =1) +
  scale_y_continuous(trans = "log2") + ylab("Drug Concentration (mg/kg)") + xlab("Passage") + theme_bw() +
  facet_wrap(vars(drug), nrow=4, ncol = 3) + 
  theme(strip.background = element_blank(), strip.text.x = element_blank())
ggsave("block2_drug.pdf", plot = last_plot(), device = "pdf", path = "~/Documents/pitt/streppneumo/lungtiters/", height = 1, width = 9)

l %>% 
  filter(drug %in% c("Pe","Ri","Va","Me")) %>%
  mutate(across(drug, factor, levels=c("Pe","Ri","Va","Me"))) %>%
  filter(gen < 15) %>%
  ggplot(mapping = aes(x=gen, y=value, shape=drug)) + geom_line() + geom_point(size =1) +
  scale_y_continuous(trans = "log2") + ylab("Drug Concentration (mg/kg)") + xlab("Passage") + theme_bw() +
  facet_wrap(vars(drug), nrow=1, ncol = 4) + 
  theme(strip.background = element_blank(), strip.text.x = element_blank())
ggsave("block3_drug.pdf", plot = last_plot(), device = "pdf", path = "~/Documents/pitt/streppneumo/lungtiters/", height = 1, width = 12)


# Plot CFU dynamics

df <- read.csv("combined.csv")

df2 <- pivot_longer(df, cols = !gen, names_to = "lineage") %>%
  mutate(log10 = log10(value+1)) %>% 
  separate(lineage, into = c("drug", "immune", "rep"), sep = c(2,4))

df3 <- df2

columns <- c("L1", "L2", "L3")

df2 <- df2 %>%
  pivot_wider(id_cols = c(gen, drug, immune), names_from = rep, values_from = log10) %>%
  mutate(mean = rowMeans(select(., all_of(columns))),
         sd = rowSds(as.matrix(select(., all_of(columns)))))

df2$drug <- recode(df2$drug, "Az" = "Azithromycin", "Ci" = "Ciprofloxacin", "Im"="Imipenem",
                 "Ce"="Cefepime","Le"="Levofloxacin", "Li"="Linezolid",
                 "Pe"="Penicillin","Me"="Meropenem", "Ri"="Rifampicin","Va"="Vancomycin")

df2 %>% 
  filter(drug %in% c("Imipenem", "Ciprofloxacin", "Azithromycin")) %>%
  mutate(across(drug, factor, levels=c("Imipenem", "Ciprofloxacin", "Azithromycin"))) %>%
  filter(gen < 15) %>%
  ggplot(mapping = aes(x=gen, y=mean, color=immune)) + 
  geom_line(mapping = aes(group = interaction(drug, immune))) +
  geom_ribbon(aes(ymin = mean-sd, ymax= mean + sd, fill = immune), alpha = 0.3, color = NA) +
  geom_point(size =1) +
  scale_y_continuous(limits = c(3,10)) + 
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) +
  ylab("CFU per Lung (log10)") +
  xlab("Passage") + 
  theme_bw() +
  facet_wrap(vars(drug))
ggsave("block1.pdf", plot = last_plot(), device = "pdf", path = "~/Documents/pitt/streppneumo/lungtiters/", height = 3, width = 9)

df2 %>% 
  filter(drug %in% c("Cefepime", "Levofloxacin", "Linezolid")) %>%
  filter(gen < 30) %>%
  ggplot(mapping = aes(x=gen, y=mean, color=immune)) + 
  geom_line(mapping = aes(group = interaction(drug, immune))) +
  geom_ribbon(aes(ymin = mean-sd, ymax= mean + sd, fill = immune), alpha = 0.3, color = NA) +
  geom_point(size =1) +
  scale_y_continuous(limits = c(3,10)) + 
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) +
  ylab("CFU per Lung (log10)") +
  xlab("Passage") + 
  theme_bw() +
  facet_wrap(vars(drug))
ggsave("block2.pdf", plot = last_plot(), device = "pdf", path = "~/Documents/pitt/streppneumo/lungtiters/", height = 3, width = 9)

df2 %>% 
  filter(drug %in% c("Vancomycin", "Penicillin", "Meropenem", "Rifampicin")) %>%
  mutate(across(drug, factor, levels=c("Penicillin","Rifampicin","Vancomycin","Meropenem"))) %>%
  filter(gen < 15) %>%
  ggplot(mapping = aes(x=gen, y=mean, color=immune)) + 
  geom_line(mapping = aes(group = interaction(drug, immune))) +
  geom_ribbon(aes(ymin = mean-sd, ymax= mean + sd, fill = immune), alpha = 0.3, color = NA) +
  geom_point(size =1) +
  scale_y_continuous(limits = c(3,10)) + 
  scale_fill_manual(values=pal) +
  scale_color_manual(values=pal) +
  ylab("CFU per Lung (log10)") +
  xlab("Passage") + 
  theme_bw() +
  facet_wrap(vars(drug), nrow = 1, ncol = 4)
ggsave("block3.pdf", plot = last_plot(), device = "pdf", path = "~/Documents/pitt/streppneumo/lungtiters/", width = 12, height = 3)


