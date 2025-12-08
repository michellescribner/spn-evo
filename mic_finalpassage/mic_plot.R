### Visualization of change in MIC for S. pneumoniae populations propagated in antibiotics

library(tidyverse)

setwd("~/Documents/pitt/streppneumo/mic_finalpassage/")

pal <- c("#287D8EFF", "#10A53DFF", "#541352FF")

# Read in raw input (duplicate readings were averaged in excel)
mic <- read.csv("rawdata/mic_input-2022-10-31-update.csv",header=TRUE)

druglist <- c( "Imipenem","Ciprofloxacin", "Azithromycin","Cefepime", "Levofloxacin", "Linezolid", "Penicillin", "Rifampicin","Vancomycin","Meropenem")

output <- data.frame()

for (i in 1:length(druglist)) {
  n <- subset(mic, mic$drug_tested == druglist[i])
  m <- subset(n, n$lineage =="TIGR4")
  n$anc <- median(m$ave_reading)
  n$mic_change <- n$ave_reading/n$anc
  n <- subset(n, n$drug == substr(druglist[i], 1,2))
  output <- rbind(output,n)
}

output <- output %>%
  group_by(lineage, drug, immune, rep, drug_tested) %>%
  summarize(median = median(mic_change, na.rm = T),
            ymin = min(mic_change, na.rm = T), 
            y25 = quantile(mic_change, 0.25),
            ymax  = max(mic_change, na.rm = T),
            y75 = quantile(mic_change, 0.25))
output$immune_rep <- paste(output$immune, output$rep, sep="")

output$immune <- recode(output$immune, "M0" = "Macrophage Depleted", "N0" = "Neutrophil Depleted", "Nd" = "Not Depleted")
output$rep <- recode(output$rep, "L1" = "Lineage 1", "L2" = "Lineage 2", "L3" = "Lineage 3")

p <- output %>%
  mutate(across(drug_tested, factor, levels=c(druglist))) %>%
  ggplot(aes(x=immune_rep,y=median,color = immune,shape=rep)) + 
  geom_point(size=2) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), stat = "identity") +
  scale_y_continuous(trans = 'log2') + 
  scale_shape(name= "Lineage") +
  scale_color_manual(values=pal, name="Immune Treatment") +
  theme_bw() + 
  theme(panel.grid.major.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~drug_tested, ncol = 3) + 
  ylab("Fold Change MIC") + 
  xlab("") 
ggsave(p, file = "mic_change_final_passage.pdf", device = "pdf", width = 6, height = 6)

