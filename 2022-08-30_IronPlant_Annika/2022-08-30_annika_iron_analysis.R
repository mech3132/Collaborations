#!bin/bash
library(tidyverse)
### Annika Iron

dat <- read.csv("annika_iron/iron_dat_reformatted.csv")
dat


dat_summarised <- dat %>% group_by(bio_rep, sample, bact, blank) %>%
  summarise( Conc_uM_mean = mean(Conc_uM)) %>% ungroup() %>%
  pivot_wider( names_from="blank", values_from="Conc_uM_mean") %>%
  mutate(diffConc_uM = tr-Bl) %>%
  mutate(Treatment = factor(ifelse(bact=="bact", "Bacteria\nadded", "No bacteria\nadded"),levels=c("No bacteria\nadded","Bacteria\nadded") ))
gg_adj <- dat_summarised %>%
  ggplot() + 
  geom_boxplot(aes(x=Treatment, y=diffConc_uM)) +
  geom_jitter(aes(x=Treatment, y=diffConc_uM), width=0.1) +
  facet_wrap(.~sample, nrow=1) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ylab("Difference in concentration\nbetween treat-blank (uM)")
gg_adj
ggsave(filename = "annika_iron/plot_conc_blank_subtracted.png",gg_adj, height=3, width=8)


gg_blanks <- dat_summarised %>%
  ggplot() + 
  geom_boxplot(aes(x=Treatment, y=Bl)) +
  geom_jitter(aes(x=Treatment, y=Bl), width=0.1) +
  facet_wrap(.~sample, nrow=1) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  ylab("Raw value of blanks (Conc_uM)")

gg_blanks
ggsave(filename = "annika_iron/plot_BLANKS.png",gg_blanks, height=3, width=8)



gg_rawval <- dat_summarised %>%
  ggplot() + 
  geom_boxplot(aes(x=Treatment, y=tr)) +
  geom_jitter(aes(x=Treatment, y=tr), width=0.1) +
  facet_wrap(.~sample, nrow=1) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  ylab("Raw value (Conc_uM)\n(Blank not subtracted)")
gg_rawval
ggsave(filename = "annika_iron/plot_raw_conc_blank_not_subtracted.png",gg_rawval, height=3, width=8)
