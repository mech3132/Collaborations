library(tidyverse)
# library(lme4)
# library(car)
# library(randomForest)
library(survival)
library(ranger)
library(ggfortify)

dat <- read.delim("survival_data.tsv") %>%
  mutate(Treatment = ifelse(Treatment=="CMR12", "CMR12a", Treatment)) %>%
  mutate(Treatment= factor(Treatment, levels=c("Control", unique(dat$Treatment)[-which(unique(dat$Treatment)=="Control")])))
  
genes <- read.delim("PA_gene_isolate.tsv") %>%
  mutate(Gene.Annotation = gsub("/","",Gene.Annotation, fixed = TRUE)) %>%
  mutate(Gene.Annotation = gsub(" ","",Gene.Annotation, fixed = TRUE)) %>%
  mutate(Gene.Identity = gsub("/","",Gene.Identity, fixed = TRUE)) %>%
  mutate(Gene.Identity = gsub(" |,|-|[(]|[)]|[&]","",Gene.Identity)) %>%
  mutate(Gene.Identity = gsub("24diacetylphloroglucinol","",Gene.Identity))

# All isolates
allIso <- unique(dat$Treatment)
allGeneComp <- unique(genes$Gene.Annotation)
allGenes <- unique(genes$Gene.Identity)
# Make unique fly names
dat <- dat %>%
  unite(ExpDate,Sex,Individual.Fly,Treatment, col="UniqueFly", sep = "", remove=FALSE) 

## Are all rows unique?
nrow(dat)
nrow(distinct(dat))
# No?


duplicated_flies <- dat %>% mutate(dup = duplicated(UniqueFly)) %>% 
  filter(dup) %>% pull(UniqueFly)
dat %>% filter(UniqueFly %in% duplicated_flies) %>%
  arrange(ExpDate, Treatment, Individual.Fly, Days)

# Plot flies
dat %>%
  ggplot(aes(x=Days, y=UniqueFly)) +
  geom_point(aes(col=factor(Turnover.event)))



### Kaplan Meier ###


dat %>% select(Days, Turnover.event)

km <- with(dat, Surv(Days, Turnover.event))
km_fit <- survfit(Surv(Days, Turnover.event) ~ 1, data=dat)
summary(km_fit)

autoplot(km_fit)

# By sex
km_sex <- survfit(Surv(Days, Turnover.event) ~ Sex, data=dat)
summary(km_sex)
autoplot(km_sex)

# By expdate
km_exp <- survfit(Surv(Days, Turnover.event) ~ ExpDate, data=dat)
summary(km_exp)
autoplot(km_exp)

# My treatment
dat %>% select(Treatment, Days) %>% table()

km_treatment <- survfit(Surv(Days, Turnover.event) ~ Treatment , data=dat)
summary(km_treatment)
ggsurvfit(km_treatment) 

autoplot(km_treatment) 

### Cox
cox <- coxph(Surv(Days, Turnover.event) ~ Treatment + Sex, data=dat)
summary(cox)
cox$coefficients

cox_fit <- survfit(cox)
summary(cox_fit)
autoplot(cox_fit)
 
