library(tidyverse)
library(lme4)
library(car)
library(randomForest)

dat <- read.delim("survival_data.tsv")
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

### Edit to calculate turnover event as numerical "days survived"
dat %>% select(Days, ExpDate) %>% table() # Max days is 11

# First, extract flies that need to be removed from exp because didn't reach day 11
toremove <- dat %>% 
  group_by(UniqueFly, Turnover.event) %>%
  summarise(maxDay=max(Days)) %>%
  filter(maxDay<11 & Turnover.event==0) %>%
  select(UniqueFly) %>%
  mutate(KEEP=FALSE) %>% ungroup()
# Now, join and filter out these flies
dat_filt <- dat %>%
  left_join(toremove) %>%
  filter(is.na(KEEP)) %>% select(-KEEP)

# Plot flies again; should have no zeros before 11 days
dat_filt %>%
  ggplot(aes(x=Days, y=UniqueFly)) +
  geom_point(aes(col=factor(Turnover.event)))

# Great- now make code that summarizes survived/not survived + days survived
# Also CMR12 isn't in dat sheet; change to CMR12a
dat_adj <- dat_filt %>%
  mutate(survived = (Days==11)&(Turnover.event==0)) %>%
  mutate(Treatment = ifelse(Treatment=="CMR12", "CMR12a", Treatment)) %>%
  mutate(Treatment= factor(Treatment, levels=c("Control", unique(dat_filt$Treatment)[-which(unique(dat_filt$Treatment)=="Control")]))) %>%
  mutate(earlyDeath = ifelse(survived,0,12-Days)) %>%
  filter(Days==11 | (Turnover.event>0)) 
dupflies <- dat_adj %>%
  mutate(dup = duplicated(UniqueFly)) %>%
  filter(dup) %>% pull(UniqueFly)
dat_adj %>% filter(UniqueFly %in% dupflies) %>% View()

nrow(dat_adj)
length(unique(dat_adj$UniqueFly))
# Look at overall rate of survival
dat_adj %>%
  ggplot() +
  geom_bar(aes(x=Treatment, fill=survived))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
# Look at early deaths
dat_adj %>%
  filter(earlyDeath>0) %>%
  ggplot() +
  geom_boxplot(aes(x=Treatment, y=earlyDeath))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
# Is survival rate correlated with how early flies die?
dat_adj %>%
  mutate(earlyDeath = ifelse(earlyDeath==0, NA, earlyDeath)) %>%
  group_by(Treatment) %>%
  summarise(survivalRate = sum(survived)/n(), averageDaysEarlyDeath = mean(earlyDeath, na.rm=TRUE)) %>%
  ggplot() + 
  geom_point(aes(x=survivalRate, y=averageDaysEarlyDeath, col=Treatment)) +
  xlab("Survival rate (to 11 days)") +
  ylab("Average number of days\nflies died early")

## Make two versions: one survival rate and another days died early
dat_survival <- dat_adj
dat_daysearly <- dat_adj %>% filter(earlyDeath>0)

## Some flies died early in control?
dat_adj %>%
  filter(!survived, Treatment=="Control")

### Really rudimentary model with only treatment
glmer_basic <- glmer(survived ~ Treatment + (1|ExpDate) + (1|Sex), data=dat_survival, family=binomial)
summary(glmer_basic)

glmer_basic_cont <- glm(earlyDeath ~ Treatment, data=dat_daysearly, family=poisson)
summary(glmer_basic_cont)


########## Model with Gene groups as predictors

# Table of genes
colnames(genes)
gene_by_iso <- genes %>% select(DSM21245:last_col())  %>% as.matrix()
rownames(gene_by_iso) <- genes$Gene.Annotation
gene_by_iso_t <- t(gene_by_iso)
# Fill in NAs as zeros
gene_by_iso_t[is.na(gene_by_iso_t)] <- 0
# Look at covariance matrix
covar_genes<- cov(gene_by_iso_t)
# Look at correlation between genes
heatmap(covar_genes, Rowv = NA, Colv = NA)
#I actually don't think this is too  bad?

# Let's group by gene cluster
geneID_by_iso <- genes %>% select(Gene.Identity,DSM21245:last_col()) %>%
  group_by(Gene.Identity) %>%
  summarise_all(max)
geneID_by_iso_t <- t(geneID_by_iso[,-1])
colnames(geneID_by_iso_t) <- geneID_by_iso$Gene.Identity
# Fill in NAs as zeros
geneID_by_iso_t[is.na(geneID_by_iso_t)] <- 0
# Look at covariance matrix
covar_geneID<- cov(geneID_by_iso_t)
# Look at correlation between genes
heatmap(covar_geneID)


### Join to metadata
genedf <- gene_by_iso_t %>% as.data.frame() %>%
  rownames_to_column(var="Treatment") 
geneIDdf <- geneID_by_iso_t %>% as.data.frame() %>%
  rownames_to_column(var="Treatment") 
# Add a control
tempcon <- matrix(0, ncol=ncol(genedf), nrow=1) %>% as.data.frame()
colnames(tempcon) <- colnames(genedf)
tempcon$Treatment[1] <- "Control"
genedf <- full_join(tempcon, genedf)

tempcon2 <- matrix(0, ncol=ncol(geneIDdf), nrow=1) %>% as.data.frame()
colnames(tempcon2) <- colnames(geneIDdf)
tempcon2$Treatment[1] <- "Control"
geneIDdf <- full_join(tempcon2, geneIDdf)

random_forest_genes <- dat_adj %>%
  select(survived, Treatment, ExpDate, Sex) %>%
  inner_join(genedf) %>%
  select(-Treatment)

random_forest_geneID <- dat_adj %>%
  select(survived, Treatment, ExpDate, Sex) %>%
  inner_join(geneIDdf) %>%
  select(-Treatment) %>%
  drop_na()


random_forest_all 

dat_adj %>%
  mutate(present=1) %>%
  select(UniqueFly, Treatment, present) %>%
  complete(UniqueFly,Treatment) %>%
  mutate(present = ifelse(is.na(present),0,present)) 

dat_adj%>%
  pivot_wider(names_from=Treatment, values_from=present)
  

test <- dat_adj %>%
  mutate(present=1) %>%
  select(UniqueFly, ExpDate, Sex, Treatment, present)%>%
  pivot_wider(names_from=Treatment, values_from=present, values_fill=0)


select(survived, Treatment, ExpDate, Sex) %>%
  inner_join(geneIDdf) %>%
  inner_join(genedf) %>%
  drop_na()

allGenesRandomForest <- randomForest(factor(survived)~., data=random_forest_genes)
allGenesRandomForest

allGeneIDsRandomForest <- randomForest(factor(survived)~., data=random_forest_geneID)
allGeneIDsRandomForest

everythingRandomForest <- randomForest(factor(survived)~.,data=random_forest_all)
everythingRandomForest

allGeneIDsRandomForest$importance %>% as.data.frame() %>%
  arrange(-MeanDecreaseGini)

everythingRandomForest$importance %>% as.data.frame() %>%
  arrange(-MeanDecreaseGini)

unique(genes$Gene.Identity)
unique(genes$Gene.Annotation)
nrow(dat)

###### summary quick
random_forest_all %>%
  select(-ExpDate, -Sex) %>%
  pivot_longer(c(-survived), values_to="PA", names_to="genetic_component")
random_forest_all %>%
  group_by(Entolysin) %>%
  summarise(prop_survive = sum(survived)/n())
  ggplot() +
  geom_jitter(aes(x = Entolysin, y=survived))
# # # Full table of genetic traits, long formal
# gene_ann <- genes %>% select(Gene.Identity, Gene.Annotation, Function)
# 
# gen_by_iso_df <- gene_by_iso %>% as.data.frame() %>%rownames_to_column(var="Gene.Annotation")
# gene_long <- gen_by_iso_df %>%
#   pivot_longer(-Gene.Annotation, names_to="Treatment", values_to="PA") %>%
#   full_join(gene_ann) %>%
#   select(Treatment, Gene.Annotation, Gene.Identity, PA)


