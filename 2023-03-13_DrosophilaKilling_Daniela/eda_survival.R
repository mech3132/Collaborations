# library(lme4)
# library(car)
# library(randomForest)
library(survival)
library(ranger)
library(ggfortify)
library(ggsurvfit)
library(survminer)
library(MASS)
library(tidyverse)



####### Load data #######
dat <- read.delim("survival_data.tsv") 
# CMR12a and CMR12 confirme to be same thing
# Factor so control is baseline
dat <- dat %>%
  mutate(Treatment = ifelse(Treatment=="CMR12", "CMR12a", Treatment)) %>%
  mutate(Treatment= factor(Treatment, levels=c("Control", unique(dat$Treatment)[-which(unique(dat$Treatment)=="Control")])))

# Get rid of extra symbols in Gene.Identity  
genes <- read.delim("PA_gene_isolate.tsv") %>%
  mutate(Gene.Annotation = gsub("/","",Gene.Annotation, fixed = TRUE)) %>%
  mutate(Gene.Annotation = gsub(" ","",Gene.Annotation, fixed = TRUE)) %>%
  mutate(Gene.Identity = gsub("/","",Gene.Identity, fixed = TRUE)) %>%
  mutate(Gene.Identity = gsub(" |,|-|[(]|[)]|[&]","",Gene.Identity)) %>%
  mutate(Gene.Identity = gsub("24diacetylphloroglucinol","",Gene.Identity))

####### Housekeeping #######

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
# Note; need to follow up with this with Daniela

####### EDA plot #######

# Plot flies
dat %>%
  ggplot(aes(x=Days, y=UniqueFly)) +
  geom_point(aes(col=factor(Turnover.event)))


########### Make Treatment x Gene table (and x Gene ID table) ############
# Table of genes
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

####### Full join to dataset #########
dat_ALLINFO <- left_join(dat, geneIDdf) %>%
  left_join(genedf)


#### Kaplan Meier ####
# Can only do one predictor at a time
# dat %>% select(Days, Turnover.event)
km <- with(dat, Surv(Days, Turnover.event))

# This is a basic fit with no predictors
km_fit <- survfit(km ~ 1, data=dat)
summary(km_fit)
autoplot(km_fit)

# See if there is a difference by sex
km_sex <- survfit(km ~ Sex, data=dat)
summary(km_sex)
autoplot(km_sex)
ggsurvplot(km_sex, data = dat, pval = TRUE)

# By expdate
km_exp <- survfit(km ~ ExpDate, data=dat)
summary(km_exp)
autoplot(km_exp)
ggsurvplot(km_exp, data = dat, pval = TRUE)


# My treatment

km_treatment <- survfit(km ~ Treatment , data=dat)
summary(km_treatment)
# ggsurvfit(km_treatment) # Could use this one too
autoplot(km_treatment, conf.int.fill = NA) 
ggsurvplot(km_treatment, data = dat, pval = TRUE)

# Use at log-ratio test for significance testing
survdiff(km ~ Treatment, data = dat)


# Could look at survival risk for each individual isolate(Restricted mean survival estimate)
km_treatment_results <- summary(km_treatment)
km_treatment_results$table %>%
  as.data.frame() %>% rownames_to_column(var="Treatment") %>%
  mutate(Treatment=gsub("Treatment=","",Treatment, fixed=TRUE))%>%
  arrange(-rmean) %>%
  mutate(Treatment=factor(Treatment, levels=unique(Treatment))) %>%
  ggplot(aes(x=Treatment)) + 
  geom_point(aes( y=rmean)) +
  geom_errorbar(aes(ymin=rmean-`se(rmean)`, ymax=rmean+`se(rmean)`)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ylab("Restricted Mean\nSurvival Estimate")
  


####### Cox Proportional Hazard models #########
cox_tr <- coxph(Surv(Days, Turnover.event) ~ Treatment, data=dat)
cox_tr.sex <- coxph(Surv(Days, Turnover.event) ~ Treatment + Sex, data=dat)
cox_tr.sex.date <- coxph(Surv(Days, Turnover.event) ~ Treatment + Sex+ExpDate, data=dat)
# summary(cox_tr.sex)
# Check AIC values; seems model without date is better
AIC(cox_tr.sex.date)
AIC(cox_tr.sex) # Best one
AIC(cox_tr)


gg_coxtest <- ggforest(cox_tr.sex, data = dat)
ggsave(filename = "cox_hr_plot.png", 
       gg_coxtest, width=6, height=5)

#### Run cox models on full gene sets #####
allgenes_formula <- formula(paste0("Surv(Days, Turnover.event) ~",paste(allGenes, collapse=" + ")))
cox_allgenes <- coxph(allgenes_formula, data=dat_ALLINFO)
AIC(cox_allgenes) # not as good as strain model
gg_coxallgenes <- ggforest(cox_allgenes, data = dat_ALLINFO)
ggsave(filename = "cox_hr_plot_allgenes.png", 
       gg_coxallgenes, width=6, height=10)
# Extract just genes of interest
if (file.exists("bestAIC_allgenes_model.RData")) {
  bestAIC_allgenes_model <- stepAIC(cox_allgenes)
  save(bestAIC_allgenes_model, file="bestAIC_allgenes_model.RData")
} else {
  load("bestAIC_allgenes_model.RData")
}
AIC(bestAIC_allgenes_model)
gg_coxbestgenes <- ggforest(bestAIC_allgenes_model, data = dat_ALLINFO)
ggsave(filename = "cox_hr_plot_bestgenes.png", 
       gg_coxbestgenes, width=8, height=4)

## All gene components
allgenecomp_formula <- formula(paste0("Surv(Days, Turnover.event) ~",paste(allGeneComp, collapse=" + ")))
cox_ALL_comp <- coxph(allgenecomp_formula, data=dat_ALLINFO)
AIC(cox_ALL_comp)
gg_coxallgenecomps <- ggforest(cox_ALL_comp, data = dat_ALLINFO)
ggsave(filename = "cox_hr_plot_allgenecomps.png", 
       gg_coxallgenecomps, width=6, height=20)
# Extract just gene components of interest
if (!file.exists("bestAIC_allgenecomp_model.RData")) {
  bestAIC_allgenecomp_model <- stepAIC(cox_ALL_comp) # Takes way to long to run
  save(bestAIC_allgenecomp_model, file="bestAIC_allgenecomp_model.RData")
} else {
  load("bestAIC_allgenecomp_model.RData")
}
gg_coxbestgenecomp <- ggforest(bestAIC_allgenecomp_model, data = dat_ALLINFO)
ggsave(filename = "cox_hr_plot_bestgenecomp.png",
       gg_coxbestgenecomp, width=8, height=4)
# 
# ###### Try looping through singles and comparing single vs group ############
## Genes
allGenesResults <- data.frame()
for ( i in allGenes ) {
  temp_frml <- formula(paste0("Surv(Days, Turnover.event) ~",i))
  temp_cox <- summary(coxph(temp_frml, data=dat_ALLINFO))
  allGenesResults <-  rbind(allGenesResults, as.data.frame(temp_cox$coefficients))
}

singleGeneResults <- allGenesResults %>% as.data.frame() %>%
  rownames_to_column(var = "Gene.Identity")
AICGeneResults <- summary(cox_allgenes)$coefficients %>% as.data.frame() %>%
  rownames_to_column(var = "Gene.Identity")
colnames(singleGeneResults) <- c("Gene.Identity","coef_single","expcoef_single","se_single","x_single","p_single")
colnames(AICGeneResults) <- c("Gene.Identity","coef_aic","expcoef_aic","se_aic","x_aic","p_aic")

full_join(singleGeneResults,AICGeneResults) %>%
  mutate(coef_aic = ifelse(is.na(coef_aic),0,coef_aic)) %>%
  ggplot(aes(x=coef_single, y=coef_aic)) + 
  geom_point(aes(col=p_single<0.05, fill=p_aic<0.05), pch=21)
# Fascinatingly, there are some genes that are only significant when others are present!

## Gene components
allGeneCompResults <- data.frame()
for ( i in allGeneComp ) {
  temp_frml <- formula(paste0("Surv(Days, Turnover.event) ~",i))
  temp_cox <- summary(coxph(temp_frml, data=dat_ALLINFO))
  allGeneCompResults <-  rbind(allGeneCompResults, as.data.frame(temp_cox$coefficients))
}

singleGeneCompResults <- allGeneCompResults %>% as.data.frame() %>%
  rownames_to_column(var = "Genes.Annotations")
AICGeneCompResults <- summary(cox_ALL_comp)$coefficients %>% as.data.frame() %>%
  rownames_to_column(var = "Gene.Annotations")
colnames(singleGeneCompResults) <- c("Gene.Annotations","coef_single","expcoef_single","se_single","x_single","p_single")
colnames(AICGeneCompResults) <- c("Gene.Annotations","coef_aic","expcoef_aic","se_aic","x_aic","p_aic")

full_join(singleGeneCompResults,AICGeneCompResults) %>%
  mutate(coef_aic = ifelse(is.na(coef_aic),0,coef_aic)) %>%
  ggplot(aes(x=coef_single, y=coef_aic)) + 
  geom_point(aes(col=p_single<0.05, fill=p_aic<0.05), pch=21)

######## Try going from broad to narrow; look at strains first, then correlate with genes

# Previous model; use to order things
cox_tr.sex
order_of_iso_effects <- cox_tr.sex$coefficients %>% sort() %>%
  names() %>% gsub("Treatment","",.)
order_of_iso_effects <- c("Control",order_of_iso_effects[-which(order_of_iso_effects=="SexM")])

cox_allgenes
order_of_gene_effects <- cox_allgenes$coefficients %>% sort(na.last = TRUE) %>%
  names()



# Use hclust on gene identity appearances so that they are clustered by how common
# they are in strains
dm_genes <- dist(t(geneID_by_iso_t),method = "euclidean")
hclust_genes <- hclust(dm_genes)

dm_isolates <-dist(geneID_by_iso_t,method = "euclidean")
hclust_isolates <- hclust(dm_isolates)

gg_heatmap_genes <- genes %>%
  pivot_longer(-c(Gene.Identity, Gene.Annotation, Function), names_to="Treatment", values_to="PA") %>%
  # mutate(Treatment=factor(Treatment, levels=order_of_iso_effects)) %>%
  mutate(Treatment=factor(Treatment, levels=hclust_isolates$labels)) %>%
  # mutate(Gene.Identity=factor(Gene.Identity, levels=hclust_genes$labels)) %>%
  mutate(Gene.Identity=factor(Gene.Identity, levels=order_of_gene_effects)) %>%
  ggplot() +
  geom_tile(aes(x=Gene.Annotation, y=Treatment, fill=factor(PA))) +
  facet_grid(.~Gene.Identity, scales="free", drop=TRUE, space = "free_x") +
  scale_fill_manual(values=c("grey","red"))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

ggsave(filename = "heatmap_genes_PA.png", 
       gg_heatmap_genes, height=4, width=25)

### Want to do random forest, but with cox models?
