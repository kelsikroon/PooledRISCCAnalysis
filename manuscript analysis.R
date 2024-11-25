# Manuscript analysis:

# devtools::install_github("kelsikroon/PICmodel") # install latest version of PICmodel R-package 

# load packages ---------------------
library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)
library(tidyr)
library(purrr)
library(PICmodel)
library(ggrepel)
library(ggpubr)
library(scales)
source("data_preparation.R")
source("plot_functions.R")

library(rcartocolor) # rcartocolor::carto_pal(n = 7, name = "Bold")
theme_set(theme_classic())

# --------------- 1. data descriptives ---------------
combo.dat.cin2.all <- create.combo.dat('cin2plus', trim.age = F)
combo.dat.cin3.all <- create.combo.dat('cin3plus', trim.age = F)
combo.dat.ca.all <- create.combo.dat('cancer', trim.age = F)

study.names <- c('POBASCAM', 'VUSAscreen', 'IMPROVE', 'NTCCitaly', 'ULslovenia', 'SwedeScreen' )

combo.dat <- do.call(rbind, lapply(study.names, data.descript))

table(combo.dat$worst_histo, combo.dat$study)

table(combo.dat.cin2.all$Source)

# (a) age/ cytology information 
age_cyt.dat <- combo.dat.cin2.all[, c("Source", "age", "cyt")]
age_cyt.dat$Source <- factor(age_cyt.dat$Source, levels= c("pob", "vusa", "improve", "italy", "slov", "swed"))

age_cyt.dat %>% group_by(Source) %>% summarise(total = n(), age.med = median(age), age.range = paste0(round(range(age)[1]), " - ", round(range(age)[2])), 
                                               cyt.NILM = paste0(sum(cyt == 1), " (", round(sum(cyt==1)/total*100, 2), ")"),  
                                               cyt.ASCUS = paste0(sum(cyt > 1 & cyt <=3), " (", round(sum(cyt > 1 & cyt <=3)/total*100, 2), ")"), 
                                               cyt.HSIL = paste0(sum(cyt >= 4), " (", round(sum(cyt >= 4)/total*100, 2), ")"))

# (b) CIN2/CIN3/cancer cases
hist.info <- combo.dat.cin2.all %>% group_by(Source) %>% 
  summarise(N= n(), cin2 = sum(worst_histo == 7), cin3 = sum(worst_histo == 8), cancer = sum(worst_histo %in% c(9, 10, 12))) %>% 
  as.data.frame() %>% 
  mutate(perc.cin2 = paste0(cin2, " (", round(cin2/N*100, 2), ")"),
         perc.cin3 = paste0(cin3, " (", round(cin3/N*100, 2), ")"),
         perc.ca =  paste0(cancer, " (", round(cancer/N*100, 2), ")")) %>% 
  select(c(Source, perc.cin2, perc.cin3, perc.ca)) ; hist.info

# (c) genotype distribution 
genotypes <- paste0("hpv", c(16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59, 66, 68))
genotypes.study <- combo.dat.cin2.all[,c("Source", genotypes)] %>% group_by(Source) %>% 
  summarise(across(hpv16:hpv68, sum, na.rm=T)) %>% 
  mutate(round(.[2:15]/rowSums(.[2:15]), 4)) %>% 
  pivot_longer(cols = `hpv16`:`hpv68`, names_to = "genotype",values_to = "perc") #%>% mutate(across(genotype), factor)

# (Figure 1) genotype summary  ---------------------
genotypes.overall <- data.frame(perc=colSums(combo.dat.cin2.all[,genotypes], na.rm=T)/dim(combo.dat.cin2.all)[1], 
                                genotype=genotypes, Source="Overall") 
genotypes.overall$ymaxes <- genotypes.study %>% group_by(genotype) %>% summarise(max = max(perc)) %>% select(max) %>% unlist()

study.genotype.plot <- ggplot(genotypes.study, aes(x=genotype, y=perc, fill=Source) ) + 
  geom_bar(stat='identity', position = position_dodge(width=0.5), width=0.4) +
  xlab("Genotype") + ylab("Positive percentage") + 
  labs(title= "Distribution of HPV genotypes in all HPV-positive women across six European studies") + 
  theme(text = element_text(size=12), legend.position.inside=c(0.9, 0.76), legend.text = element_text(size=9)) + 
  guides(fill = guide_legend(position='inside')) + 
  scale_y_continuous(expand=c(0,0), breaks=seq(0, 0.35, 0.05), labels=scales::percent, limits=c(0, 0.34)) + 
  scale_x_discrete(labels =  paste0("HPV", c(16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59, 66, 68)))  + 
  scale_fill_manual("", values=rcartocolor::carto_pal(n = 7, name = "Bold"), 
                    labels= c("POBASCAM", "VUSA-Screen", "IMPROVE", "NTCC", "Slovenian HPV Prevalence", "Swedescreen"),
                    limits = c("pob", "vusa", "improve", "italy", "slov", "swed")) + 
  geom_errorbar(data=genotypes.overall, 
                aes(x=factor(genotype), y=perc, ymin=perc, ymax=perc), col='black') + 
  geom_text(data=genotypes.overall,
            aes(x=genotype, y= ymaxes + 0.01, label=paste0(specify_decimal(perc*100, 1), "%")))
# save 950*400 size (same plot but only in those who had CIN2+/CIN3+ or cancer?)


ggpubr::ggarrange(study.genotype.plot, single.type.plot,  ncol=1, align='v', labels=c("A", "B"))


genotypes.study[genotypes.study$genotype == 'hpv16',]
genotypes.study[genotypes.study$genotype == 'hpv31',]

combo.dat.cin2.all[,c("Source", genotypes)] %>% ungroup() %>% mutate(n.pos = rowSums(.[2:15], na.rm=T)) %>% 
  mutate(more.than.1 = ifelse(n.pos > 1, 1, 0)) %>%  group_by(Source, more.than.1) %>% 
  summarise(sum.source = sum(more.than.1, na.rm=T)) %>% print(n=29)

# --------------- 2. individual study models for each endpoint (without cytology)---------------
combo.dat.cin2 <- create.combo.dat('cin2plus', interval.type = 'hierarchical', trim.age = T)
combo.dat.cin3 <- create.combo.dat('cin3plus', interval.type = 'heirarchical', trim.age = T)
combo.dat.cancer <- create.combo.dat('cancer', interval.type = 'hierarchical', trim.age = T)

studies <- c("pob", "vusa", "improve", "italy", "slov", "swed")

study.fits <- function(endpoint.data){
  covariates <- c("age.std.40", "hpv16")
  
  temp.data <- endpoint.data[endpoint.data$age <= 40, ]
  temp.data <- temp.data %>% group_by(Source) %>% mutate(age.std.40 =  0.5*(age - mean(age))/sd(age)) # add age.std for each study individually 
  model.fits <- lapply(studies, function(x)  PICmodel.fit(covariates, c(), covariates, temp.data[temp.data$Source ==x,], include.h = F, silent = F, short.runs=30, prior.type = 't4'))
  summaries.list <- lapply(model.fits, function(x) x$summary) 
  summaries.df <- cbind(summaries.list[[1]]$param, matrix(unlist(lapply(summaries.list, function(x) paste0(specify_decimal(x$theta.hat, 2), " (", specify_decimal(x$lower, 2), ", ", specify_decimal(x$upper, 2), ")"))), ncol=6)) %>% as.data.frame()
  colnames(summaries.df) <- c("parameter", studies)
  return(summaries.df)
}
# Table 2
gtools::smartbind(data.frame(parameter= "CIN2+", pob="",vusa="",improve="",italy="",slov="",swed=""), study.fits(combo.dat.cin2), 
                  data.frame(parameter= "CIN3+", pob="",vusa="",improve="",italy="",slov="",swed=""), study.fits(combo.dat.cin3), fill=0) %>% 
  write.csv(., "~/Projects/RISCC Comparison/Results/table2_20241107.csv")

# --------------- 3. single- and double- type genotype distribution in CIN2+, CIN3+ and cancer cases---------------

table3 <- do.call(rbind, list(cin2= genotype.dist(combo.dat.cin2.all), cin3= genotype.dist(combo.dat.cin3.all), ca= genotype.dist(combo.dat.ca.all)))




# --------------- 4. model fit: no age/genotype restrictions, general intercept for both prevalence and progression---------------
grouped.types <- c("hpv16", "hpv18.45", "hpvGroup1Alpha9")
# --> model fits are commented out because they take approx 10-15 mins each to run, so quicker to load from the saved file
# (a) cin2+  
# cin2pl.fit.grouped <- PICmodel.fit(c(grouped.types), c(), c(grouped.types), combo.dat.cin2, include.h = F, silent = F, short.runs=10, include.priors=F) ; cin2pl.fit.grouped$label <- "cin2pl"
# save(cin2pl.fit.grouped, file = paste0('fit_grouped_', format(Sys.Date(), "%Y%m%d"), "_cin2pl", '.RData'))
load("fit_grouped_20241107_cin2pl.RData")

# (b) cin3+
# cin3pl.fit.grouped <- PICmodel.fit(c(grouped.types), c(), c(grouped.types), combo.dat.cin3, include.h = F, silent = F, short.runs=10, include.priors=F) ; cin3pl.fit.grouped$label <- "cin3pl"
# save(cin3pl.fit.grouped, file = paste0('fit_grouped_', format(Sys.Date(), "%Y%m%d"), "_cin3pl", '.RData'))
load("fit_grouped_20241107_cin3pl.RData")

# (c) cancer 
# cancer.fit.grouped <- PICmodel.fit(c(grouped.types), c(), c(grouped.types), combo.dat.cancer, include.h = F, silent = F,short.runs=15, include.priors=T) ; cancer.fit.grouped$label <- 'cancer'
# save(cancer.fit.grouped, file = paste0('fit_grouped_', format(Sys.Date(), "%Y%m%d"), "_cancer", '.RData'))
load("fit_grouped_20241107_cancer.RData")

# individual genotypes for CIn2+ and CIN3+
cov.genotypes <- c("hpv16", "hpv18", "hpv31", "hpv33", "hpv35", "hpv45", "hpv52", "hpv58")
# (a) cin2+  
# cin2pl.fit <- PICmodel.fit(c(cov.genotypes), c(), c(cov.genotypes), combo.dat.cin2, include.h = F, silent = F, short.runs=10, include.priors=F) ; cin2pl.fit$label <- "cin2pl"
# save(cin2pl.fit, file = paste0('fit_', format(Sys.Date(), "%Y%m%d"), "_cin2pl", '.RData'))
load("fit_20241107_cin2pl.RData")

# (b) cin3+
# cin3pl.fit <- PICmodel.fit(c(cov.genotypes), c(), c(cov.genotypes), combo.dat.cin3, include.h = F, silent = F, short.runs=10, include.priors=F) ; cin3pl.fit$label <- "cin3pl"
# save(cin3pl.fit, file = paste0('fit_', format(Sys.Date(), "%Y%m%d"), "_cin3pl", '.RData'))
load("fit_20241107_cin3pl.RData")

# Figure 2: genotype forest plot for all 3 endpoints with grouped and also individual
grouped.res <- pooled.fits.endpoint('grouped_20241107') %>% subset((label != 'CIN2+' & label!= 'CIN3+') | param !='hpv16')
genotype.values <- rbind(grouped.res, pooled.fits.endpoint('20241107'))
genotype.values$param[genotype.values$param == 'hpv18.45'] <- 'hpv18/45'
genotype.values$param <- factor(genotype.values$param, levels= c("hpv16", "hpv18/45", "hpvGroup1Alpha9", "hpv18","hpv45",  "hpv31", "hpv33", 'hpv35', "hpv52", "hpv58", "intercept"))
genotype.values <- genotype.values[genotype.values$param != 'intercept',]
genotype.forest.plot(genotype.values, legend.lab='Endpoint', dodge.height = 0.4)

# number in results section:
genotype.values %>% arrange(group) %>%  subset(param=='hpvGroup1Alpha9' & label=='Cancer') %>% 
  mutate_if(is.numeric, exp) %>% mutate_if(is.numeric, round, 1)


# Figure 3: prediction plot
cancer.results.grouped <- pred.and.plot(cancer.fit.grouped, "Cancer")
cin2pl.results <- pred.and.plot(cin2pl.fit, "CIN2+")
cin3pl.results <- pred.and.plot(cin3pl.fit, "CIN3+")

ggpubr::ggarrange(cin2pl.results[[1]], cin3pl.results[[1]], cancer.results.grouped[[1]], ncol=3, 
                  align='h', labels=c("A", "B", "C"))



# Table (Supplementary):
cin2pl.results[[2]]%>%  as.data.frame() %>% mutate_if(is.numeric, round, 1) %>% 
  mutate(baseline.risk = paste0(CR.0, " (", CR.lower95.0, " - ", CR.upper95.0, ")"), 
         fifteen.risk = paste0(CR.15, " (", CR.lower95.15, " - ", CR.upper95.15, ")")) %>% 
  select(Genotype, baseline.risk, fifteen.risk)


cin3pl.results[[2]]%>%  as.data.frame() %>% mutate_if(is.numeric, round, 1) %>% 
  mutate(baseline.risk = paste0(CR.0, " (", CR.lower95.0, " - ", CR.upper95.0, ")"), 
         fifteen.risk = paste0(CR.15, " (", CR.lower95.15, " - ", CR.upper95.15, ")")) %>% 
  select(Genotype, baseline.risk, fifteen.risk)

cancer.results.grouped[[2]]%>%  as.data.frame() %>% mutate_if(is.numeric, round, 1) %>% 
  mutate(baseline.risk = paste0(CR.0, " (", CR.lower95.0, " - ", CR.upper95.0, ")"), 
         fifteen.risk = paste0(CR.15, " (", CR.lower95.15, " - ", CR.upper95.15, ")")) %>% 
  select(Genotype, baseline.risk, fifteen.risk)
