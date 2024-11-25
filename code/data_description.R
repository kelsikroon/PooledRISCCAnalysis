# data description 
rm(list=ls())
source("data_preparation.R")

genotypes <- paste0("hpv", c(16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59, 66, 68) )

# genotype summary  ---------------------
combo.dat.cin2.all <- create.combo.dat('cin2plus', trim.age = F)
combo.dat.cin3.all <- create.combo.dat('cin3plus', trim.age = F)
combo.dat.cancer.all <- create.combo.dat('cancer', trim.age = F)

genotypes.study <- combo.dat.cin2[,c("Source", genotypes)]

genotypes.study <- genotypes.study %>% group_by(Source) %>% summarise_all(funs(sum), na.rm=T) %>% 
  mutate(round(.[2:15]/rowSums(.[2:15])*100, 4)) %>% 
  pivot_longer(cols = `hpv16`:`hpv68`, names_to = "genotype",values_to = "perc")

genotypes.study$genotype <- factor(genotypes.study$genotype)

ggplot(genotypes.study, aes(x=genotype, y=perc, fill=Source) ) + theme_classic() +
  geom_bar(stat='identity', position = position_dodge(width=0.5), width=0.4) +
  xlab("Genotype") + ylab("Percentage (%)") + scale_fill_discrete("Study") + 
  theme(text = element_text(size=12), legend.position.inside=c(0.9, 0.7)) + 
  guides(fill = guide_legend(position='inside')) 


# age/ cytology information 
age_cyt.dat <- combo.dat.cin2.all[, c("Source", "age", "cyt")]
age_cyt.dat$Source <- factor(age_cyt.dat$Source, levels= c("pob", "vusa", "improve", "italy", "slov", "swed"))

age_cyt.dat %>% group_by(Source) %>% summarise(total = n(), age.med = median(age), age.range = paste0(round(range(age)[1]), " - ", round(range(age)[2])), 
                                               cyt.NILM = paste0(sum(cyt == 1), " (", round(sum(cyt==1)/total*100, 2), ")"),  
                                               cyt.ASCUS = paste0(sum(cyt > 1 & cyt <=3), " (", round(sum(cyt > 1 & cyt <=3)/total*100, 2), ")"), 
                                               cyt.HSIL = paste0(sum(cyt >= 4), " (", round(sum(cyt >= 4)/total*100, 2), ")"))  

# CIN2/CIN3/cancer cases
hist.info <- combo.dat.cin2.all %>% group_by(Source) %>% summarise(cin2pl = sum(right < Inf), N= n()) %>% as.data.frame()
hist.info$cin3pl <- combo.dat.cin3.all %>% group_by(Source) %>% summarise(cin3pl = sum(right < Inf)) %>% as.data.frame() %>% select(cin3pl) %>% unlist() %>% unname()
hist.info$cancer <- combo.dat.cancer.all %>% group_by(Source) %>% summarise(cancer = sum(right < Inf)) %>% as.data.frame() %>% select(cancer)%>% unlist() %>% unname()

hist.info %>% mutate(cin2 = cin2pl - cin3pl, cin3 = cin3pl - cancer, perc.cin2 = round(cin2/N*100, 2),perc.cin3 = round(cin3/N*100, 2),perc.ca = round(cancer/N*100, 2)) %>% 
  mutate(perc.cin2 = paste0(cin2, " (", perc.cin2, ")"),perc.cin3 = paste0(cin3, " (", perc.cin3, ")"), perc.ca = paste0(cancer, " (", perc.ca, ")")) %>% 
  select(c(Source, perc.cin2, perc.cin3, perc.ca)) %>% subset(Source == 'vusa')
