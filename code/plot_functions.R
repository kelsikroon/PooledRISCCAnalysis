# plot and tabl3 functions



genotype.dist <- function(endpoint.dat){
  endpoint.dat <- endpoint.dat[endpoint.dat$right < Inf, startsWith(colnames(endpoint.dat), 'hpv')][,1:14]
  endpoint.dat[is.na(endpoint.dat)] <- 0
  
  genotypes.mat <- matrix(NA, 13, 13)
  for (i in 1:13){
    for (j in (i+1):14){
      genotypes.mat[i, j-1] <-  sum(endpoint.dat[[i]] ==1 & endpoint.dat[[j]] == 1 & rowSums(endpoint.dat)==2)
    }
  }
  genotypes.mat <- t(genotypes.mat)
  genotypes.mat <- cbind(genotypes.mat, rep(NA, 13))
  genotypes.mat <- rbind(colSums(endpoint.dat[which(rowSums(endpoint.dat)==1),], na.rm=T), 
                         round(colSums(endpoint.dat[which(rowSums(endpoint.dat)==1),], na.rm=T)/dim(endpoint.dat[which(rowSums(endpoint.dat)==1),])[1]*100, 2),
                         colSums(endpoint.dat[which(rowSums(endpoint.dat)==2),], na.rm=T),
                         round(colSums(endpoint.dat[which(rowSums(endpoint.dat)==2),], na.rm=T)/dim(endpoint.dat[which(rowSums(endpoint.dat)==2),])[1]*100, 2),
                         genotypes.mat)
  genotypes.mat <- data.frame(genotypes.mat)
  colnames(genotypes.mat) <-  colnames(endpoint.dat)
  genotypes.mat$overall <- c(sum(genotypes.mat[1,]), round(sum(genotypes.mat[1,])/dim(endpoint.dat)[1]*100, 2), 
                             dim(endpoint.dat[which(rowSums(endpoint.dat)==2),])[1], round(dim(endpoint.dat[which(rowSums(endpoint.dat)==2),])[1]/dim(endpoint.dat)[1]*100, 2), 
                             rep(NA, 13))
  
  rownames(genotypes.mat) <- c("single", "perc.of.single", "double", "perc.of.double", colnames(endpoint.dat)[-1])
  remove.rows <- ifelse(rowSums(genotypes.mat, na.rm=T)==0, T, F)
  genotypes.mat <- genotypes.mat[!remove.rows, ]
  genotypes.mat[is.na(genotypes.mat)] <- '-'
  
  return(genotypes.mat)
}


# create function to make a prediction based on a model fit and then plot the cumulative risk for each genotype 
pred.and.plot <- function(fit, title, plot.ci=F){
  cov.genotypes <- fit$model[[1]] #c("hpv16", "hpv18", "hpv31", "hpv33", "hpv45", "hpv52", "hpv58")
  temp.data <- data.frame(diag(length(cov.genotypes)))  # create dummy data for plotting individual risk profiles
  temp.data <- rbind(rep(0, length(cov.genotypes)), temp.data)
  
  colnames(temp.data) <- cov.genotypes
  pred <- PICmodel.predict(c(cov.genotypes), c(), c(cov.genotypes), data=temp.data, 
                           time.points= seq(0, 15, 0.5), fit=fit, include.h=F, calc.CI = T)
  pred <- cbind(do.call(rbind, pred), Genotype=rep(c("hpvOther", cov.genotypes), each=31))
  if (any(pred$Genotype == 'hpv18.45')){
    pred$Genotype[pred$Genotype == 'hpv18.45'] <- 'hpv18/45'
  }
  pred.table <- pred[pred$Time %in% c(0, 1, 5, 15),] %>% 
    reshape( idvar = 'Genotype', timevar = "Time", direction='wide') %>% 
    mutate_if(is.numeric, function(x) round(x * 100, 2))
  ymax <- ifelse(title =="Cancer", 0.06, ifelse(title =='CIN2+', 0.46, 0.38))
  nudge_ylabel <- ifelse(title=='Cancer', 0.0001, 0.005)
  
  risk.plot <- ggplot(pred, aes(x=Time, y=CR, col=Genotype, group=Genotype)) + 
    theme_classic() + ylab(paste0("Cumulative ", title ," risk")) + xlab("Time since enrolment (years)") +
    geom_line(lwd=1) + theme(legend.position = 'none') + labs(title=title) + 
    scale_x_continuous(limits=c(0, 18), breaks=seq(0, 15, 3), expand=c(0,0)) +
    scale_y_continuous(limits=c(0, ymax), expand=c(0,0), labels = scales::percent_format(accuracy=1)) + 
    geom_text_repel(data=pred[pred$Time == 15,], aes(color=Genotype, label=Genotype), segment.colour =NA,
                    hjust=0, direction='y', xlim=c(15, NA), segment.size =0,  nudge_y=nudge_ylabel, box.padding = 0.085) +
    scale_color_manual(values = rcartocolor::carto_pal(n = 12, name = "Bold"))
  
  if (plot.ci) risk.plot <- risk.plot + geom_ribbon(aes(ymin=CR.lower95, ymax=CR.upper95, fill=Genotype), alpha=0.15) + 
    scale_fill_manual(values = rcartocolor::carto_pal(n = 12, name = "Bold"))
  
  return(list(risk.plot, pred.table))
}


pooled.fits.df <- function(fits, remove.interact=T){
  for (i in 1:length(fits)){
    fits[[i]][[1]]$summary$label <-  fits[[i]][[1]]$label
  }
  fits.summaries <- lapply(fits, function(x) x[[1]]$summary)
  
  fits.df <- do.call(rbind, fits.summaries)
  
  fits.df$group <- case_when(
    str_detect(rownames(fits.df), "g") ~ "prog",
    str_detect(rownames(fits.df), "p") ~ "prev",
    str_detect(rownames(fits.df), "w") ~ "clear"
  )
  fits.df$label <- factor(fits.df$label, levels=unique(fits.df$label))
  return(fits.df)
}

#!str_detect(data$param, "\\.")
library(ggtext)
library(glue)

highlight = function(x, pat, color="black", family="") {
  ifelse(grepl(pat, x), glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
}

genotype.forest.plot <- function(data, legend.lab='Model', dodge.height=0.6){
  # for progression 
  prog.genotypes <- ggplot(data[data$group %in% c("prog")   ,], 
                           aes(y = param, col=label, label=label)) +  
    theme_classic() + ylab("") + xlab("Relative risk (95% CI)") + 
    labs(title=expression('Effect of genotype on log-progression rate'~(lambda[1]))) +
         #caption = '\n\n') + 
    scale_y_discrete(limits=rev, labels = function(x) highlight(x, "hpv16|hpv18\\/45|hpvGroup1Alpha9")) +  
    scale_x_log10() + 
    geom_vline(xintercept = 1, lty=2, col='grey') + 
    geom_hline(yintercept = 7.5, lty=3, col='lightgrey') + 
    
    #geom_hline(yintercept = c(5.5, 7.5), lty=2, col='lightgrey') + 
    geom_point(aes(x=exp(theta.hat)), position=position_dodge(width=dodge.height), shape=16, size=2.5) +
    geom_linerange(aes(xmin=exp(lower), xmax=exp(upper)), position=position_dodge(width=dodge.height)) + 
    scale_color_manual(name = legend.lab, values = rcartocolor::carto_pal(n = 7, name='Bold')[c(5, 2, 3)]) +
    guides(col = guide_legend(byrow = T)) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
          legend.spacing.y = unit(1, "cm"), axis.text.y = element_markdown(),
          legend.text = element_text(margin=margin(t=5, b=5, unit='pt'), hjust=0), legend.position = "inside", legend.position.inside = c(1, 6)) 
  
  
  # for prevalence -  prevalence genotype forest plot
  prev.genotypes <- ggplot(data[data$group %in% c("prev") ,], aes(y = param, col=label)) +  
    theme_classic() + ylab("") + xlab("Odds ratio (95% CI)") + 
    labs(title=expression('Effect of genotype on prevalent proportion'~(pi))) +
    #labs(caption= 'HPV other = HPV 35/39/51/56/59/66/68-positive\nPrevalence parameter has only genotype effects and study intercepts, unless specified') + 
    scale_y_discrete(limits=rev, labels = function(x) highlight(x, "hpv16|hpv18\\/45|hpvGroup1Alpha9"))  +
    scale_x_continuous(trans='log10') + 
    geom_vline(xintercept = 1, lty=2, col='grey') + 
    geom_hline(yintercept = 7.5, lty=3, col='lightgrey') + 
    #geom_hline(yintercept = c(5.5, 7.5), lty=2, col='lightgrey') + 
    geom_point(aes(x=exp(theta.hat)), position=position_dodge(width=dodge.height), shape=16, size=2.5) +
    geom_linerange(aes(xmin=exp(lower), xmax=exp(upper)), position=position_dodge(width=dodge.height)) + 
    scale_color_manual(name = legend.lab, values= rcartocolor::carto_pal(n = 7, name='Bold')[c(5, 2, 3)]) +
    guides(col = guide_legend(byrow = T)) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
          legend.text.align = 0, legend.spacing.y = unit(1, "cm"), axis.text.y = element_markdown(),
          legend.text = element_text(margin=margin(t=5, b=5, unit='pt')),  legend.position = "inside", 
          legend.position.inside = c(0.9, 0.22)) 
  
  #1300 x 700
  pairplot <- ggarrange(prog.genotypes, prev.genotypes, common.legend=F, nrow=1, labels=c("A", "B"))
  return(pairplot)
}


# forest plot of parameter estimates ------------
pooled.fits.endpoint <- function(file.date){
  files <- list.files(pattern=paste0("fit_", file.date))
  fits <- lapply(files, function(x) mget(load(x)))
  
  for (i in 1:length(fits)){
    fits[[i]][[1]]$summary$label <-  fits[[i]][[1]]$label #str_sub(files[i], nchar(paste0), nchar(files[i]) - 6)
  }
  fits.summaries <- lapply(fits, function(x) x[[1]]$summary)
  fits.df <- do.call(rbind, fits.summaries)
  
  fits.df$group <- case_when(
    str_detect(rownames(fits.df), "g") ~ "prog",
    str_detect(rownames(fits.df), "p") ~ "prev",
    str_detect(rownames(fits.df), "w") ~ "clear"
  )
  fits.df$label <- factor(case_when(
    fits.df$label == 'cancer' ~'Cancer',
    fits.df$label == 'cin2pl' ~ 'CIN2+',
    fits.df$label == 'cin3pl' ~ 'CIN3+'
  ), levels=c("CIN2+", "CIN3+", "Cancer"))
  return(fits.df)
}



