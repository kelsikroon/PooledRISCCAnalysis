# data preparation/organisation functions 

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

# RData is quicker to read than CSV files so run this first (only need to do it once)
save.RISCC.data <- function(study){
  folder <- paste0("~/Projects/RISCC Comparison/Analysis/Data RISCC/RISCC_DB_", study, "/", study, "_")
  patient_info <- read.csv(paste0(folder, "Participants_2024.csv"), sep=';')
  hpv <- read.csv(paste0(folder, "HPV_2024.csv"), sep=';')
  cytLong <- read.csv(paste0(folder, "Cytology_2024.csv"), sep=';')
  histLong <- read.csv(paste0(folder, "Histology_2024.csv"), sep=';')
  save(patient_info, hpv, cytLong, histLong, file= paste0('~/Projects/RISCC Comparison/Analysis/Data RISCC/RISCC_DB_', study, "_all.RData"))
}
# study.names <- c('POBASCAM', 'VUSAscreen', 'IMPROVE', 'NTCCitaly', 'ULslovenia', 'SwedeScreen' )
# for (study in study.names) save.RISCC.data(study)

data.descript <- function(study){
  load(file=paste0('~/Projects/RISCC Comparison/Analysis/Data RISCC/RISCC_DB_', study, "_all.RData"))
  
  if (study =='IMPROVE'){ # only include those with follow-up permission for the IMPROVE study
    load("~/Projects/RISCC Comparison/Analysis/improve_fup.RData")
    # only need to filter HPV data frame because next we select HPV positives from all data so 
    # this selection will follow-on through the rest of the data
    hpv <- hpv[hpv$Idwoman %in% improve.fup, ]
  }
  
  # prepare HPV positives
  hpv <- hpv[hpv$ResHPV1 == 1 & hpv$VisitHPV=='A' & hpv$Intyp==1, ] # select HPV positives with genotyping information
  if (study == 'NTCCitaly'){
    hpv <- hpv[hpv$Idwoman %in% patient_info$Idwoman[patient_info$Arm==2],] # for italy only take intervention arm
  }
  hpv$DateHPV <-  as.Date(hpv$DateHPV, format='%d-%m-%Y')
  hpv <- hpv[order(hpv$Idwoman),]
  
  # prepare patient info
  patient_info <- patient_info[trimws(patient_info$Idwoman) %in% trimws(hpv$Idwoman),] # select hpv positives from patient info
  rownames(patient_info) <- patient_info$Idwoman
  patient_info$Darecl <- as.Date(patient_info$Darecl, format='%d-%m-%Y')
  patient_info <- patient_info[order(patient_info$Idwoman),]
  
  # prepare histology data
  histLong <- histLong[histLong$Idwoman %in% hpv$Idwoman, c("Idwoman", "VisitHisto", "Datehisto", "Classhisto", "Typehisto")]
  colnames(histLong) <- c("id", "visit", "date", "res", "type")
  histLong$date <- as.Date(histLong$date, format='%d-%m-%Y')
  histLong <- histLong[order(histLong$id, histLong$date),]
  
  histLong$cin2plus <- ifelse(histLong$res >= 7 & histLong$res!=13, 1, 0)
  histLong$cin3plus <-ifelse(histLong$res >= 8 & histLong$res!=13, 1, 0)
  histLong$cancer <- ifelse(histLong$res >= 9 & histLong$res!=13, 1, 0)
  
  histLong$ue.cin2plus <- ifelse(!is.na(histLong$type) & histLong$type==4 & !histLong$cin2plus, 1, 0)
  histLong$ue.cin3plus <- ifelse(!is.na(histLong$type) & histLong$type==4 & !histLong$cin3plus, 1, 0)
  histLong$ue.cancer <- ifelse(!is.na(histLong$type) & histLong$type==4 & !histLong$cancer, 1, 0)
  
  rownames(histLong) <- NULL
  histLong <- histLong %>% tidyr::drop_na(res) # remove rows with NA as a histology result result
  
  # prepare cytology data
  cytLong <- cytLong[cytLong$Idwoman %in% hpv$Idwoman, c("Idwoman", "VisitCyto", "Datecyto", "Classcyto")] # subset HPV-positives
  colnames(cytLong) <- c("id", "visit", "date", "res")
  #cytLong <- na.omit(cytLong) # remove NA values 
  cytLong$date <- as.Date(cytLong$date, format='%d-%m-%Y')
  cytLong <- cytLong[order(cytLong$id, cytLong$date),] # order by ID and date
  baseline_cyt <- cytLong[cytLong$visit=='A',]$res # store baseline cytology results
  
  # prepare genotype data
  genotype_col_num <-  case_when(
    study == 'NTCCitaly' ~ "_2",
    study =='ULslovenia' ~ "_new",
    T ~ '_1'
  )
  # for ULslovenia the situation is a little bit different because the HPV genotype test number changed for different patients, 
  # so we merge the 3 columns that have genotyping results for different patients  
  if (study == 'ULslovenia'){
    genotypes <- paste0("HPV", c(16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59, 66, 68))
    for (i in 1:14){
      hpv[[paste0(genotypes[i], "_new")]] <- apply(hpv[, grep(genotypes[i], names(hpv))], 1, max, na.rm=T)
    }
  }
  genotype_cols <- paste0("HPV", c(16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59, 66, 68), genotype_col_num)
  genotype_data <- hpv[, genotype_cols]
  genotype_data[genotype_data ==9] <- NA # 9 is the code for "missing" genotype in the data so replace with NA
  colnames(genotype_data) <- paste0("hpv", as.character(c(16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59, 66, 68)))
  
  worst_histo <- rep(NA, length(hpv$Idwoman))
  baseline_cyts <-  rep(NA, length(hpv$Idwoman))
  for (i in 1:length(hpv$Idwoman)){ # loop through all patient ids
    baseline <- patient_info$Darecl[patient_info$Idwoman == hpv$Idwoman[i]] # store baseline date for patient i 
    baseline_cyts[i] <- ifelse(!is.na(baseline_cyt[i]), baseline_cyt[i], 99)
    histo <- histLong[histLong$id == hpv$Idwoman[i] & histLong$date >= baseline,] # create subset of histology results for patient i 
    worst_histo[i] <- ifelse(!is.na(max(histo$res)), max(histo$res), 99)
  }
  df <- cbind(data.frame(id= hpv$Idwoman, worst_histo, baseline_cyts), genotype_data, study=study)
  return(df)
}

# six studies are:
# (1) POBASCAM
# (2) VUSAscreen
# (3) IMPROVE
# (4) NTCCitaly
# (5) ULslovenia
# (6) SwedeScreen
  
# endpoints are:
# (1) 'cin2plus'
# (2) 'cin3plus'
# (3) 'cancer'

# interval.option is either:
# (1) 'hierarchical' meaning histology results are prioritised over cytology resutlts for left interval 
# (2) 'all' meaning any most recent result is taken regardless of whether it is histology/cytology 

# function to create interval data from RISCC database files 
create.interval.data <- function(study, endpoint, interval.option = 'hierarchical'){
  load(file=paste0('~/Projects/RISCC Comparison/Analysis/Data RISCC/RISCC_DB_', study, "_all.RData"))
  
  if (study =='IMPROVE'){ # only include those with follow-up permission for the IMPROVE study
    load("~/Projects/RISCC Comparison/Analysis/improve_fup.RData")
    # only need to filter HPV data frame because next we select HPV positives from all data so 
    # this selection will follow-on through the rest of the data
    hpv <- hpv[hpv$Idwoman %in% improve.fup, ]
  }
  
  # prepare HPV positives
  hpv <- hpv[hpv$ResHPV1 == 1 & hpv$VisitHPV=='A' & hpv$Intyp==1, ] # select HPV positives with genotyping information
  if (study == 'NTCCitaly'){
    hpv <- hpv[hpv$Idwoman %in% patient_info$Idwoman[patient_info$Arm==2],] # for italy only take intervention arm
  }
  hpv$DateHPV <-  as.Date(hpv$DateHPV, format='%d-%m-%Y')
  hpv <- hpv[order(hpv$Idwoman),]
  
  # prepare patient info
  patient_info <- patient_info[trimws(patient_info$Idwoman) %in% trimws(hpv$Idwoman),] # select hpv positives from patient info
  rownames(patient_info) <- patient_info$Idwoman
  patient_info$Darecl <- as.Date(patient_info$Darecl, format='%d-%m-%Y')
  patient_info <- patient_info[order(patient_info$Idwoman),]
  
  # prepare histology data
  histLong <- histLong[histLong$Idwoman %in% hpv$Idwoman, c("Idwoman", "VisitHisto", "Datehisto", "Classhisto", "Typehisto")]
  colnames(histLong) <- c("id", "visit", "date", "res", "type")
  histLong$date <- as.Date(histLong$date, format='%d-%m-%Y')
  histLong <- histLong[order(histLong$id, histLong$date),]
  
  histLong$cin2plus <- ifelse(histLong$res >= 7 & histLong$res!=13, 1, 0)
  histLong$cin3plus <-ifelse(histLong$res >= 8 & histLong$res!=13, 1, 0)
  histLong$cancer <- ifelse(histLong$res >= 9 & histLong$res!=13, 1, 0)
  
  histLong$ue.cin2plus <- ifelse(!is.na(histLong$type) & histLong$type==4 & !histLong$cin2plus, 1, 0)
  histLong$ue.cin3plus <- ifelse(!is.na(histLong$type) & histLong$type==4 & !histLong$cin3plus, 1, 0)
  histLong$ue.cancer <- ifelse(!is.na(histLong$type) & histLong$type==4 & !histLong$cancer, 1, 0)
  
  rownames(histLong) <- NULL
  histLong <- histLong %>% tidyr::drop_na(res) # remove rows with NA as a histology result result
  
  # prepare cytology data
  cytLong <- cytLong[cytLong$Idwoman %in% hpv$Idwoman, c("Idwoman", "VisitCyto", "Datecyto", "Classcyto")] # subset HPV-positives
  colnames(cytLong) <- c("id", "visit", "date", "res")
  #cytLong <- na.omit(cytLong) # remove NA values 
  cytLong$date <- as.Date(cytLong$date, format='%d-%m-%Y')
  cytLong <- cytLong[order(cytLong$id, cytLong$date),] # order by ID and date
  
  # prepare genotype data
  genotype_col_num <-  case_when(
    study == 'NTCCitaly' ~ "_2",
    study =='ULslovenia' ~ "_new",
    T ~ '_1'
  )
  # for ULslovenia the situation is a little bit different because the HPV genotype test number changed for different patients, 
  # so we merge the 3 columns that have genotyping results for different patients  
  if (study == 'ULslovenia'){
    genotypes <- paste0("HPV", c(16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59, 66, 68))
    for (i in 1:14){
      hpv[[paste0(genotypes[i], "_new")]] <- apply(hpv[, grep(genotypes[i], names(hpv))], 1, max, na.rm=T)
    }
  }
  genotype_cols <- paste0("HPV", c(16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59, 66, 68), genotype_col_num)
  genotype_data <- hpv[, genotype_cols]
  genotype_data[genotype_data ==9] <- NA # 9 is the code for "missing" genotype in the data so replace with NA
  colnames(genotype_data) <- paste0("hpv", as.character(c(16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59, 66, 68)))
  
  # prepare censoring info
  cens_dates <- unname(sapply(patient_info[patient_info$Idwoman %in% hpv$Idwoman,]$DaCens, function(x) as.Date(paste("01-",x,sep=""), "%d-%m-%Y")))
  cens_dates <- as.Date(cens_dates, origin='1970-01-01') # store censoring dates
  censoring <- patient_info[patient_info$Idwoman %in% hpv$Idwoman,]$StCens # store censoring status 
  
  # start creating intervals ----
  baseline_cyt <- cytLong[cytLong$visit=='A',]$res # store baseline cytology results
  names(baseline_cyt) <- cytLong[cytLong$visit=='A',]$id # label baseline cytology results with patient ID

  hpv_pos <- hpv[hpv$Idwoman %in% names(baseline_cyt),] # only keep those who have both baseline cyt and baseline hpv
  ids <- unique(hpv_pos$Idwoman) # create a vector of ids for each woman
  right <- rep(NA, length(ids)) # create empty vector of left and right interval
  left <- rep(NA, length(ids))
  worst_histo <- rep(NA, length(ids))
  for (i in 1:length(ids)){ # loop through all patient ids
    baseline <- patient_info$Darecl[patient_info$Idwoman == ids[i]] # store baseline date for patient i 
    histo <- histLong[histLong$id == ids[i] & histLong$date >= baseline,] # create subset of histology results for patient i 
    cyto <- cytLong[cytLong$id==ids[i] & cytLong$date >= baseline,] # create subset of cytology results for patient i 
    worst_histo[i] <- ifelse(!is.na(max(histo$res)), max(histo$res), 99)
    
    # if there are no histology results, then left interval is last Pap1 date and right interval is infinity
    if (nrow(histo)==0){
      cyt_index <- length(cyto$res)-match(1,rev(cyto$res))+1 # index of last Pap1 cytology result
      cytodate <- as.Date(ifelse(is.na(cyt_index), baseline, cyto$date[cyt_index]), origin='1970-01-01') # date of last Pap1
      left[i] <- as.numeric(difftime(cytodate, baseline)) # time in days since baseline 
      right[i] <- Inf # right is infinity because they did not develop CIN2/3/ca in follow-up
      next
    }
    
    ue.endpoint <- paste0('ue.', endpoint)
    # if they had a UE or metastasis results then make their interval [histology date, Inf)
    # ---> only if they have UE *BEFORE* CIN2/3/ca result, if they have UE *AFTER* CIN2/3/ca result then the CIN2/3/ca needs to be the right interval
    if (sum(histo[[ue.endpoint]]) !=0 ){
      # UE indicator is 1 for CIN2/3/ca and 2 for UE/metastasis so a patients results might look like:
      # 1) 000112 which means they got CIN2/3/ca before UE --> we don't censor at UE time
      # 2) 000200 which means they got UE --> so we censor at the time UE is detected because they cannot develop CIN2/3/ca now
      ue_indicator <- ifelse(histo[[ue.endpoint]], 2, ifelse(histo[[endpoint]], 1, 0)) # + histo[[endpoint]]

      # next: check if there are any CIN2/3/ca results and if CIN2/3/ca comes before UE
      # if CIN2/3/ca comes after UE then we censor at the detection of UE (otherwise we continue to find the date of the CIN2/3/ca)
      ue_first <- ifelse(sum(histo[[endpoint]]) !=0 & (match(1, ue_indicator) > match(2, ue_indicator)) | all(ue_indicator == 2), 1, 0)

      if (sum(histo[[endpoint]]) ==0 | ue_first ) { #str_detect(paste(ue_indicator, collapse=''), '12', negate=T)){
        ue_index <- min(match(1, histo[[ue.endpoint]]), na.rm=TRUE) # index when a UE result was first recorded
        ue_date <- as.Date(ifelse(is.na(ue_index), baseline, histo$date[ue_index]), origin='1970-01-01')
        left[i] <- as.numeric(difftime(ue_date, baseline))
        right[i]<- Inf # right=Inf because they wont develop CIN2/3/ca during follow-up
        next
      }
    }
    
    # if all histology results are 0 then the endpoint did not develop during follow-up so left interval is the max date 
    # between last <endpoint result or last Pap 1 result, and the right interval is infinity
    if (sum(histo[[endpoint]]) ==0){
      histodate <- histo$date[length(histo$date)]
      cytodate <- cyto$date[length(cyto$date)-match(1,rev(cyto$res))+1] # last pap1 cyt date
      if (is.na(cytodate)) {
        left[i] <- difftime(histodate, baseline)
      }else{
        left[i] <- difftime(as.Date(max(histodate, cytodate)), baseline)
      }
      right[i] <- Inf 
      next
    }
    
    # we find the index of the first time the histology result is CIN2/3/ca and then save the save of the first CIN2/3/ca result
    right_index <- match(1, histo[[endpoint]])
    first_cin <- histo$date[right_index]
    r_int <- NULL
    # check if there was a cytology result >Pap 1 within 3 months before the first CIN2/3/ca result
    # as this would be the index smear that led to histology referral
    if (as.numeric(difftime(first_cin, max(cyto[cyto$date<first_cin & cyto$res>=2 & cyto$res!=99,]$date))) < 90){
      r_int <-  max(cyto[cyto$date<first_cin &  cyto$res>=2 & cyto$res!=99,]$date) # index smear date 
    }else{
      r_int <- first_cin # histology date 
    }
    right[i] <- as.numeric(difftime(r_int, baseline, units='days'))
    
    # left intervals: 
    cyto_temp <- cyto[cyto$date < r_int-1,] # filter cytology date before right interval
    histo_temp <- histo[histo$date < r_int-1,] # filter histology date before right interval
    
    if (nrow(histo_temp)==0){
      # if there are no histology results, then the left interval is the last Pap 1 cytology result 
      cyt_index <- length(cyto_temp$res)-match(1,rev(cyto_temp$res))+1 # index of last Pap1 cytology result
      cytodate <- as.Date(ifelse(is.na(cyt_index), baseline, cyto_temp$date[cyt_index]), origin='1970-01-01') # date of last Pap1
      left[i] <- as.numeric(difftime(cytodate, baseline)) # time in days since baseline 
      next
    }else{
      # otherwise if histology is available, then the left interval is either the last Pap 1 cytology result or the last histology less than endpoint
      histodate <- histo_temp$date[length(histo_temp$date)]
      cytodate <- cyto_temp$date[length(cyto_temp$date)-match(1,rev(cyto_temp$res))+1] # last pap1 cyt date
      if (interval.option =='hierarchical' | is.na(cytodate)){
        # hierarchical means that if the results are CIN1, Pap 1, CIN3, then interval = [CIN1 date, CIN3 date] and Pap1 date is ignored
        left[i] <- difftime(histodate, baseline)
      } else{
        # if interval is not hierarchical then we just take the most recent date 
        left[i] <- difftime(as.Date(max(histodate, cytodate)), baseline)
      }
      next
    }
    
    # - check which women are censored and check the date of censoring:
    #     if the date is after a right interval then ignore the censoring
    #     if the date is before a right interval, then censor becomes the left interval and right=Inf (but check with Hans)
    if (censoring[i] %in% c(1, 2)){ # if death or migration happened
      censdate <- difftime(cens_dates[i], baseline)
      if(right[i]==Inf & left[i]< censdate){
        left[i] <- censdate
        next
      }
    }
  }
  # put it all together:
  new_data <- data.frame(left=left/365, right= right/365)
  rownames(new_data) <- ids
  
  new_data$right[new_data$right < 0.25] <- 0 # set as prevalent if within the first 3 months
  new_data$left[new_data$right ==0] <- 0
  new_data$z <- ifelse(new_data$right==0, 1, ifelse(new_data$left==0, NA, 0))
  new_data$age <- as.numeric(patient_info[patient_info$Idwoman %in% ids,]$Age)
  
  new_data$cyt_bmd <- ifelse(baseline_cyt > "1" & baseline_cyt < "4", 1, 0)
  new_data$cyt_severe <- ifelse(baseline_cyt >= "4", 1, 0)
  new_data$cyt_abnormal <- ifelse(new_data$cyt_bmd ==1 | new_data$cyt_severe==1, 1, 0)
  
  new_data <- cbind(new_data, genotype_data[hpv$Idwoman %in% names(baseline_cyt),])
  new_data$cyt <- baseline_cyt
  new_data$worst_histo <- worst_histo
  #new_data <- new_data[new_data$left != 0 | new_data$right < Inf,]
  new_data <- new_data[!new_data$cyt %in% c(0, 99, "") & !is.na(new_data$cyt),] # remove patients with unsatisfactory or unavailable baseline cytology
  new_data$left[new_data$left < 0 ] <-0 # left intervals that are less than zero should not happen --> set to baseline
  
  new_data <- new_data %>% tidyr::drop_na(hpv16) # if hpv16 genotype is missing then the rest are also missing --> remove it
  return(new_data)
}




# create combined data  for a specific endpoint ---------------------
create.combo.dat <- function(endpoint, trim.age=T, interval.type='hierarchical'){
  # load functions and data ---------------------
  pob.dat <- create.interval.data("POBASCAM", endpoint, interval.type)
  #vusa.dat <- create.interval.data("VUSAscreen", endpoint, interval.type)
  load("Archive/Data Preparation/082024 Intervals/vusa_int.RData")
  improve.dat <- create.interval.data("IMPROVE", endpoint, interval.type)
  italy.dat <- create.interval.data("NTCCitaly", endpoint, interval.type)
  slov.dat <- create.interval.data("ULslovenia", endpoint, interval.type)
  swed.dat <- create.interval.data("SwedeScreen", endpoint, interval.type)

  vusa_cin2pl$cyt_bmd <- ifelse(vusa_cin2pl$cyt == 'BMD', 1, 0)
  vusa_cin2pl$cyt_severe <- ifelse(vusa_cin2pl$cyt == '>BMD', 1, 0)
  vusa_cin2pl$cyt_abnormal <- ifelse(vusa_cin2pl$cyt != 'pap1', 1, 0)

  vusa_cin3pl$cyt_bmd <- ifelse(vusa_cin3pl$cyt == 'BMD', 1, 0)
  vusa_cin3pl$cyt_severe <- ifelse(vusa_cin3pl$cyt == '>BMD', 1, 0)
  vusa_cin3pl$cyt_abnormal <- ifelse(vusa_cin3pl$cyt != 'pap1', 1, 0)

  vusa_ca$cyt_bmd <- ifelse(vusa_ca$cyt == 'BMD', 1, 0)
  vusa_ca$cyt_severe <- ifelse(vusa_ca$cyt == '>BMD', 1, 0)
  vusa_ca$cyt_abnormal <- ifelse(vusa_ca$cyt != 'pap1', 1, 0)

  vusa_cin2pl$cyt <- ifelse(vusa_cin2pl$cyt == 'BMD', 2, ifelse(vusa_cin2pl$cyt =='>BMD', 4, 1))
  vusa_cin3pl$cyt <- ifelse(vusa_cin3pl$cyt == 'BMD', 2, ifelse(vusa_cin3pl$cyt == '>BMD', 4, 1))
  vusa_ca$cyt <- ifelse(vusa_ca$cyt == 'BMD', 2, ifelse(vusa_ca$cyt == '>BMD', 4, 1))

  vusa_cin2plus <- vusa_cin2pl
  vusa_cin3plus <- vusa_cin3pl
  vusa_cancer <- vusa_ca
  # data preparation ---------------------
  genotype.cols <- paste0("hpv", c(16, 18, 31, 33, 35, 39, 45, 51, 52, 56, 58, 59, 66, 68))
  col.select <- c("left", "right", "z", "age", "cyt_bmd", "cyt_severe", "cyt", "cyt_abnormal", genotype.cols)
  datas <- c("pob", "vusa", "improve", "italy", "slov", "swed")
  dfnames <- c("pob.dat", paste0("vusa_", endpoint), paste0(datas[3:6], ".dat"))
  combo.dat <- do.call(rbind, lapply(1:6, function(x) cbind(get(dfnames[x])[,col.select], Source= datas[x]))) #substr(x, 0, nchar(x) - 4))))
  combo.dat$Source <- factor(combo.dat$Source, levels =datas)
  combo.dat$hpvOther <- with(combo.dat, ifelse(hpv35 | hpv39 | hpv51 | hpv56 | hpv59 | (hpv66 & !is.na(hpv66)) | hpv68, 1, 0))
  
  if (trim.age) combo.dat <- combo.dat[combo.dat$age >=29, ]
   
  combo.dat <- combo.dat %>% group_by(Source) %>% mutate(age.std =  0.5*(age - mean(age))/sd(age)) # add age.std for each study individually 
  combo.dat$age.std.joint <- 0.5*(combo.dat$age - mean(combo.dat$age))/sd(combo.dat$age) # add overall age.std
  
  combo.dat$z[combo.dat$right == Inf & combo.dat$left ==0 & is.na(combo.dat$z)] <- 0
  
  combo.dat$group <- factor(case_when(
    is.na(combo.dat$z) ~ 'Unknown',
    combo.dat$z == 1 ~ "Prevalent",
    combo.dat$right <Inf & combo.dat$right >0 ~ "Incident",
    T ~ "No event"
  ), levels = c("Prevalent", "Incident", 'Unknown', "No event"))
  
  combo.dat$hpv18.45 <- ifelse(combo.dat$hpv18 | combo.dat$hpv45, 1, 0)
  combo.dat$hpvGroup1Alpha9 <- with(combo.dat, ifelse(hpv31 | hpv33 | hpv35 | hpv52 | hpv58, 1, 0))
  
  return(combo.dat)
}
