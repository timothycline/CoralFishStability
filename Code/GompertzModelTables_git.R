rm(list=ls())

library(here)
SIG <- T

FishI <- readRDS(here('Output','GompertzResults','Coral_Gompertz_FishComm_InterannualFish.RDS')) %>% filter(Kbeta==1 & Bt==0)
#FishI_logit <- readRDS('~/Documents/Moorea/GompertzResults/Coral_Gompertz_FishComm_InterannualFish_Logit.RDS') %>% filter(Kbeta==1 & Bt==0)
FishMatI <- matrix(NA,nrow=3,ncol=5)
row.names(FishMatI) <- c('FO','BA','FR')
colnames(FishMatI) <- c('Bio.m','Rich','Div','Lmax','TL')
for(Hab in 1:3){
  for(Covar in 1:5){
    #Covar <-1;Hab<-1
    modline <- FishI %>% filter(CommAttr==switch(Covar,'Bio.m','Rich','FishBiomDiv','Lmax','MeanTL')) %>% pull(switch(Hab,Coef_FO,Coef_BA,Coef_FR))
    errline <- FishI %>% filter(CommAttr==switch(Covar,'Bio.m','Rich','FishBiomDiv','Lmax','MeanTL')) %>% pull(switch(Hab,Err_FO,Err_BA,Err_FR))
    errline <- ifelse(is.na(errline),1,errline)
    if(SIG){if((modline < 0 & modline+1.96*errline > 0) | (modline > 0 & modline-1.96*errline < 0)){modline <- NA}}
    FishMatI[Hab,Covar] <- modline
  }
}

TrophicGroupListI <- list()
TrophicGroupListI[[1]] <- FishMatI
for(i in 1:6){
  TG <- c(c('PrimCons','Herbivore','Grazer','ScraperExcavator','Browser','Corallivore'))[i]
  FishI <- readRDS(here('Output','GompertzResults',paste0('Coral_Gompertz_',TG,'_InterannualFish.RDS'))) %>% filter(Kbeta==1 & Bt==0)
  TrophicGroupMatI <- matrix(NA,nrow=3,ncol=4)
  row.names(TrophicGroupMatI) <- c('FO','BA','FR')
  colnames(TrophicGroupMatI) <- c('Bio.m','Rich','Div','Lmax')
  for(Hab in 1:3){
    for(Covar in 1:4){
      #Covar <-1;Hab<-1
      modline <- FishI %>% filter(CommAttr==switch(Covar,'Bio.m','Rich','FishBiomDiv','Lmax')) %>% pull(switch(Hab,Coef_FO,Coef_BA,Coef_FR))
      errline <- FishI %>% filter(CommAttr==switch(Covar,'Bio.m','Rich','FishBiomDiv','Lmax')) %>% pull(switch(Hab,Err_FO,Err_BA,Err_FR))
      errline <- ifelse(is.na(errline),1,errline)
      if(SIG){if((modline < 0 & modline+1.96*errline > 0) | (modline > 0 & modline-1.96*errline < 0)){modline <- NA}}
      TrophicGroupMatI[Hab,Covar] <- modline
    }
  }
  
  TrophicGroupListI[[i+1]] <- TrophicGroupMatI
}
names(TrophicGroupListI) <- c('FishComm','PrimCons','Herbivore','Grazer','Scraper','Browser','Corallivore')


####### MACRO
FishI <- readRDS(here('Output','GompertzResults','Macro_Gompertz_FishComm_InterannualFish.RDS')) %>% filter(Kbeta==1 & Bt==0)
#FishI_logit <- readRDS('~/Documents/Moorea/GompertzResults/Coral_Gompertz_FishComm_InterannualFish_Logit.RDS') %>% filter(Kbeta==1 & Bt==0)
FishMatI <- matrix(NA,nrow=3,ncol=5)
row.names(FishMatI) <- c('FO','BA','FR')
colnames(FishMatI) <- c('Bio.m','Rich','Div','Lmax','TL')
for(Hab in 1:3){
  for(Covar in 1:5){
    #Covar <-1;Hab<-1
    modline <- FishI %>% filter(CommAttr==switch(Covar,'Bio.m','Rich','FishBiomDiv','Lmax','MeanTL')) %>% pull(switch(Hab,Coef_FO,Coef_BA,Coef_FR))
    errline <- FishI %>% filter(CommAttr==switch(Covar,'Bio.m','Rich','FishBiomDiv','Lmax','MeanTL')) %>% pull(switch(Hab,Err_FO,Err_BA,Err_FR))
    errline <- ifelse(is.na(errline),1,errline)
    if(SIG){if((modline < 0 & modline+1.96*errline > 0) | (modline > 0 & modline-1.96*errline < 0)){modline <- NA}}
    FishMatI[Hab,Covar] <- modline
  }
}

Macro_TrophicGroupListI <- list()
Macro_TrophicGroupListI[[1]] <- FishMatI
for(i in 1:6){
  TG <- c(c('PrimCons','Herbivore','Grazer','ScraperExcavator','Browser','Corallivore'))[i]
  #FishR <- readRDS(paste0('~/Documents/Moorea/GompertzResults/Coral_Gompertz_',TG,'_FishOnR.RDS') %>% filter(Kbeta==1 & Bt==0)
  FishI <- readRDS(here('Output','GompertzResults',paste0('Macro_Gompertz_',TG,'_InterannualFish.RDS'))) %>% filter(Kbeta==1 & Bt==0)
  TrophicGroupMatI <- matrix(NA,nrow=3,ncol=4)
  row.names(TrophicGroupMatI) <- c('FO','BA','FR')
  colnames(TrophicGroupMatI) <- c('Bio.m','Rich','Div','Lmax')
  for(Hab in 1:3){
    for(Covar in 1:4){
      #Covar <-1;Hab<-1
      modline <- FishI %>% filter(CommAttr==switch(Covar,'Bio.m','Rich','FishBiomDiv','Lmax')) %>% pull(switch(Hab,Coef_FO,Coef_BA,Coef_FR))
      errline <- FishI %>% filter(CommAttr==switch(Covar,'Bio.m','Rich','FishBiomDiv','Lmax')) %>% pull(switch(Hab,Err_FO,Err_BA,Err_FR))
      errline <- ifelse(is.na(errline),1,errline)
      if(SIG){if((modline < 0 & modline+1.96*errline > 0) | (modline > 0 & modline-1.96*errline < 0)){modline <- NA}}
      TrophicGroupMatI[Hab,Covar] <- modline
    }
  }
  
  Macro_TrophicGroupListI[[i+1]] <- TrophicGroupMatI
}
names(Macro_TrophicGroupListI) <- c('FishComm','PrimCons','Herbivore','Grazer','Scraper','Browser','Corallivore')


save(list=ls(),file=here('Output','Revision_ModelTables_SIG.Rdata'))


