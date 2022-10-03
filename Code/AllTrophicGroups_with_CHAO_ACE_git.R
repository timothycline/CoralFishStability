rm(list=ls())

#load required libraries
library(RColorBrewer)
library(MuMIn)
library(dplyr)
library(glmmTMB)
library(stringr)
library(ecofolio)
library(here)

#Some custom functions used below
CV<-function(x){
  return(sd(x,na.rm=T)/mean(x,na.rm=T))
}
logit<-function(x){
  return(log(x/(1-x)))
}
ilogit<-function(x){
  return(exp(x)/(1+exp(x)))
}


#Read in required datasets
for(TrophicGroup in c('FishComm','PrimCons','Herbivore','ScraperExcavator','Grazer','Browser','Corallivore')){
  #TrophicGroup <- 'Corallivore'
  load(here('Data',paste0('TSregTable_',TrophicGroup,'.Rdata')))
  
  #missing values are just 0's in the dataset (i.e., surveys were recorded but that group of fish or benthos not recorded)
  TSregTable$MacroCover[is.na(TSregTable$MacroCover)] <- 0 #no macro observed
  TSregTable <- TSregTable %>% filter(Year < 2018 & Year > 2005)
  TSregTable$Bio.m[is.na(TSregTable$Bio.m) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  TSregTable$Rich[is.na(TSregTable$Rich) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  TSregTable$CHAO[is.na(TSregTable$CHAO) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  TSregTable$ACE[is.na(TSregTable$ACE) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  TSregTable$FishBiomDiv[is.na(TSregTable$FishBiomDiv) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  TSregTable$MeanTL[is.na(TSregTable$MeanTL) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  TSregTable$Lmax[is.na(TSregTable$Lmax) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  

  TSregTable <- TSregTable %>% mutate(HabSite = paste(Habitat,Site,sep='.')) 
  ObservedSites <- TSregTable %>% group_by(HabSite) %>% summarize(ObservedBool = (sum(Bio.m)>0 & sum(Bio.m>0)>=2)) %>% filter(ObservedBool) %>% pull(HabSite)
  TSregTable <- TSregTable %>% filter(HabSite %in% ObservedSites)
  
  #Create two tables for regression analysis
  CVtab <- TSregTable %>% group_by(Habitat,Site) %>% 
    summarize(SDCoral = log(sd(logit(CoralCover+0.01),na.rm=T)),SDMacro=log(sd(logit(MacroCover+0.01),na.rm=T)),SDbiomass=log(sd(Bio.m+0.1)),CVbiomass=log(CV(Bio.m+0.1)),CVrich=CV(Rich),CVdiv=CV(FishBiomDiv),CVlmax=CV(Lmax),CVtl=CV(MeanTL),Source=unique(Source)) %>%
    mutate(HabSite=paste(Habitat,Site,sep='.'))
  #CVtab$CVbiomass
  
  Mtab <- TSregTable %>% group_by(Habitat,Site) %>% 
    summarize(MCoral = mean(logit(CoralCover+0.01),na.rm=T),MMacro = mean(logit(MacroCover+0.01),na.rm=T),Mbiomass=mean(log(Bio.m+0.1),na.rm=T),
              Mrich=mean(Rich,na.rm=T),Mchao = mean(CHAO,na.rm=T),Mace=mean(ACE,na.rm=T),
              Mdiv=mean(FishBiomDiv,na.rm=T),Mlmax=mean(Lmax,na.rm=T),Mtl=mean(MeanTL,na.rm=T),Source=unique(Source)) %>%
    mutate(HabSite=paste(Habitat,Site,sep='.'))
  
  #ModelSelection for Fish diversity stability
  lm.base<-glmmTMB(CVbiomass ~ Habitat + (1|Site), data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mrich),by=c('Habitat','Site')),REML=FALSE)
  lm.base.source<-glmmTMB(CVbiomass ~ Source + Habitat+ (1|Site), data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mrich),by=c('Habitat','Site')),REML=FALSE)
  AICc(lm.base);AICc(lm.base.source)
  summary(lm.base.source)
  r.squaredGLMM(lm.base.source)
  
  if(AICc(lm.base)>AICc(lm.base.source)){
    FishStability.Source <- TRUE
    FishStability.Base <- lm.base.source
  }else{
    FishStability.Source <- FALSE
    FishStability.Base <- lm.base
  }
 
  if(FishStability.Source){
    
    lm.rich1<-glmmTMB(CVbiomass ~ Mchao + Habitat + Source + (1|Site), data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mchao),by=c('Habitat','Site')),REML=FALSE)
    lm.rich2<-glmmTMB(CVbiomass ~ Mchao*Habitat + Source + (1|Site), data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mchao),by=c('Habitat','Site')),REML=FALSE)
    AICc(FishStability.Base);AICc(lm.rich1);AICc(lm.rich2)
    dAIC_Mrich <- round(min(AICc(lm.rich2)-AICc(FishStability.Base),AICc(lm.rich1)-AICc(FishStability.Base)),2)
    
    lm.div1<-glmmTMB(CVbiomass ~ Mdiv + Habitat + Source + (1|Site), data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mdiv),by=c('Habitat','Site')),REML=FALSE)
    lm.div2<-glmmTMB(CVbiomass ~ Mdiv * Habitat + Source + (1|Site), data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mdiv),by=c('Habitat','Site')),REML=FALSE)
    AICc(FishStability.Base);AICc(lm.div1);AICc(lm.div2)
    dAIC_Mdiv <- round(min(AICc(lm.div1)-AICc(FishStability.Base),AICc(lm.div2)-AICc(FishStability.Base)),2)
    
    lm.tl1<-glmmTMB(CVbiomass ~ Mtl + Habitat + Source + (1|Site), data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mtl),by=c('Habitat','Site')),REML=FALSE)
    lm.tl2<-glmmTMB(CVbiomass ~ Mtl*Habitat + Source + (1|Site), data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mtl),by=c('Habitat','Site')),REML=FALSE)
    AICc(FishStability.Base);AICc(lm.tl1);AICc(lm.tl2) #dAIC
    dAIC_Mtl <- round(min(AICc(lm.tl1)-AICc(FishStability.Base),AICc(lm.tl2)-AICc(FishStability.Base)),2)
    if(TrophicGroup != 'Browser'){  
      lm.lmax1<-glmmTMB(CVbiomass ~ Mlmax + Habitat + Source + (1|Site), data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mlmax),by=c('Habitat','Site')),REML=FALSE)
      lm.lmax2<-glmmTMB(CVbiomass ~ Mlmax * Habitat + Source + (1|Site), data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mlmax),by=c('Habitat','Site')),REML=FALSE)
    }else{
      lm.lmax1<-glm(CVbiomass ~ Mlmax + Habitat + Source, data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mlmax),by=c('Habitat','Site')))
      lm.lmax2<-glm(CVbiomass ~ Mlmax * Habitat + Source, data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mlmax),by=c('Habitat','Site')))
      }
    AICc(FishStability.Base);AICc(lm.lmax1);AICc(lm.lmax2) #dAIC
    dAIC_Mlmax <- round(min(AICc(lm.lmax1)-AICc(FishStability.Base),AICc(lm.lmax2)-AICc(FishStability.Base)),2)
  
    }else{
    
    lm.rich1<-glmmTMB(CVbiomass ~ Mchao + Habitat + (1|Site), data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mchao),by=c('Habitat','Site')),REML=FALSE)
    lm.rich2<-glmmTMB(CVbiomass ~ Mchao*Habitat + (1|Site), data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mchao),by=c('Habitat','Site')),REML=FALSE)
    AICc(FishStability.Base);AICc(lm.rich1);AICc(lm.rich2)
    dAIC_Mrich <- round(min(AICc(lm.rich2)-AICc(FishStability.Base),AICc(lm.rich1)-AICc(FishStability.Base)),2)
    
    lm.div1<-glmmTMB(CVbiomass ~ Mdiv + Habitat + (1|Site), data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mdiv),by=c('Habitat','Site')),REML=FALSE)
    lm.div2<-glmmTMB(CVbiomass ~ Mdiv * Habitat + (1|Site), data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mdiv),by=c('Habitat','Site')),REML=FALSE)
    AICc(FishStability.Base);AICc(lm.div1);AICc(lm.div2)
    dAIC_Mdiv <- round(min(AICc(lm.div1)-AICc(FishStability.Base),AICc(lm.div2)-AICc(FishStability.Base)),2)
    
    lm.tl1<-glmmTMB(CVbiomass ~ Mtl + Habitat+ (1|Site), data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mtl),by=c('Habitat','Site')),REML=FALSE)
    lm.tl2<-glmmTMB(CVbiomass ~ Mtl*Habitat + (1|Site), data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mtl),by=c('Habitat','Site')),REML=FALSE)
    AICc(FishStability.Base);AICc(lm.tl1);AICc(lm.tl2) #dAIC
    dAIC_Mtl <- round(min(AICc(lm.tl1)-AICc(FishStability.Base),AICc(lm.tl2)-AICc(FishStability.Base)),2)
    if(TrophicGroup != 'Browser'){  
      lm.lmax1<-glmmTMB(CVbiomass ~ Mlmax + Habitat + (1|Site), data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mlmax),by=c('Habitat','Site')),REML=FALSE)
      lm.lmax2<-glmmTMB(CVbiomass ~ Mlmax * Habitat + (1|Site), data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mlmax),by=c('Habitat','Site')),REML=FALSE)
    }else{
      lm.lmax1<-glmmTMB(CVbiomass ~ Mlmax + Habitat, data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mlmax),by=c('Habitat','Site')),REML=FALSE)
      lm.lmax2<-glmmTMB(CVbiomass ~ Mlmax * Habitat, data = select(CVtab,Habitat,Site,CVbiomass,Source) %>% left_join(select(Mtab,Habitat,Site,Mlmax),by=c('Habitat','Site')),REML=FALSE)
    }
    AICc(FishStability.Base);AICc(lm.lmax1);AICc(lm.lmax2) #dAIC
    dAIC_Mlmax <- round(min(AICc(lm.lmax1)-AICc(FishStability.Base),AICc(lm.lmax2)-AICc(FishStability.Base)),2)
  }
  
  #ModelSelection for fish diversity benthic stability
  options(na.action = na.fail)
  Merged <- CVtab  %>% left_join(Mtab,by=c('Habitat','Site','Source'))
  
  #CORAL MODELS
  lmBase_Source_Coral<-glmmTMB(SDCoral ~ Source + Habitat + (1|Site),data = Merged)
  lmBase_NoSource_Coral<-glmmTMB(SDCoral ~ Habitat + (1|Site),data = Merged)
  AICc(lmBase_NoSource_Coral);AICc(lmBase_Source_Coral)
  
  if(AICc(lmBase_NoSource_Coral)>AICc(lmBase_Source_Coral)){
    CoralStability.Source <- TRUE
    CoralStability.Base <- lmBase_Source_Coral
  }else{
    CoralStability.Source <- FALSE
    CoralStability.Base <- lmBase_NoSource_Coral
  }
  
  if(CoralStability.Source){
    
    lmBio_Coral1<- glmmTMB(SDCoral ~ Mbiomass+Habitat + Source (1|Site) ,data = Merged)
    lmBio_Coral2<- glmmTMB(SDCoral ~ Habitat*Mbiomass+ Source + (1|Site),data = Merged)
    AICc(CoralStability.Base);AICc(lmBio_Coral1);AICc(lmBio_Coral2)
    dAIC_Bio_Coral<-round(min(AICc(lmBio_Coral1)-AICc(CoralStability.Base),AICc(lmBio_Coral2)-AICc(CoralStability.Base)),2)
    
    lmRich_Coral1<- glmmTMB(SDCoral ~ Mchao + Habitat+ Source + (1|Site),data = Merged)
    lmRich_Coral2<- glmmTMB(SDCoral ~ Mchao*Habitat+ Source + (1|Site),data = Merged)
    AICc(CoralStability.Base);AICc(lmRich_Coral1);AICc(lmRich_Coral2)
    dAIC_Rich_Coral<-round(min(AICc(lmRich_Coral1)-AICc(CoralStability.Base),AICc(lmRich_Coral2)-AICc(CoralStability.Base)),2)
    
    
    lmDiv_Coral1<- glmmTMB(SDCoral ~ Source + Mdiv + Habitat+ Source + (1|Site),data = Merged)
    lmDiv_Coral2<- glmmTMB(SDCoral ~ Source + Mdiv*Habitat+ Source + (1|Site),data = Merged)
    AICc(CoralStability.Base);AICc(lmDiv_Coral1);AICc(lmDiv_Coral2)
    dAIC_Div_Coral<-round(min(AICc(lmDiv_Coral1)-AICc(CoralStability.Base),AICc(lmDiv_Coral2)-AICc(CoralStability.Base)),2)
      
    lmTL_Coral1<- glmmTMB(SDCoral ~ Source + Mtl + Habitat+ Source + (1|Site),data = Merged)
    lmTL_Coral2<- glmmTMB(SDCoral ~ Source + Mtl*Habitat+ Source + (1|Site),data = Merged)
    AICc(CoralStability.Base);AICc(lmTL_Coral1);AICc(lmTL_Coral2)
    dAIC_TL_Coral<-round(min(AICc(lmTL_Coral1)-AICc(CoralStability.Base),AICc(lmTL_Coral2)-AICc(CoralStability.Base)),2)
    
    lmLmax_Coral1<- glmmTMB(SDCoral ~  Source + Mlmax + Habitat+ Source + (1|Site),data = Merged)
    lmLmax_Coral2<- glmmTMB(SDCoral ~  Source + Mlmax * Habitat+ Source + (1|Site),data = Merged)
    AICc(CoralStability.Base);AICc(lmLmax_Coral1);AICc(lmLmax_Coral2)
    dAIC_Lmax_Coral<-round(min(AICc(lmLmax_Coral1)-AICc(CoralStability.Base),AICc(lmLmax_Coral2)-AICc(CoralStability.Base)),2)
  
  }else{
    lmBio_Coral1<- glmmTMB(SDCoral ~ Mbiomass+Habitat + (1|Site) ,data = Merged)
    lmBio_Coral2<- glmmTMB(SDCoral ~ Habitat*Mbiomass + (1|Site),data = Merged)
    AICc(CoralStability.Base);AICc(lmBio_Coral1);AICc(lmBio_Coral2)
    dAIC_Bio_Coral<-round(min(AICc(lmBio_Coral1)-AICc(CoralStability.Base),AICc(lmBio_Coral2)-AICc(CoralStability.Base)),2)
    
    lmRich_Coral1<- glmmTMB(SDCoral ~ Mchao + Habitat + (1|Site),data = Merged)
    lmRich_Coral2<- glmmTMB(SDCoral ~ Mchao*Habitat + (1|Site),data = Merged)
    AICc(CoralStability.Base);AICc(lmRich_Coral1);AICc(lmRich_Coral2)
    dAIC_Rich_Coral<-round(min(AICc(lmRich_Coral1)-AICc(CoralStability.Base),AICc(lmRich_Coral2)-AICc(CoralStability.Base)),2)
    
    
    lmDiv_Coral1<- glmmTMB(SDCoral ~ Mdiv + Habitat + (1|Site),data = Merged)
    lmDiv_Coral2<- glmmTMB(SDCoral ~ Mdiv*Habitat + (1|Site),data = Merged)
    AICc(CoralStability.Base);AICc(lmDiv_Coral1);AICc(lmDiv_Coral2)
    dAIC_Div_Coral<-round(min(AICc(lmDiv_Coral1)-AICc(CoralStability.Base),AICc(lmDiv_Coral2)-AICc(CoralStability.Base)),2)
    
    lmTL_Coral1<- glmmTMB(SDCoral ~ Mtl + Habitat + (1|Site),data = Merged)
    lmTL_Coral2<- glmmTMB(SDCoral ~ Mtl*Habitat + (1|Site),data = Merged)
    AICc(CoralStability.Base);AICc(lmTL_Coral1);AICc(lmTL_Coral2)
    dAIC_TL_Coral<-round(min(AICc(lmTL_Coral1)-AICc(CoralStability.Base),AICc(lmTL_Coral2)-AICc(CoralStability.Base)),2)
    
    lmLmax_Coral1<- glmmTMB(SDCoral ~  Mlmax + Habitat + (1|Site),data = Merged)
    lmLmax_Coral2<- glmmTMB(SDCoral ~  Mlmax * Habitat + (1|Site),data = Merged)
    AICc(CoralStability.Base);AICc(lmLmax_Coral1);AICc(lmLmax_Coral2)
    dAIC_Lmax_Coral<-round(min(AICc(lmLmax_Coral1)-AICc(CoralStability.Base),AICc(lmLmax_Coral2)-AICc(CoralStability.Base)),2)
  }
  
    
  
  #MACROALGAE MODELS
  lmBase_Source_Macro<-glmmTMB(SDMacro ~  Source + Habitat + (1|Site),data = Merged)
  lmBase_NoSource_Macro<-glmmTMB(SDMacro ~  Habitat + (1|Site),data = Merged)
  AICc(lmBase_NoSource_Macro);AICc(lmBase_Source_Macro)
  if(AICc(lmBase_NoSource_Macro)>AICc(lmBase_Source_Macro)){
    MacroStability.Source <- TRUE
    MacroStability.Base <- lmBase_Source_Macro
  }else{
    MacroStability.Source <- FALSE
    MacroStability.Base <- lmBase_NoSource_Macro
  }
  
  if(MacroStability.Source){
    lmBio_Macro1<- glmmTMB(SDMacro ~ Mbiomass + Habitat + Source + (1|Site),data = Merged)
    lmBio_Macro2 <-glmmTMB(SDMacro ~ Mbiomass*Habitat + Source + (1|Site),data = Merged)
    AICc(MacroStability.Base);AICc(lmBio_Macro1);AICc(lmBio_Macro2)
    dAIC_Bio_Macro<-round(min(AICc(lmBio_Macro1)-AICc(MacroStability.Base),AICc(lmBio_Macro2)-AICc(MacroStability.Base)),2)
    
    lmRich_Macro1 <- glmmTMB(SDMacro ~ Mchao + Habitat + Source + (1|Site),data = Merged)
    lmRich_Macro2 <- glmmTMB(SDMacro ~ Mchao*Habitat + Source + (1|Site),data = Merged)
    AICc(MacroStability.Base);AICc(lmRich_Macro1);AICc(lmRich_Macro2)
    dAIC_Rich_Macro<-round(min(AICc(lmRich_Macro1)-AICc(MacroStability.Base),AICc(lmRich_Macro2)-AICc(MacroStability.Base)),2)
    
    lmDiv_Macro1<- glmmTMB(SDMacro ~ Mdiv + Habitat + Source + (1|Site),data = Merged)
    lmDiv_Macro2<-glmmTMB(SDMacro ~ Mdiv * Habitat + Source + (1|Site),data = Merged)
    AICc(MacroStability.Base);AICc(lmDiv_Macro1);AICc(lmDiv_Macro2)
    dAIC_Div_Macro<-round(min(AICc(lmDiv_Macro1)-AICc(MacroStability.Base),AICc(lmDiv_Macro2)-AICc(MacroStability.Base)),2)
    
    lmTL_Macro1<- glmmTMB(SDMacro ~ Mtl + Habitat + Source + (1|Site),data = Merged)
    lmTL_Macro2<-glmmTMB(SDMacro ~ Mtl*Habitat + Source  + (1|Site),data = Merged)
    AICc(MacroStability.Base);AICc(lmTL_Macro1);AICc(lmTL_Macro2)
    dAIC_TL_Macro<-round(min(AICc(lmTL_Macro1)-AICc(MacroStability.Base),AICc(lmTL_Macro2)-AICc(MacroStability.Base)),2)
    
    lmLmax_Macro1<- glmmTMB(SDMacro ~ Mlmax + Habitat + Source + (1|Site),data = Merged)
    lmLmax_Macro2 <-glmmTMB(SDMacro ~ Mlmax*Habitat + Source + (1|Site),data = Merged)
    AICc(MacroStability.Base);AICc(lmLmax_Macro1);AICc(lmLmax_Macro2)
    dAIC_Lmax_Macro<-round(min(AICc(lmLmax_Macro1)-AICc(MacroStability.Base),AICc(lmLmax_Macro2)-AICc(MacroStability.Base)),2)
  
  }else{
    
    lmBio_Macro1<- glmmTMB(SDMacro ~ Mbiomass + Habitat + (1|Site),data = Merged)
    lmBio_Macro2 <-glmmTMB(SDMacro ~ Mbiomass*Habitat + (1|Site),data = Merged)
    AICc(MacroStability.Base);AICc(lmBio_Macro1);AICc(lmBio_Macro2)
    dAIC_Bio_Macro<-round(min(AICc(lmBio_Macro1)-AICc(MacroStability.Base),AICc(lmBio_Macro2)-AICc(MacroStability.Base)),2)
    
    lmRich_Macro1 <- glmmTMB(SDMacro ~ Mchao + Habitat + (1|Site),data = Merged)
    lmRich_Macro2 <- glmmTMB(SDMacro ~ Mchao*Habitat + (1|Site),data = Merged)
    AICc(MacroStability.Base);AICc(lmRich_Macro1);AICc(lmRich_Macro2)
    dAIC_Rich_Macro<-round(min(AICc(lmRich_Macro1)-AICc(MacroStability.Base),AICc(lmRich_Macro2)-AICc(MacroStability.Base)),2)
    
    lmDiv_Macro1<- glmmTMB(SDMacro ~ Mdiv + Habitat + (1|Site),data = Merged)
    lmDiv_Macro2<-glmmTMB(SDMacro ~ Mdiv * Habitat + (1|Site),data = Merged)
    AICc(MacroStability.Base);AICc(lmDiv_Macro1);AICc(lmDiv_Macro2)
    dAIC_Div_Macro<-round(min(AICc(lmDiv_Macro1)-AICc(MacroStability.Base),AICc(lmDiv_Macro2)-AICc(MacroStability.Base)),2)
    
    lmTL_Macro1<- glmmTMB(SDMacro ~ Mtl + Habitat + (1|Site),data = Merged)
    lmTL_Macro2<-glmmTMB(SDMacro ~ Mtl*Habitat  + (1|Site),data = Merged)
    AICc(MacroStability.Base);AICc(lmTL_Macro1);AICc(lmTL_Macro2)
    dAIC_TL_Macro<-round(min(AICc(lmTL_Macro1)-AICc(MacroStability.Base),AICc(lmTL_Macro2)-AICc(MacroStability.Base)),2)
    
    lmLmax_Macro1<- glmmTMB(SDMacro ~ Mlmax + Habitat + (1|Site),data = Merged)
    lmLmax_Macro2 <-glmmTMB(SDMacro ~ Mlmax*Habitat + (1|Site),data = Merged)
    AICc(MacroStability.Base);AICc(lmLmax_Macro1);AICc(lmLmax_Macro2)
    dAIC_Lmax_Macro<-round(min(AICc(lmLmax_Macro1)-AICc(MacroStability.Base),AICc(lmLmax_Macro2)-AICc(MacroStability.Base)),2)
  }
  
  
    ##### DOES Stability predict Stability
    Coral_StabilityStability_Base <- glmmTMB(SDCoral~Habitat + (1|Site),data=Merged)
    Coral_StabilityStability_Biomass1<-glmmTMB(SDCoral~ CVbiomass + Habitat + (1|Site),data=Merged)
    Coral_StabilityStability_Biomass2<-glmmTMB(SDCoral~ CVbiomass*Habitat + (1|Site),data=Merged)
    Macro_StabilityStability_Base <- glmmTMB(SDMacro~Habitat + (1|Site),data=Merged)
    Macro_StabilityStability_Biomass1<-glmmTMB(SDMacro~ CVbiomass + Habitat + (1|Site),data=Merged)
    Macro_StabilityStability_Biomass2<-glmmTMB(SDMacro~ CVbiomass*Habitat + (1|Site),data=Merged)
  
}#End Loop Over Trophic Groups

if(TrophicGroup == 'FishComm'){
  quartz(width=7,height=4.5)
  if(TRUE){
    par(oma=c(2.5,2.5,0.5,2))
    par(mar=c(0,0,0,0))
    par(lheight=0.85)
    par(mgp=c(0,0.5,0))
    
    tckl<-0
    acx<-7/12
    lcx<-9/12
    Transparency<- '55'
    bottom <- 0.1 # Bottom of main plotting region
    
    xlim_bio<-c(0,200)
    xlim_rich<-c(0,80)
    xlim_div <-c(1,3.2)
    xlim_tl <-c(2.8,3.6)
    xlim_lmax <-c(100,200)
    
    xaxis_bio<-c(0,50,100,150,200)
    xaxis_rich<-c('0','20','40','60')
    xaxis_div <-c(1,1.5,2,2.5,3)
    xaxis_tl <-c(2.9,3.1,3.3,3.5)
    xaxis_lmax <-c(100,125,150,175,200)
    
    ylim_fish<-c(-0.05,2.2)
    #range(exp(Merged$SDMacro))
    ylim_coral<-c(-0.1,2.1)
    ylim_macro<-c(0,1.6)
      
    yaxis_fish<-c(0,1,2)
    yaxis_coral<-c(0,1,2)
    yaxis_macro<-c(0,0.5,1,1.5)
    
    pos1<- -0.7
    pos2 <- -7
    pos3 <- -13.2
    pos4 <- -19.3
    pos5 <- -25.5
    dAICline<--0.9
    dAICexp<-7/12
    dAICfrac<-0.70
    
    npreds<-300
    
    BTY <-'o'
    
  }#Setup
  
  if(TRUE){
    #BluePal<-brewer.pal(9,'Blues')[c(3,6,9)]#c('lightblue','dodgerblue','darkblue')
    GreenPal<-brewer.pal(9,'Greens')[c(3,6,9)]
    
    HabPts<-apply(Mtab[,'Habitat'],1,FUN=function(x){ifelse(x=='BA',22,ifelse(x=='FR',21,24))})
    
    par(fig=c(0.2,0.4,0.666,1))
    
    histMbio<-hist((Mtab$Mbiomass),breaks=10,plot=F)
    GrPal<-colorRampPalette(brewer.pal(9,'Greens'))
    GRS<-GrPal((length(histMbio$mids)))
    MbioCol<-GRS[apply(matrix(Mtab$Mbiomass,ncol=1),1,FUN=function(x){which.min(Mod(x-histMbio$mids))})]
    MbioOutline<-'black'#MbioCol
    #MbioOutline[Mtab$Source=='LTER']<-'black'
    
    plot(exp(CVtab$CVbiomass)~Mtab$Mchao,xlim=xlim_rich,ylim=ylim_fish,ann=F,axes=F,pch=HabPts,col=MbioOutline,bg=MbioCol,xpd=T)
    axis(1,at=xaxis_rich,labels=NA,cex.axis=acx,tck=tckl)
    axis(2,at=yaxis_fish,cex.axis=acx,tck=-0.05)
    mtext(bquote(Delta~'AICc' == .(round(dAIC_Mrich,2))),1,line=dAICline,at=min(xlim_rich)+dAICfrac*diff(xlim_rich),cex=dAICexp)
    #segments(x0=10,x1=55,y0=0.971526-0.006386*10,y1=0.971526-0.006386*55,lwd=1.5,lty=2)
    if(dAIC_Mrich < -2){
      plotMod <- list(lm.rich1,lm.rich2)[[which.min(c(AICc(lm.rich1),AICc(lm.rich2)))]]
      xvar<-Mtab$Mchao
      pred_x<-seq(min(xvar),max(xvar),length=npreds)
      pred_x_FO<-data.frame(Mchao=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('FO',npreds),Site=NA)
      pred_x_FR<-data.frame(Mchao=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('FR',npreds),Site=NA)
      pred_x_BA<-data.frame(Mchao=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('BA',npreds),Site=NA)
      #list(Mdiv=rep(seq(1,3,length=100),1),Source=as.factor(c('LTER',rep('CRIOBE',100)))[-1],Habitat=as.factor(c('BA','FR',rep('FO',100)))[-c(1,2)])
      pred_y_FO<-exp(predict(plotMod,newdata=pred_x_FO,re.form=NA))
      pred_y_FR<-exp(predict(plotMod,newdata=pred_x_FR,re.form=NA))
      pred_y_BA<-exp(predict(plotMod,newdata=pred_x_BA,re.form=NA))
      
      points(pred_x,pred_y_FO,type='l',lwd=2,lty=1)
      points(pred_x,pred_y_FR,type='l',lwd=2,lty=2)
      points(pred_x,pred_y_BA,type='l',lwd=2,lty=3)
    }
    
    box(bty=BTY)
    #mtext('Fish richness',1,line=1.25,cex=lcx)
    mtext('CV fish biomass',2,line=1.25,cex=lcx)
    
    # par(fig=c(0.31,0.35,0.5,1),new=T)
    # par(mar=c(2,0,1,0))
    # pX<-histMbio$y/max(histMbio$y)
    # pY<-histMbio$x/max(histMbio$x)
    # plot(1,xlim=c(0,1),ylim=c(0,1),type='n',axes=F,ann=F,xaxs='i',yaxs='i')
    # for(i in 2:length(histMbio$x)){
    # 	polygon(x=c(0,0,pX[i-1],pX[i]),y=c(1-pY[i],1-pY[i-1],1-pY[i-1],1-pY[i]),col=GRS[i-1],border=GRS[i-1])
    # }
    
    par(fig=c(0.4,0.6,0.666,1),new=T)
    plot(exp(CVtab$CVbiomass)~Mtab$Mdiv,xlim=xlim_div,ylim=ylim_fish,ann=F,axes=F,pch=HabPts,col=MbioOutline,bg=MbioCol,xpd=T)
    #axis(1,at=c(1.5,2.0,2.5,3.0),labels=c('1.5','',2.5,''),cex.axis=acx,tck=tckl)
    axis(1,at=xaxis_div,labels=NA,cex.axis=acx,tck=tckl)
    axis(2,at=yaxis_fish,labels=NA,cex.axis=acx,tck=tckl)
    mtext(bquote(Delta~'AICc' == .(round(dAIC_Mdiv,2))),1,line=dAICline,at=min(xlim_div)+dAICfrac*diff(xlim_div),cex=dAICexp)
    #
    box(bty=BTY)
    summary(lm.div1)
    #abline(1.41,-0.285)
    #pred_x_FO<-data.frame(expand.grid(Mdiv=seq(1,3,length=100),Source=unique(Merged$Source),Habitat=unique(Merged$Habitat),Site=unique(Merged$Site)))
    if(dAIC_Mdiv < -2){
      plotMod <- list(lm.div1,lm.div2)[[which.min(c(AICc(lm.div1),AICc(lm.div2)))]]
      xvar<-Mtab$Mdiv
      pred_x<-seq(min(xvar),max(xvar),length=npreds)
      pred_x_FO<-data.frame(Mdiv=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('FO',npreds),Site=NA)
      pred_x_FR<-data.frame(Mdiv=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('FR',npreds),Site=NA)
      pred_x_BA<-data.frame(Mdiv=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('BA',npreds),Site=NA)
      #list(Mdiv=rep(seq(1,3,length=100),1),Source=as.factor(c('LTER',rep('CRIOBE',100)))[-1],Habitat=as.factor(c('BA','FR',rep('FO',100)))[-c(1,2)])
      pred_y_FO<-exp(predict(plotMod,newdata=pred_x_FO,re.form=NA))
      pred_y_FR<-exp(predict(plotMod,newdata=pred_x_FR,re.form=NA))
      pred_y_BA<-exp(predict(plotMod,newdata=pred_x_BA,re.form=NA))
      
      points(pred_x,pred_y_FO,type='l',lwd=2,lty=1)
      points(pred_x,pred_y_FR,type='l',lwd=2,lty=2)
      points(pred_x,pred_y_BA,type='l',lwd=2,lty=3)
    }
    
    #segments(x0=1.5,x1=3.2,y0=exp(1.4007-0.2821*1.5,y1=1.4007-0.2821*3.0,lwd=1.5,lty=1)
    #mtext('Fish diversity',1,line=1.25,cex=lcx)
    #mtext('CV fish biomass',2,line=1.25,cex=lcx)
    
    
    par(fig=c(0.6,0.8,0.666,1),new=T)
    plot(exp(CVtab$CVbiomass)~Mtab$Mtl,xlim=xlim_tl,ylim=ylim_fish,ann=F,axes=F,pch=HabPts,col=MbioOutline,bg=MbioCol,xpd=T)
    axis(1,at=xaxis_tl,labels=NA,cex.axis=acx,tck=tckl)
    axis(2,at=yaxis_fish,labels=NA,cex.axis=acx,tck=tckl)
    mtext(bquote(Delta~'AICc' == .(round(dAIC_Mtl,2))),1,line=dAICline,at=min(xlim_tl)+dAICfrac*diff(xlim_tl),cex=dAICexp)
    box(bty=BTY)
    if(dAIC_Mtl < -2){
      plotMod <- list(lm.tl1,lm.tl2)[[which.min(c(AICc(lm.tl1),AICc(lm.tl2)))]]
      xvar<-Mtab$Mtl
      pred_x<-seq(min(xvar),max(xvar),length=npreds)
      pred_x_FO<-data.frame(Mtl=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('FO',npreds),Site=NA)
      pred_x_FR<-data.frame(Mtl=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('FR',npreds),Site=NA)
      pred_x_BA<-data.frame(Mtl=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('BA',npreds),Site=NA)
      #list(Mdiv=rep(seq(1,3,length=100),1),Source=as.factor(c('LTER',rep('CRIOBE',100)))[-1],Habitat=as.factor(c('BA','FR',rep('FO',100)))[-c(1,2)])
      pred_y_FO<-exp(predict(plotMod,newdata=pred_x_FO,re.form=NA))
      pred_y_FR<-exp(predict(plotMod,newdata=pred_x_FR,re.form=NA))
      pred_y_BA<-exp(predict(plotMod,newdata=pred_x_BA,re.form=NA))
      
      points(pred_x,pred_y_FO,type='l',lwd=2,lty=1)
      points(pred_x,pred_y_FR,type='l',lwd=2,lty=2)
      points(pred_x,pred_y_BA,type='l',lwd=2,lty=3)
    }
    
    #mtext('Mean trophic level',1,line=1.25,cex=lcx)
    #mtext('CV fish biomass',2,line=1.25,cex=lcx)
    
    par(fig=c(0.8,1.0,0.666,1),new=T)
    plot(exp(CVtab$CVbiomass)~Mtab$Mlmax,xlim=xlim_lmax,ylim=ylim_fish,ann=F,axes=F,pch=HabPts,col=MbioOutline,bg=MbioCol)
    axis(1,at=xaxis_lmax,labels=NA,cex.axis=acx,tck=tckl)
    axis(2,at=yaxis_fish,labels=NA,cex.axis=acx,tck=tckl)
    mtext(bquote(Delta~'AICc' == .(round(dAIC_Mlmax,2))),1,line=dAICline,at=min(xlim_lmax)+dAICfrac*diff(xlim_lmax),cex=dAICexp)
    box(bty=BTY)
     
    if(dAIC_Mlmax < -2){
      plotMod <- list(lm.lmax1,lm.lmax2)[[which.min(c(AICc(lm.lmax1),AICc(lm.lmax2)))]]
      xvar<-Mtab$Mlmax
      pred_x<-seq(min(xvar),max(xvar),length=npreds)
      pred_x_FO<-data.frame(Mlmax=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('FO',npreds),Site=NA)
      pred_x_FR<-data.frame(Mlmax=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('FR',npreds),Site=NA)
      pred_x_BA<-data.frame(Mlmax=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('BA',npreds),Site=NA)
      #list(Mdiv=rep(seq(1,3,length=100),1),Source=as.factor(c('LTER',rep('CRIOBE',100)))[-1],Habitat=as.factor(c('BA','FR',rep('FO',100)))[-c(1,2)])
      pred_y_FO<-exp(predict(plotMod,newdata=pred_x_FO,re.form=NA))
      pred_y_FR<-exp(predict(plotMod,newdata=pred_x_FR,re.form=NA))
      pred_y_BA<-exp(predict(plotMod,newdata=pred_x_BA,re.form=NA))
      
      points(pred_x,pred_y_FO,type='l',lwd=2,lty=1)
      points(pred_x,pred_y_FR,type='l',lwd=2,lty=2)
      points(pred_x,pred_y_BA,type='l',lwd=2,lty=3)
    }
    
    
    # pos1<- 0
    # pos2 <- -6.75
    # pos3 <- -13
    # pos4 <- -19.6
    # pos5 <- -26.25
    
    #mtext('A',2,outer=T,cex=lcx,las=1,line=pos2,at=0.98)
    #mtext('B',2,outer=T,cex=lcx,las=1,line=pos3,at=0.98)
    #mtext('C',2,outer=T,cex=lcx,las=1,line=pos4,at=0.98)
    #mtext('D',2,outer=T,cex=lcx,las=1,line=pos5,at=0.98)
    
    # mtext('J',2,outer=T,cex=lcx,las=1,line=0,at=0.34)
    # mtext('K',2,outer=T,cex=lcx,las=1,line=-6,at=0.34)
    # mtext('L',2,outer=T,cex=lcx,las=1,line=-12.5,at=0.34)
    # mtext('M',2,outer=T,cex=lcx,las=1,line=-18.75,at=0.34)
    
    
    #dAICexp<-7/12
    #mtext(bquote(Delta~'AIC' == .(round(dAIC_Mrich,2))),3,outer=T,line=-1,at=0.28,cex=dAICexp)
    #mtext(bquote(Delta~'AIC' == .(round(dAIC_Mdiv,2))),3,outer=T,line=-1,at=0.48,cex=dAICexp)
    #mtext(bquote(Delta~'AIC' == .(round(dAIC_Mtl,2))),3,outer=T,line=-1,at=0.68,cex=dAICexp)
    #mtext(bquote(Delta~'AIC' == .(round(dAIC_lmax,2))),3,outer=T,line=-1,at=0.88,cex=dAICexp)
    
    
    
    
    
  }#Plot FISH_DIVERISYT STABILTIY
  
  if(TRUE){
    
    par(fig=c(0,0.2,0.333,0.666),new=T)
    
    histMcoral<-hist(ilogit(Mtab$MCoral),breaks=10,plot=F)
    BuPal<-colorRampPalette(brewer.pal(9,'Blues'))
    BUS<-BuPal(length(histMcoral$mids))
    MCoralCol<-BUS[apply(matrix(ilogit(Mtab$MCoral),ncol=1),1,FUN=function(x){which.min(Mod(x-histMcoral$mids))})]
    MCoralOutline<-MCoralCol
    HabPtsCoral<-apply(Mtab[,'Habitat'],1,FUN=function(x){ifelse(x=='BA',22,ifelse(x=='FR',21,24))})
    MCoralOutline <- 'black'
    #MCoralOutline[Mtab$Source=='LTER']<-'black'
    
    plot(exp(SDCoral)~exp(Mbiomass),data=Merged,xlim=xlim_bio,ylim=ylim_coral,ann=F,axes=F,pch=HabPtsCoral,col=MCoralOutline,bg=MCoralCol,xpd=T)
    axis(1,at=xaxis_bio,labels=NA,cex.axis=acx,tck=tckl)
    axis(2,at=yaxis_coral,cex.axis=acx,tck=-0.05)
    mtext(bquote(Delta~'AICc' == .(round(dAIC_Bio_Coral,2))),1,line=dAICline,at=min(xlim_bio)+dAICfrac*diff(xlim_bio),cex=dAICexp)
    if(dAIC_Bio_Coral < -2){
      plotMod <- list(lmBio_Coral1,lmBio_Coral2)[[which.min(c(AICc(lmBio_Coral1),AICc(lmBio_Coral2)))]]
      xvar<-Mtab$Mbiomass
      pred_x<-seq(min(xvar),max(xvar),length=npreds)
      pred_x_FO<-data.frame(Mbiomass=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('FO',npreds),Site=NA)
      pred_x_FR<-data.frame(Mbiomass=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('FR',npreds),Site=NA)
      pred_x_BA<-data.frame(Mbiomass=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('BA',npreds),Site=NA)
      #list(Mdiv=rep(seq(1,3,length=100),1),Source=as.factor(c('LTER',rep('CRIOBE',100)))[-1],Habitat=as.factor(c('BA','FR',rep('FO',100)))[-c(1,2)])
      pred_y_FO<-exp(predict(plotMod,newdata=pred_x_FO,re.form=NA))
      pred_y_FR<-exp(predict(plotMod,newdata=pred_x_FR,re.form=NA))
      pred_y_BA<-exp(predict(plotMod,newdata=pred_x_BA,re.form=NA))
      
      points(exp(pred_x),pred_y_FO,type='l',lwd=2,lty=1)
      points(exp(pred_x),pred_y_FR,type='l',lwd=2,lty=2)
      points(exp(pred_x),pred_y_BA,type='l',lwd=2,lty=3)
    }
    
    
    box(bty=BTY)
    #mtext('Fish biomass',1,line=1.25,cex=lcx)
    mtext('SD coral cover',2,line=1.25,cex=lcx)
    
    
    par(fig=c(0.2,0.4,0.333,0.666),new=T)
    plot(exp(SDCoral)~Mchao,data=Merged,xlim=xlim_rich,ylim=ylim_coral,ann=F,axes=F,pch=HabPtsCoral,col=MCoralOutline,bg=MCoralCol,xpd=T)
    #summary(lm2)
    axis(1,at=xaxis_rich,labels=NA,cex.axis=acx,tck=tckl)
    axis(2,at=yaxis_coral,labels=NA,cex.axis=acx,tck=tckl)
    mtext(bquote(Delta~'AICc' == .(round(dAIC_Rich_Coral,2))),1,line=dAICline,at=min(xlim_rich)+dAICfrac*diff(xlim_rich),cex=dAICexp)
     if(dAIC_Rich_Coral < -2){
      plotMod <- list(lmRich_Coral1,lmRich_Coral2)[[which.min(c(AICc(lmRich_Coral1),AICc(lmRich_Coral2)))]]
      xvar<-Mtab$Mchao
      pred_x<-seq(min(xvar),max(xvar),length=npreds)
      pred_x_FO<-data.frame(Mchao=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('FO',npreds),Site=NA)
      pred_x_FR<-data.frame(Mchao=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('FR',npreds),Site=NA)
      pred_x_BA<-data.frame(Mchao=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('BA',npreds),Site=NA)
      #list(Mdiv=rep(seq(1,3,length=100),1),Source=as.factor(c('LTER',rep('CRIOBE',100)))[-1],Habitat=as.factor(c('BA','FR',rep('FO',100)))[-c(1,2)])
      pred_y_FO<-exp(predict(plotMod,newdata=pred_x_FO,re.form=NA))
      pred_y_FR<-exp(predict(plotMod,newdata=pred_x_FR,re.form=NA))
      pred_y_BA<-exp(predict(plotMod,newdata=pred_x_BA,re.form=NA))
      
      points(pred_x,pred_y_FO,type='l',lwd=2,lty=1)
      points(pred_x,pred_y_FR,type='l',lwd=2,lty=2)
      points(pred_x,pred_y_BA,type='l',lwd=2,lty=3)
    }
    
    box(bty=BTY)
    #mtext('Fish richness',1,line=1.25,cex=lcx)
    #mtext('SD coral',2,line=1.25,cex=lcx)
    
    par(fig=c(0.4,0.6,0.333,0.666),new=T)
    plot(exp(SDCoral)~Mdiv,data=Merged,xlim=xlim_div,ylim=ylim_coral,ann=F,axes=F,pch=HabPtsCoral,col=MCoralOutline,bg=MCoralCol,xpd=T)
    axis(1,at=xaxis_div,labels=NA,cex.axis=acx,tck=tckl)
    axis(2,at=yaxis_coral,labels=NA,cex.axis=acx,tck=tckl)
    mtext(bquote(Delta~'AICc' == .(round(dAIC_Div_Coral,2))),1,line=dAICline,at=min(xlim_div)+dAICfrac*diff(xlim_div),cex=dAICexp)
    box(bty=BTY)
     if(dAIC_Div_Coral < -2){
      plotMod <- list(lmDiv_Coral1,lmDiv_Coral2)[[which.min(c(AICc(lmDiv_Coral1),AICc(lmDiv_Coral2)))]]
      xvar<-Mtab$Mdiv
      pred_x<-seq(min(xvar),max(xvar),length=npreds)
      pred_x_FO<-data.frame(Mdiv=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('FO',npreds),Site=NA)
      pred_x_FR<-data.frame(Mdiv=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('FR',npreds),Site=NA)
      pred_x_BA<-data.frame(Mdiv=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('BA',npreds),Site=NA)
      #list(Mdiv=rep(seq(1,3,length=100),1),Source=as.factor(c('LTER',rep('CRIOBE',100)))[-1],Habitat=as.factor(c('BA','FR',rep('FO',100)))[-c(1,2)])
      pred_y_FO<-exp(predict(plotMod,newdata=pred_x_FO,re.form=NA))
      pred_y_FR<-exp(predict(plotMod,newdata=pred_x_FR,re.form=NA))
      pred_y_BA<-exp(predict(plotMod,newdata=pred_x_BA,re.form=NA))
      
      points(pred_x,pred_y_FO,type='l',lwd=2,lty=1)
      points(pred_x,pred_y_FR,type='l',lwd=2,lty=2)
      points(pred_x,pred_y_BA,type='l',lwd=2,lty=3)
    }
    
    #mtext('Fish diversity',1,line=1.25,cex=lcx)
    #mtext('SD coral',2,line=1.25,cex=lcx)
    
    par(fig=c(0.6,0.8,0.333,0.666),new=T)
    plot(exp(SDCoral)~Mtl,data=Merged,xlim=xlim_tl,ylim=ylim_coral,ann=F,axes=F,pch=HabPtsCoral,col=MCoralOutline,bg=MCoralCol,xpd=T)
    axis(1,at=xaxis_tl,labels=NA,cex.axis=acx,tck=tckl)
    axis(2,at=yaxis_coral,labels=NA,cex.axis=acx,tck=tckl)
    mtext(bquote(Delta~'AICc' == .(round(dAIC_TL_Coral,2))),1,line=dAICline,at=min(xlim_tl)+dAICfrac*diff(xlim_tl),cex=dAICexp)
    box(bty=BTY)
     if(dAIC_TL_Coral < -2){
      plotMod <- list(lmTL_Coral1,lmTL_Coral2)[[which.min(c(AICc(lmTL_Coral1),AICc(lmTL_Coral2)))]]
      xvar<-Mtab$Mtl
      pred_x<-seq(min(xvar),max(xvar),length=npreds)
      pred_x_FO<-data.frame(Mtl=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('FO',npreds),Site=NA)
      pred_x_FR<-data.frame(Mtl=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('FR',npreds),Site=NA)
      pred_x_BA<-data.frame(Mtl=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('BA',npreds),Site=NA)
      #list(Mdiv=rep(seq(1,3,length=100),1),Source=as.factor(c('LTER',rep('CRIOBE',100)))[-1],Habitat=as.factor(c('BA','FR',rep('FO',100)))[-c(1,2)])
      pred_y_FO<-exp(predict(plotMod,newdata=pred_x_FO,re.form=NA))
      pred_y_FR<-exp(predict(plotMod,newdata=pred_x_FR,re.form=NA))
      pred_y_BA<-exp(predict(plotMod,newdata=pred_x_BA,re.form=NA))
      
      points(pred_x,pred_y_FO,type='l',lwd=2,lty=1)
      points(pred_x,pred_y_FR,type='l',lwd=2,lty=2)
      points(pred_x,pred_y_BA,type='l',lwd=2,lty=3)
    }
    #mtext('Mean trophic level',1,line=1.25,cex=lcx)
    
    par(fig=c(0.8,1,0.333,0.666),new=T)
    plot(exp(SDCoral)~Mlmax,data=Merged,xlim=xlim_lmax,ylim=ylim_coral,ann=F,axes=F,pch=HabPtsCoral,col=MCoralOutline,bg=MCoralCol,xpd=T)
    axis(1,at=xaxis_lmax,labels=NA,cex.axis=acx,tck=tckl)
    axis(2,at=yaxis_coral,labels=NA,cex.axis=acx,tck=tckl)
    mtext(bquote(Delta~'AICc' == .(round(dAIC_Lmax_Coral,2))),1,line=dAICline,at=min(xlim_lmax)+dAICfrac*diff(xlim_lmax),cex=dAICexp)
    box(bty=BTY)
    if(dAIC_Lmax_Coral < -2){
      plotMod <- list(lmLmax_Coral1,lmLmax_Coral2)[[which.min(c(AICc(lmLmax_Coral1),AICc(lmLmax_Coral2)))]]
      xvar<-Mtab$Mlmax
      pred_x<-seq(min(xvar),max(xvar),length=npreds)
      pred_x_FO<-data.frame(Mlmax=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('FO',npreds),Site=NA)
      pred_x_FR<-data.frame(Mlmax=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('FR',npreds),Site=NA)
      pred_x_BA<-data.frame(Mlmax=pred_x,Source=rep('CRIOBE',npreds),Habitat=rep('BA',npreds),Site=NA)
      #list(Mdiv=rep(seq(1,3,length=100),1),Source=as.factor(c('LTER',rep('CRIOBE',100)))[-1],Habitat=as.factor(c('BA','FR',rep('FO',100)))[-c(1,2)])
      pred_y_FO<-exp(predict(plotMod,newdata=pred_x_FO,re.form=NA))
      pred_y_FR<-exp(predict(plotMod,newdata=pred_x_FR,re.form=NA))
      pred_y_BA<-exp(predict(plotMod,newdata=pred_x_BA,re.form=NA))
      
      points(pred_x,pred_y_FO,type='l',lwd=2,lty=1)
      points(pred_x,pred_y_FR,type='l',lwd=2,lty=2)
      points(pred_x,pred_y_BA,type='l',lwd=2,lty=3)
    }
    #mtext('Max length (mm)',1,line=1.25,cex=lcx)
    
    # mtext('E',2,outer=T,cex=lcx,las=1,line=pos1,at=0.64)
    # mtext('F',2,outer=T,cex=lcx,las=1,line=pos2,at=0.64)
    # mtext('G',2,outer=T,cex=lcx,las=1,line=pos3,at=0.64)
    # mtext('H',2,outer=T,cex=lcx,las=1,line=pos4,at=0.64)
    # mtext('I',2,outer=T,cex=lcx,las=1,line=pos5,at=0.64)
  } #Coral Models  
    
  if(TRUE){ 
    
    par(fig=c(0,0.2,0,0.333),new=T)
    
    histMmacro<-hist(ilogit(Mtab$MMacro),breaks=10,plot=F)
    RdPal<-colorRampPalette(brewer.pal(9,'Reds'))
    RDS<-RdPal(length(histMmacro$mids))
    MMacroCol<-RDS[apply(matrix(ilogit(Mtab$MMacro),ncol=1),1,FUN=function(x){which.min(Mod(x-histMmacro$mids))})]
    MMacroOutline<-MMacroCol
    HabPtsMacro<-apply(Mtab[,'Habitat'],1,FUN=function(x){ifelse(x=='BA',22,ifelse(x=='FR',21,24))})
    MMacroOutline<-'black'
    #MMacroOutline[Mtab$Source=='LTER']<-'black'
    
    plot(exp(SDMacro)~exp(Mbiomass),data=Merged,xlim=xlim_bio,ylim=ylim_macro,ann=F,axes=F,pch=HabPtsMacro,col=MMacroOutline,bg=MMacroCol,xpd=T)
    axis(1,at=xaxis_bio,cex.axis=acx,tck=-0.05)
    axis(2,at=yaxis_macro,cex.axis=acx,tck=-0.05)
    mtext(bquote(Delta~'AICc' == .(round(dAIC_Bio_Macro,2))),1,line=dAICline,at=min(xlim_bio)+dAICfrac*diff(xlim_bio),cex=dAICexp)
    box(bty=BTY)
    if(dAIC_Bio_Macro < -2){
      plotMod <- list(lmBio_Macro1,lmBio_Macro2)[[which.min(c(AICc(lmBio_Macro1),AICc(lmBio_Macro2)))]]
      xvar<-Mtab$Mbiomass
      pred_x<-seq(min(xvar),max(xvar),length=npreds)
      pred_x_FO<-data.frame(Mbiomass=pred_x,Habitat=rep('FO',npreds),Site=NA)
      pred_x_FR<-data.frame(Mbiomass=pred_x,Habitat=rep('FR',npreds),Site=NA)
      pred_x_BA<-data.frame(Mbiomass=pred_x,Habitat=rep('BA',npreds),Site=NA)
      #list(Mdiv=rep(seq(1,3,length=100),1),Source=as.factor(c('LTER',rep('CRIOBE',100)))[-1],Habitat=as.factor(c('BA','FR',rep('FO',100)))[-c(1,2)])
      pred_y_FO<-exp(predict(plotMod,newdata=pred_x_FO,re.form=NA))
      pred_y_FR<-exp(predict(plotMod,newdata=pred_x_FR,re.form=NA))
      pred_y_BA<-exp(predict(plotMod,newdata=pred_x_BA,re.form=NA))
      
      points(exp(pred_x),pred_y_FO,type='l',lwd=2,lty=1)
      points(exp(pred_x),pred_y_FR,type='l',lwd=2,lty=2)
      points(exp(pred_x),pred_y_BA,type='l',lwd=2,lty=3)
    }
    
    mtext(bquote('Fish biomass ' (g*'*'*m^-2)),1,line=1.35,cex=lcx)
    mtext('SD macroalgae cover',2,line=1.25,cex=lcx)
    
    
    par(fig=c(0.2,0.4,0,0.333),new=T)
    plot(exp(SDMacro)~Mchao,data=Merged,xlim=xlim_rich,ylim=ylim_macro,ann=F,axes=F,pch=HabPtsMacro,col=MMacroOutline,bg=MMacroCol,xpd=T)
    axis(1,at=xaxis_rich,cex.axis=acx,tck=-0.05)
    axis(2,at=yaxis_macro,labels=NA,cex.axis=acx,tck=tckl)
    mtext(bquote(Delta~'AICc' == .(round(dAIC_Rich_Macro,2))),1,line=dAICline,at=min(xlim_rich)+dAICfrac*diff(xlim_rich),cex=dAICexp)
    box(bty=BTY)
    if(dAIC_Rich_Macro < -2){
      plotMod <- list(lmRich_Macro1,lmRich_Macro2)[[which.min(c(AICc(lmRich_Macro1),AICc(lmRich_Macro2)))]]
      xvar<-Mtab$Mchao
      pred_x<-seq(min(xvar),max(xvar),length=npreds)
      pred_x_FO<-data.frame(Mchao=pred_x,Habitat=rep('FO',npreds),Site=NA)
      pred_x_FR<-data.frame(Mchao=pred_x,Habitat=rep('FR',npreds),Site=NA)
      pred_x_BA<-data.frame(Mchao=pred_x,Habitat=rep('BA',npreds),Site=NA)
      #list(Mdiv=rep(seq(1,3,length=100),1),Source=as.factor(c('LTER',rep('CRIOBE',100)))[-1],Habitat=as.factor(c('BA','FR',rep('FO',100)))[-c(1,2)])
      pred_y_FO<-exp(predict(plotMod,newdata=pred_x_FO,re.form=NA))
      pred_y_FR<-exp(predict(plotMod,newdata=pred_x_FR,re.form=NA))
      pred_y_BA<-exp(predict(plotMod,newdata=pred_x_BA,re.form=NA))
      
      points((pred_x),pred_y_FO,type='l',lwd=2,lty=1)
      points((pred_x),pred_y_FR,type='l',lwd=2,lty=2)
      points((pred_x),pred_y_BA,type='l',lwd=2,lty=3)
    }
    mtext('Fish richness',1,line=1.25,cex=lcx)
    
    
    
    par(fig=c(0.4,0.6,0,0.333),new=T)
    plot(exp(SDMacro)~Mdiv,data=Merged,xlim=xlim_div,ylim=ylim_macro,ann=F,axes=F,pch=HabPtsMacro,col=MMacroOutline,bg=MMacroCol,xpd=T)
    axis(1,at=xaxis_div,cex.axis=acx,tck=-0.05)
    axis(2,at=yaxis_macro,labels=NA,cex.axis=acx,tck=tckl)
    mtext(bquote(Delta~'AICc' == .(round(dAIC_Div_Macro,2))),1,line=dAICline,at=min(xlim_div)+dAICfrac*diff(xlim_div),cex=dAICexp)
    box(bty=BTY)
    
    if(dAIC_Div_Macro < -2){
      plotMod <- list(lmDiv_Macro1,lmDiv_Macro2)[[which.min(c(AICc(lmDiv_Macro1),AICc(lmDiv_Macro2)))]]
      xvar<-Mtab$Mdiv
      pred_x<-seq(min(xvar),max(xvar),length=npreds)
      pred_x_FO<-data.frame(Mdiv=pred_x,Habitat=rep('FO',npreds),Site=NA)
      pred_x_FR<-data.frame(Mdiv=pred_x,Habitat=rep('FR',npreds),Site=NA)
      pred_x_BA<-data.frame(Mdiv=pred_x,Habitat=rep('BA',npreds),Site=NA)
      #list(Mdiv=rep(seq(1,3,length=100),1),Source=as.factor(c('LTER',rep('CRIOBE',100)))[-1],Habitat=as.factor(c('BA','FR',rep('FO',100)))[-c(1,2)])
      pred_y_FO<-exp(predict(plotMod,newdata=pred_x_FO,re.form=NA))
      pred_y_FR<-exp(predict(plotMod,newdata=pred_x_FR,re.form=NA))
      pred_y_BA<-exp(predict(plotMod,newdata=pred_x_BA,re.form=NA))
      
      points((pred_x),pred_y_FO,type='l',lwd=2,lty=1)
      points((pred_x),pred_y_FR,type='l',lwd=2,lty=2)
      points((pred_x),pred_y_BA,type='l',lwd=2,lty=3)
    } #plot model if significant
    
    mtext('Fish diversity',1,line=1.25,cex=lcx)
  
    
    par(fig=c(0.6,0.8,0,0.333),new=T)
    plot(exp(SDMacro)~Mtl,data=Merged,xlim=xlim_tl,ylim=ylim_macro,ann=F,axes=F,pch=HabPtsMacro,col=MMacroOutline,bg=MMacroCol,xpd=T)
    axis(1,at=xaxis_tl,cex.axis=acx,tck=-0.05)
    axis(2,at=yaxis_macro,labels=NA,cex.axis=acx,tck=tckl)
    mtext(bquote(Delta~'AICc' == .(round(dAIC_TL_Macro,2))),1,line=dAICline,at=min(xlim_tl)+dAICfrac*diff(xlim_tl),cex=dAICexp)
    box(bty=BTY)
    
    if(dAIC_TL_Macro < -2){
      plotMod <- list(lmTL_Macro1,lmTL_Macro2)[[which.min(c(AICc(lmTL_Macro1),AICc(lmTL_Macro2)))]]
      xvar<-Mtab$Mtl
      pred_x<-seq(min(xvar),max(xvar),length=npreds)
      pred_x_FO<-data.frame(Mtl=pred_x,Habitat=rep('FO',npreds),Site=NA)
      pred_x_FR<-data.frame(Mtl=pred_x,Habitat=rep('FR',npreds),Site=NA)
      pred_x_BA<-data.frame(Mtl=pred_x,Habitat=rep('BA',npreds),Site=NA)
      #list(Mdiv=rep(seq(1,3,length=100),1),Source=as.factor(c('LTER',rep('CRIOBE',100)))[-1],Habitat=as.factor(c('BA','FR',rep('FO',100)))[-c(1,2)])
      pred_y_FO<-exp(predict(plotMod,newdata=pred_x_FO,re.form=NA))
      pred_y_FR<-exp(predict(plotMod,newdata=pred_x_FR,re.form=NA))
      pred_y_BA<-exp(predict(plotMod,newdata=pred_x_BA,re.form=NA))
      
      points((pred_x),pred_y_FO,type='l',lwd=2,lty=1)
      points((pred_x),pred_y_FR,type='l',lwd=2,lty=2)
      points((pred_x),pred_y_BA,type='l',lwd=2,lty=3)
    }
    #segments(2.95,-3.01+1.224*2.95,3.4,-3.01+1.224*3.4,lwd=1.5)
    #segments(2.95,(-2.91623+0.41822)+1.19485*2.95,3.4,(-2.91623+0.41822)+1.19485*3.4,lwd=1.5)
    mtext('Mean trophic level',1,line=1.25,cex=lcx)
    
    par(fig=c(0.8,1,0,0.333),new=T)
    plot(exp(SDMacro)~Mlmax,xlim=xlim_lmax,ylim=ylim_macro,data=Merged,ann=F,axes=F,pch=HabPtsMacro,col=MMacroOutline,bg=MMacroCol,xpd=T)
    axis(1,at=xaxis_lmax,cex.axis=acx,tck=-0.05)
    axis(2,at=yaxis_macro,labels=NA,cex.axis=acx,tck=tckl)
    mtext(bquote(Delta~'AICc' == .(round(dAIC_Lmax_Macro,2))),1,line=dAICline,at=min(xlim_lmax)+dAICfrac*diff(xlim_lmax),cex=dAICexp)
    if(dAIC_Lmax_Macro < -2){
      plotMod <- list(lmLmax_Macro1,lmLmax_Macro2)[[which.min(c(AICc(lmLmax_Macro1),AICc(lmLmax_Macro2)))]]
      xvar<-Mtab$Mlmax
      pred_x<-seq(min(xvar),max(xvar),length=npreds)
      pred_x_FO<-data.frame(Mlmax=pred_x,Habitat=rep('FO',npreds),Site=NA)
      pred_x_FR<-data.frame(Mlmax=pred_x,Habitat=rep('FR',npreds),Site=NA)
      pred_x_BA<-data.frame(Mlmax=pred_x,Habitat=rep('BA',npreds),Site=NA)
      #list(Mdiv=rep(seq(1,3,length=100),1),Source=as.factor(c('LTER',rep('CRIOBE',100)))[-1],Habitat=as.factor(c('BA','FR',rep('FO',100)))[-c(1,2)])
      pred_y_FO<-exp(predict(plotMod,newdata=pred_x_FO,re.form=NA))
      pred_y_FR<-exp(predict(plotMod,newdata=pred_x_FR,re.form=NA))
      pred_y_BA<-exp(predict(plotMod,newdata=pred_x_BA,re.form=NA))
      
      points((pred_x),pred_y_FO,type='l',lwd=2,lty=1)
      points((pred_x),pred_y_FR,type='l',lwd=2,lty=2)
      points((pred_x),pred_y_BA,type='l',lwd=2,lty=3)
    }
    box(bty=BTY)
    #mtext('Max length (mm)',1,line=1.25,cex=lcx)
    mtext(expression(paste(L[max],' (mm)')),1,line=1.25,cex=lcx)
    
    # mtext('J',2,outer=T,cex=lcx,las=1,line=pos1,at=0.31)
    # mtext('K',2,outer=T,cex=lcx,las=1,line=pos2,at=0.31)
    # mtext('L',2,outer=T,cex=lcx,las=1,line=pos3,at=0.31)
    # mtext('M',2,outer=T,cex=lcx,las=1,line=pos4,at=0.31)
    # mtext('N',2,outer=T,cex=lcx,las=1,line=pos5,at=0.31)
    
  } #Macro Models 
  
  if(TRUE){
    library(scales)
    par(fig=c(0,0.13,0.666,1),new=T)
    par(mar=c(0,0,0,0))
    plot(0,ylim=c(0,10),xlim=c(0,10),type='n',ann=F,axes=F)
    #legend(2.2,7,bty='n',legend=c('Fore','Back','Fringing'),pch=c(24,22,21),horiz=F,cex=lcx,border='black',xpd=T,y.intersp = 1)
    #legend(0,7,bty='n',seg.len=1,legend=c('','',''),lty=c(1,3,2),lwd=c(1,1,1),cex=lcx,horiz=F,xpd=T,y.intersp = 1)
    legend(-1.4,7,bty='n',legend=c('FO (solid)','BA (dotted)','FR (dashed)'),pch=c(24,22,21),horiz=F,cex=lcx,border='black',xpd=T,y.intersp = 1)
    
    
    par(fig=c(0.99,1,0.666,1),new=T)
    par(mar=c(0.1,0,0.1,0))
    par(mgp=c(0,0,0))
    ind1<-(length(GRS):1)/length(GRS)
    plot(1,type='n',ann=F,axes=F,xlim=c(0,1),ylim=c(0,1))
    polygon(x=c(0,0,1,1),y=c(0,1,1,0),col=rev(GRS)[1],border=rev(GRS)[1])
    for(i in 2:length(GRS)){
    	polygon(x=c(0,0,1,1),y=c(0,ind1[i],ind1[i],0),col=rev(GRS)[i],border=rev(GRS)[i])
    }
    polygon(x=c(0,0,1,1),y=c(0,1,1,0),border='black',lwd=1)
    #xlabb<-seq(0,180,length=length(histMbio$breaks))
    #xlabb[seq(2,length(histMbio$breaks),by=2)]<-''
    exp(histMbio$breaks)
    xlabb=c(0,5,20,90,250)
    axis(4,at=seq(0,1,length=length(histMbio$breaks))[c(1,3,6,9,11)],labels=xlabb,cex.axis=acx,tck=-0.5)
    #histMbio
    mtext(bquote('Fish biomass ' (g*'*'*m^-2)),4,line=0.75,cex=8/12,xpd=T)
    #mtext('Fish biomass',1,line=0.5,cex=8/12)
    
    par(fig=c(0.99,1,0.333,0.666),new=T)
    #par(mar=c(0,0,0,0))
    par(mgp=c(0,0,0))
    ind1<-(length(BUS):1)/length(BUS)
    plot(1,type='n',ann=F,axes=F,xlim=c(0,1),ylim=c(0,1))
    polygon(x=c(0,0,1,1),y=c(0,1,1,0),col=rev(BUS)[1],border=rev(BUS)[1])
    for(i in 2:length(BUS)){
    	polygon(x=c(0,0,1,1),y=c(0,ind1[i],ind1[i],0),col=rev(BUS)[i],border=rev(BUS)[i])
    }
    polygon(x=c(0,0,1,1),y=c(0,1,1,0),border='black',lwd=1)
    #axis(1,at=c(0,0.25,0.5,0.75,1),labels=c(0,'',150,'',300),cex.axis=acx,tck=-0.5)
    axis(4,at=c(0,0.25,0.5,0.75,1),labels=c(0,'',20,'',40),cex.axis=acx,tck=-0.5)
    
    mtext('Coral cover (%)',4,line=0.5,cex=8/12)
    
    par(fig=c(0.99,1,0,0.333),new=T)
    #par(mar=c(0,0,0,0))
    par(mgp=c(0,0,0))
    ind1<-(length(RDS):1)/length(RDS)
    plot(1,type='n',ann=F,axes=F,xlim=c(0,1),ylim=c(0,1))
    polygon(x=c(0,0,1,1),y=c(0,1,1,0),col=rev(RDS)[1],border=rev(RDS)[1])
    for(i in 2:length(RDS)){
    	polygon(x=c(0,0,1,1),y=c(0,ind1[i],ind1[i],0),col=rev(RDS)[i],border=rev(RDS)[i])
    }
    polygon(x=c(0,0,1,1),y=c(0,1,1,0),border='black',lwd=1)
    #axis(1,at=c(0,0.25,0.5,0.75,1),labels=c(0,'',150,'',300),cex.axis=acx,tck=-0.5)
    axis(4,at=c(0,0.25,0.5,0.75,1),labels=c(0,'',10,'',20),cex.axis=acx,tck=-0.5)
    mtext('Macroalgae cover (%)',4,line=0.5,cex=8/12)
    
  } #legends
  
}
  