rm(list=ls())

library(TMB)
library(dplyr)
library(tidyr)
library(doParallel)
library(here)

logit<-function(x){
	return(log(x/(1-x)))
}
iLogit<-function(x){
	return(exp(x)/(1+exp(x)))
}
Zscore<-function(x){
	return((x-mean(x,na.rm=T))/sd(x,na.rm=T))
}

#Custom AICc for model fits
AICc<-function(opt1,data,tAIC=F){
	k<-length(opt1$par)
	n<-length(na.omit(as.vector(data)))
	AIC<-2*k + 2*opt1$objective;AIC
	AICc <- AIC + (2*k^2+2*k)/(n-k-1)
	if((n-k)==1){return(999)}
	if(tAIC){return(AIC)
	}else{return(AICc)}
	
}

#custom plotting function
plotModFit<-function(mod){
	objx<-obj1
	pdf(paste0(modName,'.pdf'))
		for(j in 1:nrow(objx$data$obs)){
			par(mfrow=c(2,1))
			plot(2006:2017,objx$data$obs[j,],main=row.names(objx$data$obs)[j])
			points(2006:2017,c(rep(NA,objx$data$Kmax),objx$report()$pred[j,]),type='l')
			abline(v=2006+Kmax,col='red',lty=2,lwd=2)
			#points(2006:2017,Null1$obj1$report()$pred[j,],type='l',col='red',lwd=2)
			#points(2006:2016,Null1$obj1$report()$pred[j,2:12],type='l',col='dodgerblue',lwd=2)
			plot(objx$data$obs[j,2:12]~objx$data$Covars[j,1:11,1])
		}
	dev.off()
}

#extract parameter estimates and other model fit statisticas
extractModChar<-function(mod){
	#mod<-ModRuns[[1]]
	#dyn.load(TMB::dynlib("Analysis/regression_armaE"))
	pl1<-mod$env$parList()
	
	parnames<-names(mod$errPars$value)
	sds<-mod$errPars$sd
	BetaHab <- array(sds[which(parnames == 'Beta_hab')],dim=c(3,mod$modelForm$Bt+mod$modelForm$Kbeta,1))
	BetaRh <- sds[which(parnames == 'Beta_Rh')]
	#mlePars<-mod$report()
	d1<-data.frame(CommAttr=paste0(CommAttr[[mod$modelForm$i]],collapse='_'),
		K=mod$modelForm$K,
		Kbeta=mod$modelForm$Kbeta,
		Bt=mod$modelForm$Bt,
		Kmax=mod$modelForm$Kmax,
		Bool_Source=mod$modelForm$Bool_Source,
		dump=mod$modelForm$dump,
		Coef_BA=round(pl1$Beta_hab[1,,1],4),
		Coef_FO=round(pl1$Beta_hab[2,,1],4),
		Coef_FR=round(pl1$Beta_hab[3,,1],4),
		CoefR_BA=round(pl1$Beta_Rh[1],4),
		CoefR_FO=round(pl1$Beta_Rh[2],4),
		CoefR_FR=round(pl1$Beta_Rh[3],4),
		Err_BA=round(BetaHab[1,,1],4),
		Err_FO=round(BetaHab[2,,1],4),
		Err_FR=round(BetaHab[3,,1],4),
		ErrR_BA=round(BetaRh[1],4),
		ErrR_FO=round(BetaRh[2],4),
		ErrR_FR=round(BetaRh[3],4),
		AIC=round(AICc(mod$opt,data=mod$data$obs[,-c(1:mod$modelForm$Kmax)],tAIC=T),2),
		AICc=round(AICc(mod$opt,data=mod$data$obs[,-c(1:mod$modelForm$Kmax)],tAIC=F),2),
		Converge=ifelse(mod$opt$convergence==0,TRUE,FALSE),
		stringsAsFactors=FALSE
	)
	d1
}

compile(here('Code',"benthic_gompertz2.cpp"))
dyn.load(TMB::dynlib(here('Code',"benthic_gompertz2")))

#read in data
FunctionalGroupList <- c('FishComm','PrimCons','Herbivore','Grazer','ScraperExcavator','Browser','Corallivore')

for(ff in 1:length(FunctionalGroupList)){
  #ff<-1
  FunctionalGroup <- FunctionalGroupList[ff]
  load(here('Data',paste0('TSregTable_',FunctionalGroup,'.Rdata')))
  
  #Missing values need to be replaced with 0's, as they are just non-detections
  TSregTable$MacroCover[is.na(TSregTable$MacroCover)] <- 0 #no macro observed
  TSregTable <- TSregTable %>% filter(Year < 2018 & Year > 2005)
  TSregTable$Bio.m[is.na(TSregTable$Bio.m) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  TSregTable$Rich[is.na(TSregTable$Rich) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  TSregTable$CHAO[is.na(TSregTable$CHAO) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  TSregTable$ACE[is.na(TSregTable$ACE) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  TSregTable$FishBiomDiv[is.na(TSregTable$FishBiomDiv) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  TSregTable$MeanTL[is.na(TSregTable$MeanTL) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  TSregTable$Lmax[is.na(TSregTable$Lmax) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  
  
  TSregTable_mod <- TSregTable  %>% mutate(HabSite=paste0(Habitat,Site),HabSiteYear=paste0(Habitat,Site,Year))  %>% filter(Year <2018 & Year>2005)%>% select(-MeanLength,-MeanMass)
  TSregTable_mod$interpolate <- F
  
  #Zscore the data
  if(TRUE){
    uniqueHabSite <- TSregTable_mod %>% pull(HabSite) %>% unique() %>% sort()
  	TSregTable_mod$PropCover <- TSregTable_mod$MacroCover
  	TSregTable_mod$lPropCover <- (TSregTable_mod$MacroCover + 0.01) %>% log()
  	
  	TSregTable_mod$Z_lPropCover<-TSregTable_mod$lPropCover
  	TSregTable_mod$Z_lPropCover[TSregTable_mod$Source=='CRIOBE'] <- Zscore(TSregTable_mod$lPropCover[TSregTable_mod$Source=='CRIOBE'])
  	TSregTable_mod$Z_lPropCover[TSregTable_mod$Source=='LTER'] <- Zscore(TSregTable_mod$lPropCover[TSregTable_mod$Source=='LTER'])
  
  	TSregTable_mod$lBio.m<-(TSregTable_mod$Bio.m+0.1) %>% log()
  	TSregTable_mod$Z_lBio.m <- TSregTable_mod$lBio.m
  	TSregTable_mod$Z_lBio.m[TSregTable_mod$Source=='CRIOBE'] <- Zscore(log(TSregTable_mod$Bio.m[TSregTable_mod$Source=='CRIOBE']+0.01))
  	TSregTable_mod$Z_lBio.m[TSregTable_mod$Source=='LTER'] <- Zscore(log(TSregTable_mod$Bio.m[TSregTable_mod$Source=='LTER']+0.01))
  
  	TSregTable_mod$Z_CHAO<-(TSregTable_mod$CHAO)
  	TSregTable_mod$Z_CHAO[TSregTable_mod$Source=='CRIOBE'] <- Zscore(TSregTable_mod$CHAO[TSregTable_mod$Source=='CRIOBE'])
  	TSregTable_mod$Z_CHAO[TSregTable_mod$Source=='LTER'] <- Zscore(TSregTable_mod$CHAO[TSregTable_mod$Source=='LTER'])

  	TSregTable_mod$Z_FishBiomDiv<-(TSregTable_mod$FishBiomDiv)
  	TSregTable_mod$Z_FishBiomDiv[TSregTable_mod$Source=='CRIOBE'] <- Zscore(TSregTable_mod$FishBiomDiv[TSregTable_mod$Source=='CRIOBE'])
  	TSregTable_mod$Z_FishBiomDiv[TSregTable_mod$Source=='LTER'] <- Zscore(TSregTable_mod$FishBiomDiv[TSregTable_mod$Source=='LTER'])
  
  	TSregTable_mod$Z_Lmax<-TSregTable_mod$Lmax
  	TSregTable_mod$Z_Lmax[TSregTable_mod$Source=='CRIOBE'] <- Zscore(TSregTable_mod$Lmax[TSregTable_mod$Source=='CRIOBE'])
  	TSregTable_mod$Z_Lmax[TSregTable_mod$Source=='LTER'] <- Zscore(TSregTable_mod$Lmax[TSregTable_mod$Source=='LTER'])
  
  	TSregTable_mod$Z_MeanTL<-TSregTable_mod$MeanTL
  	TSregTable_mod$Z_MeanTL[TSregTable_mod$Source=='CRIOBE'] <- Zscore(TSregTable_mod$MeanTL[TSregTable_mod$Source=='CRIOBE'])
  	TSregTable_mod$Z_MeanTL[TSregTable_mod$Source=='LTER'] <- Zscore(TSregTable_mod$MeanTL[TSregTable_mod$Source=='LTER'])
  }
  
  #What Community attributes do you want to model
  CommAttr<<-list('NoCovar',
  				'Bio.m',
  				'Rich',
  				'FishBiomDiv',
  				'Lmax',
  				'MeanTL')
  
  #Create full list of possibel model forms
  #Bt= zero lag
  #K=year lag for benthic dynamic
  #Bool_Source = source effect, deprecated because we z-score within source
  #Kbeta = lag of fish effect
  #i = which fish community attributes
  modelForms <- expand.grid(Bt=0:0,K=1:1,Bool_Source=0:0,Kbeta=1:3,i=1:length(CommAttr))
  BadModelForms<-which(modelForms$Bt==0 & modelForms$Kbeta==0)
  if(length(BadModelForms)>0){
  	modelForms <- modelForms[-which(modelForms$Bt==0 & modelForms$Kbeta==0),]
  }
  modelForms$dump=FALSE
  modelForms$Kmax<-max(modelForms$K,modelForms$Kbeta)
  
  #PARALLEL MODEL CODE
  #cl<-makeCluster(4)
  #registerDoParallel(cl)
  ModRuns <- foreach(SimNum=1:nrow(modelForms),.packages=c('dplyr')) %do% {
    
    #Set variables for each run
  	Bt<-modelForms$Bt[SimNum]
  	K<-modelForms$K[SimNum]
  	Bool_Source<-modelForms$Bool_Source[SimNum]
  	Kbeta<-modelForms$Kbeta[SimNum]
  	i<-modelForms$i[SimNum]
  	dump<-modelForms$dump[SimNum]
  	Kmax<-modelForms$Kmax[SimNum]
  	
  	#build data sets for each run
  	BenthicIn <- TSregTable_mod  %>% xtabs(Z_lPropCover ~ HabSite+Year,data=.)
  	Habs<-as.numeric(factor(substr(row.names(BenthicIn),1,2)))-1
  	Sites<-as.numeric(factor(substr(row.names(BenthicIn),3,nchar(row.names(BenthicIn)))))-1
  	Source<-as.numeric(factor(substr(row.names(BenthicIn),nchar(row.names(BenthicIn))-3,nchar(row.names(BenthicIn)))))-1
  	
  	thisCommAttr<-CommAttr[[i]]
  	FishIn <- array(0,dim=c(nrow(BenthicIn),ncol(BenthicIn),length(thisCommAttr)))
  	RecoveryCovar <- rep(0,nrow(BenthicIn))
  	if('Bio.m' %in% thisCommAttr){
  	  mat1 <- TSregTable_mod  %>% xtabs(Z_lBio.m~HabSite+Year,data=.)
  	  mat2 <- TSregTable_mod  %>% xtabs(lBio.m~HabSite+Year,data=.)
  	  FishIn[,,which(thisCommAttr=='Bio.m')] <- mat1
  	  RecoveryCovar <- rowMeans(mat2)}
  	if('Rich' %in% thisCommAttr){
  	  mat1<-TSregTable_mod  %>% xtabs(Z_CHAO~HabSite+Year,data=.)
  	  mat2 <- TSregTable_mod  %>% xtabs(CHAO~HabSite+Year,data=.)
  	  FishIn[,,which(thisCommAttr=='Rich')] <- mat1
  	  RecoveryCovar <- rowMeans(mat2)}
  	if('FishBiomDiv' %in% thisCommAttr){
  	  mat1<-TSregTable_mod  %>% xtabs(Z_FishBiomDiv~HabSite+Year,data=.)
  	  mat2 <- TSregTable_mod  %>% xtabs(FishBiomDiv~HabSite+Year,data=.)
  	  FishIn[,,which(thisCommAttr=='FishBiomDiv')] <- mat1
  	  RecoveryCovar <- rowMeans(mat2)}
  	if('Lmax' %in% thisCommAttr){
  	  mat1 <- TSregTable_mod  %>% xtabs(Z_Lmax~HabSite+Year,data=.)
  	  mat2 <- TSregTable_mod  %>% xtabs(Lmax~HabSite+Year,data=.)
  	  FishIn[,,which(thisCommAttr=='Lmax')] <- mat1
  	  RecoveryCovar <- rowMeans(mat2)}
  	if('MeanTL' %in% thisCommAttr){
  	  mat1<-TSregTable_mod  %>% xtabs(Z_MeanTL~HabSite+Year,data=.)
  	  mat2 <- TSregTable_mod  %>% xtabs(MeanTL~HabSite+Year,data=.)
  	  FishIn[,,which(thisCommAttr=='MeanTL')] <- mat1
  	  RecoveryCovar <- rowMeans(mat2)}
  	
  	Ncovar<-dim(FishIn)[3]
    
  	
  	#build data for model input  
  	Nt<-nrow(BenthicIn)	
  	Nhabs<-length(unique(Habs))
  	NSource<-length(unique(Source))
  	dataIn=list('obs'=BenthicIn,'Covars'=FishIn,Habs=Habs,K=K,Kbeta=Kbeta,Kmax=Kmax,Bt=Bt,R_Covar=RecoveryCovar)
  	parametersIn=list(
  					R_hab=rep(0,Nhabs),
  					lSigma_R=0,
  					R_sh = rep(1,Nt),
  					Beta_Rh = rep(0,Nhabs),
  					Beta_hab=array(0,dim=c(Nhabs,Bt+Kbeta,Ncovar)),
  					Ne=matrix(0.0,nrow=Nt,ncol=ncol(BenthicIn)),
  					Phi_h=rep(0,Nhabs),
  					logsdPro=rep(0,Nt),
  					logsdObs=0)
  	  
  	#model specification if no covariates added
  		if(CommAttr[[i]] != 'NoCovar'){
  		  map1 <- list(Beta_Rh = factor(rep(NA,Nhabs)))
  			obj1 <- TMB::MakeADFun(dataIn,parametersIn,random=c('R_sh','Ne'),DLL="benthic_gompertz2",silent=T,map=map1)
  			lowerBounds=c(
  						  rep(-Inf,Nhabs),	#R_hab
  						  -Inf,	#lSigma_R
  						  rep(-Inf,Nhabs*(Bt+Kbeta)*Ncovar),#Beta_Rh*(Bt+Kbeta)*Ncovar
  					    rep(-Inf,Nhabs),	#Phi_h
  						  rep(-Inf,Nt), #logsdPro
  						  -Inf)		#logsdObs
  			upperBounds=c(
  						  rep(Inf,Nhabs),	#R_hab
  						  Inf,	#lSigma_R
  						  rep(Inf,Nhabs*(Bt+Kbeta)*Ncovar),#Beta_Rh*(Bt+Kbeta)*Ncovar
  					    rep(Inf,Nhabs),	#Phi_h
  						  rep(Inf,Nt), #logsdPro
  						  Inf)		#logsdObs
  			if(length(lowerBounds) != length(upperBounds) | length(lowerBounds) != length(obj1$par)) print('YOU HAVE A BOUNDS PROBLEM')
  		}else{ #model spects in covariates added
  				map1 <- list(Beta_Rh = factor(rep(NA,Nhabs)),Beta_hab = factor(array(NA,dim=c(Nhabs,Bt+Kbeta,Ncovar))))
  				obj1 <- TMB::MakeADFun(dataIn,parametersIn,random=c('R_sh','Ne'),DLL="benthic_gompertz2",silent=T,map=map1)
  				lowerBounds=c(
  						  rep(-Inf,Nhabs),	#R_hab
  						  -Inf,	#lSigma_R
  					    rep(-Inf,Nhabs),	#Phi_h
  						  rep(-Inf,Nt), #logsdPro
  						  -Inf)		#logsdObs
  			upperBounds=c(
  						  rep(Inf,Nhabs),	#R_hab
  						  Inf,	#lSigma_R
  					    rep(Inf,Nhabs),	#Phi_h
  						  rep(Inf,Nt), #logsdPro
  						  Inf)		#logsdObs
  			if(length(lowerBounds) != length(upperBounds) | length(lowerBounds) != length(obj1$par)) print('YOU HAVE A BOUNDS PROBLEM')
  		}
  	
  	obj1$control=list(reltol=1e-12,maxit=5000)
  	obj1$fn()
  	obj1$gr()
  	obj1$method='L-BFGS-B'
  	#obj1$par=opt1$par
  	obj1$lower=lowerBounds
  	obj1$upper=upperBounds
  	#system.time(opt1 <- do.call("optim",obj1))
  	opt1 <- do.call("optim",obj1)
  	opt1$objective <-opt1$value
  	
  	#opt1 <- nlminb(obj1$par,obj1$fn,obj1$gr,lower=lowerBounds,upper=upperBounds,control=list(iter.max=2000,eval.max=2000))
  
  	pl1 <- obj1$env$parList()#This contains all of your parameter estimates RAW as they come out of the optimizer
  	mlePars<-obj1$report()
  	errPars<-sdreport(obj1)
  	
  	obj1$data<-dataIn
  
  	obj1$modelForm<-modelForms[SimNum,]
  	obj1$Regressor <-paste(FunctionalGroup,paste0(CommAttr[[i]]),sep='_')
  	obj1$opt<-opt1
  	obj1$errPars <- errPars
  	obj1$mlePars <- mlePars
  	obj1
  }
  #stopCluster(cl)
  
  
  #MODEL SELECTION TABLE
  ModelSelectionTable<-lapply(ModRuns,FUN=extractModChar) %>% bind_rows()  %>% mutate(dAICc=AICc-min(AICc))
 
  saveRDS(ModelSelectionTable,file=here('Output',paste0('Macro_Gompertz_',FunctionalGroupList[ff],'_InterannualFish.RDS')))
}

