rm(list=ls())

library(dplyr)
library(MuMIn)
library(RColorBrewer)
library(here)


load(here('Data','CoralReturnRateTable.Rdata'))

for(TrophicGroup in c('FishComm','PrimCons','Herbivore','ScraperExcavator','Grazer','Browser','Corallivore')){
  #TrophicGroup <- 'FishComm'
  load(here('Data',paste0('TSregTable_',TrophicGroup,'.Rdata')))

  TSregTable$MacroCover[is.na(TSregTable$MacroCover)] <- 0 #no macro observed
  TSregTable <- TSregTable %>% filter(Year < 2018 & Year > 2005)
  TSregTable$Bio.m[is.na(TSregTable$Bio.m) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  TSregTable$Rich[is.na(TSregTable$Rich) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  TSregTable$CHAO[is.na(TSregTable$CHAO) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  TSregTable$ACE[is.na(TSregTable$ACE) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  TSregTable$FishBiomDiv[is.na(TSregTable$FishBiomDiv) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  TSregTable$MeanTL[is.na(TSregTable$MeanTL) & !(TSregTable$Year %in% c(2005,2018))] <- 0
  TSregTable$Lmax[is.na(TSregTable$Lmax) & !(TSregTable$Year %in% c(2005,2018))] <- 0

  CV<-function(x){
    return(sd(x,na.rm=T)/mean(x,na.rm=T))
  }
  logit<-function(x){
  	return(log(x/(1-x)))
  }
  
  TSregTable <- TSregTable %>% mutate(HabSite = paste(Habitat,Site,sep='.')) 
  ObservedSites <- TSregTable %>% group_by(HabSite) %>% summarize(ObservedBool = (sum(Bio.m)>0 & sum(Bio.m>0)>=2)) %>% filter(ObservedBool) %>% pull(HabSite)
  TSregTable <- TSregTable %>% filter(HabSite %in% ObservedSites)


  Mtab<-TSregTable %>% group_by(Habitat,Site) %>% summarize(MCoral = mean(logit(CoralCover+0.01),na.rm=T),MCrich=mean(CoralRichness,na.rm=T),MCdiv=mean(CoralDiversity,na.rm=T),MMacro = mean(logit(MacroCover+0.01),na.rm=T),Mbiomass=mean(log(Bio.m+0.01),na.rm=T),Mchao=mean(CHAO,na.rm=T), Mrich=mean(Rich,na.rm=T),Mdiv=mean(FishBiomDiv,na.rm=T),Mlmax=mean(Lmax,na.rm=T),Mtl=mean(MeanTL,na.rm=T),Source=unique(Source))
  Merged <- CoralReturnRateTable  %>% left_join(Mtab,by=c('Habitat','Site','Source'))  %>% filter(Habitat=='FO')
  Merged$PreCover <- as.numeric(Merged$PreCover)
  Merged$PostCover <- as.numeric(Merged$PostCover)

  Merged$lPreCover <- logit(as.numeric(Merged$PreCover)/100)
  Merged$lPostCover <- logit(as.numeric(Merged$PostCover)/100)


  lm.base<-lm(RR1017 ~ 1,data=Merged)
  lm.source<-lm(RR1017 ~ Source,data=Merged)
  AICc(lm.base);AICc(lm.source)

  #Merged[which(Merged$RR1017 > 8),]
  
  lm.biomass <- lm(RR1017 ~ Mbiomass,data=Merged)
  AICc(lm.base);AICc(lm.biomass)
  dAICc_biomass<-min(AICc(lm.biomass))-AICc(lm.base)
  
  lm.richness <- lm(RR1017 ~ Mchao,data=Merged)
  AICc(lm.base);AICc(lm.richness)
  dAICc_richness<-min(AICc(lm.richness))-AICc(lm.base)
  
  lm.diversity <- lm(RR1017 ~ Mdiv,data=Merged) #kinda sorta
  AICc(lm.base);AICc(lm.diversity)
  dAICc_diversity<-min(AICc(lm.diversity))-AICc(lm.base)
  lm.diversity
  
  lm.tl <- lm(RR1017 ~ Mtl,data=Merged)
  AICc(lm.base);AICc(lm.tl)
  dAICc_tl<-min(AICc(lm.tl))-AICc(lm.base)
  
  lm.lmax <- lm(RR1017 ~ Mlmax,data=Merged)
  AICc(lm.base);AICc(lm.lmax)
  dAICc_lmax<-min(AICc(lm.lmax))-AICc(lm.base)

  ###### PreCover
  lm.biomass.cover <- lm(RR1017 ~ Mbiomass + lPreCover,data=Merged)
  AICc(lm.base);AICc(lm.biomass.cover)
  dAICc_biomass.cover<-min(AICc(lm.biomass.cover))-AICc(lm.base)
  
  lm.richness.cover <- lm(RR1017 ~ Mchao* lPreCover,data=Merged)
  AICc(lm.base);AICc(lm.richness.cover)
  dAICc_richness.cover<-min(AICc(lm.richness.cover))-AICc(lm.base)
  
  lm.diversity.cover <- lm(RR1017 ~ Mdiv * lPreCover,data=Merged) #kinda sorta
  AICc(lm.base);AICc(lm.diversity.cover)
  dAICc_diversity.cover<-min(AICc(lm.diversity.cover))-AICc(lm.base)
  
  lm.tl.cover <- lm(RR1017 ~ Mtl * lPreCover,data=Merged)
  AICc(lm.base);AICc(lm.tl.cover)
  dAICc_tl.cover<-min(AICc(lm.tl.cover))-AICc(lm.base)
  
  lm.lmax.cover <- lm(RR1017 ~ Mlmax * lPreCover,data=Merged)
  AICc(lm.base);AICc(lm.lmax.cover)
  dAICc_lmax.cover<-min(AICc(lm.lmax.cover))-AICc(lm.base)
  
  AICc(lm.base);AICc(lm.biomass);AICc(lm.biomass.cover)

}#


if(TrophicGroup == 'FishComm'){
quartz(width=7,height=1.85)
pos1<- -0.9
  pos2 <- -7
  pos3 <- -13.1
  pos4 <- -19.2
  pos5 <- -25.3
if(TRUE){
  
  
  par(oma=c(2.25,2.5,0.5,2.5))
  par(mar=rep(0,4))
  par(lheight=0.85)
  par(mgp=c(0,0.5,0))
  
  tckl<-0
  acx<-7/12
  lcx<-9/12
  Transparency<- '55'
  
  bottom<-0.25
  range(exp(Merged$Mbiomass))
  range(Merged$Mchao)
  range(Merged$Mdiv)
  range(Merged$Mtl)
  range(Merged$Mlmax)
  
  xlim_bio<-c(0,180)
  xlim_rich<-c(30,80)
  xlim_div <-c(2.2,3)
  xlim_tl <-c(3,3.4)
  xlim_lmax <-c(140,190)
  
  xaxis_bio<-c(0,60,120,180)
  xaxis_rich<-c(30,40,50,60,70)
  xaxis_div <-c(2.2,2.4,2.6,2.8)
  xaxis_tl <-c(3,3.1,3.2,3.3)
  xaxis_lmax <-c(140,150,160,170,180)
  
  ylim_return<-c(0,10)
    
  yaxis_return<-c(0,2,4,6,8,10)
  
  
  
  
  dAICline<--1
  dAICexp<-7/12
  dAICfrac<-0.70
  
  npreds<-300
  
  BTY <-'l'
}

if(TRUE){
BluePal<-brewer.pal(9,'Blues')[c(3,6,9)]#c('lightblue','dodgerblue','darkblue')
GreenPal<-brewer.pal(9,'Greens')[c(3,6,9)]

histMcoral<-hist(Merged$PreCover,breaks=10,plot=F)
BuPal<-colorRampPalette(brewer.pal(9,'Blues'))
BUS<-BuPal(length(histMcoral$mids))
MCoralCol<-BUS[apply(matrix(Merged$PreCover,ncol=1),1,FUN=function(x){which.min(Mod(x-histMcoral$mids))})]
MCoralOutline<-'black'
#MCoralOutline[Mtab$Source=='LTER']<-'black'

par(fig=c(0,0.2,0,1))
plot(RR1017~exp(Mbiomass),data=Merged  %>% select(RR1017,Mbiomass) %>% na.omit(),ylim=ylim_return,xlim=xlim_bio,ann=F,axes=F,pch=24,bg=MCoralCol,col=MCoralOutline)
axis(1,at=xaxis_bio,labels=c('0','60','120',''),cex.axis=acx,tck=-0.05)
axis(2,at=yaxis_return,cex.axis=acx,tck=-0.05)
#box(bty='l')
box()
#mtext('Fish biomass (g*m2)',1,line=1.25,cex=lcx)
#mtext(bquote(Delta~'AICc' == .(round(dAICc_biomass,2))),3,outer=T,line=-1.1,at=0.07,cex=dAICcexp)
mtext(bquote(Delta~'AICc' == .(round(dAICc_biomass,2))),1,line=dAICline,at=min(xlim_bio)+dAICfrac*diff(xlim_bio),cex=dAICexp)
  
mtext(bquote('Fish biomass ' (g*'*'*m^-2)),1,line=1.35,cex=lcx)
mtext(bquote('Coral return rate '('%'*'*'*year^-1)),2,line=1.1,at=4,cex=lcx)

par(fig=c(0.2,0.4,0,1),new=T)
plot(RR1017~Mchao,data=Merged  %>% select(RR1017,Mchao) %>% na.omit(),ylim=ylim_return,xlim=xlim_rich,ann=F,axes=F,pch=24,bg=MCoralCol,col=MCoralOutline)
axis(1,at=xaxis_rich,cex.axis=acx,tck=-0.05)
axis(2,at=yaxis_return,labels=NA,cex.axis=acx,tck=tckl)
mtext(bquote(Delta~'AICc' == .(round(dAICc_richness,2))),1,line=dAICline,at=min(xlim_rich)+dAICfrac*diff(xlim_rich),cex=dAICexp)
  
box()
mtext('Fish richness',1,line=1.2,cex=lcx)

par(fig=c(0.4,0.6,0,1),new=T)
plot(RR1017~Mdiv,data=Merged  %>% select(RR1017,Mdiv) %>% na.omit(),ylim=ylim_return,xlim=xlim_div,ann=F,axes=F,pch=24,bg=MCoralCol,col=MCoralOutline)
axis(1,at=xaxis_div,cex.axis=acx,tck=-0.05)
axis(2,at=yaxis_return,labels=NA,cex.axis=acx,tck=tckl)
box()
mtext(bquote(Delta~'AICc' == .(round(dAICc_diversity,2))),1,line=dAICline,at=min(xlim_div)+dAICfrac*diff(xlim_div),cex=dAICexp)
mtext('Fish diversity',1,line=1.2,cex=lcx)

par(fig=c(0.6,0.8,0,1),new=T)
plot(RR1017~Mtl,data=Merged  %>% select(RR1017,Mtl) %>% na.omit(),ylim=ylim_return,xlim=xlim_tl,ann=F,axes=F,pch=24,bg=MCoralCol,col=MCoralOutline)
axis(1,at=xaxis_tl,cex.axis=acx,tck=-0.05)
axis(2,at=yaxis_return,labels=NA,cex.axis=acx,tck=tckl)
mtext(bquote(Delta~'AICc' == .(round(dAICc_tl,2))),1,line=dAICline,at=min(xlim_tl)+dAICfrac*diff(xlim_tl),cex=dAICexp)
  
box()
mtext('Mean trophic level',1,line=1.2,cex=lcx)

par(fig=c(0.8,1,0,1),new=T)
plot(RR1017~Mlmax,data=Merged  %>% select(RR1017,Mlmax) %>% na.omit(),ylim=ylim_return,xlim=xlim_lmax,ann=F,axes=F,pch=24,bg=MCoralCol,col=MCoralOutline)
axis(1,at=xaxis_lmax,cex.axis=acx,tck=-0.05)
axis(2,at=yaxis_return,labels=NA,cex.axis=acx,tck=tckl)
mtext(bquote(Delta~'AICc' == .(round(dAICc_lmax,2))),1,line=dAICline,at=min(xlim_lmax)+dAICfrac*diff(xlim_lmax),cex=dAICexp)
box()
summary(lm.lmax)
#segments(x0=142,x1=181,y0=-10.09203 + 0.08843*142, y1=-10.09203 + 0.08843*181,lty=2,lwd=1.5)
mtext(expression(paste(L[max],' (mm)')),1,line=1.35,cex=lcx)

#mtext('A',2,outer=T,cex=lcx,las=1,line=pos1,at=0.94)
#mtext('B',2,outer=T,cex=lcx,las=1,line=pos2,at=0.94)
#mtext('C',2,outer=T,cex=lcx,las=1,line=pos3,at=0.94)
#mtext('D',2,outer=T,cex=lcx,las=1,line=pos4,at=0.94)
#mtext('E',2,outer=T,cex=lcx,las=1,line=pos5,at=0.94)


  par(fig=c(0.99,1,0,1),new=T)
  par(mar=c(0.5,0,0.5,0))
  par(mgp=c(0,0,0))
  ind1<-(length(BUS):1)/length(BUS)
  plot(1,type='n',ann=F,axes=F,xlim=c(0,1),ylim=c(0,1))
  polygon(x=c(0,0,1,1),y=c(0,1,1,0),col=rev(BUS)[1],border=rev(BUS)[1])
  for(i in 2:length(BUS)){
  	polygon(x=c(0,0,1,1),y=c(0,ind1[i],ind1[i],0),col=rev(BUS)[i],border=rev(BUS)[i])
  }
  polygon(x=c(0,0,1,1),y=c(0,1,1,0),border='black',lwd=1)
  #axis(1,at=c(0,0.25,0.5,0.75,1),labels=c(0,'',150,'',300),cex.axis=acx,tck=-0.5)
  axis(4,at=c(0,0.25,0.5,0.75,1),labels=c(20,'',34,'',48),cex.axis=acx,tck=-0.5)
  
  mtext('Pre-disturbance \n coral cover (%)',4,line=1.25,cex=lcx)

}
