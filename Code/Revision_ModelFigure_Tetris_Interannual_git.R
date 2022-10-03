rm(list=ls())
#ModelFigure
library(dplyr)
library(RColorBrewer)
library(TMB)
library(here)

load(here('Output','Revision_ModelTables_SIG.Rdata'))
TrophicGroupListI_SIG <- TrophicGroupListI
Macro_TrophicGroupListI_SIG <- Macro_TrophicGroupListI
load(here('Output','Revision_ModelTables.Rdata'))
#load('ModelFigureSetup.Rdata')


BiomassCol<-brewer.pal(9,'Greens')[7]
RichnessCol<-brewer.pal(9,'Greys')[7]
DiversityCol<-brewer.pal(9,'Purples')[7]
TrophicCol<-brewer.pal(9,'Oranges')[6]
LmaxCol<-brewer.pal(9,'Blues')[7]

Cols<-c(BiomassCol,RichnessCol,DiversityCol,LmaxCol,TrophicCol)
ColorNames<-c('Bio.m','Rich','FishBiomDiv','Lmax','MeanTL')
CovarNames <- ColorNames
TrophicGroups<- c('FishComm','PrimCons','Herbivore','Grazer','ScraperExcavator','Browser','Corallivore')


quartz(width=3.5,height=5.75)
par(mar=c(0.5,1.2,0.5,1.2))
par(oma=c(1,1,1,1))
gridLWD<-1

for(ff in 1:7){
  
  if(ff!=1) par(new=T)
  
  par(fig=c(0,0.5,(7-ff)/7,(7-(ff-1))/7))
  if(ff==1){
    plot(0,type='n',xlim=c(0,5),ylim=c(0,3),ann=F,axes=F,yaxs='i',xaxs='i')
  }else{
    plot(0,type='n',xlim=c(0,5),ylim=c(0,3),ann=F,axes=F,yaxs='i',xaxs='i')
  }
  
  ccmax <- ifelse(ff==1,5,4)
  for(cc in 1:ccmax){
      ThisColor<-Cols[cc]
      
      col1<-ifelse(!is.na(TrophicGroupListI_SIG[[ff]][1,cc]),paste0(ThisColor),'#FFFFFF00')
      col2<-ifelse(!is.na(TrophicGroupListI_SIG[[ff]][2,cc]),paste0(ThisColor),'#FFFFFF00')
      col3<-ifelse(!is.na(TrophicGroupListI_SIG[[ff]][3,cc]),paste0(ThisColor),'#FFFFFF00')
          
      polygon(x=c(cc-1,cc-1,cc,cc),y=c(2,3,3,2),col=col1,border=col1)#FO
      polygon(x=c(cc-1,cc-1,cc,cc),y=c(1,2,2,1),col=col2,border=col2)#BA
      polygon(x=c(cc-1,cc-1,cc,cc),y=c(0,1,1,0),col=col3,border=col3)#FR
      
      #text(cc-0.5,2.5,labels=round(HabEff[2],1),col=ifelse(col2!='#FFFFFF00','white',col2),cex=6/12,font=2)
      #text(cc-0.5,1.5,labels=round(HabEff[1],1),col=ifelse(col1!='#FFFFFF00','white',col1),cex=6/12,font=2)
      #text(cc-0.5,0.5,labels=round(HabEff[3],1),col=ifelse(col3!='#FFFFFF00','white',col3),cex=6/12,font=2)
      text(cc-0.5,2.5,labels=round(TrophicGroupListI[[ff]][1,cc],1),col=ifelse(col1=='#FFFFFF00','grey40','white'),cex=6/12,font=2)
      text(cc-0.5,1.5,labels=round(TrophicGroupListI[[ff]][2,cc],1),col=ifelse(col2=='#FFFFFF00','grey40','white'),cex=6/12,font=2)
      text(cc-0.5,0.5,labels=round(TrophicGroupListI[[ff]][3,cc],1),col=ifelse(col3=='#FFFFFF00','grey40','white'),cex=6/12,font=2)
  }
  segments(x0=0,x1=ccmax,y0=0,y1=0,lwd=gridLWD)
  segments(x0=0,x1=ccmax,y0=1,y1=1,lwd=gridLWD)
  segments(x0=0,x1=ccmax,y0=2,y1=2,lwd=gridLWD)
  segments(x0=0,x1=ccmax,y0=3,y1=3,lwd=gridLWD)
  
  segments(x0=0,x1=0,y0=0,y1=3,lwd=gridLWD)
  segments(x0=1,x1=1,y0=0,y1=3,lwd=gridLWD)
  segments(x0=2,x1=2,y0=0,y1=3,lwd=gridLWD)
  segments(x0=3,x1=3,y0=0,y1=3,lwd=gridLWD)
  segments(x0=4,x1=4,y0=0,y1=3,lwd=gridLWD)
  if(ff==1) segments(x0=5,x1=5,y0=0,y1=3,lwd=gridLWD)
}
  
#mtext('Effect on coral',3)


#####MACROMODELS

for(ff in 1:7){
  
  par(new=T)
  
  par(fig=c(0.5,1,(7-ff)/7,(7-(ff-1))/7))
  if(ff==1){
    plot(0,type='n',xlim=c(0,5),ylim=c(0,3),ann=F,axes=F,yaxs='i',xaxs='i')
  }else{
    plot(0,type='n',xlim=c(0,5),ylim=c(0,3),ann=F,axes=F,yaxs='i',xaxs='i')
  }
  
  ccmax <- ifelse(ff==1,5,4)
  for(cc in 1:ccmax){
    ThisColor<-Cols[cc]
    
    col1<-ifelse(!is.na(Macro_TrophicGroupListI_SIG[[ff]][1,cc]),paste0(ThisColor),'#FFFFFF00')
    col2<-ifelse(!is.na(Macro_TrophicGroupListI_SIG[[ff]][2,cc]),paste0(ThisColor),'#FFFFFF00')
    col3<-ifelse(!is.na(Macro_TrophicGroupListI_SIG[[ff]][3,cc]),paste0(ThisColor),'#FFFFFF00')
    
    polygon(x=c(cc-1,cc-1,cc,cc),y=c(2,3,3,2),col=col1,border=col1)#FO
    polygon(x=c(cc-1,cc-1,cc,cc),y=c(1,2,2,1),col=col2,border=col2)#BA
    polygon(x=c(cc-1,cc-1,cc,cc),y=c(0,1,1,0),col=col3,border=col3)#FR
    
    #text(cc-0.5,2.5,labels=round(HabEff[2],1),col=ifelse(col2!='#FFFFFF00','white',col2),cex=6/12,font=2)
    #text(cc-0.5,1.5,labels=round(HabEff[1],1),col=ifelse(col1!='#FFFFFF00','white',col1),cex=6/12,font=2)
    #text(cc-0.5,0.5,labels=round(HabEff[3],1),col=ifelse(col3!='#FFFFFF00','white',col3),cex=6/12,font=2)
    text(cc-0.5,2.5,labels=round(Macro_TrophicGroupListI[[ff]][1,cc],1),col=ifelse(col1=='#FFFFFF00','grey40','white'),cex=6/12,font=2)
    text(cc-0.5,1.5,labels=round(Macro_TrophicGroupListI[[ff]][2,cc],1),col=ifelse(col2=='#FFFFFF00','grey40','white'),cex=6/12,font=2)
    text(cc-0.5,0.5,labels=round(Macro_TrophicGroupListI[[ff]][3,cc],1),col=ifelse(col3=='#FFFFFF00','grey40','white'),cex=6/12,font=2)
  }
  segments(x0=0,x1=ccmax,y0=0,y1=0,lwd=gridLWD)
  segments(x0=0,x1=ccmax,y0=1,y1=1,lwd=gridLWD)
  segments(x0=0,x1=ccmax,y0=2,y1=2,lwd=gridLWD)
  segments(x0=0,x1=ccmax,y0=3,y1=3,lwd=gridLWD)
  
  segments(x0=0,x1=0,y0=0,y1=3,lwd=gridLWD)
  segments(x0=1,x1=1,y0=0,y1=3,lwd=gridLWD)
  segments(x0=2,x1=2,y0=0,y1=3,lwd=gridLWD)
  segments(x0=3,x1=3,y0=0,y1=3,lwd=gridLWD)
  segments(x0=4,x1=4,y0=0,y1=3,lwd=gridLWD)
  if(ff==1) segments(x0=5,x1=5,y0=0,y1=3,lwd=gridLWD)
}
