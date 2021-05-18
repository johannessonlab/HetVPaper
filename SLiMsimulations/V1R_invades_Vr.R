# The goal of this code is to run slim simulations with variable selfing rates, variable prezygotic isolation parameters and variable degree of balancing selection between V and V1 to see how those affect the composition of the population.
#The basic scheme of a population run is to start with Vr and full outcross, introduce a small fraction of V1R individual, then sample genotype frequencies.


#-----------------------------------------------------
#Inspired from Ben Haller's code, from Slim-Extra's repository (cite:https://github.com/MesserLab/SLiM-Extras)
#IMPORTANT: update the next two values for the path to the slim software (slim_path) and slim script (script_path) for your own computer.
slim_path<-"/usr/local/bin/slim"
script_path<-"/Users/ivian/Desktop/dikaryotic_V_Reversed_story.slim"


#-----------------------------------------------------
#I need to pass values for the following parameters to the slim code:
#selfing rates
#the level of prezygotic isolation due to V/V1 interactions (from 0 to 0.6, 0% to 60%)
#the intensity of the negative frequency dependence between V and V1, wmax2, from 0 to 1.

#this is done with the following function by Ben Haller
doOneRun<-function(x)
{
	rep<-x[1]
	self1<-x[2]
	self2<-x[3]
	prezygV<-x[4]
	wmax2<-x[5]
	samplerate<-x[6]
	script<-x[7]
	system2(slim_path,args=c("-d",paste0("rep=",rep),"-d",paste0("selfing1=",self1),"-d",paste0("selfing2=",self2),"-d",paste0("prezygV=",prezygV),"-d",paste0("wmax2=",wmax2),"-d",paste0("samplerate=",samplerate),"-s",rep,shQuote(script)),stdout=F,stderr=F)
}


#-----------------------------------------------------
#The next function Loops the slim runs over the parameter combinations.

for (wm2 in c(0,0.2,0.5,1)) #different intensities of balancing selection between V and V1
{for (prez in c(0,0.2,0.4)) #different levels of prezygotic isolation
{for (s2 in c(0,0.25,0.5,0.75,0.90)) #different selfing rates
{for (i in 0:49) #number of replicated runs
{
	doOneRun(c(i,0,s2,prez,wm2,20,script_path))
	}}}}
	

#-----------------------------------------------------
#Data formatting (data that has been created by the SLiM code)
#IMPORTANT:modify the "stwd" line to the directory defined in the SLiM code to store outputs
rep<-50 #number of replicated simulations (see above in the loop that launches slim, variable i)
totgen<-1500. #total number of generation...see slim code
samplerate<-20 #genotype frequency sampling rate, see slim code
Arun<-totgen/samplerate
setwd("/Users/ivian/Desktop/SlimOuts") # directory where outputs from the SLiM code have been stored
datR<-read.table("output_RUNS_V1R_inv_Vr.txt")
colnames(datR)<-c("batch","balancing","prezyg","selfing","rep","gen","V1R","Vr","V1r","VR")
str(datR)
datRA<-subset(datR,datR$batch=="A")  # A and B correspond to two simulations with same parameters run in parallel, see slim code
summary(datRA)
datRB<-subset(datR,datR$batch=="B")
summary(datRB)
    
 
#--------------------------------------
#analysis of selfing phase generations 900-1000
selfdataRB<-data.frame(batch=double(),balancing=double(),prezyg=double(),selfing=double(),V1R=double(),Vr=double(),V1r=double(),VR=double())
str(selfdataRB)

for (bal in c(0,0.2,0.5,1)) #different intensities of balancing selection between V and V1
{for (pre in c(0,0.2,0.4)) #different levels of prezygotic isolation
{for (sel in c(0,0.25,0.5,0.75,0.90)) #different selfing rates
{for (i in 0:49) #number of replicated runs
{	
sub<-datRB[which(datRB$balancing==bal & datRB$prezyg==pre & datRB$selfing==sel ),]
selfdataRB<-rbind(selfdataRB,data.frame(batch="B",balancing=bal,prezyg=pre,selfing=sel,rep=i,V1R=mean(sub[(45+(Arun*i)):(50+(Arun*i)),7]),Vr=mean(sub[(45+(Arun*i)):(50+(Arun*i)),8]),V1r=mean(sub[(45+(Arun*i)):(50+(Arun*i)),9]),VR=mean(sub[(45+(Arun*i)):(50+(Arun*i)),10])))
}}}}
summary(selfdataRB)
selfdataRB$balancing<-as.factor(selfdataRB$balancing)



selfdataRA<-data.frame(batch=double(),balancing=double(),prezyg=double(),selfing=double(),V1R=double(),Vr=double(),V1r=double(),VR=double())
str(selfdataRA)

for (bal in c(0,0.2,0.5,1)) #different intensities of balancing selection between V and V1
{for (pre in c(0,0.2,0.4)) #different levels of prezygotic isolation
{for (sel in c(0,0.25,0.5,0.75,0.90)) #different selfing rates
{for (i in 0:49) #number of replicated runs
{	
sub<-datRA[which(datRA$balancing==bal & datRA$prezyg==pre & datRA$selfing==sel ),]
selfdataRA<-rbind(selfdataRA,data.frame(batch="A",balancing=bal,prezyg=pre,selfing=sel,rep=i,V1R=mean(sub[(45+(Arun*i)):(50+(Arun*i)),7]),Vr=mean(sub[(45+(Arun*i)):(50+(Arun*i)),8]),V1r=mean(sub[(45+(Arun*i)):(50+(Arun*i)),9]),VR=mean(sub[(45+(Arun*i)):(50+(Arun*i)),10])))
}}}}
summary(selfdataRA)
selfdataRA$balancing<-as.factor(selfdataRA$balancing)

#Full data A+B
selfdataR<-data.frame(batch=double(),balancing=double(),prezyg=double(),selfing=double(),V1R=double(),Vr=double(),V1r=double(),VR=double())
str(selfdataR)
selfdataR<-rbind(selfdataRA,selfdataRB)
selfdataR$balancing<-as.factor(selfdataR$balancing)

selfdata000<-subset(selfdataR,selfdataR$self==0)
selfdata025<-subset(selfdataR,selfdataR$self==0.25)
selfdata050<-subset(selfdataR,selfdataR$self==0.5)
selfdata075<-subset(selfdataR,selfdataR$self==0.75)
selfdata090<-subset(selfdataR,selfdataR$self==0.90)
selfdata099<-subset(selfdataR,selfdataR$self==0.99)
 
library(ggplot2)
library(cowplot)


	#-----Color fill is genotype and contour is blancing selection
g1<-ggplot(selfdata000, aes(x=prezyg, y=V1R,color=balancing,group=interaction(prezyg,balancing))) +
      geom_boxplot(fill="darkseagreen3",width=0.05, lwd=0.5,position=position_dodge(0.13),outlier.colour="darkseagreen3", outlier.shape=16,
             outlier.size=1)+
      scale_color_manual(name="Strength of balancing selection",values=c("orange", "chocolate3", "red3","red4"))+
      coord_cartesian(ylim=c(0,1))+
      xlab("Level of prezygotic isolation")+
      ylab("Haplotype frequency")+
      ggtitle("Selfing rate 00%")+
       theme_minimal_hgrid(12)
g1<-g1+geom_boxplot( aes(x=prezyg+0.005,y=Vr,color=balancing,group=interaction(prezyg,balancing)), fill= "#B66DE9",width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="#B66DE9", outlier.shape=16,
             outlier.size=1)
g1<-g1+geom_boxplot( aes(x=prezyg-0.001,y=V1r,color=balancing,group=interaction(prezyg,balancing)), fill= "gray", alpha=0.85,width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="gray", outlier.shape=16,
             outlier.size=1)

g2<-ggplot(selfdata025, aes(x=prezyg, y=V1R,color=balancing,group=interaction(prezyg,balancing))) +
      geom_boxplot(fill="darkseagreen3",width=0.05, lwd=0.5,position=position_dodge(0.13),outlier.colour="darkseagreen3", outlier.shape=16,
             outlier.size=1)+
      scale_color_manual(name="Strength of balancing selection",values=c("orange", "chocolate3", "red3","red4"))+
      coord_cartesian(ylim=c(0,1))+
      xlab("Level of prezygotic isolation")+
      ylab("Haplotype frequency")+
      ggtitle("Selfing rate 25%")+
       theme_minimal_hgrid(12)
g2<-g2+geom_boxplot( aes(x=prezyg+0.005,y=Vr,color=balancing,group=interaction(prezyg,balancing)), fill= "#B66DE9",width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="#B66DE9", outlier.shape=16,
             outlier.size=1)
g2<-g2+geom_boxplot( aes(x=prezyg-0.001,y=V1r,color=balancing,group=interaction(prezyg,balancing)), fill= "gray", alpha=0.85,width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="gray", outlier.shape=16,
             outlier.size=1)

g3<-ggplot(selfdata050, aes(x=prezyg, y=V1R,color=balancing,group=interaction(prezyg,balancing))) +
      geom_boxplot(fill="darkseagreen3",width=0.05, lwd=0.5,position=position_dodge(0.13),outlier.colour="darkseagreen3", outlier.shape=16,
             outlier.size=1)+
      scale_color_manual(name="Strength of balancing selection",values=c("orange", "chocolate3", "red3","red4"))+
      coord_cartesian(ylim=c(0,1))+
      xlab("Level of prezygotic isolation")+
      ylab("Haplotype frequency")+
      ggtitle("Selfing rate 50%")+
       theme_minimal_hgrid(12)
g3<-g3+geom_boxplot( aes(x=prezyg+0.005,y=Vr,color=balancing,group=interaction(prezyg,balancing)), fill= "#B66DE9",width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="#B66DE9", outlier.shape=16,
             outlier.size=1)
g3<-g3+geom_boxplot( aes(x=prezyg-0.001,y=V1r,color=balancing,group=interaction(prezyg,balancing)), fill= "gray", alpha=0.85,width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="gray", outlier.shape=16,
             outlier.size=1)

g4<-ggplot(selfdata075, aes(x=prezyg, y=V1R,color=balancing,group=interaction(prezyg,balancing))) +
      geom_boxplot(fill="darkseagreen3",width=0.05, lwd=0.5,position=position_dodge(0.13),outlier.colour="darkseagreen3", outlier.shape=16,
             outlier.size=1)+
      scale_color_manual(name="Strength of balancing selection",values=c("orange", "chocolate3", "red3","red4"))+
      coord_cartesian(ylim=c(0,1))+
      xlab("Level of prezygotic isolation")+
      ylab("Haplotype frequency")+
      ggtitle("Selfing rate 75%")+
      theme_minimal_hgrid(12)
g4<-g4+geom_boxplot( aes(x=prezyg+0.005,y=Vr,color=balancing,group=interaction(prezyg,balancing)), fill= "#B66DE9",width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="#B66DE9", outlier.shape=16,
             outlier.size=1)
g4<-g4+geom_boxplot( aes(x=prezyg-0.001,y=V1r,color=balancing,group=interaction(prezyg,balancing)), fill= "gray", alpha=0.85,width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="gray", outlier.shape=16,
             outlier.size=1)


g5<-ggplot(selfdata090, aes(x=prezyg, y=V1R,color=balancing,group=interaction(prezyg,balancing))) +
      geom_boxplot(fill="darkseagreen3",width=0.05, lwd=0.2,position=position_dodge(0.13),outlier.colour="darkseagreen3", outlier.shape=16,
             outlier.size=1)+
      scale_color_manual(name="Strength of balancing selection",values=c("orange", "chocolate3", "red3","red4"))+
      coord_cartesian(ylim=c(0,1))+
      xlab("Level of prezygotic isolation")+
      ylab("Haplotype frequency")+
      ggtitle("Selfing rate 90%")+
      theme_minimal_hgrid(12)
g5<-g5+geom_boxplot( aes(x=prezyg+0.005,y=Vr,color=balancing,group=interaction(prezyg,balancing)), fill= "#B66DE9",width=0.05,lwd=0.2,position=position_dodge(0.13),outlier.colour="#B66DE9", outlier.shape=16,
             outlier.size=1)
g5<-g5+geom_boxplot( aes(x=prezyg-0.001,y=V1r,color=balancing,group=interaction(prezyg,balancing)), fill= "gray", alpha=0.85,width=0.05,lwd=0.2,position=position_dodge(0.13),outlier.colour="gray", outlier.shape=16,
             outlier.size=1)


g6<-ggplot(selfdata099, aes(x=prezyg, y=V1R,color=balancing,group=interaction(prezyg,balancing))) +
      geom_boxplot(fill="darkseagreen3",width=0.05, lwd=0.2,position=position_dodge(0.13),outlier.colour="darkseagreen3", outlier.shape=16,
             outlier.size=1)+
      scale_color_manual(name="Strength of balancing selection",values=c("orange", "chocolate3", "red3","red4"))+
      coord_cartesian(ylim=c(0,1))+
      xlab("Level of prezygotic isolation")+
      ylab("Haplotype frequency")+
      ggtitle("Selfing rate 99%")+
    theme_minimal_hgrid(12)
g6<-g6+geom_boxplot( aes(x=prezyg+0.005,y=Vr,color=balancing,group=interaction(prezyg,balancing)), fill= "#B66DE9",width=0.05,lwd=0.2,position=position_dodge(0.13),outlier.colour="#B66DE9", outlier.shape=16,
             outlier.size=1)
g6<-g6+geom_boxplot( aes(x=prezyg-0.001,y=V1r,color=balancing,group=interaction(prezyg,balancing)), fill= "gray", alpha=0.85,width=0.05,lwd=0.2,position=position_dodge(0.13),outlier.colour="gray", outlier.shape=16,
             outlier.size=1)


plot_grid(g1+theme (legend.position="none"),g2+theme (legend.position="none"),g3+theme (legend.position="bottom"),align="vh",ncol=1)
plot_grid(g4+theme (legend.position="none"),g5+theme (legend.position="none"),g6+theme (legend.position="bottom"),align="vh",ncol=1)

#--------------------------------------
#analysis of the after selfing is LOST phase. generations 1300-1400, 300 gneration after selfing is lost.
latedata<-data.frame(balancing=double(),prezyg=double(),selfing=double(),V1R=double(),Vr=double(),V1r=double(),VR=double())
str(latedata)

for (bal in c(0,0.2,0.5,1)) #different intensities of balancing selection between V and V1
{for (pre in c(0,0.2,0.4,0.6)) #different levels of prezygotic isolation
{for (sel in c(0,0.25,0.5,0.75,0.90,0.99)) #different selfing rates
{for (i in 0:99) #number of replicated runs
{	
sub<-datA[which(datA$balancing==bal & datA$prezyg==pre & datA$selfing==sel ),]
latedata<-rbind(latedata,data.frame(balancing=bal,prezyg=pre,selfing=sel,rep=i,V1R=mean(sub[(65+(Arun*i)):(70+(Arun*i)),7]),Vr=mean(sub[(65+(Arun*i)):(70+(Arun*i)),8]),V1r=mean(sub[(65+(Arun*i)):(70+(Arun*i)),9]),VR=mean(sub[(65+(Arun*i)):(70+(Arun*i)),10])))
}}}}
summary(latedata)
latedata$balancing<-as.factor(latedata$balancing)

latedata000<-subset(latedata,latedata$self==0)
latedata025<-subset(latedata,latedata$self==0.25)
latedata050<-subset(latedata,latedata$self==0.5)
latedata075<-subset(latedata,latedata$self==0.75)
latedata090<-subset(latedata,latedata$self==0.90)
latedata099<-subset(latedata,latedata$self==0.99)

lg1<-ggplot(latedata000, aes(x=prezyg, y=V1R,color=balancing,group=interaction(prezyg,balancing))) +
      geom_boxplot(fill="darkseagreen3",width=0.05, lwd=0.5,position=position_dodge(0.13),outlier.colour="darkseagreen3", outlier.shape=16,
             outlier.size=1)+
      scale_color_manual(name="Strength of balancing selection",values=c("orange", "chocolate3", "red3","red4"))+
      coord_cartesian(ylim=c(0,1))+
      xlab("Level of prezygotic isolation")+
      ylab("Haplotype frequency")+
      ggtitle("Selfing rate 00%- After 00%")+
      theme(legend.position="bottom",plot.title=element_text(hjust=0.5))
lg1<-lg1+geom_boxplot( aes(x=prezyg+0.005,y=Vr,color=balancing,group=interaction(prezyg,balancing)), fill= "#B66DE9",width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="#B66DE9", outlier.shape=16,
             outlier.size=1)
lg1<-lg1+geom_boxplot( aes(x=prezyg-0.001,y=V1r,color=balancing,group=interaction(prezyg,balancing)), fill= "gray", alpha=0.85,width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="gray", outlier.shape=16,
             outlier.size=1)
plot(lg1)

lg2<-ggplot(latedata025, aes(x=prezyg, y=V1R,color=balancing,group=interaction(prezyg,balancing))) +
      geom_boxplot(fill="darkseagreen3",width=0.05, lwd=0.5,position=position_dodge(0.13),outlier.colour="darkseagreen3", outlier.shape=16,
             outlier.size=1)+
      scale_color_manual(name="Strength of balancing selection",values=c("orange", "chocolate3", "red3","red4"))+
      coord_cartesian(ylim=c(0,1))+
      xlab("Level of prezygotic isolation")+
      ylab("Haplotype frequency")+
      ggtitle("Selfing rate 00%- After 25%")+
      theme(legend.position="bottom",plot.title=element_text(hjust=0.5))
lg2<-lg2+geom_boxplot( aes(x=prezyg+0.005,y=Vr,color=balancing,group=interaction(prezyg,balancing)), fill= "#B66DE9",width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="#B66DE9", outlier.shape=16,
             outlier.size=1)
lg2<-lg2+geom_boxplot( aes(x=prezyg-0.001,y=V1r,color=balancing,group=interaction(prezyg,balancing)), fill= "gray", alpha=0.85,width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="gray", outlier.shape=16,
             outlier.size=1)
plot(lg2)

lg3<-ggplot(latedata050, aes(x=prezyg, y=V1R,color=balancing,group=interaction(prezyg,balancing))) +
      geom_boxplot(fill="darkseagreen3",width=0.05, lwd=0.5,position=position_dodge(0.13),outlier.colour="darkseagreen3", outlier.shape=16,
             outlier.size=1)+
      scale_color_manual(name="Strength of balancing selection",values=c("orange", "chocolate3", "red3","red4"))+
      coord_cartesian(ylim=c(0,1))+
      xlab("Level of prezygotic isolation")+
      ylab("Haplotype frequency")+
      ggtitle("Selfing rate 00%- After 50%")+
      theme(legend.position="bottom",plot.title=element_text(hjust=0.5))
lg3<-lg3+geom_boxplot( aes(x=prezyg+0.005,y=Vr,color=balancing,group=interaction(prezyg,balancing)), fill= "#B66DE9",width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="#B66DE9", outlier.shape=16,
             outlier.size=1)
lg3<-lg3+geom_boxplot( aes(x=prezyg-0.001,y=V1r,color=balancing,group=interaction(prezyg,balancing)), fill= "gray", alpha=0.85,width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="gray", outlier.shape=16,
             outlier.size=1)
plot(lg3)

lg4<-ggplot(latedata075, aes(x=prezyg, y=V1R,color=balancing,group=interaction(prezyg,balancing))) +
      geom_boxplot(fill="darkseagreen3",width=0.05, lwd=0.5,position=position_dodge(0.13),outlier.colour="darkseagreen3", outlier.shape=16,
             outlier.size=1)+
      scale_color_manual(name="Strength of balancing selection",values=c("orange", "chocolate3", "red3","red4"))+
      coord_cartesian(ylim=c(0,1))+
      xlab("Level of prezygotic isolation")+
      ylab("Haplotype frequency")+
      ggtitle("Selfing rate 00%- After 75%")+
      theme(legend.position="bottom",plot.title=element_text(hjust=0.5))
lg4<-lg4+geom_boxplot( aes(x=prezyg+0.005,y=Vr,color=balancing,group=interaction(prezyg,balancing)), fill= "#B66DE9",width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="#B66DE9", outlier.shape=16,
             outlier.size=1)
lg4<-lg4+geom_boxplot( aes(x=prezyg-0.001,y=V1r,color=balancing,group=interaction(prezyg,balancing)), fill= "gray", alpha=0.85,width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="gray", outlier.shape=16,
             outlier.size=1)
plot(lg4)

lg5<-ggplot(latedata090, aes(x=prezyg, y=V1R,color=balancing,group=interaction(prezyg,balancing))) +
      geom_boxplot(fill="darkseagreen3",width=0.05, lwd=0.5,position=position_dodge(0.13),outlier.colour="darkseagreen3", outlier.shape=16,
             outlier.size=1)+
      scale_color_manual(name="Strength of balancing selection",values=c("orange", "chocolate3", "red3","red4"))+
      coord_cartesian(ylim=c(0,1))+
      xlab("Level of prezygotic isolation")+
      ylab("Haplotype frequency")+
      ggtitle("Selfing rate 00%- After 90%")+
      theme(legend.position="bottom",plot.title=element_text(hjust=0.5))
lg5<-lg5+geom_boxplot( aes(x=prezyg+0.005,y=Vr,color=balancing,group=interaction(prezyg,balancing)), fill= "#B66DE9",width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="#B66DE9", outlier.shape=16,
             outlier.size=1)
lg5<-g5+geom_boxplot( aes(x=prezyg-0.001,y=V1r,color=balancing,group=interaction(prezyg,balancing)), fill= "gray", alpha=0.85,width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="gray", outlier.shape=16,
             outlier.size=1)
plot(lg5)

g6<-ggplot(selfdata099, aes(x=prezyg, y=V1R,color=balancing,group=interaction(prezyg,balancing))) +
      geom_boxplot(fill="darkseagreen3",width=0.05, lwd=0.5,position=position_dodge(0.13),outlier.colour="darkseagreen3", outlier.shape=16,
             outlier.size=1)+
      scale_color_manual(name="Strength of balancing selection",values=c("orange", "chocolate3", "red3","red4"))+
      coord_cartesian(ylim=c(0,1))+
      xlab("Level of prezygotic isolation")+
      ylab("Haplotype frequency")+
      ggtitle("Selfing rate 99%")+
      theme(legend.position="bottom",plot.title=element_text(hjust=0.5))
g6<-g6+geom_boxplot( aes(x=prezyg+0.005,y=Vr,color=balancing,group=interaction(prezyg,balancing)), fill= "#B66DE9",width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="#B66DE9", outlier.shape=16,
             outlier.size=1)
g6<-g6+geom_boxplot( aes(x=prezyg-0.001,y=V1r,color=balancing,group=interaction(prezyg,balancing)), fill= "gray", alpha=0.85,width=0.05,lwd=0.5,position=position_dodge(0.13),outlier.colour="gray", outlier.shape=16,
             outlier.size=1)
plot(g6)