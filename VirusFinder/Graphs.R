p<- ggplot(Virus_SummaryCount, aes(x=LTNP_RP, y=Count)) +
    geom_dotplot(binaxis = 'y', stackdir = 'center', stackratio = 1.5, dotsize = 0.3)
p
# Change colors
p+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
# Add box plots
ggplot(Virus_SummaryCount, aes(x=ICTV_FamilyName, y=Count, fill=LTNP_RP)) +
  geom_boxplot(fill="white")+
  geom_dotplot(binaxis = 'y', stackdir = 'center', stackratio = 1.5, dotsize = 0.8)
##Filter out Microviridae
VirusFilter<- Virus_SummaryCount[Virus_SummaryCount$ICTV_Family_Name!= "Microviridae",]
VirusFilter

Plot<- ggplot(VirusFilter, aes(x=ICTV_Family_Name, y=Count, fill=LTNP_RP)) +
  geom_boxplot(position=position_dodge(0.8))+
  geom_dotplot(binaxis='y', stackdir='center', stackratio = 1.5, dotsize = 0.5, position=position_dodge(0.8))

Plot
Plot2<- ggplot(VirusFilter, aes(x=ICTV_Family_Name, fill=LTNP_RP)) +
  #geom_bar(stat="identity", width=0.8, position=position_dodge(0.8))+
  geom_bar(aes(weight = Count), width=0.8, position=position_dodge(0.8))
#geom_boxplot(position=position_dodge(0.8)) +
#geom_dotplot(binaxis='y', stackdir='center', stackratio = 1.5, dotsize = 0.5, position=position_dodge(0.8))

Plot2
Plot3<- ggplot(VirusFilter, aes(x=ICTV_Family_Name, y=Count, fill=LTNP_RP)) +
  geom_bar(stat="identity", aes(weight = Count), width=0.8, position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', stackratio = 1.5, dotsize = 0.5, position=position_dodge(0.8))

Plot3

Plot+scale_fill_brewer(palette="Pastel2")+xlab('Virus Family Groups')+ylab('Virus Species Count')



##Select Anelloviridae
Anello <- VirusFilter[VirusFilter$ICTV_Family_Name=="Anelloviridae",]
ggplot(Anello, aes(x=ICTV_Family_Name,  y=Count, fill=LTNP_RP)) +
  geom_boxplot(position=position_dodge(0.8))+
  geom_dotplot(binaxis='y', stackdir='center', stackratio = 1.5, dotsize = 0.5, position=position_dodge(0.8))
Plot
Plot
