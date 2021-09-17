library(ggplot2)

library(readxl)
ODDS_RATIO_R_plot <- read_excel("C:/Users/savan/OneDrive/Cafgen/Analysis/Viral_Reads/ODDS_RATIO_R_plot.xlsx", sheet = "Sheet2")

RR_data <- ODDS_RATIO_R_plot
#RR_data$UpperLimit[RR_data$UpperLimit > 2] = 2
#RR_data$LowerLimit[RR_data$LowerLimit < 0.5] = 0.5
p = ggplot(data=RR_data,
           aes(x = Group,y = OddsRatio, ymin = LowerLimit, ymax = UpperLimit ))+
  geom_pointrange(aes(col=Group))+
  geom_hline(aes(fill=Group),yintercept =1, linetype=2)+
  xlab('Virus Family')+ ylab("Odds Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(ymin=LowerLimit, ymax=UpperLimit,col=Group),width=0.5,cex=1)+ 
  facet_wrap(~Virus_Family,strip.position="left",nrow=9,scales = "free_y") +
  theme_minimal() +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip()
p
##LogOdds
ODDS_RATIO_R_plot <- read_excel("C:/Users/savan/OneDrive/Cafgen/Analysis/Viral_Reads/ODDS_RATIO_R_plot.xlsx", 
                                  +     sheet = "The rest (2)")
RR_data <- ODDS_RATIO_R_plot
positions <- c("HiSeq_2500", "Novaseq_6000", "Botswana", "Uganda", "Female", "Male", "Combined")
positions_family <- c("Anelloviridae", "Microviridae", "Adenoviridae", "Baculoviridae", "Herpesviridae", "Paramyxoviridae", "Poxviridae", "Retroviridae", "Siphoviridae", "Togaviridae")
RR_data$Group <- factor(RR_data$Group, levels = positions)
RR_data$Virus_Family <- factor(RR_data$Virus_Family, levels = positions_family)
#RR_data$UpperLimit[RR_data$UpperLimit > 2] = 2
#RR_data$LowerLimit[RR_data$LowerLimit < 0.5] = 0.5
library(ggplot2)
p <- ggplot(data=RR_data,
           aes(x = Group,y = LogOddsRatio, ymin = LogLowerLimit, ymax = LogUpperLimit ))+
  geom_pointrange(aes(col=Group))+
  geom_hline(aes(fill=Group),yintercept =0, linetype=2)+
  xlab('Virus Family')+ ylab("Log Odds Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(ymin=LogLowerLimit, ymax=LogUpperLimit,col=Group),width=0.5,cex=1)+ 
  facet_wrap(~Virus_Family,strip.position="top", nrow=5) +
  #scale_color_manual(values=c("darkslateblue", "chartreuse4", "firebrick1")) +
  #scale_color_brewer(palette="Dark2") +
  theme_minimal() +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        #strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"),
        plot.margin = margin(10, 10, 10, 10, "pt")) +
  coord_flip() +
  guides(color = guide_legend(reverse = TRUE))
p
###export pdf dimensions 5.6 x 10, portrait
install.packages('svglite')
library(svglite)
ggsave(file="LogOddsRatioPlot.svg", plot=p, width=5.8, height=5.6)