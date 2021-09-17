library(readxl)
read_excel("C:/Users/savan/OneDrive/Cafgen/Sample_Phenotypes/PhenotypeCombinedBatch1and2.v2.0_WESSamplesOnly.xlsx")

CombinedVirusList_Rnd1_and_Rnd2_V2_0 <- read_excel("C:/Users/savan/OneDrive/Cafgen/Analysis/Viral_Reads/CombinedVirusList_Rnd1_and_Rnd2_V2.0.xlsx", 
                                                   +     sheet = "WithPhenoData (2)")
data <- PhenotypeCombinedBatch1and2_v2_0_WESSamplesOnly 
p <- ggplot(data, aes(x=RP_or_LTNP, y=ART_duration_at_Collection_recalculated)) +
#p + geom_boxplot()
  #stat_boxplot(aes(x= RP_or_LTNP, y=ART_duration_at_Collection_recalculated , geom='errorbar', width = 0.5, color="#FE420F", size = 1) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.5, color="#29BF12", size = 1) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.5, color="#29BF12", size = 1) +
  geom_boxplot(aes(x= RP_or_LTNP, y=ART_duration_at_Collection_recalculated), size = .75, outlier.shape = NA, color="#29BF12", fill="#29BF12", alpha=0.6) +
  geom_jitter(shape=16, size=1, color ="#2A2725", alpha=0.8, position=position_jitter(0.2)) +
  labs(title="Boxplot of Age",x="Virus Family", y = "Percent Identity") +
  #scale_color_brewer(palette="Dark2") +
  theme(plot.title=element_text(size=16, face="bold", hjust=0.5),
        axis.text.y=element_text(face="bold", hjust = 1),
        axis.text.x=element_text(face="bold", angle = -45, hjust = 0),
        axis.title=element_text(size=12,face="bold"),
        plot.margin = margin(10, 10, 0, 10, "pt"))
p


###plot Age
p <- ggplot(data, aes(x=RP_or_LTNP, y=Age_at_sampleCollection)) +
  #p + geom_boxplot()
  #stat_boxplot(aes(x= RP_or_LTNP, y=Age_at_sampleCollection, geom='errorbar', width = 0.5, color="#FE420F", size = 1) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.5, color="#574AE2", size = 1) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.5, color="#574AE2", size = 1) +
  geom_boxplot(aes(x= RP_or_LTNP, y=Age_at_sampleCollection), size = .75, outlier.shape = NA, color="#574AE2", fill="#574AE2", alpha=0.6) +
  geom_jitter(shape=16, size=1, color ="#2A2725", alpha=0.8, position=position_jitter(0.2)) +
  labs(title="Boxplot of Age",x="Virus Family", y = "Percent Identity") +
  #scale_color_brewer(palette="Dark2") +
  theme(plot.title=element_text(size=16, face="bold", hjust=0.5),
        axis.text.y=element_text(face="bold", hjust = 1),
        axis.text.x=element_text(face="bold", angle = -45, hjust = 0),
        axis.title=element_text(size=12,face="bold"),
        plot.margin = margin(10, 10, 0, 10, "pt"))
p