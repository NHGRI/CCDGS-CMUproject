library(readxl)
CombinedVirusList_Rnd1_and_Rnd2_V2_0 <- read_excel("C:/Users/savan/OneDrive/Cafgen/Analysis/Viral_Reads/CombinedVirusList_Rnd1_and_Rnd2_V2.0.xlsx", 
                                                     +     sheet = "WithPhenoData (2)")
View(CombinedVirusList_Rnd1_and_Rnd2_V2_0)
virus_summary <- CombinedVirusList_Rnd1_and_Rnd2_V2_0
p <- ggplot(virus_summary, aes(x=ICTV_NAME, y=Identities)) +
  #stat_boxplot(aes(x=ICTV_NAME, y=Identities), geom='errorbar', width = 0.5, color="#FE420F", size = 1) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.5, color="#FE420F", size = 1) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.5, color="#FE420F", size = 1) +
  geom_boxplot(aes(x=ICTV_NAME, y=Identities), size = .75, outlier.shape = NA, color="#FE420F", fill="#FE420F", alpha=0.6) +
  geom_jitter(shape=16, size=1, color ="#2A2725", position=position_jitter(0.2)) +
  labs(title="Boxplot of Percent Identities",x="Virus Family", y = "Percent Identity") +
  theme_minimal() +
  #scale_color_brewer(palette="Dark2") +
  theme(plot.title=element_text(size=16, face="bold", hjust=0.5),
        axis.text.y=element_text(face="bold", hjust = 1),
        axis.text.x=element_text(face="bold", angle = -45, hjust = 0),
        axis.title=element_text(size=12,face="bold"),
        plot.margin = margin(10, 10, 0, 10, "pt"))
p


# Rotate the box plot
#p + coord_flip()
# Notched box plot
#ggplot(ToothGrowth, aes(x=dose, y=len)) + 
  #geom_boxplot(notch=TRUE)
# Change outlier, color, shape and size
