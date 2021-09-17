library(readxl)
#Dataset <- read.table("C:/Users/savan/OneDrive/Cafgen/Analysis/Viral_Reads/VirusList_Rnd1_and_Rnd2_WIDE.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
GLM.1 <- glm(LTNP ~ Anelloviridae + COUNTRY_BWR + Gender_F + 
   Total_Virus_Contig_Count + Batch, family=binomial(logit), data=Dataset)
summary(GLM.1)
exp(coef(GLM.1))  # Exponentiated coefficients ("odds ratios")
##Calculate CI
Confint(GLM.1, level = 0.95)
