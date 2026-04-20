# this is the script for plotting PCA results from WES data

setwd("/Users/hany4/Documents/Projects/Hypertension/COEH")
workdir = getwd()

list.files("./PCA", all.files = T)

# read in the PED data
ped1kg = read.table("./PCA/integrated_call_samples_v3.20130502.ALL.panel", 
                    header = TRUE, fill = T, skip = 0, sep = "\t")
ped1kg = ped1kg[, c(1:5)]
names(ped1kg) = c("IID", "pop", "super_pop", "gender", "Relationship" )
head(ped1kg)
dim(ped1kg)

pedcoeh = read.table("./PCA/COEH_202206.ped", header = TRUE, skip = 0, sep = "\t")
pedcoeh$super_pop = c("COEH")
head(pedcoeh)
dim(pedcoeh)

# rbind two ped files
ped = rbind(ped1kg, pedcoeh)
head(ped)
dim(ped)

# read in MDS data
# mds = read.delim("./PCA/MergedCOEH.Merged1KG.matched.merged.mds", header = T, sep = "")
mds = read.delim("./PCA/MergedCOEH_maf20.Merged1KG_maf20.matched.merged.mds", header = T, sep = "")
head(mds)
dim(mds)

# merge MDS and PED
mds.ped = merge(mds, ped, by = "IID", all = T)
head(mds.ped)
dim(mds.ped)
table(mds.ped$super_pop)
table(mds.ped$pop)
class(mds.ped$super_pop)
mds.ped$super_pop = factor(mds.ped$super_pop, levels = c("AFR", "AMR", "EAS", "EUR", "SAS", "COEH"))

mds.ped[(mds.ped$IID == "PH10190"),]$pop = c("PH")
mds.ped[(mds.ped$IID == "PH10190"),]$super_pop = c("COEH")
mds.ped[(mds.ped$IID == "PH10243"),]$pop = c("PH")
mds.ped[(mds.ped$IID == "PH10243"),]$super_pop = c("COEH")

#### plot merged data - grey 1000Genome data
library(ggplot2)
library(ggrepel)

mds.ped$label = c("")
mds.ped[(mds.ped$IID == "PH10433"),]$label = c("Y")
mds.ped[(mds.ped$IID == "PH10435"),]$label = c("Y")
mds.ped[(mds.ped$IID == "PH10149"),]$label = c("Y")
mds.ped[(mds.ped$IID == "PH10237"),]$label = c("Y")
table(mds.ped$label)
table(mds.ped$super_pop)

# mds.ped[(mds.ped$pop != "PH" & mds.ped$super_pop == "AFR"),]$label = c("African")
# mds.ped[(mds.ped$pop != "PH" & mds.ped$super_pop == "AMR"),]$label = c("Admixed American")
# mds.ped[(mds.ped$pop != "PH" & mds.ped$super_pop == "EAS"),]$label = c("East Asian")
# mds.ped[(mds.ped$pop != "PH" & mds.ped$super_pop == "EUR"),]$label = c("European")
# mds.ped[(mds.ped$pop != "PH" & mds.ped$super_pop == "SAS"),]$label = c("South Asian")


ggplot(mds.ped, aes(C1, C4, color = super_pop)) +
  geom_point(shape = 20, size = 2, alpha = 0.8) +
  geom_point(data = mds.ped[grepl("PH", mds.ped$IID),], color = "black", shape = 20, size = 2) +
  geom_point(data = mds.ped[mds.ped$IID == "PH10433",], color = "dark green", shape = 20, size = 4) +
  geom_point(data = mds.ped[mds.ped$IID == "PH10149",], color = "blue", shape = 20, size = 4) +
  geom_point(data = mds.ped[mds.ped$IID == "PH10237",], color = "purple", shape = 20, size = 4) +
  geom_text_repel(aes(label = IID), data = ~ subset(., IID == "PH10433"), color = "dark green", box.padding = 1, nudge_x = 0.02, nudge_y = -0.03) +
  geom_text_repel(aes(label = IID), data = ~ subset(., IID == "PH10149"), color = "blue", box.padding = 1, nudge_y = -0.02) +
  geom_text_repel(aes(label = IID), data = ~ subset(., IID == "PH10237"), color = "purple", box.padding = 1,  nudge_y = -0.02) +
  theme_classic() +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 28, hjust = 0.5),
        legend.position="bottom") +
  xlab("PC1") +
  ylab("PC2")

ggsave("./PCA/wesPCA_maf20_c1c4.png")

