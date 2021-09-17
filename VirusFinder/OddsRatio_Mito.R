if(!require(rcompanion)){install.packages("rcompanion")}
if(!require(dplyr)){install.packages("dplyr")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(grid)){install.packages("grid")}
if(!require(pwr)){install.packages("pwr")}
library(readxl)

##import file
library(readxl)
Statistics_Rnd1_and_Rnd2_V2_0 <- read_excel("C:/Users/savan/OneDrive/Cafgen/Analysis/Viral_Reads/Statistics_Rnd1_and_Rnd2_V2.0.xlsx", 
                                              +     sheet = "UniqueMaleStats")
view(Statistics_Rnd1_and_Rnd2_V2_0)

ODDS_RATIO2 <- data.frame(Mito_CAFGEN_haplogroups, row.names = "Haplogroup")


row_fisher <- function(row, alt = 'two.sided', cnf = 0.95) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           p_val = f$p.value,
           or = f$estimate[[1]],
           or_ll = f$conf.int[1],
           or_ul = f$conf.int[2]))
}
#test_df <- data.frame(matrix(sample.int(1000, size = 4000, replace = TRUE), ncol = 4))
#colnames(test_df) <- c('a', 'b', 'c', 'd')
p <- data.frame(t(apply(ODDS_RATIO2, 1, row_fisher)))
print(p)

write.csv(p, file = "FemaleOddsRatio.csv")


# fisher test on each row
row_fisher <- function(row, alt = 'two.sided', cnf = 0.95) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           p_val = f$p.value,
           or = f$estimate[[1]],
           or_ll = f$conf.int[1],
           or_ul = f$conf.int[2]))
}

# generate sample data for test
df_test <- tbl_df(matrix(sample.int(1000, size = 4000, replace = TRUE), ncol = 4)) %>%
  rename(a = V1, b = V2, c = V3, d = V4)

# run
p <- df_test %>%
  parApply(cl, ., 1, row_fisher) %>%
  t() %>%
  tbl_df()
print(p)


##Calculating fishers exact test
Uganda <- c(93,69,298,352)
Uganda1 <- matrix(Uganda, nrow = 2, ncol = 2, byrow = TRUE)
fisher.test(Uganda1, conf.int = TRUE, or = 1)
