data <- ggplot2dotplot
library(ggplot2)
p <- ggplot(data, aes(fill=LTNP_or_RP, x=Family)) + 
  geom_bar(position="dodge")
p