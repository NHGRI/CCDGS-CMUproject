##Country stratification
ByCountry <-
  array(c(11, 178, 1, 260, 17, 185, 7, 154),
        dim = c(2, 2, 2),
        dimnames = list(
          Anelloviridae = c("Virus", "No_Virus"),
          Phenotype = c("LTNP", "RP"),
          Level = c("BWR", "UGR")))
ByCountry
mantelhaen.test(ByCountry, exact = TRUE, alternative = "greater")
#two-sided
mantelhaen.test(ByCountry, exact = TRUE, alternative = "two.sided", conf.int)
##Gender Stratification
ByGender <-
  array(c(13, 150, 5, 224, 15, 213, 3, 190),
        dim = c(2, 2, 2),
        dimnames = list(
          Anelloviridae = c("Virus", "No_Virus"),
          Phenotype = c("LTNP", "RP"),
          Level = c("Male", "Female")))
ByGender
## One sided
mantelhaen.test(ByGender, exact = TRUE, alternative = "greater")
##Two sided
mantelhaen.test(ByGender, exact = TRUE)

##TEmplate form https://stat.ethz.ch/R-manual/R-devel/library/stats/html/mantelhaen.test.html
Rabbits <-
  array(c(0, 0, 6, 5,
          3, 0, 3, 6,
          6, 2, 0, 4,
          5, 6, 1, 0,
          2, 5, 0, 0),
        dim = c(2, 2, 5),
        dimnames = list(
          Delay = c("None", "1.5h"),
          Response = c("Cured", "Died"),
          Penicillin.Level = c("1/8", "1/4", "1/2", "1", "4")))

ByBatch <-
  array(c(15, 126, 4, 169, 13, 237, 4, 245),
        dim = c(2, 2, 2),
        dimnames = list(
          Anelloviridae = c("Virus", "No_Virus"),
          Phenotype = c("LTNP", "RP"),
          Level = c("Round_1", "Round_2")))
ByBatch
mantelhaen.test(ByBatch, exact = TRUE, alternative = "greater")
#two-sided



Rabbits
Rabbits
## Classical Mantel-Haenszel test
mantelhaen.test(Rabbits)
## => p = 0.047, some evidence for higher cure rate of immediate
##               injection
## Exact conditional test
mantelhaen.test(Rabbits, exact = TRUE)
## => p - 0.040
## Exact conditional test for one-sided alternative of a higher
## cure rate for immediate injection
mantelhaen.test(Rabbits, exact = TRUE, alternative = "greater")
## => p = 0.020

## UC Berkeley Student Admissions
mantelhaen.test(UCBAdmissions)
## No evidence for association between admission and gender
## when adjusted for department.  However,
apply(UCBAdmissions, 3, function(x) (x[1,1]*x[2,2])/(x[1,2]*x[2,1]))
## This suggests that the assumption of homogeneous (conditional)
## odds ratios may be violated.  The traditional approach would be
## using the Woolf test for interaction:
woolf <- function(x) {
  x <- x + 1 / 2
  k <- dim(x)[3]
  or <- apply(x, 3, function(x) (x[1,1]*x[2,2])/(x[1,2]*x[2,1]))
  w <-  apply(x, 3, function(x) 1 / sum(1 / x))
  1 - pchisq(sum(w * (log(or) - weighted.mean(log(or), w)) ^ 2), k - 1)
}
woolf(UCBAdmissions)
## => p = 0.003, indicating that there is significant heterogeneity.
## (And hence the Mantel-Haenszel test cannot be used.)

## Agresti (2002), p. 287f and p. 297.
## Job Satisfaction example.
Satisfaction <-
  as.table(array(c(1, 2, 0, 0, 3, 3, 1, 2,
                   11, 17, 8, 4, 2, 3, 5, 2,
                   1, 0, 0, 0, 1, 3, 0, 1,
                   2, 5, 7, 9, 1, 1, 3, 6),
                 dim = c(4, 4, 2),
                 dimnames =
                   list(Income =
                          c("<5000", "5000-15000",
                            "15000-25000", ">25000"),
                        "Job Satisfaction" =
                          c("V_D", "L_S", "M_S", "V_S"),
                        Gender = c("Female", "Male"))))
## (Satisfaction categories abbreviated for convenience.)
ftable(. ~ Gender + Income, Satisfaction)
## Table 7.8 in Agresti (2002), p. 288.
mantelhaen.test(Satisfaction)
## See Table 7.12 in Agresti (2002), p. 297.
