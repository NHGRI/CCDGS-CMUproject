# partial correlaiton 

# INPUT 
#   Y matrix of dependent variables, row as dependent varaible, col as subject 
#     dimension    p*N  
#   X covariates row as subjects and col as covarates 
#     dimensiton   N*f   

# Linear regression as matrix system 
#     b = (X'X)^-1X'Y'
Residuals <- function(X,Y)
{
library(MASS)  
NormalEquationB <-  ginv(t(X)%*%X)%*%t(X)%*%t(Y)
Residual <- Y-t(X%*%NormalEquationB)
}