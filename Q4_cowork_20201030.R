##################################################### purpose: simulate mixed data containing poisson, binary and ordinal
## function : gen.PoisBinOrd
## https://www.rdocumentation.org/packages/PoisBinOrd/versions/1.4.2/topics/gen.PoisBinOrd
library(PoisBinOrd)
set.seed(20200824)
n = 100 # sample size = Number of variates
n.P = 2 # Number of Poisson variables
n.B = 2 # Number of binary variables
n.O = 2 # Number of ordinal variables
lambda.vec = sample(1:10, size = 2, replace=FALSE)# Rate vector for Poisson variables
prop.vec=runif(n=2)                               # Probability vector for binary variables
prop.list=list(c(0.3,0.6,0.7),c(0.2,0.3,0.5))     # A list of probability vectors for ordinal variables. (cdf-cutoff)

### Create overall correlation matrix, SIGMA
M = c(-0.05, 0.26, 0.14, 0.09, 0.14, 0.12, 0.13, -0.02, 0.17, 0.11,-0.04, 0.19, 0.10, 0.35, 0.39)
N = diag(6)  
N[lower.tri(N)] = M  # lower triangle 
corr.mat = N+t(N)    # symmetric matrix
diag(corr.mat) = 1   
corr.mat

# Validates the specified correlation matrix
validation.corr(n.P,n.B,n.O,corr.vec=NULL,corr.mat)

# checks if there are range violations among pp, po, pb, bb, bo, oo combinations
correlation.bound.check(n.P,n.B,n.O,lambda.vec,prop.vec,prop.list, corr.vec=NULL,corr.mat)

# print correlation limits among Poisson, binary, and ordinal variables and within themselves.
 correlation.limits(n.P,n.B,n.O,lambda.vec,prop.vec,prop.list)

# Computes the final intermediate correlation matrix, SIGMA*
final.corr.mat=overall.corr.mat(n.P,n.B,n.O,
                                lambda.vec,
                                prop.vec, 
                                prop.list,
                                corr.vec=NULL,corr.mat)
final.corr.mat

# checks if there are range violations among pp, po, pb, bb, bo, oo combinations
# correlation.bound.check(n.P,n.B,n.O,lambda.vec,prop.vec,prop.list, corr.vec=NULL,final.corr.mat)

# Generate mixed data
mymixdata = gen.PoisBinOrd(n,n.P,n.B,n.O,
                           lambda.vec,
                           prop.vec,
                           prop.list, 
                           final.corr.mat)
View(mymixdata)

# weak correlation
W_Corr = runif(n=15, min = 0.001,max = 0.399)
W_PN_Indicator = sample(c(-1,1), size = 15, replace = TRUE)
W_Corr_Final = W_Corr*W_PN_Indicator
W_Corr_Final
# strong correlation
S_Corr = runif(n=15, min = 0.600,max = 1.000)
S_PN_Indicator = sample(c(-1,1), size = 15, replace = TRUE)
S_Corr_Final = S_Corr*S_PN_Indicator
S_Corr_Final