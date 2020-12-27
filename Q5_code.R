# ----------- Title: Simulation Code for PoisNonNor--------------# 

# Library packages
require(PoisNonNor)
require(moments)
require(boot)
require(matrixcalc)
source("./PoisBinOrdNonNor.R")

# Package Problem 
# Forget to delete RNG_P_NN


# ----------- Functions to calculate estimates and CI--------------# 
# Function to cumpute point estimates for mean, variance, skewness and kurtosis for non-nomarl distribution;
non.normal.est = function(data){
    c(mean(data), var(data), skewness(data), kurtosis(data)-3)
}
# Function to compute 95% CI
# 1.Poisson 95% CI (variance stabilization)
poi.ci.f = function(est, n.obs){
    (c(-1.96, 1.96)*sqrt(1/(4*n.obs)) + sqrt(est))^2
}
# 2. Mean Variance Skewness (Kurtosis-3) 95% CI (bootstrap percentile)
nonnor_f = function(data, indices){
    d = data[indices]
    c(mean(d), var(d), skewness(d), kurtosis(d)-3)
}
nonnor.ci.f = function(data){
    re = boot(data = data, statistic = nonnor_f, R = 1000)
    mean.ci = boot.ci(re, type = "perc", index = 1)$percent[4:5]
    var.ci = boot.ci(re, type = "perc", index = 2)$percent[4:5]
    skew.ci = boot.ci(re, type = "perc", index = 3)$percent[4:5]
    kur.ci = boot.ci(re, type = "perc", index = 4)$percent[4:5]
    c(mean.ci, var.ci, skew.ci, kur.ci)
}
# 3.Correlation 95% CI (Fisher Transformation)
cor.ci.f = function(rho, n.obs){
    x = 1/2*log((1 + rho)/(1 - rho)) + c(-1.96, 1.96)*sqrt(1/n.obs)
    (exp(2*x) - 1)/(exp(2*x) + 1)
}






# ----------- Generate random correlation matrix within boundaries--------------# 
# Generate one correlation is psd within bounds and intermediate matrix is psd
gen.corr = function(para.list){
    pois.list = para.list[[1]]
    nonn.list = para.list[[2]]
    non_m = matrix(unlist(nonn.list), nrow = 3, byrow = TRUE)
    pmat = Param.fleishman(non_m[,c(3,4)])
    
    lamvec = unlist(pois.list)
    # Calculates the approximate upper and lower correlation bounds
    re = lower.upper.cors(no.pois = 3, no.bin = 0, no.ord = 0, no.nonn = 3,
                          pois.list = pois.list, nonn.list = nonn.list)
    L = re$min[upper.tri(re$min)]
    U = re$max[upper.tri(re$max)]
    stat = F
    while (stat == F) {
        rho = numeric(length(L))
        for (i in 1:length(L)) {
            rho[i] = runif(1, L[i], U[i])
        }
        cor.mat = matrix(0, nrow = 6, ncol = 6)
        cor.mat[upper.tri(cor.mat)] = rho
        cor.mat = cor.mat + t(cor.mat) + diag(1,6)
        diag(cor.mat) = 1
        if (is.positive.definite(cor.mat)){
            inter_cor = intercor.all(cor.mat, pmat, lamvec)
            stat = is.positive.definite(inter_cor)
        }
    }
    return(cor.mat)
}



# ----------- Create a simmulation function--------------# 

PoisNonNor.sim = function(n.obs, n.sim, cor.mat, para.list) {
    
    pois.list = para.list[[1]]
    nonn.list = para.list[[2]]
    pmat = matrix(unlist(nonn.list), nrow = 3, byrow = TRUE)
    lamvec = unlist(pois.list)
    rmat = pmat[,c(3,4)]
    mean.vec = pmat[,1]
    variance.vec = pmat[,2]
    
    TV = c(unlist(pois.list), unlist(nonn.list), cor.mat[upper.tri(cor.mat)])
    
    data_all = RNG.P.NN(lamvec = lamvec, cmat = cor.mat,
                    rmat = rmat, norow = n.obs*n.sim, 
                    mean.vec = mean.vec, variance.vec = variance.vec)
    
    
    est_m = matrix(nrow = n.sim, ncol = 30)
    cover_m = matrix(nrow = n.sim, ncol = 30)
    for (i in 1:n.sim) {
        # Generate data
        print(i)
        start_time <- Sys.time()
        data = data_all[(n.obs*(i-1)+1):(i*n.obs), ]
        
        # Separate data into pois and non-normal
        data.poi = data[, 1:3]
        data.nonor = data[, 4:6]
        
        # Caculate estimates and 95% CI for non-ordinal parameters
        est.poi = apply(data.poi, 2, mean) # Poisson(lambda) 
        ci.poi = sapply(est.poi, poi.ci.f, n.obs) # CI for Poisson
    
        # Calculate estimate and 95% CI for non-normal paramters
        est.nonor = sapply(c(list(data.nonor[,1]), list(data.nonor[,2]), 
                             list(data.nonor[,3])), non.normal.est, simplify = F)
        
        ci.nonnor = sapply(c(list(data.nonor[,1]), list(data.nonor[,2]),
                             list(data.nonor[,3])), nonnor.ci.f, simplify = F)
        
        # Caculate estimates and 95% for all correlation
        est.cor = cor(data)
        est.cor.vec = est.cor[upper.tri(est.cor)] # only upper traingle
        ci.cor = sapply(est.cor.vec, cor.ci.f, n.obs)
        
        # Organize Results
        est_m[i, ] = c(est.poi, unlist(est.nonor), est.cor.vec)
        est.ci = c(as.vector(ci.poi), unlist(ci.nonnor), as.vector(ci.cor))
        
        # Test whether CI cover TV 
        cover_vec = numeric(30)
        for (j in 1:30) {
            L = est.ci[2*j - 1]
            U = est.ci[2*j]
            if ((TV[j] >= L) && (TV[j] <= U)) {
                cover_vec[j] = 1
            }else cover_vec[j] = 0
        }
        cover_m[i, ] = cover_vec
        
        end_time <- Sys.time()
        t = end_time - start_time
        print(t)
    }
    
    # Summarize Result
    AE = apply(est_m, 2, mean)
    sd = apply(est_m, 2, sd)
    RB = (AE - TV)*100/TV
    SB = abs(AE - TV)/sd*100
    RMSE = sqrt(apply((est_m - TV)^2, 2, mean))
    CR = apply(cover_m, 2, mean)
    
    reuslt_table = data.frame(
        TV = round(TV,4),
        AE = round(AE, 4),
        SD = round(sd, 4),
        RB = round(RB, 4),
        SB = round(SB, 4),
        RMSE = round(RMSE, 4),
        CR = round(CR, 4)
    )
    
    name_dist = c("Poi1", "Poi2", "Poi3", 
                  "Non-nor1", "Non-nor2", "Non-nor3")
    names_cor = NULL
    for (i in 1:5) {
        for (j in (i+1):6) {
            names_cor = c(names_cor, paste(name_dist[i], name_dist[j]))
        }
    }
    names = c("lambda1", "lambda2", "lambda3",
              "mean1","var1", "skew1", "kur1",
              "mean2", "var2", "skew2", "kur2",
              "mean3", "var3", "skew3", "kur3",
              names_cor)
    
    row.names(reuslt_table) = names
    return(reuslt_table)
}


# ----------- Setting 1. CorA --------------# 
# Poisson (1, 4, 10)
# Non-normal Beta(4,2) chi-squared (32) Laplace(1,3) 
# Sample size 100 and 10000
set.seed(1)
para.list.1 = list(pois.list = list(2, 5, 10),
                   nonn.list = list(c(0.6667, 0.0317, -0.4677, -0.3750),
                                    c(32, 64, 0.5, 0.375), c(1, 2, 0, 3)))
cor.mat = gen.corr(para.list.1)
re1 = PoisNonNor.sim(100, 1000, cor.mat, para.list.1)
write.csv(re1, "result_s1.csv")

re2 = PoisNonNor.sim(10000, 1000, cor.mat, para.list.1)
write.csv(re2, "result_s2.csv")

# ----------- Setting 2. CorB --------------# 
# Poisson (2, 5, 10)
# Non-normal Beta(4,2) chi-squared (32) Laplace(1,3) 
# Sample size 100 and 10000
set.seed(2)

cor.mat = gen.corr(para.list.1)
re3 = PoisNonNor.sim(100, 1000, cor.mat, para.list.1)
write.csv(re3, "result_s3.csv")

re4 = PoisNonNor.sim(10000, 1000, cor.mat, para.list.1)
write.csv(re4, "result_s4.csv")

# ----------- Setting 3. CorC --------------# 
# Poisson (10, 20, 30)
# Normal mixture distribution  chi-squared (16) 
# Exponential (2) -> (0.5,0.25,2,6)

set.seed(1)
para.list.2 = list(pois.list = list(10, 20, 30),
                   nonn.list = list(c(2, 2, 0, -.9582),
                                    c(16, 32, 0.7071, 0.75), c(0.5,0.25,2,6)))
cor.mat = gen.corr(para.list.2)
re5 = PoisNonNor.sim(100, 1000, cor.mat, para.list.2)
write.csv(re5, "result_s5.csv")

re6 = PoisNonNor.sim(10000, 1000, cor.mat, para.list.2)
write.csv(re6, "result_s6.csv")


# ----------- Setting 4. CorD --------------# 
# Poisson (10, 20, 30)
# Normal mixture distribution  chi-squared (16) 
# Exponential (2) -> (0.5,0.25,2,6)

set.seed(4)
cor.mat = gen.corr(para.list.2)
re7 = PoisNonNor.sim(100, 1000, cor.mat, para.list.2)
write.csv(re7, "result_s7.csv")

re8 = PoisNonNor.sim(10000, 1000, cor.mat, para.list.2)
write.csv(re8, "result_s8.csv")
