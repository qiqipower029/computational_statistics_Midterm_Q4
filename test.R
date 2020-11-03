library(tidyverse)
library(PoisBinOrd)
set.seed(5)

para = list(n.P = 2, n.O =2, n.B = 2, 
             lambda.vec = c(0.02, 0.7), prop.vec = c(0.2, 0.3), prop.list = list(c(0.25, 0.50, 0.75), c(0.20, 0.45, 0.75)), 
             corr.mat = create_corr(para=para)$corr.mat)

# Computes the final intermediate correlation matrix, SIGMA*
final.corr.mat = overall.corr.mat(para$n.P,para$n.B,para$n.O,para$lambda.vec,para$prop.vec,para$prop.list, corr.vec=NULL, para$corr.mat)
para$final.corr.mat  = final.corr.mat

simulation.function = function(para, samplesize, n.sim) {
mymixdata = list(); para_esti=list(); para_true = list(); para_LCI = list(); para_UCI = list()
 
 for (i in 1:n.sim) {
    # Generate mixed data
    mymixdata[[i]] = gen.PoisBinOrd(n = samplesize, n.P = para$n.P, n.B = para$n.B, n.O = para$n.O,
                               lambda.vec = para$lambda.vec, prop.vec = para$prop.vec, prop.list = para$prop.list, 
                               final.corr.mat = para$final.corr.mat)
    # Poisson
    para_poi     = apply(mymixdata[[i]], MARGIN = 2, mean)[1:2] # Poisson Mean = Poisson Variance
    para_poi_LCI = para_poi - 1.96*sqrt(para_poi / samplesize)
    para_poi_UCI = para_poi + 1.96*sqrt(para_poi / samplesize)
    
    # Binary
    para_bin     = apply(mymixdata[[i]], MARGIN = 2, mean)[3:4]
    para_bin_LCI = para_bin - 1.96*sqrt(para_bin*(1-para_bin) / samplesize)
    para_bin_UCI = para_bin + 1.96*sqrt(para_bin*(1-para_bin) / samplesize)
    
    # Ordinal
    para_ordi = c(  cumsum( table(mymixdata[[i]][,5]))[1:3], cumsum( table(mymixdata[[i]][,6]))[1:3]  ) / samplesize # cumu proportion
    para_ordi_prop0 = c(   table(mymixdata[[i]][,5])[1], table(mymixdata[[i]][,6])[1]  ) / samplesize
    para_ordi_prop1 = c(   table(mymixdata[[i]][,5])[2], table(mymixdata[[i]][,6])[2]  ) / samplesize
    para_ordi_prop2 = c(   table(mymixdata[[i]][,5])[3], table(mymixdata[[i]][,6])[3]  ) / samplesize
    para_ordi_prop3 = c(   table(mymixdata[[i]][,5])[4], table(mymixdata[[i]][,6])[4]  ) / samplesize
    
    para_ordi_LCI1 = para_ordi[c(1,4)] - 1.96*sqrt(  para_ordi_prop0*(1-para_ordi_prop0)  /samplesize )
    para_ordi_LCI2 = para_ordi[c(2,5)] - 
                     1.96*sqrt(  (  para_ordi_prop0*(1-para_ordi_prop0) + para_ordi_prop1*(1-para_ordi_prop1) - 2*para_ordi_prop0*para_ordi_prop1  ) /samplesize )
    para_ordi_LCI3 = para_ordi[c(3,6)] - 
                     1.96*sqrt(  (     para_ordi_prop0*(1-para_ordi_prop0) +   para_ordi_prop1*(1-para_ordi_prop1) +   para_ordi_prop2*(1-para_ordi_prop2) 
                                   - 2*para_ordi_prop0*para_ordi_prop1     - 2*para_ordi_prop0*para_ordi_prop2     - 2*para_ordi_prop1*para_ordi_prop2) /samplesize )
    
    
    para_ordi_UCI1 = para_ordi[c(1,4)] + 1.96*sqrt(  para_ordi_prop0*(1-para_ordi_prop0)  /samplesize )
    para_ordi_UCI2 = para_ordi[c(2,5)] + 
                     1.96*sqrt(  (  para_ordi_prop0*(1-para_ordi_prop0) + para_ordi_prop1*(1-para_ordi_prop1) - 2*para_ordi_prop0*para_ordi_prop1  ) /samplesize )
    para_ordi_UCI3 = para_ordi[c(3,6)] + 
                     1.96*sqrt(  (     para_ordi_prop0*(1-para_ordi_prop0) +   para_ordi_prop1*(1-para_ordi_prop1) +   para_ordi_prop2*(1-para_ordi_prop2) 
                                   - 2*para_ordi_prop0*para_ordi_prop1     - 2*para_ordi_prop0*para_ordi_prop2     - 2*para_ordi_prop1*para_ordi_prop2) /samplesize )
    # Correlation
    corr = cor(mymixdata[[i]])
    to.upper = function(X) X[lower.tri(X,diag=TRUE)]
    corr_f = to.upper(corr)[c(2:6, 8:11, 13:15, 17:18, 20)]
    
    # Create function to calculate CI for rhos
    
    CI_low_cal_rho  = function(r, samplesize) { CL = 0.5*log((1+r)/(1-r)) - 1.96*sqrt(1/samplesize)  ; CI_low  = (exp(2*CL)-1)/(exp(2*CL)+1)
      return(CI_low)}
    
    CI_high_cal_rho = function(r, samplesize) { CU = 0.5*log((1+r)/(1-r)) + 1.96*sqrt(1/samplesize)  ; CI_high = (exp(2*CU)-1)/(exp(2*CU)+1)
      return(CI_high)
    }
    
    # Output
    para_esti[[i]] = as.vector(c(para_poi,para_bin,  para_ordi, corr_f))
    para_true[[i]] = c(para$lambda.vec, para$prop.vec, unlist(para$prop.list),  to.upper(para$final.corr.mat)[c(2:6, 8:11, 13:15, 17:18, 20)])
    para_LCI[[i]] = c(para_poi_LCI, para_bin_LCI, para_ordi_LCI1, para_ordi_LCI2, para_ordi_LCI3,  CI_low_cal_rho(corr_f , samplesize) )
    para_UCI[[i]] = c(para_poi_UCI, para_bin_UCI, para_ordi_UCI1, para_ordi_UCI2, para_ordi_UCI3, CI_high_cal_rho(corr_f , samplesize) )
  }

# Point Estimation
  Para_True = matrix( unlist(para_true), nrow=n.sim, byrow=TRUE)
  Para_Esti = matrix( unlist(para_esti), nrow=n.sim, byrow=TRUE)
  theta_true = apply( Para_True, MARGIN =2, mean)
  theta_esti_mean = apply( Para_Esti, MARGIN =2, mean)
  theta_esti_sd = apply( Para_Esti, MARGIN =2, sd)

# Accuracy and Precision 
  RB = (  (theta_esti_mean - theta_true)/theta_true        )*100
  SB = (  abs(theta_esti_mean - theta_true)/theta_esti_sd  )*100
  RMSE = sqrt(   apply( (Para_Esti - Para_True)^2, MARGIN =2, mean )   )
  
# Confidence Interval 
  LCI = matrix( unlist(para_LCI), nrow=n.sim, byrow=TRUE)
  UCI = matrix( unlist(para_UCI), nrow=n.sim, byrow=TRUE)
  para_true_forCI = c(para$lambda.vec, para$prop.vec, unlist(para$prop.list)[c(1,4,2,5,3,6)], to.upper(para$final.corr.mat)[c(2:6, 8:11, 13:15, 17:18, 20)])
  Para_True_forCI = matrix(rep(para_true_forCI,n.sim), nrow=n.sim, byrow=TRUE)
  
  CR = apply(  ifelse( ( Para_True_forCI - LCI > 0) & (UCI - Para_True_forCI > 0),1,0 )   ,  MARGIN =2, sum)  / n.sim
  CR = CR[c(1:4, 5,8, 6,9, 7,10, 11:25)]

  inference = data.frame(Para_True      = round(theta_true,4),
              Para_Esti_Mean = round(theta_esti_mean,4), 
              Para_Esti_SD   = round(theta_esti_sd,4),
              RB             = round(RB,4),            
              SB             = round(SB,4), 
              RMSE           = round(RMSE,4))
  
  
  return(list(result=inference, CR=CR) )
          
}

scenario_1 = simulation.function(para=para, samplesize= 100, n.sim=100)
scenario_1$result
scenario_1$CR

# 1. result: only one table (inference+ CR)
# 2. scenario 1 (n=100)/(n=10000) -- output-- csv
# 3. Poisson CI
# 4. para1, (para2, para3, para4)-- wait for Jun 

