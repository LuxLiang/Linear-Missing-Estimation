# . Simulations for normal data, 20% missing on average;
# . Sourced functions below should be in the same directory;
# . Sample size n=250, MC=1000 repeated simulations, random seed=1234;
# . X~standardized Gamma(5,1), Z~N_3(0,I_3), noise~standardized t(5);
# . Estimated seletion probabilities are bounded below from 1e-4 to avoid
#    invalid 'NaN' results;
# . Return:
#       . Beta_: a list storing all results of the called function
#       . bias_: the averaged bias of MC simulations
#       . epsd_: the empirical SE of MC simulations
#       . asysd_: the averaged asymptotic SE over MC simulations
#       . cover_: 95% coverages of beta




source("F_MAR.R")
source("CC_MAR.R")
source("IPW_MAR.R")
source("AIPW_MAR.R")



beta=c(0,0.5,1,-1,-0.5)
MC=1000

#alpha<-c(2.2,-0.9,-0.7,0,0) # 0.2
alpha<-c(0.5,-1,-0.5,0,0) # 0.4
n=250

####Full data analysis under MAR####
Beta0 <- full_MAR(xd="g",inno="t",alpha=alpha,n=n,seed=777)
bias0 <- Beta0$bias
epsd0 <- Beta0$epsd
#asysd0 <- (Beta0$asysd)/sqrt(MC*250)
cover0 <- Beta0$cover


####CC under MAR####
Beta1 <- CC_MAR(xd="g",inno="t",alpha=alpha,n=n,seed=777)
bias1 <- Beta1$bias
epsd1 <- Beta1$epsd
#asysd1 <- (Beta1$asysd)/sqrt(MC*250)
cover1 <- Beta1$cover


####IPW under MAR####
Beta2 <- ip_noaug_MAR(xd="g",inno="t",alpha=alpha,n=n,seed=777,trim=1e-4)
bias2 <- Beta2$bias
epsd2 <- Beta2$epsd
#asysd2 <- (Beta2$asysd)/sqrt(MC*250)
cover2 <- Beta2$cover


####PIPA under MAR####
Beta3 <- ip_aug_MAR(xd="g",inno="t",alpha=alpha,n=n,seed=777,trim=1e-4)
bias3 <- Beta3$bias
epsd3 <- Beta3$epsd
cover3 <- Beta3$cover




#save.image(file="xxxgamma.RData")





