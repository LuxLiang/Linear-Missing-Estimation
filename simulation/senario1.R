# . Simulations for normal data, 20% missing;
# . Sample size n=250, MC=1000 repeated simulations, random seed=777;
# . X~N(0,1), U~N_3(0,I_3), noise~N(0,1);
# . Estimated seletion probabilities are bounded below from 1e-4 to avoid
#    invalid 'NaN' results;
# . Obtain:
#     . Beta_: a list storing all results of the called function
#     . bias_: the averaged bias of MC simulations
#     . epsd_: the empirical SE of MC simulations
#     . asysd_: the averaged asymptotic SE over MC simulations
#     . cover_: 95% coverages of beta



source("F_MAR.R")
source("CC_MAR.R")
source("IPW_MAR.R")
source("AIPW_MAR.R")




beta=c(0,0.5,1,-1,-0.5)
MC=1000

#alpha<-c(3,-0.1,-0.1,0,0)
#alpha<-c(2.2,-0.9,-0.7,0,0) # 0.2
alpha<-c(0.5,-1,-0.5,0,0) # 0.4
#alpha<-c(-0.5,-0.5,-0.5,0,0)
n=250

####Full data analysis under MAR####
Beta0 <- full_MAR(xd="n",inno="n",alpha=alpha,n=n,seed=777)
bias0 <- Beta0$bias #偏差
epsd0 <- Beta0$epsd # 经验标准误 
#asysd0 <- (Beta0$asysd)/sqrt(MC*250) # 平均渐进标准误
cover0 <- Beta0$cover # 覆盖概率


####CC under MAR####
Beta1 <- CC_MAR(xd="n",inno="n",alpha=alpha,n=n,seed=777)
bias1 <- Beta1$bias
epsd1 <- Beta1$epsd
#asysd1 <- (Beta1$asysd)/sqrt(MC*250)
cover1 <- Beta1$cover
########## 置信区间图 ###########
library(ggplot2)
lower = Beta1$lower_list # MC=1000
upper = Beta1$upper_list # MC=1000
Beta_real<-rep(c(0,0.5,1,-1,-0.5),c(1000,1000,1000,1000,1000))
x<-rep(1:1000,each=1,times=5)
dat_plot <- data.frame(cbind(x, Beta_real, lower, upper))
names(dat_plot) <- c("x", "y", "lower", "upper")
dat_plot$group <- rep(c("Beta_0", "Beta_1","Beta_2","Beta_3","Beta_4"), each = 1000)
ggplot(dat_plot, aes(x = x, color = group, fill = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.8) +  # alpha 修改透明度
  geom_line(aes(y = Beta_real))




####IPW under MAR####
Beta2 <- ip_noaug_MAR(xd="n",inno="n",alpha=alpha,n=n,seed=777,trim=1e-4)
bias2 <- Beta2$bias
epsd2 <- Beta2$epsd
#asysd2 <- (Beta2$asysd)/sqrt(MC*250)
cover2 <- Beta2$cover


####AIPW-SIM under MAR####
Beta3 <- ip_aug_MAR(xd="n",inno="n",alpha=alpha,n=n,seed=777,trim=1e-4)
bias3 <- Beta3$bias
epsd3 <- Beta3$epsd
#asysd3 <- (Beta3$asysd)/sqrt(MC*250)
cover3 <- Beta3$cover




###### save #####
#save.image(file="xxxn_250.RData")





