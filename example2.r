#------Test for the case 
#-------y = 6.1 * ((xx[,1])^{2}+1)^{-1}+ 4.1 * xx[,3]*xx[,4] + -3.7 * xx[,5] - 4.2 * xx[,2] -------------------------------



SD <- 2.25;
Kn <- floor(sqrt(n));
q_ <- 5; #number of relevent variables
n <- 400; 	#sample size
p <- 5000; #number of variables
p_ <- p - q_;
nor <- (3 / (4 * q_))^{1/2};


#xx_re are the relevent variables; xx_ir are the irelevent ones
xx_re <- array(rnorm(q_ * n), c(n, q_));
xx_ir <- array(rnorm(p_ * n, sd = sqrt(SD)), c(n, p_)) + nor * (xx_re[,1] + xx_re[,3] + xx_re[,4] + xx_re[,5] + xx_re[,2]);# Alternatively

#Irelevent varialbes and yy_test are heavily correlated
xx <- cbind(xx_re, xx_ir);
#xx[,1]~ xx[,5] are relevent varialbes
yy_test <- matrix(6.1  * ((xx[,1])^{2}+1)^{-1}+ 4.1 * xx[,3]*xx[,4] + -3.7 * xx[,5] - 4.2 * xx[,2] + rnorm(n));




#-----------------------------
#--Begining of the process ---
results1 <- OHT_group(yy_test, xx, Kn, order_ = 4);


subset_1 <- rep(0, p);
subset_1[results1[[3]][1:10]] <- 1;# Get the stopped grouped OGA path. Only the first ten are considered!.

y_hat_1 <- residual(yy_test, xx, subset_ = subset_1, flag = 2);
y_residual_1 <- y_hat_1[[1]];

result_subset_1 <- indexing(y_hat_1[[2]]);#The first time subsetted indexes result!





#--Second round:----
results2 <- OHT_group(y_residual_1, xx, Kn, order_ = 4);

subset_2 <- rep(0, p);
subset_2[results2[[3]][1:10]] <- 1;# Get the stopped grouped OGA path


y_hat_2 <- residual(yy_test, xx, subset_ = (subset_2 + result_subset_1), flag = 2);
y_residual_2 <- y_hat_2[[1]];

result_subset_2 <- indexing(y_hat_2[[2]]);#The second time subsetted indexes result!
sum(y_hat_1[[1]]^{2});
sum(y_hat_2[[1]]^{2});#results2 only capture irrelevent variable(s).




#--Round three we begin the second order correlation

results3 <- OHT_group(y_residual_2^{2}, xx, Kn);

subset_3 <- rep(0, p);
subset_3[results3[[3]][1:10]] <- 1;# Get the stopped  grouped OGA path


y_hat_3 <- residual(yy_test, xx, (subset_3 + result_subset_2), flag = 2); #All of the previous significant variables! AND yy_test!!!
y_residual_3 <- y_hat_3[[1]];

result_subset_3 <- indexing(y_hat_3[[2]]);#The third time subsetted indexes result!





#--Foruth!!
results4 <- OHT_group(y_residual_3^{2}, xx, Kn);

subset_4 <- rep(0, p);
subset_4[results4[[3]][1:10]] <- 1;# Get the stopped  grouped OGA path

y_hat_4 <- residual(yy_test, xx, (subset_4 + result_subset_3), flag = 2, order_ = 4); #All of the previous significant variables! AND yy_test!!!
y_residual_4 <- y_hat_4[[1]];

result_subset_4 <- indexing(y_hat_4[[2]]);#The third time subsetted indexes result!
sum(y_hat_3[[1]]^{2});
sum(y_hat_4[[1]]^{2});#results2 only capture irrelevent variable(s).





#Have a look at how's the preformance before the 'final tuning'
E <- residual(yy_test,xx, result_subset_4, flag = 2, order_ = 4); 
sum(E[[1]]^{2}) / n;


#Final tuning for better fitting.
regressors <- power_expansion(y_hat_4[[3]], order_ = 5)[,-1];
regressors <- regressors[,colnames(unique(round(regressors, 6), MARGIN = 2))];#There might be some repeating columns.

result <- OHT(yy_test, regressors, Kn = (length(regressors[1,])));


re_ <- as.logical(result[[2]]);
colnames(regressors)[re_];
re_xx <- cbind(1, regressors[,re_]);

sum((yy_test - re_xx %*% solve(t(re_xx) %*% re_xx) %*% t(re_xx) %*% yy_test)^{2}) / n;#Sample size nomalized sum of square error
#The results are the relevent variables. Numbers follows a dash mean power.


#Remark:
#Numerically too large! Besides that everythin is fine. Delete variables with
#high power can temporarily resolve this problem.


#Remark:
#The case 
#----  y = 5.1  * xx[,2] * ((xx[,1])^{2}+1)^{-1}+ 4.1 * xx[,3]*xx[,4] + -3.7 * xx[,5] - 7.2 * xx[,2] 
#----  can be quite chanllenging