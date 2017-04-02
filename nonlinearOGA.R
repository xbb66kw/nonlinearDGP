rm(list=ls()); # clear all variables
graphics.off() #
library(MASS)# for function 'ginv'
#library(mvtnorm) 
#library(xtable);
#library("profvis");
library(Rcpp);
library(inline);

###############################################################








################
################
##Rcpp::::;###

src_2 <- '
 arma::mat y = as<arma::mat>(group_x);
 arma::colvec U = as<arma::vec>(u);
 arma::colvec result = y * y.t() * U;
 return(Rcpp::List::create(result));
'

Multiplier_C <- cxxfunction(signature(group_x = "numeric", u = "numeric"), body = src_2, plugin = "RcppArmadillo");


src_3 <- '
 arma::mat y = as<arma::mat>(group_x);
 
 
 int order_ = Rcpp::as<int>(order__);
 arma::mat result;
 //y.cols(1,1) % y.cols(2,2) % y.cols(0,0);

 y.insert_cols(0, arma::ones<arma::vec>(y.n_rows));
 
   if (order_ == 4) {
      for (int i = 0; i < y.n_cols; i ++) {
          for (int j = i; j < y.n_cols; j ++) {
              for (int k = j; k < y.n_cols; k ++) {
                  for (int h = k; h < y.n_cols; h ++) {
                      result.insert_cols(result.n_cols, y.cols(i,i) % y.cols(j,j) % y.cols(k,k) % y.cols(h,h));
                  }
              }
          }
      }
   }
 //arma::colvec result = 0;
 
 result = result.cols(1, (result.n_cols-1));
 return(Rcpp::List::create(result));
'

Expander_C <- cxxfunction(signature(group_x = "numeric", order__ = "numeric"), body = src_3, plugin = "RcppArmadillo");





#####Naming function################
####################################
fn <- function(i, K, r) {
    results <- list();
    
    if (r > 1) {
        r <- r - 1;
        item_ <- 1;
        for (j in i:K) {
            
            mid_ <- fn (j, K, r);
            
            for (S in mid_) {
                results[[item_]] <- c(j, S);
                item_ <- item_ + 1;
            }
        }
    } else if (r == 1) {
        item_ <- 1;
        for (j in i:K) {
            results[[item_]] <- j;
            item_ <- item_ + 1;
        }
    }
    
    return(results);

}

naming <- function(name_list, order_ = 4) {
    name_list <- c('NA', name_list);
    K <- length(name_list);
    r <- order_;
    results <- fn(1, K, r);
    output_ <- '';
    for (elem in results) {
        output_ <- c(output_, paste(name_list[elem], collapse = ','));
    }
    output_ <- output_[c(-1,-2)];
    return(output_);
}





OHT = function(y, X, Kn) {
    n <- nrow(X) ; p <- ncol(X); y <- y - mean(y); for(i in 1:p) {X[,i] <- X[,i] - mean(X[,i]);}     #Normalized
    
    #-OGA ------------------------------------------------------------
    jhat <- ehat1 <- integer(Kn); Xbase <- array(0,c(n,Kn));  u <- y; xnorms <- sqrt(colSums(X^2));
    for(k in 1:Kn) {  # Gram-Schmidt orthogonalization to avoid the inverse matrix calculation in Xbase
        SSE <- colSums( abs(t(u) %*% X) ) / xnorms; #colSums changes matrix with one row to a vector.
        
        if(k > 1) { SSE[jhat[1:(k-1)]] <- 0;} #Numerically not zero, manually set it to zero.
        
        jhat[k] <- which(abs(SSE) == max(abs(SSE)));
        
        if(k == 1) { 
            qk <- X[,jhat[1]] / sqrt(sum(t(X[,jhat[1]]) %*% X[,jhat[1]]));
            Xbase[,1] <- qk;
        } else { 
            rq <- X[,jhat[k]] - Xbase[,1:(k-1)] %*% t(Xbase[,1:(k-1)]) %*% X[,jhat[k]] ;
            qk <- rq / sqrt(sum(t(rq)%*%rq)); 
            Xbase[,k] <- qk;
        }
        
        u <- u - qk %*% t(qk) %*% u; 
        ehat1[k] <- abs(det( t(u) %*% u )) / ncol(y);#For HDIC
    }
    
    # OGA + HDIC
    
    HDBIC  <- n * log(ehat1[1:Kn] / (n * ncol(y))) + ncol(y) * (1 : Kn) * log(p) * log(n);
    mHDBIC  <- min(HDBIC); 
    j_hdbic <- jhat[1:which(HDBIC == mHDBIC)];
    
    # OGA + HDBIC + Trim  
    j_trim <- jremain <- j_hdbic;
    for( k in 1:length(j_hdbic)) {
        j1 <- jremain[-k]; 
        ehat1 <- array(0, c(n, ncol(y)));
        
        for(i in 1:ncol(y)) { ehat1[,i] <- y[,i] - mean(y[,i]);}     
        
        if(length(j1) > 0) { ehat1 <- y-X[,j1] %*% ginv(t(X[,j1]) %*% X[,j1]) %*% t(X[,j1]) %*% y;}
        
        HDBIC_trim <-  n * log( abs(det(t(ehat1) %*% ehat1) / ncol(y)) / n ) + ncol(y) * (length(j_hdbic) - 1) * log(n) * log(p);
        
        if(HDBIC_trim < mHDBIC) { j_trim[k] <- 0;}
    }
    
    j_trim <- j_trim[-which(j_trim==0)];
    
    if(length(j_trim)==0) { j_trim <- j_hdbic;}
    
    #----------------------------------------------------------------------------------------------
    trim_m <- integer(p);  
    trim_m[j_trim] <- 1;
    #OGA_HDBIC_Trim=c( sum(trim_m[which(beta0.rowsum!=0)])  , sum(trim_m[-which(beta0.rowsum!=0)])  )
    #----------------------------------
    #return(list( "relevant selected number"= sum(trim_m[which(beta0.rowsum!=0)]) ,  "overfitting number"= sum(trim_m[which(beta0.rowsum==0)])  , jhat      ))
    return(list(  jhat, trim_m));
    # relevant selected number : how many variables in relevant index are the selected model
    # overfitting number       : how many variables in irrelevant index are the selected model
}



####ENd of ex



##################
###########################
###########################
###########################
###########################Grouped OGA------------------


grouping <- function(x, order_ = 4) {
 M <- x;
 if (order_ > 1) {
  for (i in 2:order_) {
   M <- cbind(M, x^{i});
  }
 }
 return(M);
}


orthogonalizing <- function(x, base = NULL) {#return orthognalized x, wrt to base.
    M <- 0;
 
    if(!is.null(base)) { # assume base is orthogonalized itself
  
        M <- x[,1] - base %*% solve(t(base) %*% base) %*% t(base) %*% x[,1];
        M[,1] <- M[,1] / sqrt(t(M[,1]) %*% M[,1]);
        if (length(x[1,]) > 1) {
            for (k in 2:length(x[1,])) {
                #qk <- x[,k] - base %*% solve(t(base) %*% base) %*% t(base) %*% x[,k];
        
                qk <- x[,k] - base  %*% t(base) %*% x[,k]; #assume base is orthogonalized itself
        
                #M <- cbind(M, qk - M[,(1:k-1)] %*% solve(t(M[,(1:k-1)]) %*% M[,(1:k-1)]) %*% t(M[,(1:k-1)]) %*% qk);
        
                M <- cbind(M, qk - M[,(1:k-1)] %*% t(M[,(1:k-1)]) %*% qk);
                M[,k] <- M[,k] - mean(M[,k]);
                M[,k] <- M[,k] / sqrt(t(M[,k]) %*% M[,k]); #normalized
        
            }
        }
    } else {
        M <- Orthogonalizing_C(x)[[1]];
    }
    
    return(M);
}

src_orthogonalizing <- '
    arma::mat X = as<arma::mat>(x);
    X.cols(0,0) = X.cols(0,0) - arma::ones(X.n_rows,1) * mean(X.cols(0,0));
    arma::mat M = X.cols(0,0) / sqrt(arma::as_scalar((X.cols(0,0).t() * X.cols(0,0))));
    
    if (X.n_cols > 1) {
        for (int k = 1; k < X.n_cols; k ++) {
            M.insert_cols(M.n_cols, X.cols(k,k) - M * solve(M, X.cols(k,k)));
            
            M.cols(k,k) = M.cols(k,k) - arma::ones(M.n_rows,1) * mean(M.cols(k,k));
            M.cols(k,k) = M.cols(k,k) / sqrt(as_scalar(M.cols(k,k).t() * M.cols(k,k)));
        }
    }
    
    return(Rcpp::List::create(M));
'

Orthogonalizing_C <- cxxfunction(signature(x = 'numeric'), body = src_orthogonalizing, plugin = 'RcppArmadillo');

cp <- 1;
OHT_group <- function(y, X, Kn, order_ = 4) {
 
 n <- nrow(X); 
 p <- ncol(X);
 y <- y - mean(y); 
 ####MRIC
 ehat <- array(0, c(n, Kn));
 
 #######
 
 

   #beta0.rowsum = rowSums(beta0)   
 #-OGA ------------------------------------------------------------
 jhat <- ehat1 <- integer(Kn); 
 Xbase <- array(0, c(n, Kn * order_));  
 u <- y;
 #xnorms <- sqrt(colSums(X^2));
 #X <- A;
 #p <- 4
 rq <- apply(X, FUN = grouping, MARGIN = 2);
 GroupList <- as.list(rep(0, p))#construct list of groups;
 for (i in 1:length(rq[1,])) {
  GroupList[[i]] <- array(rq[,i], c(n, order_));
 }
 GroupList <- lapply(GroupList, FUN = orthogonalizing);
 
 for(k in 1:Kn){  # Gram-Schmidt orthogonalization to avoid the inverse matrix calculation in Xbase
  #SSE <- colSums( abs(t(u)%*%X) )/xnorms; #colSums changes matrix with one row to a vector.
  SSE_results <- rep(0, p);
  
  for (i in 1:p) {
   group_x <- GroupList[[i]];
   
   #U <- group_x %*% solve(t(group_x) %*% group_x) %*% t(group_x) %*% u; # projected vector
   #U <- group_x %*% t(group_x) %*% u; # projected vector
   if (length(dim(group_x)) <= 1) {
    U <- Multiplier_C(matrix(group_x), u)[[1]];
   } else {
    U <- Multiplier_C(group_x, u)[[1]];
   }
   
   #U <-  t(group_x) %*% u; # projected vector#test ######Yeah....this is the speed problem. Can be fixed by Rcpp.
   SSE_results[i] <- sum(U^{2});
  }
  
  
  if(k > 1) {
   SSE_results[jhat[1:(k-1)]] <- 0;
  } #Numerically not zero, manually set it to zero.
  
  jhat[k] <- which(SSE_results == max(SSE_results));
  
  if(k == 1) { 
   if (order_ == 1) { 
    qk <- matrix(GroupList[[jhat[1]]]);
   } else {
    qk <- orthogonalizing(GroupList[[jhat[1]]]);
   }
   #X[,jhat[1]] / sqrt( sum(t(X[,jhat[1]])%*%X[,jhat[1]]) ); 
   Xbase[,(1:order_)] <- qk;
  } else {
   m1 <- (k-1) * order_ + 1;
   m2 <- k * order_;
   
   if (order_ == 1) { 
    qk <- orthogonalizing(matrix(GroupList[[jhat[k]]]), Xbase[,1:(m1-1)]);
   } else {
    qk <- orthogonalizing(GroupList[[jhat[k]]], Xbase[,1:(m1-1)]);
   }
   
   Xbase[,m1:m2] <- qk;
  }
  u <- u - qk %*% t(qk) %*% u; 
  ####MRIC
  ehat[,k] <- u^{2};
  
  ######
  ehat1[k] <- abs(det( t(u) %*% u )) / ncol(y);#For HDIC
 }

 
 
 
 
 ###HDBIC
 ######
 HDBIC  <- n * log(  ehat1[1:Kn] /( n * ncol(y) ) ) + cp * (1 : Kn) * log(p) * log(n);
 mHDBIC  <- min(HDBIC); 
 j_hdbic <- jhat[1:which(HDBIC==mHDBIC)];

 # OGA + HDBIC + Trim  
 j_trim <- jremain <- j_hdbic;
 ehat2 <- array( 0, c(n, ncol(y)));
 
 for( k in 1:length(j_hdbic)) {
  j1 <- jremain[-k]; 
  
  if ( length(j1) > 0 ) {
   X_grouped <- do.call(cbind, GroupList[j1]);
   ehat2 <- y - X_grouped %*% ginv(t(X_grouped) %*% X_grouped) %*% t(X_grouped) %*% y; 
  }   
   
  HDBIC_trim <- n * log( abs(det(t(ehat2) %*% ehat2)/ncol(y))/n ) + ncol(y) * ( length(j_hdbic)-1 ) * log(n) * log(p);
   
  if( HDBIC_trim < mHDBIC ) {
   j_trim[k] <- 0;
  }
 }

 j_trim <- j_trim[-which(j_trim==0)]; 

 if(length(j_trim)==0) {
  j_trim <- j_hdbic;
 }

 #----------------------------------------------------------------------------------------------
 trim_m <- integer(p);  
 trim_m[j_trim] <- 1;
  #OGA_HDBIC_Trim=c( sum(trim_m[which(beta0.rowsum!=0)])  , sum(trim_m[-which(beta0.rowsum!=0)])  )
 #----------------------------------
  #return(list( "relevant selected number"= sum(trim_m[which(beta0.rowsum!=0)]) ,  "overfitting number"= sum(trim_m[which(beta0.rowsum==0)])  , jhat      ))
 return(list(  jhat, trim_m, j_hdbic));
  # relevant selected number : how many variables in relevant index are the selected model
 # overfitting number       : how many variables in irrelevant index are the selected model
}



















############################
###########################
#Residual finder!
residual <- function(y, x, subset_ = rep(1,length(x[1,])), order_ = 4, flag = 1) {
                xx <- x[,as.logical(subset_)];
    
    if(!is.matrix(xx)) {
        xx <- matrix(xx);
    }
    
    xx <- expander(xx, subset_ = as.logical(subset_), order_ = order_); # Possible warning for length too long.
    
    Kn <- sqrt(length(xx[1,]));
    sub_results <- OHT(y, xx, Kn);
    sub_ <- sub_results[[2]];
    
    
    
    
    name_ <- NULL #Since extract colnames from a vector results in NULL
    if (sum(as.logical(sub_)) == 1) {
        name_ <- colnames(xx)[as.logical(sub_)]
    }
    
    #OHT will defaulty regerss on 1.
    xxx <- cbind(xx[,as.logical(sub_)], 1);
    xx <- xx[,as.logical(sub_)];
    
    Y <- y - xxx %*% solve(t(xxx) %*% xxx) %*% t(xxx) %*% y;
    if (flag == 1) {
        return(Y);
    } else if (flag == 2){
        if (is.null(name_)) {
            return(list(Y, colnames(xx), xx));
        } else {
            return(list(Y, name_, xx));
        }
        
    }
}
###########################
# Expanding function:
expander <- function(xx, order_ = 4, subset_ = rep(1,length(xx[1,]))) {
                
    #print(subset_);
    if(length(xx[1,]) > 15) { 
        print("Too large a matrix, set order_ to 3");
        order_ = 3;
    }
    
    subset_ <- as.logical(subset_);
    
    subset_xx <- xx;
    
    M <- Expander_C(subset_xx, order_)[[1]];
    name_list <- which(subset_);
    
    colnames(M) <- naming(name_list, order_);
     
    return(M);
    }

#########################
#indexing function
indexing <- function(x) {
    process_ <- strsplit(paste(x, collapse = ','), ',')[[1]];
    process_ <- process_[which(process_ != 'NA')];
    sub_ <- unique(process_);
    result_subset <- rep(0, p);
    result_subset[as.integer(sub_)] <- 1;#The subsetted indexes result!
    return(result_subset);
}


power_expansion <- function(X, order_ = 7) {
   
    results <- rep(1, n);
    for (j in 1:length(X[1,])) {
        
        li_ <- strsplit(colnames(X)[j], ',')[[1]];
        li_ <- paste(li_[!(li_ == 'NA')], collapse = ',');
        
        x_ <- X[,j];
        ori_ <- li_;
        xx_ <- X[,j];
        
        for (k in 2:order_) {
            li_ <- cbind(li_, paste(ori_, '_', k));
            xx_ <- cbind(xx_, x_^{k});
        }
        
        colnames(xx_) <- li_;
        
        results <- cbind(results, xx_);
    }
   
    return(results);
}



