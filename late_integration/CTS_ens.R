program_block_li <- function(prop1,prop2,path_dataset) { 

  # prop1 = l_pred$prop1
  # prop2 = l_pred$prop2

 CTS_EnsDeconv_LS <- function(What, mygamma = 0.001, maxit = 1000,
                                 mysigma = 0.4, mybeta = 0.1) {
      #What D*n*K
      #W n*K
      i <- 1
      n <- dim(What[[1]])[1]
      K <- dim(What[[1]])[2]
      W  = list()
      W[[i]] <- Reduce(`+`, What) / length(What)
      W_noK <- W[[i]][,-K, drop = FALSE]
      a<- rep(1,K)
      
      Delta_CT <-  rep(100,(K-1))
      pren <- norm(Delta_CT, type="2")
      e <- diff_grad<-diff<-100
      difff<-new <- ind<-sum(sapply(What, function(sce) {
        tett = sum(sapply(1:K, function(k) {
          a[k]*norm(sce[,k]-W[[i]][,k],type = "2")
        }))
        return(tett)
      }))
      while ((diff > 1e-5) && (i<maxit)) {
        pre <- new
        pren <- norm(Delta_CT, type="2")  
        Delta_CT <- sapply(1:(K-1), function(sk) {
          subgrad <- lapply(What, function(sce) {
            grad0 <- sce[,K] - (rep(1,n)-rowSums(W_noK))
            grad1 <- sce[,sk] - W_noK[,sk]
            res <- a[K]*(1/norm(grad0,type = "2"))*grad0 - a[sk]*(1/norm(grad1,type = "2"))*grad1
            return(res)
          })
          Delta <- Reduce("+", subgrad)
          return(Delta)
        })
        mk = 0
        mygamma = mybeta^mk
        Z <- W_noK -mygamma*Delta_CT
        
        bvec <-  c(-1,rep(0,(K-1)))
        Amat <- t(rbind(rep(-1,(K-1)),diag(1,(K-1),(K-1))))
        if(K>2) {
          W_noK_step <- t(sapply(1:n, function(j){
            quadprog::solve.QP(Dmat =diag(1,(K-1),(K-1)), dvec = Z[j,], Amat = Amat, bvec = bvec)$solution
          }))
        } else {
          W_noK_step <- t(t(sapply(1:n, function(j) {
            quadprog::solve.QP(Dmat =diag(1,(K-1),(K-1)), dvec = Z[j,], Amat = Amat, bvec = bvec)$solution
          })))
        }
        int_diff <-  myf(W_noK,What_f = What, K_f = K, a=a)-myf(W_noK_step, What_f = What, K_f = K, a=a)
        int_diff2 <- mysigma*sum(Delta_CT*(W_noK-W_noK_step))
        while (int_diff<int_diff2) {
          mk <- mk+1
          mygamma = mybeta^mk
          Z <- W_noK -mygamma*Delta_CT
          bvec <- c(-1,rep(0,(K-1)))
          Amat <- t(rbind(rep(-1,(K-1)),diag(1,(K-1),(K-1))))
          if(K>2) {
            W_noK_step <- t(sapply(1:n, function(j) {
              quadprog::solve.QP(Dmat =diag(1,(K-1),(K-1)), dvec = Z[j,], Amat = Amat, bvec = bvec)$solution
            }))
          } else {
            W_noK_step <- t(t(sapply(1:n, function(j) {
              quadprog::solve.QP(Dmat =diag(1,(K-1),(K-1)), dvec = Z[j,], Amat = Amat, bvec = bvec)$solution
            })))
          }
          int_diff <-  myf(W_noK,What_f = What, K_f = K, a=a)-myf(W_noK_step, What_f = What, K_f = K, a=a)
          int_diff2 <- mysigma*sum(Delta_CT*(W_noK-W_noK_step))
        }
        W_noK <- W_noK_step
        i <- i+1
        W[[i]] <- cbind(W_noK, 1-rowSums(W_noK))
        new <- sum(sapply(What, function(sce) {
          tett = sum(sapply(1:K,function(k) {
            a[k]*norm(sce[,k]-W[[i]][,k],type = "2")
          }))
          return(tett)
        }))
        diff_grad <- c(diff_grad, norm(Delta_CT, type="2"))
        diff <- abs(pre - new ) 
        difff <- c(difff,diff)
        ind <- c(ind,new)
        e = mean(abs(c (as.matrix(pre-W_noK))))
      }
      W_fin =  W[[which.min(ind)]]
      rownames(W_fin) = rownames(What[[1]])
      colnames(W_fin) = colnames(What[[1]])
      return(list(W = W_fin, difff = difff))
    }
    myf <- function(W_noK, What_f = What, K_f = K, a) {
      W <-  cbind(W_noK, 1-rowSums(W_noK))
      new <- sum(sapply(What_f, function(sce) {
        tett = sum(sapply(1:K_f,function(k) {
          a[k]*norm(sce[,k]-W[,k],type = "2")
        }))
        return(tett)
      }))
      return(new)
    }
    
    list_proportions <- list("block1" = t(prop1),
                             "block2" = t(prop2))
    prop <- t(CTS_EnsDeconv_LS(What = list_proportions, maxit = 1000)$W)
    if (any(colSums(prop) != 1)) { # Sum To One 
      prop <- sapply(1:ncol(prop), function(col)
        prop[,col] / sum(prop[,col]))
    }

  return(prop)
}
