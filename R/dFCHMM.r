pseudo_resid = function(dat,res){
  pSig = res$Sigma.hat
  pmu = res$mu.hat
  num_state = dim(pmu)[1];num_state
  NSUB=dim(dat)[3]
  TIME = dim(dat)[1]
  resu1 <-matrix(NA,TIME,NSUB)
  for (sub in 1:NSUB)for (t in 1:TIME){
    for(i in 1:num_state) pSig[,,i]=res$Sigma.hat[,,i]*res$u.all[t,i,sub]
    for (i in 1:num_state) pmu[i,]=res$mu.hat[i,]*res$u.all[t,i,sub]
    wmu = colSums(pmu)
    varw = rowSums(pSig, dims = 2)-wmu%*%t(wmu)
    resu1[t,sub]=mahalanobis(dat[t,,sub], center=wmu, cov = varw) 
  }
  return(resu1)
}

hmm_viterbi <- function(res){
  n <- dim(res$u.all)[1]
  nsub <- dim(res$u.all)[3]
  iv <- matrix(NA, nsub, n)
  for (Sub in 1:nsub){
    probs <- drop(res$u.all[,,Sub])
    xi <- matrix(NA,n,dim(res$u.all)[2])
    foo <- res$delta.hat[Sub,]*probs[1,]
    xi[1,] <- foo/sum(foo)
    for (i in 2:n){
      if(is.na(dim(res$A.hat)[3])){
        foo <- apply(xi[i-1,]*res$A.hat, 2, max)*probs[i,]
      }
      else{
        foo <- apply(xi[i-1,]*res$A.hat[,,Sub], 2, max)*probs[i,]
      }
      xi[i,]<- foo/sum(foo)
    }
    iv[Sub,n] <- which.max(xi[n,])
    for (i in (n-1):1){
      if(is.na(dim(res$A.hat)[3])){
        iv[Sub, i]<- which.max(res$A.hat[,iv[Sub,i+1]]*xi[i,])
      }
      else{
        iv[Sub, i]<- which.max(res$A.hat[,iv[Sub,i+1],Sub]*xi[i,])
      }
    }
  }
  return(iv)
}

get.metastate = function(state.sequence,metastates){
  # Convert states to metastates
  num.meta = length(metastates)
  metavec = rep(NA,num.meta)
  for(m in 1:num.meta){
    metavec[metastates[[m]]] = m
  }
  N = nrow(state.sequence)
  meta.hsmm = matrix(NA,nrow=N,ncol=ncol(state.sequence))
  for(n in 1:N){
    for(time in 1:ncol(state.sequence)){
      meta.hsmm[n,time] = metavec[state.sequence[n,time]]
    }
  }
  return(meta.hsmm)
}

relabel.sequence = function(state.sequence,new.order){
  new.sequence = matrix(NA,nrow=nrow(state.sequence),ncol=ncol(state.sequence))
  for(n in 1:nrow(state.sequence)){
    new.sequence[n,] = new.order[state.sequence[n,]]
  }
  return(new.sequence)
}

fit.hmm = function(Xt,num.states, ATTEMPTS=5,tol=1e-4,maxiter=2500,verbose=FALSE){
  hmm.attempts = hmm_attempts(Xt,num.states,ATTEMPTS=ATTEMPTS,tol=tol,maxiter=maxiter,verbose=verbose)
  return(hmm.attempts[[lapply(hmm.attempts,extract.likelihood) %>% 
                         unlist() %>% 
                         which.max()]])
}

hmm_attempts = function(Xt,num.states, ATTEMPTS=5,tol=1e-4,maxiter=2500,verbose=FALSE){
  P = dim(Xt)[2]
  result.em.attempts = list()
  for(attempt in 1:ATTEMPTS){
    if(verbose){
      print(sprintf('Working on HMM, attempt = ...%g',attempt))
    }
    r.summary = summary(c(apply(Xt,3,cor)))[c(2,5)]
    sigma.summary = summary(c(apply(Xt,3,function(x){sqrt(diag(cov(x)))})))
    R = runif(num.states,r.summary[1],r.summary[2])
    sigma = runif(num.states,sigma.summary[1],sigma.summary[2])
    Sigma.init = array(NA,dim=c(P,P,num.states))
    for(s in 1:num.states){
      Rmat = matrix(R[s],nrow=P,ncol=P)
      diag(Rmat) = 1
      Sigma.init[,,s] = diag(rep(sigma[s],P))%*%Rmat%*%diag(rep(sigma[s],P))
      # Add a nugget for numerical stability
      lambda.min =  min(eigen(Sigma.init[,,s])$values)
      if(lambda.min < .1){
        Sigma.init[,,s] = Sigma.init[,,s] + diag(lambda.min + 0.1,P) 
      }
    }
    
    result.em.attempts[[attempt]] = tryCatch({hmm_cpp(Xt,
                                                       num_states = num.states,
                                                       Sigma = Sigma.init,
                                                       tol = tol,
                                                       maxiter = maxiter,
                                                       verbose=verbose)
                                              },
                                             error = function(cond){
                                               message(cond)
                                               return(result.em.attempts)
                                             })
  }
  return(result.em.attempts)
}

get.duration = function(s.hat,number_states){
  duration.mean = matrix(NA,nrow=nrow(s.hat),ncol=number_states)
  duration.sd = matrix(NA,nrow=nrow(s.hat),ncol=number_states)
  proportion = matrix(NA,nrow=nrow(s.hat),ncol=number_states)
  for(n in 1:nrow(s.hat)){
    runlength = rle(s.hat[n,])
    for(s in 1:number_states){
      duration.mean[n,s] = mean(runlength$lengths[which(runlength$values==s)])
      duration.sd[n,s] = sd(runlength$lengths[which(runlength$values==s)])
      proportion[n,s] = mean(s.hat[n,]==s)
    }
  }
  return(list(
    duration.mean = duration.mean,
    duration.sd = duration.sd,
    proportion = proportion)
  )
}

get.switching = function(s.hat){
  return(apply(s.hat,1,function(x)sum(diff(x) != 0)))
}

get.distances.Sigma = function(result,Sigma){
  num.states = result$A.hat%>%nrow()
  dist.mat = matrix(NA,nrow=num.states,ncol=num.states)
  for(j in 1:num.states){
    for(k in 1:num.states){
      dist.mat[j,k] = HS.norm(result$Sigma.hat[,,j]-Sigma[[k]])
    }
  }
  return(dist.mat)
}

permute.orders = function(dist.mat){
  num.states = ncol(dist.mat)
  states.pool = 1:num.states
  new.order = rep(NA,num.states)
  for(j in 1:ncol(dist.mat)){
    new.order[j] = which.min(dist.mat[j,])
    states.pool = states.pool[-which(states.pool == new.order[j])]
    dist.mat[,setdiff(1:num.states,states.pool)] = NA
  }
  return(new.order)
}

extract.likelihood = function(results){
  tryCatch(
    colSums(results$llk)[length(colSums(results$llk))],
    error = function(e){return(NA)})
}

HS.norm = function(A,B=NULL){
  if(is.null(B)){
    return(sqrt(sum(c(A^2))))
  }
  else{
    return(sqrt(sum(diag(t(A)%*%B))))
  }
}

KL.divergence = function(Sigma1,Sigma2,Omega1=NULL,Omega2=NULL){
  if(is.null(Omega1)){
    Omega1 = solve(Sigma1)
  }
  if(is.null(Omega2)){
    Omega2 = solve(Sigma2)
  }
  return(sum(diag((Sigma1-Sigma2)%*%(Omega2-Omega1))))
}

get.distances = function(Sigma,Omega=NULL){
  # Sigma is a list of covariance matrices
  # Omega is a list of precision matrices
  K = dim(Sigma)[3]
  if(is.null(Omega)){
    Omega = vector("list",length = K)
    for(k in 1:K){
      Omega[[k]] = solve(Sigma[,,k])
    }
  }
  dist.mat = matrix(999,nrow=K,ncol=K)
  for(j in 1:(K-1)){
    for(k in (j+1):K){
      dist.mat[j,k] = KL.divergence(Sigma[,,j],Sigma[,,k],Omega[[j]],Omega[[k]])
      dist.mat[k,j] = dist.mat[j,k]
    }
  }
  return(dist.mat)
}

IC.hmm = function(result){
  llk = extract.likelihood(result)
  num.states = nrow(result$A.hat)
  P = dim(result$Sigma.hat)[1]
  TIME = dim(result$u.all)[1]
  N = dim(result$u.all)[3]
  # num.params = num.states*(choose(P,2)+2*P) + num.states^2 # Common TPM
  num.params = num.states*(choose(P,2)+2*P) + N*num.states^2 # Subject-specific TPM
  ijhclss = matrix(NA,N,TIME) #alternatively
  iv = hmm_viterbi(result)
  for (i in 1:N) for (j in 1:TIME){
    z = iv[i,j]
    ijhclss[i,j]=-log(result$u.all[j,z,i])
  }
  hclss= sum(ijhclss)
  # hclss = -sum(c(result$u.all)*log(result$u.all),na.rm = T)
  return(list(AIC = -2*llk+2*num.params,
              BIC = -2*llk + num.params*log(TIME*N),
              ICL = -2*llk + num.params*log(TIME*N) + 2*hclss))
}

