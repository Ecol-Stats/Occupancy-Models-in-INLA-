inla.Occupancy_detCov <- function(X_det){
  
  if(class(X_det)=="list"){
    if(length(X_det)>10){
      warning("exceeded number of detection covariates, numerical issues may occur")
    }
    
    if(lapply(X_det, ncol)|>unlist()|>unique()|>length()>2){
      stop("inconsistent number of visits in provided detection covariates")
    }
    if(length(lapply(X_det, nrow) |> unlist() |> unique())>1){
      stop("inconsistent number of sites in provided detection covariates")
    }
    K<- lapply(X_det, ncol) |> unlist() |> max() # Max num of visits
    M<- lapply(X_det, nrow) |> unlist() |> unique() # Number of sites
    P <- length(X_det)
    
    if(lapply(X_det, ncol)|>unlist()|>unique()|>length()==2 & 
       1 %in% lapply(X_det, ncol)|>unlist()|>unique()){
      warning(paste("At least one covariate of dimension [",M,",1] has been provided, values for this covariate will be repeated over the max numver of visits",sep=""))
      for(l in which(lapply(X_det, ncol) |> unlist() < K)){
        X_det[[l]] <- do.call("cbind",replicate(K,X_det[[l]]))
        
      }
    }
    covariates <- do.call("cbind", lapply(1:K, function(i) {
      do.call("cbind", lapply(X_det, function(mat) mat[, i]))
    }))
    
  }
  
  if(is.data.frame(X_det)|is.matrix(X_det)){
    K<- ncol(X_det)
    M<- nrow(X_det)
    P <- 1
    covariates <- as.matrix(X_det)
  }
  
  X_mat <- matrix(NA,nrow=M,ncol=K*(P+1))
  X_mat[,seq(1,(K*(P+1)),by=(P+1))]<-1 # add Intercept at the begining of each visit-specific covariate matrix
  X_mat[, which(!(1:(K*(P+1)) %in% seq(1,(K*(P+1)),by=(P+1))))] <- covariates
  return(X_mat)
  
}