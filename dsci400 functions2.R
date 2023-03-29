pred.list2 <- list()
MMSE.list2 <- list()
actual.list2 <- list()


gamma <- 0.98^(1:dim(matrix.data)[1])
for (i in 1:length(gamma)) {
  if(gamma[i] < 10^-3){
    gamma[i] <- 0
  }
}




# traditional covarience matric 

for (stocks in stock_picks) {
  
  matrix.data <- get(paste0(stocks, sep = ".", "Hankel.Matrix"))[-1]
  
  q <- ncol(matrix.data) # because Q=M
  
  normailzed.matrix <- data.frame(matrix(nrow =nrow(matrix.data), ncol = q)) # generate empty dataframe that will hold all normalized values
  
  for (i in 1:nrow(matrix.data)) {
    for (j in 1:ncol(matrix.data)) {
      normailzed.matrix[i,j] <- matrix.data[i,j]/matrix.data[i,q]
    }
  }
  
  normailzed.centered.matrix <- data.frame(matrix(nrow = nrow(matrix.data), ncol = ncol(matrix.data)))
  
  for (i in 1:ncol(normailzed.matrix)) {
    normailzed.centered.matrix[i] <- normailzed.matrix[i] - mean(normailzed.matrix[[i]])
  }
  
  
  normailzed.centered.matrix <- normailzed.centered.matrix[-q]
  
  matrix.X <- as.matrix(normailzed.centered.matrix)
  
  
  gamma.diag <- gamma * diag(nrow(matrix.data))
  gamma.const <- (1-0.98)/(1-0.98^(nrow(matrix.data)-1))
  cov.sample <- gamma.const * (t(matrix.X) %*% gamma.diag %*% matrix.X)
  
  matrix.y <- cov.sample[,1:(ncol(cov.sample)/2)]
  matrix.z <- cov.sample[,-(1:ncol(cov.sample)/2)]
  
  cov.yy <- cov(matrix.y)
  cov.zz <- cov(matrix.z)
  
  cov.yz <- cov(matrix.y,matrix.z)
  cov.zy <- cov(matrix.z, matrix.y)
  
  
  
  eigenvectors.yy <- eigen(cov.yy)
  
  
  matrixYML <- eigenvectors.yy$vectors # this is our V matrix 
  
  matrixYML <- as.matrix(matrixYML)
  
  
  G <- solve(t(matrixYML[1:5,1:2])  %*%  matrixYML[1:5,1:2]) %*% t( matrixYML[1:5,1:2]) # equation 8 
  y <- cbind(as.numeric(cov.sample[1,1:5]))
  sigma.zw <- cov.zy[,1:5] %*% t(G) # equation 10
  sigma.ww <- G %*% cov.yy[,1:5] %*% t(G) # equation 11
  w <- G %*% y
  
  AA <- sigma.zw %*% solve(sigma.ww) %*% w
  
  
  
  
  
  
  pred.list2[[stocks]] <- as.numeric((AA + c(mean(normailzed.matrix[[6]]),
                                            mean(normailzed.matrix[[7]]),
                                            mean(normailzed.matrix[[8]]),
                                            mean(normailzed.matrix[[9]]),
                                            mean(normailzed.matrix[[10]])))*matrix.data[1,10])
  
  actual.list2[[stocks]] <- as.numeric(matrix.data[1,6:10]) 
  
  MMSE.list2[[stocks]] <- (as.numeric(matrix.data[1,6:10]) - as.numeric((AA + c(mean(normailzed.matrix[[6]]),
                                                                               mean(normailzed.matrix[[7]]),
                                                                               mean(normailzed.matrix[[8]]),
                                                                               mean(normailzed.matrix[[9]]),
                                                                               mean(normailzed.matrix[[10]])))*matrix.data[1,10]))^2
  
  
  
}


###################



# traditional covarience matric 

for (stocks in stock_picks) {
  
  matrix.data <- get(paste0(stocks, sep = ".", "Hankel.Matrix"))[-1]
  
  q <- ncol(matrix.data) # because Q=M
  
  normailzed.matrix <- data.frame(matrix(nrow =nrow(matrix.data), ncol = q)) # generate empty dataframe that will hold all normalized values
  
  for (i in 1:nrow(matrix.data)) {
    for (j in 1:ncol(matrix.data)) {
      normailzed.matrix[i,j] <- matrix.data[i,j]/matrix.data[i,q]
    }
  }
  
  normailzed.centered.matrix <- data.frame(matrix(nrow = nrow(matrix.data), ncol = ncol(matrix.data)))
  
  for (i in 1:ncol(normailzed.matrix)) {
    normailzed.centered.matrix[i] <- normailzed.matrix[i] - mean(normailzed.matrix[[i]])
  }
  
  
  normailzed.centered.matrix <- normailzed.centered.matrix[-q]
  
  matrix.X <- as.matrix(normailzed.centered.matrix)
  
  
  gamma.diag <- gamma * diag(nrow(matrix.data))
  gamma.const <- (1-0.98)/(1-0.98^(nrow(matrix.data)-1))
  cov.sample <- gamma.const * (t(matrix.X) %*% gamma.diag %*% matrix.X)
  

  
  matrix.y <- cov.sample[,1:(ncol(cov.sample)/2)]
  matrix.z <- cov.sample[,-(1:ncol(cov.sample)/2)]
  
  cov.yy <- cov.sample[1:5,1:5]
  cov.zz <- cov.sample[6:10,6:10]
  
  cov.yz <- cov.sample[1:5,6:10]
  cov.zy <- cov.sample[1:5,1:5]
  
  
  
  eigenvectors.yy <- eigen(cov.yy)
  
  
  matrixYML <- eigenvectors.yy$vectors # this is our V matrix 
  
  matrixYML <- as.matrix(matrixYML)
  
  
  G <- solve(t(matrixYML[1:5,1:2])  %*%  matrixYML[1:5,1:2]) %*% t( matrixYML[1:5,1:2]) # equation 8 
  y <- cbind(as.numeric(cov.sample[1,1:5]))
  sigma.zw <- cov.zy[,1:5] %*% t(G) # equation 10
  sigma.ww <- G %*% cov.yy[,1:5] %*% t(G) # equation 11
  w <- G %*% y
  
  AA <- sigma.zw %*% solve(sigma.ww) %*% w
  
  
  
  
  
  
  pred.list2[[stocks]] <- as.numeric((AA + c(mean(normailzed.matrix[[6]]),
                                             mean(normailzed.matrix[[7]]),
                                             mean(normailzed.matrix[[8]]),
                                             mean(normailzed.matrix[[9]]),
                                             mean(normailzed.matrix[[10]])))*matrix.data[1,10])
  
  actual.list2[[stocks]] <- as.numeric(matrix.data[1,6:10]) 
  
  MMSE.list2[[stocks]] <- (as.numeric(matrix.data[1,6:10]) - as.numeric((AA + c(mean(normailzed.matrix[[6]]),
                                                                                mean(normailzed.matrix[[7]]),
                                                                                mean(normailzed.matrix[[8]]),
                                                                                mean(normailzed.matrix[[9]]),
                                                                                mean(normailzed.matrix[[10]])))*matrix.data[1,10]))^2
  
  
  
}


###############
##############33
################3




  
  matrix.data <- get(paste0(stocks, sep = ".", "Hankel.Matrix"))[-1]
  
  q <- ncol(matrix.data) # because Q=M
  
  matrix.data <- get(paste0(stocks, sep = ".", "Hankel.Matrix"))[-1]
  
  q <- ncol(matrix.data) # because Q=M
  
  normailzed.matrix <- data.frame(matrix(nrow =nrow(matrix.data), ncol = q)) # generate empty dataframe that will hold all normalized values
  
  for (i in 1:nrow(matrix.data)) {
    for (j in 1:ncol(matrix.data)) {
      normailzed.matrix[i,j] <- matrix.data[i,j]/matrix.data[i,q]
    }
  }
  
  normailzed.centered.matrix <- data.frame(matrix(nrow = nrow(matrix.data), ncol = ncol(matrix.data)))
  
  for (i in 1:ncol(normailzed.matrix)) {
    normailzed.centered.matrix[i] <- normailzed.matrix[i] - mean(normailzed.matrix[[i]])
  }
  
  
  normailzed.centered.matrix <- normailzed.centered.matrix[-q]
  
  matrix.X <- as.matrix(normailzed.centered.matrix)
  
  
  gamma.diag <- gamma * diag(nrow(matrix.data))
  gamma.const <- (1-0.98)/(1-0.98^(nrow(matrix.data)-1))
  cov.sample <- gamma.const * (t(matrix.X) %*% gamma.diag %*% matrix.X)
  
  matrix.y <- cov.sample[,1:(ncol(cov.sample)/2)]
  matrix.z <- cov.sample[,-(1:ncol(cov.sample)/2)]
  
  cov.yy <- cov(matrix.y)
  cov.zz <- cov(matrix.z)
  
  cov.yz <- cov(matrix.y,matrix.z)
  cov.zy <- cov(matrix.z, matrix.y)
  
  
  
  eigenvectors.yy <- eigen(cov.yy)
  
  
  matrixYML <- eigenvectors.yy$vectors # this is our V matrix 
  
  matrixYML <- as.matrix(matrixYML)
  
  
  G <- solve(t(matrixYML[1:M,1:L])  %*%  matrixYML[1:M,1:L]) %*% t( matrixYML[1:M,1:L]) # equation 8 
  y <- cbind(as.numeric(cov.sample[1,1:M]))
  sigma.zw <- cov.zy[,1:M] %*% t(G) # equation 10
  sigma.ww <- G %*% cov.yy[,1:M] %*% t(G) # equation 11
  w <- G %*% y
  
  AA <- sigma.zw %*% solve(sigma.ww) %*% w
  
  
  
  
  
  
  pred.list2[[stocks]] <- as.numeric((AA + as.numeric(sapply(normailzed.matrix[(M+1):(N)], mean)))*matrix.data[1,N])
  
  actual.list2[[stocks]] <- as.numeric(matrix.data[1,(M+1):(N)]) 
  
  MMSE.list2[[stocks]] <- (as.numeric(matrix.data[1,(M+1):(N)]) - as.numeric((AA + as.numeric(sapply(normailzed.matrix[(M+1):(N)],
                                                                                                     mean)))*matrix.data[1,N]))^2
  
  print(paste0(stocks, sep = " ", "DONE"))


