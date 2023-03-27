pred.list <- list()
MMSE.list <- list()
actual.list <- list()

# traditional covarience matric 

for (stocks in stock_picks) {
  
  matrix.data <- get(paste0(stocks, sep = ".", "Hankel.Matrix"))[-1]
  
  matrix.y <- matrix.data[,1:(length(matrix.data)/2)]
  matrix.z <- matrix.data[,-(1:length(matrix.data)/2)]
  
  cov.yy <- cov(matrix.y)
  cov.zz <- cov(matrix.z)
  
  cov.yz <- cov(matrix.y,matrix.z)
  cov.zy <- cov(matrix.z, matrix.y)
  
  
  
  eigenvectors.yy <- eigen(cov.yy)
  
  
  matrixYML <- eigenvectors.yy$vectors # this is our V matrix 
  
  matrixYML <- as.matrix(matrixYML)
  
  
  G <- solve(t(matrixYML[1:5,1:2])  %*%  matrixYML[1:5,1:2]) %*% t( matrixYML[1:5,1:2]) # equation 8 
  y <- cbind(as.numeric(matrix.data[1,1:5]))
  sigma.zw <- cov.zy[,1:5] %*% t(G) # equation 10
  sigma.ww <- G %*% cov.yy[,1:5] %*% t(G) # equation 11
  w <- G %*% y
  
  
  
  pred.list[[stocks]] <- as.numeric(sigma.zw %*% solve(sigma.ww) %*% w)
  
  actual.list[[stocks]] <- as.numeric(matrix.data[1,6:11])
  
  MMSE.list[[stocks]] <- (as.numeric(matrix.data[1,6:11]) - as.numeric(sigma.zw %*% solve(sigma.ww) %*% w))^2
  
  
  
}

##########################################
############################################
###############################################



for (stocks in stock_picks) {
  
  matrix.data <- get(paste0(stocks, sep = ".", "Hankel.Matrix"))[-1]
  
  for (i in 1:nrow(matrix.data)) {
    
  
  matrix.y <- matrix.data[,1:(length(matrix.data)/2)]
  matrix.z <- matrix.data[,-(1:length(matrix.data)/2)]
  
  cov.yy <- cov(matrix.y)
  cov.zz <- cov(matrix.z)
  
  cov.yz <- cov(matrix.y,matrix.z)
  cov.zy <- cov(matrix.z, matrix.y)
  
  
  
  eigenvectors.yy <- eigen(cov.yy)
  
  
  matrixYML <- eigenvectors.yy$vectors # this is our V matrix 
  
  matrixYML <- as.matrix(matrixYML)
  
  
  G <- solve(t(matrixYML[1:5,1:2])  %*%  matrixYML[1:5,1:2]) %*% t( matrixYML[1:5,1:2]) # equation 8 
  y <- cbind(as.numeric(matrix.data[i,1:5]))
  sigma.zw <- cov.zy[,1:5] %*% t(G) # equation 10
  sigma.ww <- G %*% cov.yy[,1:5] %*% t(G) # equation 11
  w <- G %*% y
  
  
  
  pred.list[[stocks]][[i]] <- as.numeric(sigma.zw %*% solve(sigma.ww) %*% w)
  
  actual.list[[stocks]][[i]] <- as.numeric(matrix.data[i,6:11])
  
  MMSE.list[[stocks]][[i]] <- (as.numeric(matrix.data[i,6:11]) - as.numeric(sigma.zw %*% solve(sigma.ww) %*% w))^2
  
  }
  
  
  
}

pred.list[["AAPL"]][[5]]
