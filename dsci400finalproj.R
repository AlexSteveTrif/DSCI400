require(tidyverse)
require(quantmod)
require(tidyquant)
require(ggplot2)
require(dplyr)
require(readxl)

# desktop wd 

setwd("C:/Users/Alex/OneDrive/Documents/oc/2nd semester/DSCI400/final project")

# laptop wd 

setwd("C:/Users/alext/OneDrive/Documents/oc/2nd semester/DSCI400/final project")



Stocks <- read_excel("Stocks.xlsx")
stock_picks <- Stocks$Symbol # we can add stocks to this excel sheet! 


# weekly pull function (no excel sheet)



get_stocks_Daily_Pull <- function(symbols, start_date, end_date, metric, ...){
  lags <- c(...)
  stock_weeklyClose_df_list <- list()
  stock_weekly_Return_list <- list()
  weekly_close_lag <- list()
  for (symbol in symbols) {
    stock <- getSymbols.yahoo(symbol, src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)
    stock_weekly <- to.daily(stock, indexAt = "last", OHLC = FALSE)
    colnames(stock_weekly) <- paste0(symbol, sep = ".", c("Open", "High", "Low", "Close", "Volume", "Adjusted"))
    stock_weeklyClose <- stock_weekly[, paste0(symbol, sep = ".", metric)]
    stock_weeklyClose_df <- data.frame(Close = coredata(stock_weeklyClose))
    
    stock_weeklyClose_df_list[[symbol]] <- stock_weeklyClose_df
    
  }
  stocks_weekly_close <- do.call(cbind, sapply(stock_weeklyClose_df_list, "[",1))
  stocks_weekly_close <- data.frame(stocks_weekly_close)
  colnames(stocks_weekly_close) <- symbols
  # adding date columns in front
  date <- getSymbols.yahoo(symbols[1], src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)
  getdate <- to.daily(date, indexAt = "last", OHLC = FALSE)
  colnames(getdate) <- paste0(symbols[1], sep = ".", c("Open", "High", "Low", "Close", "Volume", "Adjusted"))
  get.date <- getdate[, paste0(symbols[1], ".Close")]
  dates <- data.frame(Date = index(get.date))
  dates$Date <- as.character(dates$Date)

  stocks_weekly_close <- cbind(dates, stocks_weekly_close)
  assign(paste0(metric, sep = ".", "daily"), stocks_weekly_close, envir = .GlobalEnv)
  stocks_weekly_close <- arrange(stocks_weekly_close, desc(Date))
  for (symbol in symbols) {
    stk <- stocks_weekly_close %>% select(c(symbol))
    stk_lag <- data.frame(stk)
    for (i in lags){
      stk_lag <- cbind(stk_lag, lag(stk, n = i))
    }
    stk_lag <- cbind(stocks_weekly_close$Date,stk_lag)
    colnames(stk_lag) <- c("date", "day0", paste0("day",lags))
    stk_lag <- na.omit(stk_lag)
    stk_lag <- arrange(stk_lag, date)
    weekly_close_lag[[symbol]] <- stk_lag
    
    assign(paste0(symbol, sep = ".", "Hankel.Matrix"), stk_lag, envir = .GlobalEnv)
    # i want to tag on weekly returns to the lag data to generate a covarience matrix 
    
    
  }
  names(weekly_close_lag) <- symbols
  
  for (i in 1:length(weekly_close_lag)) {
    write.csv(weekly_close_lag[[i]], paste0(metric, sep = ".", "Hankel.Matrix" ,names(weekly_close_lag[i]), ".csv"), row.names = FALSE)}
  

  rm(list = paste0(metric, sep = ".", "daily"), envir = .GlobalEnv)
  

}


get_stocks_Daily_Pull(stock_picks,  "2020-12-01",
                      "2023-02-10", "Adjusted",
                      1:10)

1:10

# linear regression model 
##################################

matrix.data <- AAPL.Hankel.Matrix[-1]

matrix.data

linear.model <- lm(day10  ~ day0 + day1 + day2 + day3 + day4, data = matrix.data)

summ <- summary(linear.model)

summ$r.squared
summ$coefficients


predict(linear.model, newdata = matrix.data[1,1:5])




as.numeric(matrix.data[1,1:5])



matrix.data[1,11]


####################33333


matrix.data <- AAPL.Hankel.Matrix[-1]

require(GGally)

ggpairs(matrix.data)



# decreasing exponential weights 

dim(matrix.data)[1] # this is our K 

gamma <- 0.98^(1:dim(matrix.data)[1])
for (i in 1:length(gamma)) {
  if(gamma[i] < 10^-3){
    gamma[i] <- 0
  }
}

gamma


matrix.y <- matrix.data[,1:(length(matrix.data)/2)]
matrix.z <- matrix.data[,-(1:length(matrix.data)/2)]

cov.yy <- cov(matrix.y)
cov.zz <- cov(matrix.z)

cov.yz <- cov(matrix.y,matrix.z)
cov.zy <- cov(matrix.z, matrix.y)

cov.mat.data <- cov(matrix.data) # xx


# eigenvalues are calculated from the covariance matrix by solving the characteristic equation

# solution to the characteristic equation yields the eigenvalues for the covariance matrix

eigenvectors.yy <- eigen(cov.yy)

eigenvectors.yy$values #eigenvalues

eigenvectors.yy$vectors # eigenvectors


eigenvalues.yy <- diag(5) * eigenvectors.yy$values



matrixYML <- eigenvectors.yy$vectors # this is our V matrix 



matrixYML <- as.matrix(matrixYML)

matrixYML
 

round(matrixYML %*% t(matrixYML))     # test to see if this matrix is unitary

# this is our identity matrix, so this matrix is unitary 


matrixYML %*% eigenvalues.yy %*% t(Conj(matrixYML)) # equation 4 

as.matrix(cov.yy)


solve(diag(5)) # note that the solution to an identity matrix is the identity matrix  



# so the equation to G (8) should just be equal to matrixYML (AKA as eigenvectors)

solve(round(t(matrixYML)%*%  matrixYML ,1)) %*% matrixYML == matrixYML # equation 8 


# now lets consider the first L column (first 2 principle components)

matrixYML[,1]

matrixYML[,2]

# note their dot product 

matrixYML[,1] %*% matrixYML[,2]

# its zero! These vectors are orthogonal and form a basis 




# now lets consider the first M  rows 

G <- solve(t(matrixYML[1:5,1:2])  %*%  matrixYML[1:5,1:2]) %*% t( matrixYML[1:5,1:2]) # equation 8 

G

y <- cbind(as.numeric(matrix.data[1,1:5]))

dim(y)


sigma.zw <- cov.zy[,1:5] %*% t(G) # equation 10

sigma.ww <- G %*% cov.yy[,1:5] %*% t(G) # equation 11 

w <- G %*% y

w

matrix.data[1,5:11] #is our estimator 

sigma.zw %*% solve(sigma.ww) %*% w # this is our estimator 

as.numeric(matrix.data[1,5:11]) - as.numeric(sigma.zw %*% solve(sigma.ww) %*% w)




matrix.data[1,6:11] - sigma.zw %*% solve(sigma.ww) %*% w


dim(sigma.zw)
dim(sigma.ww)
dim(w)


AA <- sigma.zw %*% solve(sigma.ww) %*% w # this is our estimator 

dim(sigma.zw)
dim(sigma.ww)
dim(w)

AA







sum((matrix.data[1,6:11] - AA)^2) # sum of squared errors from using covarience matrix without exponential weighting 



################################
################################
################################
################################

################################  Normalize Each Row???

matrix.data



matrix.data <- AAPL.Hankel.Matrix[-1]
q <- ncol(matrix.data) # because Q=M
q

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


hist(normailzed.centered.matrix$X3, breaks = 25)




matrix.y <- normailzed.centered.matrix[,1:(length(matrix.data)/2)]
matrix.z <- normailzed.centered.matrix[,-(1:length(matrix.data)/2)]

cov.yy <- cov(matrix.y)
cov.zz <- cov(matrix.z)

cov.yz <- cov(matrix.y,matrix.z)
cov.zy <- cov(matrix.z, matrix.y)

cov.mat.data <- cov(matrix.data) # xx


# eigenvalues are calculated from the covariance matrix by solving the characteristic equation

# solution to the characteristic equation yields the eigenvalues for the covariance matrix

eigenvectors.yy <- eigen(cov.yy)

eigenvectors.yy$values #eigenvalues

eigenvectors.yy$vectors # eigenvectors


eigenvalues.yy <- diag(5) * eigenvectors.yy$values



matrixYML <- eigenvectors.yy$vectors # this is our V matrix 



matrixYML <- as.matrix(matrixYML)

matrixYML


round(matrixYML %*% t(matrixYML))     # test to see if this matrix is unitary

# this is our identity matrix, so this matrix is unitary 


matrixYML %*% eigenvalues.yy %*% t(Conj(matrixYML)) # equation 4 

as.matrix(cov.yy)


solve(diag(5)) # note that the solution to an identity matrix is the identity matrix  



# so the equation to G (8) should just be equal to matrixYML (AKA as eigenvectors)

solve(round(t(matrixYML)%*%  matrixYML ,1)) %*% matrixYML == matrixYML # equation 8 


# now lets consider the first L column (first 2 principle components)

matrixYML[,1]

matrixYML[,2]

# note their dot product 

matrixYML[,1] %*% matrixYML[,2]

# its zero! These vectors are orthogonal and form a basis 




# now lets consider the first M  rows and first L columns 


G <- solve(t(matrixYML[1:5,1:2])  %*%  matrixYML[1:5,1:2]) %*% t( matrixYML[1:5,1:2]) # equation 8 

G

dim(G)[2] # y vector should have same number of rows

y <- cbind(as.numeric(matrix.data[1,1:5]))

dim(y)




sigma.zw <- cov.zy[,1:5] %*% t(G) # equation 10

sigma.ww <- G %*% cov.yy[,1:5] %*% t(G) # equation 11 




w <- G %*% y

dim(w)




w

AA <- sigma.zw %*% solve(sigma.ww) %*% w # this is our estimator 

dim(sigma.zw)
dim(sigma.ww)
dim(w)

AA


(AA + c(mean(normailzed.matrix[[6]]),
       mean(normailzed.matrix[[7]]),
       mean(normailzed.matrix[[8]]),
       mean(normailzed.matrix[[9]]),
       mean(normailzed.matrix[[10]])))*matrix.data[1,11] # these are our predictions converted back to stock prices 




(matrix.data[1,6:11] - (AA + c(mean(normailzed.matrix[[6]]),
                             mean(normailzed.matrix[[7]]),
                             mean(normailzed.matrix[[8]]),
                             mean(normailzed.matrix[[9]]),
                             mean(normailzed.matrix[[10]])))*matrix.data[1,11])^2 # these are our predictions converted back to stock prices 


(matrix.data[1,6:11] - AA)

check scree plot 





# trying to use sample matrix 

################################################
##################################################
################################################
##################################################
################################################
##################################################
################################################
##################################################
################################################
##################################################
################################################
##################################################
################################################
##################################################


matrixYML



matrix.X <- as.matrix(normailzed.centered.matrix)

dim(matrix.X)
dim(t(matrix.X))


gamma <- 0.98^(1:dim(matrix.data)[1])
for (i in 1:length(gamma)) {
  if(gamma[i] < 10^-3){
    gamma[i] <- 0
  }
}

gamma.diag <- gamma * diag(542)

gamma.diag
dim(gamma.diag)

gamma.const <- (1-0.98)/(1-0.98^541)
gamma.const

cov.sample <- gamma.const * (t(matrix.X) %*% gamma.diag %*% matrix.X)


t(matrix.X) %*% gamma.diag %*% matrix.X


cov.sample




matrix.y <- cov.sample[,1:(length(matrix.data)/2)]
matrix.z <- cov.sample[,-(1:length(matrix.data)/2)]

cov.yy <- cov(matrix.y)
cov.zz <- cov(matrix.z)

cov.yz <- cov(matrix.y,matrix.z)
cov.zy <- cov(matrix.z, matrix.y)

cov.mat.data <- cov(matrix.data) # xx


# eigenvalues are calculated from the covariance matrix by solving the characteristic equation

# solution to the characteristic equation yields the eigenvalues for the covariance matrix

eigenvectors.yy <- eigen(cov.yy)

eigenvectors.yy$values #eigenvalues

eigenvectors.yy$vectors # eigenvectors


eigenvalues.yy <- diag(5) * eigenvectors.yy$values



matrixYML <- eigenvectors.yy$vectors # this is our V matrix 



matrixYML <- as.matrix(matrixYML)

matrixYML


round(matrixYML %*% t(matrixYML))     # test to see if this matrix is unitary

# this is our identity matrix, so this matrix is unitary 


matrixYML %*% eigenvalues.yy %*% t(Conj(matrixYML)) # equation 4 

as.matrix(cov.yy)


solve(diag(5)) # note that the solution to an identity matrix is the identity matrix  



# so the equation to G (8) should just be equal to matrixYML (AKA as eigenvectors)

solve(round(t(matrixYML)%*%  matrixYML ,1)) %*% matrixYML == matrixYML # equation 8 


# now lets consider the first L column (first 2 principle components)

matrixYML[,1]

matrixYML[,2]

# note their dot product 

matrixYML[,1] %*% matrixYML[,2]

# its zero! These vectors are orthogonal and form a basis 




# now lets consider the first M  rows 

G <- solve(t(matrixYML[1:5,1:2])  %*%  matrixYML[1:5,1:2]) %*% t( matrixYML[1:5,1:2]) # equation 8 

G

dim(G)[2] # y vector should have same number of rows

y <- cbind(as.numeric(normailzed.centered.matrix[1,1:5]))

dim(y)



sigma.zw <- cov.zy[,1:5] %*% t(G) # equation 10

sigma.ww <- G %*% cov.yy[,1:5] %*% t(G) # equation 11 




w <- G %*% y

dim(w)




w

AA <- sigma.zw %*% solve(sigma.ww) %*% w # this is our estimator 

dim(sigma.zw)
dim(sigma.ww)
dim(w)

AA


(AA + c(mean(normailzed.matrix[[6]]),
        mean(normailzed.matrix[[7]]),
        mean(normailzed.matrix[[8]]),
        mean(normailzed.matrix[[9]]),
        mean(normailzed.matrix[[10]])))*matrix.data[1,11] # these are our predictions converted back to stock prices 




(matrix.data[1,6:11] - (AA + c(mean(normailzed.matrix[[6]]),
                               mean(normailzed.matrix[[7]]),
                               mean(normailzed.matrix[[8]]),
                               mean(normailzed.matrix[[9]]),
                               mean(normailzed.matrix[[10]])))*matrix.data[1,11])^2 # these are our predictions converted back to stock prices 















