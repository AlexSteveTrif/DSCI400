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




# test 
stock_picks <- Stocks$Symbol[1:5]

stock_picks
# weekly pull function (no excel sheet)






# N is our 2x the number of trading days we wish to predict 

# L is the number of principle components we want to consider 

# M is the length of predction period and is equal to N/2




Stocks <- read_excel("Stocks.xlsx")
stock_picks <- Stocks$Symbol




func <- function(N, M, L) {
  
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
                        1:N)
  
  
  pred.list <- list()
  MMSE.list <- list()
  actual.list <- list()
  d.stat.list <- list()
  
  pred.list2 <- list()
  MMSE.list2 <- list()
  actual.list2 <- list()
  d.stat.list2 <- list()
  
  
  
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
    
    
    G <- solve(t(matrixYML[1:M,1:L])  %*%  matrixYML[1:M,1:L]) %*% t( matrixYML[1:M,1:L]) # equation 8 
    y <- cbind(as.numeric(matrix.data[1,1:M]))
    sigma.zw <- cov.zy[,1:M] %*% t(G) # equation 10
    sigma.ww <- G %*% cov.yy[,1:M] %*% t(G) # equation 11
    w <- G %*% y
    
    
    
    pred.list[[stocks]] <- as.numeric(sigma.zw %*% solve(sigma.ww) %*% w)
    
    actual.list[[stocks]] <- as.numeric(matrix.data[1,(M+1):(N+1)])
    
    MMSE.list[[stocks]] <- (as.numeric(matrix.data[1,(M+1):(N+1)]) - as.numeric(sigma.zw %*% solve(sigma.ww) %*% w))^2
    
    d.stat.list[[stocks]] <- mean(ifelse(((matrix.data[1,(M+2):(N+1)] - matrix.data[1,(M+1):(N)])*(matrix.data[1,(M+2):(N+1)] - as.numeric(sigma.zw %*% solve(sigma.ww) %*% w)))>0, 1,0))

    
    
    
  }
  
  
  
  
  
  gamma <- 0.98^(1:dim(matrix.data)[1])
  for (i in 1:length(gamma)) {
    if(gamma[i] < 10^-3){
      gamma[i] <- 0
    }
  }
  
  for (stocks in stock_picks) {
    
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
    
    d.stat.list2[[stocks]] <- mean(ifelse(((matrix.data[1,(M+2):(N+1)] - matrix.data[1,(M+1):(N)])*(matrix.data[1,(M+2):(N+1)] - as.numeric((AA + as.numeric(sapply(normailzed.matrix[(M+1):(N)], mean)))*matrix.data[1,N])))>0, 1,0))

    
    print(paste0(stocks, sep = " ", "DONE"))
  }
  
  
  Pred.Tables <- list()
  Pred.Tables2 <- list()
  
  
  for (stocks in stock_picks) {
    
    Pred.data <- get(paste0(stocks, sep = ".", "Hankel.Matrix"))
    
    pred.table <- data.frame(date = as.Date(Pred.data[(M+1):(N+1),1]),
                             actual = as.numeric(Pred.data[(M+1):(N+1),2]), pred = pred.list[[paste0(stocks)]])
    
    # assign(paste0(stocks, sep = ".", "pred.table"), pred.table, envir = .GlobalEnv)
    
    Pred.Tables[[stocks]] <- pred.table
    
    
    
    
    
    
    
  }
  
  
  for (stocks in stock_picks) {
    
    Pred.data <- get(paste0(stocks, sep = ".", "Hankel.Matrix"))
    
    pred.table <- data.frame(date = as.Date(Pred.data[(M+1):(N),1]),
                             actual = as.numeric(Pred.data[(M+1):(N),2]), pred = pred.list2[[paste0(stocks)]])
    
    # assign(paste0(stocks, sep = ".", "pred.table2"), pred.table, envir = .GlobalEnv)
    
    Pred.Tables2[[stocks]] <- pred.table
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  assign("pred.list", pred.list, envir = .GlobalEnv)
  assign("actual.list", actual.list, envir = .GlobalEnv)
  assign("MSE.list", MMSE.list, envir = .GlobalEnv)
  assign("d.stat.list", d.stat.list, envir = .GlobalEnv)
  
  assign("pred.list2", pred.list2, envir = .GlobalEnv)
  assign("actual.list2", actual.list2, envir = .GlobalEnv)
  assign("MSE.list2", MMSE.list2, envir = .GlobalEnv)
  assign("d.stat.list2", d.stat.list2, envir = .GlobalEnv)
  
  assign("Pred.Tables", Pred.Tables, envir = .GlobalEnv)
  assign("Pred.Tables2", Pred.Tables2, envir = .GlobalEnv)
  
  
  
  for (stocks in stock_picks) {
    rm(list = paste0(stocks, sep = ".", "Hankel.Matrix"), envir = .GlobalEnv)
  }
  
}


func(10,5,2)


func(100, 50, 2)


func(30, 15, 2)







plot.gen <- function(stocks){
  
  
  pred.table <- Pred.Tables[[stocks]]
  
  ggplot(data = pred.table) + 
    geom_line(aes(x = date, y = actual), color = "maroon") + 
    geom_line(aes(x = date, y = pred), color = "goldenrod") + 
    xlab("Date") +  ylab("Stock Price")
  
}

plot.gen2 <- function(stocks){
  
  
  pred.table <- Pred.Tables2[[stocks]] 
  
  ggplot(data = pred.table) + 
    geom_line(aes(x = date, y = actual), color = "maroon") + 
    geom_line(aes(x = date, y = pred), color = "goldenrod") + 
    xlab("Date") +  ylab("Stock Price")
  
}

# gold is our prediction, red is our actual 



plot.gen("MSFT")
plot.gen2("MSFT")



plot.gen("AAPL")
plot.gen2("AAPL")


plot.gen("GOOG")
plot.gen2("GOOG")



 
















