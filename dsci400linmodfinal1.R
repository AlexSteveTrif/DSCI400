


lin.pred.list <- list()
lin.MMSE.list <- list()
lin.list <- list()


for (stock in stock_picks) {
  
  
  matrix.data <- get(paste0(stock, sep = ".", "Hankel.Matrix"))[-1]
  
  linear.model <- lm(day10 ~ day0 + day1 + day2 + day3 + day4, data = matrix.data)
  
  lin.list[[stock]] <- summary(linear.model)
  
  
  for (i in 1:nrow(matrix.data)) {
  
  lin.pred.list[[stock]][[i]] <- predict(linear.model, newdata = matrix.data[i,1:5])
  
  lin.MMSE.list[[stock]][[i]] <- (matrix.data[i,11] - predict(linear.model, newdata = matrix.data[i,1:5]))^2
  
  }
  
}



#################3333
#################3333
#################3333
#################3333#################3333
#################3333
#################3333
#################3333#################3333
#################3333
#################3333
#################3333#################3333
#################3333
#################3333
#################3333#################3333
#################3333
#################3333
#################3333#################3333
#################3333
#################3333
#################3333#################3333
#################3333
#################3333
#################3333


lin.pred.list <- list()
lin.MMSE.list <- list()
lin.list <- list()


for (stock in 1:length(stock_picks)) {
  matrix.data <- get(paste0(stock_picks[stock], sep = ".", "Hankel.Matrix"))[-1]
  
  linear.model <- lm(day10 + day9 ~ day0 + day1 + day2 + day3 + day4, data = matrix.data)
  
  lin.list[[stock]] <- summary(linear.model)
  
  lin.pred.list[[stock]] <- predict(linear.model, newdata = matrix.data[1,1:5])
  
  lin.MMSE.list[[stock]] <- (matrix.data[1,11] - predict(linear.model, newdata = matrix.data[1,1:5]))^2
  
}

#################3333
#################3333
#################3333
#################3333


