setwd("D:/rcode/treehouse")
AAPL = read.csv("AAPL.csv",header = TRUE)
IBM = read.csv("IBM.csv",header = TRUE)
GOOG = read.csv("GOOG.csv",header = TRUE)
BP = read.csv("BP.csv",header = TRUE)
XOM = read.csv("XOM.csv",header = TRUE)
COST = read.csv("COST.csv",header = TRUE)
GS = read.csv("GS.csv",header = TRUE)
AdjPrices = cbind(AAPL[,"Adj.Close"],IBM[,"Adj.Close"],GOOG[,"Adj.Close"],
          BP[,"Adj.Close"],XOM[,"Adj.Close"],COST[,"Adj.Close"],
          GS[,"Adj.Close"])
r = log(AdjPrices[2:nrow(AdjPrices),]/AdjPrices[1:(nrow(AdjPrices)-1),])
w = c(0.15,0.2,0.2,0.15,0.1,0.15,0.05)
r_w = apply(r*w,1,sum)
r_w_asc = sort(r_w)
#(a) Use Historical Simulation
VaR95 = r_w_asc[round(0.05*length(r_w_asc))]
CVaR95 = mean(r_w_asc[1:round(0.05*length(r_w_asc))])
#(b)
#density(r_w)
ExMean = mean(r_w) #Expected mean
std = sd(r_w) #Standard deviation
#Plot the density distribution
library(ggplot2)
p<-ggplot(data.frame(r_w), aes(x = r_w))
p + geom_density(color = 'black', fill = 'gray')
  #+geom_vline(aes(xintercept = 0.005, color='red'),linetype='dashed')
#cor(r)
#Assume t-distribution
df = 251-2-1
lambda = std/sqrt(df/(df-2))
t_VaR95 = ExMean + lambda*qt(0.05,df) #non-centrality delta is omitted
f_t = function(x){
  f = gamma((df+1)/2)/gamma(df/2)*sqrt(lambda/pi/df)*(1+lambda*(x-ExMean)^2/2)^(-(df+1)/2)
  return(f)
}
t_CVaR95 = ExMean + lambda*(dt(qt(0.05,df),df)/0.05*(df+(qt(0.05,df))^2)/(df-1)) 
#t_CVaR95 = ExMean + lambda*(f_t(qt(0.05,df))/0.05*((df+(qt(0.05,df))^2)/(df-1)))

#(d)
#install.packages("lpSolve")
library(lpSolve)
dim1=AAPL[2:252,"Date"]
dim2=c('AAPL','IBM','GOOG','BP','XOM','COST','GS')
WEIGHTS=array(dim=c(251,7),dimnames = list(dim1,dim2))
OptimalReturn = array(dim=c(251,1))
for(i in 1:251){
  f.obj <- r[i,]
  f.con <- matrix(c(rep(1,7)),nrow=1)
  f.dir <- c ( "=")
  f.rhs <- c (1)
  lp.result <- lp ( "max" , f.obj, f.con, f.dir, f.rhs)
  #lp.result
  WEIGHTS[i,] = lp.result$solution
  OptimalReturn[i] = lp.result$objval
}
#Success: the objective function on 12/31/2016 is 0.0345337