library(xts)
library(quantmod)
library(PerformanceAnalytics)
library(ggplot2)
library(tidyverse)


# Test

# set begin-end date and stock namelist
begin_date <- "2014-01-01"
end_date <- "2019-12-31"

data_part <- 0.75


indexes<- c("USDBRL=X", "^BVSP", "^IXIC", "^DJI", "^GSPC")

getSymbols(indexes, src = "yahoo", from = begin_date, to = end_date) # S&P 500


stock_namelist <- c("NATU3.SA", "MRVE3.SA", "VALE3.SA",
                    "BTOW3.SA", "BRKM5.SA", "BRFS3.SA", "BBDC3.SA",
                    "BBSE3.SA", "BRAP4.SA", "BRML3.SA", "CIEL3.SA",
                    "CSAN3.SA", "CVCB3.SA", "ECOR3.SA", "EGIE3.SA",
                    "FLRY3.SA", "HYPE3.SA", "JBSS3.SA", "LAME4.SA",
                    "LREN3.SA", "TIMP3.SA", "VVAR3.SA",
                    "PETR4.SA", "ABEV3.SA", "B3SA3.SA", "BBAS3.SA",
                    "BBDC4.SA", "CMIG4.SA", "CSNA3.SA", "CYRE3.SA",
                    "EMBR3.SA", "GGBR3.SA", "GGBR4.SA", "GOAU4.SA",
                    "ITSA4.SA", "ITUB4.SA", "PCAR4.SA", "USIM5.SA")

# download data from YahooFinance
prices <- xts()

for (stock_index in 1:length(stock_namelist))
  prices <- cbind(prices, Ad(
    getSymbols(
      stock_namelist[stock_index],
      from = begin_date,
      to = end_date,
      auto.assign = FALSE)
  )
)

prices <- prices %>% na.omit()


colnames(prices) <- stock_namelist
indexClass(prices) <- "Date"


# compute log-returns and linear returns
X_log <- diff(log(prices))[-1]
X_log[is.na(X_log)] <- 0

lag(prices)

X_lin <- (prices / lag(prices) - 1)[-1]
X_lin[is.na(X_lin)] <- 0

N <- ncol(X_log)  # number of stocks
T <- nrow(X_log)  # number of days


##### Plot of Normalized Prices #####
plot(
  prices / rep(prices[1,], each = nrow(prices)),
  col = rainbow10equal,
  legend.loc = "topleft",
  main = "Normalized prices"
)


# split data into training and set data
T_trn <- round(data_part * T)

X_log_trn <- X_log[1:T_trn,]
X_lin_trn <- X_lin[1:T_trn,]

X_log_tst <- X_log[(T_trn + 1):T,]
X_lin_tst <- X_lin[(T_trn + 1):T,]


mu <- colMeans(X_log_trn)

Sigma <- cov(X_log_trn)


###################### MARKOWITZ PORTFOLIO OPTIMIZATION ######################
library(CVXR)

portolioMarkowitz <- function(mu, Sigma, lmd = 0.5) {
  w <- Variable(nrow(Sigma))
  
  prob <- Problem(Maximize(t(mu) %*% w - lmd * quad_form(w, Sigma)),
                  constraints = list(w >= 0, w <= 0.1, sum(w) == 1))
  
  result <- solve(prob)
  
  return(as.vector(result$getValue(w)))
}

portolioGMVP <- function(Sigma) {
  w <- Variable(nrow(Sigma))
  
  prob <- Problem(Minimize(quad_form(w, Sigma)),
                  constraints = list(w >= 0, w <= 0.1, sum(w) == 1))
  
  result <- solve(prob)
  
  return(as.vector(result$getValue(w)))
}

w_Markowitz <- (portolioMarkowitz(mu, Sigma))
w_GMVP <- (portolioGMVP(Sigma))




# combine portfolios
w_all <- cbind("GMVP" = w_GMVP,
               "Markowitz" = w_Markowitz)

rownames(w_all) <- colnames(X_lin)

# compute returns of all portfolios
ret_all <- xts(X_lin %*% w_all, index(X_lin))

ret_all_trn <- ret_all[1:T_trn,]
ret_all_tst <- ret_all[-c(1:T_trn),]



###################### DOWNSIDE RISK PORTFOLIO OPTIMIZATION ######################
library(CVXR)

portfolioDR <- function(X, lmd = 0.5, alpha = 2) {
  T <- nrow(X)
  N <- ncol(X)
  X <- as.matrix(X)
  mu <- colMeans(X)
  w <- Variable(N)
  
  prob <- Problem(Maximize(t(w) %*% mu - (lmd / T) * 
                             sum(pos(t(mu) %*% w - X %*% w)) ^ alpha),
                  constraints = list(w >= 0, w <= 0.1, sum(w) == 1))
  
  result <- solve(prob)
  
  return(as.vector(result$getValue(w)))
}

w_DR_alpha1 <- portfolioDR(X_log_trn, alpha = 1)
w_DR_alpha2 <- portfolioDR(X_log_trn, alpha = 2)
w_DR_alpha3 <- portfolioDR(X_log_trn, alpha = 3)

# combine portfolios
w_all <- cbind(w_all,
               "DR-alpha-1" = w_DR_alpha1,
               "DR-alpha-2" = w_DR_alpha2,
               "DR-alpha-3" = w_DR_alpha3)

# compute returns of all portfolios
ret_all <- xts(X_lin %*% w_all, index(X_lin))
ret_all_trn <- ret_all[1:T_trn,]
ret_all_tst <- ret_all[-c(1:T_trn),]





###################### CVaR PORTFOLIO OPTIMIZATION ######################

portolioCVaR <- function(X, lmd = 0.5, alpha = 0.95) {
  T <- nrow(X)
  N <- ncol(X)
  X <- as.matrix(X)
  mu <- colMeans(X)
  
  # variables
  w <- Variable(N)
  z <- Variable(T)
  zeta <- Variable(1)
  # problem
  prob <-
    Problem(
      Maximize(t(w) %*% mu - lmd * zeta - (lmd / (T * (
        1 - alpha
      ))) * sum(z)),
      constraints = list(z >= 0, z >= -X %*% w - zeta,
                         w >= 0, w <= 0.1, sum(w) == 1)
    )
  
  result <- solve(prob)
  
  return(as.vector(result$getValue(w)))
}

w_CVaR095 <- portolioCVaR(X_log_trn, alpha = 0.95)
w_CVaR099 <- portolioCVaR(X_log_trn, alpha = 0.99)

# combine portfolios
w_all <- cbind(w_all,
               "CVaR-alpha-0.95" = w_CVaR095,
               "CVaR-alpha-0.99" = w_CVaR099)

# compute returns of all portfolios
ret_all <- xts(X_lin %*% w_all, index(X_lin))
ret_all_trn <- ret_all[1:T_trn,]
ret_all_tst <- ret_all[-c(1:T_trn),]



###################### MAX DRAWDOWN PORTFOLIO OPTIMIZATION ######################

portfolioMaxDD <- function(X, c = 0.2) {
  T <- nrow(X)
  N <- ncol(X)
  X <- as.matrix(X)
  X_cum <- apply(X, MARGIN = 2, FUN = cumsum)
  mu <- colMeans(X)
  # variables
  w <- Variable(N)
  u <- Variable(T)
  # problem
  prob <- Problem(
    Maximize(t(w) %*% mu),
    constraints = list(w >= 0, w <= 0.1, sum(w) == 1,
                       u <= X_cum %*% w + c,
                       u >= X_cum %*% w,
                       u[-1] >= u[-T])
  )
  result <- solve(prob)
  return(as.vector(result$getValue(w)))
}

w_MaxDD_c018 <- portfolioMaxDD(X_log_trn, c = 0.18)
w_MaxDD_c021 <- portfolioMaxDD(X_log_trn, c = 0.21)
w_MaxDD_c024 <- portfolioMaxDD(X_log_trn, c = 0.24)


# combine portfolios
w_all <- cbind(
  w_all,
  "Max-DD-c-018" = w_MaxDD_c018,
  "Max-DD-c-021" = w_MaxDD_c021,
  "Max-DD-c-024" = w_MaxDD_c024
)

# compute returns of all portfolios
ret_all <- xts(X_lin %*% w_all, index(X_lin))
ret_all_trn <- ret_all[1:T_trn,]
ret_all_tst <- ret_all[-c(1:T_trn),]





###################### AVERAGE DRAWDOWN PORTFOLIO OPTIMIZATION ######################

portfolioAveDD <- function(X, c = 0.2) {
  T <- nrow(X)
  N <- ncol(X)
  X <- as.matrix(X)
  X_cum <- apply(X, MARGIN = 2, FUN = cumsum)
  mu <- colMeans(X)
  # variables
  w <- Variable(N)
  u <- Variable(T)
  # problem
  prob <- Problem(Maximize(t(w) %*% mu),
                  constraints = list(w >= 0, w <= 0.1,sum(w) == 1,
                                     mean(u) <= mean(X_cum %*% w) + c,
                                     u >= X_cum %*% w,
                                     u[-1] >= u[-T]))
  
  result <- solve(prob)
  
  return(as.vector(result$getValue(w)))
}

w_AveDD_c004 <- portfolioAveDD(X_log_trn, c = 0.04)
w_AveDD_c006 <- portfolioAveDD(X_log_trn, c = 0.06)
w_AveDD_c008 <- portfolioAveDD(X_log_trn, c = 0.08)


# combine portfolios
w_all <- cbind(
  w_all,
  "Ave-DD-c-004" = w_AveDD_c004,
  "Ave-DD-c-006" = w_AveDD_c006,
  "Ave-DD-c-008" = w_AveDD_c008
)

# compute returns of all portfolios
ret_all <- xts(X_lin %*% w_all, index(X_lin))
ret_all_trn <- ret_all[1:T_trn,]
ret_all_tst <- ret_all[-c(1:T_trn),]






###################### CDaR PORTFOLIO OPTIMIZATION ######################

portfolioCDaR <- function(X, c = 0.1, alpha = 0.95) {
  T <- nrow(X)
  N <- ncol(X)
  X <- as.matrix(X)
  X_cum <- apply(X, MARGIN = 2, FUN = cumsum)
  mu <- colMeans(X)
  # variables
  w <- Variable(N)
  z <- Variable(T)
  zeta <- Variable(1)
  u <- Variable(T)
  # problem
  prob <- Problem(
    Maximize(t(w) %*% mu),
    constraints = list(
      w >= 0, w <= 0.1,
      sum(w) == 1,
      zeta + (1 / (T * (1 - alpha))) * sum(z) <= c,
      z >= 0,
      z >= u - X_cum %*% w - zeta,
      u >= X_cum %*% w,
      u[-1] >= u[-T]
    )
  )
  result <- solve(prob)
  
  return(as.vector(result$getValue(w)))
}

w_CDaR095_c014 <- portfolioCDaR(X_log_trn, c = 0.14, alpha = 0.95)
w_CDaR099_c016 <- portfolioCDaR(X_log_trn, c = 0.16, alpha = 0.99)
w_CDaR095_c016 <- portfolioCDaR(X_log_trn, c = 0.16, alpha = 0.95)
w_CDaR099_c018 <- portfolioCDaR(X_log_trn, c = 0.18, alpha = 0.99)


# combine portfolios
w_all <- cbind(w_all,
               "CDaR095-c-014" = w_CDaR095_c014,
               "CDaR099-c-016" = w_CDaR099_c016,
               "CDaR095-c-016" = w_CDaR095_c016,
               "CDaR099-c-018" = w_CDaR099_c018)

# compute returns of all portfolios
ret_all <- xts(X_lin %*% w_all, index(X_lin))
ret_all_trn <- ret_all[1:T_trn,]
ret_all_tst <- ret_all[-c(1:T_trn),]

# performance
t(table.AnnualizedReturns(ret_all_trn))

t(table.AnnualizedReturns(ret_all_tst))




port_to_be <- w_AveDD_c008

test <- t(as.list(round(port_to_be, 4)))
colnames(test) <- stock_namelist

df <- data.frame(t((test)))
df <- subset(df, t..test.. >0)

assets <- rownames(df)


ggplot(data = df, aes(x = assets, y = as.numeric(df$t..test..) , fill = assets)) +
  geom_bar(stat = "identity",
           position = position_dodge(),
           colour = "black") +
  geom_text(
    aes(label = sprintf("%.02f %%", as.numeric(df$t..test..) * 100)),
    position = position_dodge(width = 0.9),
    vjust = -0.25, check_overlap = TRUE) +
  ggtitle("Max Drawdown Optimal Weights - C 0.24") + theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Assets", y = "Weight (%)")








###### cumulative PnL over time ######

BVSP_ret<- (BVSP$BVSP.Close / lag(BVSP$BVSP.Close) - 1)[-1]

final_perf <- cbind(BVSP_ret, ret_all$`CDaR095-c-016`)


{
  chart.CumReturns(final_perf,
                   main = "Performance of different portfolios",
                   wealth.index = TRUE, legend.loc = "topleft",
                   colorset = rich8equal)
  addEventLines(xts("training", index(X_lin[T_trn])),
                    srt = 90, pos = 2, lwd = 2,
                col = "darkblue")
}



final_perf <- final_perf[index(final_perf) >= '2019-01-01' ]

{
  chart.CumReturns(final_perf,
                   main = "Performance of different portfolios",
                   wealth.index = TRUE, legend.loc = "topleft",
                   colorset = rich8equal)
  addEventLines(xts("training", index(X_lin[T_trn])),
                srt = 90, pos = 2, lwd = 2,
                col = "darkblue")
}



chart.Drawdown(final_perf, main = "Drawdown of portfolios", 
               legend.loc = "bottomleft", colorset = rich8equal)



test <- t(as.list(round(port_to_be, 4) * 2000))
colnames(test) <- stock_namelist

df <- data.frame(t((test)))
df <- subset(df, t..test.. >0)

df
