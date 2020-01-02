loaded_pack <- c("readr", "bindr", "reshape2", "plyr",
                 "tidyverse", "cutpointr", "compiler", 
                 "caTools", "caret", "xgboost", "RCurl",
                 "smotefamily", "smooth", "Mcomp",
                 "quantmod", "lubridate", "XML", "furrr")

lapply(loaded_pack, FUN = function(X) {do.call("require", list(X))})


############################## DATA EXTRATION FROM YAHOO ############################## 
min_date = "2008-01-01"
max_date = Sys.Date()



opt_sample <- 0.75

indexes<- c("USDBRL=X", "^BVSP", "^IXIC", "^DJI", "^GSPC")

symb<- c("NATU3.SA", "MRVE3.SA", "VALE3.SA",
         "BTOW3.SA", "BRKM5.SA", "BRFS3.SA", "BBDC3.SA",
         "BBSE3.SA", "BRAP4.SA", "BRML3.SA", "CIEL3.SA",
         "CSAN3.SA", "CVCB3.SA", "ECOR3.SA", "EGIE3.SA",
         "FLRY3.SA", "HYPE3.SA", "JBSS3.SA", "LAME4.SA",
         "LREN3.SA", "TIMP3.SA", "VVAR3.SA",
         "PETR4.SA", "ABEV3.SA", "B3SA3.SA", "BBAS3.SA",
         "BBDC4.SA", "CMIG4.SA", "CSNA3.SA", "CYRE3.SA",
         "EMBR3.SA", "GGBR3.SA", "GGBR4.SA", "GOAU4.SA",
         "ITSA4.SA", "ITUB4.SA", "PCAR4.SA", "USIM5.SA")

getSymbols(indexes, src = "yahoo", from = min_date, to = max_date, periodicity = "monthly") # S&P 500

USDBRL <- as_tibble(data.frame(date = index(`USDBRL=X`), coredata(`USDBRL=X`)))
BVSP <- as_tibble(data.frame(date = index(BVSP), coredata(BVSP)))
NASDAQ <- as_tibble(data.frame(date = index(IXIC), coredata(IXIC)))
DOWJONES <- as_tibble(data.frame(date = index(DJI), coredata(DJI)))
S_P <- as_tibble(data.frame(date = index(GSPC), coredata(GSPC)))


getSymbols(symb, src = "yahoo", from = min_date, to = max_date, periodicity = "monthly")



#######################################################################
volat <- function(x) {
  volat <- (sd(x) / mean(x))
  return(volat)
}

return_new <- function(x, y) {
  y <- replace(y, y == 0, mean(y))
  return_new <- 1 - (mean(x) / mean(y))
  return(return_new)
}

risk_port <- function(x) {
  Data <-  as_tibble(x)
  colnames(Data)[1] <- "Open"
  colnames(Data)[4] <- "Close"
  
  descr <- summary(Data$Close)
  
  dist <- Data$Open / Data$Close
  vol <- round(volat((dist)) * 100, 2)
  
  retorno <- round(return_new(Data$Open, Data$Close) * 100, 2)
  
  Expct_Ret <- round(150000 * (return_new(Data$Open, Data$Close)), 2)
  
  VaR <-
    round(abs(return_new(Data$Open, Data$Close) - 1.282 * (volat(dist))) *
            150000, 2)
  
  rel_r_r <- round(abs(1 - (VaR / Expct_Ret)), 2)
  
  print(paste("Volatilidade do Ativo:", vol))
  print(paste("Retorno Esperado:", retorno))
  print(paste("relacao risco x retorno:", rel_r_r))
  print(paste("Repasse Esperado:", Expct_Ret * 0.6))
}


##### Faixas:
##### 1 - 150 mil
##### 2 - 275 mil
##### 3 - 400 mil


ma_7 <- function(x) {
  x <- as.numeric(lag(x, n = 1L))
  y <- as.numeric((rollmean(x, 7, align = "right", fill = NA)))
  y[is.na(y)] <- 0
  
  return(y)
}

ma_21 <- function(x) {
  x <- as.numeric(lag(x, n = 1L))
  y <- as.numeric((rollmean(x, 21, align = "right", fill = NA)))
  y[is.na(y)] <- 0
  
  return(y)
}

macd <- function(x) {
  x <- as.numeric(lag(x, n = 1L))
  c <- as.numeric((rollmean(x, 12, align = "right", fill = NA)))
  d <- as.numeric((rollmean(x, 21, align = "right", fill = NA)))
  macd <- c - d
  macd[is.na(macd)] <- 0
  
  return(macd)
}

upper_band <- function(x) {
  x <- as.numeric(x)
  x <-  lag(x, n = 1L)
  up <- as.numeric((rollmean(x, 21, align = "right", fill = NA)))
  sd <- sd(na.omit(up))
  up <- up + (sd * 2)
  
  return(up)
}

lower_band <- function(x) {
  x <- as.numeric(x)
  x <-  lag(x, n = 1L)
  low <- as.numeric((rollmean(x, 21, align = "right", fill = NA)))
  sd <- sd(na.omit(low))
  low <- low - (sd * 2)
  
  return(low)
}

trend_period <- function(x, y) {
  x <- as.numeric(lag(x, n = 1L))
  x[is.na(x)] <- 0
  x <- if_else(x > lag(x, n = 1L), 1, 0)
  
  test <- as.numeric((rollsum(x, y, align = "right", fill = NA)))
  test[is.na(test)] <- 0
  
  return(test)
}

rev_mm <- function(x) {
  x <-  round(ma_7(lag(x, n = 1L)), 0)
  y <-  round(ma_21(lag(x, n = 1L)), 0)
  
  rev_mm <- as.numeric(x > y)
  rev_mm[is.na(rev_mm)] <- 0
  return(rev_mm)
}

rev_macd <- function(x) {
  x <- as.numeric(lag(x, n = 1L))
  x <-  macd(x)
  x <- as.numeric((rollmean(x, 12, align = "right", fill = NA)))
  y <- as.numeric((rollmean(x, 21, align = "right", fill = NA)))
  
  rev_mm <- as.numeric(x > y)
  rev_mm[is.na(rev_mm)] <- 0
  return(rev_mm)
}

engolf_alta <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  opening <- if_else(round(lag(x, n = 2L), 1) > round(lag(x, n = 1L), 1), 1, 0)
  opening[is.na(opening)] <- 0
  
  closing <- if_else(round(lag(y, n = 2L), 1) > round(lag(y, n = 1L), 1), 1, 0)
  closing[is.na(closing)] <- 0
  
  engolf_alta <- opening * closing
  
  return(engolf_alta)
  
}

kicking <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  kicking <- if_else(round(lag(x, n = 2L), 1) > round(lag(y, n = 1L), 1), 1, 0)
  kicking[is.na(kicking)] <- 0
  
  return(kicking)
}

discret_df <- function(x) {
  
  Data <- as_tibble(data.frame(date = index(x), coredata(x)))
  colnames(Data)[1] <- "date"
  colnames(Data)[2] <- "Open"
  colnames(Data)[3] <- "High"
  colnames(Data)[4] <- "Low"
  colnames(Data)[5] <- "Close"
  colnames(Data)[6] <- "Volume"
  
  Data <- Data %>% dplyr::left_join(USDBRL,  by = "date")
  Data <- Data %>% dplyr::left_join(BVSP,  by = "date")
  Data <- Data %>% dplyr::left_join(NASDAQ,  by = "date")
  Data <- Data %>% dplyr::left_join(DOWJONES,  by = "date")
  Data <- Data %>% dplyr::left_join(S_P,  by = "date")
  
  Data <- Data %>%
    mutate(
      class = if_else(Close > Open*1.005, 1, 0),
      close_maior_max = as.numeric(round(lag(Close, n = 2L), 0) > round(lag(High, n = 1L), 0)),
      close_maior_min = as.numeric(round(lag(Close, n = 2L), 0) > round(lag(Low, n = 1L), 0)),
      close_alta = as.numeric(lag(Close, n = 1L) > lag(Open, n = 1L)),
      max_igual_close = as.numeric(lag(Close, n = 1L) == lag(High, n = 1L)),
      min_igual_Open = as.numeric(lag(Open, n = 1L) == lag(Low, n = 1L)),
      doji = as.numeric(round(lag(Close, n = 1L), 1) ==  round(lag(Open, n = 1L), 1)),
      engolf_alta = engolf_alta(Open, Close),
      kicking = kicking(Open, Close),
      
      tend_2_p = if_else(trend_period(Close, 2) > 1, 1, 0),
      tend_4_p = if_else(trend_period(Close, 4) > 1, 1, 0),
      tend_9_p = if_else(trend_period(Close, 9) > 1, 1, 0),
      tend_17_p = if_else(trend_period(Close, 17) > 1, 1, 0),
      
      max_upper_band = if_else(upper_band(lag(Close, n = 1L)) < lag(High, n = 1L), 1, 0),
      close_maior_upper_band = if_else(upper_band(lag(Close, n = 1L)) < lag(Close, n = 1L), 1, 0),
      close_menor_lower_band = if_else(lower_band(lag(Close, n = 1L)) > lag(Close, n = 1L), 1, 0),
      min_lower_band = if_else(lower_band(lag(Close, n = 1L)) < lag(Low, n = 1L), 1, 0),
      pos_macd = if_else(macd(Close) > 0, 1, 0),
      rev_macd = rev_macd(Close),
      rev_mm = rev_mm(Close),
      Open_maior_mm_7 = if_else(ma_7(lag(Close, n = 1L)) < lag(Open, n = 1L), 1, 0),
      Low_maior_mm_7 = if_else(ma_7(lag(Close, n = 1L)) < lag(Low, n = 1L), 1, 0),
      High_maior_mm_7 = if_else(ma_7(lag(Close, n = 1L)) < lag(High, n = 1L), 1, 0),
      Close_maior_mm_7 = if_else(ma_7(lag(Close, n = 1L)) < lag(Close, n = 1L), 1, 0),
      Open_maior_mm_21 = if_else(ma_21(lag(Close, n = 1L)) < lag(Open, n = 1L), 1, 0),
      Low_maior_mm_21 = if_else(ma_21(lag(Close, n = 1L)) < lag(Low, n = 1L), 1, 0),
      High_maior_mm_21 = if_else(ma_21(lag(Close, n = 1L)) < lag(High, n = 1L), 1, 0),
      Close_maior_mm_21 = if_else(ma_21(lag(Close, n = 1L)) < lag(Close, n = 1L), 1, 0),
      mm_7_vol = if_else(ma_7(lag(Volume, n = 1L)) > lag(Volume, n = 1L), 1, 0),
      mm_21_vol = if_else(ma_21(lag(Volume, n = 1L)) > lag(Volume, n = 1L), 1, 0),
      
      dol_close_maior_max = as.numeric(round(lag(USDBRL.X.Close, n = 2L), 0) > 
                                         round(lag(USDBRL.X.High, n = 1L), 0)),
      dol_close_maior_min = as.numeric(round(lag(USDBRL.X.Close, n = 2L), 0) > 
                                         round(lag(USDBRL.X.Low, n = 1L), 0)),
      dol_close_alta = as.numeric((lag(USDBRL.X.Close, n = 1L)) > 
                                    (lag(USDBRL.X.Open, n = 1L))),
      dol_close_alta_2 = as.numeric((lag(USDBRL.X.Close, n = 1L)) > 
                                      (lag(USDBRL.X.Close, n = 2L))),
      dol_tend_2_p = if_else(trend_period(USDBRL.X.Close, 2) > 1, 1, 0),
      dol_tend_4_p = if_else(trend_period(USDBRL.X.Close, 4) > 1, 1, 0),
      dol_tend_9_p = if_else(trend_period(USDBRL.X.Close, 9) > 1, 1, 0),
      dol_Close_maior_mm_21 = if_else(ma_21(lag(Close, n = 1L)) < 
                                        lag(Close, n = 1L), 1, 0),
      dol_pos_macd = if_else(macd(Close) > 0, 1, 0),
      dol_rev_macd = rev_macd(Close),
      
      ibov_close_maior_max = as.numeric(round(lag(BVSP.Close, n = 2L), 0) > 
                                          round(lag(BVSP.High, n = 1L), 0)),
      ibov_close_maior_min = as.numeric(round(lag(BVSP.Close, n = 2L), 0) > 
                                          round(lag(BVSP.Low, n = 1L), 0)),
      ibov_close_alta = as.numeric(lag(BVSP.Close, n = 1L) > lag(BVSP.Open, n = 1L)),
      ibov_close_alta_2 = as.numeric(lag(BVSP.Close, n = 1L) > lag(BVSP.Close, n = 2L)),
      ibov_tend_2_p = if_else(trend_period(BVSP.Close, 2) > 1, 1, 0),
      ibov_tend_4_p = if_else(trend_period(BVSP.Close, 4) > 1, 1, 0),
      ibov_tend_9_p = if_else(trend_period(BVSP.Close, 9) > 1, 1, 0),
      ibov_Close_maior_mm_21 = if_else(ma_21(lag(Close, n = 1L)) < 
                                         lag(Close, n = 1L), 1, 0),
      ibov_pos_macd = if_else(macd(Close) > 0, 1, 0),
      ibov_rev_macd = rev_macd(Close),
      
      nasdq_close_maior_max = as.numeric(round(lag(IXIC.Close, n = 2L), 0) > 
                                           round(lag(IXIC.High, n = 1L), 0)),
      nasdq_close_maior_min = as.numeric(round(lag(IXIC.Close, n = 2L), 0) > 
                                           round(lag(IXIC.Low, n = 1L), 0)),
      nasdq_close_alta = as.numeric(lag(IXIC.Close, n = 1L) > 
                                      lag(IXIC.Open, n = 1L)),
      nasdq_close_alta_2 = as.numeric(lag(IXIC.Close, n = 1L) > 
                                        lag(IXIC.Close, n = 2L)),
      nasdq_tend_2_p = if_else(trend_period(IXIC.Close, 2) > 1, 1, 0),
      nasdq_tend_4_p = if_else(trend_period(IXIC.Close, 4) > 1, 1, 0),
      nasdq_tend_9_p = if_else(trend_period(IXIC.Close, 9) > 1, 1, 0),
      nasdq_Close_maior_mm_21 = if_else(ma_21(lag(Close, n = 1L)) < 
                                          lag(Close, n = 1L), 1, 0),
      nasdq_pos_macd = if_else(macd(Close) > 0, 1, 0),
      nasdq_rev_macd = rev_macd(Close),
      
      dwjon_close_maior_max = as.numeric(round(lag(DJI.Close, n = 2L), 0) > 
                                           round(lag(DJI.High, n = 1L), 0)),
      dwjon_close_maior_min = as.numeric(round(lag(DJI.Close, n = 2L), 0) > 
                                           round(lag(DJI.Low, n = 1L), 0)),
      dwjon_close_alta = as.numeric(lag(DJI.Close, n = 1L) > lag(DJI.Open, n = 1L)),
      dwjon_close_alta_2 = as.numeric(lag(DJI.Close, n = 1L) > lag(DJI.Close, n = 2L)),
      dwjon_tend_2_p = if_else(trend_period(DJI.Close, 2) > 1, 1, 0),
      dwjon_tend_4_p = if_else(trend_period(DJI.Close, 4) > 1, 1, 0),
      dwjon_tend_9_p = if_else(trend_period(DJI.Close, 9) > 1, 1, 0),
      dwjon_Close_maior_mm_21 = if_else(ma_21(lag(Close, n = 1L)) < lag(Close, n = 1L), 1, 0),
      dwjon_pos_macd = if_else(macd(Close) > 0, 1, 0),
      dwjon_rev_macd = rev_macd(Close),
      
      s_p_close_maior_max = as.numeric(round(lag(GSPC.Close, n = 2L), 0) > 
                                         round(lag(GSPC.High, n = 1L), 0)),
      s_p_close_maior_min = as.numeric(round(lag(GSPC.Close, n = 2L), 0) > 
                                         round(lag(GSPC.Low, n = 1L), 0)),
      s_p_close_alta = as.numeric(lag(GSPC.Close, n = 1L) > lag(GSPC.Open, n = 1L)),
      s_p_close_alta_2 = as.numeric(lag(GSPC.Close, n = 1L) > lag(GSPC.Close, n = 2L)),
      s_p_tend_2_p = if_else(trend_period(GSPC.Close, 2) > 1, 1, 0),
      s_p_tend_4_p = if_else(trend_period(GSPC.Close, 4) > 1, 1, 0),
      s_p_tend_9_p = if_else(trend_period(GSPC.Close, 9) > 1, 1, 0),
      s_p_Close_maior_mm_21 = if_else(ma_21(lag(Close, n = 1L)) < lag(Close, n = 1L), 1, 0),
      s_p_pos_macd = if_else(macd(Close) > 0, 1, 0),
      s_p_rev_macd = rev_macd(Close),
      day = lubridate::wday(as.Date(lag(date, n=1L)))

    ) %>%
    dplyr::select(class:day)
  Data[is.na(Data)] <- 0
  return(Data)
}

########################################################################################


Data <- rbind(discret_df(BBDC3.SA), discret_df(BBSE3.SA),
              discret_df(BRAP4.SA), discret_df(BRML3.SA),
              discret_df(CIEL3.SA), discret_df(CSAN3.SA),
              discret_df(CVCB3.SA), discret_df(ECOR3.SA),
              discret_df(EGIE3.SA), discret_df(FLRY3.SA),
              discret_df(HYPE3.SA), discret_df(JBSS3.SA),
              discret_df(LAME4.SA), discret_df(LREN3.SA),
              discret_df(SANB11.SA), discret_df(TIMP3.SA),
              discret_df(ABEV3.SA), discret_df(B3SA3.SA),
              discret_df(BRFS3.SA), discret_df(BRKM5.SA),
              discret_df(BTOW3.SA), discret_df(NATU3.SA),
              discret_df(MRVE3.SA), discret_df(MGLU3.SA),
              discret_df(BBAS3.SA), discret_df(BBDC4.SA),
              discret_df(CMIG4.SA), discret_df(CSNA3.SA),
              discret_df(CYRE3.SA), discret_df(EMBR3.SA),
              discret_df(GGBR3.SA), discret_df(GGBR4.SA),
              discret_df(ITSA4.SA), discret_df(ITUB4.SA),
              discret_df(PCAR4.SA), discret_df(PETR4.SA),
              discret_df(RAIL3.SA), discret_df(USIM5.SA),
              discret_df(VALE3.SA), discret_df(GOAU4.SA),
              discret_df(SUZB3.SA), discret_df(VVAR3.SA)
)


########################### XGBOOST TEST ###########################

set.seed(1)

train_data <- Data[sample(1:nrow(Data), nrow(Data) * opt_sample, replace = FALSE), ]
test_data <- Data[sample(1:nrow(Data), nrow(Data) * (1 - opt_sample), replace = FALSE), ]

dtrain <- xgb.DMatrix(data = as.matrix(train_data)[, -1], 
                      label = as.matrix(train_data)[, 1])

dtest <- xgb.DMatrix(data = as.matrix(test_data)[, -1], 
                     label = as.matrix(test_data)[, 1])

bstSparse <- xgb.train(data = dtrain, max_depth = 8,
                       eta = 1, nthread = 4, lambda = 5,
                       nrounds = 50, nfold = 25,
                       objective = "binary:logistic")

a <- as.matrix(predict(bstSparse, dtest))
data_result <- dplyr::as_tibble(cbind(test_data$class, a))

b <- cutpointr(data = data_result, x = a, 
               class = V1, pos_class = 1,
               direction = ">=",
               method = maximize_metric,
               metric = youden)
summary(b)

mat <- xgb.importance (feature_names = colnames(test_data), model = bstSparse)

#xgb.plot.importance (importance_matrix = mat[1:10])




################ XGBoost - Para Feature Importance e Eliminar Multicolinearidade #####################
#start 2h47 - multisession - end



tec_prob <- function(x) {
  
  #plan("multisession")
  
  Data <- as.matrix(discret_df(x))
  
  dtest <- xgb.DMatrix(data = as.matrix(Data)[, -1], label = as.matrix(Data)[, 1])
  
  a <- as.matrix(predict(bstSparse, dtest)  * 100) - 5 ##### ADICIONO CERTO GRAU DE INCERTEZA
  
  b <- as.character(as.data.frame(strsplit(names(x)[1], split = '.', fixed = TRUE))[1, 1])
  
  return(paste(b, "- Probabilidade de Alta:", round(last(a), 2)))
  
}


prob_diarias <- rbind(tec_prob(ABEV3.SA),tec_prob(B3SA3.SA),
                      tec_prob(BBAS3.SA),tec_prob(BBDC4.SA),
                      tec_prob(CSNA3.SA),tec_prob(CMIG4.SA),
                      tec_prob(CYRE3.SA),tec_prob(EMBR3.SA),
                      tec_prob(GGBR3.SA),tec_prob(GGBR4.SA),
                      tec_prob(ITSA4.SA),tec_prob(ITUB4.SA),
                      tec_prob(PCAR4.SA),tec_prob(PETR4.SA),
                      tec_prob(RAIL3.SA),tec_prob(USIM5.SA),
                      tec_prob(VALE3.SA),tec_prob(GOAU4.SA),
                      tec_prob(BIDI4.SA)
                      )

View(prob_diarias)
