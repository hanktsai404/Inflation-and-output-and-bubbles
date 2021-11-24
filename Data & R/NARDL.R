# A NARDL implementation on bubbles, inflation and output
# Created by Hank Tsai
# 07-01-2021

# Libraries
library(ggplot2)
library(ggthemes)
library(scales)
library(fUnitRoots)
library(vars)
library(ARDL)
library(stats)
library(lmtest)
library(sandwich)
library(nlWaldTest)
library(car)
library(moments)

# Lag function
lag = function(seq, time)
{
  n = length(seq)
  for(i in c(1:time))
  {
    seq = c(NA, seq[-n])
  }
  return(seq)
}

df = read.csv("data_preprocessed.csv")
df$found = df$CAPE - df$BUBBLE

summary(df)


# Identify pos and neg change
dinf = df$INF - lag(df$INF, 1)
dipi = df$IPI - lag(df$IPI, 1)
dinf_pos = c()
dinf_neg = c()
inf_pos = c()
inf_neg = c()
dipi_pos = c()
dipi_neg = c()
ipi_pos = c()
ipi_neg = c()

for(i in c(1: length(dinf)))
{
  dinf_pos = c(dinf_pos, max(c(dinf[i], 0)))
  dinf_neg = c(dinf_neg, min(c(dinf[i], 0)))
  inf_pos = c(inf_pos, sum(dinf_pos, na.rm = TRUE))
  inf_neg = c(inf_neg, sum(dinf_neg, na.rm = TRUE))
  dipi_pos = c(dipi_pos, max(c(dipi[i], 0)))
  dipi_neg = c(dipi_neg,min(c(dipi[i], 0)))
  ipi_pos = c(ipi_pos, sum(dipi_pos, na.rm = TRUE))
  ipi_neg = c(ipi_neg, sum(dipi_neg, na.rm = TRUE))
}
df$dinf = dinf
df$dinf_pos = dinf_pos
df$dinf_neg = dinf_neg
df$inf_pos = inf_pos
df$inf_neg = inf_neg
df$dipi = dipi
df$dipi_pos = dipi_pos
df$dipi_neg = dipi_neg
df$ipi_pos = ipi_pos
df$ipi_neg = ipi_neg

# Graphs
ggplot(df, aes(x = as.Date(Date), group = 1)) + geom_line(aes(y = CAPE, colour = "CAPE"), size = 1) +
  geom_line(aes(y = THRES, colour = "Threshold"), size = 1) + labs(x = "Time", y = "CAPE") +
  theme_economist_white() + scale_colour_economist() + ggtitle("CAPE") +
  scale_x_date(date_breaks = "2 year", date_labels = "%Y") +
  scale_colour_manual("", values = c("CAPE" = "#00AEAE", "Threshold" = "red"))

ggplot(df, aes(x = as.Date(Date), y = BUBBLE), group = 1) +geom_line(size = 1) +
  labs(x = "Time", y = "Bubbles") + theme_economist_white() + ggtitle("Bubbles") +
  scale_x_date(date_breaks = "2 year", date_labels = "%Y")

# adf test
adfTest(df$INF) #I(0)
adfTest(df$IPI) #I(1)
adfTest(df$BUBBLE) #I(1) DF statistic is below zero
trans_bubble = df$BUBBLE - min(df$BUBBLE) + 1
df$dbub = log(trans_bubble) - log(lag(trans_bubble, 1))
adfTest(df$dbub)

# lag selection: ardl AIC and BIC grid search
best_mod_aic = uecm(dbub ~ inf_pos + inf_neg + ipi_pos + ipi_neg, data = df, order = c(1,1,1,1,1))
best_mod_bic = uecm(dbub ~ inf_pos + inf_neg + ipi_pos + ipi_neg, data = df, order = c(1,1,1,1,1))
best_p_aic = 1
best_q_aic = 1
best_p_bic = 1
best_q_bic = 1
for(p in c(1:12)){
  for(q in c(1:12)){
    mod = uecm(dbub ~ inf_pos + inf_neg + ipi_pos + ipi_neg, data = df, order = c(p,q,q,q,q))
    if(AIC(mod)<AIC(best_mod_aic)){
      best_mod_aic = mod
      best_p_aic = p
      best_q_aic = q
    }
    if(BIC(mod)<BIC(best_mod_bic)){
      best_mod_bic = mod
      best_p_bic = p
      best_q_bic = q
    }
  }
}
best_p_aic
best_q_aic
best_p_bic
best_q_bic
# Optimal lag: (6,8) by AIC

# NARDL
mod = uecm(dbub ~ inf_pos + inf_neg + ipi_pos + ipi_neg, data = df, order = c(6,8,8,8,8))
summary(mod) # inf up bub down ipi up bub down
hac = NeweyWest(mod)
coeftest(mod, vcov. = hac)
bounds_f_test(mod, case = 3, alpha = 0.05, vcov_matrix = hac)


# Wald test: inf
nlWaldtest(mod, texts = c("(b[3]/b[2]) - (b[4]/b[2]) = 0"), Vcov = hac) # No asymmetric in long inf
nlWaldtest(mod, texts = 
             c("b[12] + b[13] + b[14] + b[15] + b[16] + b[17] + b[18] + b[19] - b[20] - b[21] - b[22] - b[23] - b[24] - b[25] - b[26] - b[27] = 0"), Vcov = hac) # No asymmetric in short inf

# Wald test: ipi
nlWaldtest(mod, texts = c("(b[5]/b[2]) - (b[6]/b[2]) = 0"), Vcov = hac) # No asymmetric in long ipi
nlWaldtest(mod, texts = 
             c("b[28] + b[29] + b[30] + b[31] + b[32] + b[33] +b[34] + b[35] - b[36] - b[37] - b[38] - b[39] - b[40] - b[41] - b[42] - b[43] = 0"), Vcov = hac) # No asymmetric in short ipi

# Diagnostics
bptest(mod) # Heteroskedasticity problem detected -> fixed with HAC std error
dwtest(mod) # No serial correlation

# Fitted graph
df$fit = c(NA, NA, NA, NA, NA, NA, NA, NA, mod$fitted.values)
df$res1 = c(NA, NA, NA, NA, NA, NA, NA, NA, mod$residuals)
df$ddbub = df$dbub - lag(df$dbub, 1)
ggplot(df, aes(x = as.Date(Date), group = 1)) + geom_line(aes(y = ddbub, colour = "Actual bubble growth"), size = 0.8) +
  geom_line(aes(y = fit, colour = "fit"), size = 1) + labs(x = "Time", y = "growth") +
  theme_economist_white() + scale_colour_economist() + ggtitle("Fitted values of NARDL") +
  scale_x_date(date_breaks = "2 year", date_labels = "%Y") +
  scale_colour_manual("", values = c("Actual bubble growth" = "#2828FF", "fit" = "#FF9224"))

ggplot(df, aes(x = as.Date(Date), y = res1), group = 1) + geom_line() + labs(x = "Time", y = "Residuals") +
  theme_economist_white() + scale_colour_economist() + ggtitle("Residuals") +
  scale_x_date(date_breaks = "2 year", date_labels = "%Y")





# Model for foundamental
# lag selection
df$dfound = log(df$found) - log(lag(df$found, 1))
best_mod_aic = uecm(dfound ~ inf_pos + inf_neg + ipi_pos + ipi_neg, data = df, order = c(1,1,1,1,1))
best_mod_bic = uecm(dfound ~ inf_pos + inf_neg + ipi_pos + ipi_neg, data = df, order = c(1,1,1,1,1))
best_p_aic = 1
best_q_aic = 1
best_p_bic = 1
best_q_bic = 1
aic_arr = c()
for(p in c(1:12)){
  for(q in c(1:12)){
    mod = uecm(dfound ~ inf_pos + inf_neg + ipi_pos + ipi_neg, data = df, order = c(p,q,q,q,q))
    if(AIC(mod)<AIC(best_mod_aic)){
      best_mod_aic = mod
      best_p_aic = p
      best_q_aic = q
    }
    if(BIC(mod)<BIC(best_mod_bic)){
      best_mod_bic = mod
      best_p_bic = p
      best_q_bic = q
    }
  }
}
best_p_aic
best_q_aic
best_p_bic
best_q_bic

# NARDL
adfTest(df$found)
adfTest(df$dfound)

mod = uecm(dfound ~ inf_pos + inf_neg + ipi_pos + ipi_neg, data = df, order = c(6,8,8,8,8))
bptest(mod)
dwtest(mod)  # heteroskedasticity, no serial correlation
summary(mod)  # inf up found up ipi up inf up
hac = NeweyWest(mod)
coeftest(mod, vcov. = hac)
bounds_f_test(mod, case = 3, alpha = 0.05, vcov_matrix = hac)  # Cointegration detected

# Wald test: inf
nlWaldtest(mod, texts = c("(b[3]/b[2]) - (b[4]/b[2]) = 0"), Vcov = hac) # Asymmetric in long inf
nlWaldtest(mod, texts = 
             c("b[12] + b[13] + b[14] + b[15] + b[16] + b[17] + b[18] + b[19] - b[20] - b[21] - b[22] - b[23] - b[24] - b[25] - b[26] - b[27] = 0"), Vcov = hac) # No asymmetric in short inf

# Wald test: ipi
nlWaldtest(mod, texts = c("(b[5]/b[2]) - (b[6]/b[2]) = 0"), Vcov = hac) # Asymmetric in long ipi
nlWaldtest(mod, texts = 
             c("b[28] + b[29] + b[30] + b[31] + b[32] + b[33] +b[34] + b[35] - b[36] - b[37] - b[38] - b[39] - b[40] - b[41] - b[42] - b[43] = 0"), Vcov = hac) # Asymmetric in short ipi


df$fit = c(NA, NA, NA, NA, NA, NA, NA, NA, mod$fitted.values)
df$res2 = c(NA, NA, NA, NA, NA, NA, NA, NA, mod$residuals)
df$ddfound = df$dfound - lag(df$dfound, 1)

ggplot(df, aes(x = as.Date(Date), group = 1)) + geom_line(aes(y = ddfound, colour = "Actual fundamental growth"), size = 0.8) +
  geom_line(aes(y = fit, colour = "fit"), size = 1) + labs(x = "Time", y = "growth") +
  theme_economist_white() + scale_colour_economist() + ggtitle("Fitted values of NARDL") +
  scale_x_date(date_breaks = "2 year", date_labels = "%Y") +
  scale_colour_manual("", values = c("Actual fundamental growth" = "#2828FF", "fit" = "#FF9224"))

ggplot(df, aes(x = as.Date(Date), y = res2), group = 1) + geom_line() + labs(x = "Time", y = "Residuals") +
  theme_economist_white() + scale_colour_economist() + ggtitle("Residuals") +
  scale_x_date(date_breaks = "2 year", date_labels = "%Y")

df$res1_std = (df$res1 - mean(df$res1, na.rm = TRUE))/sd(df$res1, na.rm = TRUE)
df$res2_std = (df$res2 - mean(df$res2, na.rm = TRUE))/sd(df$res2, na.rm = TRUE)
ggplot(df, aes(x = as.Date(Date), group = 1)) + geom_line(aes(y = res1_std, colour = "Residuals for bubbles"), size = 0.8) +
  geom_line(aes(y = res2_std, colour = "Residuals for foundamental"), size = 1) + labs(x = "Time", y = "Residuals") +
  theme_economist_white() + scale_colour_economist() + ggtitle("Standardized Residuals") +
  scale_x_date(date_breaks = "2 year", date_labels = "%Y") +
  scale_colour_manual("", values = c("Residuals for bubbles" = "#2828FF", "Residuals for foundamental" = "#FF9224"))


ggplot(df, aes(x = as.Date(Date), group = 1)) + geom_line(aes(y = dbub, colour = "Bubble growth"), size = 1) +
  geom_line(aes(y = dfound, colour = "Fundamental Growth"), size = 1) + geom_line(aes(y = INF, colour = "Inflation"), size = 1) +
  labs(x = "Time", y = "Growth") +
  theme_economist_white() + scale_colour_economist()+
  scale_x_date(date_breaks = "2 year", date_labels = "%Y") +
  scale_colour_manual("", values = c("Bubble growth" = "#00AEAE", "Fundamental Growth" = "red", "Inflation" = "orange"))



