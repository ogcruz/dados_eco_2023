# Importando o banco 
chuva <- read.csv2("https://gitlab.procc.fiocruz.br/oswaldo/eco2019/raw/master/exemplos/rain_ALL.csv")

# Colocando no formato ts a chuva da AP 1.0 '
chuva.ts  <- ts(chuva$AP_1.0, start = c(1981,1), frequency=12) 

# Gráfico da série temporal
plot(chuva.ts, main="Pluviosidade AP 1.0")

# Verificando a alto correlação da chuva no período
acf(chuva.ts, main="Autocorrelação da chuva na AP 1.0")

# Testando a autocorrelação 
# H0: rho_{h} = 0 
# H1: rho_{h}\neq 0
Box.test(chuva.ts, lag=20, type="Ljung-Box")

# Decompondo a série
# Via Moving Averages
plot(decompose(chuva.ts))

# Via stl (Seasonal Decomposition Of Time Series By Loess)
plot(stl(chuva.ts, s.window="periodic"))

# Testando a estacionariedade
# Teste de Dickey-Fuller
# H0: A série temporal não é Estacionária
# H1: A série temporais é Estacionária
library(tseries)
adf.test(chuva.ts)

# Avaliando a sazonalidade
boxplot(chuva.ts ~ cycle(chuva.ts))  
monthplot(chuva.ts) 
library(forecast)
seasonplot(chuva.ts,col=rainbow(6),lwd=2)  
ggseasonplot(chuva.ts,col=rainbow(6),lwd=2)  

# Utilizando trace = T, será possível verificar todo o processo de criação e teste dos modelos
modelo1 <-  auto.arima(chuva.ts, trace = F, allowdrift=F)
summary(modelo1)
tsdiag(modelo1)

# Neste modelo, será feito uma busca maior para uma solução "mais ótimizada"
modelo2 <-  auto.arima(chuva.ts, trace = F, stepwise = F, approximation = F, parallel = TRUE)
summary(modelo2)
tsdiag(modelo2)

# Utilizando trace = T, será possível verificar todo o processo de criação e teste dos modelos
modelo3 <-  auto.arima(log(chuva.ts), trace = F, allowdrift=F)
summary(modelo3)
tsdiag(modelo3)

# Neste modelo, será feito uma busca maior para uma solução "mais ótimizada"
modelo4 <-  auto.arima(log(chuva.ts), trace = F, stepwise = F, approximation = F,parallel = TRUE)
summary(modelo4)
tsdiag(modelo4)

# Teste de Shapiro-Wilk
shapiro.test(modelo1$residuals)
shapiro.test(modelo2$residuals)
shapiro.test(modelo3$residuals)
shapiro.test(modelo4$residuals)

# Comparando o ajuste dos modelo dado original
plot(chuva.ts)
lines(fitted(modelo1),col="red")
lines(fitted(modelo2),col="blue")
legend(1980, 500, legend=c("Dados", "Modelo 1", "Modelo 2"),
       col=c("black", "red", "blue"), lty=1:3, cex=1,
       title="Legenda", text.font=4, bg='lightblue')

AIC(modelo1)
AIC(modelo2)

# Comparando o ajuste dos modelos transformados pelo log
plot(log(chuva.ts))
lines(fitted(modelo3),col="red")
lines(fitted(modelo4),col="blue")
legend(1980, 6.3, legend=c("Dados", "Modelo 3", "Modelo 4"),
       col=c("black", "red", "blue"), lty=1:3, cex=1,
       title="Legenda", text.font=4, bg='lightblue')

AIC(modelo3)
AIC(modelo4)

# Comparando as predições de ambos os modelos
# Estimando as previsões
prev1 = forecast(modelo1, h=24)
plot(prev1) 

prev2 = forecast(modelo2, h=24)
plot(prev2) 

prev3 = forecast(modelo3, h=24)
plot(prev3) 

prev4 = forecast(modelo4, h=24)
plot(prev3) 

# Comparar os modelos com a escala original
plot(prev1, main='Comparação Modelo 1 vs Modelo 2')
lines(prev2$mean, col="red")
legend("topleft", legend=c("Dados", "Modelo 1", "Modelo 2"),
       col=c("black", "red", "blue"), lty=1:3, cex=1,
       title="Legenda", text.font=4, bg='lightblue')

# Comparar os modelos com a escala log
plot(prev3, main='Comparação Modelo 3 vs Modelo 4')
lines(prev4$mean, col="red")
legend("topleft", legend=c("Dados", "Modelo 3", "Modelo 4"),
       col=c("black", "red", "blue"), lty=1:3, cex=1,
       title="Legenda", text.font=4, bg='lightblue')
