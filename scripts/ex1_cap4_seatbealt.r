data(Seatbelts)
class(Seatbelts)
Seatbelts

# Analisando a série multivariada 
plot(Seatbelts)
acf(Seatbelts)

# Convertendo o dataset Seatbelts em formato data.frame
# floor : arredondamento para os anos 
# cycle : ciclos (nesse caso são em meses)
# labels=month.abb : Aparecer os meses por nomes
Seatbelts.df <- data.frame(ano=floor(time(Seatbelts)),
                           mes=factor(cycle(Seatbelts), labels=month.abb),
                           mes2=cycle(Seatbelts), Seatbelts)
names(Seatbelts.df)

# Começando pelo modelo ARIMA
# Colocando no formato ts somente a variavel DriversKilled
drivers.ts  <- ts(Seatbelts.df$DriversKilled,start = c(1964,1), frequency=12) 

drivers.ts

plot(drivers.ts)

acf(drivers.ts, main="DriversKilled")

# Testando a autocorrelação 
# H0: rho_{h} = 0 
# H1: rho_{h}\neq 0
Box.test(drivers.ts, lag=20, type="Ljung-Box")

# Decompondo a série
plot(decompose(drivers.ts))
plot(stl(drivers.ts, s.window="periodic"))

# Testando a estacionariedade
# Teste de Dickey-Fuller
# H0: A série temporal não é Estacionária
# H1: A série temporais é Estacionária
library(tseries)
adf.test(drivers.ts)

# Avaliando a sazonalidade
boxplot(drivers.ts ~ cycle(drivers.ts))  
monthplot(drivers.ts) 
library(forecast)
seasonplot(drivers.ts,col=rainbow(6),lwd=2)  
ggseasonplot(drivers.ts,col=rainbow(6),lwd=2)  

# tilizando trace = T, será possível verificar todo o processo de criação e teste dos modelos
modelo1 <-  auto.arima(drivers.ts, trace = F, allowdrift=F)
modelo1
tsdiag(modelo1)

# Neste modelo, será feito uma busca maior para uma solução "mais ótimizada"
modelo2 <-  auto.arima(drivers.ts, trace = F, stepwise = F, approximation = F,parallel = TRUE)
modelo2
tsdiag(modelo2)

# Teste de normalidade dos resíduos Shapiro-Wilk
shapiro.test(modelo1$residuals)
shapiro.test(modelo2$residuals)


# Comparando o ajuste dos modelos
plot(drivers.ts)
lines(fitted(modelo1),col="red")
lines(fitted(modelo2),col="blue")
legend(1976, 200, legend=c("Dados", "Modelo 1", "Modelo 2"),
       col=c("black", "red", "blue"), lty=1:3, cex=1,
       title="Legenda", text.font=4, bg='lightblue')

# Comparando as predições de ambos os modelos
# Estimando as previsões
prev1 = forecast(modelo1, h=24)
plot(prev1) 

prev2 = forecast(modelo2, h=24)
plot(prev2) 

# Comparar os modelos
plot(prev1, main='Comparação Modelo 1 vs Modelo 2')
lines(prev2$mean, col="red")
legend("topleft", legend=c("Modelo 1", "Modelo 2"),
       col=c("blue", "red"), lty=1:3, cex=1,
       title="Legenda", text.font=4, bg='lightblue')

# Correlação cruzada Motoristas de carro mortos x Preço da gasolina
plot(Seatbelts.df$PetrolPrice, Seatbelts.df$DriversKilled, xlab="Preço da gasolina", ylab="Motoristas de carro mortos", pch=19)
abline(lm(Seatbelts.df$DriversKilled ~ Seatbelts.df$PetrolPrice), col="red", lwd=2) # regression line (y~x)
lines(lowess(Seatbelts.df$DriversKilled ~ Seatbelts.df$PetrolPrice), col="blue", lwd=2) # lowess line (x,y)

cor.test(Seatbelts.df$DriversKilled, Seatbelts.df$PetrolPrice)
plot(ccf(Seatbelts.df$DriversKilled, Seatbelts.df$PetrolPrice))

# Modelo GAM
library(mgcv)

hist(Seatbelts.df$DriversKilled)
summary(Seatbelts.df$DriversKilled)
var(Seatbelts.df$DriversKilled)


# modelo poisson vazio só com as componentes temporais
modGAM1 <- gam(DriversKilled ~ s(ano) + s(mes2), family = poisson(), data=Seatbelts.df)
summary(modGAM1)
plot(modGAM1)
acf(modGAM1$residuals)
Box.test(modGAM1$residuals, type="Ljung-Box")
shapiro.test(modGAM1$residuals)

# modelo poisson full com as variáveis explicativas e as componentes temporais
modGAM2 <- gam(DriversKilled ~  s(kms) + s(rear) + s(PetrolPrice) + s(ano) + s(mes2), family = poisson(), data=Seatbelts.df)
summary(modGAM2)
plot(modGAM2)
acf(modGAM2$residuals)
Box.test(modGAM2$residuals, type="Ljung-Box")
shapiro.test(modGAM2$residuals)

# modelo binomial-negativo vazio só com as componentes temporais
modGAM3 <- gam(DriversKilled ~ s(ano) + s(mes2), family = nb(), data=Seatbelts.df)
summary(modGAM3)
plot(modGAM3)
acf(modGAM3$residuals)
Box.test(modGAM3$residuals, type="Ljung-Box")
shapiro.test(modGAM3$residuals)

# modelo binomial-negativo full com as variáveis explicativas e as componentes temporais
modGAM4 <- gam(DriversKilled ~ s(kms) + s(rear) + s(PetrolPrice) + s(ano) + s(mes2), family = nb(), data=Seatbelts.df)
summary(modGAM4)
plot(modGAM4)
acf(modGAM4$residuals)
Box.test(modGAM4$residuals, type="Ljung-Box")
shapiro.test(modGAM4$residuals)

AIC(modGAM1)
AIC(modGAM2)
AIC(modGAM3)
AIC(modGAM4)





