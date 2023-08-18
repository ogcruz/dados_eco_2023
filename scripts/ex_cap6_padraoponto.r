library(splancs) # class .spp
library(spatstat) # class .ppp
library(gpclib)
library(maptools)

# lendo os bancos
# local <- 'https://gitlab.procc.fiocruz.br/oswaldo/eco2019/raw/master/dados/'
local <- '/home/tassinari/Documentos/cursos_ecologicos_2019/Bookdown/dados/'

homic <- read.table(paste0(local,"homic2.txt"), col.names=c("x","y"))
suic <- read.table(paste0(local,"suic.dat"), col.names=c("x","y"))
acid <- read.table(paste0(local,"acid2.txt"), col.names=c("x","y"))


##################################################################################################
# 1- Testar a CSR e explorar os Kernels para as causas de suicídio e acidentes de carro em POA/RS
##################################################################################################

# Construindo as caixas (contornos) segundo a distribuição de cada conjunto de pontos
homic.caixa <- as.owin(list(xrange=range(homic$x),yrange=range(homic$y)))
acid.caixa <- as.owin(list(xrange=range(acid$x),yrange=range(acid$y)))
suic.caixa <- as.owin(list(xrange=range(suic$x),yrange=range(suic$y)))

# Convertendo para a class ppp
homic.ppp <- as.ppp(homic.spp, homic.caixa)
acid.ppp <- as.ppp(acid.spp, acid.caixa)
suic.ppp <- as.ppp(suic.spp, caixa.suic)

# Construindo os quadrantes com as respectivas contagens
homicQ <- quadratcount(homic.ppp, nx = 4, ny = 4)
acidQ <- quadratcount(acid.ppp, nx = 4, ny = 4)
suicQ <- quadratcount(suic.ppp, nx = 4, ny = 4)


par(mfrow=c(1, 3))
plot(homicQ, main="Homicidio")
plot(homic.ppp, add = TRUE)

plot(acidQ, main="Acidente")
plot(acid.ppp, add = TRUE)

plot(suicQ, main="Suicidio")
plot(suic.ppp, add = TRUE)

# Testando a Completa Aletoriedade Espacial (CSR - complete spatial randomness)
# H0: Os pontos estão distribuidos aleatoriamente no espaço
# H1: Os pontos podem formar clusters ou estão dispersos no espaço
quadrat.test(homicQ)
quadrat.test(acidQ)
quadrat.test(suicQ)


# Fazendo os kernels
# Plotando os mapas

# contorno de porto alegre
homic.spp <- as.points(homic)
acid.spp <- as.points(acid)
suic.spp <- as.points(suic)
contorno.poa <- as.points(read.table(paste0(local,"contpoa.dat"),col.names=c("x","y")))

polymap(contorno.poa)
pointmap(homic.spp, pch=19, add=T,cex=0.5)
pointmap(suic.spp, pch=19, add=T, col=2,cex=0.5)
pointmap(acid.spp, pch=19, add=T, col=3,cex=0.5)

# Escolhendo a largura de banda de 200
homicKer <- kernel2d(homic.spp,contorno.poa,h0=200,nx=100,ny=100)
suicKer <- kernel2d(suic.spp,contorno.poa,h0=200,nx=100,ny=100)
acidKer <- kernel2d(acid.spp,contorno.poa,h0=200,nx=100,ny=100)

par(mfrow=c(1, 3))
polymap(contorno.poa,axes=FALSE)
image(homicKer, col=gray(32:0/32),add=TRUE)
pointmap(homic.spp, add=T, col=4,pch=19,cex=0.5)
title("Homicidios")

polymap(contorno.poa,axes=FALSE)
image(suicKer, col=gray(32:0/32),add=TRUE)
pointmap(suic.spp, add=T, col=4,pch=19,cex=0.5)
title("Suicidios")

polymap(contorno.poa,axes=FALSE)
image(acidKer, col=gray(32:0/32),add=TRUE)
pointmap(acid.spp, add=T, col=4,pch=19,cex=0.5)
title("Acidentes")

# Passando as linhas de contorno

par(mfrow=c(1, 3))
polymap(contorno.poa,axes=FALSE)
image(homicKer, col=gray(32:0/32),  main="Aleatorio", add=TRUE)
contour(homicKer,add=T)
title("Homicidios")

polymap(contorno.poa,axes=FALSE)
image(suicKer, col=gray(32:0/32), main="Regular", add=TRUE)
contour(suicKer, add=T)
title("Suicidios")

polymap(contorno.poa,axes=FALSE)
image(acidKer, col=gray(32:0/32), main="Cluster", add=TRUE)
contour(acidKer, add=T)
title("Acidentes")

# Estimando a largura de banda otima

homic.Mse2d <- mse2d(homic.spp, contorno.poa, nsmse=20,range=0.2)
suic.Mse2d <- mse2d(suic.spp, contorno.poa, nsmse=20,range=0.2)
acid.Mse2d <- mse2d(acid.spp, contorno.poa, nsmse=20,range=0.2)

par(mfrow=c(1, 3))
plot(homic.Mse2d$h, homic.Mse2d$mse, type="l", xlab="h",  ylab="Mse Homicidios")
plot(suic.Mse2d$h, suic.Mse2d$mse, type="l", xlab="h",  ylab="Mse Suicidios")
plot(acid.Mse2d$h, acid.Mse2d$mse, type="l", xlab="h",  ylab="Mse Acidentes")

##################################################################################################
# 2 - Fazer a razão de kernel entre suicídios e acidentes de carro em POA/RS
##################################################################################################

suic_acid.ratio <- kernrat(suic.spp, acid.spp, contorno.poa, h1=300, h2=600, nx=100, ny=100)

polymap(contorno.poa,axes=FALSE)
image(suic_acid.ratio, col=gray(32:0/32), add=TRUE)
pointmap(suic.spp, add=T, col=2,pch=19,cex=0.7)
pointmap(acid.spp, add=T, col=4,pch=19,cex=0.7)
title("Razao Suicidios/Acidentes")
legend("topleft", legend=c("Suicidios", "Acidentes"),
       col=c("red", "blue"),  pch = 19, cex=1,
       title="Legenda", text.font=4, bg='lightblue')



##################################################################################################
# 3 - Inspecionar o processo pontual de segunta ordem, utilizando as funções K, G e L para as causas de suicídio e acidentes de carro em POA/RS.
##################################################################################################

# Funcao K 
par(mfrow=c(1, 3))
plot(envelope(Y = homic.ppp, fun = Kest, nsim = 99), main="Homicidios")
plot(envelope(Y = suic.ppp, fun = Kest, nsim = 99), main="Suicidios")
plot(envelope(Y = acid.ppp, fun = Kest, nsim = 99), main="Acidentes")

# Funcao G
par(mfrow=c(1, 3))
plot(envelope(Y = homic.ppp, fun = Gest, nsim = 99), main="Homicidios")
plot(envelope(Y = suic.ppp, fun = Gest, nsim = 99), main="Suicidios")
plot(envelope(Y = acid.ppp, fun = Gest, nsim = 99), main="Acidentes")

# Funcao L
par(mfrow=c(1, 3))
plot(envelope(Y = homic.ppp, fun = Lest, nsim = 99), main="Homicidios")
plot(envelope(Y = suic.ppp, fun = Lest, nsim = 99), main="Suicidios")
plot(envelope(Y = acid.ppp, fun = Lest, nsim = 99), main="Acidentes")



