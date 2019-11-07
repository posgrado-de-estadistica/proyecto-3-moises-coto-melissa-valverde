##Bibliotecas
library(sp)
library(rgdal)
library(raster)
library(spatstat)
library(dplyr)
library(maptools)
library(lattice)
library(RColorBrewer)
library(sp)
library(lattice)
library(georob)
library(multcomp)
library(car)
library(gstat)
library(raster)

setwd("C:/Users/mvalverde/Documents/Proyecto_3/")

base1 <- read.csv2("C:/Users/mvalverde/Documents/Proyecto_3/SISMOS_CR_LIS_UCR_DF.csv", header=TRUE)

##Creando SpatialPointsDataFrame

sp_point <- cbind(base1$Lon, base1$Lat)
colnames(sp_point) <- c("Lon","Lat")
head(sp_point)
proj <- CRS("+proj=utm +zone=16 +datum=WGS84")
data.sp <- SpatialPointsDataFrame(coords=sp_point,base1,proj4string=proj)

##Verificando los duplicados 

if((cbind(base1$Lon,base1$Lat) %>% duplicated() %>% sum())==0)
{
  print("No hay duplicados")
} else {
  print("Hay duplicados")
}

##Definiendo la ventana de analisis

costarica <- readRDS("CRI_adm2.rds")
margen<-
  xcord<-(-c(85.60,84.92,84.66,83.88,83.71,83.69,83.43,83.09,83.03,82.52,82.65,
             82.83,82.82,82.91,82.91,82.67,82.88,82.78,82.97,82.83,82.92,
             83.01,83.18,83.21,83.38,83.25,83.27,83.59,83.77,83.64,84.08,
             84.73,85.05,84.85,85.11,85.43,85.69,85.91,85.69,85.98))

ycord<-c(11.27,10.99,11.13,10.76,10.82,10.99,10.45,10.03,10.02,9.55,
         9.47,9.59,9.49,9.45,9.12,8.95,8.74,8.52,8.32,8.04,8.02,8.24,
         8.37,8.59,8.66,8.49,8.34,8.43,8.62,8.99,9.33,9.53,10.06,
         9.85,9.53,9.82,9.90,10.32,10.80,10.91)

wndw <- owin(poly = list(x = rev(xcord),
                         y = rev(ycord)))

##Creando el objeto tipo ppp

base<- ppp(data.sp@coords[,1],data.sp@coords[,2],window=wndw)
names(base)

sismos<- ppp(data.sp@coords[,1],data.sp@coords[,2],window=wndw)
names(sismos)

##Graficando directamente el objeto ppp

par(mfrow=c(1,1))
plot(sismos)
plot(sismos,pch=16,main="", cex=.6,axes=TRUE,cols="black",cex.main=0.5,xlim=c(-85.5,-83),ylim=c(8,11.5),border=F)
mtext("Distribución de los sismos en Costa Rica", side =3, line = 2, outer = FALSE, at = NA,
      adj = 0.5, padj = 1, cex = 1.5, col = 1, font = 4)
plot(costarica,xlim=c(-85.5,-83),ylim=c(8,11.5),add=T)

##Gráficos para visualizar la densidad

plot(sismos,pch=16,main="", cex=.6,axes=TRUE,cols="green",cex.main=0.5,xlim=c(-85.5,-83),ylim=c(8,11.5),border=F)
plot(costarica,xlim=c(-85.5,-83),ylim=c(8,11.5),add=T)
plot(density(sismos))
plot(costarica,xlim=c(-85.5,-83),ylim=c(8,11.5),add=T)
contour(density(sismos))
plot(costarica,xlim=c(-85.5,-83),ylim=c(8,11.5),add=T)

plot(log(Magnitud_Mw)~log(Profundidad_km), data.sp@data, pch=as.integer(Tipo), col=Tipo)
plot(log(Magnitud_Mw)~log(Acelaracion_cms2), data.sp@data, pch=as.integer(Tipo),col=Tipo)

plot(log(Magnitud_Mw)~1, data.sp@data, pch=as.integer(Tipo),col=Tipo)
plot(log(Profundidad_km)~1, data.sp@data, pch=as.integer(Tipo), col=Tipo)
plot(log(Acelaracion_cms2)~1, data.sp@data, pch=as.integer(Tipo),col=Tipo)

xyplot(log(Magnitud_Mw)~log(Acelaracion_cms2) | Tipo, data.sp@data, groups=Tipo, panel=function(x, y, ...){ panel.xyplot(x, y, ...)
  panel.loess(x, y, ...)
}, auto.key=TRUE)

xyplot(log(Magnitud_Mw)~log(Profundidad_km) | Tipo, data.sp@data, groups=Tipo, panel=function(x, y, ...){ panel.xyplot(x, y, ...)
  panel.loess(x, y, ...)
}, auto.key=TRUE)


## Analisis exploratorio y modelamiento
### variogramas

## Magnitud

r.lm1<- lm(log(Magnitud_Mw)~1,data.sp@data)
summary(r.lm1)

op <- par(mfrow=c(2, 2)); plot(r.lm1); par(op)

v <- variogram(log(Magnitud_Mw)~1,data=data.sp, width=1)
plot(v)


plot(r.sv1 <- sample.variogram(residuals(r.lm1), locations=data.sp@data[, c("Lat","Lon")],
                               lag.dist.def=0.5, max.lag=10,
                               estimator="matheron"), type="l",
     main="Variograma de los residuos del modelo log(Magnitud_Mw)~1 (Cressie)")
lines(r.sv.spher <- fit.variogram.model(r.sv1, variogram.mode="RMexp",
                                        param=c(variance=0.02, nugget=0.05, scale=1.2)),col="red")

plot(variogram(log(Magnitud_Mw)~1,data=data.sp,  width=1,alpha = c(0, 45, 90, 135)))


## Ajuste modleo lineal guassiano

r.georob.m0.spher.reml1 <- georob(log(Magnitud_Mw)~1, data.sp, locations=~Lat+Lon,
                                  variogram.model="RMexp", param=c(variance=0.02, nugget=0.05, scale=1.2),
                                  tuning.psi=1000)
summary(r.georob.m0.spher.reml1)

r.georob.m0.spher.ml1 <- update(r.georob.m0.spher.reml1,
                                control=control.georob(ml.method="ML"))

summary(r.georob.m0.spher.ml1)

extractAIC(r.georob.m0.spher.reml1, REML=TRUE)

extractAIC(r.georob.m0.spher.ml1)

#############################

# Profundidad

r.lm1<- lm(log(Profundidad_km)~1,data.sp@data)
summary(r.lm1)

v <- variogram(log(Profundidad_km)~1,data=data.sp, width=1)
plot(v)


plot(r.sv1 <- sample.variogram(residuals(r.lm1), locations=data.sp@data[, c("Lat","Lon")],
                               lag.dist.def=0.5, max.lag=10,
                               estimator="matheron"), type="l",
     main="Variograma de los residuos del modelo log(Profundidad_km)~1 (Cressie)")
lines(r.sv.spher <- fit.variogram.model(r.sv1, variogram.mode="RMexp",
                                        param=c(variance=4, nugget=0.05, scale=40)),col="red")

fit.variogram(v, vgm(c("Exp", "Mat", "Sph")))

plot(variogram(log(Profundidad_km)~1,data=data.sp,  width=1,alpha = c(0, 45, 90, 135)))

summary(r.sv.spher)

## Ajuste modleo lineal guassiano

r.georob.m0.spher.reml2 <- georob(log(Profundidad_km)~1, data.sp, locations=~Lat+Lon,
                                  variogram.model="RMexp", param=c(variance=4, nugget=0.05, scale=40),
                                  tuning.psi=1000)
summary(r.georob.m0.spher.reml2)

r.georob.m0.spher.ml2 <- update(r.georob.m0.spher.reml2,
                                control=control.georob(ml.method="ML"))

summary(r.georob.m0.spher.ml2)

extractAIC(r.georob.m0.spher.reml2, REML=TRUE)

extractAIC(r.georob.m0.spher.ml2)

#########################

## Aceleración

r.lm1<- lm(log(Acelaracion_cms2)~1,data.sp@data)
summary(r.lm1)

v <- variogram(log(Acelaracion_cms2)~1,data=data.sp, width=1)
plot(v)


plot(r.sv1 <- sample.variogram(residuals(r.lm1), locations=data.sp@data[, c("Lat","Lon")],
                               lag.dist.def=0.5, max.lag=10,
                               estimator="matheron"), type="l",
     main="Variograma de los residuos del modelo log(Acelaracion_cms2)~1 (Cressie)")
lines(r.sv.spher <- fit.variogram.model(r.sv1, variogram.mode="RMspheric",
                                        param=c(variance=0.7, nugget=0.05, scale=2)),col="red")

fit.variogram(v, vgm(c("Exp", "Mat", "Sph","Gau")))

plot(variogram(log(Profundidad_km)~1,data=data.sp,  width=1,alpha = c(0, 45, 90, 135)))

summary(r.sv.spher)

## Ajuste modleo lineal guassiano

r.georob.m0.spher.reml3 <- georob(log(Profundidad_km)~1, data.sp, locations=~Lat+Lon,
                                  variogram.model="RMspheric", param=c(variance=0.7, nugget=0.05, scale=2),
                                  tuning.psi=1000)
summary(r.georob.m0.spher.reml3)

r.georob.m0.spher.ml3 <- update(r.georob.m0.spher.reml3,
                                control=control.georob(ml.method="ML"))

summary(r.georob.m0.spher.ml3)

extractAIC(r.georob.m0.spher.reml3, REML=TRUE)

extractAIC(r.georob.m0.spher.ml3)

#############################

## Magnitud y aceleración

r.lm2<- lm(log(Magnitud_Mw)~log(Acelaracion_cms2),data.sp@data)
summary(r.lm2)
op <- par(mfrow=c(2, 2)); plot(r.lm2); par(op)

v <- variogram(log(Magnitud_Mw)~log(Acelaracion_cms2),data=data.sp, width=1)
plot(v)

fit.variogram(v, vgm(c("Exp", "Mat", "Sph","Gau")))

plot(r.sv2 <- sample.variogram(residuals(r.lm2), locations=data.sp@data[, c("Lat","Lon")],
                               lag.dist.def=0.5, max.lag=10,
                               estimator="matheron"), type="l",
     main="Variograma de los residuos del modelo log(Magnitud_Mw)~log(Acelaracion_cms2) (Cressie)")
lines(r.sv.spher <- fit.variogram.model(r.sv2, variogram.mode="RMgauss",
                                        param=c(variance=29, nugget=0.005, scale=126)),col="red")

plot(variogram(log(Magnitud_Mw)~log(Acelaracion_cms2),data=data.sp,  width=1,alpha = c(0, 45, 90, 135)))

summary(r.sv.spher)

## Ajuste modleo lineal guassiano

r.georob.m0.spher.reml4 <- georob(log(Magnitud_Mw)~log(Acelaracion_cms2), data.sp, locations=~Lat+Lon,
                                  variogram.model="RMgauss", param=c(variance=29, nugget=0.005, scale=126),
                                  tuning.psi=1000)
summary(r.georob.m0.spher.reml4)

r.georob.m0.spher.ml4 <- update(r.georob.m0.spher.reml4,
                                control=control.georob(ml.method="ML"))

summary(r.georob.m0.spher.ml4)

extractAIC(r.georob.m0.spher.reml4, REML=TRUE)

extractAIC(r.georob.m0.spher.ml4)


################################

## Kriging Univariado

vt <- variogram(f, meuse)
v.fit1 <- fit.variogram(v, vgm(0.02, "Exp", 1.27, 0.05))
vt.fit <- fit.variogram(vt, vgm(1, "Exp", 300, 1))
vt.fit

r<-raster(data.sp,nrows=100,ncols=100)
coordinates(g) <- c("x", "y")
g<-as(r,"SpatialGrid")
plot(g)


## Ordinario
lz.sk <- krige(log(Magnitud_Mw)~1,data.sp, g, v.fit1, beta = 5.9)
summary(lz.sk )

summary(lz.sk$var1.pred)


## Universal

lz.ok <- krige(log(Magnitud_Mw)~1,data.sp, g, v.fit1)
summary(lz.sk )

#respuestas multivariadas CoKriging

cor(as.data.frame(data.sp@data)[c("Profundidad_km", "Magnitud_Mw", "Acelaracion_cms2")])

ck <- gstat(NULL, "logProf", log(Profundidad_km)~1, data=data.sp)
ck <- gstat(ck, "logAce", log(Acelaracion_cms2)~1, data=data.sp)
ck <- gstat(ck, "logMag", log(Magnitud_Mw)~1, data=data.sp)
ck

vm <- variogram(ck, width=1)
plot(vm)

pal = function(n = 9) brewer.pal(n, "Reds")
vm.fit <- fit.lmc(vm, ck, vgm(1, "Sph", 800, 1))
cok.maps <- predict(vm.fit,g)
cok.maps$cov.logProf.logAce
summary(cok.maps)
names(cok.maps)
print(spplot.vcov(cok.maps, cuts=6, col.regions=pal(7)))