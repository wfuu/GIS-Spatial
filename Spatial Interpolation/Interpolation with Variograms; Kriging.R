
########Read .shp file into R#############
# setwd()

library(sp)
library(rgdal)
BMcD <- readOGR(".", "BMcD")
BMcD$Fldf <- factor(BMcD$Fldf)
names(BMcD)
spTransform(BMcD, CRS("+proj=laea +lat_0=51.998807 +lon_0=5.692291"))


#######Plot the coordinates##############
plot(BMcD)
bubble(BMcD, "Zn")


########Explore statistical data#######################
boxplot(Zn ~ Fldf, BMcD, width=table(BMcD$Fldf), col="grey")


######Create SptialPixelsDataFrame#################
BMcD_grid <- as(readGDAL("BMcD_fldf.txt"), "SpatialPixelsDataFrame")
names(BMcD_grid) <- "Fldf"
BMcD_grid$Fldf <- as.factor(BMcD_grid$Fldf)



######Plot points and grid overlay################
pts = list("sp.points", BMcD, pch = 4, col = "white")
spplot(BMcD_grid, "Fldf", col.regions=1:3, sp.layout=list(pts))


##########Set parameters for interpolated plots#######
bluepal <- colorRampPalette(c("azure1", "steelblue4"))
brks <- c(0,130,155,195,250,330,450,630,890,1270,1850)
cols <- bluepal(length(brks)-1)
scols <- c("green", "red")



#######Examine residuals and significance differences in Zinc (ANOVA)######
library(ipred)
res <- errorest(Zn ~ 1, data = as(BMcD, "data.frame"), model=lm, est.para=control.errorest(k=nrow(BMcD), random=FALSE, predictions=TRUE))
round(res$error, 2)
fres <- lm(Zn ~ Fldf, data=BMcD)
anova(fres)
eres <- errorest(Zn ~ Fldf, data = as(BMcD, "data.frame"), model=lm, est.para=control.errorest(k=nrow(BMcD), random=FALSE, predictions=TRUE))
round(eres$error, 2)


####Estimate an aspatail model###########
library(maptools)
BMcD_grid$lm_pred <- predict(fres, newdata=BMcD_grid)
image(BMcD_grid, "lm_pred", breaks=brks, col=cols)
title("Flood frequency model interpolation")
pe <- BMcD$Zn-eres$predictions
symbols(coordinates(BMcD), circles=sqrt(abs(pe)), fg="black", bg=scols[(pe < 0)+1], inches=FALSE, add=TRUE)
legend("topleft", fill=cols, legend=leglabs(brks), bty="n", cex=0.8)


####Thin Plate Spline interpolation#######
library(fields)
pe_tps <- numeric(nrow(BMcD))
cBMcD <- coordinates(BMcD)
for (i in seq(along=pe_tps)) {
  tpsi <- Tps(cBMcD[-i,], BMcD$Zn[-i])
  pri <- predict(tpsi, cBMcD[i,,drop=FALSE])
  pe_tps[i] <- BMcD$Zn[i]-pri
}
round(sqrt(mean(pe_tps^2)), 2)
tps <- Tps(coordinates(BMcD), BMcD$Zn)


#####Plot spline results########
BMcD_grid$spl_pred <- predict(tps, coordinates(BMcD_grid))
image(BMcD_grid, "spl_pred", breaks=brks, col=cols)
title("Thin plate spline model")
symbols(coordinates(BMcD), circles=sqrt(abs(pe_tps)), fg="black", bg=scols[(pe_tps < 0)+1], inches=FALSE, add=TRUE)
legend("topleft", fill=cols, legend=leglabs(brks), bty="n", cex=0.8)


#####Estimate spatial dependence in data##########
library(gstat)
cvgm <- variogram(Zn~1, data=BMcD, width=100, cutoff=1000)
efitted <- fit.variogram(cvgm, vgm(psill=1, model="Exp", range=100, nugget=1))
efitted


###Plot semivariogram#######
plot(cvgm, model=efitted, plot.numbers=TRUE, col="black")


####Interpolate using original kriging############
OK_fit <- gstat(id="OK_fit", formula = Zn ~ 1, data = BMcD, model=efitted)
pe <- gstat.cv(OK_fit, debug.level=0, random=FALSE)$residual
round(sqrt(mean(pe^2)), 2)
z <- predict(OK_fit, newdata=BMcD_grid, debug.level=0)
BMcD_grid$OK_pred <- z$OK_fit.pred
BMcD_grid$OK_se <- sqrt(z$OK_fit.var)


####Plot OK model##################
image(BMcD_grid, "OK_pred", breaks=brks, col=cols)
title("Fitted exponential OK model")
symbols(coordinates(BMcD), circles=sqrt(abs(pe)), fg="black", bg=scols[(pe < 0)+1], inches=FALSE, add=TRUE)
legend("topleft", fill=cols, legend=leglabs(brks), bty="n", cex=0.8)


#########Add flood frequency to the semivariance estimation############
cvgm <- variogram(Zn~Fldf, data=BMcD, width=100, cutoff=1000)
uefitted <- fit.variogram(cvgm, vgm(psill=1, model="Exp", range=100, nugget=1))
uefitted


#####Plot universial kriging model#########
plot(cvgm, model=uefitted, plot.numbers=TRUE, col="black")


####Interpolate using the flood frequency UK model########
UK_fit <- gstat(id="UK_fit", formula = Zn ~ Fldf, data = BMcD, model=uefitted)
pe_UK <- gstat.cv(UK_fit, debug.level=0, random=FALSE)$residual
round(sqrt(mean(pe_UK^2)), 2)
z <- predict(UK_fit, newdata=BMcD_grid, debug.level=0)
BMcD_grid$UK_pred <- z$UK_fit.pred
BMcD_grid$UK_se <- sqrt(z$UK_fit.var)


######Plot UK model#########
image(BMcD_grid, "UK_pred", breaks=brks, col=cols)
title("Flood frequency UK model")
symbols(coordinates(BMcD), circles=sqrt(abs(pe_UK)), fg="black", bg=scols[(pe_UK < 0)+1], inches=FALSE, add=TRUE)
legend("topleft", fill=cols, legend=leglabs(brks), bty="n", cex=0.8)





#####Plot all interpolation for comparison######
pts = list("sp.points", BMcD, pch = 4, col = "black", cex=0.5)
spplot(BMcD_grid, c("lm_pred", "spl_pred", "OK_pred", "UK_pred"), at=brks, col.regions=cols, sp.layout=list(pts))


#####Write raster to TIF#####
writeGDAL(BMcD_grid["UK_pred"], "UK_pred.tif")




