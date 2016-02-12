
#setwd()

install.packages('maps')
install.packages('maptools')
install.packages('rgdal')
install.packages('RColorBrewer')
install.packages('classInt')

#PROJECTIONS

library(maps) 

oldpar<-par()

world <- map("world", res=0)
str(world)
head(world$names)
plot(world)

states <- map("state", res=0)
str(states)
head(states$names)
plot(states)

#Convert object to SpatialLines 

library(maptools)

spworld <- map2SpatialLines(world, proj4string = CRS("+proj=longlat"))
spstates <- map2SpatialLines(states, proj4string = CRS("+proj=longlat"))

str(spworld, max.level=2)
str(spstates,max.level=2)

plot(spworld)
plot(spstates)

#Compare with LAEA projection
library(rgdal)

world.laea <- spTransform(spworld, CRS("+proj=laea +lat_0=0 +lon_0=0"))
states.laea <- spTransform(spstates, CRS("+proj=laea +lat_0=43.0758 +lon_0=-89.3976"))

par(mfrow = c(2, 2), pty = "s", cex.axis = 0.5)
plot(spworld, axes = T)
title(main = "Longitude and\nLatitude")
plot(world.laea, axes = T)
title(main = "Lambert Azimuthal\nEqual Area")
plot(spstates, axes = T)
title(main = "Longitude and\nLatitude")
plot(states.laea, axes = T)
title(main = "Lambert Azimuthal\nEqual Area")

#SPATIAL REFERENCING

par(oldpar)

map.states <- map("state", plot = T, fill = T, res=0)

list.names.states <- strsplit(map.states$names,":")
tail(list.names.states)

map.IDs <- sapply(list.names.states, function(x) x[1])
tail(map.IDs)

#Convert object to SpatialPolygons
polystates <- map2SpatialPolygons(map.states, IDs = map.IDs,proj4string = CRS("+proj=longlat"))

summary(polystates)

plot(polystates)

#Project object with LAEA
states.laea <- spTransform(polystates, CRS("+proj=laea +lat_0=43.0758 +lon_0=-89.3976"))
plot(states.laea)

#Merge spatial polygon ID to file
sp.IDs <- sapply(slot(states.laea, "polygons"), function(x) slot(x,"ID"))
tail(sp.IDs)

sat <- read.csv("sat.csv", stringsAsFactors = F,row.names = 1)
head(sat)

#Link data
states.sat <- SpatialPolygonsDataFrame(polystates,sat)
summary(states.sat)

states.sat.laea <- spTransform(states.sat, CRS("+proj=laea +lat_0=43.0758 +lon_0=-89.3976"))
plot(states.sat.laea)
summary(states.sat.laea)

#MAPPING
par(mfrow = c(1,1), pty = "s", cex.axis = 0.5)
library(RColorBrewer)
display.brewer.all()

library(classInt)

#Visualize distribution of SAT verbal
plotvar <- states.sat.laea$verbal
nclr <- 5
plotclr <- brewer.pal(nclr, "Greys")
plotclr
class <- classIntervals(plotvar, nclr, style = "quantile")
class
colcode <- findColours(class, plotclr, digits = 3)
colcode
plot(states.sat.laea, col = colcode)

#More visualization
plotclr <- brewer.pal(nclr, "Purples")
class <- classIntervals(plotvar, nclr, style = "quantile")
colcode <- findColours(class, plotclr, digits = 3)
plot(states.sat.laea, col = colcode, border = "grey",axes = T)
title(main = "SAT math scores in 1999")
legend("bottomleft", legend = names(attr(colcode,"table")), fill = attr(colcode, "palette"))

#Write to .shp file for further visualization inside QGIS 
writeSpatialShape(states.sat.laea, "sat")