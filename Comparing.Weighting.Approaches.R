# setwd()

library(maptools)
library(rgdal)
library(spdep)
library(RANN)

# sids<-readOGR(dsn="path",layer="sids2")

# we can see that the sids object is a SpatialPolygonsDataFrame
class(sids)
names(sids)

# create an initial contiguity-based weight matrix 
sids_nbq<-poly2nb(sids) # Queen based
sids_nbr<-poly2nb(sids, queen=FALSE) # Rook based

# set the coordinates of the sids SpatialPolygonsDataFrame for visualization and future weight matrix creation
coords<-coordinates(sids) # get centroids

plot(sids) 
plot(sids_nbq, coords, add=T) # all regions that share a border and a vertex are treated as one node
plot(sids) # clear the plot
plot(sids_nbr, coords, add=T) # regions that only share a border are treated as neighbors

# nearest neighbors approach to create connectivity matrices
IDs<-row.names(as(sids, "data.frame"))
sids_kn1<-knn2nb(knearneigh(coords, k=1), row.names=IDs)
sids_kn2<-knn2nb(knearneigh(coords, k=2), row.names=IDs)
sids_kn4<-knn2nb(knearneigh(coords, k=4), row.names=IDs)

plot(sids)
plot(sids_kn1, coords, add=T)
plot(sids)
plot(sids_kn2, coords, add=T)
plot(sids)
plot(sids_kn4, coords, add=T)

# distance-based matrices
dist<-unlist(nbdists(sids_kn1, coords))
summary(dist)

# set the max distance
max_k1<-max(dist)

# create different weight matrices
sids_kd1<-dnearneigh(coords, d1=0, d2=0.75*max_k1, row.names=IDs)
sids_kd2<-dnearneigh(coords, d1=0, d2=1*max_k1, row.names=IDs)
sids_kd3<-dnearneigh(coords, d1=0, d2=1.5*max_k1, row.names=IDs)

plot(sids)
plot(sids_kd1, coords, add=T)
plot(sids)
plot(sids_kd2, coords, add=T)
plot(sids)
plot(sids_kd3, coords, add=T)

# weight presentation
sids_nbq_w<-nb2listw(sids_nbq) # standard format
sids_nbq_w

sids_nbq_wb<-nb2listw(sids_nbq, style="B") # binary format
sids_nbq_wb

# inverse distance weighting
dist<-nbdists(sids_nbq, coordinates(sids))
idw<-lapply(dist, function(x) 1/(x/1000))
sids_nbq_idwb<-nb2listw(sids_nbq, glist=idw, style="B")
summary(unlist(sids_nbq_idwb$weights))

# Moran's I for global spatial dependence
moran.test(sids$SIDR79, listw=sids_nbq_w) # a Moran's I of zero indicates non spatial dependence

# simulate test of obtained Moran's I value (where it might fall)
set.seed(1234)
perm<-moran.mc(sids$SIDR79,listw=sids_nbq_w,nsim=999)
perm

mean(perm$res[1:999])

var(perm$res[1:999])

summary(perm$res[1:999])

hist(perm$res, freq=TRUE, breaks=20, xlab="Simulated Moran's I")

# local indicators of spatial association (LISA)
# examine the data for the geographic locations of counties that disproportionately 
# contribute to that autocorrelation via local clustering and spatial dependence
fips <- order(sids$FIPSNO)
nclocI <- localmoran(sids$SIDR79, sids_nbq_w)
printCoefmat(data.frame(nclocI[fips,],row.names=sids$FIPSNO[fips]), check.names=FALSE)

# plot the global Moran's I
nci <- moran.plot(sids$SIDR79,sids_nbq_w,labels=as.character(sids$NAME),xlim=c(-1,6.5),ylim=c(-1,4.5),xlab="SIDS Rate",ylab="SL SIDS Rate")

# create a map to show geographic distribution of the clusters associated with the LISA results
infl <- apply(nci$is.inf, 1, any)
x <- sids$SIDR79
lhx <- cut(x, breaks=c(min(x), mean(x), max(x)), labels=c("L","H"), include.lowest=TRUE)
wx <- lag(sids_nbq_w, sids$SIDR79)
lhwx <- cut(wx, breaks=c(min(wx), mean(wx), max(wx)),labels=c("L", "H"), include.lowest=TRUE)
lhlh <- interaction(lhx, lhwx, infl, drop=TRUE)
cols <- rep(1, length(lhlh))
cols[lhlh == "H.L.TRUE"] <- 2
cols[lhlh == "L.H.TRUE"] <- 3
cols[lhlh == "H.H.TRUE"] <- 4
plot(sids, col=grey.colors(4, 0.95, 0.55, 2.2)[cols])
legend("topright", legend=c("None", "HL", "LH", "HH"),fill=grey.colors(4, 0.95, 0.55, 2.2), bty="n", cex=0.8,y.intersp=0.8)

# from this we see that there is spatial dependence and spatial clustering
# need to control for these in subsequent analyses