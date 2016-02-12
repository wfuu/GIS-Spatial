## This script is a near replication of the satellite imagery data visualization implemented by
## the Commerce Data Usability Project at http://commercedataservice.github.io/tutorial_viirs_part1/

## Section 1: Leading off

library(doParallel)
library(foreach)
library(raster)
library(sp)
library(rgdal)
library(ggmap)
library(plotly)
##################

## Section 2: Loading data (VIIRS Day/Night Band Satellite Data)
## Step 1: Download raster file of VIIRS DNB monthly composite that excludes stray light 
## Step 2: Set directory path 
imagery <- "/Users/afu/Documents/School/QMSS II/Projects/NOAA replication/imagery"

## Step 3: Obtain a list of TIF files, load in the first file in list
tifs <- list.files(imagery, pattern = "\\.tif")
rast <- raster(paste(imagery,"/", tifs[1], sep=""))

## Step 4: Specify WGS84 as the projection of the raster file
wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
projection(rast) <- CRS(wgs84)

## Step 5: Download MSA shapefiles and use the same WGS84 projection. Use a function to download
shape_direct <- function(url, shp) {
        library(rgdal)
        temp = tempfile()
        download.file(url, temp) ## download the URL target to the temp file
        unzip(temp,exdir=getwd()) ## unzip that file
        return(readOGR(paste(shp,".shp",sep=""),shp))
}
msa <- shape_direct(url="http://www2.census.gov/geo/tiger/GENZ2014/shp/cb_2014_us_cbsa_20m.zip", 
                    shp= "cb_2014_us_cbsa_20m")
projection(msa) <- CRS(wgs84)

## Step 6: Load population data by MSA, preprocess (subset, reorder, convert) the data to only MSAs
msa_pop <- read.csv("http://www.census.gov/popest/data/metro/totals/2014/files/CBSA-EST2014-alldata.csv")
msa_pop <- msa_pop[msa_pop$LSAD=="Metropolitan Statistical Area",]
msa_pop <- msa_pop[order(msa_pop$POPESTIMATE2014),]
msa_pop$NAME <- as.character(msa_pop$NAME) 
##################

## Section 3: Mapping (only the top 10 most populous cities)
cities <- c("New York, NY", "Los Angeles, CA","Chicago, IL", "Houston, TX",
            "Philadelphia, PA", "Phoenix, AZ", "San Antonio, TX", "San Diego, CA",     
            "Dallas, TX", "San Jose, CA")

## Step 1: Set graph layout
par(mai=c(0,0,0,0), mfrow = c(2,5), bg='#001a4d', bty='n')

## Step 2: Loop through data to create a montage of city maps
coords <- data.frame()
for(i in 1:length(cities)) {
        
        ## Geocode each city name by coords
        temp_coord <- geocode(cities[i], source = "google")
        coords <- rbind(coords,temp_coord)
        
        ## Specify the spatial boundary/frame and crop the raster file by it
        e <- extent(temp_coord$lon - 1, temp_coord$lon + 1,
                    temp_coord$lat - 0.25, temp_coord$lat + 0.25)
        rc <- crop(rast, e)    
        
        ## Run k-means to find natural intervals within the radiance distribution
        ## Varying the number of clusters affects the radiance distribution for all cities
        sampled <- as.vector(rc)
        clusters <- 20 
        clust <- kmeans(sampled,clusters)$cluster
        combined <- as.data.frame(cbind(sampled,clust))
        brk <- sort(aggregate(combined[,1], list(combined[,2]), max)[,2])
        
        ## Plotting
        plot(rc, breaks=brk, col=colorRampPalette(c("#001a4d","#0066FF", "yellow"))(clusters), 
             legend=F,yaxt='n',xaxt='n',frame = F, asp= 2)
        text(temp_coord$lon ,temp_coord$lat + 0.15,
             substr(cities[i],1,regexpr(",",cities[i])-1), 
             col="white", cex=1.25)
        
        rm(combined)
}

## Step 3: Process raster data into vector data (for further analysis)
## Use a function to allow parallel processing
masq <- function(shp,rast,i){
        
        ## Extract one polygon based on index value i
        polygon <- shp[i,]           ## Extract one polygon
        extent <- extent(polygon)    ## Extract the polygon extent 
        
        ## Raster extract
        outer <- crop(rast, extent)  ## Extract raster by polygon extent
        inner <- mask(outer,polygon) ## Keeps values from raster extract that are within polygon
        
        ## Convert cropped raster into a vector
        ## Specify coordinates
        coords <- expand.grid(seq(extent@xmin,extent@xmax,(extent@xmax-extent@xmin)/(ncol(inner)-1)),
                              seq(extent@ymin,extent@ymax,(extent@ymax-extent@ymin)/(nrow(inner)-1)))
        ## Convert raster into vector
        data <- as.vector(inner)
        
        ## Package data in neat dataframe
        data <- cbind(as.character(shp@data$CBSAFP[i]),coords, data) 
        colnames(data)<-c("GEOID","lon","lat","avg_rad")  
        data <- data[!is.na(data$avg_rad),] 
        
        return(data)
}
##################

## Section 4: Comparing measures through visualizations
## 1. Radiance level by city (select 5)
## MSAs by GEOID
msa_list <- c(16180,19140,45820,42540,35620)
radiances <- data.frame() 

## Loop MSAs and organize their radiance levels
for(i in msa_list){
        print(i)
        
        # Extract MSA i polygon
        shp_temp <- msa[msa@data$GEOID==i,]
        
        # Extract MSA abbreviated name
        if(regexpr("-",as.character(shp_temp@data$NAME)[1])[1]==-1){
                loc = as.character(substr(as.character(shp_temp@data$NAME)[1],1,regexpr(",",as.character(shp_temp@data$NAME)[1])-1))
        } else{
                loc = as.character(substr(as.character(shp_temp@data$NAME)[1],1,regexpr("-",as.character(shp_temp@data$NAME)[1])-1))
        }
        
        # Extract the radiances, append to radiances placeholder
        rad <- masq(shp_temp,rast,1)$avg_rad 
        temp <- data.frame(loc = as.character(paste(loc,"(TNL = ",round(sum(rad),0),")",sep="")), avg_rad = rad) 
        radiances <- rbind(radiances,temp)
}

## Plotting
ggplot(radiances, aes(x=log(avg_rad))) +
        geom_histogram(position="identity", alpha=0.4) +
        facet_grid(. ~ loc)

## Styling: Remove all axes labels for style
x <- list(
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
)
y <- list(
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
) 

## Initiate a plotly graph without axes (interactive graph)
ggplotly()  %>% layout(xaxis=x, yaxis=y)

## 2. Total nighttime light by population
## Set up comparisons
registerDoParallel(cores=2)
extract <- foreach(i=1:nrow(msa@data), .combine=rbind) %dopar% {
        data <- masq(msa,rast,i)
        data.frame(GEOID = data$GEOID[1],sum = sum(data$avg_rad))
}
extract$GEOID <- as.numeric(as.character(extract$GEOID))

## Join in population data
joined <- merge(extract, msa_pop[,c("CBSA","NAME","POPESTIMATE2014")],by.x="GEOID",by.y="CBSA")
colnames(joined) <- c("GEOID","TNL","MSA","Population")

## Plotting
plot_ly(joined, 
        x = log(TNL), 
        y = log(Population), 
        text = paste("MSA: ", MSA),
        mode = "markers", 
        color = TNL,colors="PuOr")  %>% 
        layout(title="Total Nighttime Light vs. Population", showlegend = F)

