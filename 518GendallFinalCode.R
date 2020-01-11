#####################
#PACKAGES FOR FINAL PROJECT
#####################

#Install Libraries
install.packages("raster")
install.packages("ncdf4")
install.packages("sp")
install.packages("spatialEco")
install.packages("spgwr")
install.packages("lmtest")
install.packages("gvlma")
#####
#Load Libraries
library(ncdf4) #reading netcdf files
library(rgdal) #geospatial files

library(sp) #spatial data analysis tools
library(spdep) #gwr and weighting schemes
library(spgwr) #gwr
library(spatialEco) #spatial analysis and modelling
library(GISTools) #spatial data analysis tools
library(raster) #rasters and spatial data analysis tools
library(ape) #spatial autocorrelation
library(gstat) #spatial autocorrelation (variogram)

library(plyr) #data manipulation
library(dplyr) #data manipulation

library(maptools) #tools for mapping spatial onjects
library(tmap) #mapping
library(BAMMtools) #pallete explorer?
library(shinyjs) #pallete explorer?

library(gtable) #table

library("grid") #grid graphic
library("gridExtra") #extension for grid graphics

library(lmtest) #for BP test for heteroskedasticity of OLS residuals
library(gvlma) #tests OLS assumptions are met or not
#########setting working directory
dir <- "~/Desktop/GEOG518/FinalProject"
setwd(dir)





#############################
#DATA PREP & DESCRIPTIVE STATS
#############################

##############
#Kelp data prep
##############

#Load kelp shapefile into usable format (from Washington State DNR, department of natural ressources)
kelp <- readOGR(dsn = "./Data" , layer = "2014kelp")
#readOGR function is used to read open source geospatial data into a suitable spatial vector object

#shows Spatial Polygons Data Frame
View(kelp)

####################
#summary statistics of kelp SST area
summary(kelp$Shape_Area) 

#individually computed summary descriptive stats for table:
#mean
meanBed <- mean(kelp$Shape_Area, na.rm = TRUE) #Use na.rm = TRUE to ignore NA values in calculation
meanBed <- round(meanBed, digits = 3)
meanBed

#standard deviation
sdBed <- sd(kelp$Shape_Area, na.rm = TRUE) #Calculate the SD, ignoring NA values
sdBed <- round(sdBed, digits = 3)
sdBed

#Mode
modeBed <- as.numeric(names(sort(table(kelp$Shape_Area), decreasing = TRUE))[1]) #make frequency table of bed area variable and sort it in desending order and extract the first row (Most Frequent)
modeBed
#Bed area is unique to each polygon, continous variable and non repeat

#Median
medBed <- median(kelp$Shape_Area, na.rm = TRUE)
medBed <- round(medBed, digits = 3)
medBed

#skew
skewBed <- e1071::skewness(kelp$Shape_Area, na.rm = TRUE)
skewBed <- round(skewBed, digits = 3)
skewBed

#Kurtosis
kurtBed <- e1071::kurtosis(kelp$Shape_Area, na.rm = TRUE)
kurtBed <- round(kurtBed, digits = 3)
kurtBed

#CoV 
CoVBed <- (sdBed / meanBed) * 100
CoVBed <- round(CoVBed, digits = 3)
CoVBed

#Normal distribution test
normBed_PVAL <- shapiro.test(kelp$Shape_Area)$p.value
normBed_PVAL <- round(normBed_PVAL, digits = 62)
normBed_PVAL

#Normal distribution test on log transformed bed area
normBedL_PVAL <- shapiro.test(log10(kelp$Shape_Area))$p.value
normBedL_PVAL <- round(normBedL_PVAL, digits = 28)
normBedL_PVAL


####################
#Histogram of bed area
hist(kelp$Shape_Area, breaks = 100, main = " ", xlab = "Kelp Bed Area (m^2)", las = 1, ylim= c(0, 1200), xlim= c(0, 14000000)) 

png("./Outputs/2014KelpBasicHist.png") #Create an object to print the table to
hist(kelp$Shape_Area, breaks = 100, main = " ", xlab = "Kelp Bed Area (m^2)", las = 1, ylim= c(0, 1200), xlim= c(0, 14000000)) 
dev.off() #Print table

#Histogram of log transformed kelp bed area
hist(log10(kelp$Shape_Area), breaks = 100, main= " ", xlab = "Log of Kelp Bed Area (m^2)", las = 1) #Base R style

png("./Outputs/2014KelpLogHist.png") #Create an object to print the table to
hist(log10(kelp$Shape_Area), breaks = 100, main= " ", xlab = "Log of Kelp Bed Area (m^2)", las = 1) #Base R style
dev.off() #Print table

######################
#transforming kelp to the same projection as SST for next section

#transform kelp projection from NAD 83 HARN (EPSG 2926) to WGS 84
kelpP <- spTransform(kelp, CRS("+init=epsg:4326"))

################################
#SST data prep
################################

#loading netcdf SST 4km world file for august 2014 (downloaded from NASA Earth Data from MODIS-terra)
nc.brick <-brick("./Data/T20140322014059.L3m_MO_SST_sst_4km.nc")
dim(nc.brick)
#plot(nc.brick)

#make Netcdf file into dataframe
SST.df <- as.data.frame(nc.brick[[1]], xy=T)

#viewing data
# head(SST.df)
# View(SST.df)

#Subset
###############
SST.df <- subset(SST.df, x >= -125 & x <= -122)
SST.df <- subset(SST.df, y >= 47.5 & y <= 48.5)

#remove NA values
SST.df <- na.omit(SST.df)

#maked dataframe into a csv
write.csv(SST.df, file.choose())

#get average SST for the month of August 2014 for the study area around Washington State
mean(SST.df$layer, na.rm=T)

#############
#SST dataframe to spatial points dataframe
#############

#setting coordinates using x y columns in dataframe 
Coords <- SST.df[,c("x","y")]

#creating spatialpointdataframe using coordinates, dataframe of SST for data and the CRS WGS 84 (because of lat long and netcdf file format)
SST <- SpatialPointsDataFrame(coords = Coords, data = SST.df, proj4string = crs("+init=epsg:4326"))

################
#map of SST
################
map_SST <-  tmap_mode("plot") +
  tm_shape(SST) + 
  tm_symbols(size=.05, 
             col= "layer",
             palette = "-RdBu") +
  tm_shape(kelpP) + 
  tm_borders(col = NA, lwd = .4, lty = "solid", alpha = NA,
             group = NA) +
  tm_scale_bar()

map_SST



######################
#Spatial interpolation of SST to enchance resolution (IDW)
######################

###IDW OR KRIGING? Look at variogram & tried kriging multiple times R crashed everytime
head(SST)

###Look at semivariogram 
f.0 <- as.formula(layer ~ 1) 
#No trend, 0 order polynomial

var.smpl <- variogram(f.0, SST, cloud = FALSE) 

dat.fit  <- fit.variogram(var.smpl, fit.ranges = TRUE, fit.sills = TRUE,
                          vgm(model="Exp")) #gaussian gau, spherical sph, exponential exp
plot(var.smpl, dat.fit)


#IDW
# Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(SST, "regular", n=100000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object

proj4string(grd) <- proj4string(SST)
P.idw <- gstat::idw(layer ~ 1, SST, newdata=grd, idp=10)
r       <- raster(P.idw)

#saving a raster vs a shapefile
writeRaster(r, "raster_SST_test4.tif", format = "GTiff")
writeOGR(SST, dsn = ".", layer = "SSTTest", driver = "ESRI Shapefile")

##ploting the raster made with idw two ways
plot(r)

#mapping IDW with tmap
map_SSTIDW<- tmap_mode("plot") +
  tm_shape(r) + 
  tm_raster(n=10,palette = "-RdBu",
            title="SST") +
  tm_legend(legend.outside=TRUE) +
  tm_shape(kelpSST) + 
  tm_borders(col = NA, lwd = .4, lty = "solid", alpha = NA,
             group = NA) +
  tm_scale_bar()

map_SSTIDW

####################
#Validating SST IDW

# Leave-one-out validation routine
IDW.out <- vector(length = length(SST))
for (i in 1:length(SST)) {
  IDW.out[i] <- gstat::idw(layer ~ 1, SST[-i,], SST[i,], idp=10)$var1.pred
}
#change idp for changng the power

# Plot the differences
OP <- par(pty="s", mar=c(4,3,0,0))
par(mar = c(5.1,4.1,3.1,2.1))
plot(IDW.out ~ SST$layer, asp=1, xlab="Observed", ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5))
abline(lm(IDW.out ~ SST$layer), col="red", lw=2,lty=2)
abline(0,1)
par(OP)
rtmean <- sqrt( sum((IDW.out - SST$layer)^2) / length(SST))

####################
##### Table of Rootmean square error for IDW power decision

#2 is 0.224180
#5 is 0.159108
#8 is 0.152287
#9 is 0.152049
#10 is 0.15199657 (10 it is)
#11 is 0.1520064
#10.5 is 0.15199758

Power = c(2,5,8,9,10,11)
RMSE = c(0.2242, 0.1591, 0.1523, 0.1521, 0.1519, 0.1520) 


data.for.tableIDW = data.frame(Power,RMSE)

#Make table 3
tableIDW <- tableGrob(data.for.tableIDW, rows = rownames(data.for.tableIDW)) #make a table "Graphical Object" (GrOb) 
tCaption <- textGrob("IDW", gp = gpar(fontsize = 09))
padding <- unit(5, "mm") 

tableIDW <- gtable_add_rows(tableIDW, 
                          heights = grobHeight(tCaption) + padding, 
                          pos = 0) 

tableIDW <- gtable_add_grob(tableIDW,
                          tCaption, t = 1, l = 2, r = ncol(data.for.tableIDW) + 1)
grid.arrange(tableIDW, newpage = TRUE)

#export table to png
png("./Outputs/IDW.png") #Create an object to print the table to
grid.arrange(tableIDW, newpage = TRUE)
dev.off() #Print table





###################################
#Putting kelp & SST together
##################################

plot(r) #ploting the raster
plot(kelpP, add=T)

#joining the kelp with SST
kelpSST <- extract(r, kelpP, fun=mean, sp = TRUE)

#################################
#Interactive map of Kelp polygons with associated temperature values
map_kelpSST <-  tmap_mode("view") +
  tm_shape(kelpSST) + 
  tm_polygons(col = "var1.pred", 
              style= "jenks",
              title = "Temperature ˚C", 
              palette = "-RdBu")+
  tm_scale_bar()

map_kelpSST

tmap_save(map_kelpSST, filename = "kelpSST_map.html")

####################################
#summary statistics of kelp bed SST

summary(kelpSST@data$var1.pred) 

#####################################
#individually computed summary descriptive stats for table

#mean
meanSST <- mean(kelpSST@data$var1.pred, na.rm = TRUE) #Use na.rm = TRUE to ignore NA values in calculation
meanSST <- round(meanSST, digits = 3)
meanSST

#standard deviation
sdSST <- sd(kelpSST@data$var1.pred, na.rm = TRUE) #Calculate the SD, ignoring NA values
sdSST <- round(sdSST, digits = 3)
sdSST

#Mode
modeSST <- as.numeric(names(sort(table(kelpSST@data$var1.pred), decreasing = TRUE))[1]) #make frequency table of SST area variable and sort it in desending order and extract the first row (Most Frequent)
modeSST
#SST area is unique to each polygon, continous variable and non repeat

#Median
medSST <- median(kelpSST@data$var1.pred, na.rm = TRUE)
medSST <- round(medSST, digits = 3)
medSST

#skew
skewSST <- e1071::skewness(kelpSST@data$var1.pred, na.rm = TRUE)
skewSST <- round(skewSST, digits = 3)
skewSST

#Kurtosis
kurtSST <- e1071::kurtosis(kelpSST@data$var1.pred, na.rm = TRUE)
kurtSST <- round(kurtSST, digits = 3)
kurtSST

#CoV 
CoVSST <- (sdSST / meanSST) * 100
CoVSST <- round(CoVSST, digits = 3)
CoVSST

#Normal distribution test
normSST_PVAL <- shapiro.test(kelpSST@data$var1.pred)$p.value
normSST_PVAL <- round(normSST_PVAL, digits = 54)
normSST_PVAL

############
#Making a descriptive stats table
##########

#Create a table of descriptive stats
Variable = c("Kelp Bed Area ft^2" , "Sea Surface Temperature ˚C")
Mean = c(meanBed, meanSST) #Create an object for the means
SD = c(sdBed, sdSST) #Create an object for the standard deviations
Median = c(medBed, medSST) #Create an object for the medians
Mode <- c(modeBed, modeSST) #Create an object for the modes
Skewness <- c(skewBed, skewSST) #Create an object for the skewness
Kurtosis <- c(kurtBed, kurtSST) #Create an object for the kurtosis
CoV <- c(CoVBed, CoVSST) #Create an object for the CoV
Normality <- c(normBed_PVAL, normSST_PVAL) #Create an object for the normality PVALUE
L.Normality <- c(normBedL_PVAL, "-") #Create an object for the normality for log transformed bed area PVALUE

data.for.table1 = data.frame(Variable, Mean, SD, Median, Skewness, Kurtosis, CoV, Normality, L.Normality)

#Make table 1
table1 <- tableGrob(data.for.table1, rows = rownames(data.for.table1)) #make a table "Graphical Object" (GrOb) 
t1Caption <- textGrob("Table 1: Descriptive Statistics", gp = gpar(fontsize = 09))
padding <- unit(5, "mm") 

table1 <- gtable_add_rows(table1, 
                          heights = grobHeight(t1Caption) + padding, 
                          pos = 0) 

table1 <- gtable_add_grob(table1,
                          t1Caption, t = 1, l = 2, r = ncol(data.for.table1) + 1)
grid.arrange(table1, newpage = TRUE)

#export table to png
png("./Outputs/DescriptiveStats.png") #Create an object to print the table to
grid.arrange(table1, newpage = TRUE)
dev.off() #Print table


############
#Histogram of SST
############

par(mar= c(4.1,1,1,1))
hist(kelpSST@data$var1.pred, breaks = 100, main = " ", xlab = "Kelp Sea Surface Temperature ˚C") #Base R style

png("./Outputs/2014SSTBasicHist.png") #Create an object to print the table to
hist(kelpSST@data$var1.pred, breaks = 100, main = " ", xlab = "Kelp Sea Surface Temperature ˚C") #Base R style
dev.off() #Print table







####################################
#LINEAR REGRESSION KELP & SST
#####################################

#regular linear regression
plot(kelpSST$Shape_Area~kelpSST$var1.pred)

#Perform a linear regression on the two variables. You should decide which one is dependent.
lm.model <- lm(kelpSST$Shape_Area~kelpSST$var1.pred)
#Add the regression model to the plot you created
abline(lm.model)
#Get the summary of the results
summary(lm.model)

gvlma(lm.model)

AIC(lm.model)

############
#log transformed linear regression
###########

plot(log10(kelpSST$Shape_Area)~kelpSST$var1.pred)
#Perform a linear regression on the two variables. You should decide which one is dependent.
lm.model <- lm(log10(kelpSST$Shape_Area)~kelpSST$var1.pred)
#Add the regression model to the plot you created
abline(lm.model)
#Get the summary of the results
summary(lm.model)

gvlma(lm.model)

AIC(lm.model)


####################################################
#TESTING LINEAR REGRESSION ASSUMPTIONS FOR RESIDUALS
#####################################################

#mapping and plotting residuals 
####################
#You want to determine if the model residuals are spatially clustered. 
#First obtain the residuals from the model
model.resids <- as.data.frame(residuals.lm(lm.model))
#Then add the residuals to your spatialpolygon dataframe
kelpSST$residuals <- residuals.lm(lm.model)
#Observe the result to make sure it looks correct
head(kelpSST)

#interactive map of residuals 
map_kelpresids <-  tmap_mode("view") +
  tm_shape(kelpSST) + 
  tm_polygons(col = "residuals", 
              style= "jenks",
              title = "Residuals", 
              palette = "magma")+
  tm_scale_bar()

map_kelpresids

tmap_save(map_kelpresids, filename = "LRresiduals_map.html")

#####plot the residuals 
plot(kelpSST$residuals~kelpSST$var1.pred)
lm.model <- lm(kelpSST$residuals~kelpSST$var1.pred)
abline(lm.model)

plot(log10(kelpSST$Shape_Area)~kelpSST$residuals)
lm.model <- lm(log10(kelpSST$Shape_Area)~kelpSST$residuals)
abline(lm.model)


#testing residuals for heteroskedasticity with Breusch Pagan Test (using the variance function of residuals and a chi-squared test)
############################
qchisq(.95, df = 1) 

bptest(lm.model)


#testing residuals for normality
####################
#Histogram of residuals
hist(kelpSST@data$residuals, breaks = 100, main = " ", xlab = "Kelp Residuals") #Base R style

png("./Outputs/2014RESBasicHist.png") #Create an object to print the table to
hist(kelpSST@data$residuals, breaks = 100, main = " ", xlab = "Kelp Residuals") #Base R style
dev.off() #Print table

#Normal distribution test for residualts
normRES_PVAL <- shapiro.test(kelpSST@data$residuals)$p.value
normRES_PVAL <- round(normRES_PVAL, digits = 25)
normRES_PVAL









#########################################
#TESTING FOR AUTOCORRELATION (Moran's I & Variograms) OF VARIABLES & RESIDUALS
#########################################

####Make into dataframe with lat & long (DONE LOWER DOWN IN GWR)
kelpSST.coords <- sp::coordinates(kelpSST)
kelpSST$X <- kelpSST.coords[,1]
kelpSST$Y <- kelpSST.coords[,2]
head(kelpSST)

#####Making a matrix of inverse distance weights with points (not polys)
Distances <- as.matrix(dist(cbind(kelpSST$X, kelpSST$Y)))

#inversing
Distances.inv <- 1/Distances

#transforming diagonal (the distance of one point to itself) back to zero after the inversing
diag(Distances.inv) <- 0

#viewing first 5 rows and column so computer doesn't crash
Distances.inv[1:5, 1:5]

#running moran's I
Moran.I(kelpSST@data$Shape_Area, Distances.inv)

Moran.I(kelpSST@data$var1.pred, Distances.inv)

Moran.I(kelpSST@data$residuals, Distances.inv)

#################
#table for autocorrelation results
##################

Variable = c("Kelp Bed Area ft^2" , "Sea Surface Temperature ˚C", "OLS Residuals", "Log OLS Residuals ")
Observed = c(0.0183, 0.662, 0.0170, 0.0986) 
SD = c(0.00484, 0.00533,0.00483, 0.00533) 
Expected = c(-0.000692, -0.000692,-0.000692, -0.000692) 
PVAL <- c("8.710e-05","0","2.605e-4","0") 

data.for.table3 = data.frame(Variable, Observed, SD, Expected, PVAL)

#Make table 3
table3 <- tableGrob(data.for.table3, rows = rownames(data.for.table3)) #make a table "Graphical Object" (GrOb) 
t3Caption <- textGrob("Table 1: Moran's I results", gp = gpar(fontsize = 09))
padding <- unit(5, "mm") 

table3 <- gtable_add_rows(table3, 
                          heights = grobHeight(t3Caption) + padding, 
                          pos = 0) 

table3 <- gtable_add_grob(table3,
                          t3Caption, t = 1, l = 2, r = ncol(data.for.table3) + 1)
grid.arrange(table3, newpage = TRUE)

#export table to png
png("./Outputs/Moran.png") #Create an object to print the table to
grid.arrange(table3, newpage = TRUE)
dev.off() #Print table

###################
#Semivariograms looking at autocorrelation 
#SST
head(kelpSST)

f.0 <- as.formula(var1.pred ~ 1) 
#No trend, 0 order polynomial

var.smpl <- variogram(f.0, kelpSST, cloud = FALSE) 

dat.fit  <- fit.variogram(var.smpl, fit.ranges = TRUE, fit.sills = TRUE,
                          vgm(model="Exp")) #gaussian gau, spherical sph, exponential exp
plot(var.smpl, dat.fit)

#####################
#VARIOGRAM for AREA

f.0 <- as.formula(log10(Shape_Area) ~ 1) 
#No trend, 0 order polynomial

var.smpl <- variogram(f.0, kelpSST, cloud = FALSE) 

dat.fit  <- fit.variogram(var.smpl, fit.ranges = TRUE, fit.sills = TRUE,
                          vgm(model="Sph")) #gaussian gau, spherical sph, exponential exp
plot(var.smpl, dat.fit)


#####################
#VARIOGRAM for RESIDUALS

f.0 <- as.formula(residuals ~ 1) 
#No trend, 0 order polynomial

var.smpl <- variogram(f.0, kelpSST, cloud = FALSE) 

dat.fit  <- fit.variogram(var.smpl, fit.ranges = TRUE, fit.sills = TRUE,
                          vgm(model="Sph")) #gaussian gau, spherical sph, exponential exp
plot(var.smpl, dat.fit)












###NOT USED IN PAPER
#######################################
#GEOGRAPHIC WEIGHTED REGRESSION
#######################################

#Let's say you are continuing with your data from the regression analysis. 

#The first thing you need to do is to add the polygon coordinates to the spatialpolygondataframe.
#You can obtain the coordinates using the "coordinates" function from the sp library
kelpSST.coords <- sp::coordinates(kelpSST)
#Observe the result
head(kelpSST.coords)
#Now add the coordinates back to the spatialpolygondataframe
kelpSST$X <- kelpSST.coords[,1]
kelpSST$Y <- kelpSST.coords[,2]
head(kelpSST)

####FIXED
###Determine the bandwidth for GWR: this will take a while
GWRbandwidth <- gwr.sel(log10(kelpSST@data$Shape_Area)~kelpSST@data$var1.pred, 
                        data=kelpSST, coords=cbind(kelpSST$X,kelpSST$Y),adapt=F) 

###Perform GWR on the two variables with the bandwidth determined above
###This will take a looooooong while
gwr.model = gwr(log10(kelpSST$Shape_Area)~kelpSST$var1.pred, 
                data=kelpSST, coords=cbind(kelpSST$X,kelpSST$Y), 
                bandwidth = GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

#Print the results of the model
gwr.model

######ADAPTIVE
###Determine the bandwidth for GWR: this will take a while
GWRbandwidth <- gwr.sel(log10(kelpSST@data$Shape_Area)~kelpSST@data$var1.pred, 
                        data=kelpSST, coords=cbind(kelpSST$X,kelpSST$Y),adapt=T) 

###Perform GWR on the two variables with the bandwidth determined above
###This will take a looooooong while
gwr.model = gwr(log10(kelpSST$Shape_Area)~kelpSST$var1.pred, 
                data=kelpSST, coords=cbind(kelpSST$X,kelpSST$Y), 
                adapt = GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 


#Look at the results in detail
results<-as.data.frame(gwr.model$SDF)
head(results)

#Now for the magic. Let's add our local r-square values to the map
kelpSST$localr <- results$localR2

# #Create choropleth map of r-square values
# local.r.square <- kelpSST$localr
# shades <- auto.shading(local.r.square, n=6, cols = brewer.pal(6, 'Oranges'))
# choropleth(kelpSST, local.r.square, shades) #map the data with associated colours
# choro.legend(4769000, 1248500, shades) #add a legend (you might need to change the location)

map_kelprsq <-  tmap_mode("view") +
  tm_shape(kelpSST) + 
  tm_polygons(col = "localr", 
              style= "jenks",
              title = "R squared", 
              palette = "magma")+
  tm_scale_bar()

map_kelprsq

tmap_save(map_kelprsq, filename = "gwr_rsquared.html")

#Time for more magic. Let's map the coefficients
kelpSST$coeff <- results$kelpSST.var1.pred
head(kelpSST)

#Create choropleth map of the coefficients
# local.coefficient <- pm.income.poly$coeff
# shades <- auto.shading(local.coefficient, n=6, cols = brewer.pal(6, 'Oranges'))
# choropleth(income.tracts, local.coefficient, shades) #map the data with associated colours
# choro.legend(3864000, 1965000, shades) #add a legend (you might need to change the location)

map_kelpcoef <-  tmap_mode("view") +
  tm_shape(kelpSST) + 
  tm_polygons(col = "coeff", 
              style= "jenks",
              title = "Coefficient",
              midpoint = 0, 
              palette = "magma")+
  tm_scale_bar()

map_kelpcoef

tmap_save(map_kelpcoef, filename = "gwr_coefficient.html")

