# Choose Agro_Data.csv
Agro<-read.csv("data/Agro_Data.csv")
xy <- Agro[2:3]
df <- Agro[-2:-3]
AgroSP <- SpatialPointsDataFrame(coords=xy, data=df)
plot(AgroSP,pch=1,cex=0.5)

# Names of the variables, soil pH("pH"), Soil organic matter ("OM"), Phosphorus ("P"), Potassium ("K"), Magnesium ("Mg"), Calcium ("Ca"), Aluminum ("Al"), Nitrate("NO_3"), Electrical conductivity("EC"), Elevation (meters) ("Elev"). 


# Choose Agro_Mask.tif
Agro_mask<-readGDAL("data/Agro_mask.tif")
proj4string(AgroSP)<-proj4string(Agro_mask)
Agro_mask@data[Agro_mask@data==2]=1
Agro_mask@data[Agro_mask@data==0]=NA

spplot(AgroSP['Ca'])
spplot(AgroSP['K'])

lx<-diff(range(coordinates(AgroSP)[,1])) #Approximate field width
ly<-diff(range(coordinates(AgroSP)[,2])) #Approximate field length

# You can continue the analysis here... 

## Variogram cloud
vario.claud <- variogram(Sand~1, ArboSP, cloud=TRUE, curoff=350)
