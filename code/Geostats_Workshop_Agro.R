# Choose Agro_Data.csv
Agro<-read.csv("data/Agro_Data.csv")
xy <- Agro[2:3]
df <- Agro[-2:-3]
AgroSP <- SpatialPointsDataFrame(coords=xy, data=df)
plot(AgroSP,pch=1,cex=0.5)

# Names of the variables, soil pH("pH"), Soil organic matter ("OM"), Phosphorus ("P"), Potassium ("K"), Magnesium ("Mg"), Calcium ("Ca"), Aluminum ("Al"), Nitrate("NO_3"), Electrical conductivity("EC"), Elevation (meters) ("Elev"). 


# Choose Agro_Mask.tif
Agro_mask<-readGDAL("data/Agro_mask.tif") # read raster
proj4string(AgroSP)<-proj4string(Agro_mask) # adjust projection
Agro_mask@data[Agro_mask@data==2]=1 # set binary data
Agro_mask@data[Agro_mask@data==0]=NA # replace zeros by NA

spplot(AgroSP['Ca']) # plot distribution of Ca values
spplot(AgroSP['K']) # plot distribution of K values

lx<-diff(range(coordinates(AgroSP)[,1])) #Approximate field width
ly<-diff(range(coordinates(AgroSP)[,2])) #Approximate field length


##============================================================================
## PART 2 - EXPERIMENTAL VARIOGRAM AND AUTOCORRELATION MEASURES
##============================================================================

# Semi-variogram for % Sand. 
# Variogram cloud
vario.SandCloud <- variogram(Ca~1,AgroSP,cloud=TRUE, cutoff=350)
plot(vario.SandCloud) # plot data variation at different distances

# Average at each lag ('classical' variogram)
vario.Sand <- variogram(Ca~1,AgroSP, cutoff=350, width=50)
plot(vario.Sand)

# Variogram for total Calcium
vario.tC <- variogram(Ca~1,AgroSP, cutoff=350, width=50)
plot(vario.tC)

# Choose the file GearyC.R
source("code/GearyC.R")
# Choose the file MoranI.R
source("code/MoranI.R")

# Geary C
corgramC<-GearyC(cbind(AgroSP$X,AgroSP$Y),AgroSP$Ca,seq(0,350,50))
plot(corgramC) # Il y a d'autocorrelation espaciale dès le quatrieme classe, 173.98m

# Moran's I
corgramI<-MoranI(cbind(AgroSP$X,AgroSP$Y),AgroSP$Ca,seq(0,350,50))
plot(corgramI)


##====================================================
## PART 3 - UNCONDITIONAL SIMULATIONS
##====================================================

# Unconditional simulations
# Pure nugget (no spatial structure)
v <- vgm(1,'Nug',0) # create nugget model; the nugget represents the small-scale variability of the data, and it's the y-intercept of the variogram
g.dummy <- gstat(formula = z~1, locations=AgroSP, dummy = TRUE, beta = 0,model = v, nmax = 20)
g.prdNug<-predict(g.dummy,Agro_mask,nsim=9)
spplot(g.prdNug)

# Spherical model with 173.98 m range. 
v <- vgm(1, "Sph", 173.98)
g.dummy <- gstat(formula = z~1, locations=AgroSP, dummy = TRUE, beta = 0,model = v, nmax = 20)
g.prdSph<-predict(g.dummy,Agro_mask,nsim=9)
spplot(g.prdSph)

# Gaussian model with 173.98 m range. 
v <- vgm(1, "Gau", 173.98, nug = 0.0001)
# Note: we add a small nugget effect here for mathematical reasons (the variance-covariance matrix would not be invertible due to covariances at small distances being too similar).
g.dummy <- gstat(formula = z~1, dummy = TRUE, beta = 0,model = v, nmax = 20)
g.prdGau<-predict(g.dummy,Agro_mask,nsim=9)
spplot(g.prdGau)

# VARIOGRAM MODELING
vario.Sand <- variogram(Ca~1, AgroSP, cutoff=350, width=50)
plot(vario.Sand)
v <- vgm(500, "Sph", 200, nug = 250)
model.Sand = fit.variogram(vario.Sand, model = v)
plot(vario.Sand,model=model.Sand)
# To give an idea of the 'goodness of fit': 
attr(model.Sand,'SSErr') # Sums of squared differences between the variogram model and the experimental variogram. 

# Calcium
vario.Ca <- variogram(Ca ~ 1, AgroSP,cutoff=350, width=50)
plot(vario.Ca)
v <- vgm(1500000, "Sph", 200, nug = 250)
model.Ca = fit.variogram(vario.Ca, model = v)
plot(vario.Ca, model=model.Ca)


##===========================================================
## PART 4 - INTERPOLATION - KRIGING - CONDITIONAL SIMULATIONS
##===========================================================

# Inverse Distance Weighting
idwr <- idw(pH ~ 1, AgroSP, Agro_mask,idp=2)
spplot(idwr['var1.pred'])

# Thin plate spline
tps <- Tps(xy, Agro$Ca)
ras<-raster(Agro_mask)
spline.Ca <- interpolate(ras, tps)
spline.Ca <- mask(spline.Ca, ras)
spplot(spline.Ca)

# Simple Kriging
vario.Ca <- variogram(Ca~1,AgroSP, cutoff=350, width=50)
v <- vgm(500, "Sph", 200, nug = 250)
model.Ca = fit.variogram(vario.Ca, model = v)
kri <- krige(Ca ~ 1, AgroSP, Agro_mask, model = model.Ca,beta=1)
spplot(kri['var1.pred'])

# Ordinary Kriging
kri <- krige(Ca ~ 1, AgroSP, Agro_mask, model = model.Ca,maxdist=178)
spplot(kri['var1.pred']) # Map of kriging estimates
spplot(kri['var1.var']) # Map of kriging variances

# Save the kriged estimates to a TIF file in your working directory
writeGDAL(kri['var1.pred'],'CaKri.tif')

# Leave-one-out Cross-validation
crossval.Ca <- krige(Ca ~ 1, AgroSP[-1,], AgroSP[1,], model = model.Ca, maxdist=200)
for (i in 2:nrow(AgroSP)){
k <- krige(Ca ~ 1, AgroSP[-i,], AgroSP[i,], model = model.Ca, maxdist=200)
crossval.Ca=rbind(crossval.Ca,k);
}
plot(AgroSP$Ca,crossval.Ca$var1.pred)
cor(AgroSP$Ca,crossval.Ca$var1.pred) #Correlaton coefficient between estimated and actual values


# CONDITIONAL GAUSSIAN SIMULATIONS
#Calcium
vario.Ca <- variogram(Ca ~ 1, AgroSP,cutoff=350, width=50)
plot(vario.Ca)
v <- vgm(500, "Sph", 200, nug = 250)
model.Ca = fit.variogram(vario.Ca, model = v)

Ca.sim <- krige(formula = Ca~1, ArboSP, Arbo_mask, model = model.Ca,
                nmax = 15, nsim = 9)
spplot(Ca.sim)

# pH
vario.pH <- variogram(pH ~ 1, ArboSP,cutoff=350, width=50)
plot(vario.pH)
v <- vgm(500, "Gau", 200, nug = 250)
model.pH = fit.variogram(vario.pH, model = v)

pH.sim <- krige(formula = pH~1, ArboSP, Arbo_mask, model = model.pH,
                nmax = 35, nsim = 9)
spplot(pH.sim)



