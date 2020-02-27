## BEFORE YOU BEGIN
install.packages(c('gstat','lattice','rgdal','spdep','vegan','SpatialPack','ggplot2','raster','fields'))

# Note for Linux users. If the installation of rdgal does not work, enter this command in the terminal:
# > sudo apt-get install libgdal-dev gdal-bin proj-bin libproj-dev

## PART I - LOAD NECESSARY PACKAGES

#Load each package
library(gstat)
library(lattice)
library(rgdal)
library(spdep)
library(vegan)
library(SpatialPack)
library(raster)
library(fields)


# Select the file 'Arboretum_Data.csv' on your computer. 
Arbo<-read.csv(file.choose())

# Variables names/acronyms
# Coordinates("X","Y"), "ID", Total carbon 0-15 cm of soil("tC_015"),Total carbon 15-30 cm of soil ("tC_1530"), Total carbon in forest floow ("FFtC"), Total soil carbon("tC"), Sand %("Sand"), Silt %("Silt"), Clay %("Clay"), Soil pH("pH"), Potassium ("K"), Phosphorus("P"), Magnesium ("Mg"), Calcium ("Ca"), N mineralization("Min"), Basal areas for: Red maple("RM"), Sugar maple("SuM"), Silver maple("SiM"), American beech("BE"), Red oak("RO"), White birch("WB"), Gray birch("GB"), Yellow birch("YB"), Basswood("BW"), Bitternut hickory ("BH"), Shagbark Hickory ("SH"), American elm("AE"), White ash("WA"), Black ash("BA"), Ironwood("IW"), Musclewood("MW"), Hemlock("HE") 

# Select the Arbo_mask.tif file on your computer
Arbo_mask<-readGDAL(file.choose())
xy <- Arbo[1:2]
df <- Arbo[-1:-2]
ArboSP <- SpatialPointsDataFrame(coords=xy, data=df)
proj4string(ArboSP)<-proj4string(Arbo_mask)

# Specify the theme to be used in package SP
sp.theme(set = TRUE, regions = list(col = terrain.colors(100)))
sp.theme(set = TRUE, regions = list(col = colorRampPalette(c("black","brown","orange","yellow","white"))(50)))

# Visualize some variables 
spplot(ArboSP['Ca'])
spplot(ArboSP['tC'])
spplot(ArboSP['pH'])


##============================================================================
## PART 2 - EXPERIMENTAL VARIOGRAM AND AUTOCORRELATION MEASURES
##============================================================================

# Semi-variogram for % Sand. 
# Variogram cloud
vario.SandCloud <- variogram(Sand~1,ArboSP,cloud=TRUE, cutoff=350)
plot(vario.SandCloud)

# Average at each lag ('classical' variogram)
vario.Sand <- variogram(Sand~1,ArboSP, cutoff=350, width=50)
plot(vario.Sand)

# Variogram for total carbon
vario.tC <- variogram(tC ~ 1, ArboSP, cutoff=350, width=50)
plot(vario.tC)

# Choose the file GearyC.R
source(file.choose())
# Choose the file MoranI.R
source(file.choose())

# C de Geary
corgramC<-GearyC(cbind(ArboSP$X,ArboSP$Y),ArboSP$Sand,seq(0,350,50))
ggplot(corgramC, aes(x=distances, y=C)) + geom_errorbar(aes(ymin=C-(1.96*sqrt(varC)), ymax=C+(1.96*sqrt(varC))), width=5)+geom_point()+geom_hline(aes(yintercept=1))+xlim(0,350)+ylim(0,1.2)+ggtitle("Geary's C")

# I de Moran
corgramI<-MoranI(cbind(ArboSP$X,ArboSP$Y),ArboSP$Sand,seq(0,350,50))
ggplot(corgramI, aes(x=distances, y=I)) + geom_errorbar(aes(ymin=I-(1.96*sqrt(varI)), ymax=I+(1.96*sqrt(varI))), width=5)+geom_point()+geom_hline(aes(yintercept=0))+xlim(0,350)+ylim(-0.2,1)+ggtitle("Moran's I")


##====================================================
## PART 3 - UNCONDITIONAL SIMULATIONS
##====================================================

# Unconditional simulations
# Pure nugget (no spatial structure)
v <- vgm(1,'Nug',0)
g.dummy <- gstat(formula = z~1, locations=ArboSP, dummy = TRUE, beta = 0,model = v, nmax = 20)
g.prdNug<-predict(g.dummy,Arbo_mask,nsim=9)
spplot(g.prdNug)

# Spherical model with 200 m range. 
v <- vgm(1, "Sph", 200)
g.dummy <- gstat(formula = z~1, locations=ArboSP, dummy = TRUE, beta = 0,model = v, nmax = 20)
g.prdSph<-predict(g.dummy,Arbo_mask,nsim=9)
spplot(g.prdSph)

# Gaussian model with 200 m range. 
v <- vgm(1, "Gau", 200, nug = 0.0001)
# Note: we add a small nugget effect here for mathematical reasons (the variance-covariance matrix would not be invertible due to covariances at small distances being too similar).
g.dummy <- gstat(formula = z~1, dummy = TRUE, beta = 0,model = v, nmax = 20)
g.prdGau<-predict(g.dummy,Arbo_mask,nsim=9)
spplot(g.prdGau)

# VARIOGRAM MODELING
vario.Sand <- variogram(Sand~1, ArboSP, cutoff=350, width=50)
plot(vario.Sand)
v <- vgm(500, "Sph", 200, nug = 250)
model.Sand = fit.variogram(vario.Sand, model = v)
plot(vario.Sand,model=model.Sand)
# To give an idea of the 'goodness of fit': 
attr(model.Sand,'SSErr') # Sums of squared differences between the variogram model and the experimental variogram. 

# Variogram for Sugar maple
vario.SuM <- variogram(SuM ~ 1, ArboSP, cutoff=350, width=50)
plot(vario.SuM)
v <- vgm(500, "Sph", 200, nug = 250) # 'Initial' values can be guessed by visually inspecting the experimental variogram
model.SuM = fit.variogram(vario.SuM, model = v) # Loot at model.SuM to see estimated values. 
plot(vario.SuM, model=model.SuM)

# Total carbon
vario.tC <- variogram(tC ~ 1, ArboSP, cutoff=350, width=50)
plot(vario.tC)
v <- vgm(1, "Nug", 0)
# Pure nugget effect with no spatial component
model.tC = fit.variogram(vario.tC, model = v)
plot(vario.tC, model=model.tC)

# Calcium
vario.Ca <- variogram(Ca ~ 1, ArboSP,cutoff=350, width=50)
plot(vario.Ca)
v <- vgm(1500000, "Sph", 200, nug = 250)
model.Ca = fit.variogram(vario.Ca, model = v)
plot(vario.Ca, model=model.Ca)

##===========================================================
## PART 4 - INTERPOLATION - KRIGING - CONDITIONAL SIMULATIONS
##===========================================================

# Inverse Distance Weighting
idwr <- idw(pH ~ 1, ArboSP, Arbo_mask,idp=2)
spplot(idwr['var1.pred'])

# Thin plate spline
tps <- Tps(xy, Arbo$pH)
ras<-raster(Arbo_mask)
spline.pH <- interpolate(ras, tps)
spline.pH <- mask(spline.pH, ras)
spplot(spline.pH)

# Simple Kriging
vario.pH <- variogram(pH~1,ArboSP, cutoff=350, width=50)
v <- vgm(500, "Sph", 200, nug = 250)
model.pH = fit.variogram(vario.pH, model = v)
kri <- krige(pH ~ 1, ArboSP, Arbo_mask, model = model.pH,beta=1)
spplot(kri['var1.pred'])

# Ordinary Kriging
kri <- krige(pH ~ 1, ArboSP, Arbo_mask, model = model.pH,maxdist=200)
spplot(kri['var1.pred']) # Map of kriging estimates
spplot(kri['var1.var']) # Map of kriging variances

# Save the kriged estimates to a TIF file in your working directory
writeGDAL(kri['var1.pred'],'phKri.tif')

# Leave-one-out Cross-validation
crossval.pH <- krige(pH ~ 1, ArboSP[-1,], ArboSP[1,], model = model.pH, maxdist=200)
for (i in 2:nrow(ArboSP)){
k <- krige(pH ~ 1, ArboSP[-i,], ArboSP[i,], model = model.pH, maxdist=200)
crossval.pH=rbind(crossval.pH,k);
}
plot(ArboSP$pH,crossval.pH$var1.pred)
cor(ArboSP$pH,crossval.pH$var1.pred) #Correlaton coefficient between estimated and actual values


# CONDITIONAL GAUSSIAN SIMULATIONS
#Calcium
vario.Ca <- variogram(Ca ~ 1, ArboSP,cutoff=350, width=50)
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

##====================================================
## PART 5 - BIVARIATE ANALYSIS
##====================================================

rm(Arbo.g)
# Note: we use the scale() fonction here to standardize the data, which is common practice in multivariate analysis. 
Arbo.g <- gstat(id="Wa", formula=scale(WA)~1, data=ArboSP, nmax = 10)
Arbo.g <- gstat(Arbo.g, "SuM", scale(SuM)~1, ArboSP, nmax = 10)
Xmodel <- vgm(0.3, "Sph", 100, nug = 0.1)
Xmodel <- vgm(0.7, "Sph", 250,add.to=Xmodel)
Xvario <- variogram(Arbo.g, cutoff=300,40)
XArbo.fit = fit.lmc(Xvario, Arbo.g,Xmodel,fit.lmc=TRUE)

# Cross-variogram
plot(Xvario, XArbo.fit)
xv=Xvario$id=='Wa.SuM'
v1=Xvario$id=='Wa'
v2=Xvario$id=='SuM'

# Calculation of hulls of perfect correlations
hulls=sqrt(Xvario[v1,3]*Xvario[v2,3])
hulls=cbind(hulls,-hulls)

plot(Xvario[xv,2],Xvario[xv,3],ylim=range(hulls))
par(new=FALSE)
lines(Xvario[v1,2],hulls[,1],col='red')
lines(Xvario[v1,2],hulls[,2],col='red')

# Codispersion coefficient (values between 0 et 1 at each lag)
codcoef<-Xvario[xv,3]/sqrt(Xvario[v1,3]*Xvario[v2,3])
plot(Xvario[v1,2],codcoef)

# Correlation with Dutilleul's modified T-test
# Sugar maple and pH
corr1<-modified.ttest(ArboSP$SuM,ArboSP$pH, coordinates(ArboSP))
print(corr1) # Note the really low effective sample size (adjusted degrees of freedom), due to the strong spatial correlation in the variables

# Sugar maple and total soil carbon
corr2<-modified.ttest(ArboSP$SuM,ArboSP$tC, coordinates(ArboSP))
print(corr2) # Weak spatial structure for Carbon = effective sample size (adjusted degrees of freedom) near the actual sample size value. 