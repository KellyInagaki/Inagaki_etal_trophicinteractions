
#.....................................................................................#
# Trophic interactions will expand geographically, but be less intense as oceans warm
# Inagaki, K.Y., Pennino, M.G., Floeter, S.R., Hay, M.E., Longo, G.O.
#.....................................................................................#

#.....................................................................................#
##### Load required packages #####
library(sp)
library(INLA)
library(geoR)
library(dismo)
library(rgeos)
library(rgdal)
library(stats)
library(spdep)
library(Hmisc)
library(raster)
library(GGally)
library(fields)
library(ggplot2)
library(maptools)
library(corrplot)
library(sdmpredictors)

#.....................................................................................#
##### Load dataframe with biomass and scaled environmental data #####
load("BiomassData.RData")
summary(sdmdata)

#.....................................................................................#
##### Exploration of the dataset #####
#.......Step 1: Check correlation among explicative variables ........#
matrix<-rcorr(as.matrix(sdmdata[,c(2:9)]), type = "spearman")

# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


corrplot(matrix$r, type="lower", tl.col = "black",method="number",p.mat = matrix$P, sig.level = 0.05)# in Corrplot "X" are no significant variables. We look at correlation among variables

#...... Step 2: Check multicollinearity among variables .......#
source("HighstatLib.r")
corvif(sdmdata[,2:9])

#.....................................................................................#
##### ESTIMATION #####
#...... Built the mesh .......#
sea=readOGR(".","Mar_Simplificado")
bound<-inla.sp2segment(sea)
mesh<-inla.mesh.2d(loc=as.matrix(sdmdata[,10:11]),
                   boundary=bound,max.edge=c(4,9),offset=c(5,3), cutoff=2)


##### Definition of the spde #####
# Total diameter
size = min(c(diff(range(sdmdata$lon)), diff(range(sdmdata$lat))))
# Mean range
range0 = size / 2

spde <- inla.spde2.pcmatern(mesh,
                            prior.sigma = c(1, 0.3), prior.range = c(range0, 0.3))

# Matrix which link data with the mesh
A.est <- inla.spde.make.A(mesh, loc=cbind(sdmdata$lon, sdmdata$lat))

#  inla.stack to stimate
stk.est<-inla.stack(data=list(y=sdmdata$biomass),
                    A=list(A.est, 1),
                    effects=list(spatial=1:spde$n.spde,
                                 data.frame(beta0=1, sdmdata[,2:9])),
                    tag='est')
head(sdmdata)

##### Fitting the model #####
formula.1 <- y~-1 + beta0 +  SST + pH+ f(spatial,model=spde)
model.est <- inla(formula.1, 
                  data=inla.stack.data(stk.est), family="gaussian" ,
                  control.compute=list(cpo=TRUE, waic=TRUE, return.marginals=TRUE), 
                  control.predictor=list(A=inla.stack.A(stk.est), compute=TRUE, 
                                         quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975)),
                  control.inla=list(strategy = "laplace"),
                  num.threads = 3,
                  verbose=T)

summary(model.est) 

##### How probable are these variables? #####
1-inla.pmarginal(0, model.est$marginals.fixed$pH)
1-inla.pmarginal(0, model.est$marginals.fixed$SST)

##### Plot posteriors #####
par(mfrow=c(2,2))
plot(model.est$marginals.fixed$beta0, type="l",main="Posterior distribution of Intercept")
abline(v=0, col="red", lwd=2)
plot(model.est$marginals.fixed$pH, type="l",main="Posterior distribution of pH")
abline(v=0, col="red", lwd=2)
plot(model.est$marginals.fixed$SST, type="l",main="Posterior distribution of SST")
abline(v=0, col="red", lwd=2)

##### Posterior distribution hyperpars #####
# Check the range is smaller than the offset
spde.result = inla.spde2.result(model.est, "spatial", spde, do.transform=TRUE)

range<-inla.emarginal(function(x) x, spde.result$marginals.range.nominal[[1]]) 

# Check the range is smaller than the offset
range < max(diff(range(sdmdata[,10])), diff(range(sdmdata[,11])))*0.40 #Yes!!!

# Plot 
par(mfrow=c(2,2), mar=c(3,3,1,0.5)) 
plot(spde.result$marginals.range.nominal[[1]], type='l') 

##### Interpolate the posterior mean and sd #####
# plot in a grid m X m 
# Customize the grid to predict
cat_rec=readOGR(".","zona")
coast <- gDifference(cat_rec,sea)
sea2=gDifference(cat_rec,coast)

bbox(sea2)
(dxy <- apply(bbox(sea2),1, diff))
(r <- dxy[1]/dxy[2])
m<-150
proj.grid.mat <- 
  inla.mesh.projector(mesh, 
                      xlim=bbox(sea2)[1,],
                      ylim=bbox(sea2)[2,] ,
                      dims=c(r, 1)*m)

plot(sea2)
points(proj.grid.mat$lattice$loc, pch=20, cex=0.5)

# clean (set NA to the values outside boundary)
ov <- over(SpatialPoints(proj.grid.mat$lattice$loc, sea2@proj4string), sea2)

# check grid points inside the map
i.map <- is.na(ov)

# Plot the points where we will predict
par(mar=c(0,0,0,0))
plot(cat_rec)
points(proj.grid.mat$lattice$loc[!i.map,], col="red", cex=0.2)
points(proj.grid.mat$lattice$loc[i.map,], col="blue", cex=0.2)

# consider only those inside map
proj.grid.mat$lattice$loc[i.map, ]

# Project the values of the mean and sd of the spatial effect 
mean.g <- inla.mesh.project(proj.grid.mat, model.est$summary.random$spatial$mean)
sd.g <- inla.mesh.project(proj.grid.mat, model.est$summary.random$spatial$sd)
quantile_0.025 <- inla.mesh.project(proj.grid.mat, model.est$summary.random$spatial$`0.025quant`)
quantile_0.975 <- inla.mesh.project(proj.grid.mat, model.est$summary.random$spatial$`0.975quant`)
sd.g[i.map] <- mean.g[i.map] <- quantile_0.025[i.map] <- quantile_0.975[i.map] <- NA

# Plot mean spatial effect
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           mean.g, axes=TRUE, xlab="Longitude", ylab="Latitude")
plot(wrld_simpl,add=T, axes=TRUE,col='dark grey')

#####Plot SD of spatial effect #####
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           sd.g, axes=TRUE, xlab="Longitude", ylab="Latitude")
plot(wrld_simpl,add=T, axes=TRUE,col='dark grey')

##### Plot quantile_0.025 of spatial effect #####
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           quantile_0.025 , axes=TRUE, xlab="Longitude", ylab="Latitude")
plot(wrld_simpl,add=T, axes=TRUE,col='dark grey')

##### Plot quantile_0.975 of spatial effect #####
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           quantile_0.975, axes=TRUE, xlab="Longitude", ylab="Latitude")
plot(wrld_simpl,add=T, axes=TRUE,col='dark grey')

##### PREDICTION #####
##### Change dataframe for future scenarios
#####Load environmental predictors for predictions
load("envi_pred.RData")

# Matrix which link the mesh with coordinates to predict
A.pred <- inla.spde.make.A(mesh, loc=proj.grid.mat$lattice$loc[!i.map, ])

# Stack to predict
stk.pred <- inla.stack(data=list(y=NA),
                       A=list(A.pred, 1), 
                       effects=list(spatial=1:spde$n.spde,
                                    data.frame(beta0 = 1, 
                                              envi_pred)),
                       tag='pred')

stk <- inla.stack(stk.est, stk.pred)

# Run prediction model
model.pred <- inla(formula.1, 
                   data=inla.stack.data(stk), family="gaussian",
                   control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=1), #link:link is a vector of
                   #length given by the size of the response variable with values 1 if the corresponding
                   #data is missing and NA otherwise
                   control.inla=list(strategy = "laplace"), # Strategy
                   control.mode=list(theta=model.est$mode$theta, restart=TRUE), #Mode 
                   control.results=list(return.marginals.random=FALSE,
                                        return.marginals.predictor=FALSE), # Avoid some marginals
                   num.threads = 3,
                   verbose=T)


###### Plot predictions #####
# index for the prediction data
idx <- inla.stack.index(stk, 'pred')$data

summary(model.pred$summary.fitted.val$mean[idx])

# Organize probabilities into a matrix to visualize
prob.mean <- prob.sd <- prob.0.025<- prob.0.975 <- matrix(NA, proj.grid.mat$lattice$dims[1],
                                                          proj.grid.mat$lattice$dims[2])
prob.mean[!i.map] <- c((model.pred$summary.fitted.val$mean[idx]))
prob.sd[!i.map] <- c(model.pred$summary.fitted.val$sd[idx])
prob.0.025[!i.map] <- c(model.pred$summary.fitted.val$`0.025quant`[idx])
prob.0.975[!i.map] <- c(model.pred$summary.fitted.val$`0.975quant`[idx])

#Plot biomass mean 
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
          prob.mean, axes=TRUE, xlab=("Longitude"), ylab="Latitude")
plot(wrld_simpl,add=T, axes=TRUE,col='dark grey')

#Plot biomass sd
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           prob.sd, axes=TRUE, xlab=("Longitude"), ylab="Latitude")
plot(wrld_simpl,add=T, axes=TRUE,col='dark grey')

#Plot biomass quantile_0.025 
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           prob.0.025, axes=TRUE, xlab=("Longitude"), ylab="Latitude")
plot(wrld_simpl,add=T, axes=TRUE,col='dark grey')

#Plot biomass quantile_0.975 
image.plot(proj.grid.mat$x, 
           proj.grid.mat$y,
           prob.0.975, axes=TRUE, xlab=("Longitude"), ylab="Latitude")
plot(wrld_simpl,add=T, axes=TRUE,col='dark grey')

