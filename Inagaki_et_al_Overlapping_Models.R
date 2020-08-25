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
library(tidyverse)
library(lubridate)
library(gcookbook)


#Importing data
oc_Acanthuridae=raster("aca_oc_atual_cut.asc")
bm_Acanthuridae=raster("aca_bm_atual_cut.asc")
fp_Acanthuridae=raster("aca_fp_atual_cut.asc")

# overlapping
data(wrld_simpl)
par(mfrow=c(1,4), mar=c(5,4,4,5))
summary(oc_Acanthuridae)
plot((oc_Acanthuridae/0.99999999),col=tim.colors(100)[1:100],main="atual Acanthuridae occurrence", axes=T)
plot(wrld_simpl,xlim=c(-88,-34), ylim=c(-28,35), add=T, axes=TRUE,col='dark grey')
box()

summary(bm_Acanthuridae)
plot((bm_Acanthuridae),col=tim.colors(100)[1:100],main="atual Acanthuridae Biomass", axes=T)
plot(wrld_simpl,xlim=c(-88,-34), ylim=c(-28,35), add=T, axes=TRUE,col='dark grey')
box()

summary(fp_Acanthuridae)
plot((fp_Acanthuridae),col=tim.colors(100)[1:100],main="atual Acanthuridae feeding pressure", axes=T)
plot(wrld_simpl,xlim=c(-88,-34), ylim=c(-28,35), add=T, axes=TRUE,col='dark grey')
box()


multiplicado_aca<-(oc_Acanthuridae/0.9999999)* (bm_Acanthuridae) * (fp_Acanthuridae)
summary(multiplicado_aca)
plot((multiplicado_aca/0.4952246),col=tim.colors(100)[1:100],main="Intensity of trophic interactions", axes=T)
plot(wrld_simpl,xlim=c(-88,-34), ylim=c(-28,35), add=T, axes=TRUE,col='dark grey')
box()


# Save .ascii
writeRaster(multiplicado_aca, filename="Acanthuridae_atual_all_projection_escaled.asc", 
            format="ascii", overwrite=TRUE)
multiplicado_aca <- raster("Acanthuridae_atual_all_projection_escaled.asc")



##### 2050 #####
# importing data
oc_Acanthuridae_205085max=raster("aca_oc_2050_cut.asc")
bm_Acanthuridae_205085max=raster("aca_bm_2050_cut.asc")
fp_Acanthuridae_205085max=raster("aca_fp_2050_cut.asc")

# overlapping
data(wrld_simpl)
par(mfrow=c(1,4), mar=c(5,4,4,5))

summary(oc_Acanthuridae_205085max)
plot((oc_Acanthuridae_205085max),col=tim.colors(100)[1:100],main="205085max Acanthuridae occurrence", axes=T)
plot(wrld_simpl,xlim=c(-88,-34), ylim=c(-28,35), add=T, axes=TRUE,col='dark grey')
box()

summary(bm_Acanthuridae_205085max)
plot((bm_Acanthuridae_205085max),col=tim.colors(100)[1:100],main="205085max Acanthuridae Biomass", axes=T)
plot(wrld_simpl,xlim=c(-88,-34), ylim=c(-28,35), add=T, axes=TRUE,col='dark grey')
box()

summary(fp_Acanthuridae_205085max)
plot((fp_Acanthuridae_205085max),col=tim.colors(100)[1:100],main="205085max Acanthuridae feeding pressure", axes=T)
plot(wrld_simpl,xlim=c(-88,-34), ylim=c(-28,35), add=T, axes=TRUE,col='dark grey')
box()


multiplicado_aca_205085max<-(oc_Acanthuridae_205085max) * (bm_Acanthuridae_205085max) * (fp_Acanthuridae_205085max)
summary(multiplicado_aca_205085max)
plot((multiplicado_aca_205085max/0.4908844),col=tim.colors(100)[1:100],main="Intensity of trophic interactions", axes=T)
plot(wrld_simpl,xlim=c(-88,-34), ylim=c(-28,35), add=T, axes=TRUE,col='dark grey')
box()


# Save .ascii
writeRaster(multiplicado_aca_205085max, filename="Acanthuridae_205085max_all_projection_escaled.asc", 
            format="ascii", overwrite=TRUE)
multiplicado_aca_205085max <- raster("Acanthuridae_205085max_all_projection_escaled.asc")


##### 2100 #####
# importing data
oc_Acanthuridae_210085max=raster("aca_oc_2100_cut.asc")
bm_Acanthuridae_210085max=raster("aca_bm_2100_cut.asc")
fp_Acanthuridae_210085max=raster("aca_fp_2100_cut.asc")

# overlapping
data(wrld_simpl)
par(mfrow=c(1,4), mar=c(5,4,4,5))

summary(oc_Acanthuridae_210085max)
plot((oc_Acanthuridae_210085max),col=tim.colors(100)[1:100],main="210085max Acanthuridae occurrence", axes=T)
plot(wrld_simpl,xlim=c(-88,-34), ylim=c(-28,35), add=T, axes=TRUE,col='dark grey')
box()

summary(bm_Acanthuridae_210085max)
plot((bm_Acanthuridae_210085max),col=tim.colors(100)[1:100],main="210085max Acanthuridae Biomass", axes=T)
plot(wrld_simpl,xlim=c(-88,-34), ylim=c(-28,35), add=T, axes=TRUE,col='dark grey')
box()

summary(fp_Acanthuridae_210085max)
plot((fp_Acanthuridae_210085max),col=tim.colors(100)[1:100],main="210085max Acanthuridae feeding pressure", axes=T)
plot(wrld_simpl,xlim=c(-88,-34), ylim=c(-28,35), add=T, axes=TRUE,col='dark grey')
box()


multiplicado_aca_210085max<-(oc_Acanthuridae_210085max) * (bm_Acanthuridae_210085max) * (fp_Acanthuridae_210085max)
summary(multiplicado_aca_210085max)
plot((multiplicado_aca_210085max/0.4991983),col=tim.colors(100)[1:100],main="Intensity of trophic interactions", axes=T)
plot(wrld_simpl,xlim=c(-88,-34), ylim=c(-28,35), add=T, axes=TRUE,col='dark grey')
box()

# Save .ascii
writeRaster(multiplicado_aca_210085max, filename="Acanthuridae_210085max_all_projection_escaled.asc", 
            format="ascii", overwrite=TRUE)
multiplicado_aca_210085max <- raster("Acanthuridae_210085max_all_projection_escaled.asc")


##### Extracting values #####
# atual
aca_atual <- raster("Acanthuridae_atual_all_projection_escaled.asc")
aca_atual = as.data.frame(cbind(coordinates(aca_atual),getValues(aca_atual)))
summary(aca_atual)
head(aca_atual)

to.remove <- which(!complete.cases(aca_atual))
aca_atual <- aca_atual[-to.remove,]
summary(aca_atual)
write.table(aca_atual,"aca_atual.txt")

# 2050
aca_2050 <- raster("Acanthuridae_205085max_all_projection_escaled.asc")
aca_2050 = as.data.frame(cbind(coordinates(aca_2050),getValues(aca_2050)))
summary(aca_2050)
head(aca_2050)

to.remove <- which(!complete.cases(aca_2050))
aca_2050 <- aca_2050[-to.remove,]
summary(aca_2050)

write.table(aca_2050,"aca_2050.txt")

# 2100
aca_2100 <- raster("Acanthuridae_210085max_all_projection_escaled.asc")
aca_2100 = as.data.frame(cbind(coordinates(aca_2100),getValues(aca_2100)))
summary(aca_2100)
head(aca_2100)

to.remove <- which(!complete.cases(aca_2100))
aca_2100 <- aca_2100[-to.remove,]
summary(aca_2100)

write.table(aca_2100,"aca_2100.txt")



##### Calculating differences #####
# 2050
acadif_2050 <- aca_2050$V3 - aca_atual$V3
acadif_2050 = cbind(aca_2050$x, aca_2050$y, acadif_2050)
acadif_2050 = as.data.frame(acadif_2050)
colnames(acadif_2050) = c("long", "lat", "dif")
acadif_2050 = aggregate(acadif_2050$dif ~ acadif_2050$lat, FUN = mean )
head(acadif_2050)
colnames(acadif_2050) = c("lat", "dif")
summary(acadif_2050)
pad_2050 <- acadif_2050$dif/0.0544787
acadif_2050 = cbind(acadif_2050$lat,pad_2050)
summary(acadif_2050)
acadif_2050 = as.data.frame(acadif_2050)
colnames(acadif_2050) = c("lat", "dif")
write.table(acadif_2050,"aca_dif_2050.txt")

# 2100
acadif_2100 <- aca_2100$V3 - aca_atual$V3
acadif_2100 = cbind(aca_2100$x, aca_2100$y, acadif_2100)
acadif_2100 = as.data.frame(acadif_2100)
colnames(acadif_2100) = c("long", "lat", "dif")
acadif_2100 = aggregate(acadif_2100$dif ~ acadif_2100$lat, FUN = mean )
head(acadif_2100)
colnames(acadif_2100) = c("lat", "dif")
summary(acadif_2100)
pad_2100 <- acadif_2100$dif/0.051946
acadif_2100 = cbind(acadif_2100$lat,pad_2100)
summary(acadif_2100)
acadif_2100 = as.data.frame(acadif_2100)
colnames(acadif_2100) = c("lat", "dif2100")
write.table(acadif_2100,"aca_dif_2100.txt")

##### PLOT #####
# 2050
aca_2050dif <- ggplot(acadif_2050, aes(x=lat, y=dif)) +
  geom_area(stat="identity") +
  geom_hline(yintercept=0, color="red") +
  coord_flip() +
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(),
        axis.line = element_blank())
aca_2050dif


# 2100
aca_2100dif <- ggplot(acadif_2100, aes(x=lat, y=dif2100)) +
  geom_area(stat="identity") +
  geom_hline(yintercept=0, color="red") +
  coord_flip() +
  theme(panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank(),
        axis.line = element_blank())
aca_2100dif

