rm(list=ls())

install.packages("mice")

library(sp)
library(mice)
library(raster)
library(gdal)
library(albers)
library(sf)
library(terra)

mydir <- "~/Desktop/data"

sr_orig <- read.delim("plantSr.txt")
sr_final_lat_long <- read.csv("real.csv")

sr_orig <-sr_orig[,c("Site", "Latitude", "Longitude", "X87Sr86Sr")]
sr_final <- sr_final_lat_long[,c("Latitude", "Longitude")]

albers<- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

df <- st_as_sf(x = sr_orig,                         
               coords = c("Longitude", "Latitude"),
               crs="+init=epsg:4326")

sr_proj <- st_transform(df, crs = albers)


df <- st_as_sf(x = sr_final,                         
               coords = c("Longitude", "Latitude"),
               crs="+init=epsg:4326")

sr_final <- st_transform(df, crs = albers)



##NOTE: This chunk only needs to be run ONCE, change to eval=FALSE for subsequent runs to save time.

Ras_Path <- "raster/raster"
OutRas_Path <- "raster/rasterNA"

#{r eval=FALSE, message=FALSE, warning=FALSE}
files<-list.files(path=paste(mydir,Ras_Path,sep="/"), pattern="\\.tif$")
outfiles <- paste0("p_", files)

eusa<-as.vector(c(-10000000, -4000000, 2000000, 7000000))

r.age<-raster("~/Desktop/data/raster/rasterNA/p_basement_age_reproj.tif")
r.ai<-raster("~/Desktop/data/raster/rasterNA/p_ai_reproj.tif")
r.bulk<-raster("~/Desktop/data/raster/rasterNA/p_rbulk_reproj.tif")
r.cec<-raster("~/Desktop/data/raster/rasterNA/p_rcec_reproj.tif")
r.clay<-raster("~/Desktop/data/raster/rasterNA/p_rclay_reproj.tif")
r.dust<-raster("~/Desktop/data/raster/rasterNA/p_dust_reproj.tif")
r.elevation<-raster("~/Desktop/data/raster/rasterNA/p_elevation_reproj.tif")
r.GUM<-raster("~/Desktop/data/raster/rasterNA/p_gum_mask3.tif")
r.lc<-raster("~/Desktop/data/raster/rasterNA/p_lc_reproj.tif")
r.litho<-raster("~/Desktop/data/raster/rasterNA/p_litho.tif")
r.m1<-raster("~/Desktop/data/raster/rasterNA/p_rm1_reproj.tif")
r.map<-raster("~/Desktop/data/raster/rasterNA/p_map_reproj.tif")
r.maxage_geol<-raster("~/Desktop/data/raster/rasterNA/p_agemax.tif")
r.meanage_geol<-raster("~/Desktop/data/raster/rasterNA/p_agemean.tif")
r.minage_geol<-raster("~/Desktop/data/raster/rasterNA/p_agemin.tif")
r.nfert<-raster("~/Desktop/data/raster/rasterNA/p_nfert_reproj.tif")
r.orc<-raster("~/Desktop/data/raster/rasterNA/p_rorc_reproj.tif")
r.pet<-raster("~/Desktop/data/raster/rasterNA/p_pet_reproj.tif")
r.ph<-raster("~/Desktop/data/raster/rasterNA/p_rph_reproj.tif")
r.phkcl<-raster("~/Desktop/data/raster/rasterNA/p_phkcl_reproj.tif")
r.salt<-raster("~/Desktop/data/raster/rasterNA/p_salt_reproj.tif")
r.sand<-raster("~/Desktop/data/raster/rasterNA/p_sand_reproj.tif")
r.silt<-raster("~/Desktop/data/raster/rasterNA/p_silt_reproj.tif")
r.soilthick<-raster("~/Desktop/data/raster/rasterNA/p_soilthick_reproj.tif")
r.srsrq1<-raster("~/Desktop/data/raster/rasterNA/p_srsrq1.tif")
r.srsrq3<-raster("~/Desktop/data/raster/rasterNA/p_srsrq3.tif")
r.ssa<-raster("~/Desktop/data/raster/rasterNA/p_ssa.tif")
r.ssaw<-raster("~/Desktop/data/raster/rasterNA/p_ssaw.tif")
r.xx<-raster("~/Desktop/data/raster/rasterNA/p_xx.tif")


agexy<-raster::extract(r.age, sr_proj, method='bilinear', na.rm=TRUE) 
aixy<-raster::extract(r.ai, sr_proj, mmethod='bilinear', na.rm=TRUE)
bulkxy<-raster::extract(r.bulk,sr_proj, method='bilinear', na.rm=TRUE)
cecxy<-raster::extract(r.cec,sr_proj, method='bilinear', na.rm=TRUE) 
clayxy<-raster::extract(r.clay,sr_proj, method='bilinear', na.rm=TRUE)
dustxy<-raster::extract(r.dust, sr_proj, method='bilinear', na.rm=TRUE)
elevationxy<-raster::extract(r.elevation, sr_proj, method='bilinear', na.rm=TRUE)
GUMxy<-raster::extract(r.GUM,sr_proj, method='simple', na.rm=TRUE) #b/c categorical
lcxy=raster::extract(r.lc,sr_proj, method='bilinear', na.rm=TRUE)
lithoxy<-raster::extract(r.litho,sr_proj, method='bilinear', na.rm=TRUE)
m1xy<-raster::extract(r.m1, sr_proj, method='bilinear', na.rm=TRUE) 
mapxy<-raster::extract(r.map, sr_proj, method='bilinear', na.rm=TRUE)
maxage_geolxy<-raster::extract(r.maxage_geol,sr_proj, method='bilinear', na.rm=TRUE)
meanage_geolxy<-raster::extract(r.meanage_geol,sr_proj, method='bilinear', na.rm=TRUE)
minage_geolxy<-raster::extract(r.minage_geol,sr_proj, method='bilinear', na.rm=TRUE)
nfertxy<-raster::extract(r.nfert,sr_proj, method='bilinear', na.rm=TRUE)
orcxy<-raster::extract(r.orc,sr_proj, method='bilinear', na.rm=TRUE)
petxy<-raster::extract(r.pet, sr_proj, method='bilinear', na.rm=TRUE)
phxy<-raster::extract(r.ph,sr_proj, method='bilinear', na.rm=TRUE)
phkclxy<-raster::extract(r.phkcl,sr_proj, method='bilinear', na.rm=TRUE)
saltxy<-raster::extract(r.salt, sr_proj, method='bilinear', na.rm=TRUE)
sandxy<-raster::extract(r.sand,sr_proj, method='bilinear', na.rm=TRUE)
siltxy<-raster::extract(r.silt,sr_proj, method='bilinear', na.rm=TRUE)
soilthickxy<-raster::extract(r.soilthick,sr_proj, method='bilinear', na.rm=TRUE)
srsrq1xy<-raster::extract(r.srsrq1,sr_proj, method='bilinear', na.rm=TRUE)
srsrq3xy<-raster::extract(r.srsrq3,sr_proj, method='bilinear', na.rm=TRUE)
ssaxy<-raster::extract(r.ssa,sr_proj, method='bilinear', na.rm=TRUE)
ssawxy<-raster::extract(r.ssaw,sr_proj, method='bilinear', na.rm=TRUE)
xxxy<-raster::extract(r.xx,sr_proj, method='bilinear', na.rm=TRUE)
sr_proj_xy <- data.frame(sr_proj, sr_orig[c("Latitude","Longitude")], agexy, aixy, bulkxy,
                         cecxy,clayxy,dustxy,elevationxy,GUMxy,
                         lcxy,lithoxy,m1xy,mapxy,maxage_geolxy,
                         meanage_geolxy,minage_geolxy, nfertxy,
                         orcxy,petxy,phxy,phkclxy,saltxy,sandxy,
                         siltxy,soilthickxy,srsrq1xy,srsrq3xy,
                         ssaxy, ssawxy,xxxy)
colnames(sr_proj_xy)<-c("Site","X87Sr86Sr","geometry","Latitude","Longitude",
                        "age","ai","bulk","cec",
                        "clay","dust","elevation","GUM",
                        "lc","litho","m1","map","maxage_geol",
                        "meanage_geol","minage_geol","nfert",
                        "orc","pet","ph","phkcl","salt",
                        "sand","silt","soilthick","srsrq1",
                        "srsrq3","ssa","ssaw","xx") 
sr_proj_xy<-subset(sr_proj_xy,select=c("Site","X87Sr86Sr","Latitude","Longitude",
                                       "age","ai","bulk","cec",
                                       "clay","dust","elevation","GUM",
                                       "lc","litho","m1","map","maxage_geol",
                                       "meanage_geol","minage_geol","nfert",
                                       "orc","pet","ph","phkcl","salt",
                                       "sand","silt","soilthick","srsrq1",
                                       "srsrq3","ssa","ssaw","xx")) 

imputed <- mice(sr_proj_xy, m = 1, method = "rf", remove.collinear=FALSE)
completed <- complete(imputed)

write.csv(completed, "~/Desktop/data/final_data.csv", row.names=FALSE)


### New Points

agexy<-raster::extract(r.age, sr_final, method='bilinear', na.rm=TRUE) 
aixy<-raster::extract(r.ai, sr_final, mmethod='bilinear', na.rm=TRUE)
bulkxy<-raster::extract(r.bulk,sr_final, method='bilinear', na.rm=TRUE)
cecxy<-raster::extract(r.cec,sr_final, method='bilinear', na.rm=TRUE) 
clayxy<-raster::extract(r.clay,sr_final, method='bilinear', na.rm=TRUE)
dustxy<-raster::extract(r.dust, sr_final, method='bilinear', na.rm=TRUE)
elevationxy<-raster::extract(r.elevation, sr_final, method='bilinear', na.rm=TRUE)
GUMxy<-raster::extract(r.GUM,sr_final, method='simple', na.rm=TRUE) #b/c categorical
lcxy=raster::extract(r.lc,sr_final, method='bilinear', na.rm=TRUE)
lithoxy<-raster::extract(r.litho,sr_final, method='bilinear', na.rm=TRUE)
m1xy<-raster::extract(r.m1, sr_final, method='bilinear', na.rm=TRUE) 
mapxy<-raster::extract(r.map, sr_final, method='bilinear', na.rm=TRUE)
maxage_geolxy<-raster::extract(r.maxage_geol,sr_final, method='bilinear', na.rm=TRUE)
meanage_geolxy<-raster::extract(r.meanage_geol,sr_final, method='bilinear', na.rm=TRUE)
minage_geolxy<-raster::extract(r.minage_geol,sr_final, method='bilinear', na.rm=TRUE)
nfertxy<-raster::extract(r.nfert,sr_final, method='bilinear', na.rm=TRUE)
orcxy<-raster::extract(r.orc,sr_final, method='bilinear', na.rm=TRUE)
petxy<-raster::extract(r.pet, sr_final, method='bilinear', na.rm=TRUE)
phxy<-raster::extract(r.ph,sr_final, method='bilinear', na.rm=TRUE)
phkclxy<-raster::extract(r.phkcl,sr_final, method='bilinear', na.rm=TRUE)
saltxy<-raster::extract(r.salt, sr_final, method='bilinear', na.rm=TRUE)
sandxy<-raster::extract(r.sand,sr_final, method='bilinear', na.rm=TRUE)
siltxy<-raster::extract(r.silt,sr_final, method='bilinear', na.rm=TRUE)
soilthickxy<-raster::extract(r.soilthick,sr_final, method='bilinear', na.rm=TRUE)
srsrq1xy<-raster::extract(r.srsrq1,sr_final, method='bilinear', na.rm=TRUE)
srsrq3xy<-raster::extract(r.srsrq3,sr_final, method='bilinear', na.rm=TRUE)
ssaxy<-raster::extract(r.ssa,sr_final, method='bilinear', na.rm=TRUE)
ssawxy<-raster::extract(r.ssaw,sr_final, method='bilinear', na.rm=TRUE)
xxxy<-raster::extract(r.xx,sr_final, method='bilinear', na.rm=TRUE)
sr_final_xy <- data.frame(sr_final, sr_final_lat_long[c("Latitude","Longitude")], agexy, aixy, bulkxy,
                         cecxy,clayxy,dustxy,elevationxy,GUMxy,
                         lcxy,lithoxy,m1xy,mapxy,maxage_geolxy,
                         meanage_geolxy,minage_geolxy, nfertxy,
                         orcxy,petxy,phxy,phkclxy,saltxy,sandxy,
                         siltxy,soilthickxy,srsrq1xy,srsrq3xy,
                         ssaxy, ssawxy,xxxy)
colnames(sr_final_xy)<-c("geometry","Latitude","Longitude",
                        "age","ai","bulk","cec",
                        "clay","dust","elevation","GUM",
                        "lc","litho","m1","map","maxage_geol",
                        "meanage_geol","minage_geol","nfert",
                        "orc","pet","ph","phkcl","salt",
                        "sand","silt","soilthick","srsrq1",
                        "srsrq3","ssa","ssaw","xx") 
sr_final_xy<-subset(sr_final_xy,select=c("Latitude","Longitude",
                                       "age","ai","bulk","cec",
                                       "clay","dust","elevation","GUM",
                                       "lc","litho","m1","map","maxage_geol",
                                       "meanage_geol","minage_geol","nfert",
                                       "orc","pet","ph","phkcl","salt",
                                       "sand","silt","soilthick","srsrq1",
                                       "srsrq3","ssa","ssaw","xx")) 

nas <- data.frame(rowSums(is.na(sr_final_xy)))
colnames(nas) <- "nas"

sr_final_xy$nas <- nas$nas

high_na_df <- sr_final_xy %>% filter(nas >= 10)

final_no_na <- sr_final_xy %>% filter(nas < 1)

imputed <- mice(final_no_na, m = 1, method = "rf", remove.collinear=FALSE)
completed <- complete(imputed)

write.csv(completed, "~/Desktop/data/final_data_new_points.csv", row.names=FALSE)
