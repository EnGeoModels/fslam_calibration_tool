#==============================
# R-FSLAM (calibration module)
#==============================

library(sf)
library(dplyr)
library(cvms)
library(caret)
library(nsga2R)
library(plotly)
library(psych)
library(pROC)

# Fixed parameters
g <- 9.81 # gravity (m/s^2)
dw <- 1000 # water density (kg/m^3)
b <- 5 # cell size (m)
ds <- 2000 # soil density (Kg/m^3)

#x <- c(Cs,Cr,phi_min,phi_max,Pa,z,Pe,k,por,CN)

# Input data for landslide (MORLE) and no-landslide (no-MORLE) conditions
points <- read.csv('.../MORLE_input_data/points_sample.csv')
soil <- read.csv('.../MORLE_input_data/soil.csv')
hmtu <- read.csv('.../MORLE_input_data/hmtu.csv')

lower_b <-  c(0.5,0.5,-5,-5,0,1,1,1,-0.3,1) # Parameter Lower boundaries 
upper_b <-  c(5,5,0,0,1,4.5,1,1,0,1) # Parameter Upper boundaries 

# FSLAM function
fslam <- function(x){
  # Editing parameters
  
  nrow_soil <- nrow(soil)
  nrow_hmtu <- nrow(hmtu)
  
  soil$h[2:nrow_soil] <- soil$h[2:nrow_soil] %>% as.numeric() + x[6]
  
  soil$Cmax[2:nrow_soil] <- (soil$Cmax[2:nrow_soil] %>% as.numeric)*x[1]
  hmtu$Cr_max <- c('(kPa)',(hmtu$Cr_max[2:(nrow_hmtu-1)] %>% as.numeric)*x[2],999)
  
  soil$phimin[2:nrow_soil] <- (soil$phimin[2:nrow(soil)] %>% as.numeric)+x[3]
  soil$phimax[2:nrow_soil] <-(soil$phimax[2:nrow(soil)] %>% as.numeric)+x[4]
  
  A <- (hmtu$A[2:nrow_hmtu] %>% as.numeric)*x[10]
  B <- (hmtu$B[2:nrow_hmtu] %>% as.numeric)*x[10]
  C <- (hmtu$C[2:nrow_hmtu] %>% as.numeric)*x[10]
  D <- (hmtu$D[2:nrow_hmtu] %>% as.numeric)*x[10]
  
  hmtu$A[2:nrow_hmtu] <- ifelse(A>100,100,A)
  hmtu$B[2:nrow_hmtu] <- ifelse(B>100,100,B)
  hmtu$C[2:nrow_hmtu] <- ifelse(C>100,100,C)
  hmtu$D[2:nrow_hmtu] <- ifelse(D>100,100,D)
  
  options(scipen = 999)
  soil$Ks[2:nrow_soil] <- (soil$Ks[2:nrow(soil)] %>% as.numeric)*x[8]
  
  soil$porosity[2:nrow_soil] <-  (soil$porosity[2:nrow_soil] %>% as.numeric)+x[9]
  
  # Defining mean and sd of stochastic variables (phi, Cr, Cs)
  
  df_phi <- data.frame(phi_min=soil$phimin[2:nrow(soil)] %>% as.numeric,
                       phi_max=soil$phimax[2:nrow(soil)] %>% as.numeric) %>% 
    mutate(phi_min_rad=phi_min*pi/180, phi_max_rad=phi_max*pi/180,
           tan_phi_min=tan(phi_min_rad),tan_phi_max=tan(phi_max_rad))
  
  
  tanphi <- (df_phi$tan_phi_max + df_phi$tan_phi_min)/2
  
  tanphisd <- (df_phi$tan_phi_max - df_phi$tan_phi_min)/4
  
  Cs <- apply(data.frame(Cs_min=(soil$Cmin[2:nrow(soil)] %>% as.numeric)*1000,
                         Cs_max=(soil$Cmax[2:nrow(soil)] %>% as.numeric)*1000),1,mean)
  
  Cssd <- ((soil$Cmax[2:nrow(soil)] %>% as.numeric)*1000 - (soil$Cmin[2:nrow(soil)] %>% as.numeric)*1000)/4
  
  Cr <-  apply(data.frame(Cr_min = (hmtu$Cr_min[2:nrow(hmtu)] %>% as.numeric)*1000,
                          Cr_max = (hmtu$Cr_max[2:nrow(hmtu)] %>% as.numeric)*1000),1,mean)
  
  Crsd <- ((hmtu$Cr_max[2:nrow(hmtu)] %>% as.numeric)*1000 - (hmtu$Cr_min[2:nrow(hmtu)] %>% as.numeric)*1000)/4
  
  # Assigning values:
  
  # Non-fixed parameters (soil and hmtu)
  
  points$Cs_m <- NA
  for(i in 1:length(Cs)){
    
    if(which(points$soil==i) %>% length == 0){
      next    
    } else {
      points[which(points$soil==i),]$Cs_m <- Cs[i]
    }
    
  } # effective cohesion
  
  points$Cs_sd <- NA
  for(i in 1:length(Cssd)){
    
    if(which(points$soil==i) %>% length == 0){
      next
    } else {
      points[which(points$soil==i),]$Cs_sd <- Cssd[i]
      
    }} 
  
  points$Cr_m <- NA
  for(i in 1:length(Cr)){
    
    if(which(points$lucl==i) %>% length == 0){
      next
    } else {
      points[which(points$lucl==i),]$Cr_m <- Cr[i]
    }
    
  } # apparent cohesion
  
  points$Cr_sd <- NA
  for(i in 1:length(Crsd)){
    
    if(which(points$lucl==i) %>% length == 0){
      next
    } else {
      points[which(points$lucl==i),]$Cr_sd <- Crsd[i]
      
    }}
  
  points$tanphi_m <- NA
  for(i in 1:length(tanphi)){
    
    if(which(points$soil==i) %>% length == 0){
      next
    } else {
      points[which(points$soil==i),]$tanphi_m <- tanphi[i]
      
    }} # tangent of internal friction angle
  
  points$tanphi_sd <- NA
  for(i in 1:length(tanphisd)){
    
    if(which(points$soil==i) %>% length == 0){
      next
    } else {
      points[which(points$soil==i),]$tanphi_sd <- tanphisd[i]
      
    }} 
  
  points$K <- NA
  for(i in 1:(nrow_soil-1)){
    
    if(which(points$soil==i) %>% length == 0){
      next
    } else {
      points[which(points$soil==i),]$K <- (soil$Ks[2:nrow_soil] %>% as.numeric)[i]
      
    }} # permeability
  
  points$z <- NA
  for(i in 1:(nrow_soil-1)){
    
    if(which(points$soil==i) %>% length == 0){
      next
    } else {
      points[which(points$soil==i),]$z <- (soil$h[2:nrow_soil] %>% as.numeric)[i]
      
    }} # soil depth
  
  points$n <- NA
  for(i in 1:(nrow_soil-1)){
    
    if(which(points$soil==i) %>% length == 0){
      next
    } else {
      points[which(points$soil==i),]$n <- (soil$porosity[2:nrow_soil] %>% as.numeric)[i]
      
    }} # porosity
  
  points$hsg <- NA
  for(i in 1:(nrow_soil-1)){
    
    if(which(points$soil==i) %>% length == 0){
      next
    } else {
      points[which(points$soil==i),]$hsg <- (soil$hsg[2:nrow_soil])[i]
      
    }} # hydrologic soil group
  
  
  # Updating Pa and Pe values
  points$Pa <- points$Pa*x[5] # effective infiltration - antecedent recharge (mm/day)
  points$Pe <- points$Pe*x[7] # Rainfall event (mm)
  
  # Preliminary updated table
  points_soil <- points
  
  #summary(points_soil)
  
  # curve number
  CN <- points_soil$CN
  
  # Estimating qe (recharge from event)
  Ia <- (5080/CN-50.8)
  Ia <- ifelse(Ia<0,yes = 0, no = Ia)
  
  qe <-  sapply(1:nrow(points_soil), function(x){
    #x=1223
    qe <- ifelse(points_soil$Pe[x] >= Ia[x], yes = points_soil$Pe[x]-
                   (points_soil$Pe[x]-Ia[x])^2/(points_soil$Pe[x]+4*Ia[x]),
                 no = points_soil$Pe[x])
    
  })
  
  points_final <- points_soil %>% mutate(C=Cs_m+Cr_m,CN,qe,K_mmd=K*3600*24*10^3,
                                         
                                         # Estimating water table
                                         ha = (area/b)*Pa/(K_mmd*z*sin(slope)*cos(slope)),
                                         ha_filter = ifelse(ha>1,1,ha)*z, # water table after Pa (m)
                                         h = ha_filter + qe/1000/n,
                                         h_filter = ifelse(h>z,z,h),# total water table in m (Pa+Pe)
                                         he = h_filter - ha_filter,
                            
                                         
                                         A = z*ds*g*sin(2*slope)/2,
                                         
                                         # SF post antecedent recharge
                                         D_ini = tan(slope)/(1-(ha_filter/z)*(dw/ds)),
                                         #D_ini = ds*tan(slope)/(ds-ha_filter/z*dw),
                                         SF_mean_ini = (tanphi_m/D_ini + C/A),
                                         #SF_sd_ini = sqrt(tanphi_sd^2/D_ini^2+(Cr_sd^2+Cs_sd^2)/A^2), # in paper
                                         SF_sd_ini = sqrt(A^2*tanphi_sd^2+D_ini^2*(Cs_sd^2+Cr_sd^2))/(D_ini*A), # in paper
                                         
                                         # SF post event
                                         D_fin = tan(slope)/(1-(h_filter/z)*(dw/ds)),
                                         #D_fin = ds*tan(slope)/(ds-h_filter/z*dw),
                                         SF_mean_fin = (tanphi_m/D_fin + C/A),
                                         #SF_sd_fin = sqrt(tanphi_sd^2/D_fin^2+(Cr_sd^2+Cs_sd^2)/A^2), # in paper
                                         SF_sd_fin = sqrt(A^2*tanphi_sd^2+D_fin^2*(Cs_sd^2+Cr_sd^2))/(D_fin*A), # in paper
                                         
                                         # Landslide prediction (Safe Factor)
                                         SF_ini_pred = ifelse(SF_mean_ini > 1,0,1),
                                         SF_fin_pred = ifelse(SF_mean_fin > 1,0,1),
                                         
                                         # Landslide prediction (Probability of Failure)
                                         POF_ini = pnorm(1, SF_mean_ini, SF_sd_ini),
                                         POF_fin = pnorm(1, SF_mean_fin, SF_sd_fin),
                                         POF_ini_pred = ifelse(POF_ini > 0.5,1,0),
                                         POF_fin_pred = ifelse(POF_fin > 0.5,1,0)
  )
  
  
  #hist(points_final$h_filter/points_final$z)
  #summary(points_final)
  points_final <- points_final %>% mutate(ratio_ha_z = (ha_filter/z) %>% round(2),
                                          ratio_he_z = (he/z) %>% round(2),
                                          ratio_ha_h = (ha_filter/h_filter) %>% round(2),
                                          ratio_he_h = (he/h_filter) %>% round(2))
  
  return(points_final)
} # FSLAM model function

# Multi-objective calibration
fslam_fit_multi <- function(x){
  
  points_final <- fslam(x)
  
      # Metrics (post event)
  cm <- data.frame(target = points_final$des, prediction = points_final$POF_fin_pred)
  eval <- confusionMatrix(data = factor(cm$prediction, levels = c(1,0)),
                          reference = factor(cm$target, levels = c(1,0)))
  acc_fin <- eval$overall[1] %>% as.numeric()
  #auc_fin <- auc(points_final$des, points_final$POF_fin_pred) %>% as.numeric

  
  # Metrics (pre event)
  cm <- data.frame(target = rep(0,nrow(points_final)), prediction = points_final$POF_ini_pred)
  #cm <- data.frame(target = points_final$des, prediction = points_final$POF_ini_pred)
  eval <- confusionMatrix(data = factor(cm$prediction, levels = c(1,0)),
                          reference = factor(cm$target, levels = c(1,0)))
  acc_ini <- eval$overall[1] %>% as.numeric()
  #auc_ini<- auc(points_final$des, points_final$POF_ini_pred) %>% as.numeric
  
  return(c(-acc_ini,-acc_fin))
}

# Single-objective calibration
fslam_fit_single <- function(x){
  
  points_final <- fslam(x)
  
  # Metrics (post event)
  cm <- data.frame(target = points_final$des, prediction = points_final$POF_fin_pred)
  eval <- confusionMatrix(data = factor(cm$prediction, levels = c(1,0)),
                          reference = factor(cm$target, levels = c(1,0)))
  acc_fin <- eval$overall[1] %>% as.numeric()
  #auc_fin <- auc(points_final$des, points_final$POF_fin_pred) %>% as.numeric
  
  names(acc_fin) <- 'accuracy'
  return(-acc_fin)
  
} 

# Function to get only accuracy values

get_acc <- function(x){
  
  points_final <- fslam(x)
  
  # Metrics (pre event)
  cm <- data.frame(target = rep(0,nrow(points_final)), prediction = points_final$POF_ini_pred)
  eval <- confusionMatrix(data = factor(cm$prediction, levels = c(1,0)),
                          reference = factor(cm$target, levels = c(1,0)))
  #acc_ini <- eval$overall[[1]]
  stat_ini <- c(eval$byClass[1],eval$byClass[2],eval$byClass[11])
  
  cm <- data.frame(target = points_final$des, prediction = points_final$POF_fin_pred)
  eval <- confusionMatrix(data = factor(cm$prediction, levels = c(1,0)),
                          reference = factor(cm$target, levels = c(1,0)))
  #acc_fin <- eval$overall[[1]]
  stat_fin <- c(eval$byClass[1],eval$byClass[2],eval$byClass[11])
  
  return(data.frame(stat_ini = stat_ini, stat_fin = stat_fin))
  
}

#####################################
# Multi-objective calibration results
######################################

results <- nsga2R(fn = fslam_fit, varNo = 10, objDim = 2, generations = 10, popSize = 100,  
                  lowerBounds = lower_b, upperBounds = upper_b)


df <- data.frame(ID = 1:nrow(results$objectives), 
                 acc_fin =  round(results$objectives[,2]*-1,3),
                 acc_ini =  round(results$objectives[,1]*-1,3))

# Plotting Pareto curve
plot_ly(data = df, x = ~auc_ini,y = ~auc_fin)

# Select best set of parameters by filteting objetive function

dft <- cbind(round(results$parameters,3),round(results$objectives*-1,2))%>% as.data.frame() 
colnames(dft) <- c('Cs','Cr','phi_min','phi_max','Pa','z','Pe','k','por','CN','acc_ini','acc_fin')

df_s <- dft %>% filter(acc_fin >= 0.60 & acc_ini >= 0.85) 
print(df_s) # show the filtered parameter set

#####################################
# Single-objective calibration results
######################################

p <- data.frame(name = c('Cs','Cr','phi_min','phi_max','Pa','z','Pe','k','por','CN'),
                initial = c(1,1,0,0,1,0,1,1,0,1),
                min = lower_b,
                max = upper_b)
results <- dds(fn = fslam_fit_single, p, m = 100, r = 0.2, plot = T)
print(results$f_best*-1) # Print optimal accuracy
x <- results$p_best %>% round(3) # Optimal parameter set
print(x)

# Map visualization and accuracy

fslam_map <- function(x,
                      #var,
                      crs,
                      fun,
                      export.shp = F, 
                      path = NULL){
  
  points_final <- fslam(x)
  points_final <- points_final %>% drop_na()
  
  # Convert data frame to sf object
  points_sf <- sf::st_as_sf(x = points_final, 
                        coords = c("x", "y"),
                        crs = crs)
  
  # interactive map:
  
  #library(mapview)
  
    mapviewOptions(fgb = FALSE)
    map_fin <- mapview(points_sf, zcol = c('POF_fin_pred'), legend = T, alpha = 0, 
                       col.regions = c('blue','red'),
                       na.color = 'white',
                       layer.name = 'model prediction (PoF final)',
                       cex = 4)
    map_ini <- mapview(points_sf, zcol = c('POF_ini_pred'), legend = T, alpha = 0, 
                       col.regions = c('blue','red'),
                       na.color = 'white',
                       layer.name = 'model prediction (PoF antecedent)',
                       cex = 4) 
    #print(map_ini+map_fin)
    
    mapviewOptions(fgb = FALSE)
    
    map_ha <- mapview(points_sf, zcol = c('ratio_ha_h'), legend = T, alpha = 0, 
                      na.color = 'white',
                      layer.name = 'ha/h',
                      cex = 4) 
    map_he <- mapview(points_sf, zcol = c('ratio_he_h'), legend = T, alpha = 0, 
                      na.color = 'white',
                      layer.name = 'he/h',
                      cex = 4) 
    
    #print(map_ha + map_he)

    mapviewOptions(fgb = FALSE)
    map_rha <- mapview(points_sf, zcol = c('ratio_ha_z'), legend = T, alpha = 0, 
                       na.color = 'white',
                       layer.name = 'ha/z',
                       cex = 4) 
    
    map_rhe <- mapview(points_sf, zcol = c('ratio_he_z'), legend = T, alpha = 0, 
                       na.color = 'white',
                       layer.name = 'he/z',
                       cex = 4) 
    
    print(map_ini + map_fin + map_ha + map_he + map_rha + map_rhe)

  
  acc <- get_acc(x)
  acc <- c(acc$stat_ini[2], acc$stat_fin[3])
  names(acc) <- c('Acc_pre_event','Acc_post_event')
  
  print(acc*-1)
  
  if(export.shp == T){
    
    writeOGR(points_sf %>% as('Spatial'), dsn = path, 
             layer = "points_result", 
             driver = "ESRI Shapefile",
             overwrite_layer = T)
    
    #setwd(path)
    #write.csv(points_final, file = 'Final_Table.csv',row.names = F)
    
  } 
  
} # Visualization tool


# Choose row number (e.g. 1) to apply that parameter set 
x <- df_s %>% slice(1) %>% as.numeric 

# Plot the interactive map
fslam_map(x, 
          crs = "+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs", # set projection
          export.shp = F, # shape file will be exported?
          path = '.../best_calibrated.shp') # Indicate shapefile path to be saved, only if export.shp = T

