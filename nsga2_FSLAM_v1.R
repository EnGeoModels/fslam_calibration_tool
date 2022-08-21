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

# For MORLE conditions
fslam_fit <- function(x){
  
  # Editing parameters
  
  nrow_soil <- nrow(soil)
  nrow_hmtu <- nrow(hmtu)
  
  soil$h[2:nrow_soil] <- x[6]
  
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
    
    points[which(points$soil==i),]$Cs_m <- Cs[i]
    
 } # effective cohesion
 
 points$Cs_sd <- NA
 for(i in 1:length(Cssd)){
   
   points[which(points$soil==i),]$Cs_sd <- Cssd[i]
   
 } 
 
 points$Cr_m <- NA
 for(i in 1:length(Cr)){
   #i <- 5
   points[which(points$lucl==i),]$Cr_m <- Cr[i]
   
 } # apparent cohesion
 
 points$Cr_sd <- NA
 for(i in 1:length(Crsd)){
   
   points[which(points$lucl==i),]$Cr_sd <- Crsd[i]
   
 } 

 points$tanphi_m <- NA
 for(i in 1:length(tanphi)){
   
   points[which(points$soil==i),]$tanphi_m <- tanphi[i]
   
 } # tangent of internal friction angle
 
 points$tanphi_sd <- NA
 for(i in 1:length(tanphisd)){
   
   points[which(points$soil==i),]$tanphi_sd <- tanphisd[i]
   
 } 
 
 points$K <- NA
 for(i in 1:(nrow_soil-1)){
   
   points[which(points$soil==i),]$K <- (soil$Ks[2:nrow_soil] %>% as.numeric)[i]
   
 } # permeability
 
 points$z <- NA
 for(i in 1:(nrow_soil-1)){
   
   points[which(points$soil==i),]$z <- (soil$h[2:nrow_soil] %>% as.numeric)[i]
   
 } # soil depth
 
 points$n <- NA
 for(i in 1:(nrow_soil-1)){
   
   points[which(points$soil==i),]$n <- (soil$porosity[2:nrow_soil] %>% as.numeric)[i]
   
 } # porosity
 
 points$hsg <- NA
 for(i in 1:(nrow_soil-1)){
   
   points[which(points$soil==i),]$hsg <- (soil$hsg[2:nrow_soil])[i]
   
 } # hydrologic soil group
 
  # Updating Pa and Pe values
  points$Pa <- points$Pa*x[5] # effective infiltration (mm/day)
  points$Pe <- points$Pe*x[7] # Rainfall event (mm/day)
  
  # Preliminary updated table
  points_soil <- points

  #summary(points_soil)
 
  # curve number
  CN <- points_soil$CN
  
  # Estimating qe (recharge from event)
  Ia <- (5080/CN-50.8)
  
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
                                         he = (ha_filter+qe/1000/n),
                                         he_filter = ifelse(he>z,z,he),
                                         h = he_filter, # total water table in m (Pa+Pe)
                                         
                                         A = z*ds*g*sin(2*slope)/2,
                                         
                                         # SF post antecedent recharge
                                         D_ini = tan(slope)/(1-(ha_filter/z)*(dw/ds)),
                                         SF_mean_ini = (tanphi/D_ini + C/A),
                                         SF_sd_ini = sqrt(tanphi_sd^2/D_ini^2+(Cr_sd^2+Cs_sd^2)/A^2), # in paper
                                         
                                         # SF post event
                                         D_fin = tan(slope)/(1-(h/z)*(dw/ds)),
                                         SF_mean_fin = (tanphi/D_fin + C/A),
                                         SF_sd_fin = sqrt(tanphi_sd^2/D_fin^2+(Cr_sd^2+Cs_sd^2)/A^2), # in paper
                                         
                                         # Landslide prediction (Safe Factor)
                                         SF_ini_pred = ifelse(SF_mean_ini > 1,0,1),
                                         SF_fin_pred = ifelse(SF_mean_fin > 1,0,1),
                                         
                                         # Landslide prediction (Probability of Failure)
                                         POF_ini = pnorm(1, SF_mean_ini, SF_sd_ini),
                                         POF_fin = pnorm(1, SF_mean_fin, SF_sd_fin),
                                         POF_ini_pred = ifelse(POF_ini > 0.5,1,0),
                                         POF_fin_pred = ifelse(POF_fin > 0.5,1,0)
  )
  
  
  
  #summary(points_final)
  points_final <- points_final %>% dplyr::select(x=x, y=y, des, ha = ha_filter, h,
                                                 SF_ini = SF_mean_ini, SF_fin = SF_mean_fin,
                                                 SF_ini_pred, SF_fin_pred, POF_ini, POF_fin,
                                                 POF_ini_pred, POF_fin_pred)
  
  # Metrics (post event)
  cm <- data.frame(target = points_final$des, prediction = points_final$POF_fin_pred)
  eval <- confusionMatrix(data = factor(cm$prediction, levels = c(1,0)),
                          reference = factor(cm$target, levels = c(1,0)))
  acc_fin <- eval$overall[1] %>% as.numeric()
  #auc_fin <- auc(points_final$des, points_final$POF_fin_pred) %>% as.numeric

  
  # Metrics (pre event)
  cm <- data.frame(target = rep(0,nrow(points_final)), prediction = points_final$POF_ini_pred)
  eval <- confusionMatrix(data = factor(cm$prediction, levels = c(1,0)),
                          reference = factor(cm$target, levels = c(1,0)))
  acc_ini <- eval$overall[1] %>% as.numeric()
  #auc_ini<- auc(points_final$des, points_final$POF_ini_pred) %>% as.numeric
  
  return(c(-acc_ini,-acc_fin))

}

# For no-MORLE conditions
fslam_fit_no_morle <- function(x){
  
  # Editing parameters
  
  nrow_soil <- nrow(soil)
  nrow_hmtu <- nrow(hmtu)
  
  soil$h[2:nrow_soil] <- x[6]
  
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
    
    points[which(points$soil==i),]$Cs_m <- Cs[i]
    
  } # effective cohesion
  
  points$Cs_sd <- NA
  for(i in 1:length(Cssd)){
    
    points[which(points$soil==i),]$Cs_sd <- Cssd[i]
    
  } 
  
  points$Cr_m <- NA
  for(i in 1:length(Cr)){
    #i <- 5
    points[which(points$lucl==i),]$Cr_m <- Cr[i]
    
  } # apparent cohesion
  
  points$Cr_sd <- NA
  for(i in 1:length(Crsd)){
    
    points[which(points$lucl==i),]$Cr_sd <- Crsd[i]
    
  } 
  
  points$tanphi_m <- NA
  for(i in 1:length(tanphi)){
    
    points[which(points$soil==i),]$tanphi_m <- tanphi[i]
    
  } # tangent of internal friction angle
  
  points$tanphi_sd <- NA
  for(i in 1:length(tanphisd)){
    
    points[which(points$soil==i),]$tanphi_sd <- tanphisd[i]
    
  } 
  
  points$K <- NA
  for(i in 1:(nrow_soil-1)){
    
    points[which(points$soil==i),]$K <- (soil$Ks[2:nrow_soil] %>% as.numeric)[i]
    
  } # permeability
  
  points$z <- NA
  for(i in 1:(nrow_soil-1)){
    
    points[which(points$soil==i),]$z <- (soil$h[2:nrow_soil] %>% as.numeric)[i]
    
  } # soil depth
  
  points$n <- NA
  for(i in 1:(nrow_soil-1)){
    
    points[which(points$soil==i),]$n <- (soil$porosity[2:nrow_soil] %>% as.numeric)[i]
    
  } # porosity
  
  points$hsg <- NA
  for(i in 1:(nrow_soil-1)){
    
    points[which(points$soil==i),]$hsg <- (soil$hsg[2:nrow_soil])[i]
    
  } # hydrologic soil group
  
  # Updating Pa and Pe values
  points$Pa <- points$Pa*x[5] # effective infiltration (mm/day)
  points$Pe <- points$Pe*x[7] # Rainfall event (mm/day)
  
  # Preliminary updated table
  points_soil <- points
  
  #summary(points_soil)
  
  # curve number
  CN <- points_soil$CN
  
  # Estimating qe (recharge from event)
  Ia <- (5080/CN-50.8)
  
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
                                         he = (ha_filter+qe/1000/n),
                                         he_filter = ifelse(he>z,z,he),
                                         h = he_filter, # total water table in m (Pa+Pe)
                                         
                                         A = z*ds*g*sin(2*slope)/2,
                                         
                                         # SF post antecedent recharge
                                         D_ini = tan(slope)/(1-(ha_filter/z)*(dw/ds)),
                                         SF_mean_ini = (tanphi/D_ini + C/A),
                                         SF_sd_ini = sqrt(tanphi_sd^2/D_ini^2+(Cr_sd^2+Cs_sd^2)/A^2), # in paper
                                         
                                         # SF post event
                                         D_fin = tan(slope)/(1-(h/z)*(dw/ds)),
                                         SF_mean_fin = (tanphi/D_fin + C/A),
                                         SF_sd_fin = sqrt(tanphi_sd^2/D_fin^2+(Cr_sd^2+Cs_sd^2)/A^2), # in paper
                                         
                                         # Landslide prediction (Safe Factor)
                                         SF_ini_pred = ifelse(SF_mean_ini > 1,0,1),
                                         SF_fin_pred = ifelse(SF_mean_fin > 1,0,1),
                                         
                                         # Landslide prediction (Probability of Failure)
                                         POF_ini = pnorm(1, SF_mean_ini, SF_sd_ini),
                                         POF_fin = pnorm(1, SF_mean_fin, SF_sd_fin),
                                         POF_ini_pred = ifelse(POF_ini > 0.5,1,0),
                                         POF_fin_pred = ifelse(POF_fin > 0.5,1,0)
  )
  
  
  
  #summary(points_final)
  points_final <- points_final %>% dplyr::select(x=x, y=y, des, ha = ha_filter, h,
                                                 SF_ini = SF_mean_ini, SF_fin = SF_mean_fin,
                                                 SF_ini_pred, SF_fin_pred, POF_ini, POF_fin,
                                                 POF_ini_pred, POF_fin_pred)
  
  # Metrics (post event)
  cm <- data.frame(target = rep(0,nrow(points_final)), prediction = points_final$POF_fin_pred)
  eval <- confusionMatrix(data = factor(cm$prediction, levels = c(1,0)),
                          reference = factor(cm$target, levels = c(1,0)))
  acc_fin <- eval$overall[1] %>% as.numeric()
  auc_fin <- auc(points_final$des, points_final$POF_fin_pred) %>% as.numeric
  
  
  # Metrics (pre event)
  cm <- data.frame(target = rep(0,nrow(points_final)), prediction = points_final$POF_ini_pred)
  eval <- confusionMatrix(data = factor(cm$prediction, levels = c(1,0)),
                          reference = factor(cm$target, levels = c(1,0)))
  acc_ini <- eval$overall[1] %>% as.numeric()
  auc_ini<- auc(points_final$des, points_final$POF_ini_pred) %>% as.numeric
  
  return(c(-acc_ini,-acc_fin))
  
}

# Multi-objective calibration results (Change function in fn accordingly to assess MORLE or no-MORLE) 

results <- nsga2R(fn = fslam_fit, varNo = 10, objDim = 2, generations = 10, popSize = 100,  
                  lowerBounds = lower_b, upperBounds = upper_b)


df <- data.frame(ID = 1:nrow(results$objectives), 
                 acc_fin =  round(results$objectives[,2]*-1,3),
                 acc_ini =  round(results$objectives[,1]*-1,3))

# Plotting Pareto curve
plot_ly(data = df, x = ~auc_ini,y = ~auc_fin)

# Select best set of parameters by filteting objetive function
id <- df %>% filter(acc_fin >= 0.6 & acc_ini >= 0.9) %>% dplyr::select(ID)
id <- id$ID
r_df <- round(results$parameters,3)[c(id),] %>% as.data.frame() 
colnames(r_df) <- c('Cs','Cr','phi_min','phi_max','Pa','z','Pe','k','por','CN')
r_df <- r_df %>% mutate(df %>% filter(acc_fin >= 0.6 & acc_ini >= 0.9)) %>% dplyr::select(-ID)
r_df # show the filtered parameter set


# Map visualization and accuracy

fslam_map <- function(x,var,
                      export.shp = F, 
                      path = NULL){
  
  # Editing parameters
  
  nrow_soil <- nrow(soil)
  nrow_hmtu <- nrow(hmtu)
  
  soil$h[2:nrow_soil] <- x[6]
  
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
    
    points[which(points$soil==i),]$Cs_m <- Cs[i]
    
  } # effective cohesion
  
  points$Cs_sd <- NA
  for(i in 1:length(Cssd)){
    
    points[which(points$soil==i),]$Cs_sd <- Cssd[i]
    
  } 
  
  points$Cr_m <- NA
  for(i in 1:length(Cr)){
    #i <- 5
    points[which(points$lucl==i),]$Cr_m <- Cr[i]
    
  } # apparent cohesion
  
  points$Cr_sd <- NA
  for(i in 1:length(Crsd)){
    
    points[which(points$lucl==i),]$Cr_sd <- Crsd[i]
    
  } 
  
  points$tanphi_m <- NA
  for(i in 1:length(tanphi)){
    
    points[which(points$soil==i),]$tanphi_m <- tanphi[i]
    
  } # tangent of internal friction angle
  
  points$tanphi_sd <- NA
  for(i in 1:length(tanphisd)){
    
    points[which(points$soil==i),]$tanphi_sd <- tanphisd[i]
    
  } 
  
  points$K <- NA
  for(i in 1:(nrow_soil-1)){
    
    points[which(points$soil==i),]$K <- (soil$Ks[2:nrow_soil] %>% as.numeric)[i]
    
  } # permeability
  
  points$z <- NA
  for(i in 1:(nrow_soil-1)){
    
    points[which(points$soil==i),]$z <- (soil$h[2:nrow_soil] %>% as.numeric)[i]
    
  } # soil depth
  
  points$n <- NA
  for(i in 1:(nrow_soil-1)){
    
    points[which(points$soil==i),]$n <- (soil$porosity[2:nrow_soil] %>% as.numeric)[i]
    
  } # porosity
  
  points$hsg <- NA
  for(i in 1:(nrow_soil-1)){
    
    points[which(points$soil==i),]$hsg <- (soil$hsg[2:nrow_soil])[i]
    
  } # hydrologic soil group
  
  # Updating Pa and Pe values
  points$Pa <- points$Pa*x[5] # effective infiltration (mm/day)
  points$Pe <- points$Pe*x[7] # Rainfall event (mm/day)
  
  # Preliminary updated table
  points_soil <- points
  
  #summary(points_soil)
  
  # curve number
  CN <- points_soil$CN
  
  # Estimating qe (recharge from event)
  Ia <- (5080/CN-50.8)
  
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
                                         he = (ha_filter+qe/1000/n),
                                         he_filter = ifelse(he>z,z,he),
                                         h = he_filter, # total water table in m (Pa+Pe)
                                         
                                         A = z*ds*g*sin(2*slope)/2,
                                         
                                         # SF post antecedent recharge
                                         D_ini = tan(slope)/(1-(ha_filter/z)*(dw/ds)),
                                         SF_mean_ini = (tanphi/D_ini + C/A),
                                         SF_sd_ini = sqrt(tanphi_sd^2/D_ini^2+(Cr_sd^2+Cs_sd^2)/A^2), # in paper
                                         
                                         # SF post event
                                         D_fin = tan(slope)/(1-(h/z)*(dw/ds)),
                                         SF_mean_fin = (tanphi/D_fin + C/A),
                                         SF_sd_fin = sqrt(tanphi_sd^2/D_fin^2+(Cr_sd^2+Cs_sd^2)/A^2), # in paper
                                         
                                         # Landslide prediction (Safe Factor)
                                         SF_ini_pred = ifelse(SF_mean_ini > 1,0,1),
                                         SF_fin_pred = ifelse(SF_mean_fin > 1,0,1),
                                         
                                         # Landslide prediction (Probability of Failure)
                                         POF_ini = pnorm(1, SF_mean_ini, SF_sd_ini),
                                         POF_fin = pnorm(1, SF_mean_fin, SF_sd_fin),
                                         POF_ini_pred = ifelse(POF_ini > 0.5,1,0),
                                         POF_fin_pred = ifelse(POF_fin > 0.5,1,0)
  )
  
  
  #summary(points_final)
  points_final <- points_final %>% dplyr::select(x=x, y=y, des, ha = ha_filter, h,
                                                 SF_ini = SF_mean_ini, SF_fin = SF_mean_fin,
                                                 SF_ini_pred, SF_fin_pred, POF_ini, POF_fin,
                                                 POF_ini_pred, POF_fin_pred)
  
  # Convert data frame to sf object
  points_sf <- st_as_sf(x = points_final, 
                        coords = c("x", "y"),
                        crs = "+proj=utm +zone=51 +datum=WGS84 +units=m +no_defs")
  
  # interactive map:
  
  library(mapview)
  
  if(var=='POF_fin_pred'){
    mapviewOptions(fgb = FALSE)
    map <- mapview(points_sf, zcol = c('POF_fin_pred'), legend = T, alpha = 0, 
                   col.regions = c('blue','red'),
                   na.color = 'white',
                   layer.name = 'model prediction (PoF final)',
                   cex = 4) 
    print(map)
    
  } else if(var=='POF_ini_pred') {
    mapviewOptions(fgb = FALSE)
    map <- mapview(points_sf, zcol = c('POF_ini_pred'), legend = T, alpha = 0, 
                   col.regions = c('blue','red'),
                   na.color = 'white',
                   layer.name = 'model prediction (PoF initial)',
                   cex = 4) 
    print(map)
  }
  
  # Metrics (post event)
  cm <- data.frame(target = points_final$des, prediction = points_final$POF_fin_pred)
  eval <- confusionMatrix(data = factor(cm$prediction, levels = c(1,0)),
                          reference = factor(cm$target, levels = c(1,0)))
  acc_fin <- eval$overall[1] %>% as.numeric()
  #auc_fin <- auc(points_final$des, points_final$POF_fin_pred) %>% as.numeric
  
  # Metrics (pre event)
  cm <- data.frame(target = rep(0,nrow(points_final)), prediction = points_final$POF_ini_pred)
  eval <- confusionMatrix(data = factor(cm$prediction, levels = c(1,0)),
                          reference = factor(cm$target, levels = c(1,0)))
  acc_ini <- eval$overall[1] %>% as.numeric()
  #auc_ini<- auc(points_final$des, points_final$POF_ini_pred) %>% as.numeric
  
  acc <- c(acc_ini,acc_fin)
  
  names(acc) <- c('Acc_pre_event','Acc_post_event')
  print(acc)
  
  if(export.shp == T){
    
    writeOGR(points_sf %>% as('Spatial'), dsn = path, 
             layer = "points_result", 
             driver = "ESRI Shapefile",
             overwrite_layer = T)
    
    setwd(path)
    write.csv(points_final, file = 'Final_Table.csv',row.names = F)
    
  } 
  
}

# Choose row number (e.g. 18) to apply that parameter set 
x <- r_df %>% slice(1) %>% as.numeric

# Choose conditions to be plotted (POF_fin_pred or POF_ini_pred')
var <- 'POF_ini_pred'

fslam_map(x,var, 
          export.shp = F, # shapefile will not be exported
          path = 'C:/disco_D/tesis_maestria/06_FSLAM_inR/shp') # Indicate shapefile path to be saved, only if export.shp = T

