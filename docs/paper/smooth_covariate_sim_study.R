#' ---
#' title: "Occupancy models in INLA simulation study"
#' ---

#' First we load the libraris and helper functions
library(INLA)
library(inlabru)
library(fmesher)
library(tidyverse)
library(spatstat)
library(sf)
library(terra)
library(ggplot2)
library(gt)
library(dplyr)
library(viridis)
library(viridisLite)
library(scico)
library(patchwork)


inla.Occupancy_detCov <- function(X_det){
  
  if(class(X_det)=="list"){
    if(length(X_det)>10){
      warning("exceeded number of detection covariates, numerical issues may occur")
    }
    
    if(lapply(X_det, ncol)%>%unlist()%>%unique()%>%length()>2){
      stop("inconsistent number of visits in provided detection covariates")
    }
    if(length(lapply(X_det, nrow) %>% unlist() %>% unique())>1){
      stop("inconsistent number of sites in provided detection covariates")
    }
    K<- lapply(X_det, ncol) %>% unlist() %>% max() # Max num of visits
    M<- lapply(X_det, nrow) %>% unlist() %>% unique() # Number of sites
    P <- length(X_det)
    
    if(lapply(X_det, ncol)%>%unlist()%>%unique()%>%length()==2 & 
       1 %in% lapply(X_det, ncol)%>%unlist()%>%unique()){
      warning(paste("At least one covariate of dimension [",M,",1] has been provided, values for this covariate will be repeated over the max numver of visits",sep=""))
      for(l in which(lapply(X_det, ncol) %>% unlist() < K)){
        X_det[[l]] <- do.call("cbind",replicate(K,X_det[[l]]))
        
      }
    }
    covariates <- do.call("cbind", lapply(1:K, function(i) {
      do.call("cbind", lapply(X_det, function(mat) mat[, i]))
    }))
    
  }
  
  if(is.data.frame(X_det)|is.matrix(X_det)){
    K<- ncol(X_det)
    M<- nrow(X_det)
    P <- 1
    covariates <- as.matrix(X_det)
  }
  
  X_mat <- matrix(NA,nrow=M,ncol=K*(P+1))
  X_mat[,seq(1,(K*(P+1)),by=(P+1))]<-1 # add Intercept at the begining of each visit-specific covariate matrix
  X_mat[, which(!(1:(K*(P+1)) %in% seq(1,(K*(P+1)),by=(P+1))))] <- covariates
  return(X_mat)
  
}

hill_4P <- function(x,Int_occ,a,b,c){
  y = Int_occ + (a*x^b)/(c^b + x^b)
  y
}

#' Define spatial domain

win <- owin(c(0,300), c(0,300))
npix <- 1000
Domain <- rast(nrows=npix, ncols=npix,
               xmax=win$xrange[2],xmin=win$xrange[1],
               ymax = win$yrange[2],ymin=win$yrange[1])

values(Domain) <- 1:ncell(Domain)
xy <- crds(Domain)

#' Define a regular grid
cell_size = 3
customGrid <- st_make_grid(Domain,cellsize = c(cell_size,cell_size)) %>% 
  st_cast("MULTIPOLYGON") %>%
  st_sf() %>%
  mutate(cellid = row_number())

#' number of cells
ncells <- nrow(customGrid)


#' Spatial boundary
boundary_sf = st_bbox(c(xmin = 0, xmax = 300, ymax = 0, ymin = 300)) %>%
  st_as_sfc()

#' Create a fine mesh for simulating the random field
mesh_sim = fm_mesh_2d(loc.domain = st_coordinates(boundary_sf)[,1:2],
                      offset = c(-0.1, -.2),
                      max.edge = c(4, 50))

#' Matern model
matern_sim <- inla.spde2.pcmatern(mesh_sim,
                                  prior.range = c(100, 0.5),
                                  prior.sigma = c(1, 0.5))


range_spde = 100
sigma_spde = 1

#' Precision matrix
Q = inla.spde.precision(matern_sim, theta = c(log(range_spde),
                                              log(sigma_spde)))
#' Simulate three spatial fields
seed = 12345
sim_field = inla.qsample(n = 2, Q = Q, seed = seed)


#' Obtain the centroid of each cell
coord_grid  = st_coordinates(customGrid %>% st_centroid())
#' A projector matrix
A_proj = inla.spde.make.A(mesh_sim, loc = coord_grid)

#' Spatial components
omega_s = (A_proj %*% sim_field)[,1] # spatial random field
x_s = (A_proj %*% sim_field)[,2] # spatial environmental covariate

#' plot linear predictor (withour spatial effect)
#+  fig.width=5, fig.height=5
ggplot(data.frame(x=x_s,y =  hill_4P(x_s,Int_occ = qlogis(0.4),a = 5,b=3,c=-3)  ),aes(x=x,y=plogis(y)))+ geom_line()

#' create rasters
x_rast = rast(data.frame(x = coord_grid[,1], y = coord_grid[,2],x_s))
#+ fig.width=5, fig.height=5
plot(x_rast)

#' Occupancy probabilities
psi <- inla.link.logit(hill_4P(x_s,Int_occ = qlogis(0.4),a = 5,b=3,c=-3) + omega_s , inverse = T)
psi_rast = rast(data.frame(x = coord_grid[,1], y = coord_grid[,2],psi))

#+  fig.width=5, fig.height=5
plot(psi_rast)
#' True occupancy state

set.seed(seed)
z <- rbinom(ncells, size = 1, prob = psi)


#' number of cells/sites in the sample
nsites = round(ncells*.40)
site_id = sample(1:ncells, size=nsites, replace=FALSE) # cell id
#' add an indicator of whether a cell is in the sample or not
customGrid$sample <- ifelse(customGrid$cellid%in%site_id,1,0)


#' Number of vistis
K= 3
#' Observational process model coeficcients
alpha <- c(NA,NA)
alpha[1] <- qlogis(0.6) # Base line detection probability
alpha[2] <- 1 # detection covariate g1 effect


#' Detection probabilities and observed occurrences
y <- p.mat <- matrix(NA,nrow = nsites,ncol=3) # Create empty matrix to store the results

g2 <- array(runif(n = nsites * K, -1, 1), dim = c(nsites, K)) # detection covariate

# loop over visits
for(j in 1:K){
  p.mat[,j] <- inla.link.logit(alpha[1] + 
                                 alpha[2]*g2[, j], inverse = T)
  y[,j] <- rbinom(n = nsites,size = 1,prob = p.mat[,j]*z[site_id] )
}

#' Get data sets

Occ_data_1 <- customGrid %>%
  st_centroid() %>%
  filter(sample==1) %>%
  dplyr::select(-c('sample'))

Obs_data <- data.frame(y = y, # detectio/non-detection data
                       g2 = g2, # survey level covariate
                       cellid = site_id)

SSOM <- left_join(Occ_data_1,Obs_data,by = "cellid") 

SSOM = SSOM %>%
  dplyr::select(-cellid) %>% 
  mutate(terra::extract(x_rast,st_coordinates(SSOM)))

#' append coordinates as columns

SSOM[,c('x.loc','y.loc')] <- st_coordinates(SSOM) 
SSOM <- SSOM %>% st_drop_geometry(SSOM)


#' Get occurrence data, coordinates, detection covariates and group site-level covariates
#' 
Y_mat <- SSOM %>% dplyr::select(num_range("y.",1:3)) %>% as.matrix()
XY_coords <- SSOM %>% dplyr::select(c(x.loc, y.loc)) %>% as.matrix()
X_occ <- SSOM %>% dplyr::select(x_s) %>%
  mutate(group_xs = inla.group(x_s),
         Int_occ =1, # Baseline occupancy porbability intercept
         spatial_field = rep(NA,nrow(SSOM))) # A matrix is supplied in the f() function (see online supplementary material)
X_det <-  SSOM %>% dplyr::select(num_range("g2.",1:3)) %>% 
  inla.Occupancy_detCov() # detection helper function (see details on online supplementary material) 


#' Define SPDE model
mesh = fm_mesh_2d(loc.domain = st_coordinates(boundary_sf)[,1:2],
                  offset = c(-0.1, -.2),
                  max.edge = c(15, 30))
matern <- inla.spde2.pcmatern(mesh,
                              prior.range = c(100, 0.5),
                              prior.sigma = c(1, 0.5))

#' A projector matrix
A_sp <- inla.spde.make.A(mesh = mesh,loc = XY_coords)

#' Data list
data_list = as.list(X_occ)
data_list$Y = Y_mat
data_list$X = X_det

#' model formula
formula_occ <- inla.mdata(Y,X) ~  -1 + Int_occ +  f(group_xs, model = "rw2") +  f(spatial_field, model=matern,  A.local = A_sp)
#' Run the model
model_occ <- inla(formula_occ,    # model formula
                    data= data_list,        # data 
                    family= 'occupancy', # model likelihood
                    # priors
                    control.fixed =  list(prec = 1/2.72, prec.intercept = 1/2.72),
                    # compute WAIC and DIC
                    control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                    verbose = FALSE,
                    # choose link functions for:
                    # (i) the state process (control.link)
                    # (ii) the observation process (link.simple)
                    control.family = list(control.link = list(model = "logit"),
                                          link.simple = "logit",
                                          # priors for hyperparameters
                                          hyper = list(
                                            beta1 = list(param = c(0,1), initial = -1),
                                            beta2 = list(param = c(0,1/2.72))
                                          )
                    ))


#' Plot model results
#' 
p2 = model_ssom2$summary.random$group_xs %>% 
  ggplot() +
  geom_line(aes(ID, y = mean)) +
  geom_ribbon(aes(ID, ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.5) +
  geom_line(data = data.frame(x = x_s, 
                              y =   hill_4P(x_s,Int_occ = 0,a = 5,b=3,c=-3)),
            aes(x,y), color = "red")
p2


