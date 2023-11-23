library(INLA)
library(inlabru)
library(fmesher)
library(tidyverse)
library(sf)
library(terra)
library(dplyr)

###########################################################################
# Simulate data -----------------------------------------------------------
###########################################################################


# Define spatial domain
win <- owin(c(0,300), c(0,300))
npix <- 1000
Domain <- rast(nrows=npix, ncols=npix,
               xmax=win$xrange[2],xmin=win$xrange[1],
               ymax = win$yrange[2],ymin=win$yrange[1])
values(Domain) <- 1:ncell(Domain)
xy <- crds(Domain)

# Define regular grid
cell_size = 3
customGrid <- st_make_grid(Domain,cellsize = c(cell_size,cell_size)) %>% 
  st_cast("MULTIPOLYGON") %>%
  st_sf() %>%
  mutate(cellid = row_number())

# number of cells
ncells <- nrow(customGrid)

# Spatial boundary
boundary_sf = st_bbox(c(xmin = 0, xmax = 300, ymax = 0, ymin = 300)) |>
  st_as_sfc()
# Create a fine mesh
mesh_sim = fm_mesh_2d(loc.domain = st_coordinates(boundary_sf)[,1:2],
                      offset = c(-0.1, -.2),
                      max.edge = c(4, 50))
# Matern model
matern_sim <- inla.spde2.pcmatern(mesh_sim,
                                  prior.range = c(100, 0.5),
                                  prior.sigma = c(1, 0.5))

range_spde = 100
sigma_spde = 1

# Precision matrix
Q = inla.spde.precision(matern_sim, theta = c(log(range_spde),
                                              log(sigma_spde)))
# Simulate three spatial fields
seed = 12345
sim_field = inla.qsample(n = 3, Q = Q, seed = seed)

# Obtain the centroid of each cell
coord_grid  = st_coordinates(customGrid |> st_centroid())
# A matrix
A_proj = inla.spde.make.A(mesh_sim, loc = coord_grid)

# Spatial covariates 
x_s = (A_proj %*% sim_field)[,2] # spatial environmental covariate
g_s = (A_proj %*% sim_field)[,3] # spatial detection covariate

# create rasters
customGrid$xs <- x_s

x_covariate = rast(data.frame(x = coord_grid[,1], y = coord_grid[,2],x_s))
g_covariate = rast(data.frame(x = coord_grid[,1], y = coord_grid[,2],g_s)) 


beta <- c(NA,NA,NA)
beta[1] <- qlogis(0.35) # Base line occupancy probability
beta[2] <- -1.75 # environmental covariate effect
beta[3] <- -1 # environmental covariate quadratic effect

# Occupancy probabilities
psi <- inla.link.logit(beta[1] -scale(x_s**2) , inverse = T)
#psi_rast = rast(data.frame(x = coord_grid[,1], y = coord_grid[,2],qlogis(psi))) 
plot(x_s,qlogis(psi))
# True occupancy state

set.seed(seed)
z <- rbinom(ncells, size = 1, prob = psi)

# number of cells/sites in the sample
nsites = round(ncells*.20)
site_id = sample(1:ncells, size=nsites, replace=FALSE) # cell id
# add an indicator of whether a cell is in the sample or not
customGrid$sample <- ifelse(customGrid$cellid%in%site_id,1,0)

min_nvisits = 1 # minimum number of visits
max_nvisits = 5 # maximum number of visits
# Probabiliies of drawing 1 thru 5 visits per site
probs = rep(1/length(min_nvisits:max_nvisits),length(min_nvisits:max_nvisits))
# Number of visits
nvisits = sample(min_nvisits:max_nvisits,nsites, prob = probs, replace = T)

# Observational process model coeficcients
alpha <- c(NA,NA)
alpha[1] <- qlogis(0.6) # Base line detection probability
alpha[2] <- 1 # detection covariate effect

# Detection probabilities
p <- inla.link.logit(alpha[1] + alpha[2]*g_s[site_id], inverse = T)
y <- rbinom(n = nsites,size = nvisits,prob = p*z[site_id] )


Occ_data_1 <- customGrid |>
  st_centroid() |>
  filter(sample==1) |>
  dplyr::select(-c('sample'))

y_counts = data.frame(y = y , cellid = site_id, nvisits = nvisits)

SSOM <- left_join(Occ_data_1,y_counts,by = "cellid")

# ggplot()+tidyterra::geom_spatraster(data=psi_rast)+
#   geom_sf(data=SSOM,aes(colour=factor(ifelse(y>0,1,0))))+scale_fill_viridis()


###########################################################################
# Analysis ----------------------------------------------------------------
###########################################################################

 # scale covariate (the whole raster)

x_covariate <- x_covariate %>% scale()


# grouped covariate values for smoothing

x_covariate$group_xs <- inla.group(values(x_covariate),n=25)
plot(x_covariate)

# Extract the covariate values

# evaluate covariates at each cell 

SSOM = SSOM %>% dplyr::select(-cellid) %>% 
  mutate(terra::extract(x_covariate$x_s,st_coordinates(SSOM)),
         terra::extract(g_covariate,st_coordinates(SSOM)),
         terra::extract(x_covariate$group_xs,st_coordinates(SSOM)))

val_xs = sort(unique(SSOM$group_xs)) 


mesh = fm_mesh_2d(loc.domain = st_coordinates(boundary_sf)[,1:2],
                  offset = c(-0.1, -.2),
                  max.edge = c(15, 30))
matern <- inla.spde2.pcmatern(mesh,
                              prior.range = c(100, 0.5),  
                              prior.sigma = c(1, 0.5))


# projector matrix A
A_sp <- inla.spde.make.A(mesh = mesh, 
                         loc = st_coordinates(SSOM))
# index set
iset_sp <- inla.spde.make.index(name = "spatial_field", matern$n.spde)

# build the stack
stk <- inla.stack(data=list(Ycounts = SSOM$y, # observed occurrences
                            Ncounts = SSOM$nvisits, # number of visits
                            det_cov = SSOM$g_s, # detection covariate
                            Int_det = rep(1,length(SSOM$y))), # Det Intercept
                  A=list(A_sp,1),  # the A matrix; the 1 is included to make the list(covariates)
                  effects=list(c(list(Int_occ=1), # Occ Intercept
                                 iset_sp),  #the spatial index
                               # the covariates
                               list(occ_cov = SSOM$x_s,
                                    group_xs = SSOM$group_xs)), 
                  #this is a quick name so yo can call upon easily
                  tag='ssom')


formula_SI <- inla.mdata(cbind(Ycounts,Ncounts),Int_det,det_cov) ~  -1 + 
    f(group_xs,model = "rw2", values = val_xs,scale.model = TRUE,constr = FALSE) 

model_SI <- inla(formula_SI, # model formula
                 data=inla.stack.data(stk), # data stack
                 family= '0binomialS', # model likelihood
                 # priors 
                 control.fixed =  list(prec = 1/2.72, prec.intercept = 1/2.72),
                 # matrix of predictors
                 control.predictor=list(A=inla.stack.A(stk),compute=TRUE),
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
                                         beta2 = list(param = c(0,1/2.72)))
                 )
)

model_SI$summary.fixed

adjust_int = function(x)
  x-log(cell_size^2)


data.frame(inla.tmarginal(
  adjust_int, model_SI$marginals.fixed$Int_occ),
  par = "beta[0]", true.value = beta[1]) %>%
ggplot(aes(x,y,colour = par)) +
  geom_line() +
  geom_vline(aes(xintercept = true.value), linewidth = 0.6) +
  facet_wrap(~par,labeller = label_parsed,scales = 'free_x')+
  theme(legend.position = 0)


sample = inla.posterior.sample(1000, model_SI)

ff = function(...)
   group_xs
out = inla.posterior.sample.eval(ff, sample)

data.frame(model_SI$summary.random$group_xs) %>%
  mutate(m1 = apply(out,1,mean),
         q1 = apply(out,1,quantile, 0.025),
         q2 = apply(out,1,quantile,0.975)) %>%
  mutate(eta.true = (beta[1] -scale(ID**2)) ) %>%
  ggplot() +
   geom_ribbon(aes(ID, ymin = -q1, ymax = -q2),
               alpha = 0.25,fill="tomato")+
  geom_line(aes(ID, -m1,color="estimated")) + 
  geom_line(aes(ID, eta.true,color="true"))
