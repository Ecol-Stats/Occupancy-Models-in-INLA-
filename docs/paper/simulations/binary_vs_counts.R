library(tidyverse)
library(INLA)
library(inlabru)


# Simulate Data -----------------------------------------------------------

# Choose sample sizes and prepare observed data array y
set.seed(2023)                  # So we all get same data set
M <- 400                      # Number of sites
J <- 4                       # Number of presence/absence measurements
y <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data
x <- rnorm(M)

beta <- c(qlogis(0.4),1.5)
# Parameter values
psi <- plogis(beta[1] + beta[2]*x)  # Probability of occupancy or presence
p <- 0.3                    # Probability of detection

# Generate presence/absence data (the truth)
z <- rbinom(n = M, size = 1, prob = psi)  

# Generate detection/nondetection data (i.e. presence/absence measurements)
for(j in 1:J){
  y[,j] <- rbinom(n = M, size = 1, prob = z*p)
}


# Binary Data  ----------------------------------------------------


data_bin = data.frame(site = rep(1:M, J),
                      n = 1,
                      y = c(y),
                      x_cov = rep(x, J),
                      int_detection =1,
                      visit = rep(1:J,each=M)
)


formula1 = inla.mdata(cbind(y,n),int_detection) ~  x_cov 


model_bin <- inla(formula1, data=data_bin, family= '0binomialS',verbose = FALSE,
               control.compute = list( config = TRUE,dic  = T, waic = T),
               control.fixed = list(prec.intercept = 1/2.72,prec = 1/2.72),
               control.family = list(control.link = list(model = "logit"),
                                     link.simple = "logit",
                                     hyper = list(beta1 = list(param = c(0,1),
                                                               initial = 0),
                                                  beta2 = list(param = c(0,1)))))

#inla.rerun(model_bin)

model_bin$summary.hyperpar

plot(inla.tmarginal(function(x) -x ,model_bin$marginals.fixed$`(Intercept)`),type="l")

# Counts data analysis ----------------------------------------------------


data_wide <- data_bin %>% pivot_wider(names_from = visit,values_from = y)
data_wide <- data_wide %>% mutate(y = apply(data_wide[,5:8],1,sum),n=4)
data_wide$n <- J


formula2 = inla.mdata(cbind(y,n),int_detection) ~ x_cov 


model_counts <- inla(formula2, data=data_wide, family= '0binomialS',verbose = FALSE,
               control.compute = list( config = TRUE,dic  = T, waic = T),
               control.fixed = list(prec.intercept = 1/2.72,prec = 1/2.72),
               control.family = list(control.link = list(model = "logit"),
                                     link.simple = "logit",
                                     hyper = list(beta1 = list(param = c(0,1),
                                                               initial = 0),
                                                  beta2 = list(param = c(0,1)))))


model_counts$summary.hyperpar

# plot results ------------------------------------------------------------

results =  data.frame(rbind(data.frame(inla.tmarginal(function(x) -x ,model_counts$marginals.fixed$`(Intercept)`),
                          model="counts"),
      data.frame(inla.tmarginal(function(x) -x ,model_bin$marginals.fixed$`(Intercept)`),
                 model="binary")), par = "beta[0]", true.value = beta[1])

      
results = rbind(results,
                data.frame(rbind(data.frame(inla.tmarginal(function(x) -x ,model_counts$marginals.fixed$x_cov),
                                           model="counts"),
                                data.frame(inla.tmarginal(function(x) -x ,model_bin$marginals.fixed$x_cov),
                                           model="binary")), par = "beta[1]", true.value = beta[2])
                )

results = rbind(results,
                data.frame(rbind(data.frame(model_counts$marginals.hyperpar$`beta1 for 0binomialS observations`,
                                            model="counts"),
                                 data.frame(model_bin$marginals.hyperpar$`beta1 for 0binomialS observations`,
                                            model="binary")), par = "p", true.value = qlogis(p))
)



ggplot(data = results, aes(x,y,colour = model)) +
  geom_line() +
  geom_vline(aes(xintercept = true.value), linewidth = 0.6) +
    facet_wrap(~par,labeller = label_parsed,scales = 'free_x')




