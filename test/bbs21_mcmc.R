rm(list = ls())
library(tidyverse)
# library(spBayes)

ci_beta <- function(mat){
  out <- array(dim = c(ncol(mat), 5))
  for(i in 1:ncol(mat)){
    out[i, ] = c(mean(mat[, i]), sd(mat[, i]), 
                 quantile(mat[, i], probs = c(0.025, 0.5, 0.975)))
  }
  colnames(out) <- c("Mean", "SD", "2.5%", "50%", "97.5%")
  rownames(out) <- paste("beta", 1:ncol(mat) - 1, sep = '')
  return(out)
}

bbs21 <- read.csv("../data/BBS21.csv")
bbs21 <- bbs21 %>%
  filter(Latitude < 60) %>%
  filter(!(RouteDataID == 6376050)) %>%
  filter(BirdCount < 500)

ids <- 1:nrow(bbs21)
y <- as.numeric(bbs21[ids, "BirdCount"])
X <- as.matrix(bbs21[, c("Longitude", "Latitude", "NCar", "Noise")])
# X <- cbind(rep(1, length(y)), X)
S <- as.matrix(bbs21[ids, c("Longitude", "Latitude")])

# beta.starting <- coefficients(glm(y~X, family="poisson"))
# beta.tuning <- t(chol(vcov(glm(y~X, family="poisson"))))
# n.batch <- 500
# batch.length <- 50
# n.samples <- n.batch*batch.length

# t1 <- Sys.time()
# m.1 <- spGLM(y~X, family="poisson", coords=S,
#              starting=list("beta"=beta.starting, "phi"=0.06,"sigma.sq"=1, "w"=0),
#              tuning=list("beta"=beta.tuning, "phi"=0.5, "sigma.sq"=0.5, "w"=0.5),
#              priors=list("beta.Flat", "phi.Unif"=c(0.03, 1500), "sigma.sq.IG"=c(2, 1)),
#              amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
#              cov.model="exponential", verbose=TRUE, n.report=10)
# t2 <- Sys.time()
# print(t2-t1)
# burn.in <- 0.9*n.samples
# sub.samps <- burn.in:n.samples
# 
# w.hat <- m.1$p.w.samples[,sub.samps]
# beta.hat <- m.1$p.beta.theta.samples[sub.samps, 1:5]

# write.table(beta.hat, 
#             file = "bbs21_beta_hat.txt", 
#             col.names = FALSE, row.names = FALSE)
# write.table(w.hat, 
#             file = "bbs21_z_hat.txt", 
#             col.names = FALSE, row.names = FALSE)

w.hat <- read.table("bbs21_z_hat.txt")
beta.hat <- read.table("bbs21_beta_hat.txt")

print(ci_beta(beta.hat))

y.hat <- apply(exp(cbind(rep(1, nrow(X)), X) %*% t(beta.hat) + w.hat), 
               2, function(x){rpois(dim(X)[1], x)})
y.hat.mu <- apply(y.hat, 1, median)

bbs21$yhat <- log(y.hat.mu)
bbs21$postmedian_z <- apply(w.hat, 1, median)

library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(sp)
library(raster)
library(viridis)
library(MBA)

surf <- mba.surf(bbs21[ , c("Longitude","Latitude", "yhat")], 
                 no.X = 200, no.Y = 200, h = 5, m = 1, n = 1, 
                 extend=FALSE)$xyz.est

surf_rast <- raster(surf)
usa <- ne_countries(country = "United States of America", 
                    scale = "medium", returnclass = "sf")
# us_rast <- crop(surf_rast, usa)
us_rast <- raster::extract(surf_rast, usa, cellnumbers = TRUE)
full_grid <- expand.grid(surf$x, rev(surf$y))
us_df <- data.frame(full_grid[us_rast[[1]][, "cell"], ], 
                    z = us_rast[[1]][, "value"])
colnames(us_df) <- c("x", "y", "z")
us_df <- na.omit(us_df)
us_df$z <- exp(us_df$z)

world <- ne_countries(scale = "medium", returnclass = "sf")
brks <- quantile(bbs21$yhat, c(0, 0.2, 0.4, 0.6, 0.8, 1))

blue.red <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
col.br <- colorRampPalette(blue.red)

ggplot(data = world) +
  geom_sf(fill = "antiquewhite") +
  coord_sf(xlim = c(-149.40764, -39.89025), 
           ylim = c(17.06898, 76.91902), expand = FALSE) +
  # geom_raster(data = surf_df, aes(x = x, y = y, fill = z), alpha = 0.75) +
  geom_raster(data = us_df, aes(x = x, y = y, fill = z), alpha = 0.9) +
  # scale_fill_viridis(option = "plasma", trans = "log", direction = -1, 
  #                    label = function(x) sprintf("%.1f", x)) +
  scale_fill_gradientn(colours = col.br(20), trans = "log",
                       label = function(x) sprintf("%.1f", x)) +
  xlab("Longitude") +
  ylab("Latitude") +
  # labs(fill = latex2exp::TeX('y(s)$')) +
  labs(fill = "Predicted\nAvian\nCount") +
  theme_bw() +
  theme(panel.grid.major = element_line(color = gray(.8), 
                                        linetype = "dotted",
                                        linewidth = 0.25), 
        panel.background = element_rect(fill = "aliceblue"),
        legend.position = c(0.92, 0.15),
        legend.background = element_blank(),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_text(size = 8),
        legend.title.align = 0.5,
        legend.text = element_text(size = 8),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        aspect.ratio = 1)




