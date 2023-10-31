# library(MASS)
# library(dplyr)
# library(ggplot2)
# library(akima)

n.grid <- 1000
# Grid <- expand.grid(x.easting = seq(0, 1, length.out = n.grid),
#                     x.northing = seq(0, 1, length.out = n.grid))
Grid <- data.frame(x = c(0,0,1,1,runif(n.grid - 4, 0, 1)),
                   y = c(0,1,0,1,runif(n.grid - 4, 0, 1)))
n <- nrow(Grid)

# Spatial field
distance <- as.matrix(dist(Grid))
phi <- 5
# omega <- MASS::mvrnorm(n     = 1,
#                        mu    = rep(0,n),
#                        Sigma = 0.4 * exp(- phi * distance))
V <- exp(- phi * distance)
V <- Rfast::cholesky(V)
omega <- rnorm(n = n, mean = 0, sd = 0.02)
omega <- crossprod(V, omega)

# omega_obj <- list(x = Grid$x.easting, y = Grid$x.northing,
#                   z = omega)
# fields_loc <- cbind(runif(100),
#                     runif(100))
# fields_omega <- interp.surface(omega_obj, fields_loc)
# image.plot(omega_obj)
# quilt.plot(fields_loc, fields_omega) 

# ggplot(Grid %>% mutate(z = omega), aes(x = x, y = y)) +
#   geom_raster(aes(fill = z)) +
#   # geom_tile(aes(fill = z)) +
#   scale_fill_viridis_c(option = "inferno", direction = -1) +
#   # scale_fill_gradient2(midpoint = 2) +
#   xlab("Easting") +
#   ylab("Northing") +
#   theme_bw() +
#   theme(axis.line = element_line(color='black'),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         aspect.ratio = 1)

akima_omega <- interp(x = Grid$x, y = Grid$y, z = omega,
                      xo = seq(0, 1, length.out = 1000),
                      yo = seq(0, 1, length.out = 1000),
                      linear = TRUE)
akima_omega <- interp2xyz(akima_omega, data.frame = TRUE)

ggplot(akima_omega, aes(x = x, y = y)) +
  geom_raster(aes(fill = z)) +
  scale_fill_viridis_c(option = "inferno", direction = -1) +
  # scale_fill_gradient2(midpoint = 2) +
  xlab("Easting") +
  ylab("Northing") +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        aspect.ratio = 1)

