# rm(list = ls())
# 
# library(Rfast)

n <- 7000
Grid <- data.frame(x = runif(n, 0, 1),
                   y = runif(n, 0, 1))

distance <- as.matrix(dist(Grid))
phi <- 3
Rphi <- exp(- phi * distance)

# t11 <- Sys.time()
# chol_base <- chol(Rphi)
# t12 <- Sys.time()
# cat("Base:\n")
# print(t12 - t11)

t21 <- Sys.time()
chol_Rfast <- Rfast::cholesky(Rphi, parallel = TRUE)
t22 <- Sys.time()
cat("Rfast (parallel = TRUE):\n")
print(t22 - t21)

t31 <- Sys.time()
chol_Rfast2 <- Rfast::cholesky(Rphi, parallel = FALSE)
t32 <- Sys.time()
cat("Rfast (parallel = FALSE):\n")
print(t32 - t31)
