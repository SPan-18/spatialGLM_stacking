# test posterior contraction: psi3
N <- 1000
y <- 3.5
b <- 1
z <- rpois(N, y)

prior_alpha <- 10
prior_kappa <- 1
y_prior <- rDY(n.samples = 1000, 
               alpha = prior_alpha, kappa = prior_kappa, psi = "psi3")
post_alpha <- sum(z) + prior_alpha
post_kappa <- N * b + prior_kappa
y_post <- rDY(n.samples = 1000, 
              alpha = post_alpha, kappa = post_kappa, psi = "psi3")

h2 <- density(exp(y_prior))
h3 <- density(exp(y_post))
plot(h2, col = "green", lty = 2, ylim = c(0, 2.5))
lines(h3, col = "red", lty = 2)

summary(exp(y_post))
y

# test posterior contraction: psi2
N <- 100
p <- 0.2
b <- 1
z <- rbinom(N, 1, p)

prior_alpha <- 1
prior_kappa <- 2
y_prior <- rDY(n.samples = 1000, 
               alpha = prior_alpha, kappa = prior_kappa, psi = "psi2")
post_alpha <- sum(z) + prior_alpha
post_kappa <- N * b + prior_kappa
y_post <- rDY(n.samples = 1000, 
              alpha = post_alpha, kappa = post_kappa, psi = "psi2")

logit.inv <- function(x) return(exp(x)/(1 + exp(x)))

h2 <- density(logit.inv(y_prior))
h3 <- density(logit.inv(y_post))
plot(h2, col = "green", lty = 2, ylim = c(0, 10))
lines(h3, col = "red", lty = 2)

summary(logit.inv(y_post))




