varc <- c(0.2, 0.5, 0.8)
s <- 3
# rho <- varc/(varc + pi^2/3) #logit
rho <- varc/(varc + 1) #probit
print(rho)
time <- 5
n <- 50
m <- 192

# DSDo
m * (n*rho * time /(2*(1 - rho + n*rho*time)))
m * (1 - rho + n*rho*s)/(2*(1-rho + n* rho *time))
m * (1 - rho + n*rho*(time-s))/ (2*(1-rho + n* rho *time))


# USDo
m* ((1 - rho + n*rho)/(2*(1 - rho + n * rho * time)))
m* n*rho/(1 - rho + n*rho*time)

