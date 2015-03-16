#!/usr/bin/env Rscript

set.seed(666)
identity <- F
exact <- F

library("MASS")

lpr_mvn <- function (value, grad=FALSE, inv.cov)
{
  g <- - inv.cov %*% value
  r <- as.vector (t(value) %*% g / 2)
  if (grad)
  { attr(r,"grad") <- as.vector(g)
    if (is.na(as.vector(g))) write("ERRORG", stderr())
  }
  if (is.na(as.vector(r))) write("ERROR", stderr())

  r
}

logmvdnorm <- function (x, mu, Sigma)
{
    d <- x - mu
    - length(x)/2 * log(2 * pi) - log(det(Sigma))/2 - d %*% solve(Sigma) %*% d / 2
}

grad_logmvdnorm <- function (x, mu, Sigma) - c(solve(Sigma) %*% (x - mu))

HMC <- function (U, grad_U, epsilon, L, mass, current_q)
{
    q <- current_q
    p <- mvrnorm(1, rep.int(0, length(q)), mass) # multivariate normal variates
    current_p <- p
    # Make a half step for momentum at the beginning
    p <- p - epsilon * grad_U(q) / 2
    # Alternate full steps for position and momentum
    for (i in 1:L)
    {
        # Make a full step for the position
        q <- q + epsilon * p
        # Make a full step for the momentum, except at end of trajectory
        if (i < L) p <- p - epsilon * grad_U(q)
    }
    # Make a half step for momentum at the end.
    p <- p - epsilon * grad_U(q) / 2
    # Negate momentum at end of trajectory to make the proposal symmetric
    p <- -p
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_U <- U(current_q)
    current_K <- -logmvdnorm(current_p, rep.int(0, length(q)), mass)
    proposed_U <- U(q)
    proposed_K <- -logmvdnorm(p, rep.int(0, length(q)), mass)
    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    if (log(runif(1)) < (current_U-proposed_U+current_K-proposed_K))
    {
        return (list(q = q, a = T))  # accept
    }
    else
    {
        return (list(q = current_q, a = F)) # reject
    }
}

exact_HMC <- function (t, mu, Sigma, q)
{
    a <- mvrnorm(1, rep.int(0, length(q)), Sigma)
    b <- q - mu
    a * sin(t) + b * cos(t) + mu
}

mu <- c(0.07447646, 0.9947984)
Sigma <- matrix(c(1.24064E-04, 0.0000562561, 0.0000562561, 0.0026316385), c(2, 2))
if (identity) {
    mass <- diag(2)
    epsilon <- 1/50
    L <- 75
    t <- pi/2
} else {
    mass <- solve(Sigma)
    epsilon <- 1/5000
    L <- 20
    t <- pi / 2
}
# U <- function(q) -logmvdnorm(q, mu, Sigma) # + 4404.4148
# grad_U <- function(q) -grad_logmvdnorm(q, mu, Sigma)
U <- function(q) lpr_mvn(q - mu, grad = F, solve(Sigma)) # + 4404.4148
grad_U <- function(q) attr(lpr_mvn(q - mu, grad = T, solve(Sigma)), "grad")
M <- 10000
q <- rep.int(1, 2)
accept <- 0
write(paste(c("state", "U", "x", "y"), collapse="\t"), stdout())
write(paste(c(0, -U(q), q), collapse="\t"), stdout())
for (i in 1:M) {
    if (exact) {
        q <- exact_HMC(t, mu, Sigma, q)
        accept <- accept + 1
    } else {
        h <- HMC(U, grad_U, epsilon, L, mass, q)
        q <- h$q
        if (h$a) accept <- accept + 1
    }
    write(paste(c(i, -U(q), q), collapse="\t"), stdout())
}
write(paste("accept: ", accept / M), stderr())
warnings()
