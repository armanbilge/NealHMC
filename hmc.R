#!/usr/bin/env Rscript

library("MASS")

logmvdnorm <- function (x, mu, Sigma)
{
    d <- x - mu
    return -length(x)/2 * log(2 * pi) - log(det(Sigma))/2 - d * solve(Sigma) * d / 2
}

HMC <- function (U, grad_U, epsilon, L, mass, current_q)
{
    q <- current_q
    p <- mvrnorm(length(q), 0, mass) # multivariate normal variates
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
    if (log(runif(1)) < current_U-proposed_U+current_K-proposed_K)
    {
        return (q)  # accept
    }
    else
    {
        return (current_q)  # reject
    }
}
