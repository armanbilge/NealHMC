#!/usr/bin/env Rscript

# set.seed(666)
treeModel <- T
identity <- F
exact <- F
heatmap <- T
java <- T
parallel <- T

library(MASS)

if (parallel) {
    library(foreach)
    library(doParallel)
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl, cores = detectCores() - 1)
}

if (java) {
    library(rJava)
    .jinit("beast.jar")
    Double <- J("java.lang.Double")
    ArrayList <- J("java.util.ArrayList")
    TraceCorrelation <- J("dr.inference.trace.TraceCorrelation")
    TraceType <- J("dr.inference.trace.TraceFactory")$TraceType
} else {
    library(coda)
}

lpr_mvn <- function (value, grad=FALSE, inv.cov)
{
  g <- - inv.cov %*% value
  r <- as.vector (t(value) %*% g / 2)
  if (grad)
  { attr(r,"grad") <- as.vector(g)
  }

  r
}

logmvdnorm <- function (x, mu, Sigma)
{
    d <- x - mu
    - length(x)/2 * log(2 * pi) - log(det(Sigma))/2 - t(d) %*% solve(Sigma) %*% d / 2
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
    # write(current_U-proposed_U+current_K-proposed_K, stderr())
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

sample <- function (U, grad_U, epsilon, L, mass, M) {
    q <- rep.int(1, dim(mass)[1])
    accept <- 0
    X = matrix(nrow = M+1, ncol = length(q)+2)
    X[1,] <- c(0, -U(q), q)
    for (i in 1:M) {
        if (exact) {
            q <- exact_HMC(t, mu, Sigma, q)
            accept <- accept + 1
        } else {
            h <- HMC(U, grad_U, epsilon, L, mass, q)
            q <- h$q
            if (h$a) accept <- accept + 1
        }
        X[i+1,] <- c(i, -U(q), q)
    }
    list(X = X, a = accept/M)
}

if (heatmap) {
    mu <- c(0.07447646, 0.9947984)
    Sigma <- matrix(c(1.24064E-04, 0.0000562561, 0.0000562561, 0.0026316385), c(2, 2))
    epsilon <- sqrt(min(eigen(Sigma)$values))
    L <- 4
    M <- 10000
    U <- function(q) -logmvdnorm(q, mu, Sigma)
    grad_U <- function(q) -grad_logmvdnorm(q, mu, Sigma)
    n <- 8
    m <- 8
    ESS <- matrix(ncol = n+m+1, nrow = n+m+1)
    for (i in -n:m) {
        data <- foreach (j = -n:m, .packages = c("MASS")) %dopar%  {
            sample(U, grad_U, epsilon, L, diag(c(2^i, 2^j)), M)$X
            # write.table(X, paste("x", i, j, ".log", sep = ""), quote = F, sep = "\t", col.names = c("state", "U", "x", "y"), row.names = F)
        }
        for (j in 1:length(data)) {
            X <- data[[j]]
            if (java) {
                al <- new(ArrayList, as.integer(M*9/10+1))
                for (d in X[(M/10+1):(M+1),2]) {
                    al$add(new(Double, d))
                }
                tr <- new(TraceCorrelation, al, TraceType$DOUBLE, as.integer(1))
                ess <- .jsimplify(tr$getESS())
            } else {
                mcmc <- as.mcmc(X[(M/10):(M+1),], start = M/10)
                ess <- effectiveSize(mcmc)[2]
            }
            ESS[i+n+1,j] <- ess
        }
    }
    print(ESS)
    heatmap(ESS, col = heat.colors(256), Rowv=NA, Colv=NA, labRow=as.character(-n:m), labCol=as.character(-n:m))
} else {
    t <- pi / 2
    if (treeModel) {
        mu <- c(0.07447646, 0.9947984)
        Sigma <- matrix(c(1.24064E-04, 0.0000562561, 0.0000562561, 0.0026316385), c(2, 2))
        if (identity) {
            mass <- diag(2)
            epsilon <- 1/50
            L <- 75
        } else {
            mass <- solve(Sigma)
            epsilon <- 1/10000
            L <- 8
        }
    } else {
        mu <- c(1, 1)
        Sigma <- matrix(c(0.01, 0, 0, 0.00001), c(2, 2))
        if (identity) {
            mass <- diag(2)
            epsilon <- 1/200
            L <- 40
        } else {
            mass <- Sigma
            epsilon <- 1/30000
            L <- 40
        }
    }

    U <- function(q) -logmvdnorm(q, mu, Sigma) # + 4404.4148
    grad_U <- function(q) -grad_logmvdnorm(q, mu, Sigma)
    # U <- function(q) - length(q)/2 * log(2 * pi) - log(det(Sigma))/2 - lpr_mvn(q - mu, grad = F, solve(Sigma)) # + 4404.4148
    # grad_U <- function(q) -attr(lpr_mvn(q - mu, grad = T, solve(Sigma)), "grad")
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
}

if (parallel) {
    stopCluster(cl)
}
