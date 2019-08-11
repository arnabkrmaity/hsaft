#'
#'
#' This function extends the main function \code{\link{hsaft}} to create correlation among groups.
#'
#'
#'
#' @references Stephanie van der Pas, James Scott, Antik Chakraborty and Anirban Bhattacharya (2016). horseshoe:
#' Implementation of the Horseshoe Prior. R package version 0.1.0.
#' https://CRAN.R-project.org/package=horseshoe
#'
#' Arnab Kumar Maity, Anirban Bhattacharya, Bani K. Mallick, and Veerabhadran Baladandayuthapani (2017).
#' Joint Bayesian Estimation and Variable Selection for TCPA Protein Expression Data
#'
#'
#'
#'@param ct Response, a \eqn{n*2} matrix with first column as response and second column as right censored indicator,
#'1 is event time and 0 is right censored.
#'@param X Matrix of covariates, dimension \eqn{n*p}.
#'@param method.tau Method for handling \eqn{\tau}. Select "truncatedCauchy" for full
#' Bayes with the Cauchy prior truncated to [1/p, 1], "halfCauchy" for full Bayes with
#' the half-Cauchy prior, or "fixed" to use a fixed value (an empirical Bayes estimate,
#' for example).
#'@param tau  Use this argument to pass the (estimated) value of \eqn{\tau} in case "fixed"
#' is selected for method.tau. Not necessary when method.tau is equal to"halfCauchy" or
#' "truncatedCauchy". The default (tau = 1) is not suitable for most purposes and should be replaced.
#'@param method.sigma Select "Jeffreys" for full Bayes with Jeffrey's prior on the error
#'variance \eqn{\sigma^2}, or "fixed" to use a fixed value (an empirical Bayes
#'estimate, for example).
#'@param Sigma2 A fixed value for the error variance \eqn{\sigma^2}. Not necessary
#'when method.sigma is equal to "Jeffreys". Use this argument to pass the (estimated)
#'value of Sigma2 in case "fixed" is selected for method.sigma. The default (Sigma2 = 1)
#'is not suitable for most purposes and should be replaced.
#'@param burn Number of burn-in MCMC samples. Default is 1000.
#'@param nmc Number of posterior draws to be saved. Default is 5000.
#'@param thin Thinning parameter of the chain. Default is 1 (no thinning).
#'@param alpha Level for the credible intervals. For example, alpha = 0.05 results in
#'95\% credible intervals.
#'@param r number of groups.
#'@param n.seq a vector of sample sizes for all groups.
#'@param pk number of covariates in each group.
#'
#'
#'
#'@return \item{SurvivalHat}{Predictive survival probability.}
#'\item{LogTimeHat}{Predictive log time.}
#'\item{BetaHat}{Posterior mean of Beta, a \eqn{p} by 1 vector.}
#' \item{LeftCI}{The left bounds of the credible intervals.}
#' \item{RightCI}{The right bounds of the credible intervals.}
#' \item{BetaMedian}{Posterior median of Beta, a \eqn{p} by 1 vector.}
#' \item{Sigma2Hat}{Posterior mean of error variance \eqn{\sigma^2}. If method.sigma =
#' "fixed" is used, this value will be equal to the user-selected value of Sigma2
#' passed to the function.}
#' \item{TauHat}{Posterior mean of global scale parameter tau, a positive scalar.
#' If method.tau = "fixed" is used, this value will be equal to the user-selected value
#' of tau passed to the function.}
#' \item{BetaSamples}{Posterior samples of Beta.}
#' \item{TauSamples}{Posterior samples of tau.}
#' \item{Sigma2Samples}{Posterior samples of Sigma2.}
#' \item{BHat}{Posterior samples of b which is the mean of \eqn{\beta}.}
#' \item{LikelihoodSamples}{Posterior Samples of likelihood.}
#'
#'
#'
#' @examples
#' \dontrun{# Examples for hsaftgroupcorr function
#' burnin <- 500   # number of burnin
#' nmc    <- 1000  # number of Markov Chain samples
#' y.sd   <- 1     # standard deviation of the data
#' p      <- 80    # number of covariates
#' r      <- 5     # number of groups
#' p      <- 80    # number of covariate in each group
#' n1     <- 40    # sample size of 1st group
#' n2     <- 50    # sample size of 2nd group
#' n3     <- 70    # sample size of 3rd group
#' n4     <- 100   # sample size of 4th group
#' n5     <- 120   # sample size of 5th group
#' n      <- sum(c(n1, n2, n3, n4, n5))  # total sample size
#' n.seq  <- c(n1, n2, n3, n4, n5)
#' Beta   <- matrix(smoothmest::rdoublex(p * r), nrow = r, ncol = p, byrow = TRUE)
#' # from double exponential distribution
#' beta   <- as.vector(t(Beta))  # vectorize Beta
#' x1     <- mvtnorm::rmvnorm(n1, mean = rep(0, p))
#' x2     <- mvtnorm::rmvnorm(n2, mean = rep(0, p))
#' x3     <- mvtnorm::rmvnorm(n3, mean = rep(0, p))
#' x4     <- mvtnorm::rmvnorm(n4, mean = rep(0, p))
#' x5     <- mvtnorm::rmvnorm(n5, mean = rep(0, p))  # from multivariate normal distribution
#' y.mu1  <- x1 %*% Beta[1, ]
#' y.mu2  <- x2 %*% Beta[2, ]
#' y.mu3  <- x3 %*% Beta[3, ]
#' y.mu4  <- x4 %*% Beta[4, ]
#' y.mu5  <- x5 %*% Beta[5, ]
#' y1     <- stats::rnorm(n1, mean = y.mu1, sd = y.sd)
#' y2     <- stats::rnorm(n2, mean = y.mu2, sd = y.sd)
#' y3     <- stats::rnorm(n3, mean = y.mu3, sd = y.sd)
#' y4     <- stats::rnorm(n4, mean = y.mu4, sd = y.sd)
#' y5     <- stats::rnorm(n5, mean = y.mu5, sd = y.sd)
#' y      <- c(y1, y2, y3, y4, y5)
#' x      <- Matrix::bdiag(x1, x2, x3, x4, x5)
#' X      <- as.matrix(x)
#' y      <- as.numeric(as.matrix(y))  # from normal distribution
#' T      <- exp(y)   # AFT model
#' C      <- rgamma(n, shape = 1.75, scale = 3)  # censoring time
#' time   <- pmin(T, C)  # observed time is min of censored and true
#' status = time == T   # set to 1 if event is observed
#' ct     <- as.matrix(cbind(time = time, status = status))  # censored time

#' posterior.fit <- hsaftgroupcorr(ct, X, method.tau = "truncatedCauchy", method.sigma = "Jeffreys",
#'                                 burn = burnin, nmc = nmc,
#'                                 r = r, n.seq = n.seq, pk = p)
#' summary(posterior.fit$BetaHat)
#'}
#'
#' @export

# 20 November 2016
# We modify this code to update \beta for blockwise for each tumor
# Introduction of Correlation among groups
# calculate likelihood for lpml
# No c
# inclusion of n_k
# Inclusion of f(n_k)
# compute predictive log(survival time)


hsaftgroupcorr <- function(ct, X, method.tau = c("fixed", "truncatedCauchy","halfCauchy"), tau = 1,
                           method.sigma = c("fixed", "Jeffreys"), Sigma2 = 1,
                           burn = 1000, nmc = 5000, thin = 1, alpha = 0.05,
                           r, n.seq, pk)
{

  method.tau = match.arg(method.tau)

  method.sigma = match.arg(method.sigma)

  ptm=proc.time()
  niter =burn+nmc
  effsamp=(niter -burn)/thin
  n=nrow(X)
  p=ncol(X)

  time        <- ct[, 1]
  status      <- ct[, 2]
  censored.id <- which(status == 0)
  n.censored  <- length(censored.id)  # number of censored observations
  X.censored   <- X[censored.id, ]
  y <- logtime <- log(time)   # for coding convenience, since the whole code is written with y

  ## parameters ##
  beta           <- rep(0, p)
  Beta           <- matrix(0, nrow = pk, ncol = r);
  lambda         <- rep(1,p);
  Lambda         <- matrix(1, nrow = pk, ncol = r)
  sigma_sq       <- Sigma2;
  b              <- rep(0, pk)
  B              <- rep(0, p)
  sigma.b.square <- rep(1, pk)
  # f.n.seq        <- ifelse(n.seq < mean(n.seq), 1/n.seq, n.seq/sum(n.seq))
  f.n.seq        <- n.seq/n.seq
  betan          <- rep(0, p)

  ## output ##
  betaout       <- matrix(0, p, effsamp)
  lambdaout     <- matrix(0, p, effsamp)
  tauout        <- rep(0, effsamp)
  sigmaSqout    <- rep(1, effsamp)
  bout          <- matrix(0, pk, effsamp)
  likelihoodout <- matrix(0, n, effsamp)
  predsurvout   <- matrix(0, n, effsamp)
  logtimeout    <- matrix(0, n, effsamp)



  ## start Gibb's sampling ##
  for(i in 1:niter)
  {

    mean.impute <- X.censored %*% beta
    sd.impute   <- sqrt(sigma_sq)
    ## update censored data ##
    time.censored <- msm::rtnorm(n.censored, mean = mean.impute, sd = sd.impute, lower = logtime[censored.id])
    # truncated at log(time) for censored data
    y[censored.id] <- time.censored

    mean  <- X %*% beta
    sd    <- sqrt(sigma_sq)
    predictive.survivor <- stats::pnorm(mean/sd, lower.tail = FALSE)

    ## update beta ##
    sum.nk <- 0
    sum.pk <- 0

    sum.nk1 <- sum.nk + 1
    sum.nk  <- sum.nk + n.seq[1]
    sum.pk1 <- sum.pk + 1
    sum.pk  <- sum.pk + pk

    lambdak <- lambda[sum.pk1:sum.pk]
    Xk <- X[sum.nk1:sum.nk, sum.pk1:sum.pk]
    yk <- y[sum.nk1:sum.nk]


    Dk     <- tau^2 * (1/f.n.seq[1]) * diag(lambdak^2)
    D.invk <- chol2inv(chol(Dk))
    D.inv  <- D.invk
    Ak     <- t(Xk) %*% Xk + D.invk
    A.invk <- chol2inv(chol(Ak))

    beta.mean  <- A.invk %*% (t(Xk) %*% yk + D.invk %*% B[sum.pk1:sum.pk])
    beta.Sigma <- sigma_sq * A.invk
    betak      <- t(mvtnorm::rmvnorm(1, mean = beta.mean, sigma = beta.Sigma))

    betan[sum.pk1:sum.pk] <- betak * sqrt(f.n.seq[1])
    beta[sum.pk1:sum.pk]  <- betak

    for(k in 2:r)
    {
      sum.nk1 <- sum.nk + 1
      sum.nk  <- sum.nk + n.seq[k]
      sum.pk1 <- sum.pk + 1
      sum.pk  <- sum.pk + pk

      lambdak <- lambda[sum.pk1:sum.pk]
      Xk <- X[sum.nk1:sum.nk, sum.pk1:sum.pk]
      yk <- y[sum.nk1:sum.nk]

      Dk     <- tau^2 * (1/f.n.seq[k]) * diag(lambdak^2)
      D.invk <- chol2inv(chol(Dk))
      D.inv  <- Matrix::bdiag(D.inv, D.invk)
      Ak     <- t(Xk) %*% Xk + D.invk
      A.invk <- chol2inv(chol(Ak))

      beta.mean  <- A.invk %*% (t(Xk) %*% yk + D.invk %*% B[sum.pk1:sum.pk])
      beta.Sigma <- sigma_sq * A.invk
      betak      <- t(mvtnorm::rmvnorm(1, mean = beta.mean, sigma = beta.Sigma))

      betan[sum.pk1:sum.pk] <- betak * sqrt(f.n.seq[k])
      beta[sum.pk1:sum.pk]  <- betak
    }


    Beta  <- matrix(beta, nrow = pk, ncol = r, byrow = FALSE)
    D.inv <- as.matrix(D.inv)


    ## update lambda_j's in a block using slice sampling ##
    eta = 1/(lambda^2)
    upsi = stats::runif(p,0,1/(1+eta))
    tempps = betan^2/(2*sigma_sq*tau^2)
    ub = (1-upsi)/upsi
    # now sample eta from exp(tempv) truncated between 0 & upsi/(1-upsi)
    Fub = 1 - exp(-tempps*ub) # exp cdf at ub
    Fub[Fub < (1e-4)] = 1e-4;  # for numerical stability
    up = stats::runif(p,0,Fub)
    eta = -log(1-up)/tempps
    lambda = 1/sqrt(eta);
    Lambda <- matrix(lambda, nrow = pk, ncol = r, byrow = FALSE)

    ## update tau ##
    ## Only if prior on tau is used
    if(method.tau == "halfCauchy"){
      tempt = sum((betan/lambda)^2)/(2*sigma_sq)
      et = 1/tau^2
      utau = stats::runif(1,0,1/(1+et))
      ubt = (1-utau)/utau
      Fubt = stats::pgamma(ubt,(p+1)/2,scale=1/tempt)
      Fubt = max(Fubt,1e-8) # for numerical stability
      ut = stats::runif(1,0,Fubt)
      et = stats::qgamma(ut,(p+1)/2,scale=1/tempt)
      tau = 1/sqrt(et)
    }#end if

    if(method.tau == "truncatedCauchy"){
      tempt = sum((betan/lambda)^2)/(2*sigma_sq)
      et = 1/tau^2
      utau = stats::runif(1,0,1/(1+et))
      ubt_1=1
      ubt_2 = min((1-utau)/utau,p^2)
      Fubt_1 = stats::pgamma(ubt_1,(p+1)/2,scale=1/tempt)
      Fubt_2 = stats::pgamma(ubt_2,(p+1)/2,scale=1/tempt)
      Fubt_2 = max(Fubt_2, 1e-8) # for numerical stability
      ut = stats::runif(1,Fubt_1,Fubt_2)
      et = stats::qgamma(ut,(p+1)/2,scale=1/tempt)
      tau = 1/sqrt(et)
    }

    ## update sigma_sq ##
    if(method.sigma == "Jeffreys"){

      E_1  <- max(t(y-X%*%beta)%*%(y-X%*%beta),(1e-10))
      E_21 <- max(sum((beta)^2/((tau*lambda))^2),(1e-10))
      E_22 <- max(t(beta - B) %*% D.inv %*% (beta - B), 1e-10)
      E_2  <- min(c(E_21, E_22))  # for numerical stability

      sigma_sq= 1/stats::rgamma(1, (n + p)/2, scale = 2/(E_1+E_2))
    }


    ## update bP for proteins ##
    for(j in 1:pk)
    {
      denominator <- sigma_sq * tau^2
      tau.b       <- (1/denominator) * sum(f.n.seq/Lambda[j, ]^2) + (1/sigma.b.square[j])
      tau.b.inv   <- 1/tau.b
      b.mean      <- (tau.b.inv/denominator) * sum((Beta[j, ])/(Lambda[j, ]^2/f.n.seq))
      b[j]       <- stats::rnorm(1, mean = b.mean, sd = sqrt(tau.b.inv))
    }
    B <- rep(b, k)

    # likelihood
    likelihood <- stats::dnorm(y, mean = X %*% beta, sd = sqrt(sigma_sq))
    logt       <- X %*% beta


    if (i%%500 == 0)
    {
      print(i)
    }



    if(i > burn && i%%thin== 0)
    {
      betaout[, (i - burn)/thin]       <- beta
      lambdaout[, (i - burn)/thin]     <- lambda
      tauout[(i - burn)/thin]          <- tau
      sigmaSqout[(i-burn)/thin]        <- sigma_sq
      bout[, (i - burn)/thin]          <- b
      likelihoodout[ ,(i - burn)/thin] <- likelihood
      predsurvout[ ,(i - burn)/thin]   <- predictive.survivor
      logtimeout[, (i - burn)/thin]    <- logt
    }
  }


  pMean   <- apply(betaout, 1, mean)
  pMedian <- apply(betaout, 1, stats::median)
  pLambda <- apply(lambdaout, 1, mean)
  pSigma  <- mean(sigmaSqout)
  pTau    <- mean(tauout)
  pB      <- apply(bout, 1, mean)
  pPS     <- apply(predsurvout, 1, mean)
  pLogtime<- apply(logtimeout, 1, mean)


  #construct credible sets
  left  <- floor(alpha*effsamp/2)
  right <- ceiling((1-alpha/2)*effsamp)

  BetaSort     <- apply(betaout, 1, sort, decreasing = F)
  left.points  <- BetaSort[left, ]
  right.points <- BetaSort[right, ]

  result=list("SurvivalHat" = pPS, "LogTimeHat" = pLogtime, "BetaHat"=pMean, "LeftCI" = left.points,
              "RightCI" = right.points,"BetaMedian"=pMedian, "LambdaHat" = pLambda,
              "Sigma2Hat"=pSigma,"TauHat"=pTau,"BetaSamples"=betaout,
              "TauSamples" = tauout, "Sigma2Samples" = sigmaSqout, "BHat" = pB,
              "LikelihoodSamples" = likelihoodout)
  return(result)
}

