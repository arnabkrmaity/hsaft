#' Function to implement the horseshoe shrinkage prior in Bayesian survival regression
#'
#'
#' This function employs the algorithm provided by van der Pas et. al. (2016) for
#' log normal Accelerated Failure Rate (AFT) model to fit survival regression. The censored observations are updated
#' according to the data augmentation of approach of Tanner and Wong (1984).
#'
#'  The model is:
#'  \eqn{t_i} is response,
#'  \eqn{c_i} is censored time,
#'  \eqn{t_i^* = \min_(t_i, c_i)} is observed time,
#'  \eqn{w_i} is censored data, so \eqn{w_i = \log t_i^*} if \eqn{t_i} is event time and
#'  \eqn{w_i = \log t_i^*} if \eqn{t_i} is right censored
#'  \eqn{\log t_i=X\beta+\epsilon, \epsilon \sim N(0,\sigma^2)}
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
#' \item{LikelihoodSamples}{Posterior Samples of likelihood.}
#'
#'
#'
#' @examples \dontrun{
#' burnin <- 500   # number of burnin
#' nmc    <- 1000  # number of Markov Chain samples
#' y.sd   <- 1     # standard deviation of the data
#' p      <- 80    # number of covariates
#' n      <- 40    # number of samples
#' beta   <- as.vector(smoothmest::rdoublex(p))  # from double exponential distribution
#' x      <- mvtnorm::rmvnorm(n, mean = rep(0, p))  # from multivariate normal distribution
#' y.mu   <- x %*% beta  # mean of the data
#' y      <- as.numeric(stats::rnorm(n, mean = y.mu, sd = y.sd))  # from normal distribution
#' T      <- exp(y)   # AFT model
#' C      <- rgamma(n, shape = 1.75, scale = 3)  # censoring time
#' time   <- pmin(T, C)  # observed time is min of censored and true
#' status = time == T   # set to 1 if event is observed
#' ct     <- as.matrix(cbind(time = time, status = status))  # censored time
#'
#' posterior.fit <- hsaft(ct, x, method.tau = "truncatedCauchy", method.sigma = "Jeffreys",
#'                        burn = burnin, nmc = nmc)
#' summary(posterior.fit$BetaHat)
#'}
#'
#' @export

# 20 November 2016
# We modify this code to update \beta for blockwise for each tumor
# calculate likelihood for lpml
# compute predictive log(survival time)
# Sample beta using C++


hsaft <- function(ct, X, method.tau = c("fixed", "truncatedCauchy","halfCauchy"), tau = 1,
                  method.sigma = c("fixed", "Jeffreys"), Sigma2 = 1,
                  burn = 1000, nmc = 5000, thin = 1, alpha = 0.05)
{

  method.tau = match.arg(method.tau)

  method.sigma = match.arg(method.sigma)

  ptm=proc.time()
  N=burn+nmc
  effsamp=(N-burn)/thin
  n=nrow(X)
  p=ncol(X)

  time         <- ct[, 1]
  status       <- ct[, 2]
  censored.id  <- which(status == 0)
  n.censored   <- length(censored.id)  # number of censored observations
  X.censored   <- X[censored.id, ]
  y <- logtime <- log(time)   # for coding convenience, since the whole code is written with y


  ## parameters ##
  Beta=rep(0,p); lambda=rep(1,p);
  sigma_sq = Sigma2;

  ## output ##
  betaout       <- matrix(0, p, effsamp)
  lambdaout     <- matrix(0, p, effsamp)
  tauout        <- rep(0, effsamp)
  sigmaSqout    <- rep(1, effsamp)
  likelihoodout <- matrix(0, n, effsamp)
  predsurvout   <- matrix(0, n, effsamp)
  logtimeout    <- matrix(0, n, effsamp)




  ## which algo to use ##
  if(p>n)
  {
    algo=1
  } else {
    algo=2
  }

  ## matrices ##
  I_n=diag(n)
  l0=rep(0,p)
  l1=rep(1,n)
  l2=rep(1,p)
  if(algo==2)
  {
    Q_star=t(X)%*%X
  }


  ## start Gibb's sampling ##
  for(i in 1:N)
  {
    mean.impute <- X.censored %*% Beta
    sd.impute   <- sqrt(sigma_sq)
    ## update censored data ##
    time.censored <- msm::rtnorm(n.censored, mean = mean.impute, sd = sd.impute, lower = logtime[censored.id])
    # truncated at log(time) for censored data
    y[censored.id] <- time.censored

    mean  <- X %*% Beta
    sd    <- sqrt(sigma_sq)
    predictive.survivor <- stats::pnorm(mean/sd, lower.tail = FALSE)


    ## update beta ##
    if(algo==1)
    {
      lambda_star=tau*lambda
      U=as.numeric(lambda_star^2)*t(X)
      ## step 1 ##
      u=stats::rnorm(l2,l0,lambda_star)
      v=X%*%u + stats::rnorm(n)
      ## step 2 ##
      v_star=solve((X%*%U+I_n),((y/sqrt(sigma_sq))-v))
      Beta=sqrt(sigma_sq)*(u+U%*%v_star)
    }
    else if(algo==2)
    {
      lambda_star=tau*lambda
      L=chol((1/sigma_sq)*(Q_star+diag(1/as.numeric(lambda_star^2),p,p)))
      v=solve(t(L),t(t(y)%*%X)/sigma_sq)
      mu=solve(L,v)
      u=solve(L,stats::rnorm(p))
      Beta=mu+u
    }

    ## update lambda_j's in a block using slice sampling ##
    eta = 1/(lambda^2)
    upsi = stats::runif(p,0,1/(1+eta))
    tempps = Beta^2/(2*sigma_sq*tau^2)
    ub = (1-upsi)/upsi
    # now sample eta from exp(tempv) truncated between 0 & upsi/(1-upsi)
    Fub = 1 - exp(-tempps*ub) # exp cdf at ub
    Fub[Fub < (1e-4)] = 1e-4;  # for numerical stability
    up = stats::runif(p,0,Fub)
    eta = -log(1-up)/tempps
    lambda = 1/sqrt(eta);

    ## update tau ##
    ## Only if prior on tau is used
    if(method.tau == "halfCauchy"){
      tempt = sum((Beta/lambda)^2)/(2*sigma_sq)
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
      tempt = sum((Beta/lambda)^2)/(2*sigma_sq)
      et = 1/tau^2
      utau = stats::runif(1,0,1/(1+et))
      ubt_1=1
      ubt_2 = min((1-utau)/utau,p^2)
      Fubt_1 = stats::pgamma(ubt_1,(p+1)/2,scale=1/tempt)
      Fubt_2 = stats::pgamma(ubt_2,(p+1)/2,scale=1/tempt)
      #Fubt = max(Fubt,1e-8) # for numerical stability
      ut = stats::runif(1,Fubt_1,Fubt_2)
      et = stats::qgamma(ut,(p+1)/2,scale=1/tempt)
      tau = 1/sqrt(et)
    }

    ## update sigma_sq ##
    if(method.sigma == "Jeffreys"){
      if(algo==1)
      {
        E_1=max(t(y-X%*%Beta)%*%(y-X%*%Beta),(1e-10))
        E_2=max(sum(Beta^2/((tau*lambda))^2),(1e-10))
      } else {
        E_1=max(t(y-X%*%Beta)%*%(y-X%*%Beta),1e-8)
        E_2=max(sum(Beta^2/((tau*lambda))^2),1e-8)
      }
      sigma_sq= 1/stats::rgamma(1, (n + p)/2, scale = 2/(E_1+E_2))
    }

    # likelihood
    likelihood <- stats::dnorm(y, mean = X %*% Beta, sd = sqrt(sigma_sq))
    logt       <- X %*% Beta

    if (i%%500 == 0)
    {
      print(i)
    }



    if(i > burn && i%%thin== 0)
    {
      betaout[ ,(i-burn)/thin]         <- Beta
      lambdaout[ ,(i-burn)/thin]       <- lambda
      tauout[(i - burn)/thin]          <- tau
      sigmaSqout[(i - burn)/thin]      <- sigma_sq
      likelihoodout[ ,(i - burn)/thin] <- likelihood
      predsurvout[ ,(i - burn)/thin]   <- predictive.survivor
      logtimeout[, (i - burn)/thin]    <- logt
    }
  }


  pMean=apply(betaout,1,mean)
  pMedian=apply(betaout,1,stats::median)
  pLambda=apply(lambdaout,1,mean)
  pSigma=mean(sigmaSqout)
  pTau=mean(tauout)
  pPS <- apply(predsurvout, 1, mean)
  pLogtime<- apply(logtimeout, 1, mean)


  #construct credible sets
  left <- floor(alpha*effsamp/2)
  right <- ceiling((1-alpha/2)*effsamp)

  BetaSort <- apply(betaout, 1, sort, decreasing = F)
  left.points <- BetaSort[left, ]
  right.points <- BetaSort[right, ]

  result=list("SurvivalHat" = pPS, "LogTimeHat" = pLogtime, "BetaHat"=pMean, "LeftCI" = left.points,
              "RightCI" = right.points,"BetaMedian"=pMedian, "LambdaHat" = pLambda,
              "Sigma2Hat"=pSigma,"TauHat"=pTau,"BetaSamples"=betaout,
              "TauSamples" = tauout, "Sigma2Samples" = sigmaSqout, "LikelihoodSamples" = likelihoodout)
  return(result)
}

