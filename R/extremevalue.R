#' Return level estimation - Only intended for developer use
#'
#' Estimation of return level for a given threshold value using Peaks Over Threshold model.
#' The function is called by \code{\link{KRDetect.outliers.EV}} and is not intended for use by regular users of the package.
#'
#' @param r a return period.
#' @param u a threshold value.
#' @param sigma_u a scale parameter of Generalized Pareto distribution.
#' @param xi a shape parameter of Generalized Pareto distribution.
#' @param lambda_u a relative frequency of the number of threshold value exceedances.
#' @param theta_gomes an extremal index.
#' @details This function computes the estimate of return level for a given threshold value using Peaks Over Threshold model.
#' The function is exported for developer use only. It does not perform any checks on inputs since it is only convenience function used within \code{\link{KRDetect.outliers.EV}}.
#' @return A numeric value of return lever corresponding to return period \code{r}
#' @references Coles S (2001). An Introduction to Statistical Modeling of Extreme Values. 3 edition. London: Springer. ISBN 1-85233-459-2.
#'
#' Pickands J (1975). Statistical inference using extreme order statistics. The Annals of Statistics, 3(1), 119-131.
return.level.est <- function(r,
                             u,
                             sigma_u,
                             xi,
                             lambda_u,
                             theta_gomes) {
  z_r = u + (sigma_u / xi) * ((lambda_u ^ (-1) * (1 - (1 - r ^ (-1)) ^ (theta_gomes ^ (-1)))) ^ (-xi) - 1)

  return(z_r)
}

#' Extremal index estimation - Only intended for developer use
#'
#' Estimation of an extremal index using the block maxima approach suggested by (Gomes, 1993).
#' The function is called by \code{\link{KRDetect.outliers.EV}} and is not intended for use by regular users of the package.
#'
#' @param x a numeric vector of observations.
#' @param block.length a numeric value giving the length of blocks.
#' @details This function computes the estimate of extremal index suggested by (Gomes, 1993).
#' The function is exported for developer use only. It does not perform any checks on inputs since it is only convenience function used within \code{\link{KRDetect.outliers.EV}}.
#' @return A numeric value of an extremal index
#' @references Gomes M (1993). On the estimation of parameter of rare events in environmental time series. In Statistics for the Environment, volume 2 of Water Related Issues, pp. 225-241. Wiley.
#'
#' Heffernan JE, Stephenson AG (2016). ismev: An Introduction to Statistical Modeling of Extreme Values. R package version 1.41, URL http://CRAN.R-project.org/package=ismev.
#' @importFrom ismev gev.fit
extremal.index.gomes <- function(x,
                                 block.length) {
  blocks = seq(1, length(x), block.length)
  k = length(blocks)

  #-------------------------------
  # stationary data
  #-------------------------------
  max.stat = c(rep(0, k))
  for (i in 1:k) {
    if (i < k) {
      max.stat[i] = max(x[c(blocks[i]:(blocks[(i + 1)] - 1))])
    } else {
      max.stat[i] = max(x[c(blocks[i]:length(x))])
    }
  }

  #-------------------------------
  # permuted data
  #-------------------------------
  max.perm = c(rep(0, k))
  set.seed(1)
  x.permut = sample(x)

  for (i in 1:k) {
    if (i <  k) {
      max.perm[i] = max(x.permut[c(blocks[i]:(blocks[(i + 1)] - 1))])
    } else {
      max.perm[i] = max(x.permut[c(blocks[i]:length(x))])
    }
  }

  #---------------------------------------------------------------
  # gev.fit for stationary data
  #---------------------------------------------------------------
  gev.model = gev.fit(max.stat, show = FALSE)
  stat.sigma = gev.model$mle[2]
  stat.mu = gev.model$mle[1]

  #---------------------------------------------------------------
  # gev.fit for permuted data
  #---------------------------------------------------------------
  gev.model = gev.fit(max.perm, show = FALSE)
  perm.sigma = gev.model$mle[2]
  perm.mu = gev.model$mle[1]

  #-------------------------------
  # extremal index estimation
  #-------------------------------
  xi = (perm.sigma - stat.sigma) / (perm.mu - stat.mu)
  extremal.index = min(1, (perm.sigma / stat.sigma) ^ (-1 / xi))

  return(extremal.index)
}

#' Identification of outliers using extreme value theory
#'
#' Identification of outliers in environmental data using semiparametric method based on kernel smoothing and extreme value theory (Holesovsky et al., 2018). The outliers are identified as observations whose values are exceeded on average once a given period that is specified by the user.
#'
#' @param x a numeric vector of observations.
#' @param perform.smoothing a logical value specifying if data smoothing is performed. If \code{TRUE} (default), data are smoothed.
#' @param bandwidth.type a character string specifying the type of bandwidth, must be \code{"local"} (default) or \code{"global"}.
#' @param bandwidth.value a local bandwidth array (for \code{bandwidth.type = "local"}) or global bandwidth value (for \code{bandwidth.type = "global"}) for kernel regression estimation. If \code{bandwidth.type = "NULL"} (default) a data-adaptive local plug-in (Herrmann, 1997) (for \code{bandwidth.type = "local"}) or data-adaptive global plug-in (Gasser et al., 1991) (for \code{bandwidth.type = "global"}) bandwidth is used instead.
#' @param extremal.index.min a numeric value giving the extremal index for identification of outliers with extremely low value. If \code{extremal.index.min = NULL} (default), estimate of (Gomes, 1993) is used.
#' @param extremal.index.max a numeric value giving the extremal index for identification of outliers with extremely high value. If \code{extremal.index.max = NULL} (default), estimate of (Gomes, 1993) is used.
#' @param block.length a numeric value giving the length of blocks for estimation of extremal index. Default is \eqn{round(sqrt(length(na.omit(x))))}.
#' @param threshold.min a threshold value for residuals with low values, that is used to estimate the parameters of Generalized Pareto distribution. If \code{threshold.min = NULL} (default), threshold is estimated as 90\% quantile of smoothing residuals.
#' @param threshold.max a threshold value for residuals with high values, that is used to estimate the parameters of Generalized Pareto distribution. If \code{threshold.max = NULL} (default), threshold is estimated as 90\% quantile of smoothing residuals.
#' @param return.period a positive numeric value giving return period. Default is \code{r = 120}, which means that observations whose values are exceeded on average once every 120 observations are detected as outliers.
#' @details This function identifies outliers in time series using two-step procedure (Holesovsky et al., 2018). The procedure consists of kernel smoothing and extreme value estimation of high threshold exceedances for smoothing residuals.
#' Outliers with both extremely high and extremely low values are identified.
#' Crucial for the method is the choice of return period - parameter defining the criterion for outliers detection.
#' The outliers with extremely high values are detected as observations whose values are exceeded on average once a given return.period of observations. Analogous, the outliers with extremely low values are identified.
#' @return A list is returned with elements:
#' \item{method.type}{a character string giving the type of method used for outlier idetification}
#' \item{x}{a numeric vector of observations}
#' \item{index}{a numeric vector of index design points assigned to individual observations}
#' \item{smoothed}{a numeric vector of estimates of the kernel regression function (smoothed data)}
#' \item{sigma_u.min}{a numeric value giving scale parameter of Generalised Pareto distribution used for identification of outliers with extremely low value}
#' \item{sigma_u.max}{a numeric value giving scale parameter of Generalised Pareto distribution used for identification of outliers with extremely high value}
#' \item{xi.min}{a numeric value giving shape parameter of Generalised Pareto distribution used for identification of outliers with extremely low value}
#' \item{xi.max}{a numeric value giving shape parameter of Generalised Pareto distribution used for identification of outliers with extremely high value}
#' \item{lambda_u.min}{a numeric value giving relative frequency of the number of threshold value exceedances and identification of outliers with extremely low value}
#' \item{lambda_u.max}{a numeric value giving relative frequency of the number of threshold value exceedances and identification of outliers with extremely high value}
#' \item{extremal.index.min}{a numeric value giving extremal index used for identification of outliers with extremely low value}
#' \item{extremal.index.max}{a numeric value giving extremal index used for identification of outliers with extremely high value}
#' \item{threshold.min}{a numeric value giving threshold used in Peaks Over Threshold model and identification of outlier with extremely low value}
#' \item{threshold.max}{a numeric value giving threshold used in Peaks Over Threshold model and identification of outlier with extremely high value}
#' \item{return.level.min}{a numeric value giving return level used for identification of outliers with extremely low value}
#' \item{return.level.max}{a numeric value giving return level used for identification of outliers with extremely high value}
#' \item{outlier.min}{a logical vector specyfing the identified outliers with extremely low value. \code{TRUE} means that corresponding observation from vector \code{x} is detected as outlier}
#' \item{outlier.max}{a logical vector specyfing the identified outliers with extremely high value. \code{TRUE} means that corresponding observation from vector \code{x} is detected as outlier}
#' \item{outlier}{a logical vector specyfing the identified outliers with both extremely low and extremely high value. \code{TRUE} means that corresponding observation from vector \code{x} is detected as outlier}
#' @references Holesovsky J, Campulova M, Michalek J (2018). Semiparametric Outlier Detection in Nonstationary Times Series: Case Study for Atmospheric Pollution in Brno, Czech Republic. Atmospheric Pollution Research, 9(1).
#'
#' Theo Gasser, Alois Kneip & Walter Koehler (1991) A flexible and fast method for automatic smoothing. Journal of the American Statistical Association 86, 643-652. https://doi.org/10.2307/2290393
#'
#' E. Herrmann (1997) Local bandwidth choice in kernel regression estimation. Journal of Graphical and Computational Statistics 6, 35-54.
#'
#' Herrmann E, Maechler M (2013). lokern: Kernel Regression Smoothing with Local or Global Plug-in Bandwidth. R package version 1.1-5, URL http://CRAN.R-project.org/package=lokern.
#'
#' Gomes M (1993). On the estimation of parameter of rare events in environmental time series. In Statistics for the Environment, volume 2 of Water Related Issues, pp. 225-241. Wiley.
#'
#' Heffernan JE, Stephenson AG (2016). ismev: An Introduction to Statistical Modeling of Extreme Values. R package version 1.41, URL http://CRAN.R-project.org/package=ismev.
#'
#' Coles S (2001). An Introduction to Statistical Modeling of Extreme Values. 3 edition. London: Springer. ISBN 1-85233-459-2.
#'
#' Pickands J (1975). Statistical inference using extreme order statistics. The Annals of Statistics, 3(1), 119-131.
#' @examples data("mydata", package = "openair")
#' x = mydata$o3[format(mydata$date, "%m %Y") == "12 2002"]
#' result = KRDetect.outliers.EV(x)
#' KRDetect.outliers.plot(result)
#' @importFrom ismev gpd.fit
#' @importFrom stats na.omit quantile
#' @export
KRDetect.outliers.EV <- function(x,
                                 perform.smoothing = TRUE,
                                 bandwidth.type = "local",
                                 bandwidth.value = NULL,
                                 extremal.index.min = NULL,
                                 extremal.index.max = NULL,
                                 block.length = round(sqrt(length(na.omit(x)))),
                                 threshold.min = NULL,
                                 threshold.max = NULL,
                                 return.period = 120) {
  if (!bandwidth.type %in% c("local", "global")) {
    stop("Invalid bandwidth type")
  }

  data.input = data.frame(x)
  data.input$index = c(1:dim(data.input)[1])
  data = na.omit(data.input)

  if (dim(data)[1] < 4) {
    stop("Not enough data")
  }

  data$smoothed = NA

  if (perform.smoothing) {
    data$smoothed = smoothing(data$index, data$x, bandwidth.type, bandwidth.value)
    data$residuals = data$x - data$smoothed
  } else {
    data$residuals = data$x
  }

  #-------------------------------
  # extremal index estimation
  #-------------------------------
  if (is.null(extremal.index.min)) {
    extremal.index.min = extremal.index.gomes(x = (-data$residuals), block.length = block.length)
  }

  if (is.null(extremal.index.max)) {
    extremal.index.max = extremal.index.gomes(x = data$residuals, block.length = block.length)
  }

  #-------------------------------
  # threshld value computing
  #-------------------------------
  if (is.null(threshold.min)) {
    threshold.min = quantile(-data$residuals, 0.9)
  }

  if (is.null(threshold.max)) {
    threshold.max = quantile(data$residuals, 0.9)
  }

  #------------------------------------------
  # Generalised Pareto distribution fit
  #------------------------------------------
  gpd.model.min = gpd.fit(-data$residuals, threshold.min, show = FALSE)
  sigma_u.min = gpd.model.min$mle[1]
  xi.min = gpd.model.min$mle[2]
  lambda_u.min = gpd.model.min$nexc/length(data$residuals)

  gpd.model.max = gpd.fit(data$residuals, threshold.max, show = FALSE)
  sigma_u.max = gpd.model.max$mle[1]
  xi.max = gpd.model.max$mle[2]
  lambda_u.max = gpd.model.max$nexc/length(data$residuals)

  #-------------------------------
  # return lever estimation
  #-------------------------------
  return.level.min = return.level.est(r = return.period, u = threshold.min, sigma_u = sigma_u.min, xi = xi.min, lambda_u = lambda_u.min, theta_gomes = extremal.index.min)
  return.level.max = return.level.est(r = return.period, u = threshold.max, sigma_u = sigma_u.max, xi = xi.max, lambda_u = lambda_u.max, theta_gomes = extremal.index.max)

  #-------------------------------
  # outlier identification
  #-------------------------------
  data$outliers.min = FALSE
  data$outliers.min[-data$residuals > return.level.min] = TRUE

  data$outliers.max = FALSE
  data$outliers.max[data$residuals > return.level.max] = TRUE

  data$outliers = FALSE
  data$outliers[(-data$residuals > return.level.min) | (data$residuals > return.level.max)] = TRUE

  data.output = merge(data.input, data, by = c("index", "x"), all.x = TRUE)

  result = list(method.type = "EV",
                x = data.output$x,
                index = data.output$index,
                smoothed = data.output$smoothed,
                sigma_u.min = sigma_u.min,
                sigma_u.max = sigma_u.max,
                xi.min = xi.min,
                xi.max = xi.max,
                lambda_u.min = lambda_u.min,
                lambda_u.max = lambda_u.max,
                extremal.index.min = extremal.index.min,
                extremal.index.max = extremal.index.max,
                threshold.min = threshold.min,
                threshold.max = threshold.max,
                return.level.min = return.level.min,
                return.level.max = return.level.max,
                outlier.min = data.output$outliers.min,
                outlier.max = data.output$outliers.max,
                outlier = data.output$outliers)

  return(result)
}
