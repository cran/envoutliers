#' Kernel regression smoothing - Only intended for developer use
#'
#' Nonparametric estimation of regression function using kernel regression with local or global data-adaptive plug-in bandwidth.
#' The function is called by \code{\link{KRDetect.outliers.changepoint}}, \code{\link{KRDetect.outliers.controlchart}} and \code{\link{KRDetect.outliers.EV}}. This function is not intended for use by regular users of the package.
#'
#' @param x a numeric vector of design points.
#' @param y a numeric vector of data values of the same length as \code{x}.
#' @param bandwidth.type a character string specifying the type of bandwidth, must be \code{"local"} or \code{"global"}.
#' @param bandwidth.value a local bandwidth array (for \code{bandwidth.type = "local"}) or global bandwidth value (for \code{bandwidth.type = "global"}) for kernel regression estimation. If \code{bandwidth.type = "NULL"} (default) a data-adaptive local plug-in (Herrmann, 1997) (for \code{bandwidth.type = "local"}) or data-adaptive global plug-in (Gasser et al., 1991) (for \code{bandwidth.type = "global"}) bandwidth is used instead.
#' @details This function computes the estimate of kernel regression function using a local or global data-adaptive plug-in algorithm.
#' The function is exported for developer use only. It does not perform any checks on inputs since it is only convenience function for computing estimates of the kernel regression function using functions implemented in package \pkg{lokern}.
#' @return A numeric vector of estimates of the kernel regression function (smoothed data).
#' @references Gasser T, Kneip A, Kohler W (1991). A flexible and fast method for automatic smoothing. Journal of the American Statistical Association, 86, 643–652.
#'
#' Herrmann E (1997). Local bandwidth choice in kernel regression estimation. Journal of Computational and Graphical Statistics, 6(1), 35–54.
#'
#' Eva Herrmann; Packaged for R and enhanced by Martin Maechler (2016). lokern: Kernel Regression Smoothing with Local or Global Plug-in Bandwidth. R package version 1.1-8. https://CRAN.R-project.org/package=lokern
#' @importFrom lokern glkerns lokerns
smoothing <- function(x,
                      y,
                      bandwidth.type,
                      bandwidth.value) {

  if (bandwidth.type == "local" & is.null(bandwidth.value)) {
    data.smoothed = lokerns(x = x, y = y, x.out = x, deriv = 0)$est
  }
  if (bandwidth.type == "local" & !is.null(bandwidth.value)) {
    data.smoothed = lokerns(x = x, y = y, x.out = x, deriv = 0, bandwidth = bandwidth.value)$est
  }
  if (bandwidth.type == "global" & is.null(bandwidth.value)) {
    data.smoothed = glkerns(x = x, y = y, x.out = x, deriv = 0)$est
  }
  if (bandwidth.type == "global" & !is.null(bandwidth.value)) {
    data.smoothed = glkerns(x = x, y = y, x.out = x, deriv = 0, bandwidth = bandwidth.value)$est
  }

  return(data.smoothed)
}

#' Changepoint outlier detection plot - Only intended for developer use
#'
#' Plot of results obtained using function \code{\link{KRDetect.outliers.changepoint}} for identification of outliers using changepoint analysis.
#' The function is called by \code{\link{KRDetect.outliers.plot}} and is not intended for use by regular users of the package.
#'
#' @param x a list obtained as an output of function \code{\link{KRDetect.outliers.changepoint}} for identification of outliers using changepoint analysis.
#' @param segments a logical variable specifying if vertical lines representing individual segments are plotted.
#' @details This function plots the results obtained using function \code{\link{KRDetect.outliers.changepoint}} identificating outliers using changepoint analysis based method.
#' The function is exported for developer use only. It does not perform any checks on inputs since it is only convenience function for plotting results obtained using functions implemented in package \pkg{envoutliers}.
#' @importFrom graphics plot points segments
KRDetect.changepoint.plot <- function(x,
                                      segments) {
  plot(x$index, x$x, cex = 0.8, cex.axis = 0.7, cex.lab = 0.7, xlab = "x", ylab = "data values", mgp=c(2, 1, 0))
  points(x$index, x$smoothed, type = "l", col = "darkgrey")
  points(x$index[x$outlier], x$x[x$outlier], col = "red", cex = 0.8, pch = 16)

  if (segments) {
    for (i in 1:max(x$changepoints, na.rm = TRUE)) {
      segments(min(x$index[x$changepoints == i], na.rm = TRUE), min(x$x, na.rm = TRUE),  min(x$index[x$changepoints == i], na.rm = TRUE), max(x$x, na.rm = TRUE), lty = 2 )
    }
  }
}

#' Control chart outliers detection plot - Only intended for developer use
#'
#' Plot of results obtained using function \code{\link{KRDetect.outliers.controlchart}} for identification of outliers using control charts.
#' The function is called by \code{\link{KRDetect.outliers.plot}} and is not intended for use by regular users of the package.
#'
#' @param x a list obtained as an output of function \code{\link{KRDetect.outliers.controlchart}} for identification of outliers using control charts.
#' @param all a logical variable specifying if individual graphs for outliers detected using control chart x, R and s are plotted together with graph visualising outliers detected based on at least 1 control chart.
#' If \code{all = FALSE}, only one graph visualising outliers detected based on at least 1 control chart is plotted.
#' @details This function plots the results obtained using function \code{\link{KRDetect.outliers.controlchart}} identificating outliers using control charts.
#' The function is exported for developer use only. It does not perform any checks on inputs since it is only convenience function for plotting results obtained using functions implemented in package \pkg{envoutliers}.
#' @importFrom graphics close.screen par plot points screen split.screen
KRDetect.controlchart.plot <- function(x,
                                       all) {
  par.old = par(no.readonly = TRUE)
  on.exit(par(par.old))

  if (all) {
    split.screen(figs = c(2,2) )
    screen(1)
    par(mar = c(1.5,1.5,1,0))  #DLHP

    plot(x$index, x$x, cex = 0.6, cex.axis = 0.7, cex.main = 0.7, cex.lab = 0.7, main = "Outliers detected using control chart x", xlab = "x", ylab = "data values", mgp=c(1, 0.5, 0))
    points(x$index, x$smoothed, type = "l", col = "darkgrey")
    points(x$index[x$outlier.x], x$x[x$outlier.x], col = "red", cex = 0.8, pch = 16)

    screen(2)
    par(mar = c(1.5,1.5,1,0))  #DLHP

    plot(x$index, x$x, cex = 0.6, cex.axis = 0.7, cex.main = 0.7, cex.lab = 0.7, main = "Outliers detected using control chart R", xlab = "x", ylab = "data values", mgp=c(1, 0.5, 0))
    points(x$index, x$smoothed, type = "l", col = "darkgrey")
    points(x$index[x$outlier.R], x$x[x$outlier.R], col = "red", cex = 0.8, pch = 16)

    screen(3)
    par(mar = c(1.5,1.5,1,0))  #DLHP

    plot(x$index, x$x, cex = 0.6, cex.axis = 0.7, cex.main = 0.7, cex.lab = 0.7, main = "Outliers detected using control chart s", xlab = "x", ylab = "data values", mgp=c(1, 0.5, 0))
    points(x$index, x$smoothed, type = "l", col = "darkgrey")
    points(x$index[x$outlier.s], x$x[x$outlier.s], col = "red", cex = 0.8, pch = 16)

    screen(4)
    par(mar = c(1.5,1.5,1,0))  #DLHP

    plot(x$index, x$x, cex = 0.6, cex.axis = 0.7, cex.main = 0.7, cex.lab = 0.7, main = "Outliers - all", xlab = "x", ylab = "data values", mgp=c(1, 0.5, 0))
    points(x$index, x$smoothed, type = "l", col = "darkgrey")
    points(x$index[x$outlier], x$x[x$outlier], col = "red", cex = 0.8, pch = 16)

    close.screen(all.screens = TRUE)
  } else {
    split.screen(figs = c(1,1) )
    screen(1)
    plot(x$index, x$x, cex = 0.8, cex.axis = 0.7, cex.main = 0.7, cex.lab = 0.7, main = "Outliers - all", xlab = "x", ylab = "data values", mgp=c(2, 1, 0))
    points(x$index, x$smoothed, type = "l", col = "darkgrey")
    points(x$index[x$outlier], x$x[x$outlier], col = "red", cex = 0.8, pch = 16)
  }


}

#' Extreme value outlier detection plot - Only intended for developer use
#'
#' Plot of results obtained using function \code{\link{KRDetect.outliers.EV}} for identification of outliers using extreme value theory.
#' The function is called by \code{\link{KRDetect.outliers.plot}} and is not intended for use by regular users of the package.
#'
#' @param x a list obtained as an output of function \code{\link{KRDetect.outliers.EV}} for identification of outliers using extreme value theory.
#' @param all a logical variable specifying if individual graphs for outliers with extremely low and extremely high value are plotted together with graph visualising outliers with both extremely low and extremely high value.
#' If \code{all = FALSE}, only one graph visualising outliers with both extremely low and extremely high value is plotted.
#' @details This function plots the results obtained using function \code{\link{KRDetect.outliers.EV}} identificating outliers using changepoint analysis based method.
#' The function is exported for developer use only. It does not perform any checks on inputs since it is only convenience function for plotting results obtained using functions implemented in package \pkg{envoutliers}.
#' @importFrom graphics par plot points screen split.screen
KRDetect.EV.plot <- function(x, all) {
  par.old = par(no.readonly = TRUE)
  on.exit(par(par.old))

  if (all) {
    split.screen(figs = c(2,1) )
    screen(1)
    split.screen(figs = c(1,2) )
    screen(3)
    par(mar = c(1.5,1.5,1,0))
    plot(x$index, x$x, cex = 0.6, cex.axis = 0.7, cex.main = 0.7, cex.lab = 0.7, main = "Outliers max", xlab = "x", ylab = "data values", mgp=c(1, 0.5, 0))
    points(x$index, x$smoothed, type = "l", col = "darkgrey")
    points(x$index[x$outlier.max], x$x[x$outlier.max], col = "red", cex = 0.8, pch = 16)

    screen(4)
    par(mar = c(1.5,1.5,1,0))
    plot(x$index, x$x, cex = 0.6, cex.axis = 0.7, cex.main = 0.7, cex.lab = 0.7, main = "Outliers min", xlab = "x", ylab = "data values", mgp=c(1, 0.5, 0))
    points(x$index, x$smoothed, type = "l", col = "darkgrey")
    points(x$index[x$outlier.min], x$x[x$outlier.min], col = "red", cex = 0.8, pch = 16)

    screen(2)
    par(mar = c(1.5,1.5,1,0))
    plot(x$index, x$x, cex = 0.6, cex.axis = 0.7, cex.main = 0.7, cex.lab = 0.7, main = "Outliers - all", xlab = "x", ylab = "data values", mgp=c(1, 0.5, 0))
    points(x$index, x$smoothed, type = "l", col = "darkgrey")
    points(x$index[x$outlier], x$x[x$outlier], col = "red", cex = 0.8, pch = 16)
  } else {
    split.screen(figs = c(1,1) )
    screen(1)
    par(mar = c(3,2.5,1,0))
    plot(x$index, x$x, cex = 0.8, cex.axis = 0.7, cex.lab = 0.7, xlab = "x", ylab = "data values", mgp=c(2, 1, 0))
    points(x$index, x$smoothed, type = "l", col = "darkgrey")
    points(x$index[x$outlier], x$x[x$outlier], col = "red", cex = 0.8, pch = 16)
  }
}

#' Outlier detection plot
#'
#' Plot of results obtained using functions \code{\link{KRDetect.outliers.changepoint}}, \code{\link{KRDetect.outliers.controlchart}} and \code{\link{KRDetect.outliers.EV}} for identification of outliers.
#' The function graphically visualizes results obtained using functions for outlier detection implemented in package \pkg{envoutliers}.
#'
#' @param x a list obtained as an output of function \code{\link{KRDetect.outliers.changepoint}}, \code{\link{KRDetect.outliers.controlchart}} or \code{\link{KRDetect.outliers.EV}} for identification of outliers.
#' @param segments a logical variable specifying if vertical lines representing individual segments are plotted. Only required for results obtained using \code{\link{KRDetect.outliers.changepoint}} function. Default is \code{segments = TRUE}.
#' @param all a logical variable.
#' Only required for results obtained using functions \code{\link{KRDetect.outliers.controlchart}} and \code{\link{KRDetect.outliers.EV}}. Default is \code{all = TRUE}.
#' In case of results obtained using function \code{\link{KRDetect.outliers.controlchart}} specifying if individual graphs for outliers detected using control chart x, R and s are plotted together with graph visualising outliers detected based on at least 1 control chart.
#' If \code{all = FALSE}, only one graph visualising outliers detected based on at least 1 control chart is plotted.
#' In case of results obtained using function \code{\link{KRDetect.outliers.EV}} specifying if individual graphs for outliers with extremely low and extremely high value are plotted together with graph visualising outliers with both extremely low and extremely high value.
#' If \code{all = FALSE}, only one graph visualising outliers with both extremely low and extremely high value is plotted.
#' @details This function plots the results obtained using function \code{\link{KRDetect.outliers.changepoint}}, \code{\link{KRDetect.outliers.controlchart}} or \code{\link{KRDetect.outliers.EV}} implemented in package \pkg{envoutliers} and identificating outliers.
#' @examples data("mydata", package = "openair")
#' x = mydata$o3[format(mydata$date, "%m %Y") == "12 2002"]
#' result = KRDetect.outliers.EV(x)
#' KRDetect.outliers.plot(result)
#' @export
KRDetect.outliers.plot <- function(x,
                                   segments = TRUE,
                                   all = TRUE) {
  if (x$method.type == "changepoint") {
    KRDetect.changepoint.plot(x,segments)
  } else if (x$method.type == "controlchart") {
    KRDetect.controlchart.plot(x, all)
  } else if (x$method.type == "EV") {
    KRDetect.EV.plot(x, all)
  }
}
