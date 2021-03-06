#' Outlier detection plot
#'
#' This function is deprecated. Use \code{\link{plot.KRDetect}} instead.
#'
#' Plot of results obrained using functions KRDetect.outliers.changepoint, KRDetect.outliers.controlchart and KRDetect.outliers.EV for identification of outliers.
#' The function graphically visualizes results obtained using functions for outlier detection implemented in package envoutliers.
#'
#' @param x a list obtained as the output of function KRDetect.outliers.changepoint, KRDetect.outliers.controlchart or KRDetect.outliers.EV for identification of outliers.
#' @param all a logical variable,
#' in case of results obtained using function KRDetect.outliers.controlchart specifying if individual graphs for outliers detected using control chart x, R and s are plotted together with graph visualising outliers detected based on at least 1 control chart.
#' If all = FALSE, only one graph visualising outliers detected based on at least 1 control chart is plotted.
#' in case of results obtained using function KRDetect.outliers.EV specifying if individual graphs for outliers with extremely low and extremely high value are plotted together with graph visualising outliers with both extremely low and extremely high value.
#' If all = FALSE, only one graph visualising outliers with both extremely low and extremely high value is plotted.
#' Only required for results obtained using functions KRDetect.outliers.controlchart and KRDetect.outliers.EV. Default is all = TRUE.
#' @param segments a logical variable specifying if vertical lines representing individual segments are plotted. Only required for results obtained using KRDetect.outliers.changepoint function. Default is segments = TRUE.
#' @details This function plots the results obtained using function KRDetect.outliers.changepoint, KRDetect.outliers.controlchart or KRDetect.outliers.EV implemented in package envoutliers and identificating outliers.
#' @export
KRDetect.outliers.plot <- function(x, all, segments) {
  .Deprecated("plot.KRDetect")
}
