#' neuroca - multivariate component analysis methods for neuroimaging data
#'
#' @author Bradley Buchsbaum <brad.buchsbaum@gmail.com>
#' @docType package
#' @import assertthat
#' @importFrom stats cor coef lsfit median model.matrix optim predict quantile resid residuals rnorm sd t.test terms var as.formula
#' @importFrom utils combn
#' @importFrom purrr map map_dbl
#' @name neuroca
#' @useDynLib neuroca
#' @importFrom Rcpp sourceCpp
NULL