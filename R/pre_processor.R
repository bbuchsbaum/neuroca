# 
# 
# compose_scaler <- function(X, ...) {
#   
# }
# 
# dual_scaler <- function(X, center_col=TRUE, scale_col=FALSE, center_row=FALSE, scale_row=FALSE) {
#   cs <- col_scaler(X, center=center_col, scale=scale_col)
#   rs <- row_scaler(pre_process(cs, X), center=center_row, scale=scale_row)
#   
#   
#   f <- function(newX) {
#     
#     print("dual_scaler: pp")
#     pre_process(rs, pre_process(cs, newX))
#   }
#   
#   structure(list(pre_process=f), class=c("dual_scaler", "pre_processor"))
#   
# }
# 
# 
# row_scaler <- function(X, center=TRUE, scale=TRUE, lazy=FALSE) {
#   
#   center_vec <- if (center) {
#     rowMeans(X)
#   }
#   
#   scale_vec <- if (scale) {
#     matrixStats::rowSds(X)
#   }
#   
#   f <- if (center & scale) {
#     function(newX) {
#       print("row scale")
#       apply(newX, 2, function(x) (x - center_vec)/scale_vec)
#     }
#   } else if (center & !scale) {
#     function(newX) apply(newX, 2, function(x) (x - center_vec))
#   } else if (scale && ! center) {
#     function(newX) apply(newX, 2, function(x) (x/scale_vec))
#   } else {
#     function(newX) newX
#   }
#     
#   #Xp <- apply(X, 1, function(x) x - meanvec)
#   structure(list(center_vec=center_vec, scale_vec=scale_vec, pre_process=f), class=c("row_scaler", "pre_processor"))
# }
# 
# col_scaler <- function(X, center=TRUE, scale=TRUE) {
#   
#   center_vec <- if (center) {
#     colMeans(X)
#   }
#   
#   scale_vec <- if (scale) {
#     matrixStats::colSds(X)
#   }
#   
#   f <- if (center & scale) {
#     function(newX) {
#       print("col scale")
#       t(apply(newX, 1, function(x) (x - center_vec)/scale_vec))
#     }
#   } else if (center & !scale) {
#     function(newX) t(apply(newX, 1, function(x) (x - center_vec)))
#   } else if (scale && ! center) {
#     function(newX) t(apply(newX, 1, function(x) (x/scale_vec)))
#   } else {
#     function(newX) newX
#   }
#   
#   #Xp <- apply(X, 1, function(x) x - meanvec)
#   structure(list(center_vec=center_vec, scale_vec=scale_vec, pre_process=f), class=c("col_scaler", "pre_processor"))
# }
# 
# 
# 
# 
# pre_process.row_scaler <- function(obj, Xnew, ...) { obj$pre_process(Xnew) }
# 
# pre_process.col_scaler <- function(obj, Xnew, ...) { obj$pre_process(Xnew) }
# 
# pre_process.dual_scaler <- function(obj, Xnew, ...) { obj$pre_process(Xnew) }
# 
# 

