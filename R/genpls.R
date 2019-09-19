# genpls <- function(X, Y, A=NULL, M=NULL, ncomp=min(dim(X)), 
#                    preproc=center(), deflation=FALSE, stripped=FALSE) {
#   
#   Y <- as.matrix(Y)
#   pcon <- prep_constraints(X, A, M)
#   A <- pcon$A
#   M <- pcon$M
#   
#   nobj <- dim(X)[1] # n in paper
#   npred <- dim(X)[2] # p in paper
#   nresp <- dim(Y)[2]
#   
#   V <- R <- matrix(0, nrow = npred, ncol = ncomp)
#   tQ <- matrix(0, nrow = ncomp, ncol = nresp) # Y loadings; transposed
#   B <- array(0, dim = c(npred, nresp, ncomp))
#   
#   if (!stripped) {
#     P <- R
#     U <- TT <- matrix(0, nrow = nobj, ncol = ncomp)
#     fitted <- array(0, dim = c(nobj, nresp, ncomp))
#   }
#   
#   
#   procres <- prep(preproc, X)
#   Xp <- procres$Xp
#   
#   S <- crossprod(X,Y)
#   #CP_i <- CP
#   
#   
#   
#   for (i in 1:ncomp) {
#     if (ncol(S) > 3) {
#       sres <- RSpectra::svds(S,1)
#       u_i <- sres$v
#       v_i <- sres$u
#     } else {
#       u_i <- matrix(rnorm(nrow(Y)), nrow(Y), 1)
#       v_i <- matrix(rnorm(ncol(X)), ncol(x), 1)
#     }
#     
#     conv <- 1e-5
#     v_old <- v_i
#     u_old <- u_i
#     
#     z <- matrix(0, nrow(Y), ncomp)
#     while (TRUE) {
#       u_i <- crossprod(CP_i, v_i) 
#       u_i <- u_i / sqrt(sum(u_i^2))
#       
#       v_i <- CP_i %*% u_i
#       v_i <- v_i / sqrt(sum(v_i^2))
#       
#       d <- abs(u_i - u_old) + abs(v_i - v_old)
#       print(d)
#       if (d < conv) {
#         break
#       }
#     }
# 
#     z[,i] <- X %*% v_i
#     
#     CP_i <- P %*% CP_i
#     R_i <- 
#     
#     
#     
#   }
#   
# }