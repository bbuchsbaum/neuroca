

do_null_sim <- function(type, nreps=1000, nboot=1000, nvars=200) {
  res <- replicate(nreps, {
    print("go")
    X <- matrix(rnorm(100*nvars), 100, nvars)
    Y <- factor(rep(letters[1:4], length.out=100))
    S <- factor(rep(1:10, each=10))
  
    bres <- bada(Y,X,S, center=TRUE)
    boot <- bootstrap(bres, nboot=nboot, type=type)
      list(z1=boot$zboot_scores)
  }, simplify=FALSE)
  
  z1 <- do.call(rbind, lapply(res, "[[", "z1"))
  data.frame(type=type, FPR1=sum(abs(z1[,1]) > 1.96)/nrow(z1), FPR2=sum(abs(z1[,2]) > 1.96)/nrow(z1), FPR3=sum(abs(z1[,3]) > 1.96)/nrow(z1))
}
  
sim1 <- do_null_sim(type="projection")
sim2 <- do_null_sim(type="unrotated")
sim3 <- do_null_sim(type="split_half")





