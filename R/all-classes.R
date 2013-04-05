#' @rdname plsfit-class
#' @export
setClass("plsfit", 
         representation(Y="factor", X="matrix", d="numeric",brainSaliences="matrix",designSaliences="matrix",brainScores="matrix", designScores="matrix", nsigcomp="numeric",
                        Centroids="matrix", GlobalCentroid="numeric", svd="list",
                        ZBootDesign="matrix", ZBootBrain="matrix", strata="factor"),
         prototype(d=numeric(), brainSaliences = matrix(), designSaliences=matrix(), brainScores=matrix(), designScores=matrix(), nsigcomp=as.integer(0), ZBootBrain=matrix(), ZBootDesign=matrix(), strata=factor())
)


