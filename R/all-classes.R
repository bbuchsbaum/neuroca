#' @rdname plsfit-class
#' @export
setClass("plsfit", 
         representation(d="numeric",brainSaliences="matrix",designSaliences="matrix",brainScores="matrix", designScores="matrix", nsigcomp="integer",
                        ZBootDesign="matrix", ZBootBrain="matrix", strata="factor"),
         prototype(d=numeric(), brainSaliences = matrix(), designSaliences=matrix(), brainScores=matrix(), designScores=matrix(), nsigcomp=as.integer(0), ZBootBrain=matrix(), ZBootDesign=matrix(), strata=factor())
)