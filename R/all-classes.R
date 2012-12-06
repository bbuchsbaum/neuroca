#' @rdname plsfit-class
#' @export
setClass("plsfit", 
         representation(d="numeric",brainSaliences="matrix",designSaliences="matrix",brainScores="matrix", designScores="matrix", nsigcomp="integer",
                        ZDesign="matrix", ZBrain="matrix")
         prototype(d=numeric(), brainSaliences = matrix(), designSaliences=matrix(), brainScores=matrix(), designScores=matrix(), nsigcomp=NA, ZBrain=matrix(), ZDesign=matrix())
)