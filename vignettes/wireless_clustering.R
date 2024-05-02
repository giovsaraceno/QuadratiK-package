## ----srr-tags, eval = FALSE, echo = FALSE-------------------------------------
#  #' srr tags
#  #'
#  #'
#  #' @srrstats {G1.5} application to the wireless data set in the
#  #'                  associated paper

## -----------------------------------------------------------------------------
library(QuadratiK)
head(wireless)

## -----------------------------------------------------------------------------
wire <- wireless[,-8]
labels <- wireless[,8]
wire_norm <- wire/sqrt(rowSums(wire^2))
set.seed(2468)
res_pk <- pkbc(as.matrix(wire_norm),3:6)

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(2468)
#  res_validation <- validation(res_pk, true_label = labels)
#  res_validation$IGP
#  round(res_validation$metrics, 8)

## ----fig.width=6, fig.height=8------------------------------------------------
summary_clust <- extract_stats(res_pk,4)

## -----------------------------------------------------------------------------
summary_clust

