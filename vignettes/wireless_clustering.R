## ----srr-tags, eval = FALSE, echo = FALSE-------------------------------------
#  #' srr tags
#  #'
#  #'
#  #' @srrstats {G1.5} application to the wireless data set in the
#  #'                  associated paper

## ----message=FALSE------------------------------------------------------------
library(QuadratiK)
head(wireless)

## -----------------------------------------------------------------------------
wire <- wireless[,-8]
labels <- wireless[,8]
wire_norm <- wire/sqrt(rowSums(wire^2))
set.seed(2468)
res_pk <- pkbc(as.matrix(wire_norm),3:5)

## -----------------------------------------------------------------------------
set.seed(2468)
res_validation <- pkbc_validation(res_pk, true_label = labels)
res_validation$IGP
round(res_validation$metrics, 5)

## ----eval=FALSE---------------------------------------------------------------
#  help(pkbc_validation)

## ----fig.show = "hold", fig.width=6, fig.height=6-----------------------------
plot(res_pk)

## ----fig.show = "hold", fig.width=8, fig.height=5-----------------------------
plot(res_pk, k=4, true_label = labels)

## -----------------------------------------------------------------------------
summary_clust <- stats_clusters(res_pk,4)

## -----------------------------------------------------------------------------
summary_clust

