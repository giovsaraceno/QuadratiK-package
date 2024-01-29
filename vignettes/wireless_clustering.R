## -----------------------------------------------------------------------------
library(QuadratiK)
head(wireless)

## -----------------------------------------------------------------------------
wire <- wireless[,-8]
labels <- wireless[,8]
wire_norm <- wire/sqrt(rowSums(wire^2))
set.seed(2468)
res_pk <- pkbc(as.matrix(wire_norm),3:6)

