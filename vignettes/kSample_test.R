## ----srr-tags, eval = FALSE, echo = FALSE-------------------------------------
# #' srr tags
# #'
# #'
# #' @srrstats {G1.5} k-sample test example in the associated paper

## ----eval=FALSE---------------------------------------------------------------
# help(kb.test)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(mvtnorm)
library(QuadratiK)
library(ggplot2)
sizes <- rep(50,3)
eps <- 1
set.seed(2468)
x1 <- rmvnorm(sizes[1], mean = c(0,sqrt(3)*eps/3))
x2 <- rmvnorm(sizes[2], mean = c(-eps/2,-sqrt(3)*eps/6))
x3 <- rmvnorm(sizes[3], mean = c(eps/2,-sqrt(3)*eps/6))
x <- rbind(x1, x2, x3)
y <- as.factor(rep(c(1, 2, 3), times = sizes))

## ----fig.width=6, fig.height=4------------------------------------------------
ggplot(data.frame(x = x, y = y), aes(x = x[,1], y = x[,2], color = y)) +
  geom_point(size = 2) +
  labs(title = "Generated Points", x = "X1", y = "X2") +
  theme_minimal()

## ----fig.width=6, fig.height=4------------------------------------------------
set.seed(2468)
h_k <- select_h(x = x, y = y, alternative = "location")

## -----------------------------------------------------------------------------
h_k$h_sel

## -----------------------------------------------------------------------------
set.seed(2468)
k_test <- kb.test(x = x, y = y, h = h_k$h_sel)
show(k_test)

## -----------------------------------------------------------------------------
summary_ktest <- summary(k_test)
summary_ktest$summary_tables

## ----eval=FALSE---------------------------------------------------------------
# k_test_h <- kb.test(x = x, y = y)

## ----eval=FALSE---------------------------------------------------------------
# help(select_h)

