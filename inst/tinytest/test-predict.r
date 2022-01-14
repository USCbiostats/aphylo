# context("Prediction functions")

# test_that("Prediction works", {
  
# Jmorr tree
tree  <- matrix(c(1, 2, 1, 3, 3, 4, 3, 5), ncol=2, byrow = TRUE)
tree  <- as.phylo(tree)
X     <- c(0, 0, 0)
atree <- new_aphylo(tip.annotation = X, tree = tree)

# Model parameters
Pi  <- .3
psi <- c(.01, .02)
mu  <- c(.2, .1)
eta <- c(1, 1)

# John's spreadsheet
ans0 <- c(
  0.009268032907,
  0.006404403815,
  0.006831065254,
  0.006041731783,
  0.006041731783
)[with(atree$tree, c(tip.label, node.label))]

# Computing brute-force and using the pre-order
# ans1 <- predict_brute_force(atree, psi, mu_d = mu, mu_s = mu, Pi)
# saveRDS(ans1, "inst/tinytest/test-predict-ans1.rds")
ans1 <- readRDS("test-predict-ans1.rds")
ans2 <- predict_pre_order(atree, psi, mu_d = mu, mu_s = mu, eta, Pi, loo = FALSE)


expect_equivalent(ans0, ans1$posterior)
expect_equivalent(ans1$posterior, ans2[,1])

# # Alternative take
# p_up <- matrix(ncol=2, nrow = length(ans0))
# p_pred <- double(nrow(p_up))
# l  <- LogLike(atree, psi = psi, mu_d = mu, mu_s = mu, eta = eta, Pi = Pi)
# for (i in rev(atree$pseq)) {
#   
#   # If the root note
#   if (Ntip(atree) == (i - 1)) {
#     
#     p_up[i, 1] = 1 - Pi
#     p_up[i, 2] = Pi
#     
#   }
#   
#   # Computing the offsprings
#   for (o in atree$offspring[[i]]) {
#     p_up[o, 1] = p_up[i, 1] * (1 - mu[1]) + p_up[i, 2] * mu[2]
#     p_up[o, 2] = p_up[i, 1] * mu[1]       + p_up[i, 2] * (1 - mu[2])
#   }
#  
#   p_pred[i] <-  p_up[i, 2] * l$Pr[[1]][i, 2]/
#     (p_up[i, 1] * l$Pr[[1]][i, 1] + p_up[i, 2] * l$Pr[[1]][i, 2])
# }
# 
# p_pred
# ans0

# Random test
set.seed(122331)
atree <- raphylo(4, psi = psi, mu_d = mu, mu_s = mu, eta = eta, Pi = Pi)

# ans0 <- predict_brute_force(atree, psi, mu_d = mu, mu_s = mu, Pi)
# saveRDS(ans0, file = "inst/tinytest/test-predict-ans0.rds")
ans0 <- readRDS("test-predict-ans0.rds")
ans1 <- predict_pre_order(atree, psi, mu_d = mu, mu_s = mu, eta, Pi, loo = FALSE)

expect_equivalent(ans0$posterior, ans1[,1])

# })

# test_that("Calling the prediction function works", {

set.seed(137245)

x <- raphylo(10)
x_obs <- rdrop_annotations(x, .5)
res   <- suppressWarnings({
  aphylo_mcmc(
    x_obs ~ psi + mu_d + Pi, priors = bprior(),
    control = list(nsteps = 2e3, burnin = 0)
    )
  })

ans0 <- predict_pre_order(
  x    = x_obs,
  psi  = res$par[c("psi0", "psi1")],
  mu_d = res$par[c("mu_d0", "mu_d1")],
  mu_s = res$par[c("mu_d0", "mu_d1")],
  eta  = c(1,1)/2, #res$par[c("eta0", "eta1")],
  Pi   = res$par["Pi"], loo = FALSE
  )

ans1 <- predict(res, ids = list(1:Nnode(x_obs, internal.only = FALSE)), loo = FALSE)

expect_silent(plot(prediction_score(res)))
expect_equivalent(ans0, ans1)
  
# Predictions on new data ------------------------------------------------------

set.seed(88)
res <- suppressMessages(
  suppressWarnings(
    aphylo_mcmc(
      x_obs ~ psi + mu_d + mu_s + eta + Pi, priors = bprior(),
      control = list(nsteps = 2e3, burnin = 0))
    )
  )
z <- raphylo(20, node.type = sample.int(2, 19, TRUE)-1)

ans0_a <- predict(res, loo = FALSE)
ans1_a <- predict_pre_order(
  res$dat,
  psi = res$par[c("psi0", "psi1")],
  eta = res$par[c("eta0", "eta1")],
  mu_s = res$par[c("mu_s0", "mu_s1")],
  mu_d = res$par[c("mu_d0", "mu_d1")],
  Pi  = res$par["Pi"], loo = FALSE
)

ans0_b <- predict(res, newdata = z, loo = FALSE)
ans1_b <- predict_pre_order(
  z,
  psi = res$par[c("psi0", "psi1")],
  eta = res$par[c("eta0", "eta1")],
  mu_s = res$par[c("mu_s0", "mu_s1")],
  mu_d = res$par[c("mu_d0", "mu_d1")],
  Pi  = res$par["Pi"], loo = FALSE
  )

expect_equivalent(ans0_a, ans1_a)
expect_equivalent(ans0_b, ans1_b)

expect_true(nrow(ans0_a) != nrow(ans0_b))

# Prediction score function ----------------------------------------------------

# test_that("Best vs Worse Prediction score", {
  
set.seed(123)
y <- sample(c(0,1), 20, replace = TRUE)

# Perfect prediction score
ans0 <- prediction_score(cbind(y), cbind(y))
expect_equivalent(ans0$obs[1], 1)

# Worse
ans1 <- prediction_score(cbind(y), 1 - cbind(y))
expect_equivalent(ans1$obs[1], 0)

# })

# test_that("Random prediction score", {

set.seed(123)
y <- cbind(sample(c(0,1), 20, replace = TRUE, prob = c(.8, .2)))
a <- .3

p0 <- mean(aphylo:::predict_random(1, y, diag(20), alpha0 = a, alpha1 = 1-a))
p1 <- aphylo:::prediction_score_rand(y, diag(20), alpha0 = a, alpha1 = 1-a)

expect_equivalent(p0, p1, tol = 1e-1)
  
# Mutliphylo -------------------------------------------------------------------
set.seed(192318)
net <- c(raphylo(20), raphylo(10))
ans <- aphylo_mcmc(net ~ mu_d + mu_s + Pi, control = list(nsteps = 2e3, burnin = 0))

pred <- prediction_score(ans)
expect_output(pred, "Prediction score")

# Prediction of multiple samples
set.seed(1123);p_all <- predict(ans, nsamples = 10)
set.seed(1123);p_1 <- predict(ans, which.tree = 1, nsamples = 10)
expect_equivalent(p_all[[1]], p_1[[1]])

# Predicting all, including interior nodes -------------------------------------
pred0 <- predict(ans, ids = lapply(Nnode(ans, internal.only = FALSE), seq_len))
pred1 <- predict(ans, ids = lapply(Nnode(ans, internal.only = FALSE), function(a) {
  seq_len(a - 5) + 5
  }))

expect_equivalent(rep(NA_real_, 10), unlist(lapply(pred1, "[", i = 1:5)))
expect_true(!any(is.na(unlist(pred0))))
expect_equal(
  lapply(pred0, "[", -c(1:5)),
  lapply(pred1, "[", -c(1:5))
)
