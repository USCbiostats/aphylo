context("Formulae")

test_that("Formulas create the right model", {
  
  set.seed(0772)
  x <- rdrop_annotations(raphylo(50), .5)

  suppressWarnings({  
  m_mu <- aphylo_formula(x ~ mu_d)
  m_mu_psi0 <- aphylo_formula(x ~ mu_d + psi)
  m_mu_psi1 <- aphylo_formula(x ~ psi)
  m_mu_psi2 <- aphylo_formula(x ~ psi + mu_d)
  })
  
  expect_equal(deparse(m_mu_psi0)[-1], deparse(m_mu_psi1)[-1])
  expect_equal(deparse(m_mu_psi0)[-1], deparse(m_mu_psi2)[-1])
  
})

test_that("Errors are caught", {
  
  set.seed(071)
  x <- rdrop_annotations(raphylo(50), .5)
  
  expect_error(suppressWarnings(aphylo_formula(y ~ psi)), "be found")
  y <- 1
  expect_error(suppressWarnings(aphylo_formula(y ~ psi)), "should be either")
  
  expect_error(aphylo_formula(x ~ mu_d, c(mu_d0=1, mu_d1=1, psi1=0)), "overspecified")
  expect_error(aphylo_formula(x ~ mu_d, c(mu_d0=1)), "missing")
  
  expect_error(suppressWarnings(aphylo_formula(x~mu_d(9))), "Arguments passed to")
  expect_error(suppressWarnings(aphylo_formula(x~mu_d + I(x))), "supported")
  
  
})

test_that("Multiple chains formulas work", {
  
  set.seed(121)
  
  x <- raphylo(30)
  p <- matrix(runif(7*2), nrow=2, dimnames = NULL)
  
  expect_warning(
    ans <- aphylo_mcmc(
      x~psi+mu_d+eta+Pi, params = p,
      control = list(
        nchains=2, nsteps=500, burnin=10,
        conv_checker = NULL)),
    "matched by position"
  )
  
  expect_s3_class(ans, "aphylo_estimates")
  
})

