context("Formulae")

test_that("Formulas create the right model", {
  
  set.seed(0772)
  x <- rdrop_annotations(sim_annotated_tree(100), .5)
  
  m_mu <- aphylo_formula(x ~ mu)
  m_mu_psi0 <- aphylo_formula(x ~ mu + psi)
  m_mu_psi1 <- aphylo_formula(x ~ psi)
  m_mu_psi2 <- aphylo_formula(x ~ psi + mu)
  
  expect_equal(deparse(m_mu_psi0)[-1], deparse(m_mu_psi1)[-1])
  expect_equal(deparse(m_mu_psi0)[-1], deparse(m_mu_psi2)[-1])
  
})

test_that("Errors are caught", {
  
  set.seed(071)
  x <- rdrop_annotations(sim_annotated_tree(100), .5)
  
  expect_error(aphylo_formula(y ~ psi), "be found")
  y <- 1
  expect_error(aphylo_formula(y ~ psi), "should be an")
  
  expect_error(aphylo_formula(x ~ mu, c(mu0=1, mu1=1, psi1=0)), "overspecified")
  expect_error(aphylo_formula(x ~ mu, c(mu0=1)), "missing")
  
  expect_error(aphylo_formula(x~mu(9)), "Arguments passed to")
  expect_error(aphylo_formula(x~mu + I(x)), "supported")
  
  
})