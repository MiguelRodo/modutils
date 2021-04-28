test_that("test_lr works", {
  set.seed(2106); n_id <- 30
  test_tbl <- tibble::tibble(id = rep(as.character(1:n_id), each = 2),
                             grp = rep(as.character(1:2), each = n_id),
                             y_re = rep(rnorm(n_id), each = 2),
                             x = rnorm(n_id * 2),
                             y = x^4 * ifelse(grp == "1", 3, 0) + y_re + rnorm(20, sd = 0.3))
  test_tbl$y <- test_tbl$y +  abs(min(test_tbl$y)) + 0.001
  test_tbl$y <- test_tbl$y/max(test_tbl$y + 0.001)

  # check that it runs for each model type
  # ------------------------

  # fit model
  mod_lm_fix <- lm(y ~ -1, data = test_tbl)
  mod_lm_int <- lm(y ~ 1, data = test_tbl)
  mod_lm <- lm(y ~ splines::ns(x, 3), data = test_tbl)
  mod_lmer <- lme4::lmer(y ~ x + (1|grp), data = test_tbl)
  mod_glmer <- suppressMessages(suppressWarnings(lme4::glmer(y ~ x + (1|grp), data = test_tbl,
                                                             family = 'binomial')))

  # check associated extractor functions
  # ------------------------

  # check that we extract log-likelihoods
  expect_true(is.numeric(get_ll(mod_lm)))
  expect_true(is.numeric(get_ll(mod_lmer)))
  expect_true(is.numeric(get_ll(mod_glmer)))

  # check that we extract estimates properly
  expect_true(tibble::is_tibble(get_coef(mod_lm)))
  expect_true(tibble::is_tibble(get_coef(mod_lm_fix)))
  expect_identical(ncol(get_coef(mod_lm_fix)), 5L)
  expect_identical(nrow(get_coef(mod_lm_fix)), 0L)
  expect_identical(nrow(get_coef(mod_lmer)), 2L)
  expect_identical(colnames(get_coef(mod_glmer)), c("var", "est", "se", "test_stat", "p"))

  # check that degrees of freedom are correct
  expect_identical(get_dof(mod_lm_int), 1L)
  expect_identical(get_dof(mod_lm), 4L)
  expect_identical(get_dof(mod_lmer), 2L)
  expect_identical(get_dof(mod_glmer), 2L)

  # check lr test output
  # ------------------------

  # fit simpler models
  mod_lm_lin <- lm(y ~ x, data = test_tbl)
  mod_lmer_int <- lme4::lmer(y ~ (1|grp), data = test_tbl)
  mod_glmer_int <- suppressMessages(suppressWarnings(lme4::glmer(y ~ (1|grp), data = test_tbl,
                                                             family = 'binomial')))
  # check
  test_lr(mod_glmer, mod_glmer_int)
  test_lr(mod_lmer, mod_lmer_int)
  expect_identical(test_lr(mod_lm, mod_lm_lin)$var,
                   "splines::ns(x, 3)1; splines::ns(x, 3)2; splines::ns(x, 3)3")
  expect_identical(test_lr(mod_lm, list("x" = mod_lm_lin))$var, "x")
  expect_identical(test_lr(mod_lm, mod_lm_lin, "x")$var, "x")

  # compare to wald output
  wald_tbl_lm_lin <- test_wald(mod_lm_lin, list("x"))
  lr_tbl_lm_lin <- test_lr(mod_lm_lin, mod_lm_int)
  expect_true(abs(wald_tbl_lm_lin$p - lr_tbl_lm_lin$p) < 0.06)

  set.seed(2106)
  test_tbl_norm <- tibble::tibble(x = rnorm(20, sd = 0.3))
  test_tbl_norm$y <- test_tbl_norm$x^3 + rnorm(20, sd = 0.2)
  test_tbl_norm$y_ind <- rnorm(20, sd = 50)
  test_tbl_norm$x_ind <- rnorm(20, sd = 30)

  # check where there is a relationship
  mod_lm_norm <- lm(y ~ splines::ns(x, 3), data = test_tbl_norm)
  mod_lm_norm_int <- lm(y ~ 1, data = test_tbl_norm)
  wald_tbl_lm_spline <- test_wald(mod_lm_norm, list("x"), match_condn = "any")
  lr_tbl_lm_spline <- test_lr(mod_lm_norm, mod_lm_norm_int, "x")
  expect_true(abs(wald_tbl_lm_spline$p - lr_tbl_lm_spline$p) < 0.05)

  # check where there's no relationship
  mod_lm_norm_ind <- lm(y_ind ~ splines::ns(x_ind, 3), data = test_tbl_norm)
  mod_lm_norm_int_ind <- lm(y_ind ~ 1, data = test_tbl_norm)
  wald_tbl_lm_spline <- test_wald(mod_lm_norm_ind, list("x"), match_condn = "any")
  lr_tbl_lm_spline <- test_lr(mod_lm_norm_ind, mod_lm_norm_int_ind, "x")
  expect_true(abs(wald_tbl_lm_spline$p - lr_tbl_lm_spline$p) < 0.3)

})
