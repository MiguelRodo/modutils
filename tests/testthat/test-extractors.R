test_that("multiplication works", {
  set.seed(2106); n_id <- 30
  test_tbl <- tibble::tibble(id = rep(as.character(1:n_id), each = 2),
                             grp = rep(as.character(1:2), each = n_id),
                             grp_add = rep(as.character(1:4), each = n_id/2),
                             y_re = rep(rnorm(n_id), each = 2),
                             x = rnorm(n_id * 2),
                             y = x^4 * ifelse(grp == "1", 3, 0) + y_re + rnorm(20, sd = 0.3))
  test_tbl$y <- test_tbl$y +  abs(min(test_tbl$y)) + 0.001
  test_tbl$y <- test_tbl$y/max(test_tbl$y + 0.001)
  test_tbl$y_01 <- test_tbl$y/diff(range(test_tbl$y))
  test_tbl$y_01 <- pmin(test_tbl$y_01, pmax(test_tbl$y_01, 0.99), 0.01)

  # check that it runs for each model type
  # ------------------------

  # fit model
  mod_lm_fix <- lm(y ~ -1, data = test_tbl)
  mod_lm_int <- lm(y ~ 1, data = test_tbl)
  mod_lm_cont <- lm(y ~ x, data = test_tbl)
  mod_lm <- lm(y ~ splines::ns(x, 3), data = test_tbl)
  mod_lmer <- lme4::lmer(y ~ x + (1|grp), data = test_tbl)
  mod_lmer_re2 <- lme4::lmer(y ~ x + (1|grp) + (x|grp_add), data = test_tbl)
  mod_lmertest <- lmerTest::lmer(y ~ x + (1|grp), data = test_tbl)
  mod_glmer <- suppressMessages(suppressWarnings(lme4::glmer(y ~ x + (1|grp), data = test_tbl,
                                                             family = 'binomial')))
  mod_beta <- betareg::betareg(y ~ x, data = test_tbl)

  # check that we extract log-likelihoods
  expect_true(is.numeric(ll(mod_lm)))
  expect_true(is.numeric(ll(mod_lmer)))
  expect_true(is.numeric(ll(mod_glmer)))
  expect_true(is.numeric(ll(mod_beta)))

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
  expect_identical(get_dof(mod_lmer_re2), 2L)
  expect_identical(get_dof(mod_glmer), 2L)

  # check that we extract the residuals
  expect_identical(class(get_resid(mod_lm)), c("tbl_df", "tbl", "data.frame"))
  expect_identical(colnames(get_resid(mod_lmer)), c(".resid_raw", ".resid_std"))
  expect_identical(nrow(get_resid(mod_glmer)), 60L)

  # check that we extract the predictions
  expect_identical(class(get_pred(mod_lm)), c("tbl_df", "tbl", "data.frame"))
  expect_identical(colnames(get_pred(mod_lmer)), c(".pred_resp", ".pred_lin", ".pred_resp_no_re", ".pred_lin_no_re"))
  expect_identical(nrow(get_pred(mod_glmer)), 60L)

  # check that we get the standard errors out

})
