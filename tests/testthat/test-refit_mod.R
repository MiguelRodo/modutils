test_that("refit_mod works", {

  test_tbl <- data.frame(x = (1:10)^2, grp = rep(1:5, each = 2), dur = 1:10)
  test_tbl$y <- (test_tbl$x + test_tbl$grp * 2)

  # ===================
  # Mixed-effects models - lme4
  # ===================

  on.exit(try(detach("package:lme4", unload = TRUE), silent = TRUE))
  on.exit(try(detach("package:lme4", unload = TRUE), silent = TRUE),
          add = TRUE)
  suppressMessages(library(lme4))

  # Negative binomial
  # ------------------

  # model formula as chr
  mod_fml_as_chr <- modutils::get_mod_fml_as_chr(var_dep = "y",
                                                 var_exp = "x",
                                                 var_offset = "dur",
                                                 var_re = "grp")[2]

  # model call as chr
  mod_expr_as_chr <- paste0("try(lme4::glmer.nb(", mod_fml_as_chr, ", data = test_tbl))")
  # fit model
  mod_orig <- suppressWarnings(eval(rlang::parse_expr(mod_expr_as_chr)))

  # refit
  mod_refit_tolpwrss <- suppressWarnings(refit_mod(mod_expr_as_chr, pkg = 'lme4', family = 'nb'))
  expect_true(class(mod_refit_tolpwrss) != 'try-error')

  # refit with fixed theta
  theta <- lme4::getME(mod_orig, "glmer.nb.theta")
  mod_refit_theta_fix <- suppressWarnings(refit_mod(mod_expr_as_chr, pkg = 'lme4', family = 'nb',
                                                    theta_fix = theta/2))
  expect_equal(theta/2, lme4::getME(mod_refit_theta_fix, "glmer.nb.theta"))

})
