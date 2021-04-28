test_that("test_wald works", {

  # test data
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
  mod_lm <- lm(y ~ splines::ns(x, 3), data = test_tbl)
  mod_lmer <- lme4::lmer(y ~ x + (1|grp), data = test_tbl)
  mod_glmer <- suppressMessages(suppressWarnings(lme4::glmer(y ~ x + (1|grp), data = test_tbl,
                                            family = 'binomial')))

  # check that a tibble is returned
  expect_identical(class(test_wald(fit = mod_lm, var = list('0')))[1],
                   "tbl_df")
  expect_identical(class(test_wald(fit = mod_lmer, var = list('0')))[1],
                   "tbl_df")
  expect_identical(class(test_wald(fit = mod_glmer, var = list('0')))[1],
                   "tbl_df")

  # check that continuous var matches default output
  # and that the p-value is the same as if you square root
  # the stat and compare to N(0,1) dbn
  # ------------------------

  # fit models
  mod_lm_x <- lm(y ~ x, data = test_tbl)
  mod_lmer_x <- lme4::lmer(y ~ x + (1|grp), data = test_tbl)
  mod_glmer_x <- suppressMessages(suppressWarnings(lme4::glmer(y ~ x + (1|grp), data = test_tbl,
                                                             family = 'binomial')))

  # lm
  lm_stat_x <- coefficients(summary(mod_lm_x))["x","t value"]
  lm_tbl_x_wald <- test_wald(fit = mod_lm_x, var = list("x"))
  expect_equal(lm_stat_x, sqrt(lm_tbl_x_wald$stat[1]))
  lm_p_x_norm <- 2 * pnorm(sqrt(lm_tbl_x_wald$stat[1]), lower.tail = FALSE)
  expect_equal(lm_p_x_norm, lm_tbl_x_wald$p[1])

  # lmer
  lmer_stat_x <- coefficients(summary(mod_lmer_x))["x","t value"]
  lmer_tbl_x_wald <- test_wald(fit = mod_lmer_x, var = list("x"))
  expect_equal(lmer_stat_x, sqrt(lmer_tbl_x_wald$stat))
  expect_identical(1L, lmer_tbl_x_wald$df)
  lmer_p_x_norm <- 2 * pnorm(sqrt(lmer_tbl_x_wald$stat), lower.tail = FALSE)
  expect_equal(lmer_p_x_norm, lmer_tbl_x_wald$p[1])

  # glmer
  glmer_stat_x <- coefficients(summary(mod_glmer_x))["x","z value"]
  glmer_tbl_x_wald <- test_wald(fit = mod_glmer_x, var = list("x"))
  expect_equal(abs(glmer_stat_x), sqrt(glmer_tbl_x_wald$stat[1])) #
  expect_identical(1L, glmer_tbl_x_wald$df)
  glmer_p_x_norm <- 2 * pnorm(sqrt(glmer_tbl_x_wald$stat), lower.tail = FALSE)
  expect_equal(glmer_p_x_norm, glmer_tbl_x_wald$p[1])

  # check that iintercept stat matches default output
  # and that the p-value is the same as if you square root
  # the stat and compare to N(0,1) dbn
  # ------------------------

  # fit models
  set.seed(3)
  test_tbl$y_straight <- test_tbl$x + rnorm(nrow(test_tbl), sd = 0.001)
  mod_lm_x_lin <- lm(y_straight ~ x, data = test_tbl)
  mod_lm_x_spline_df_1 <- lm(y_straight ~ splines::ns(x), data = test_tbl)
  mod_lm_x_spline_df_3 <- lm(y_straight ~ splines::ns(x, df = 3), data = test_tbl)


  # lm
  lm_stat_x_lin <- coefficients(summary(mod_lm_x_lin))["x","t value"]
  lm_tbl_x_wald_1 <- test_wald(fit = mod_lm_x_spline_df_1, var = list("x"),
                                  match_condn = 'any')
  expect_equal(lm_stat_x_lin, sqrt(lm_tbl_x_wald_1$stat))
  lm_tbl_x_wald_3 <- test_wald(fit = mod_lm_x_spline_df_3, var = list("x"),
                                  match_condn = 'any')
  # numbers are huge and models are different, so a difference of less than five
  # percent in the actual stat seems fine
  expect_true(abs(diff(c(lm_stat_x_lin, sqrt(lm_tbl_x_wald_3$stat))))/lm_stat_x_lin < 0.05)
  expect_identical(1L, lm_tbl_x_wald_1$df)
  expect_identical(3L, lm_tbl_x_wald_3$df)

  # check that we get warnings as we should
  # -------------------

  # error if match_condn = 'exact'
  expect_error(test_wald(fit = mod_lm_x_spline_df_1, var = list("x"),
                              match_condn = 'exact'))

  # check that we can match when we use non-alphanumeric characters
  expect_identical(lm_tbl_x_wald_3,
                   test_wald(fit = mod_lm_x_spline_df_3, var = list("splines::ns(x"),
                                  match_condn = 'start'))


})

test_that(".bracket_non_an_chr works", {
  expect_identical(.bracket_non_an_chr("abc"), "abc")
  # check that non-an chr's aren't being removed
  expect_false(identical(.bracket_non_an_chr("ab(c"), "abc"))
  # check that we replace one
  expect_identical(.bracket_non_an_chr("ab:c"), "ab[[:]]c")
  expect_false(identical(.bracket_non_an_chr("ab(c"), "ab(c"))
  # check that we replace all
  expect_identical(.bracket_non_an_chr("!ab(c@"), "[[!]]ab[[(]]c[[@]]")
  # check that sequences are not bracketed inside
  expect_identical(.bracket_non_an_chr("!ab((*#c@@"), "[[!]]ab[[(]]{2}[[*]][[#]]c[[@]]{2}")
  # check that it works on vector input, too
  expect_identical(.bracket_non_an_chr(c("!ab(c@", "abc")), c("[[!]]ab[[(]]c[[@]]", "abc"))

  # check that we match when we expect to match
  str_out <- .bracket_non_an_chr("@@")
  expect_true(stringr::str_detect("@@", str_out))
  expect_false(stringr::str_detect("@", str_out))
})
