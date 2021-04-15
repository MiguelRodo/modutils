test_that("multiplication works", {

  out_list <- list(null = "resp ~ 1 + (1|SubjectID) + offset(log(n_cell)) + Sex + DaysSinceEntry",
                   full = "resp ~ 1 + (1|SubjectID) + offset(log(n_cell)) + Sex + DaysSinceEntry + Progressor + QFT + splines::ns(timeFromMaxTimeToTB, knots = c(1, 3, 5))",
                   Progressor = "resp ~ 1 + (1|SubjectID) + offset(log(n_cell)) + Sex + DaysSinceEntry + QFT + splines::ns(timeFromMaxTimeToTB, knots = c(1, 3, 5))",
                   QFT = "resp ~ 1 + (1|SubjectID) + offset(log(n_cell)) + Sex + DaysSinceEntry + Progressor + splines::ns(timeFromMaxTimeToTB, knots = c(1, 3, 5))",
                   timeFromMaxTimeToTB = "resp ~ 1 + (1|SubjectID) + offset(log(n_cell)) + Sex + DaysSinceEntry + Progressor + QFT",
                   `Progressor~p~QFT` = "resp ~ 1 + (1|SubjectID) + offset(log(n_cell)) + Sex + DaysSinceEntry + splines::ns(timeFromMaxTimeToTB, knots = c(1, 3, 5))",
                   `Progressor~p~timeFromMaxTimeToTB` = "resp ~ 1 + (1|SubjectID) + offset(log(n_cell)) + Sex + DaysSinceEntry + QFT",
                   `QFT~p~timeFromMaxTimeToTB` = "resp ~ 1 + (1|SubjectID) + offset(log(n_cell)) + Sex + DaysSinceEntry + Progressor")

  # full
  expect_true(identical(create_mod_fml_as_chr(var_dep = "resp",
                                  var_exp = c("Progressor", "QFT"),
                                  var_exp_spline = var_exp_spline <- list("timeFromMaxTimeToTB" = list("fn" = "ns",
                                                                                                       "params" = list("knots" = c(1, 3, 5)))),
                                  var_re = "SubjectID",
                                  var_offset = "log(n_cell)",
                                  var_conf = c("Sex", "DaysSinceEntry"),
                                  rhs_text = NULL,
                                  int = TRUE,
                                  sub_exp = 2),
                        out_list))

  # removing only one explanatory variable
  expect_true(identical(create_mod_fml_as_chr(var_dep = "resp",
                                  var_exp = c("Progressor", "QFT"),
                                  var_exp_spline = var_exp_spline <- list("timeFromMaxTimeToTB" = list("fn" = "ns",
                                                                                                       "params" = list("knots" = c(1, 3, 5)))),
                                  var_re = "SubjectID",
                                  var_offset = "log(n_cell)",
                                  var_conf = c("Sex", "DaysSinceEntry"),
                                  rhs_text = NULL,
                                  int = TRUE,
                                  sub_exp = 1),
                        out_list[-c(6:8)]))

  # removing no explanatory variables
  expect_true(identical(create_mod_fml_as_chr(var_dep = "resp",
                                              var_exp = c("Progressor", "QFT"),
                                              var_exp_spline = var_exp_spline <- list("timeFromMaxTimeToTB" = list("fn" = "ns",
                                                                                                                   "params" = list("knots" = c(1, 3, 5)))),
                                              var_re = "SubjectID",
                                              var_offset = "log(n_cell)",
                                              var_conf = c("Sex", "DaysSinceEntry"),
                                              rhs_text = NULL,
                                              int = TRUE,
                                              sub_exp = 0),
                        out_list[1:2]))

  # check adding rhs_text works
  out_list_rhs <- list(null = "resp ~ 1 + (1|SubjectID) + offset(log(n_cell)) + Sex + DaysSinceEntry + BMI + splines::ns(Age)",
                       full = "resp ~ 1 + (1|SubjectID) + offset(log(n_cell)) + Sex + DaysSinceEntry + BMI + splines::ns(Age) + Progressor + QFT + splines::ns(timeFromMaxTimeToTB, knots = c(1, 3, 5))")
  expect_true(identical(create_mod_fml_as_chr(var_dep = "resp",
                                              var_exp = c("Progressor", "QFT"),
                                              var_exp_spline = var_exp_spline <- list("timeFromMaxTimeToTB" = list("fn" = "ns",
                                                                                                                   "params" = list("knots" = c(1, 3, 5)))),
                                              var_re = "SubjectID",
                                              var_offset = "log(n_cell)",
                                              var_conf = c("Sex", "DaysSinceEntry"),
                                              rhs_text = "BMI + splines::ns(Age)",
                                              int = TRUE,
                                              sub_exp = 0),
                        out_list_rhs))

  # checking that we can add no explanatory variables
  expect_true(identical(create_mod_fml_as_chr(var_dep = "resp",
                                              var_re = "SubjectID",
                                              var_offset = "log(n_cell)",
                                              var_conf = c("Sex", "DaysSinceEntry"),
                                              rhs_text = "BMI + splines::ns(Age)",
                                              int = TRUE,
                                              sub_exp = 0),
                        out_list_rhs[1]))
})
