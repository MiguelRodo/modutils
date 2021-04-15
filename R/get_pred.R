if(FALSE){
  library(lme4)
  library(merTools)
  library(tibble)
  library(magrittr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  theme_set(theme_cowplot())



  set.seed(2106); n_id <- 30
  test_tbl <- tibble::tibble(id = rep(as.character(1:n_id), each = 2),
                             grp = rep(as.character(1:2), each = n_id),
                             y_re = rep(rnorm(n_id), each = 2),
                             x = rnorm(n_id * 2),
                             z = rnorm(n_id * 2),
                             y = x^4 * ifelse(grp == "1", 3, 0) + z +  y_re + rnorm(20, sd = 0.3)) %>%
    mutate(y = y + abs(min(y)) + 0.001) %>%
    mutate(y = y/max(y + 0.001))


  get_link_inv_fn <- function(mod) UseMethod("get_link_inv_fn")
  get_link_inv_fn.lm <- function(mod) function(x) x
  get_link_inv_fn.lmerMod <- function(mod) function(x) x
  get_link_inv_fn.glmerMod <- function(mod){
    link <- summary(mod)$link
    switch(link,
           "inverse" = function(x) 1/x,
           "log" = exp,
           "logit" = function(x) 1/(1+exp(-x)),
           stop(paste0("link ", link, " not recognised in get_link_inv_fn.glmer")))
  }


  get_fml <- function(mod_sum) UseMethod("get_fml")
  get_fml.summary.merMod <- function(mod_sum) mod_sum$call$formula


  get_mod_cov <- function(mod_sum){
    fml <- get_fml(mod_sum = mod_sum)
    fml_comps_rhs <- stringr::str_split(fml, pattern = "[[~]]|[[+]]")[[3]]
    cov_vec <- purrr::map(fml_comps_rhs, function(x){
      if(stringr::str_detect(x, "[[|]]")) return(NULL)
      stringr::str_trim(x)
    }) %>%
      purrr::compact() %>%
      as.character()
    # remove references to removing or specifying intercept
    cov_vec[!cov_vec %in% c("1", "-1")]

  }

  get_coef <- function(mod_sum) UseMethod("get_coef")
  get_coef.summary.merMod <- function(mod_sum) coefficients(mod_sum)[,"Estimate"]
  get_coef.summary.lm <- function(mod_sum) coefficients(mod_sum)
  #' @title Get linear predictor
  get_lp <- function(cov, coef, data){
    # get vector to add to
    lp_vec <- rep(0, nrow(data))

    # add intercept
    if(any(grepl("Intercept", names(coef)))){
      lp_vec <- lp_vec + coef[[grep("Intercept", names(coef))]]
    }

    # add other covariates
    for(x in cov){
      # add covariate if it appears
      # exactly in coefficients, i.e. continuous covariates
      if(x %in% names(coef)){
        if(!x %in% colnames(data)){
          stop(paste0("covariate ", x, " not found in colnames(data) in get_lp"))
        }
        lp_vec <- lp_vec + data[[x]] * coef[[grep(x, names(coef))]]
      }

      # if there are bs or ns splines specified as their own terms, will need to do that
      if(stringr::str_detect(x, "ns[[(]]|bs[[(]]")){
        spline_command_loc <- stringr::str_locate(x, "ns[[(]]|bs[[(]]")[,"end"][[1]]
        x_after_spline_command <-
          first_comma_after_bs_or_ns_loc <- str_locate(stringr::str_sub(x,
                                                                        bs_or_ns_loc+ 1))
      }

      # get indices that match to current covariate
      # - may be more than one as current covariate has levels
      cov_ind_vec <- grep(x, names(coef))

      # stop is covariate not found in coefficient names
      if(length(cov_ind_vec) == 0){
        stop(paste0("covariate ", x, " not found in names(coef) in get_lp"))
      }

      # for each individaul that has right level, add corresponding coefficient
      for(cov_ind in cov_ind_vec){
        coef_curr <- names(coef)[cov_ind] # current coefficient name
        cov_level <- stringr::str_remove(coef_curr, x) # level of covariate
        lp_vec <- lp_vec + coef[[cov_ind]] * (data[[x]] == cov_level) # add coefficient for all that match right level
      }
    }

    lp_vec
  }

  # y - outcome
  # x - spline term
  # z - extra term
  # id - random effect term
  # grp - interacts with splines

  mod_lm <- lm(y ~ x, data = test_tbl)
  mod_lmer <- lme4::lmer(y ~ x + (1|id), data = test_tbl)
  mod_glmer <- lme4::glmer(y ~ x + grp + (1|id), data = test_tbl,
                           family = 'Gamma')
  mod_sum_lmer <- summary(mod_lmer)

  get_pred <- function(mod, data = NULL, type = "response"){

    # preparation
    # ------------------

    # get model summary info
    mod_sum <- summary(mod)

    # get components needed for calculation
    cov_vec <- get_mod_cov(mod_sum = mod_sum) # non-re and non-int model covariates
    coef_vec <- get_coef(mod_sum = mod_sum) # model coefficients

    # prepare data used for prediction
    if(is.null(data)) data <- mod@frame
    if(class(data)[1] == 'matrix') data <- as.data.frame(data)

    # get predictions
    # --------------------

    # get linear predictors
    lp_vec <- get_lp(cov = cov_vec, coef = coef_vec, data = data)
    if(type == 'link') return(lp_vec) # return if that is right scale

    # get predictions on response scale
    inv_link_fn <- get_link_inv_fn(mod = mod)
    inv_link_fn(lp_vec)

  }


  get_pred(mod = mod_glmer, data = test_tbl)



  plot_tbl <- test_tbl %>%
    mutate(pred_lm = predict(mod_lm),
           pred_lmer = predict(mod_lmer),
           #pred_lmer_man = get_pred(mod = mod_lmer,
           #                         data = test_tbl),
           pred_glmer_man = get_pred(mod = mod_glmer,
                                     data = test_tbl),
           pred_glmer = predict(mod_glmer, type = 'response'),
           pred_glmer_newdata = predict(mod_glmer, newdata = test_tbl %>%
                                          mutate(id = "1"),
                                        type = 'resp')) %>%
    pivot_longer(-(id:y),
                 names_to = "type",
                 values_to = "pred")

  ggplot(plot_tbl %>%
           filter(type == 'pred_glmer_newdata')) +
    geom_point(aes(x, y, col = grp)) +
    geom_line(aes(x, y = pred, col = grp, linetype = type), size = 2, alpha = 0.5) +
    scale_colour_brewer(palette = "Set1")



}
