#' @title Perform likelihood ratio test
#'
#' @param mod_large Model in which other model(s) are nested.
#' @param mod_small list or model. If a model, then coerced to a
#' list. Names of list elements are taken to be the names of the
#' variables missing in \code{mod_small} that are included in \code{mod_large}.
#' @param var character vector. If provided, then it
#' gives the names of the variables omitted in \code{mod_small}. If omitted,
#' then the names of the variables missing in \code{mod_small} are automatically
#' determined. Default is \code{NULL}.
#'
#'@return A tibble with the following columns:
#'\tabular{rrrrr}{
#'var \tab names of variables missing in mod small \cr
#'test \tab "lr" \cr
#'stat \tab likelihood ratio \cr
#'df \tab degrees of freedom \cr
#'p \tab p-value \cr
#'}
#'
#'@export
#'
#'@examples
#'  set.seed(2106); n_id <- 30
#' test_tbl <- tibble::tibble(id = rep(as.character(1:n_id), each = 2),
#'                           grp = rep(as.character(1:2), each = n_id),
#'                           y_re = rep(rnorm(n_id), each = 2),
#'                           x = rnorm(n_id * 2),
#'                           y = x^4 * ifelse(grp == "1", 3, 0) + y_re + rnorm(20, sd = 0.3))
#'test_tbl$y <- test_tbl$y +  abs(min(test_tbl$y)) + 0.001
#'test_tbl$y <- test_tbl$y/max(test_tbl$y + 0.001)
#'
#'mod_lm <- lm(y ~ splines::ns(x, 3), data = test_tbl)
#'mod_lm_int <- lm(y ~ 1, data = test_tbl)
#'test_lr(mod_lm, mod_lm_int, 'x')
test_lr <- function(mod_large, mod_small, var = NULL){

  if(!class(mod_small) == 'list'){
    mod_small <- list(mod_small)
  }
  if(is.null(var)) var <- names(mod_small)

  purrr::map_df(seq_along(mod_small), function(i){
    .test_lr(mod_large = mod_large, mod_small = mod_small[[i]],
             var = var[i])
  })
}

.test_lr <- function(mod_large, mod_small, var = NULL){
  # preparation
  # ----------------

  # extract log-likelihoods
  ll1 <- get_ll(mod_large)
  ll0 <- get_ll(mod_small)

  # extract degrees of freedom
  df1 <- get_dof(mod_large)
  df0 <- get_dof(mod_small)
  df <- df1 - df0
  if(df1 <= df0) stop("mod_large has fewer degrees of freedom than mod_small")

  # get coefficients that differ between the two
  if(is.null(var) || var == ""){
    coef1 <- get_coef(mod_large)$var
    coef0 <- get_coef(mod_small)$var
    var <- setdiff(coef1, coef0)
  }

  # perform test
  # ----------------

  # calculate test statistic
  lr_stat <- 2 * (ll1 - ll0)
  df_test <- df1 - df0
  lr_p <- pchisq(lr_stat, df = df_test, lower.tail = FALSE)

  tibble::tibble(var = ifelse(is.null(var), NA_character_,
                              paste0(var, collapse = "; ")),
                 test = "lr",
                 stat = lr_stat,
                 df = df_test,
                 p = lr_p)
}

#' @title Generic function to extract log-likelihoods
get_ll <- function(mod) UseMethod("get_ll")
get_ll.lm <- function(mod) logLik(mod)[1]
get_ll.lmerMod <- function(mod){
  summary(mod)$logLik[1]
}
get_ll.glmerMod <- function(mod){
  summary(mod)$logLik[1]
}

#' @title Generic function to extract estimates
get_coef <- function(mod) UseMethod("get_coef")
get_coef.lm <- function(mod){
  coef_mat <- coef(summary(mod))
  tibble::tibble(var = rownames(coef_mat),
                 est = coef_mat[,"Estimate"],
                 se = coef_mat[,"Std. Error"],
                 test_stat = coef_mat[,"t value"],
                 p = coef_mat[,"Pr(>|t|)"])
}
get_coef.lmerMod <- function(mod){
  coef_mat <- coef(summary(mod))
  tibble::tibble(var = rownames(coef_mat),
                 est = coef_mat[,"Estimate"],
                 se = coef_mat[,"Std. Error"],
                 test_stat = coef_mat[,"t value"],
                 p = NA)
}
get_coef.glmerMod <- function(mod){
  coef_mat <- coef(summary(mod))
  tibble::tibble(var = rownames(coef_mat),
                 est = coef_mat[,"Estimate"],
                 se = coef_mat[,"Std. Error"],
                 test_stat = coef_mat[,"z value"],
                 p = coef_mat[,"Pr(>|z|)"])
}



#' @title Generic function to extract fixed-effects degrees of freedom
get_dof <- function(mod) UseMethod("get_dof")
get_dof.lm <- function(mod) length(mod$coefficients)
get_dof.lmerMod <- function(mod) ncol(coef(mod)$grp)
get_dof.glmerMod <- function(mod) ncol(coef(mod)$grp)
