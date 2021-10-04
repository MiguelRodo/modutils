#' @title Extract log-likelihoods
#'
#' @description A wrapper around \code{logLik} essentially. Implemented
#' just in case there is a model that does not have a \code{logLik} generic.
#'
#' @param mod All objects that have a method for them in \code{broom::glance} or
#' \code{broom.mixed::glance} or \code{modutils:::.;;} or \code{stats::logLik}.
#' Run \code{library(broom.mixed); utils::methods(broom:glance)}
#' or \code{library(stats); utils::methods(stats::logLik)} to see a list.
#'@param get_ll function. A function that takes a model as input and returns
#'the log-likelihood. Provided so that the user need not write a method for
#'a package if the log-likelihood function is not available.
#' @return Numeric vector of length one containing the log-likelihood.
#'
#' @export
#'
#' @examples
#' set.seed(3)
#' test_tbl <- data.frame(x = 1:10)
#' test_tbl$y <- 2 * test_tbl$x + rnorm(10)
#' mod_lm <- lm(y ~ x, data = test_tbl)
#' ll(mod_lm)
ll <- function(mod, get_ll = NULL){
  if(!is.null(get_ll)){
    ll <- get_ll(mod)
    if(!is.numeric(ll)) stop(paste0("Function get_ll supplied to ll returns an object of class ",
                                    paste0(class(mod), collapse = "/"), "when it should be a numeric vector of length 1."))
    if(length(ll) != 1) stop(paste0("Function get_ll supplied to ll returns an object of length ",
                                    length(ll), "insted of length 1."))
    return(ll)
  }
  ll <-  ll <- try(logLik(mod)[[1]], silent = TRUE)
  if(class(ll) != 'try-error'){
    return(ll)
  } else{
    ll <- try(broom.mixed::glance(mod)$logLik, silent = TRUE)
  }
  if(class(ll) != 'try-error'){
    return(ll)
  } else try(modutils::.ll(mod), silent = TRUE)


  if(class(ll) == 'try-error'){
    stop(paste0("ll not available for models of class",
                paste0(class(mod), collapse = "/"), ". Add method for such a model for generic modutils::.ll or pass such a function to the get_ll parameter of ll"))
  }
  ll
}

#' @title Extract log-likelihood
.ll <- function(mod) UseMethod(".ll")


#' @title Extract model estimates
#' @export
get_coef <- function(mod) UseMethod("get_coef")
#' @export
get_coef.lm <- function(mod){
  coef_mat <- coef(summary(mod))
  tibble::tibble(var = rownames(coef_mat),
                 est = coef_mat[,"Estimate"],
                 se = coef_mat[,"Std. Error"],
                 test_stat = coef_mat[,"t value"],
                 p = coef_mat[,"Pr(>|t|)"])
}
#' @export
get_coef.lmerMod <- function(mod){
  coef_mat <- coef(summary(mod))
  tibble::tibble(var = rownames(coef_mat),
                 est = coef_mat[,"Estimate"],
                 se = coef_mat[,"Std. Error"],
                 test_stat = coef_mat[,"t value"],
                 p = NA)
}
#' @export
get_coef.glmerMod <- function(mod){
  coef_mat <- coef(summary(mod))
  zt_val <- ifelse(any(stringr::str_detect(colnames(coef_mat), "t value")),
                   "t value",
                   "z value")
  zt_p <- ifelse(any(stringr::str_detect(colnames(coef_mat), "Pr([[>]][[|]]t[[|]])")),
                   "Pr(>|t|)",
                   "Pr(>|z|)")
  tibble::tibble(var = rownames(coef_mat),
                 est = coef_mat[,"Estimate"],
                 se = coef_mat[,"Std. Error"],
                 test_stat = coef_mat[,zt_val],
                 p = coef_mat[,zt_p])
}

#' @export
get_coef.glmmTMB <- function(mod){
  coef_mat <- coef(summary(mod))$cond
  tibble::tibble(var = rownames(coef_mat),
                 est = coef_mat[,"Estimate"],
                 se = coef_mat[,"Std. Error"],
                 test_stat = coef_mat[,"z value"],
                 p = coef_mat[,"Pr(>|z|)"])
}


#' @title Extract fixed-effects degrees of freedom
#' @export
get_dof <- function(mod) UseMethod("get_dof")
#' @export
get_dof.lm <- function(mod) length(mod$coefficients)
#' @export
get_dof.lmerMod <- function(mod) ncol(coef(mod)[[1]])
#' @export
get_dof.glmerMod <- function(mod) ncol(coef(mod)[[1]])
#' @export
get_dof.glmmTMB <- function(mod) nrow(coef(summary(mod))$cond)

#' @title Extract residuals
#'
#'
#'
#' @param mod object of class 'lm'.
#'
#' @details
#' Notes:
#' broom::augment.lm gives .resid_std column as standardised deviance residuals (standardised by standard error calculated using all samples)
#' deviance is good as it's more N(0,1) (when standardised, I would guess) than Pearson (can't remember what that is now, maybe the working?)

#'
#' @return A tibble with an entry for each observation for each of the following columns:
#' \tabular{rr}{
#' .resid_raw \tab raw residuals (predicted - actual) \cr
#' .resid_std \tab standardised residuals ((predicted - actual)(standard deviation of response)) \cr
#' }
#' @export
get_resid <- function(mod) UseMethod("get_resid")
#' @export
get_resid.lm <- function(mod) tibble::tibble(.resid_raw = residuals(mod, type = 'deviance'),
                                             .resid_std = residuals(mod, type = 'deviance')/sigma(mod))
#' @export
get_resid.lmerMod <- function(mod) tibble::tibble(.resid_raw = residuals(mod, type = 'deviance', scaled = FALSE),
                                                  .resid_std = residuals(mod, type = 'deviance', scaled = TRUE))
#' @export
get_resid.glmerMod <- function(mod) tibble::tibble(.resid_raw = residuals(mod, type = 'deviance', scaled = FALSE),
                                                   .resid_std = residuals(mod, type = 'deviance', scaled = TRUE))
#' @export
get_resid.default <- function(mod){
  broom_augment_tbl <- try(broom.mixed::augment(mod),
                           silent = TRUE)
  if(class(broom_augment_tbl)[1] == 'try-error'){
    resid_vec_raw <- try(residuals(mod), silent = TRUE)
    if(class(resid_vec_raw)[1] == 'try-error'){
      return(tibble::tibble(
        .resid_raw = NA_real_,
        .resid_std = NA_real_
      )[-1,])
    }
    tibble(.resid_raw = resid_vec_raw,
           .resid_std = resid_vec_raw/sd(resid_vec_raw))
  }
  if(!'.std.resid' %in% colnames(broom_augment_tbl)){
    broom_augment_tbl %<>%
      mutate(.resid_std = .resid/sd(.resid))
    try(broom_augment_tbl %<>% rename(.resid_raw = .resid),
        silent = TRUE)
  } else{
    broom_augment_tbl %<>%
      rename(.resid_std = .std.resid)
  }

  broom_augment_tbl %>%
    select(.resid_raw, .resid_std)
}

#' @title Extract predicted values
#'
#' @param mod model object of class 'lm', 'lmerMod' or 'glmerMod'.
#'
#' @return
#' A tibble with column \code{pred}.
#' @export
get_pred <- function(mod, ...) UseMethod("get_pred")

#' @export
get_pred.lm <- function(mod, ...) tibble::tibble(.pred_resp = predict(mod, ...),
                                                 .pred_lin = predict(mod, ...),
                                                 .pred_resp_no_re =  predict(mod, ...),
                                                 .pred_lin_no_re =  predict(mod, ...))

#' @export
get_pred.lmerMod <- function(mod, ...) tibble::tibble(.pred_resp = predict(mod, type = 'response', ...),
                                                      .pred_lin = predict(mod, type = 'link', ...),
                                                      .pred_resp_no_re = predict(mod, type = 'response', re.form = ~0, ...),
                                                      .pred_lin_no_re = predict(mod, type = 'link', re.form = ~0, ...))

#' @export
get_pred.glmerMod <- function(mod, ...) tibble::tibble(.pred_resp = predict(mod, type = 'response', ...),
                                                        .pred_lin = predict(mod, type = 'link', ...),
                                                        .pred_resp_no_re = predict(mod, type = 'response', re.form = ~0, ...),
                                                        .pred_lin_no_re = predict(mod, type = 'link', re.form = ~0, ...))
#' @export
get_pred.glmmTMB <- function(mod, ...) tibble::tibble(.pred_resp = predict(mod, type = 'response', ...),
                                                      .pred_lin = predict(mod, type = 'link', ...),
                                                      .pred_resp_no_re = predict(mod, type = 'response', re.form = ~0, ...),
                                                      .pred_lin_no_re = predict(mod, type = 'link', re.form = ~0, ...))

#' @export
get_pred.default <- function(mod, ...){
  return(tibble::tibble(
    .pred_resp = NA_real_,
    .pred_lin = NA_real_,
    .pred_resp_no_re = NA_real_,
    .pred_resp_no_re = NA_real_
  )[-1,])
  broom_augment_tbl <- try(broom.mixed::augment(mod),
                           silent = TRUE)
  if(class(broom_augment_tbl)[1] == 'try-error'){
    pred_vec <- try(predict(mod), silent = TRUE)
    if(class(pred_vec)[1] == 'try-error'){
      return(tibble::tibble(
        .pred_resp = NA_real_,
        .pred_lin = NA_real_,
        .pred_resp_no_re = NA_real_,
        .pred_resp_no_re = NA_real_
      )[-1,])
    }
    tibble(.pred = resid_vec_raw,
           .resid_std = resid_vec_raw/sd(resid_vec_raw))
  }
  if(!'.resid_std' %in% colnames(broom_augment_tbl)){
    broom_augment_tbl %<>%
      mutate(.resid_std = .resid/sd(.resid))
    try(broom_augment_tbl %<>% rename(.resid_raw = .resid),
        silent = TRUE)
  }

  broom_augment_tbl %>%
    select(.resid_raw, .resid_std)
}

#' @title Extract the standard error of prediction
#'
#' @param mod model object of class 'lm', 'lmerMod' or 'glmerMod'.
get_se <- function(fit, data){

  # preparation
  # ==================

  # checks for availability
  if(missing(data)) stop("data missing in get_se")
  if(missing(fit)) stop("fit missing in get_se")

  # ensure data is a matrix or a data frame
  if(!class(data) %in% c("matrix", "data.frame")){
    stop(paste0("data has class ", paste0(class(data), collapse = " "), "instead of data.frame or matrix in get_se"))
  }


}

get_est_and_vcov <- function(fit){
  UseMethod("get_est_and_vcov")
}

get_est_and_vcov.glmerMod <- function(fit){
  est_vec <- coefficients(summary(fit))[,"Estimate"]
  vcov_mat <- as.matrix(vcov(fit))
  list('est' = est_vec, 'vcov' = vcov_mat)
}

get_est_and_vcov.glmmTMB <- function(fit){
  est_vec <- coefficients(summary(fit))$cond[,"Estimate"]
  vcov_mat <- vcov(fit)$cond
  list('est' = est_vec, 'vcov' = vcov_mat)
}

get_est_and_vcov.lmerMod <- function(fit){
  est_vec <- coefficients(summary(fit))[,"Estimate"]
  vcov_mat <- as.matrix(vcov(fit))
  list('est' = est_vec, 'vcov' = vcov_mat)
}

get_est_and_vcov.lm <- function(fit){
  est_vec <- coefficients(fit)
  vcov_mat <- vcov(fit)
  list('est' = est_vec, 'vcov' = vcov_mat)
}

get_est_and_vcov.list <- function(fit){
  est_vec <- fit$est
  vcov_mat <- fit$vcov
  list('est' = est_vec, 'vcov' = vcov_mat)
}
