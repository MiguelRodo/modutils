#' @title Extract log-likelihoods
#' @export
get_ll <- function(mod) UseMethod("get_ll")
#' @export
get_ll.lm <- function(mod) logLik(mod)[1]
#' @export
get_ll.lmerMod <- function(mod){
  summary(mod)$logLik[1]
}
#' @export
get_ll.glmerMod <- function(mod){
  summary(mod)$logLik[1]
}

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
