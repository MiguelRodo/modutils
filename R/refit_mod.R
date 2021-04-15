#' @title Fit a model using "easier" optimisation settings
#'
#' @description Intended to automatically fit a model using easier settings.
#' @param mod_expr_as_chr character. Expression to run model as a character,
#' wrapped in a try call, e.g. \code{try(lm(y~x, data = data_mod))}.
#' @param family 'nb'. 'nb' denotes negative binomial.
#' @param pkg 'lme4'. Package used to fit model.
#' @param theta_fix numeric. If \code{NULL}, then negative binomial model is
#' refit using less stringent \code{tolPwrss} parameter (changing it from 1e-7 to 1e-6,
#' and then to 1e-5 if that doesn't fit).
#' @param env environment. Environment in which to fit model. Default is calling environment
#' of \code{refit_mod}.
#'
#' @details
#' For mixed-effects negative binomial models fit using the `lme4` package, if \code{is.null(theta_fix)},
#' then model is refit using less stringent \code{tolPwrss} parameter (changing it from 1e-7 to 1e-6,
#' and then to 1e-5 if that doesn't fit). If provided, then model is refit using \code{theta_fix} as
#' the fixed dispersion parameter.
#'
#' @export
#'
#' @examples
#' test_tbl <- data.frame(x = (1:10)^2, grp = rep(1:5, each = 2), dur = 1:10)
#' test_tbl$y <- (test_tbl$x + test_tbl$grp * 2)
#' suppressMessages(library(lme4))
#' mod_fml_as_chr <- modutils::get_mod_fml_as_chr(var_dep = "y", var_exp = "x", var_offset = "dur", var_re = "grp")[2]
#'
#' mod_expr_as_chr <- paste0("try(lme4::glmer.nb(", mod_fml_as_chr, ", data = test_tbl))")
#'
#' # refitting nb model with more lenient tolPwrss optimisation argument
#' mod_refit_tolpwrss <- suppressWarnings(refit_mod(mod_expr_as_chr, pkg = 'lme4', family = 'nb'))
#'
#' # refit with fixed theta
#' theta <- lme4::getME(mod_orig, "glmer.nb.theta")
#'
#' mod_refit_theta_fix <- suppressWarnings(refit_mod(mod_expr_as_chr, pkg = 'lme4', family = 'nb', theta_fix = theta/2))
#' lme4::getME(mod_refit_theta_fix, "glmer.nb.theta"))
refit_mod <- function(mod_expr_as_chr, pkg, family, theta_fix = NULL, env = rlang::caller_env()){
  force(env)
  if(pkg == 'lme4'){ # if pkg is lme4
    if(family == 'nb'){ # if fitting a negative binomial model
      if(is.null(theta_fix)){ # if no fixed theta available
        # adjust command to use more lenient tolPwrss
        mod_expr_as_chr <- paste0(stringr::str_sub(mod_expr_as_chr, end = -3),
                                  ", control = glmerControl(tolPwrss = 1e-6)))")
        mod <- eval(rlang::parse_expr(mod_expr_as_chr), envir = env)
        if(class(mod) == 'try-error'){
          # adjust command to use even more lenient tolPwrss
          mod_expr_as_chr <- paste0(stringr::str_sub(mod_expr_as_chr, end = -5),
                                    "5)))")
          mod <- eval(rlang::parse_expr(mod_expr_as_chr), envir = env)
        }
      } else{
        mod_expr_as_chr <- stringr::str_replace(mod_expr_as_chr, "glmer.nb",
                                                "glmer")
        mod_expr_as_chr <- paste0(stringr::str_sub(mod_expr_as_chr, end = -3),
                                  ", family = MASS::negative.binomial(theta = ",
                                  theta_fix, ")))")
        mod <- eval(rlang::parse_expr(mod_expr_as_chr), envir = env)
      }
      return(mod)
    }
  }
}
