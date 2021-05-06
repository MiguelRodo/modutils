#' @title Create model formula as character from components
#'
#' @param var_dep character. Name of dependent variable. Goes on left-hand side of formula.
#' @param int logical. If \code{TRUE}, then formula includes an intercept. Default is \code{TRUE}.
#' @param var_re character vector. Specifies individual random effects in the \code{lme4} notation..
#' If named, then the name of each
#' element specifies the variable that has a random effect for each level of the corresponding element. If unnamed,
#' then it is assumed that the intercept has a random effect. For example, \code{var_re = c("time" = "x1",
#' "x2")} adds \code{(time|x1) + (1|x2)} to the formula. If blank, then no random effect terms are added.
#' @param var_offset character vector. Text to be added as \code{"offset(<var_offset>)"},
#' e.g. \code{var_offset = "log(days)"} adds \code{offset(log(days))} to the model. If blank,
#' then no offset term is added. That that this can be a character vector, but it's probably
#' the case that \code{lme4::{g}lmer}, at least, requires that there is only one offset term.
#' @param var_conf character vector. Specifies confounders to be added.
#' @param rhs_text character. Text simply to be concatenated to the model formula. Note that "+"
#' is added before adding this.
#' @param sub_exp non-negative integer. Specifies maximum number of explanatory variables
#' to be excluded from the model formulae. If for example \code{sub_exp = 1}, then
#' model formula as strings are returned where each explanatory variable is excluded on its own (but
#' all other explanatory variables of interest and other terms (such as confounders and random effects))
#' are kept. If greater than one, then all possible combinations of explanatory variables to exclude
#' are excluded, along with all possible combinations of explanatory variables to exclude of size less than this.
#' @param var_exp character vector. Specifies explanatory variables to be included as is in model.
#' @param var_exp_spline named list of named list(s).. For each element, the name denotes then name of the explanatory variable
#' for which a spline term is required. Each element is a named list, where the following names denote the following things.
#' The name \code{"splines_pkg_fn"} specifies the spline function to be used from the \code{splines} pkg, e.g.
#' \code{splines_pkg_fn = "ns"} means that the \code{ns} function from the \code{splines} pkg is used. The
#' \code{knots}-named element specifies the knots to be used. No other named elements are supported at present
#' (but can be easily added).
#'
#' @details This was written for the following reasons:
#' - So that if we want to fit multiple, similar models (for example, the null model, the full model
#' and all models less one variable) the formulae could be obtained quickly.
#' - So that if we want to link model formula to character arguments in functions, it can be done easily.
#' This is useful if, for example, the same function that fits the model is also used to create
#' plots of model output, as the key argument (the variables for which we want plots) can be
#' used to specify the model formula and then be used to create the plots.
#'
#' @return A named list, where the name specifies the type of formula and the
#' corresponding element is the model formula as a stringr. The full model formula (with all explanatory variables)
#' has name \code{'full'}. The null model (with zero explanatory variables of interest but
#' including all random effects, offsets, confounders and any additional text specified by \code{rhs_text}) has name
#' null. The full model less one or more explanatory variables is denoted by the name of the excluded
#' explanatory variables, separated by \code{~p~} (the 'p' denoted plus).
#'
#' @importFrom magrittr %<>%
#'
#' @export
get_mod_fml_as_chr <- function(var_dep, var_exp = NULL, var_exp_spline = NULL, var_re = NULL,
                               var_offset = NULL, var_conf = NULL,
                               rhs_text = NULL,
                               int = TRUE, sub_exp = 0){

  # set dependent variable
  lhs <- var_dep

  # specify intercept or not
  rhs_null <- ifelse(int, "1", "-1")

  # random effects
  # -----------------
  if(!is.null(var_re)){
    # set names equal to "1" if not specified already
    if(is.null(names(var_re))) names(var_re) <- rep("1", length(var_re))
    # add each term as a random effect in turn
    for(i in seq_along(var_re)){
      # covariate whose coefficient changes for each level of var_re
      re_cov <- ifelse(names(var_re)[i] == "", "1", names(var_re)[i])
      # random effect grouping factor
      re_grp <- var_re[[i]]
      # add to model formula as chr
      rhs_null <- paste0(rhs_null, " + (", re_cov, "|", re_grp, ")")
    }
  }

  # offset terms
  # -------------------
  for(i in seq_along(var_offset)){
    rhs_null <- paste0(rhs_null, " + offset(", var_offset[i], ")")
  }

  # confounders
  # -------------------
  if(!is.null(var_conf)){
    rhs_null <- paste0(rhs_null, " + ", paste0(var_conf, collapse = " + "))
  }

  # additional text
  # -------------------
  if(!is.null(rhs_text)){
    rhs_null <- paste0(rhs_null, " + ", rhs_text)
  }

  if(is.null(var_exp) && is.null(var_exp_spline)) return(list(null = paste0(lhs, " ~ ", rhs_null)))

  # explanatory variables
  # -------------------

  # non-splines
  rhs_comp_exp_list <- switch(as.character(is.null(var_exp)),
                              "TRUE" = list(),
                              "FALSE" = setNames(as.list(var_exp), var_exp))

  # splines
  for(i in seq_along(var_exp_spline)){
    # specify splines pkg, its fn and the variable we need a splines basis for
    rhs_spline_add <- paste0("splines::", var_exp_spline[[i]]$fn,
                             "(", names(var_exp_spline)[i])
    # add additional parameters if available
    if('params' %in% names(var_exp_spline[[i]])){
      # extract parameters
      params_curr <- var_exp_spline[[i]]$params
      # loop over parameters
      for(j in seq_along(params_curr)){
        # get parameter name
        param_nm <- names(params_curr)[j]
        # if parameter is knots, then you expect a numeric vector,
        # so add argument accordingly
        if(param_nm == 'knots'){
          param_arg <- paste0("c(", paste0(params_curr[[j]], collapse = ", "), ")")
          rhs_spline_add <- paste0(rhs_spline_add, ", ", param_nm, " = ", param_arg)
        }
        if(param_nm == "df"){
          param_arg <- params_curr[[j]]
          rhs_spline_add <- paste0(rhs_spline_add, ", ", param_nm, " = ", param_arg)
        }

      }
    }
    # add closing bracket for this spline term
    rhs_spline_add <- paste0(rhs_spline_add, ")")

    rhs_comp_exp_list <- rhs_comp_exp_list %<>% append(setNames(list(rhs_spline_add),
                                                                names(var_exp_spline)[i]))
  }

  # create list of rhs elements
  rhs_comp_exp <- paste0(rhs_comp_exp_list, collapse = " + ")
  rhs_full <- paste0(rhs_null, " + ", rhs_comp_exp)
  rhs_list <- list(null = rhs_null,
                   full = rhs_full)

  # create nested models of full model
  if(sub_exp > 0 && length(rhs_comp_exp_list) > 1){
    n_var <- length(rhs_comp_exp_list)
    for(n_exc in 1:min(sub_exp, n_var - 1)){
      #exc_mat <- gtools::permutations(n_var, n_exc)
      exc_mat <- combn(x = 1:n_var, m = n_exc)
      for(j in 1:ncol(exc_mat)){
        #rhs_comp_exp_curr <- rhs_comp_exp_list[-exc_mat[j,]]
        rhs_comp_exp_inc <- rhs_comp_exp_list[-exc_mat[,j]]
        rhs_list_add_nm <- paste0(names(rhs_comp_exp_list[exc_mat[,j]]),
                                  collapse = "~p~")
        rhs_list_add_fml <- paste0(rhs_null, " + ",
                                   paste0(rhs_comp_exp_inc, collapse = " + "))
        rhs_list %<>% append(setNames(list(rhs_list_add_fml),
                                      rhs_list_add_nm))
      }
    }
  }

  fml_as_chr_list <- purrr::map(rhs_list, function(rhs){
    paste0(lhs, " ~ ", rhs)
  }) %>%
    setNames(names(rhs_list))

  fml_as_chr_list
}
