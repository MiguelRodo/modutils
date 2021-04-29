#' @title Calculate Wald statistics
#'
#' @param fit object of class 'glmerMod' or 'list'. Methods are defined for each of these
#' to then extract the estimates and the var-cov matrix before calculating the Wald statistics.
#' If a list, then must have elements with names \code{est} and \code{vcov}, corresponding
#' to the numeric vector of estimates and the variance-covariance matrix of
#' @param var list of character vector(s). Each element is a vector of name(s) of variables for which
#' the Wald test is required. If \code{"0"},
#' then all variables are automatically included. If \code{"-1"}, then all variables except
#' intercept are included. No default.
#' @param match_condn 'exact', 'start' or 'any'. If \code{'exact'}, then
#' variables in \code{var} match coefficient names in \code{fit} only if they match exactly (e.g.
#' \code{'abc'} matches \code{'abc'} but neither \code{'a'} nor \code{'c'} do). If \code{ 'start'},
#' then they match only if they correspond to the start of coefficient names in \code{fit} (e.g.
#' \code{'abc'} and \code{'a'} match \code{'abc'} \code{'c'} does not). If \code{'any'},
#' then they match as long as they appear at some point in the coefficient name in \code{fit} (e.g.
#' \code{'abc'}, \code{'a'} and \code{'c'} all match \code{'abc'}). Default is \code{'start'}.
#' @param match_print 'all', 'inexact' or 'none'. If \code{'all'}, then
#' all matched coefficient names in \code{fit} are printed. If \code{'inexact'}, then
#' only inexact matches are printed. If \code{'none'}, then none are printed. Default is \code{'all'}.
#' @param intercept_nm character. Name for intercept coefficient. Default is \code{"(Intercept)"}. Only
#' used if \code{class(fit) = 'list'}.
#'
#' @return A tibble with the following elements:
#' \tabular{ll}{
#' test \tab 'wald' \cr
#' stat \tab test statistic \cr
#' df \tab degrees of freedom \cr
#' p \tab p-value \cr
#' }
#'
#' @export
test_wald <- function(fit, var, match_condn = 'start', match_print = 'inexact',
                           intercept_nm = "(Intercept)"){

  # preparation
  # ==================

  # checks for availability
  if(missing(var)) stop("var missing in test_wald")
  if(missing(fit)) stop("fit missing in test_wald")

  # ensure var is a list
  if(!is.list(var)){
    warning("var is coerced to a list in test_wald")
    var <- list(var)
  }
  # ensure 0 and 1 are characters
  var <- purrr::map(var, function(x){
    if(is.character(x)) return(x)
    if(x == 0){
      warning("0 as numeric coerced to 0 as character in var in test_wald")
      return("0")
    } else if(x == -1){
      warning("-1 as numeric coerced to 0 as character in var in test_wald")
      return("-1")
    }
    stop("input to var numeric but neither 0 nor -1 and so is not meaningful in var param in test_wald")
    })

  # extract estimates and vcov mat
  # -----------------
  fit_list <- get_est_and_vcov(fit = fit)

  # calculate Wald statistics
  # ==================

  purrr::map_df(var, function(var){
    .test_wald(est = fit_list$est, vcov = fit_list$vcov,
                  var = var, match_condn = match_condn,
                  match_print = match_print, intercept_nm = intercept_nm)})
}

get_est_and_vcov <- function(fit){
  UseMethod("get_est_and_vcov")
}
get_est_and_vcov.glmerMod <- function(fit){
  est_vec <- coefficients(summary(fit))[,"Estimate"]
  vcov_mat <- as.matrix(vcov(fit))
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

#' @title Calculate wald statistic from coef estimates and vcov-mat
.test_wald <- function(est, vcov, var, match_condn, match_print, intercept_nm){

  # preparation
  # ==================

  # checks for availability
  if(missing(var)) paste0("var missing in test_wald")
  if(missing(est)) paste0("est missing in test_wald")
  if(missing(vcov)) paste0("vcov missing in test_wald")

  # set var to all names if given as NULL
  # -------------------

  if(identical(var, "0")) var <- names(est)
  if(identical(var, "-1")){
    var <- names(est)
    var <- var[!var == intercept_nm]

  }

  # get indices of vcoefficients
  # ------------------------
  coef_ind_vec <- .get_var_ind(est = est, vcov = vcov,
                               var = var, match_condn = match_condn,
                               match_print = match_print)

  # get wald stats
  # ==================
  .test_wald_ind(est = est, vcov = vcov, ind = coef_ind_vec, var = var)
}

#' @title Get indices of variables to include in wald test
#'
#' @return A numeric vector o
.get_var_ind <- function(est, vcov, var, match_condn, match_print){

  # checks
  # ----------------------

  if(missing(est) || missing(vcov) || missing(var) ||
     missing(match_condn) ||  missing(match_print)) {
    stop("At least one of est, vcov, var, match_condn and match_print missing in .get_var_ind")
  }
  if(!match_print %in% c("all", "inexact", "none")){
    stop("match_print has value ", paste0(class(match_print), collapse = " "),
         " when it should a value of one of 'all', 'inexact' or 'none'")
  }
  if(!match_condn %in% c("exact", "start", "any")){
    stop("match_condn has argument ", paste0(match_condn, collapse = " "),
         " when it should be one of 'exact', 'start' or 'any' in .get_var_ind")
  }
  if(!is.numeric(est)){
    stop("est has class ", paste0(class(est), collapse = " "),
         " when it should have class numeric in .get_var_ind")
  }
  if(!is.matrix(vcov)){
    stop("vcov has class ", paste0(class(vcov), collapse = " "),
         " when it should have class matrix in .get_var_ind")
  }
  if(!is.character(var)){
    stop("var has class ", paste0(class(var), collapse = " "),
         " when it should have class character in .get_var_ind")
  }

  if(is.null(names(est))) stop("est must be a named vector in .get_var_ind")
  if(is.null(rownames(vcov))) stop("vcov must be matrix with rownames in .get_var_ind")

  # check that vcov and est names are ordered the same
  vcov_name_vec <- rownames(vcov)
  vcov_name_vec_comp <- .bracket_non_an_chr(vcov_name_vec) # double-bracket non-an chr for comparisons
  est_name_vec <- names(est)
  if(!identical(est_name_vec, vcov_name_vec)){
    stop(paste0("The coefficient estimates' names do not match the names in the var-cov matrix."))
  }

  # check that all elements specified by var are in the rownames
  var_missing_vec <- NULL
  for(i in seq_along(var)){
    var_curr <- var[i]
    var_comp <- .bracket_non_an_chr(var_curr)
    if(!any(stringr::str_detect(vcov_name_vec,
                                var_comp))){
      var_missing_vec <- c(var_missing_vec, var_curr)
    }
  }
  if(!is.null(var_missing_vec)){
    stop(paste0("The variables ",
                paste0(var_missing_vec, collapse = " "),
                " are not found in var parameter in .get_var_ind"))
  }

  # get indices to be tested
  # ----------------

  # check for each element in var
  purrr::map(var, function(var_curr){
    var_curr <- .bracket_non_an_chr(var_curr)
    # match against each coefficient name
    ind_vec_lgl <- purrr::map_lgl(vcov_name_vec_comp, function(nm){
      # match at any point
      if(match_condn == 'any') return(stringr::str_detect(nm, var_curr))
      # exact match
      if(match_condn == 'exact') return(identical(nm, var_curr))
      # match from start (checked up above that match_condn must now have value 'start')
      nm_sub <- stringr::str_sub(nm, end = stringr::str_length(var_curr))
      identical(nm_sub, var_curr)
    })
    # return matching indices
    which(ind_vec_lgl)
  }) %>%
    unlist() %>%
    unique()

}

#'@title Place double square brackets around non-alphanumeric characters
#'
#'@description Useful for string comparisons
#'
#'@param x character vector.
#'
#'@return Character, with non-alphanumeric characters now being bracketed by double square brackets.
#'
#'@examples
#'
.bracket_non_an_chr <- function(x){
  purrr::map_chr(x, function(x_ind){
    x_split_vec <- purrr::map_chr(1:stringr::str_length(x_ind),
                                  function(i) stringr::str_sub(x_ind, i, i))
    is_an_vec <- purrr::map_lgl(x_split_vec, function(chr) stringr::str_detect(chr, "\\w"))

    #x_out <- ifelse(!is_an_vec[1], paste0("[[", x_split_vec[1]), x_split_vec[1])
    # if string has length 1
    if(stringr::str_length(x_ind) == 1){
      x_out <- ifelse(!is_an_vec[1], paste0(x, "[[", x_ind, "]]"), x_ind)
      return(x_out)
    }
    # now we must count the number of times each value is repeated, if it is non-alphanumeric (otherwise it's 1)
    x_split_vec_rem <- x_split_vec
    k <- length(x_split_vec_rem)
    x_split_vec_rep <- ""
    i <- 1
    info_tbl <- tibble::tibble(chr = NULL, is_an = NULL, rep = NULL, add = NULL)
    while(k > 0){
      chr_curr <- x_split_vec_rem[1]
      is_an <- stringr::str_detect(chr_curr, "\\w")
      # reduce vector to draw from
      x_split_vec_rem <- x_split_vec_rem[-1]
      # increase indent, so that it will eventually stop
      k <- length(x_split_vec_rem)
      rep_curr <- 1
      if(is_an){
        info_tbl %<>%
          dplyr::bind_rows(tibble::tibble(chr = chr_curr, is_an = FALSE, rep = 1, add = chr_curr))
      } else if(k == 0){
        # if there are none to follow
        info_tbl %<>%
          dplyr::bind_rows(tibble::tibble(chr = chr_curr, is_an = TRUE,
                                          rep = 1, add = paste0("[[", chr_curr, "]]")))
        } else{

        # wrap in double brackets for comparison
        chr_next_comp <- paste0("[[", x_split_vec_rem[1], "]]")
        # if it's not a match to the next one, then just bind on
        if(!stringr::str_detect(chr_curr, chr_next_comp)){
          info_tbl %<>%
            dplyr::bind_rows(tibble::tibble(chr = chr_curr, is_an = TRUE, rep = 1, add =
                                              paste0("[[", chr_curr, "]]")))
        } else{
          # if a match, continue increasing rep until it isn't or we reach the end
          while(stringr::str_detect(chr_curr, chr_next_comp)){
            # increase number of reps as we have a rep
            rep_curr <- rep_curr + 1
            # remove latest "first" elem from remaining values
            x_split_vec_rem <- x_split_vec_rem[-1]
            # count number of remaining values
            k <- length(x_split_vec_rem)
            # if there are no remaining values, set next chr to an an character so we force no match
            chr_next_comp <- ifelse(k == 0, "abc", paste0("[[", x_split_vec_rem[1], "]]"))
          }
          info_tbl %<>%
            dplyr::bind_rows(tibble::tibble(chr = chr_curr, is_an = TRUE, rep = rep_curr,
                                            add = paste0("[[", chr_curr, "]]{", rep_curr, "}")))
        }
      }
    }
    paste0(info_tbl$add, collapse = "")
  })
}

#' @title Calculate Wald statistic given summary statistics
#'
#' @param est numeric vector. Estimates of parameters.
#' @param vcov numeric matrix. Variance-covarance matrix of parameters.
#' @param ind logical or numeric vector. Indicates which are the parameters
#' in \code{est} that we need to test against the null hypothesis of being zero.

#'
#' @return A tibble with the following elements:
#' \tabular{ll}{
#' stat \tab Wald test statistic \cr
#' df \tab Degrees of freedom \cr
#' p \tab P-value \cr
#' }
#'
#' @examples
#' est_vec <- c(0.2,5,1)
#' vcov_mat <- as.matrix(rbind(c(1, 0.3, 0.1), c(0.3, 4, 0.4), c(0.1, 0.4, 2.2)))
#' ind_vec_lgl <- c(TRUE, FALSE, TRUE)
#' ind_vec_num <- c(1, 3)
#' # using logical vector to indicate which parameters are to be tested
#' .test_wald(est = est_vec, vcov = vcov_mat, ind = ind_vec_lgl)
#' # using numeric vector to indicate which parameters are to be tested
#' .test_wald(est = est_vec, vcov = vcov_mat, ind = ind_vec_num)
.test_wald_ind <- function(est, vcov, ind, var){

  # check entries
  # --------------------

  if(missing(ind) || missing(vcov) || missing(est)){
    stop("at least one of est, vcov and ind params in .test_wald missing")
  }
  if(missing(var)) var <- NA
  if(is.null(var)) var <- NA


  if(!is.matrix(vcov)){
    stop(paste0("vcov in .test_wald has class ",
                paste0(class(vcov), " "),
                " when it should have class numeric"))
  }
  if(!is.numeric(vcov[1,1])) stop(paste0("vcov entries in .test_wald have class ",
                                         paste0(class(vcov[1,1]), " "),
                                         " when it should have class numeric"))
  if(!is.numeric(est[1])) stop(paste0("est entries in .test_wald have class ",
                                      paste0(class(est[1]), " "),
                                      " when it should have class numeric"))

  # get degrees of freedom and check ind class
  df <- switch(class(ind),
               "logical" = sum(ind),
               "integer" = ,
               "numeric" = length(ind),
               stop("params ind in .test_wald of class ", paste0(class(ind, collapse = " ")),
                    " when it should be logical, integer or numeric"))


  # calculate wald test stats
  # ------------------------

  vcov_mat_ind <- vcov[ind,ind]
  est_vec_ind <- est[ind]
  wald_stat <- matrix(est_vec_ind, nrow = 1) %*%
    solve(vcov_mat_ind) %*%
    matrix(est_vec_ind, ncol = 1)
  wald_stat <- as.numeric(wald_stat)
  wald_p <- pchisq(wald_stat, df = df, lower.tail = FALSE)

  tibble::tibble(var = paste0(var, collapse = "; "),
                 "test" = 'wald', "stat" = wald_stat, "df" = length(ind),  "p" = wald_p)
}
