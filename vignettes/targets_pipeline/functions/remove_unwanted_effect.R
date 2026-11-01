#-----------------------#
# Functions
#-----------------------#

#' Remove Unwanted Effects from a brmsfit Model
#'
#' This function calculates posterior residuals from a \code{brmsfit} model and combines them with
#' factor-specific fitted values (potentially excluding random effects or other parts of the model),
#' thereby producing adjusted outcomes that highlight the contribution of a specified factor or subset
#' of model terms.
#'
#' @param fit A \code{brmsfit} object, resulting from a model fitted by \code{\link[brms]{brm}}.
#' @param newdata A data frame or list containing new data. Passed to \code{\link[brms]{fitted}}
#'   to obtain factor-specific fitted values at specified covariate levels.
#' @param robust A logical value indicating whether to use robust (median-based) summaries rather
#'   than means. Defaults to \code{FALSE}.
#' @param correct_by_offset A logical value indicating whether to divide the residuals by
#'   \code{exp(offset)} (from \code{fit$data$offset}). Defaults to \code{TRUE}.
#' @param re_formula A formula specifying which random effects (if any) to include when generating
#'   fitted values. Defaults to \code{~0}, which removes random effects and thus isolates the
#'   contribution of fixed effects in the new data.
#'
#' @return A \code{tibble} containing posterior summaries of:
#'   \itemize{
#'     \item Adjusted outcomes (prefix: \code{adjusted___}): The combined values of the specified
#'     factor's fitted counts and the residuals.
#'     \item Residuals (prefix: \code{residuals___}): The model's posterior residuals, possibly
#'     normalised by the offset.
#'     \item Fitted values for the factor (prefix: \code{fitted___}): The model's fitted values based
#'     on the \code{re_formula} and provided \code{newdata}.
#'   }
#'
#' @details
#' The function proceeds as follows:
#' \enumerate{
#'   \item Extracts posterior residuals (via \code{\link[brms]{residuals}}).
#'   \item (Optionally) divides these residuals by the exponential of the offset, if \code{correct_by_offset = TRUE}.
#'   \item Obtains new fitted values from the model (via \code{\link[brms]{fitted}}), usually excluding random effects
#'         by specifying \code{re_formula = ~0}.
#'   \item Adds these residuals to the factor-specific fitted values to obtain adjusted outcomes
#'         that highlight the contribution of the factor of interest.
#'   \item Summarises all these draws (residuals, fitted values, adjusted outcomes) and returns them
#'         in a single \code{tibble}.
#' }
#'
#' This method is particularly useful for examining how a factor or other subset of the model
#' affects the outcome when other model components (e.g., random intercepts) are removed.
#' It can assist in visualising or quantifying the partial contribution of certain terms.
#'
#' @examples
#' \dontrun{
#' # Suppose 'fit' is a brmsfit model object predicting a 'counts' outcome
#' # We create a new data frame 'some_data' for which we want partial predictions
#' adjusted_results <- remove_unwanted_effect(
#'   fit,
#'   newdata = some_data,
#'   robust = TRUE,
#'   correct_by_offset = TRUE,
#'   re_formula = ~0
#' )
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_cols
#' @importFrom tibble as_tibble
#' @importFrom brms posterior_summary fitted residuals
#'
#' @export
remove_unwanted_effect = function(fit,
                                  newdata,
                                  robust = FALSE,
                                  correct_by_offset = T,
                                  re_formula = ~ 0) {
  # Calculate residuals: observed counts minus fitted values, normalised by exp(offset)
  # This places residuals on a consistent scale, making them addable to adjusted predictions later.
  fitted_residuals =   fit |> residuals(robust = robust, summary = FALSE)
  
  # Correct by offset
  if (correct_by_offset)
    fitted_residuals = fitted_residuals |>
      sweep(2, fit$data$offset |> exp(), FUN = "/")
  
  # Extract fitted values for the specified factor only, removing random effects by setting re_formula = ~0
  # 'resp = factor' focuses on the selected response variable (factor)
  fitted_values_ethnicity <- fitted(
    fit,
    newdata = newdata,
    re_formula = re_formula,
    summary = FALSE,
    offset = 0
  )
  
  # Adjusted counts are obtained by adding the factor-specific fitted values and the normalised residuals
  adjusted_counts = fitted_values_ethnicity + fitted_residuals
  
  # Summarise residuals into a tibble, prefixed to denote their source
  fitted_residuals_tbl =
    fitted_residuals |>
    posterior_summary(robust = robust) |>
    as_tibble()
  fitted_residuals_tbl |> colnames() = paste0("residuals___", fitted_residuals_tbl |> colnames())
  
  # Summarise the factor-only fitted values into a tibble, prefixed accordingly
  fitted_values_ethnicity_tbl =
    fitted_values_ethnicity |>
    posterior_summary(robust = robust) |>
    as_tibble()
  fitted_values_ethnicity_tbl |> colnames() = paste0("fitted___", fitted_values_ethnicity_tbl |> colnames())
  
  # Summarise the adjusted counts (factor + residuals) into a tibble, prefixed for clarity
  adjusted_counts_tbl =
    adjusted_counts |>
    posterior_summary(robust = robust) |>
    as_tibble()
  adjusted_counts_tbl |> colnames() = paste0("adjusted___", adjusted_counts_tbl |> colnames())
  
  # Combine all three resulting tables into one tibble
  adjusted_counts_tbl |>
    bind_cols(fitted_residuals_tbl) |>
    bind_cols(fitted_values_ethnicity_tbl)
}