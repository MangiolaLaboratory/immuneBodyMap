fit_to_age_monotonic_changes = function(fit) {
  full_bins <- c(
    "age_decade1",
    "age_decade2",
    "age_decade3",
    "age_decade4",
    "age_decade6",
    "age_decade7",
    "age_decade8",
    "age_decade9",
    "age_decade10"
  )
  
  vars <- brms::variables(fit)            # all parameter names in the model
  
  # a bin is present if its *population* coefficient exists
  present_bins <- full_bins[paste0("b_", full_bins) %in% vars]
  
  ## need at least two bins to make a split
  if (length(present_bins) < 2) {
    return(tibble())                      # empty result → nothing to contrast
  }
  
  ## helper to build one contrast for an arbitrary vector of bins
  build_contrast <- function(k, bins) {
    younger <- bins[1:k]
    older   <- bins[(k + 1):length(bins)]
    glue::glue(
      "({paste(older,   collapse = ' + ')})/{length(older)} - ",
      "({paste(younger, collapse = ' + ')})/{length(younger)} > 0"
    )
  }
  
  
  # 2. loop over split points *within the available bins*
  
  purrr::map_dfr(seq_len(length(present_bins) - 1), function(k) {
    h_txt <- build_contrast(k, present_bins)
    
    # total = population + random
    h_tot <- brms::hypothesis(fit, h_txt, scope = "coef", group = "tissue_groups")$hypothesis %>%
      tibble::as_tibble() %>%
      dplyr::transmute(
        component    = "total",
        tissue       = Group,
        split_after  = present_bins[k],
        younger_bins = paste(present_bins[1:k], collapse = ","),
        older_bins   = paste(present_bins[(k + 1):length(present_bins)], collapse = ","),
        estimate     = Estimate,
        ci_lower     = CI.Lower,
        ci_upper     = CI.Upper,
        post_prob    = Post.Prob
      )
    
    ##  fixed = population-level only
    h_fix <- brms::hypothesis(fit, h_txt, scope = "standard")$hypothesis %>%
      tibble::as_tibble() %>%
      dplyr::transmute(
        component    = "fixed",
        tissue       = "population",
        split_after  = present_bins[k],
        younger_bins = paste(present_bins[1:k], collapse = ","),
        older_bins   = paste(present_bins[(k + 1):length(present_bins)], collapse = ","),
        estimate     = Estimate,
        ci_lower     = CI.Lower,
        ci_upper     = CI.Upper,
        post_prob    = Post.Prob
      )
    
    dplyr::bind_rows(h_tot, h_fix)
  })
}