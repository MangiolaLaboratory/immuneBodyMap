fit_to_sex_differ_at_decades = function(fit, decades = 5:7) {
  
  full_bins <- paste0('age_decade', decades)
  
  vars <- brms::variables(fit)            # all parameter names in the model
  
  # a bin is present if its *population* coefficient exists
  present_bins <- full_bins[paste0("b_", full_bins) %in% vars]
  
  h_txt = present_bins %>% paste0(':sexmale')
  
  if (5 %in% decades & !"age_decade5:sexmale" %in% h_txt) {
    h_txt = h_txt %>% append('sexmale')      
  }
  
  h_txt = h_txt %>% paste0(' > 0')
  
  h_tot <- brms::hypothesis(fit, h_txt, scope = "coef", group = "tissue_groups")$hypothesis %>%
    tibble::as_tibble() %>%
    mutate(
      age_decade = str_extract(Hypothesis, "age_decade\\d+"),
      age_decade = if_else(is.na(age_decade), "age_decade5", age_decade)
    ) %>% 
    dplyr::transmute(
      component    = "total",
      tissue       = Group,
      age_decade  = age_decade,
      estimate     = Estimate,
      ci_lower     = CI.Lower,
      ci_upper     = CI.Upper,
      post_prob    = Post.Prob
    )
  
  h_fix <- brms::hypothesis(fit, h_txt, scope = "standard")$hypothesis %>%
    tibble::as_tibble() %>%
    mutate(
      age_decade = str_extract(Hypothesis, "age_decade\\d+"),
      age_decade = if_else(is.na(age_decade), "age_decade5", age_decade)
    ) %>% 
    dplyr::transmute(
      component    = "fixed",
      tissue       = "population",
      age_decade  = age_decade,
      estimate     = Estimate,
      ci_lower     = CI.Lower,
      ci_upper     = CI.Upper,
      post_prob    = Post.Prob
    )
  
  dplyr::bind_rows(h_tot, h_fix)
  
}