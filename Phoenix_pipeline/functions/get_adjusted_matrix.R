get_adjusted_matrix = function(effect_removed_df, column_adjusted) {
  column_adjusted = enquo(column_adjusted)
  
  m =
    effect_removed_df |>
    unnest(!!column_adjusted) |>
    # dplyr::filter(analysis == "observed_proportion") |>
    select(.feature, adjusted___Estimate, sample_id) |>
    pivot_wider(names_from = sample_id, values_from = adjusted___Estimate) |>
    tidybulk:::as_matrix(rownames = ".feature") |>
    as("sparseMatrix")  |>
    Matrix::Matrix(sparse = T)
  
  # Cap infinite
  max_rm_infinite =
    m |>
    _[!m |> is.infinite()] |>
    quantile(0.999)
  
  m |>
    _[m > max_rm_infinite] =
    max_rm_infinite
  
  m |>
    _[m < 0] =
    0
  
  return(m)
}