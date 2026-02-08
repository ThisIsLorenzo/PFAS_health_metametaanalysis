build_vcv_from_overlap <- function(dat_sub, ps_details,
                                   id_col = "ma_e_id",
                                   vi_col = "se") {
  dat_vcv <- dat_sub %>%
    mutate(
      estimate_id = .data[[id_col]],
      vi = (.data[[vi_col]])^2
    ) %>%
    arrange(estimate_id)
  
  ps_long <- ps_details %>%
    transmute(
      estimate_id = ma_e_id,
      primary_doi = primary_doi
    ) %>%
    distinct() %>%
    semi_join(dat_vcv %>% select(estimate_id), by = "estimate_id")
  
  i_idx <- match(ps_long$estimate_id, dat_vcv$estimate_id)
  doi_levels <- sort(unique(ps_long$primary_doi))
  j_idx <- match(ps_long$primary_doi, doi_levels)
  
  A <- Matrix::sparseMatrix(
    i = i_idx,
    j = j_idx,
    x = 1,
    dims = c(nrow(dat_vcv), length(doi_levels))
  )
  
  shared <- Matrix::tcrossprod(A)
  n_i <- Matrix::rowSums(A)
  
  shared_t <- summary(shared)
  denom <- sqrt(n_i[shared_t$i] * n_i[shared_t$j])
  rho_x <- ifelse(denom > 0, shared_t$x / denom, 0)
  
  rho <- Matrix::sparseMatrix(
    i = shared_t$i,
    j = shared_t$j,
    x = rho_x,
    dims = dim(shared)
  )
  
  diag(rho) <- ifelse(n_i > 0, 1, 0)
  
  vi <- dat_vcv$vi
  rho_t <- summary(rho)
  
  V <- Matrix::sparseMatrix(
    i = rho_t$i,
    j = rho_t$j,
    x = rho_t$x * sqrt(vi[rho_t$i] * vi[rho_t$j]),
    dims = dim(rho)
  )
  
  diag(V) <- vi
  V <- Matrix::forceSymmetric(V, uplo = "L")
  
  list(V = V, dat_vcv = dat_vcv)
}
