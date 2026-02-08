custom_meta_aggregate_logOR_new <- function(data, ps_details){
  
  n_studies <- dplyr::n_distinct(data$ma_id)
  n_es <- nrow(data)
  n_info <- tibble::tibble(n_studies = n_studies, n_es = n_es)
  
  data <- data %>%
    dplyr::mutate(
      ma_id   = factor(ma_id),
      ma_e_id = factor(ma_e_id)  # keep original ids for DOI linkage
    )
  
  # You really need at least 2 estimates AND at least 2 ma_id for robust clustering
  if (n_es >= 2 && n_studies >= 2) {
    
    # NEW: overlap-based VCV from observed shared primary studies
    vcv_obj <- build_vcv_from_overlap(
      dat_sub    = data,
      ps_details = ps_details,
      id_col     = "ma_e_id",
      vi_col     = "se"
    )
    
    VCV   <- vcv_obj$V
    data2 <- vcv_obj$dat_vcv  # reordered to match VCV
    
    suppressWarnings(
      mod <- metafor::rma.mv(
        yi = logOR,
        V  = VCV,
        random = list(~1 | ma_id, ~1 | ma_e_id),
        test = "t",
        sparse = TRUE,
        data = data2,
        control = list(iter.max = 1000, rel.tol = 1e-8)
      )
    )
    
    mod_rob <- metafor::robust(
      mod,
      cluster = ma_id,
      adjust = TRUE,
      clubSandwich = TRUE
    )
    
    tau2_raw <- sum(mod$sigma2)
    sigma2_v <- sum(1/mod$vi) * (mod$k - 1) / (sum(1/mod$vi)^2 - sum((1/mod$vi)^2))
    I2_raw <- 100 * (tau2_raw / (tau2_raw + sigma2_v))
    CV_raw <- sqrt(tau2_raw) / abs(mod$b[1])
    
    res <- tibble::tibble(
      estimate = mod_rob$beta[[1]],
      ci.lb = mod_rob$ci.lb,
      ci.ub = mod_rob$ci.ub,
      pval = mod_rob$pval,
      tval = mod_rob$zval,
      df = mod_rob$ddf,
      tau2 = tau2_raw,
      I2_logOR = I2_raw,
      CV = CV_raw
    )
    
  } else if (n_es == 1) {
    
    est <- data$logOR[1]
    se_val <- data$se[1]
    
    res <- tibble::tibble(
      estimate = est,
      ci.lb = est - qnorm(0.975) * se_val,
      ci.ub = est + qnorm(0.975) * se_val,
      pval = 2 * (1 - pnorm(abs(est / se_val))),
      tval = NA_real_,
      df = NA_real_,
      tau2 = NA_real_,
      I2_logOR = NA_real_,
      CV = NA_real_
    )
    
  } else {
    
    # If you land here, you typically have >=2 estimates but <2 ma_id.
    # Robust clustering by ma_id isn't meaningful, so return NA or a simple FE model.
    # Here I keep your original approach but using rma() as a fallback.
    suppressWarnings(
      mod <- metafor::rma(
        yi = logOR,
        vi = se^2,
        method = "EE",
        data = data
      )
    )
    
    res <- tibble::tibble(
      estimate = mod$beta[[1]],
      ci.lb = mod$ci.lb,
      ci.ub = mod$ci.ub,
      pval = mod$pval,
      tval = NA_real_,
      df = NA_real_,
      tau2 = NA_real_,
      I2_logOR = NA_real_,
      CV = NA_real_
    )
  }
  
  summary_data <- dplyr::bind_cols(res, n_info) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), \(x) round(x, 5))) %>%
    dplyr::select(estimate, ci.lb, ci.ub, pval, tval, df, tau2, I2_logOR, n_studies, n_es, CV)
  
  return(summary_data)
}
