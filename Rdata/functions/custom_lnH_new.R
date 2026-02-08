custom_meta_aggregate_lnH_new <- function(data, ps_details){
  
  n_studies_lnH <- n_distinct(data$ma_id)
  n_es_lnH <- nrow(data)
  n_info <- tibble(n_studies_lnH = n_studies_lnH, n_es_lnH = n_es_lnH)
  
  data <- data %>%
    mutate(
      ma_id   = factor(ma_id),
      ma_e_id = factor(ma_e_id)   # keep original ids for DOI linkage
    )
  
  # Need at least 2 estimates for rma.mv to be meaningful
  if (n_es_lnH >= 2 && n_studies_lnH >= 2) {
    
    # --- NEW: overlap-based VCV on this subset ---
    vcv_obj <- build_vcv_from_overlap(
      dat_sub    = data,
      ps_details = ps_details,
      id_col     = "ma_e_id",
      vi_col     = "SElnH"
    )
    
    VCV_H2 <- vcv_obj$V
    data2  <- vcv_obj$dat_vcv
    
    suppressWarnings(
      mod <- rma.mv(
        yi = lnH,
        V  = VCV_H2,
        random = list(~1 | ma_id, ~1 | ma_e_id),
        test = "t",
        sparse = TRUE,
        data = data2,
        control = list(iter.max = 1000, rel.tol = 1e-8)
      )
    )
    
    mod_rob <- robust(
      mod,
      cluster = ma_id,
      adjust = TRUE,
      clubSandwich = TRUE
    )
    
    tau2_raw <- sum(mod$sigma2)
    sigma2_v <- sum(1/mod$vi) * (mod$k - 1)/(sum(1/mod$vi)^2 - sum((1/mod$vi)^2))
    I2_raw <- 100 * (sum(mod$sigma2)/(sum(mod$sigma2) + sigma2_v))
    CV_raw <- sqrt(tau2_raw) / abs(mod$b[1])
    
    res <- tibble(
      lnH = mod_rob$beta[[1]],
      ci.lb_lnH = mod_rob$ci.lb,
      ci.ub_lnH = mod_rob$ci.ub,
      pval_lnH = mod_rob$pval,
      tval_lnH = mod_rob$zval,
      df_lnH = mod_rob$ddf,
      tau2_lnH = tau2_raw,
      I2_lnH = I2_raw,
      CV_lnH = CV_raw
    )
    
  } else if (n_es_lnH == 1) {
    
    # NOTE: your existing single-estimate branch looks inconsistent
    # (it uses logOR/se, but we are in lnH/SElnH land).
    # A consistent single-estimate lnH branch is:
    
    est <- data$lnH[1]
    se_val <- data$SElnH[1]
    
    res <- tibble(
      lnH = est,
      ci.lb_lnH = est - qnorm(0.975) * se_val,
      ci.ub_lnH = est + qnorm(0.975) * se_val,
      pval_lnH = 2 * (1 - pnorm(abs(est / se_val))),
      tval_lnH = NA_real_,
      df_lnH = NA_real_,
      tau2_lnH = NA_real_,
      I2_lnH = NA_real_,
      CV_lnH = NA_real_
    )
    
  } else {
    
    # If there are 0 rows or only 1 study, you can return NA summary
    res <- tibble(
      lnH = NA_real_,
      ci.lb_lnH = NA_real_,
      ci.ub_lnH = NA_real_,
      pval_lnH = NA_real_,
      tval_lnH = NA_real_,
      df_lnH = NA_real_,
      tau2_lnH = NA_real_,
      I2_lnH = NA_real_,
      CV_lnH = NA_real_
    )
  }
  
  bind_cols(res, n_info) %>%
    mutate(across(where(is.numeric), \(x) round(x, 5))) %>%
    select(lnH, ci.lb_lnH, ci.ub_lnH, pval_lnH, tval_lnH, df_lnH,
           tau2_lnH, I2_lnH, n_studies_lnH, n_es_lnH, CV_lnH)
}
