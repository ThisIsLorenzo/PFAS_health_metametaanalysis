custom_meta_aggregate_lnH <- function(data, rho = 0.5){
  
  n_studies_lnH <- length(unique(data$ma_id))
  n_es_lnH <- nrow(data)
  n_info <- tibble(n_studies_lnH = n_studies_lnH, n_es_lnH = n_es_lnH)
  data$ma_id <- as.factor(data$ma_id)
  
  if(n_studies_lnH >= 2){
    data$ma_e_id <- as.factor(1:nrow(data))
    
    VCV_H2 <- impute_covariance_matrix(vi = data$SElnH^2,  
                                    cluster = data$ma_id,
                                    r = rho)
    
    suppressWarnings(mod <- rma.mv(yi = lnH,  
                                   V = VCV_H2,
                                   random = list(~1 | ma_id, ~1 | ma_e_id),
                                   test = "t",
                                   sparse = TRUE,
                                   data = data,
                                   control = list(iter.max = 1000, rel.tol = 1e-8)))
    
    mod_rob <- robust(mod, 
                      cluster = ma_id, 
                      adjust = TRUE, 
                      clubSandwich = TRUE) 
    
    tau2_raw <- sum(mod$sigma2)
    sigma2_v <- sum(1/mod$vi) * (mod$k - 1)/(sum(1/mod$vi)^2 - sum((1/mod$vi)^2))
    I2_raw <- 100 * (sum(mod$sigma2)/(sum(mod$sigma2) + sigma2_v))
    CV_raw <- sqrt(tau2_raw) / abs(mod$b[1])
    
    res <- tibble(
      lnH = mod_rob$beta[[1]],
      ci.lb_lnH = mod_rob$ci.lb,
      ci.ub_lnH = mod_rob$ci.ub,
      pval_lnH = mod_rob$pval,
      tau2_lnH = tau2_raw,
      I2_lnH = I2_raw,
      CV_lnH = CV_raw
    )
    
  } else if(n_es_lnH == 1){
    
    est <- data$logOR
    se_val <- sqrt(data$se^2)
    
    res <- tibble(
      lnH = est,
      ci.lb_lnH = est - qnorm(0.975) * se_val,
      ci.ub_lnH = est + qnorm(0.975) * se_val,
      pval_lnH = 2 * (1 - pnorm(abs(est / se_val))),
      tau2_lnH = NA_real_,
      I2_lnH = NA_real_,
      CV_lnH = NA_real_
    )
    
  } else {
    
    suppressWarnings(mod <- rma(yi = logOR,
                                vi = se^2,
                                method = "EE",
                                data = data)) 
    
    res <- tibble(
      lnH = mod$beta[[1]],
      ci.lb_lnH = mod$ci.lb,
      ci.ub_lnH = mod$ci.ub,
      pval_lnH = mod$pval,
      tau2_lnH = NA_real_,
      I2_lnH = NA_real_,
      CV_lnH = NA_real_
    )
  } 
  
  summary_data <- bind_cols(res, n_info) %>%
    mutate(across(where(is.numeric), \(x) round(x, 5))) %>%
    select(lnH, ci.lb_lnH, ci.ub_lnH, pval_lnH, tau2_lnH, I2_lnH, n_studies_lnH, n_es_lnH, CV_lnH)
  
  return(summary_data)
}
