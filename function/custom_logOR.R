custom_meta_aggregate_logOR <- function(data, rho = 0.5){
  
  n_studies <- length(unique(data$ma_id))
  n_es <- nrow(data)
  n_info <- tibble(n_studies = n_studies, n_es = n_es)
  data$ma_id <- as.factor(data$ma_id)
  
  if(n_studies >= 2){
    data$ma_e_id <- as.factor(1:nrow(data))
    
    VCV <- impute_covariance_matrix(vi = data$se^2,  
                                    cluster = data$ma_id,
                                    r = rho)
    
    suppressWarnings(mod <- rma.mv(yi = logOR,  
                                   V = VCV,
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
    
  } else if(n_es == 1){
    
    est <- data$logOR
    se_val <- sqrt(data$se^2)
    
    res <- tibble(
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
    
    suppressWarnings(mod <- rma(yi = logOR,
                                vi = se^2,
                                method = "EE",
                                data = data)) 
    
    res <- tibble(
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
  
  summary_data <- bind_cols(res, n_info) %>%
    mutate(across(where(is.numeric), \(x) round(x, 5))) %>%
    select(estimate, ci.lb, ci.ub, pval, tval, df, tau2, I2_logOR, n_studies, n_es, CV)
  
  return(summary_data)
}
