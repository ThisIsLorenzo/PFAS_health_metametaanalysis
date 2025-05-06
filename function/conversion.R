# Define a conversion function with an optional baseline risk parameter
convert_to_logOR <- function(effect_measure, point_estimate, ma_l_ci, ma_u_ci, regression_type = NA, baseline_risk = NA) {
  # Ensure the effect measure is treated as character
  effect_measure <- as.character(effect_measure)
  # Ensure point estimate, lower CI, and upper CI are numeric
  point_estimate <- suppressWarnings(as.numeric(as.character(point_estimate)))
  ma_l_ci <- suppressWarnings(as.numeric(as.character(ma_l_ci)))
  ma_u_ci <- suppressWarnings(as.numeric(as.character(ma_u_ci)))
  
  # Conversion for Odds ratio and log Odds ratio
  if (effect_measure == "Odds ratio") {
    # Convert OR to logOR
    return(c(log(point_estimate), log(ma_l_ci), log(ma_u_ci)))
    
  } else if (effect_measure == "log Odds ratio") {
    # Already in logOR
    return(c(point_estimate, ma_l_ci, ma_u_ci))
    
    # Conversion for SMD (standardized mean difference)
  } else if (effect_measure == "SMD (standardized mean difference)") {
    # Following Borenstein et al. (2021): logOR = SMD * pi/sqrt(3)
    logOR_point_estimate <- point_estimate * (pi / sqrt(3))
    
    # Convert lower and upper confidence intervals using the same logic
    logOR_l_ci <- ma_l_ci * (pi / sqrt(3))
    logOR_u_ci <- ma_u_ci * (pi / sqrt(3))
    
    return(c(logOR_point_estimate, logOR_l_ci, logOR_u_ci))
    
    # Conversion for Relative risk / Risk ratio with optional baseline risk
  } else if (effect_measure == "Risk ratio") {
    if (!is.na(baseline_risk)) {
      # Using the conversion: OR = RR * (1-p0) / (1-RR*p0)
      OR <- point_estimate * (1 - baseline_risk) / (1 - point_estimate * baseline_risk)
      OR_l_ci <- ma_l_ci * (1 - baseline_risk) / (1 - ma_l_ci * baseline_risk)
      OR_u_ci <- ma_u_ci * (1 - baseline_risk) / (1 - ma_u_ci * baseline_risk)
      return(c(log(OR), log(OR_l_ci), log(OR_u_ci)))
    } else {
      # If no baseline risk is provided, fallback to the approximation:
      return(c(log(point_estimate), log(ma_l_ci), log(ma_u_ci)))
    }
    
    # Conversion for Fisher's z (often used for correlations)
  } else if (grepl("fishers z", effect_measure, ignore.case = TRUE)) {
    # Convert z to correlation, then to logOR.
    # r = tanh(z)
    # Approximation: logOR = r * pi/sqrt(3)
    r <- tanh(point_estimate)
    return(r * (pi / sqrt(3)))
    
    # Conversion for regression coefficients assumed from logistic regression:
  } else if (effect_measure %in% c("Beta (regression coefficient) value")) {
    # Check if regression_type is 'logistic'
    if (!is.na(regression_type) && regression_type == "logistic") {
      return(c(point_estimate, ma_l_ci, ma_u_ci))
    } else {
      return(rep(NA, 3)) # Not logistic regression, cannot convert
    }
  } else {
    warning("Effect measure conversion not implemented for: ", effect_measure)
    return(rep(NA, 3)) # Return NA for all three
  }
}