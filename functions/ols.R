#######################################################################
# OLS Schätzer with HC robust White Standard Errors ###################
# Author: Florian Schulz ##############################################
#######################################################################

ols = function(formula, dataframe, intercept = TRUE, var_method = "standard") {
  # Transform formula into data frame (uses stats package)
  data = model.frame(formula, dataframe)
  
  # Define dependent variable
  y = as.matrix(data[, 1])
  # ... and independet variable
  if (intercept) {
    X = as.matrix(cbind(1, data[, -1]))
    colnames(X)[1] = "(const)"
  } else { X = as.matrix(data[, -1])}
  
  # Estimate beta based on standard formula (X'X)^-1 X'y
  out = list(beta = as.numeric(solve(crossprod(X))%*%crossprod(X,y)))
  
  # Calculate fitted values
  out$fitted = as.numeric(crossprod(t(X), out$beta))
  
  # Calculate residuals (res = y - X beta)
  out$res = as.numeric(y - out$fitted)
  
  # Calculate degrees of freedom
  out$n  = nrow(X)
  k      = ncol(X)
  out$df = out$n - k
  
  # Calculate sum of squared residuals
  out$ssr = sum(out$res^2)
  
  # Calculate variance-covariance matrix
  ## Maybe TODO: implement HAC estimators (kind of difficult, since needs weights)
  if (var_method == "standard") {
    out$Vp = 1/out$df * out$ssr * solve(crossprod(X))
  } else if (var_method  == "HC") {
    out$Vp = solve(crossprod(X))%*%(crossprod(X, diag(as.numeric(out$res^2)))%*%X)%*%solve(crossprod(X))
  } else {
    warning("Use 'standard' or 'HC' as var_method, else 'standard' will be used.")
    Vp = 1/out$df * as.numeric(crossprod(out$res)) * solve(crossprod(X))
  }
  
  # Calculate standard errors of estimates
  out$se = sqrt(diag(out$Vp))
  
  # Calculate p-values
  out$p_value = 2*pt(abs(out$beta/out$se), df = out$df, lower.tail = FALSE)
  
  # Calculate sum of squared explained values
  sse = if (intercept) {
    sum((out$fitted - mean(out$fitted))^2)
  } else { sum(out$fitted^2) }
  
  # Calculate R squared
  out$r.squared = sse/(sse + out$ssr)
  
  # Calculate adjusted R squared
  out$adj.r.squared = 1 - (1 - out$r.squared) * ((out$n - as.numeric(intercept))/out$df)
  
  # Calculate F statistic
  ftest.stat = (sse/(k-as.numeric(intercept)))/(out$ssr/out$df)
  out$f.stat = c(test.stat = ftest.stat, num.par = k - as.numeric(intercept), df = out$df,
                 p.value = pf(ftest.stat, k - as.numeric(intercept), out$df, lower.tail = FALSE))
  
  return(out)
}