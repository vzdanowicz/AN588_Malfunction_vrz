z.prop.test <- function(p1, n1, p2 = NULL, n2 = NULL, p0, alternative = "two-sided", conf.level = 0.95) 
{
  
  if (is.null(p2)|is.null(n2)) {
    Z <- (p1 - p0) / sqrt(p0 * (1-p0)/n1)
    PVAL <- 1 - pnorm(Z, lower.tail = TRUE) + pnorm(z, lower.tail = FALSE)
    if ((n1 * p0 < 5) | (n1 * (1 - p0) < 5 ))
      return(warning("check distribution, does not appear normal"))
    lower <- p1 - qnorm(0.975) * sqrt(p1 * (1 - p1)/n1)
    upper <- p1 + qnorm(0.975) * sqrt(p1 * (1 - p1)/n1)
    ci <- c(lower, upper)
    one.sample.prop <- list(z.test.statistic = as.numeric(Z), p.value = as.numeric(PVAL), conf.level = ci)
    
    return(one.sample.prop)
  }
  
  else  {
    pstar <- (p1 + p2)/(n1 + n2)
    Z <- (p2 - p1)/sqrt((pstar * (1 - pstar)) * (1/n1) + (1/n2))
    PVAL <- 1 - pnorm(Z, lower.tail = TRUE) + pnorm(Z, lower.tail = FALSE)
    if ((n1 * p0 < 5) | (n1 * (1 - p0) < 5 ) | (n2 * p0 < 5) | (n2 * (1 - p0) < 5)){
      warning("check distribution, does not appear normal")
    }
    alpha = 1 - (conf.level/100)
    crit <- qnorm(1 - alpha/2)  # identify critical values
    test <- p < -crit || p > crit  # boolean test
    two.sample.prop <- list(z.test.statistic = as.numeric(Z), p.value = as.numeric(PVAL), conf.level = conf.level, critical.value = as.numeric(crit), test = test)
    
    return(two.sample.prop)
  }}





p0 <- p0[OK]
if (any((p0 <= 0) | (p0 >= 1))) 
  stop("elements of 'p0' must be in (0,1)")


  PARAMETER <- k - 1
}
else {
  PARAMETER <- k
  names(NVAL) <- names(ESTIMATE)
}
names(PARAMETER) <- "df"
x <- cbind(x, n - x)
E <- cbind(n * p, n * (1 - p))

STATISTIC <- sum((abs(x - E) - YATES)^2/E)
names(STATISTIC) <- "X-squared"
if (alternative == "two.sided") 
  PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
else {
  if (k == 1) 
    z <- sign(ESTIMATE - p) * sqrt(STATISTIC)
  else z <- sign(DELTA) * sqrt(STATISTIC)
  PVAL <- pnorm(z, lower.tail = (alternative == "less"))
}
RVAL <- list(statistic = STATISTIC, parameter = PARAMETER, 
             p.value = as.numeric(PVAL), estimate = ESTIMATE, null.value = NVAL, 
             conf.int = CINT, alternative = alternative, method = METHOD, 
             data.name = DNAME)
class(RVAL) <- "htest"
return(RVAL)
}