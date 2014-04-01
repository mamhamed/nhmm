pois_NHMM_forecast = function(R, X, m, D, lambda, theta, xrange = NULL, H = 1, ...)
{
    if (is.null(xrange))
        xrange = qpois(0.001, min(lambda)) : qpois(0.999, max(lambda))
    
    n = length(R)
    allprobs = outer(R, lambda, dpois)
    allprobs = ifelse (!is.na(allprobs), allprobs, 1)
    state_probs = compute_state_probs(m, n+H, X, theta)
    
    foo = state_probs$initial * allprobs[1,]
    sumfoo = sum(foo)
    lscale = log(sumfoo)
    foo = foo / sumfoo
    
    for (k in 2:n)
    {
        foo = foo %*% state_probs$transition[[k-1]] * allprobs[k,]
        sumfoo = sum(foo)
        lscale = lscale + log(sumfoo)
        foo = foo / sumfoo
    }
    
    xi = matrix(NA, nrow = m, ncol = H)
    for (k in 1:H)
    {
        foo = foo %*% state_probs$transition[[n+k-1]]
        xi[,k] = foo
    }
    
    allprobs = outer(xrange, lambda, dpois)
    fdists = allprobs %*% xi[, 1:H]
    
    list(xrange = xrange, fdists = fdists)
}