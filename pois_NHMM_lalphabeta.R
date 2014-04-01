pois_NHMM_lalphabeta = function(R, m, D, state_probs, lambda)
{
    #if (is.null(delta)) 
    #    delta = solve(t(diag(m) - gamma + 1), rep(1,m))
    n = length(R)
    lalpha = lbeta = matrix(NA, m, n)
    allprobs = outer(R, lambda[1:m], dpois)

    foo = state_probs$initial * allprobs[1,]    # foo is alpha_1(i) for i = 1..m
    sumfoo = sum(foo)
    lscale = log(sumfoo)
    foo = foo / sumfoo
    lalpha[,1] = log(foo) + lscale
    
    for (k in 2:n)
    {
        foo = foo %*% state_probs$transition[[k-1]] * allprobs[k,]
        sumfoo = sum(foo)
        lscale = lscale + log(sumfoo)
        foo = foo / sumfoo
        lalpha[,k] = log(foo) + lscale
    }
    
    lbeta[,n] = rep(0,m)
    foo = rep(1/m, m)
    lscale = log(m)
    for (k in (n-1):1)
    {
        foo = state_probs$transition[[k]] %*% (allprobs[k+1,] * foo)
        lbeta[,k] = log(foo) + lscale
        sumfoo = sum(foo)
        foo = foo / sumfoo
        lscale = lscale + log(sumfoo)
    }
    
    return (list(la = lalpha, lb = lbeta))
}