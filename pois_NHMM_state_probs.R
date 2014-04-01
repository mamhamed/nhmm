pois_NHMM_state_probs = function(R, X, m, D, lambda, theta)
{
    n = length(R)
    state_probs = compute_state_probs(m, n, X, theta)
    fb = pois_NHMM_lalphabeta(R, m, D, state_probs, lambda)
    la = fb$la
    lb = fb$lb
    c = max(la[,n])
    llk = c + log(sum(exp(la[,n] - c)))
    A = exp(la + lb - llk)
    A = A / apply(A, 2, sum)
    return (A)
    
#     stateprobs = matrix(NA , ncol = n, nrow = m)
#     for (i in 1:n) 
#         stateprobs[,i] = exp(la[,i] + lb[,i] - llk)
#     stateprobs
}