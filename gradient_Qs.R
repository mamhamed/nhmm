gradient_Qs = function(m, D, n, A, B, state_probs)
{
    Qs_last_term = list(1:(n-1))
    for (k in 1:(n-1))
    {
        Qs_last_term[[k]] = B[[k]] * (1 - state_probs$transition[[k]])
    }
    
    grad_gamma = matrix(1:(m*m), m, m)
    grad_rho = matrix(rep(0, (m*D)), m, D)
    for (j in 1:m)
    {
        for (i in 1:m)
        {
            grad_gamma[j,i] = sum(sapply(1:(n-1), function(k) Qs_last_term[[k]][j,i]) )
            for (d in 1:D)
                grad_rho[i,d] = grad_rho[i,d] + 
                    sum(sapply(1:(n-1), function(k) Qs_last_term[[k]][j,i] * X[k+1,d]) )
        }
    }
    
    grad_delta = A[,1] * (1 - state_probs$initial)
    
    return (list(delta = grad_delta, gamma = grad_gamma, rho = grad_rho))
}