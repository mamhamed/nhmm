compute_dQs_dnu = function(m, X, n, A, B, state_probs, theta)        
{
    # calculate first derivative d_Qs / d_nu
    Qs_last_term = list(1:(n-1))
    for (k in 1:(n-1))
    {
        Qs_last_term[[k]] = B[[k]] * (theta$gamma + 
                        matrix(rep(X[k+1,] %*% t(theta$rho), m), m, m, byrow = TRUE)) *
                        (1 - state_probs$transition[[k]])
    }
    
    Qs_first_term = A[,1] * (theta$delta + X[1,] %*% t(theta$rho) ) *
                        (1 - state_probs$initial)
    
    Qs_first_derivative = sum(Qs_first_term) +
                        sum(sapply( 1:(n-1), function(k) sum(Qs_last_term[[k]]) ))
    
    return (Qs_first_derivative)
}