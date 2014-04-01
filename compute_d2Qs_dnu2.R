compute_d2Qs_dnu2 = function(X, n, A, B, state_probs, theta)        
{
    # calculate second derivative d2_Qs / d_nu_2
    Qs_last_term = list(1:(n-1))
    for (k in 1:(n-1))
    {
        Qs_last_term[[k]] = B[[k]] * (theta$gamma + 
                    matrix(rep(X[k+1,] %*% t(theta$rho), m), m, m, byrow = TRUE))^2 *
                    (1 - state_probs$transition[[k]]) * state_probs$transition[[k]]
    }
    Qs_first_term = A[,1] * (theta$delta + X[1,] %*% t(theta$rho) )^2 *
                    (1 - state_probs$initial) * state_probs$initial
    
    Qs_second_derivative = -sum(Qs_first_term) - 
                    sum(sapply( 1:(n-1), function(k) sum(Qs_last_term[[k]]) ))
    
    return (Qs_second_derivative)
}