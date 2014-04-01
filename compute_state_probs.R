compute_state_probs = function(m, n, X, theta)
{
    transition_probs = list(1:(n-1))
    for (k in 1:(n-1))
    {
        s = theta$gamma + matrix(rep(X[k+1,] %*% t(theta$rho), m), m, m, byrow = TRUE)
        s_max = max(s)
        transition_probs[[k]] = exp(s - s_max)
        #transition_probs = transition_probs / matrixrep(rowSums(transition_probs), m), m, m)
        transition_probs[[k]] = transition_probs[[k]] / apply(transition_probs[[k]], 1, sum)
    }
    
    s = theta$delta + X[1,] %*% t(theta$rho)
    s_max = max(s)
    initial_probs = exp(s - s_max)
    initial_probs = initial_probs / sum(initial_probs)
    
    return (list(initial = initial_probs, transition = transition_probs))
}