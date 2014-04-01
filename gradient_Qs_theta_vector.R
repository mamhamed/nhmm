gradient_Qs_theta_vector = function(theta_vector, m, D, n, A, B)
{
    theta_list = list(delta = theta_vector[1:m],
                      gamma = matrix(theta_vector[m + 1:(m*m)], m, m),
                      rho = matrix(theta_vector[m + m*m + 1:(D*m)], m, D) )
    
    state_probs = compute_state_probs(m, n, X, theta_list)
    
    gradQ_list = gradient_Qs(m, D, n, A, B, state_probs)
    
    gradQ_vector = c(-gradQ_list$delta, -gradQ_list$gamma, -gradQ_list$rho)
    return (gradQ_vector)
}