compute_Qs_nu = function(nu, m, n, X, A, B, theta, phi)
{
    theta_nu = list(delta = theta$delta + nu * phi$delta,
                    gamma = theta$gamma + nu * phi$gamma,
                    rho = theta$rho + nu * phi$rho)
    
    state_probs = compute_state_probs(m, n, X, theta_nu)
    
    return (-compute_Qs(n, A, B, state_probs))
}