newton_raphson = function(m, n, X, A, B, state_probs, theta, phi,
                          maxiter_newton = 100, tol = 1e-6, epsilon = 1e-3, ...)
{
    nu = 1
    for (i in 1:maxiter_newton)
    {
        theta_nu = list(delta = theta$delta + nu * phi$delta,
                        gamma = theta$gamma + nu * phi$gamma,
                        rho = theta$rho + nu * phi$rho)
        
        state_probs = compute_state_probs(m, n, X, theta_nu)
        
        Qs_second_derivative = compute_d2Qs_dnu2(X, n, A, B, state_probs, theta) 
        
        #if (abs(Qs_second_derivative) < epsilon)
        #{
        #    print('WARNING: denominator is too small')
        #    break
        #}
        
        Qs_first_derivative = compute_dQs_dnu(X, n, A, B, state_probs, theta) 
        
        eta = 1
        Qs_eta_first_derivative = Qs_first_derivative
        while (TRUE)
        {
            nu.next = nu - eta * Qs_first_derivative / 1
            theta_eta = list(delta = theta$delta + nu.next * phi$delta,
                            gamma = theta$gamma + nu.next * phi$gamma,
                            rho = theta$rho + nu.next * phi$rho)
            
            Qs_old_first_derivative = Qs_eta_first_derivative
            Qs_eta_first_derivative = compute_dQs_dnu(X, n, A, B, state_probs, theta_eta)
            
            if ( eta != 1 &
                     Qs_eta_first_derivative * Qs_old_first_derivative < 0 ) break
            
            eta = eta / 2
        }
        
        if (abs(nu.next - nu) < tol) return (nu.next)
        
        nu = nu.next
    }
 
    cat("Warning: Unable to find solution to within desired tolerance of", tol)
    return (NA)
}