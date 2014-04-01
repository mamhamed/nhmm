NHMM_conjugate_gradient = function(X, m, D, n, A, B, state_probs, lambda, theta0,
                                        maxiter_conjugate = 30, tol = 1e-4, ...)
{
    theta = theta0
    gradQ = gradient_Qs(m, D, n, A, B, state_probs)
    phi = list(delta = -gradQ$delta, gamma = -gradQ$gamma, rho = -gradQ$rho)
    Qs = compute_Qs(n, A, B, state_probs)

    for (iter in 1:maxiter_conjugate)
    {
        #nu = newton_raphson(m, n, X, A, B, state_probs, theta, phi)
        #nu = uniroot(compute_dQs_dnu_nu, interval = c(-1, 1), 
        #             m = m, n = n, X = X, A = A, B = B, theta = theta, phi = phi)
        #nu = multiroot(compute_dQs_dnu_nu, start = 1, 
        #               m = m, n = n, X = X, A = A, B = B, theta = theta, phi = phi)
        #nl_min = nlm(compute_Qs_nu, p = 1,
        #             m = m, n = n, X = X, A = A, B = B, theta = theta, phi = phi)
        #nu = nl_min$estimate
        
        nu_optim = optim(1, compute_Qs_nu, gr=NULL, 
                         m = m, n = n, X = X, A = A, B = B, theta = theta, phi = phi, 
                         method = "BFGS")
        nu = nu_optim$par
        
        theta.next = list(delta = theta$delta + nu * phi$delta,
                          gamma = theta$gamma + nu * phi$gamma,
                          rho = theta$rho + nu * phi$rho)
        
        state_probs = compute_state_probs(m, n, X, theta.next)
        
        gradQ.next = gradient_Qs(m, D, n, A, B, state_probs)
        eta = sum( (gradQ.next$delta - gradQ$delta) * gradQ.next$delta,
                   (gradQ.next$gamma - gradQ$gamma) * gradQ.next$gamma,
                   (gradQ.next$rho - gradQ$rho) * gradQ.next$rho ) /
            sum( gradQ$delta * gradQ$delta, 
                 gradQ$gamma * gradQ$gamma,
                 gradQ$rho * gradQ$rho )
        phi.next = list(delta = gradQ$delta - eta * phi$delta, 
                        gamma = gradQ$gamma - eta * phi$gamma, 
                        rho = gradQ$rho - eta * phi$rho)
        
        Qs.next = compute_Qs(n, A, B, state_probs)
        crit = abs(Qs.next - Qs)
        
        # for debugging
        if (is.nan(crit)) {
            print('NHMM_conjugate_gradient crit is nan')
            cat('theta.next$delta', theta.next$delta, '\n')
            cat('theta.next$gamma', theta.next$gamma, '\n')
            cat('theta.next$rho', theta.next$rho, '\n')
            cat('gradQ.next$delta', gradQ.next$delta, '\n')
            cat('gradQ.next$gamma', gradQ.next$gamma, '\n')
            cat('gradQ.next$rho', gradQ.next$rho, '\n')
            cat('phi.next$delta', phi.next$delta, '\n')
            cat('phi.next$gamma', phi.next$gamma, '\n')
            cat('phi.next$rho', phi.next$rho, '\n')
            cat('Qs.next', Qs.next, '\n')
            return (NA)
        }
        if (crit < tol)
            return (theta)
        
        theta = theta.next
        gradQ = gradQ.next
        phi = phi.next
        Qs = Qs.next
    }
    print(paste("Conjugate gradient no convergence after ", maxiter_conjugate, " iterations"))
    return (theta)
}