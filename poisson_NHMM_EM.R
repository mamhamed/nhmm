poisson_NHMM_EM = function(R, X, m, D, lambda, theta, maxiter = 30, tol = 1e-2, ...)
{
    # R is observed sequence
    # X[t,] is d-dimensional time covariate sequence
    # m is number of states
    # D is dimension of time covariate Xt
    # lambda_i is state-dependent Poisson distribution rate for modeling P(Rt | St = i)
    # gamma_ji is logistic function parameter for modeling P(St = i | St-1 = j, Xt = x)
    # delta_i is logistic function parameter for modeling P(S1 = i | X1 = x)
    # rho[i,] is d-dimensional state-dependent time covariate coefficients in 
    # .. logistic function for modeling P(St = i | St-1 = j, Xt = x) and P(S1 = i | X1 = x)
    
    cat(sprintf("******* fitting HMM with %s states ******", m), "\n")
    n = length(R)
    np = m*m + m + m*D + m
    
    lambda.next = lambda
    theta.next = theta

    for (iter in 1:maxiter)
    {
        cat("******* iter", iter)
        lallprobs = outer(R, lambda[1:m], dpois, log = TRUE)
        state_probs = compute_state_probs(m, n, X, theta)
        
        # forward-backward algorithm
        fb = pois_NHMM_lalphabeta(R, m, D, state_probs, lambda)
        la = fb$la
        lb = fb$lb
        c = max(la[,n])
        llk = c + log(sum(exp(la[,n] - c)))
        
        B = list(1:(n-1))
        for (k in 1:(n-1))
        {
            B[[k]] = state_probs$transition[[k]] * 
                        exp(matrix(rep(la[,k], m), m, m) - llk +
                            matrix(rep(lb[,k+1], m), m, m, byrow = TRUE) +
                            matrix(rep(lallprobs[k+1,], m), m, m, byrow = TRUE) )
            B[[k]] = B[[k]] / sum(B[[k]])
        }

        A = exp(la + lb - llk)
        lambda.next = A %*% R / apply(A, 1, sum)
        A = A / apply(A, 2, sum)
        
        theta.next = NHMM_conjugate_gradient(X, m, D, n, A, B, state_probs, lambda, theta)
#         theta_vector_optim = optim(par = c(theta$delta, theta$gamma, theta$rho), 
#                                    fn = compute_Qs_theta_vector, 
#                                    gr = gradient_Qs_theta_vector,
#                                    m = m, D = D, n = n, A = A, B = B, method = "BFGS")
#         theta.next = list(delta = theta_vector_optim$par[1:m],
#                           gamma = matrix(theta_vector_optim$par[m + 1:(m*m)], m, m),
#                           rho = matrix(theta_vector_optim$par[m + m*m + 1:(D*m)], m, D) )
            
        crit = sum(abs(lambda - lambda.next)) + 
                sum(abs(theta$gamma - theta.next$gamma)) + 
                sum(abs(theta$delta - theta.next$delta)) + 
                sum(abs(theta$rho - theta.next$rho))
        
        # for debugging
        if (is.nan(crit)) {
            print('poisson_NHMM_EM crit is nan')
            cat('lambda.next', lambda.next, '\n')
            cat('la', la, '\n')
            cat('lb', lb, '\n')
            cat('delta.next', theta.next$delta, '\n')
            cat('gamma.next', theta.next$gamma, '\n')
            cat('rho.next', theta.next$rho, '\n')
            return (NA)
        }
        
        AIC = -2*(llk - np)
        BIC = -2*llk + np * log(n)
        print(BIC)
        if (crit < tol)
        {
            return (list(lambda = lambda, delta = theta$delta, 
                         gamma = theta$gamma, rho = theta$rho,
                         mllk = -llk, AIC = AIC, BIC = BIC))
        }
        
        lambda = lambda.next
        theta = theta.next
    }
    
    print(paste("No convergence after ", maxiter, " iterations"))
    return (list(lambda = lambda, delta = theta$delta, 
                 gamma = theta$gamma, rho = theta$rho,
                 mllk = -llk, AIC = AIC, BIC = BIC))
}