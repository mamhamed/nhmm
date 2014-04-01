if(!exists("pois_NHMM_lalphabeta", mode = "function")) source("pois_NHMM_lalphabeta.R")
if(!exists("poisson_NHMM_EM", mode = "function")) source("poisson_NHMM_EM.R")
if(!exists("pois_NHMM_state_probs", mode = "function")) source("../pois_NHMM_state_probs.R")
if(!exists("NHMM_conjugate_gradient", mode = "function")) source("NHMM_conjugate_gradient.R")
if(!exists("compute_Qs", mode = "function")) source("compute_Qs.R")
if(!exists("gradient_Qs", mode = "function")) source("gradient_Qs.R")
if(!exists("compute_state_probs", mode = "function")) source("compute_state_probs.R")
if(!exists("newton_raphson", mode = "function")) source("newton_raphson.R")
if(!exists("compute_dQs_dnu", mode = "function")) source("compute_dQs_dnu.R")
if(!exists("compute_d2Qs_dnu2", mode = "function")) source("compute_d2Qs_dnu2.R")
if(!exists("pois_NHMM_forecast", mode = "function")) source("pois_NHMM_forecast.R")

library(rootSolve)

# constants definitions
min_m = 2
max_m = 4
all_lambda0 = c(0, 1, 4, 8, 16, 32, 64)
num_truncate = 141#24*100 #num_bins
num_predict = 16

#x = pois_HMM_generate_sample(2000, m, lambda0, gamma0, delta0)
y = c(rep(c(0, 6, 15, 4, 3, 18, 0), 20), 0)
y = y + rpois(num_truncate, 1)
# portal_index = 12
# y = access_counts[[portal_index]]
# y = y[1:num_truncate+18000]
# cat('portal url is', substr(dataIn_times$times[url_indices[portal_index]], 0, 50), '\n')

D = 7
X = t(matrix(rep(diag(D), (num_truncate + num_predict - 1) %/% D + 1), 
             D, num_truncate + num_predict))

# model selection
model_result = list()
theta = list()
for (m in min_m:max_m) 
{
    start_time = Sys.time()
    cat('number of state is', m, '\n')
    #lambda0 = all_lambda0[1:m]
    lambda0 = sapply(1:m-1, function(x) 2^x)
    gamma0 = matrix(rep(1, m*m), m, m)
    delta0 = rep(1, m)
    rho0 = matrix(rep(1, m*D), m, D)
    theta0 = list(delta = delta0, gamma = gamma0, rho = rho0)
    p_nhmm = poisson_NHMM_EM(y, X[1:num_truncate,], m, D, lambda0, theta0)
    
    if (any(is.na(p_nhmm))) 
        if (length(model_result$num_state) == 0) {
            next
        } else break
    
    model_result$num_state = c(model_result$num_state, m)
    model_result$aic_list = c(model_result$aic_list, p_nhmm$AIC)
    model_result$bic_list = c(model_result$bic_list, p_nhmm$BIC)
    model_result$lambda = c(model_result$lambda, p_nhmm$lambda)
    model_result$delta = c(model_result$delta, p_nhmm$delta)
    model_result$gamma = c(model_result$gamma, p_nhmm$gamma)
    model_result$rho = c(model_result$rho, p_nhmm$rho)
    model_result$mllk = c(model_result$mllk, p_nhmm$mllk)
    finish_time = Sys.time()
    print(finish_time - start_time)
}

optimal_index = which.min(model_result$bic_list)
m = model_result$num_state[optimal_index]
lambda = model_result$lambda[(m*(m-1)/2 - min_m*(min_m-1)/2) + 1:m]
theta = list(delta = model_result$delta[(m*(m-1)/2 - min_m*(min_m-1)/2) + 1:m],
             gamma = matrix(model_result$gamma[(m*(m-1)*(2*m-1)/6 - 
                            min_m*(min_m-1)*(2*min_m-1)/6) + 1:(m*m)], m, m),
             rho = matrix(model_result$rho[D*(m*(m-1)/2 - min_m*(min_m-1)/2) + 
                                               1:(D*m)], m, D) )
stateprobs = pois_NHMM_state_probs(y, X[1:num_truncate,], m, D, lambda, theta)
forecasts = pois_NHMM_forecast(y, X, m, D, lambda, theta, H = num_predict)

# print and plot results
cat("\n")
cat("number of HMM states fitted =", model_result$num_state, "\n\n")
cat("AIC =", model_result$aic_list, "\n")
cat("BIC =", model_result$bic_list, "\n\n")
cat("chosen HMM poisson lambda =", lambda, "\n\n")
print(theta)

par(mfrow = c(1, 1))
plot(y, type = 'h',
     xlab = paste('time in multiple of', period_in_seconds / 3600, 'hours'),
     ylab = 'access count',
     main = substr(dataIn_times$times[url_indices[portal_index]], 0, 50))

if (m < 4) {
    par(mfrow = c(m, 1)) 
} else if (m < 7) {
    par(mfrow = c( (m-1)%/%2 + 1, 2))
} else 
    par(mfrow = c(2, 2))

plot_colors = c('blue', 'red', 'green', 'magenta', 'purple', 'pink', 'brown')
for (i in 1:m)
{
    #layout(1)
    plot(1:num_truncate, stateprobs[i,], col = plot_colors[i], type = 'h', 
         xlab = paste('time in multiple of', period_in_seconds / 3600, 'hours'), 
         ylab = 'probability', 
         main = sprintf('state %d, lambda = %0.3f', i, lambda[i]) )
}

if (num_predict < 4) {
    par(mfrow = c(num_predict, 1)) 
} else if (num_predict < 7) {
    par(mfrow = c( (num_predict-1)%/%2 + 1, 2))
} else 
    par(mfrow = c(2, 2))

for (i in 1:num_predict)
    plot(forecasts$xrange, forecasts$fdists[,i], 
         col = plot_colors[i %% 7 + 1], type = 'h', 
         xlab = paste('access count in', period_in_seconds / 3600, 'hours'), 
         ylab = 'probability', 
         main = sprintf('%d-th period of %d hours from last observation', 
                        i, period_in_seconds / 3600) )