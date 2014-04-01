compute_Qs = function(n, A, B, state_probs)
{
    Qs_last_term = list(1:(n-1))
    for (k in 1:(n-1))
    {
#         if ( any(is.na(state_probs$transition[[k]])) | any(state_probs$transition[[k]] <= 0 ) )
#         {
#             print(state_probs$transition[[k]])
#         }
        Qs_last_term[[k]] = B[[k]] * log(state_probs$transition[[k]])
    }
    
#     if ( any(is.na(state_probs$transition[[k]])) | any(state_probs$initial <= 0 ) )
#     {
#         print(state_probs$initial)
#     }
    Qs_first_term = A[,1] * log(state_probs$initial)
    Qs = sum(Qs_first_term) + sum(sapply( 1:(n-1), function(k) sum(Qs_last_term[[k]]) ))
    
    return (Qs)
}