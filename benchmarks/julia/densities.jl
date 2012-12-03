# Examples borrowed from: http://arxiv.org/pdf/1211.3759.pdf

# 2/3 N(x|0,1) + 1/3 N(x|0,1/10)
function kurtotic(theta)
    if theta[3] <= 0 || theta[6] <= 0
        return -Inf
    end
    llk = 0
    for i in 1:length(y)
      llk += log(theta[1] * pdf(Normal(theta[2], theta[3]), y[i]) + 
                 theta[4] * pdf(Normal(theta[5], theta[6]), y[i]))
    end
    return llk
end

function bimodal_density(x)
    lik = (1/2) * pdf(Normal(-1,2/3),x) + (1/2) * pdf(Normal(1,2/3),x)
    return log(lik)
end

# 1/2 N(x|-1,2/3) + 1/2 N(x|1,2/3)
function bimodal_loglikelihood(theta)
    if theta[3] < 0 || theta[6] < 0
        return -Inf
    end
    llk = 0
    for i in 1:length(y)
      llk += log(theta[1] * pdf(Normal(theta[2], theta[3]), y[i]) + 
                 theta[4] * pdf(Normal(theta[5], theta[6]), y[i]))
    end
    return -llk
end

# requires y and x
function logistic_regression(theta)
    llk = 0
    eta = [1 / (1 + exp(- x[i,:] * theta))[1] for i in 1:length(y)]
    for i = 1:length(y)
        llk += y[i] * log(eta[i]) + (1-y[i]) * log(1-eta[i])
    end
    return llk
end

