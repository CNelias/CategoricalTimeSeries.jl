using StatsBase
using Random


"""
Used in the computation of cramer's V and cohen's K.
Returns the lagged bivariate probability of two given categories, Pij.
Given i and j two categories, and l a lag,
Pij is the probability to have the category j at time t + l, if we have i at time t.

    inputs :
    - time series : the data to analyse
    - lag : the value of lag at which one wants Pij
    - category1 : the first category
    - category2 : the second category

    output :
    - Pij : the probability to observe j at t + l if we observe i at time t.
"""
function LaggedBivariateProbability(TimeSerie, Lag::Int, Category1, Category2)
    Pij = 0
    lagged_serie_length = length(TimeSerie) - Lag
    for i in 1:lagged_serie_length
        if (TimeSerie[i] == Category1) && (TimeSerie[i + Lag] == Category2)
            Pij += 1/(lagged_serie_length)
        end
    end
    return Pij
end

"""
If an array of lags 'Lag' is provided instead of a single value, computes Pij for each lag value and returns an array.
"""
function LaggedBivariateProbability(TimeSerie, Lags::Array{Int64,1}, Category1, Category2)
    Pij = [LaggedBivariateProbability(TimeSerie,L,Category1,Category2) for L in Lags]
    return Lags,Pij
end

"""
If no categories are explicitely provided, will a CxC symmetric probability matrix, with C the number of != categories.
the element [i,j] of this matrix corresponds to Pij at lag 'Lag'. Use this function if you want to have the contingency table at lag 'Lag'.
"""
function LaggedBivariateProbability(TimeSeries, Lag::Int64)
    categories = sort(unique(TimeSeries))
    Pij = zeros(length(categories), length(categories))
    lagged_serie_length = length(TimeSeries) - Lag
    for (idxi, i) in enumerate(categories)
        for (idxj, j) in enumerate(categories)
            for t in 1:lagged_serie_length
                if (TimeSeries[t] == i) && (TimeSeries[t + Lag] == j)
                        Pij[idxi, idxj] += 1/(lagged_serie_length)
                end
            end
        end
    end
    return Pij
end


"""
Returns the relative frequency of a given category (of type ::Any).
For example, if we have 3 category a,b and c, and a occurs 300 times in a time series of length 900,
its relative frequency  is 1/3.

"""
function relative_frequency(TimeSerie, idx::Any)
    Pi = 0
    for i in TimeSerie
        if i == unique(TimeSerie)[idx]
            Pi += 1/length(TimeSerie)
        end
    end
    return Pi
end

"""
Computes the lagged cohen coefficient κ describing how much the different categories
are correlated to each others at different lag values in the time-series. K is a measure of agreement.
input:
    -Serie : a categorical time series
    -lag: an integer lag value
output:
    -κ : cohen's coefficient
"""
function cohen_coefficient(Serie, lag = 1)
    if typeof(lag) != Int64
        throw("typeError : 'lag' value needs to be integer.")
    end
    categories = unique(Serie)
    K = 0
    pi_denominateur = 0
    rf_squared = [relative_frequency(Serie,i)^2 for i in 1:length(categories)]
    lagged_pij = LaggedBivariateProbability(Serie, lag)
    for (idxi, i) in enumerate(categories)
        K += (lagged_pij[idxi, idxi] - rf_squared[idxi])
        pi_denominateur =+ rf_squared[idxi]
    end
    K = K/(1-pi_denominateur)
    return K
end

"""
If 'lag' is provided as an array of lags, computes and returns an array of K for each elements of 'lag'.
"""
function cohen_coefficient(Serie, Lags::Array{Int64,1})
    if typeof(Lags) != Array{Int64,1}
        throw("typeError : 'lag' value needs to be integer.")
    end
    categories = unique(Serie)
    K = zeros(length(Lags))
    rf_squared = [relative_frequency(Serie,i)^2 for i in 1:length(categories)]
    for (lidx,l) in enumerate(Lags)
        k = 0
        pi_denominateur = 0
        lagged_pij = LaggedBivariateProbability(Serie, l)
        for (idxi, i) in enumerate(categories)
            k += (lagged_pij[idxi, idxi] - rf_squared[idxi])
            pi_denominateur =+ rf_squared[idxi]
        end
        K[lidx] =  k/(1-pi_denominateur)
    end
    return K
end


"""
    cramer_coefficient(Serie, Lags::Array{Int64,1})

Returns the cramer's V coefficient for the given 'Lags' values.

Cramer's v is used in categorical data analysis to test the degree of association between
a given set of categorical variable and another set of variables (example : eye color and hair color).
Here, it measures the association between the categorical values of a time series, and the values at a later lagged time.

V lies in [0,1]. 0 no association, 1 perfect association.
It is symmetric, meaning v(A,B) = v(B,A), so the information contained in the asymmetry between the variables is lost.

    input :
    - Serie : the time series containing the data
    - Lags : Array containing the lags at which V will be computed. If an int is given
                a single point value will be returned.

    output :
    - V : the values of v for the given lags.

"""
function cramer_coefficient(Serie, Lags::Array{Int64,1})
    if typeof(Lags) != Array{Int64,1}
        throw("typeError : 'Lags' needs to be an Array{Int64,1}.")
    end
    Categories = unique(Serie)
    V = zeros(length(Lags))
    d = length(Categories)-1
    rf = [relative_frequency(Serie, i) for i in 1:length(Categories)]
    for (lidx,lag) in enumerate(Lags)
        v = 0
        lagged_pij = LaggedBivariateProbability(Serie, lag)
        for i in 1:d+1
            for j in 1:d+1
                v =+ ( (lagged_pij[i,j]-rf[i]*rf[j])^2 ) / (rf[i]*rf[j])
            end
        end
        V[lidx] = sqrt(v/d)
    end
    return V
end

function cramer_coefficient(Serie, lag = 1)
    if typeof(lag) != Int64
        throw("typeError : 'lag' value needs to be integer.")
    end
    v = 0
    d = length(unique(Serie))-1
    lagged_pij = LaggedBivariateProbability(Serie, lag)
    rf = [relative_frequency(Serie,i) for i in 1:d+1]
    for i in 1:d+1
        for j in 1:d+1
            v =+ ((lagged_pij[i,j]-rf[i]*rf[j])^2)/(rf[i]*rf[j])
        end
    end
    return sqrt(v/d)
end

"""
    H(x)

Estimates the probability distribution of the categories in 'x' and computes its entropy.
It characterizes the amount of information carried by the distribution.
it tells how "disperse" the distribution is.

"""
function H(x)
    x_count = countmap(x)
    total_occurences = sum(values(x_count))
    entropy = 0.0
    for occurence in values(x_count)
        px = occurence/total_occurences
        entropy -= px*log2(px)
    end
    return entropy
end

"""
Computing the conditional entropy of y given x: H(Y|X),
Given in bits or shannons (meaning using the log2 for the measurement).
Wikipedia: https://en.wikipedia.org/wiki/Conditional_entropy

"""
function conditional_entropy(y,x)
    xy_count = countmap(collect(zip(y,x)))
    x_count = countmap(x)
    total_occurences = sum(values(xy_count))
    entropy = 0.0
    for xy in keys(xy_count)
        p_xy = xy_count[xy]/total_occurences
        p_x = x_count[xy[2]]/total_occurences
        entropy += p_xy*log2(p_x/p_xy)
    end
    return entropy
end


"""
    theils_u(x, Lags)

Measures what portion of the information associated to the values in 'x' is know at t + lag if the value 'x' is known at t.
Input :
    x : a categorical time-series
    lags : an array of lags onto U is computed.
returns:
    U : an array of U for the given values in 'Lags'.
"""
function theils_u(x, Lags)
    return [(H(x[1:end-l]) - conditional_entropy(x[1:end-l], x[l+1:end])) / H(x[1:end-l]) for l in Lags]
end

"""
    rate_evolution(Series)

A way to test the stationarity of the input categorical time-series 'Series'.
if it varies linearly, then the time series is more or less stationary.
"""
function rate_evolution(Series)
    categories = unique(Series)
    RATE = []
    for c in categories
        init = zeros(length(Series))
        for v in 1:length(Series)
            if Series[v] == c
                init[v] = 1
            end
        end
        push!(RATE,cumsum(init))
    end
    return RATE
end


"""
    bootstrap_CI(Series, lags, coefficient, n_iter = 1000)

Provides 95% a confidence interval by shuffling 'Series' 'n_iter' times, computing the values of 'coefficient'
then finding the value of the top and bottom 2,5% to get the interval.
This is done for every point in 'lags' (can be costly if 'Series' is long).

Input :
    - Series : input array of categorical data
    - lags (Array{Int64,1}) : the lag values at which the analysis is conducted
    - coef_func (**function**) : the function for which the CI needs to be computed.
            'coefficient' can be one of the following **functions** : 'cramer_coefficient, cohen_coefficient, theils_U'.
    - n_iter (Int64) : how many iterations are run for the bootstrap procedure.
returns :
    - top (Array{Float64,1}) : Array of values for the upper limit of the CI.
    - bottom (Array{Float64,1}) : Array of values for the lower limit of the CI.
"""
function bootstrap_CI(Series,  lags, coef_func, n_iter = 1000)
    critical_values = zeros(n_iter)
    bootstrap_storage = zeros(n_iter, length(lags))
    for i in 1:n_iter
        bootstrap_storage[i,:] = coef_func(shuffle(Series), lags)
    end
    top_values = sort(bootstrap_storage, dims = 1)[n_iter - div(n_iter,40),:]
    bottom_values = sort(bootstrap_storage, dims = 1)[div(n_iter,40),:]
    return top_values, bottom_values
end

export cramer_coefficient, cohen_coefficient, conditional_entropy, LaggedBivariateProbability, H, theils_u, rate_evolution
