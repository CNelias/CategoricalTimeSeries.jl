# Correlations
The study of categorical data prevents the usage of standard tools like the autocorrelation function, as they are often not defined. The following functions provide ways to study categorical serial dependences.  
Most of these methods are described in C. Weiss's book "*An Introduction to Discrete-Valued Time Series*" (2018)[1].
## Main functions
- - -
**cramer_coefficient — Function**
- - -
```Julia
cramer_coefficient(series, lags)
```
Measures average association between elements of ```series``` at time t and time t + ```lags```. Cramer's V is an unsigned measurement : its values lies in [0,1], 0 being perfect independence and 1 perfect dependence. k can be biased, for more informations, refer to [1].
> **Parameters**:

>>* **series** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array containing input categorical time-series.
>>* **lags** ([Array{Int,1}](https://docs.julialang.org/en/v1/base/arrays/)): lag values at which cramer's coefficient is computed. Alternatively, `lags` can be an integer, a single integer value will then be returned.  

> **Returns**: `V`, the value of cramer's coefficient for each value in `lags`.

- - -
**cohen_coefficient — Function**
- - -
```Julia
cohen_coefficient(series, lags)
```
Measures average association between elements of ```series``` at time t and time t + ```lags```.
Cohen's k is a signed measurement : its values lie in [-pe/(1 -pe), 1], with positive (negative) values indicating positive (negative) serial dependence at `lags`. pe is probability of agreement by chance.

> **Parameters**:

>>* **series** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array containing input categorical time-series.
>>* **lags** ([Array{Int,1}](https://docs.julialang.org/en/v1/base/arrays/)): lag values at which Cohen's coefficient is computed. Alternatively, `lags` can be an integer, a single integer value will then be returned.  

> **Returns**: `K`, the value of Cohen's coefficient for each value in `lags`.

- - -
**theils_u — Function**
- - -
```Julia
theils_u(series, Lags)
```
Measures average portion of information known about `series` at t + `lags` given that `series` is known at time t. Theil's U makes use of concepts borrowed from *information theory*
U is an unsigned measurement: its values lies in [0,1], 0 meaning no information shared and 1 complete knowledge (determinism).

> **Parameters**:

>>* **series** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array containing input categorical time-series.
>>* **lags** ([Array{Int,1}](https://docs.julialang.org/en/v1/base/arrays/)): lag values at which Theil's U is computed. Alternatively, `lags` can be an integer, a single integer value will then be returned.  

> **Returns**: `U`, the value of Theil's U for each value in `lags`.

## Confidence interval
Depending on the length of the time-series and the method used, the estimated value of serial dependence might fluctuate a lot around its true value.
It is therefore useful to relate estimations to a corresponding confidence interval to know how significant given results are. The following function provides a confidence interval via bootstrap:
- - -
**bootstrap_CI — Function**
- - -
```Julia
bootstrap_CI(series, lags, coef_func, n_iter = 1000, interval_size = 0.95)
```
Returns a top and bottom limit of the confidence interval at values of `lags`. The width of the confidence interval can be choosen (defaults to 95%). The returned confidence interval corresponds to the null hypothesis (no serial dependence), if the estimated serial dependence lies in this interval, no significant correlations can be claimed.

> **Parameters**:

>>* **series** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array containing input categorical time-series.
>>* **lags** ([Array{Int,1}](https://docs.julialang.org/en/v1/base/arrays/)): lag values at which the CI is computed.
>>* **coef_func** ([function](https://docs.julialang.org/en/v1/manual/functions/)): the function for which the CI needs to be computed.
            `coef_func` can be one of the following **functions** : `cramer_coefficient`, `cohen_coefficient` or `theils_U`.
>>* **n_iter** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): number of iterations for the bootstrap procedure. The higher, the more precise but more computationaly demanding. Defaults to 1000.
>>* **interval_size** ([Float](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): Desired size of the confidence interval. Defaults to 0.95, for a 95% confidence interval.

> **Returns**: `(top_values, bottom_values)`, the top and bottom limit for confidence interval, for each point in `lags`.

## Example
Using the Pewee [birdsong data](https://github.com/johncwok/CategoricalTimeSeries.jl/tree/main/test) (1943) one can do a serial dependence plot using Cohen's cofficient as follow :
```
using DelimitedFiles, Plots
using CategoricalTimeSeries

#reading 'pewee' time-series test folder.
data_path = joinpath(dirname(dirname(pathof(CategoricalTimeSeries))), "test", "pewee.txt")
series = readdlm(data_path,',')[1,:]
lags = collect(1:25)
v = cohen_coefficient(series, lags)
t, b = bootstrap_CI(series, lags, cramer_coefficient)
a = plot(lags, v, xlabel = "Lags", ylabel = "K", label = "Cramer's k")
plot!(a, lags, t, color = "red", label = "Limits of 95% CI"); plot!(a, lags, b, color = "red", label = "")
```

<img src=https://user-images.githubusercontent.com/34754896/136663737-f30f20bf-c42b-4979-b514-637be8b7f404.PNG width = "600">


[1] DOI : 10.1002/9781119097013
