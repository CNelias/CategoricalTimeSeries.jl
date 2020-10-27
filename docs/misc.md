# Miscalleneous
Here we regroup all additional useful functions that do not necessarily deserve a section of their own.

- - -
**rate_evolution — Function**
- - -
```
rate_evolution(series)
```
The rate of evolution is a way to test the stationarity of a categorical time-series.
If the rate evolves more or less linearly, then the time-series can reasonably be considered stationary.
It is most informative to plot the rate of evolution of each categories on the same graph for a direct visual inspection.
> **Parameters**:

>>* **series** ([Array{any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array of categorical time-series.

> **Returns**: `RATE`, Array containing, for each category, an array representing it's rate of evolution.

- - -
**LaggedBivariateProbability — Function**
- - -
```
LaggedBivariateProbability(serie, Lags::Array{Int64,1}, Category1, Category2)
```
Returns the lagged bivariate probability of two given categories, Pij.
Given i and j two categories, and l a lag (or array of lags),
Pij is the probability to have the category j at time t + l, if we have i at time t.
> **Parameters**:

>>* **serie** ([Array{any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array of categorical time-series.
>>* **category1**
>>* **category2**

> **Returns**: `pij`, Array containing, for each value in `lags`, the lagged bivariate probability.

- - -
**varcov — Function**
- - -
```Julia
varcov(ts::Array{Float64,2})
```
Computes the covariance-variance matrix of a given multivariate time-series. This can also be used for a univariate time-series but the input should still be 2-D.
> **Parameters**:

>>* **ts** ([Array{Float,2}](https://docs.julialang.org/en/v1/base/arrays/)): 2-D input array of multivariate time-series.

> **Returns**: `cov_matrix` the correpsonding covariance matrix.

- - -
**power_spectrum — Function**
- - -
```Julia
power_spectrum(x::Array{Float64,1}, window::Int, step::Int)
```
Computes an estimation of the power-spectrum of the input time-series `x`.
> **Parameters**:

>>* **x** ([Array{Float,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array of real-valued time-series.
>>* **window** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): Integer specifying the size of the window for averaging Must be shorter than length(x). Recommended value is 1/10th of length(x).
>>* **step** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): Parameters controlling the overlap between the windows. Shouldn't be biggger than div(window,2).  

> **Returns**: `pxx`, the estimated power-spectrum.
