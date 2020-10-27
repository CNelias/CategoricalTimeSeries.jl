# Spectral Envelope

The **spectral envelope** is a tool to study cyclic behaviors in categorical data. It is more informative than the traditional approach of attributing a different number to each category for power-spectral density estimation. <br/>

For each frequency in the spectrum, the **spectral envelope** finds an optimal real-numbered mapping that maximizes the normed power-spectral density at this point. Therefore, no matter what mapping is choosen for the different categories, the power-spectral density will always be bounded by the spectral envelope.

The spectral envelope was defined by David S. Stoffer in *DAVID S. STOFFER, DAVID E. TYLER, ANDREW J. MCDOUGALL, Spectral analysis for categorical time series: Scaling and the spectral envelope*.

## Main functions
- - -
**spectral_envelope — Function**
- - -
```Julia
spectral_envelope(ts; m = 3)
```
Computes the spectral envelope of an input categorical time-series.  
The degree of smoothing can be chosen by the user.

> **Parameters**:

>>* **ts** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array containing input categorical time-series.
>>* **m** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): Smoothing parameter. corresponds to how many neighboring points
        are to be involved in the smoothing (weighted average). Defaults to 3.  

> **Returns**: `(freq, se, eigev)`, with `freq` the frequencies of the power-spectrum, `se` <br/> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the values of the spectral envelope for each frequency in 'freq'.
    `eigvecs` contains <br/> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the optimal real-valued mapping for each frequency point.

- - -
**get_mappings — Function**
- - -
```
get_mapping(data, freq; m = 3)
```

Computes, for a given frequency `freq`, the optimal mappings for the categories in `data`. Scans the vincinity of `freq` to find the maximum of the spectral envelope, prints a sum up and returns the obtained mappings.
> **Parameters**:

>>* **data** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array containing input categorical time-series.
>>* **freq** ([Float](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): Frequency for which the mappings are wanted. The vincinity of 'freq' will be scaned to find maximal value of the spectral envelope.  
>>* **m** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): Smoothing parameter. corresponds to how many neighboring points
        are to be involved in the smoothing (weighted average). Defaults to 3.  

> **Returns**: `mappings`, the optimal mappings for the found maxima around 'freq'.



## Example
Applying the spectral envelope to study a [segment of DNA](https://github.com/johncwok/SpectralEnvelope.jl/tree/master/test) from the Epstein-Barr virus and plotting the results:
```
using DelimitedFiles, Plots

data = readdlm("..\\test\\DNA_data.txt")
f, se, eigvecs = spectral_envelope(data; m = 4)

plot(f, se, xlabel = "Frequency", ylabel = "Intensity", title = "test data: extract of Epstein virus DNA", label = "spectral envelope")
```
<img src=https://user-images.githubusercontent.com/34754896/91556982-eef72680-e933-11ea-85f3-fab6aea17258.PNG width = "500">

To get the associated optimal mapping for the peak at frequency 0.33:
```
mappings = get_mappings(data, 0.33)
>> position of peak: 0.33 strengh of peak: 0.6
print(mappings)
>> ["A : 0.54", "G : 0.62", "T : -0.57", "C : 0.0"]
```


## Additional functions
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

- - -
**varcov**
- - -
```Julia
varcov(ts::Array{Float64,2})
```
Computes the covariance-variance matrix of a given multivariate time-series. This can also be used for a univariate time-series but the input should still be 2-D.
> **Parameters**:

>>* **ts** ([Array{Float,2}](https://docs.julialang.org/en/v1/base/arrays/)): 2-D input array of multivariate time-series.

> **Returns**: `cov_matrix` the correpsonding covariance matrix.
