# Spectral Envelope

The **spectral envelope** is a tool to study cyclic behaviors in categorical data. It is more informative than the traditional approach of attributing a different number to each category for power-spectral density estimation. <br/>

For each frequency in the spectrum, the **spectral envelope** finds an optimal real-numbered mapping that maximizes the normed power-spectral density at this point. Therefore, no matter what mapping is choosen for the different categories, the power-spectral density will always be bounded by the spectral envelope.

The spectral envelope was defined by David S. Stoffer in *DAVID S. STOFFER, DAVID E. TYLER, ANDREW J. MCDOUGALL, Spectral analysis for categorical time series: Scaling and the spectral envelope*.\

## Main functions
- - -
**<div class="border border-black-fade bg-red-light p-2 mb-2">
  spectral_envelope
</div>**
```Julia
spectral_envelope(ts; m = 3)
```
Computes the spectral envelope of an input categorical time-series.  
The degree of smoothing can be chosen by the user.

> **Parameters**:

>>* **ts** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array containing input categorical time-series.
>>* **m** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): Smoothing parameter. corresponds to how many neighboring points
        are to be involved in the smoothing (weighted average). Defaults to 3.  

> **Returns**: `(freq, se, eigev)` with `freq` the frequencies of the power-spectrum, `se` <br/> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the values of the spectral envelope for each frequency in 'freq'.
    `eigvecs` contains <br/> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the optimal real-valued mapping for each frequency point.
- - -
**get_mappings**

`get_mapping(data, freq; m = 3)`

Computes, for a given frequency 'freq', the optimal mappings for the categories in 'data'. Scans the vincinity of 'freq' to find the maximum of the spectral envelope, prints a sum up and returns the obtained mappings.
> **Parameters**:

>>* **data** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array containing input categorical time-series.
>>* **freq** ([Float](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): Frequency for which the mappings are wanted. The vincinity of 'freq' will be scaned to find maximal value of the spectral envelope.  
>>* **m** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): Smoothing parameter. corresponds to how many neighboring points
        are to be involved in the smoothing (weighted average). Defaults to 3.  

> **Returns**: `mappings` the optimal mappings for the found maxima around 'freq'.



## Example
