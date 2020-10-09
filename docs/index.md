# Categorical Time-Series Analysis

## Introduction
**CategoricalTimeSeries.jl** is a Julia package regrouping methods of categorical time-series analysis. The term *categorical* commonly refers to two types of data: *nominal* and *ordinal*. <br/>
**Nominal**, or labeled values represent *discrete* units that have no intrinsic *order*, like common types of pet:

- Cat
- Dog
- Bird

**Ordinal** values on the other hand represent *discrete* and *ordered* units, like the size of a coffee cup:

1. Small
2. Medium
3. Large

When categorical data is layed out in function of time, one speaks of *categorical time-series*. <br/>
<br/>
Often, especially when dealing with ordinal time-series, it is enough to map the different values to a set of integers to carry the analysis. However, when it is not sufficient **CategoricalTimeSeries.jl** is here to help.
## Overview
The package lets you carry three main kind of analysis: **Spectral analysis**, **Data clustering** and **correlations analysis**. Other functionnalities are also avalaible (see misc.). Here is a quick overview of these methods, for more details go to the specific sections.
#### Spectral analysis
The standard approach to study spectral properties in categorical time-series is to map the different values to a set of numbers. While the overall shape of the spectrum is usually unaffected by this operation, peaks representing cyclic behaviors in the time-series can completely disappear depending on the choice of mapping. To tackle this issue, this package proposes an implementation of the spectral envelope method. It was developed by *David S. Stoffer*  in order to identify optimal mappings.
#### Data clustering
A categorical time-series might have many different possible values, but these values are rarely independant. Some categories are usually equivalent to one-another. For example, among the values `"male"`, `"female"`, `"house"` and `"car"`, it is evident that `"male"` and `"female"` refer to similar concepts and could be grouped in one single category. Data clustering is the task of identifying relationships of equivalency in the time-series in order to reduce it to a more efficient representation. **CategoricalTimeSeries.jl** offers an implementation of the *information bottleneck* to this end.
#### Correlation analysis
The notion of autocorrelation function is formally not defined for a categorical time-series. Yet, it might be of interest to know how interdependent the values of the time-series are. Efforts have been made to generalize the concept of linear correlations to categorical time-series. This package implements several of these methods.
#### Misc.
Additional functions are also provided in the package, allowing for instance to test for stationarity, look for change points or simply get confidence intervals.

##Installation
Like any other Julia package, **CategoricalTimeSeries.jl** can be installed with
```Julia
using Pkg
Pkg.clone("https://github.com/johncwok/CategoricalTimeSeries.jl.git")
```
To use it, you need to import it:
```Julia
using CategoricalTimeSeries
```
