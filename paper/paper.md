---
title: 'CategoricalTimeSeries.jl: A toolbox for categorical time-series analysis'
tags:
  - Julia
  - Categorical Time-series
  - Spectral analysis
  - Association measurement
  - Clustering

authors:
 - name: Corentin Nelias
   orcid: 0000-0001-6266-5575 
   affiliation: "1, 2"
affiliations:
 - name: Max Planck Institute for Dynamics and Self-Organization
   index: 1
 - name: Department of Physics, Georg-August-Universität Göttingen
   index: 2
date: 8 septembre 2021
bibliography: paper.bib
---

# Introduction

**CategoricalTimeSeries.jl** is a Julia @Julia toolbox made to help with the analysis of categorical time-series. 

The common approach to deal with categorical time-series consists in transforming the data via a mapping to obtain a real-valued time-series.
This enables the use of traditional time-series analysis methods. However, most of these methods (power-spectral density estimation, correlation coefficients, dimensionality reduction etc...)
are not invariant under general transformations and will produce different results based on the choice of mapping. 
Therefore, depending on the type of categorical data and the problem at hand, it can be desirable to have methods that work with the direct categorical values themselves. 

The purpose of **CategoricalTimeSeries.jl** is to provide such tools. The package comes with an extensive documentation available online: https://categoricaltimeseriesjl.readthedocs.io/en/latest/

# Statement of need

While several implementations of categorical time-series analysis methods are accessible online, they are written in different languages, some of which are not free (e.g. Matlab). Additionally, no implementation for methods such as the *spectral envelope* of the *random projection* (see *Overview of functionalities* below) were available online. This package centralizes and implements most of the standard methods of categorical time-series analysis in a single toolbox fully written in the Julia language. 


# Overview of functionalities

**CategoricalTimeSeries.jl** was designed to be easy to use and produce results that are simple to plot. 
Consequently, the methods implemented in the package take as input 1-D arrays of any type. 
Type conversion and pre-processing (when needed) are done inside of the methods automatically without the need for additional coding from the user.
The results are either formatted in a way that can be plotted directly with the build-in ```Plots``` pacakge, or a helper function is provided for visualization and interpretation.

We now present the main functionalities provided by **CategoricalTimeSeries.jl**:
### Spectral analysis
The spectral envelope method [@Stoffer:1998] is used to study the power-spectrum of categorical time-series. 
As stated in the introduction, the power-spectrum of a time-series is not invariant under a generic transformation. 
A wrong choice of mapping can potentially flatten certain peaks and render them unnoticeable.
For each frequency, the spectral envelope seeks the mapping that maximizes the value of the power-spectrum normalized by the total variance.
The ```spectral_envelope``` function takes as input a time-series (1-D array) and returns all the frequencies of the spectrum and the values of the intensity associated with the optimal mappings. It also returns the mappings.
For a finer study of the mappings themselves, the ```get_mappings``` function can be used.

### Association analysis
The notion of autocorrelation function is formally not defined for a categorical time-series [@Weiss:2018].
Yet, it might be of interest to know how inter-dependent the values of the time-series are. 
We implemented several coefficients generalizing the concept of linear correlations to categorical time-series.
Cramer's coefficient, Cohen's coefficient and Theil's U can respectively be computed via the ```cramer_coefficient```, ```cohen_coefficient``` and ```theils_u``` functions.
They take as input a 1-D array representing the time-series to study and an array of lags storing the lag values at which the coefficients are evaluated. For example:
```Julia
cramer_coefficient(series, lags)
```

### Motif recognition
Time-series can present repeating motifs that are worthwhile identifying. However, simple line-search algorithm are not adapted for all motifs (@Pevzner2000)
Moreover, the lack of proper distance measurement complicates the saerch in the context of categorical time-series.  
An implementation using the *random projection* method [@Buhler:2002] is used here.
The identification of potential motifs is done with the ```detect_motifs``` function.
It takes as required input a time-series (1-D array), the length of the motifs to look for, and the number of allowed errors. 
It returns an instance of the ```pattern``` structure which stores properties of the identified motif such as shapes, repetition number and positions.


### Data clustering.
If certain categories inside a time-series present functional similarities, one might wish to clusther them together into a single equivalent representation.
This reduces the total number of categories and can simplify the analysis of the time-series. To do this, we use an implementation based on the *Information bottleneck* concept [@Tishby:2000; @Strouse:2016].
After an initial bottleneck model of the structure ```IB``` is instantiated, it can be optimized with the ```IB_optimize!``` function to reveal potential clusters of categories. An overview of the results can be obtained with the ```print_results``` function. 

# Acknowledgements

The author thanks Nori Jacoby discussing and providing hinsight on the *Information bottleneck* concept.

# References
