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

# Overview

**CategoricalTimeSeries.jl** was designed to be easy to use and produce results that are simple to plot. 
Consequently, the methods implemented in the package take as input 1-D arrays of any type. 
Type conversion and pre-processing (when needed) are done inside of the methods automatically without the need for additional coding from the user.
The results are either formatted in a way that can be plotted directly with the build-in ```Plots``` pacakge, or a helper function is provided for visualization and interpretation.

We now present the main functionalities provided by **CategoricalTimeSeries.jl**:
#### Spectral analysis
The spectral envelope method (@Stoffer1998) is used to study the power-spectrum of categorical time-series. 
As stated in the introduction, the power-spectrum of a time-series is not invariant under a generic transformation. 
A wrong choice of mapping can potentially flatten certain peaks and render them unnoticeable.
For each frequency, the spectral envelope seeks the mapping that maximizes the value of the power-spectrum normalized by the total variance.
The ```spectral_envelope``` function takes as input a time-series (1-D array) and returns all the frequencies of the spectrum and the values of the intensity associated with the optimal mappings. It also returns the mappings.
For a finer study of the mappings themselves, the ```get_mappings``` function can be used.

#### Association analysis
The notion of autocorrelation function is formally not defined for a categorical time-series. 
Yet, it might be of interest to know how inter-dependent the values of the time-series are. 
We implemented several coefficients generalizing the concept of linear correlations to categorical time-series.
Cramer's coefficient, Cohen's coefficient and Theil's U can respectively be computed via the ```cramer_coefficient```, ```cohen_coefficient``` and ```theils_u``` functions.
They take as input a 1-D array representing the time-series to study and an array of lags storing the lag values at which the coefficients are evaluated. For example:
```Julia
cramer_coefficient(series, lags)
```

#### Motif recognition
Time-series can present repeating motifs that are worthwhile identifying. However, simple line-search algorithm are not adapted for all motifs (@Pevzner2000)
Moreover, the lack of proper distance measurement complicates the saerch in the context of categorical time-series.  
An implementation using the *random projection* method (@Buehler2002) is used here.
The identification of potential motifs is done with the ```detect_motifs``` function.
It takes as required input a time-series (1-D array), the length of the motifs to look for, and the number of allowed errors. 
It returns an instance of the ```pattern``` structure which stores properties of the identified motif such as shapes, repetition number and positions.


#### Data clustering.
If certain categories inside a time-series present functional similarities, one might wish to clusther them together into a single equivalent representation.
This reduces the total number of categories and can simplify the analysis of the time-series. To do this, we use an implementation based on the *Information bottleneck* concept (@Tishby2000, @Strouse2016).
After an initial bottleneck model of the structure ```IB``` is instantiated, it can be optimized with the ```IB_optimize!``` function to reveal potential clusters of categories. An overview of the results can be obtained with the ```print_results``` function. 







# Applications

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
