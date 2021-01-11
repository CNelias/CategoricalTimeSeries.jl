# Motif recognition

Time-series sometimes present **repeating motifs** (or patterns) that are worthwhile identifying. The detection of such motifs can be difficult depending on the amount of noise in the time-series. <br/>

In the case of categorical time-series, the lack of proper metric to measure distance between motifs can make their detection tricky. Improper distances like the number of differences between the two motifs is commonly used.

This package proposes a detection algorithm based on JEREMY BUHLER and MARTIN TOMPA's paper "*Finding Motifs Using Random Projections*". This algorithm although very precise is not exact. Therefore, when you are done detecting potential motifs with the `detect_motifs` function, you can refine your results with `find_motifs` for an exact search.
<br/> The main functions return instances of a class called **pattern**:
- - -
**pattern — Class**
- - -
A class storing useful information about found motifs in a time-series. An array of `pattern` instances is returned when the searching algorithm is done running.
>**Attributes**:

- **shape** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): Array containing the shape (or contour) of the first found repetition of the motif.
- **repetitions** ([Array{Array{Any,1},1}](https://docs.julialang.org/en/v1/base/arrays/)): all the different shapes from the motif's repetitions, they can vary a bit from one to the next.
- **positions** ([Array{Int,1}](https://docs.julialang.org/en/v1/base/arrays/)): the positions at which the different repetitions of the motif were found.

## Main functions
- - -
**detect_motifs — Function**
- - -
```Julia
detect_motifs(ts, w, d, t = w - d; iters = 1000, tolerance = 0.95)
```
Detects all motifs of length 'w' occuring more often than chance, being identical to each other up to 'd' differences inside of imput time-series 'ts'.
Returns an array of `pattern`, inside of which the patterns are classified by how frequently they are observed. The first elements is therefore the most frequently observed motif, and so on.
> **Parameters**:

>>* **ts** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): input time-series in which motifs are searched for.
>>* **w** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): length of motifs to look for.
>>* **d** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): allowed errors (differences) between motifs repetitions.
>>* **t = w - d** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): size of the masks to use for random projection in the detection (defaults to w - d).
>>* **iters = 1000** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): the numbers of iterations for the random projection process (defaults to 1000)
>>* **tolerance = 0.95** ([Float](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): threshold of motif identification. Defaults to 0.95

> **Returns** :
>>* **motifs** : array of `pattern` instances classified by frequency. motifs[1] is therefore the most frequent motif, motifs[2] the second most observed and so on.

- - -
**find_motifs — Function**
- - -
```Julia
find_motifs(ts, shape, d)
```
Given a motif of shape 'shape' (array{any,1}), looks for all the repetitions of it which differ only up to 'd' differences inside of the input time-series 'ts'.
Input:

>**Parameters**:
>>* **ts** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)) : time-series in which to look for motifs
>>* **shape** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): shape (aray{any,1}) of the motif to look for.
>>* **d** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): allowed errors (differences) between motifs

>* **Returns** :

>>* **motif** : an instance of `pattern` containing the found repetition of the input 'shape'.

## Plotting
To help visualize results, two simple plotting functions are provided.
- - -
**plot_motif — Function**
- - -
```Julia
plot_motif(m::pattern)
```
Plots all repetitions of an input `pattern` instance on top of each other to see how similar they are to each other.
> **Parameters**:

>>* **m** : Instance of the `pattern` class

- - -
**plot_motif — Function**
- - -
```Julia
plot_motif(m::pattern, ts)
```
Plots all repetitions of an input `pattern` instance on top of the input time-series 'ts' to better visualize their repartition in time.
> **Parameters**:

>>* **m** : Instance of the `pattern` class
>>* **ts** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): Input time-series


## Example
In an improvisation from John coltrane [segment of DNA](https://github.com/johncwok/SpectralEnvelope.jl/tree/master/test), we extract a time-series of intervals.
A spectral envelope analysis reveal a peak at period ~8, so we expect to find motifs of length ~8 and allow for 1 error between them.
After detection, we visualize the most frequent motif:
```Julia
using DelimitedFiles

data = readdlm("..\\test\\countdown.txt")
motifs = detect_motifs(data, 8, 1)
most_frequent_motif = motifs[1]
plot_motif(most_frequent_motif)
```

<img src=https://user-images.githubusercontent.com/34754896/91556982-eef72680-e933-11ea-85f3-fab6aea17258.PNG width = "500">

We notice that the motif `[.....]` seems to be the underlying shape so we do an exact search, and plot the found motif on top of each other:
```Julia
true_shape = [...]
motif = find_motifs(true_shape, data)
plot_motif(motif)
```
Now, we visualize the repetitions of the motif in the time-series:
```Julia
plot_motif(motif, ts)
```
