#Spectral Envelope

The **spectral envelope** is a tool to study cyclic behaviors in categorical data. It is more informative than the traditional approach of attributing a different number to each category for power-spectral density estimation. <br/>

For each frequency in the spectrum, the **spectral envelope** finds an optimal real-numbered mapping that maximizes the normed power-spectral density at this point. Therefore, no matter what mapping is choosen for the different categories, the power-spectral density will always be bounded by the spectral envelope.

The spectral envelope was defined by David S. Stoffer in *DAVID S. STOFFER, DAVID E. TYLER, ANDREW J. MCDOUGALL, Spectral analysis for categorical time series: Scaling and the spectral envelope*.\

## Main functions
```Julia
spectral_envelope(ts; m = 3)
```
Computes the spectral envelope of an input categorical time-series.

The degree of smoothing can be chosen by the user.

**Parameters**:

* ts: Array containing input categorical time-series.
* m: Smoothing parameter. corresponds to how many neighboring points
        are to be involved in the smoothing (weighted average). Defaults to 3.
**Returns**:

* freq: Array containing the frequencies of the power-spectrum (or spectral envelope).
* se: Values of the spectral envelope for each frequency in 'freq'.
* eigvec: Array containing the optimal real-valued mapping for each frequency point.

To use the spectral envelope, call the function ```spectral_envelope```, you can then easily plot the results and extract the mapping for a given frequency.
Here is an example with DNA data from a portion of the Epstein virus:
```Julia
using DelimitedFiles, Plots
# extracting data
data = readdlm("..\\test\\DNA_data.txt")
# spectral envelope analysis
f, se, eigvecs = spectral_envelope(data; m = 4)
# plotting the results
plot(f, se, xlabel = "Frequency", ylabel = "Intensity", title = "test data: extract of Epstein virus DNA", label = "spectral envelope")
```
<img src=https://user-images.githubusercontent.com/34754896/91556982-eef72680-e933-11ea-85f3-fab6aea17258.PNG width = "600">


To get the **optimal mappings** for a given frequency, you can use the ```get_mapping(data, freq; m = 3)```. With the previous DNA example, we see a peak at 0.33. To get the corresponding mappings:
```Julia
mappings = get_mappings(data, 0.33)
>> position of peak: 0.33 strengh of peak: 0.6
print(mappings)
>> ["A : 0.54", "G : 0.62", "T : -0.57", "C : 0.0"]
```

The function scans the vincinity of the provided goal frequency and returns the mapping for the found maxima. It also prints the positions and intensity of the peak so that you may control that you actually identified the desired peak and not a nearby sub-peak.<br/>
The codons A and G have a similar mapping, so they could potentially have similar functions : this is however not a necessity, as the spectral envelope only seeks to maximize the power-spectrum. If you want to study equivalency of categories, you should also check the results with a clustering algorithm like https://github.com/johncwok/IntegerIB.jl.git.

## Example
