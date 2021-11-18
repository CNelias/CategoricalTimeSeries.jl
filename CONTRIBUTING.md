## Issues
If you are experiencing troubles/bugs with some functionallities of the package, please open an issue and I will try to resolve it.
Alternatively, you can also open an issue to give me feedback or suggest new features.

## List of improvements
Despite all my efforts to make this package as complete as possible, I haven't had the time to implement certain improvements.
If you want to help, your contribution will be appreciated!

Here is a list of potential improvements to the existing methods:
- Implement common windowing functions for the **spectral envelope** method: The original paper post processes the results with a smoothing triangular function (of width ```m```), however the core estimation is still based on the periodogram. This effectively passes the data through a normalized boxcar function which is known to introduce bias.
As is the case in power-spectral density estimation, (e.g. Welch's method) a carefully choosen window such as the Hanning, Hamming or Blackman function could help reducing the bias and improve the accuracy of the results. An even better improvement would be to use the window functions of the slepian sequence, to allow averaging without frequency range restriction. However, introducing such windows will change the mathematical starting point of the spectral envelope and I do not know if closed-form solutions can be found for all windows.
- Implement simulated annealing for the **Information bottleneck** method: due to the probabilistic nature of this algorithm, results do not always converge to the absolute minima of the cost function. 
At the moment, the package samples n optimizations and selects the one that has the lowest cost. This is a working approach, but can be computationally costly depending on the amount of datapoints and categories. Simulated annealing could be a way to reduce computation time.
- Find a way to display the results of the **Information bottleneck** clusterings in a more visually comprehensive manner. At the moment, a DataFrame filled with 1 and 0s indicating appartenance to different clusters is used,
which can be hard to interpret when many categories are involved.
- add doc entry for ```apply_mapping```` function and correct export from apply_mappings to apply_mapping.

