## List of improvements

Despite all my efforts to make this package as complete as possible, I haven't had the time to implement certain improvements.
If you want to help, your contribution will be appreciated!
Here is a list of the things that could be added to the existing methods:
- Implement windowing for the **spectral envelope** method. The original paper post processes the results with a smoothing function, however the core estimation is still based on the periodogram.
As is the case in power-spectral density estimation, (e.g Welch's method) a carefully choosen window function could help reducing the bias and potentially improve the accuracy of the results.
- Implement simulated annealing for the **Information bottleneck** method: due to the probabilistic nature of this algorithm, results do not always converge to the absolute minima of the cost function. 
At the moment, the package samples n optimizations and selects the one that has the lowest cost function, which works but can be computationally costly. Simulated annealing could be a way to reduce computation time.
- Find a way to display the results of the **Information bottleneck** clusterings in a more visually comprehensive manner. At the moment, a DataFrame filled with 1 and 0s indicating appartenance to different clusters is used,
which can be hard to interpret when many categories are involved.

