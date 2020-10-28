module CategoricalTimeSeries

using IntegerIB
using SerialDependence
using SpectralEnvelope

export IB, search_optima!, brute_optimize!, IB_optimize!, calc_metrics, get_IB_curve, get_y, print_results
export cramer_coefficient, cohen_coefficient, conditional_entropy, LaggedBivariateProbability, H, theils_u, rate_evolution
export spectral_envelope, get_mappings, detrend, smooth, power_spectrum

end
