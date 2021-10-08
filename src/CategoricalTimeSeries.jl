module CategoricalTimeSeries

include("IntegerIB.jl")
include("SerialDependence.jl")
include("SpectralEnvelope.jl")
include("MotifRecognition.jl")

export IB, search_optima!, brute_optimize!, IB_optimize!, calc_metrics, get_IB_curve, get_y, print_results
export cramer_coefficient, cohen_coefficient, conditional_entropy, LaggedBivariateProbability, H, theils_u, rate_evolution, bootstrap_CI
export spectral_envelope, get_mappings, detrend, smooth, power_spectrum, apply_mappings
export detect_motifs, mapping, find_motifs, plot_motif, expected_matches, least_occurence_threshold

end
