using Combinatorics
using Distributions
using DataStructures
using Random

function big_binomial(n,k)
    bb = factorial(big(n))/(factorial(big(k))*factorial(big(n - k)))
    return floor(bb)
end

"""
    error_dist(x1, x2)

Returns the number of differences between two sorted patterns 'x1' and 'x2'.
"""
function error_dist(x1, x2)
    return sum(x1 .!= x2)
end

"""
Returns dictionnary containing the probability of each symbol found in ts.
"""
function symbol_probabilities(ts)
    # Sample size should be at least 10 times more than the numer of symbols.
    proba = values(counter(ts))
    proba = proba./sum(proba)
end

"""
    match_probability(ts, w, d)

Estimates and returns the expectation value of motif-pair of length 'w'
matching up to 'd' errors in the provided time-series 'ts'.
"""
function expected_matches(ts, w, d)
    p = symbol_probabilities(ts)
    p_symbol_match = p_motif_match(p, 2)
    p_symbol_mismatch = 1 - p_symbol_match
    match_proba = 0
    for i in 0:d
        match_proba += binomial(w, i)*p_symbol_match^(w-i)*p_symbol_mismatch^i
    end
    return big_binomial(length(ts), 2)*match_proba
end

"""
    masked_match(w, d, t)

Returns the probability of two random masked generated motifs of length 't'
to match despite 'd' error. 'w' is the size of the full motifs.
!! this is a probability PER ITERATION !!
inputs (Int):
    w : motif size (window size).
    d : # of allowed errors between motifs.
    t : length of projections after applying mask. Defaults to w - d.
"""
function masked_match(w, d, t)
    match_proba = binomial(w - t, d)/binomial(w, d)
    return match_proba
end

"""
    least_occurence_threshold(w, d, p; confidence = 0.95)

Computes the least number of occurences that a motif has to have in order not to be considered a product of chance.
input:
    ts : input time-series
    w : motif size (window size).
    d : # of allowed errors between motifs.
    confidence (optional) : the level of desired confidence. If  set to 0.5, the expected number of motifs chains matching
        'match_number' times by chance will be 0.5. (defaults to 0.5)
Returns:
    match_number : match number by chance allowed by 'confidence'.
"""
function least_occurence_threshold(ts, w, d; confidence = 0.5)
    p = symbol_probabilities(ts)
    E = []
    max = 20
    for r in 2:max
        p_symbol_match = p_motif_match(p, r)
        p_symbol_mismatch = 1 - p_symbol_match
        match_proba = 0
        for i in 0:d
            match_proba += binomial(w, i)*p_symbol_match^(w-i)*p_symbol_mismatch^i
        end
        push!(E, big_binomial(length(ts), r)*match_proba)
    end
    idxs = findall(x -> x < confidence, E)
    if isempty(idxs)
        @warn "Number of identical match by chance >= $max. \n Try considering longer motifs or fewer allowed errors. \n Assigning value of match by chance to $max"
        return max
    end
    return idxs[1]
end



"""
    p_motif_match(p, match_number)

probability of two motif matching 'match_number' times by chance.
"""
function p_motif_match(p, match_number)
    p_match = copy(p)
    for i in 1:match_number-1
        p_match = p_match .* p
    end
    return sum(p_match)
end
