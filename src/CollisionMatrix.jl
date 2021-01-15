using Combinatorics
using DataStructures
using Random

include("Thresholds.jl")

"""
    S_matrix(ts, window)

Returns S matrix containing all segments of length 'window' in 'ts'.
Done with a sliding window.
"""
function S_matrix(ts, window)
    S = Array{typeof(ts[1]),2}(undef, length(ts) - window, window)
    for index in 1:length(ts) - window
        S[index,:] = ts[index:index+window-1]
    end
    return S
end

"""
Returns all possible masks of length d among the possible w positions.
positive masks: the indexes represent the columns that are KEPT.
    w : motif length
    t : length of masked motifs
"""
function get_masks(w, t)
    if t > w
        error("mask length must be smaller than w")
    end
    return collect(combinations(collect(1:w), t))
end

"""
    exclusion_mask!(c_matrix, exclusion_zone)

Applies exclusion zone to collision matrix (in place).
"""
function exclusion_mask!(c_matrix, exclusion_zone)
    for column in 1:size(c_matrix)[2] - exclusion_zone
        for row in 1:column + exclusion_zone
            c_matrix[row, column] = -Inf
        end
    end
    for i in 0:exclusion_zone
        c_matrix[:, end - i] .= -Inf
    end
end

"""
A version of `findall` adapted for loops to go faster.
"""
function findall_fast(f, a::Array{T, N}, stocker) where {T, N}
    j = 1
    @inbounds for i in eachindex(a)
        @inbounds if f(a[i])
        @inbounds    stocker[j] = i
            j += 1
        end
    end
    return stocker[1:j-1]
end

"""
    update_c_matrix!(collision_matrix, masked_S)

Updates the collision matrix by adding 1 in the row where repetitions of a motif are found.
The row are the index of the repetitions, the column represent the first found motif of this type.
An exclusion zone is applied to tackle trivial neighbours motifs.
"""
function update_c_matrix!(collision_matrix, masked_S)
    masked_S_list = [masked_S[index,:] for index in 1:size(masked_S)[1]]
    total_count = counter(masked_S_list)
    # sub_count = Dict(key => total_count[key] for key in keys(total_count) if total_count[key] > 1)
    stocker = Vector{Int}(undef, length(masked_S_list))
    for motif in keys(total_count)
        if total_count[motif] > 1
            p = findall_fast(x -> x == motif, masked_S_list, stocker)
            # stocker .= 0
            for index in 1:length(p)-1
                @inbounds if abs(p[index+1] - p[index]) >= length(masked_S[1,:])
                @inbounds collision_matrix[p[index+1], p[1]] += 1
                end
            end
        end
    end
end



"""
    collision_matrix(ts, w, d)

Constructs and returns the collision matrix of time-series 'ts'.
inputs (Int):
    ts : input time-series
    w : motif size (window size)
    d : # of allowed errors between motifs
    e : exclusion zone to get rid of trivial matches.
    t : length of projections after applying mask. Defaults to w - d.
    iters : the number of cycles used to construct the collision matrix.
returns (Int):
    C : collision matrix
"""
function collision_matrix(ts, w, d, t = w - d, e = div(w,2); iters = 1000)
    S = S_matrix(ts, w)
    # pre-alocate memory
    C = zeros(length(ts) - w, length(ts) - w)
    exclusion_mask!(C, e) #apply exclusion zone to collision matrix
    masks = get_masks(w, t)
    for i in 1:iters
        mask = rand(masks)
        update_c_matrix!(C, S[:, mask]) #selecting masked version of S
    end
    return C
end
