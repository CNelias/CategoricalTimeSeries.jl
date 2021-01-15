module MotifRecognition

using Plots


include("CollisionMatrix.jl")
include("Thresholds.jl")


"""
    sort_decreasing(arr)

Returns the indices needed to sort input array in decreasing order.
"""
function sort_decreasing(arr)
    liste = copy(arr)
    sorted_indexes = Int.(zeros(length(arr)))
    for i in 1:length(liste)
        max_idx = findfirst(isequal(maximum(liste)),liste)
        sorted_indexes[i] = max_idx
        liste[max_idx] = 0
    end
    return sorted_indexes
end

"""
    A class whose instances contain all usefull informations about a given motif
found in a categorical time-series.
"""
struct pattern
    shape
    instances
    positions
end

"""
    detect_motifs(ts, w, d, t = w - d; iters = 1000, tolerance = 0.95)

Detects all motifs of length 'w' occuring more often than chance, being identical up to 'd' differences.
Input:

    w : length of motifs to look for
    d : allowed errors (differences) between motifs
    t : size of the masks to use for random projection in the detection
    iters : the numbers of iterations for the random projection process (defaults to 1000)
    tolerance : threshold of motif identification.

returns :
    motifs : a list of motif instances containing all usefull informations about the found motifs (positions, frequency, shapes ...)
"""
function detect_motifs(ts, w, d, t = w - d; iters = 800, tolerance = 0.7)
    CM = collision_matrix(ts, w, d, t; iters = iters)
    min_occurences = least_occurence_threshold(ts, w, d)
    # println("Expected similar matches by chance : $min_occurences")
    expected_matches = iters * masked_match(w, d, t) * tolerance
    occurences = sum(CM .> expected_matches, dims = 1) .+ 1
    idx = sort_decreasing(occurences[1,:]) #the indices needed to sort the final dict and list
    sorted_indices = [idx[i] for i in 1:length(idx)-1 if abs(idx[i+1] - idx[i]) > w] # removing trivial matches
    if isempty(occurences[occurences .>= min_occurences])
        println("No motif found more frequently than chance.")
    end
    list = []
    positions = []
    for j in sorted_indices
        tmp_list = []
        tmp_positions = []
        for i in 1:size(CM)[1]-div(w,2)
            if CM[i,j] > expected_matches
                push!(tmp_list, ts[i:i+w-1]); push!(tmp_positions, i)
                CM[i+1:i+div(w,2), j+1:j+div(w,2)] .= 0 # removing trivial matches
            end
        end
        if isempty(tmp_list) == false
            push!(tmp_list, ts[j:j+w-1]); push!(list, tmp_list)
            push!(tmp_positions, j); push!(positions, tmp_positions)
        end
    end
    motifs = [pattern(list[i][end], list[i], positions[i]) for i in 1:length(list)]
    return motifs
end

"""
    find_motifs(ts, w, d)

Given a motif of shape 'shape' (array{any,1}), look for all the repetitions of it which differ only up to 'd' differences.
Input:

    ts : time-series in which to look for motifs
    shape : shape (aray{any,1}) of the motif to look for.
    d : allowed errors (differences) between motifs

returns :
    motif : an instance of 'pattern' containing the found repetition of the input 'shape'.

"""
function find_motifs(ts, shape, d)
    w = length(shape)
    ts_copy = copy(ts)
    found_motifs = []
    positions = []
    for err in 0:d
        for i in 1:length(ts)-w
            if error_dist(ts_copy[i:i+w-1], shape) <= err
                push!(found_motifs, ts_copy[i:i+w-1]); push!(positions, i)
                ts_copy[i:i+w-1] .= 0
            end
        end
    end
    return pattern(shape, found_motifs, positions)
end

"""
    plot_motif(motif, ts; subplots = false)

Plots all instances of the given motif. If the corresponding time series 'ts'
is provided, plots them on top of it, preserving the time-orderding. Otherwise,
plots all instances of 'motif' on top of each other to facilitate their comparison.
if subplots is given false, will only plot the motifs on top of the time-series.
"""
function plot_motif(motif, ts; subplots = true)
    len = length(motif.shape)
    p2 = plot(1:length(ts), ts, color = "grey", xlabel = "Time", ylabel = "Value", label = "", title = "motif positions")
    for (idx, p) in enumerate(motif.positions)
        plot!(p2, p:p+len-1, ts[p:p+len-1], lw = 3, label = "")
    end
    if subplots
        p1 = plot_motif(motif)
        a = plot(p1, p2, layout = (2,1))
    else
        a = p2
    end
    display(a)
    return a
end

function plot_motif(motif)
    a = plot(motif.shape, lw = 3, yticks = minimum(motif.shape):maximum(motif.shape), label = "#1", title = "Occurences of input motif")
    for (idx, m) in enumerate(motif.instances[1:end-1])
        idx += 1
        plot!(a, m, lw = 3, label = "#$idx")
    end
    display(a)
    return a
end

"""
    mapping(ts)
returns an integer mapping of input time-series 'ts'.
"""
function mapping(ts)
    categories = unique(ts)
    arr = copy(ts)
    for (idx, elt) in arr
        for (pos, c) in categories
            if elt == c
                arr[idx] = pos
            end
        end
    end
    return arr
end

 export detect_motifs, mapping, find_motifs, plot_motif, expected_matches, least_occurence_threshold


end
