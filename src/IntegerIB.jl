module IntegerIB

using Random
using DataFrames
import Base.merge!

"""
    Computes the KL divergence for two 1D probability distributions.
"""
function KL_single(P, Q)
    return sum(P.*log2.(P./Q))
end

"""
    Computes the KL divergence for a collection of distributions in matrix form.
Given P of shape MxN and Q of shape MxL, returns a matrix of shape NxL.
"""
function KL(P, Q)
    kl = zeros(size(P)[2], size(Q)[2])
    for i in 1:size(P)[2]
        for j in 1:size(Q)[2]
            kl[i,j] = KL_single(P[:,i],Q[:,j])
        end
    end
    return kl
end


"""
    Computes the entropy of a distribution. If distribution is of shape (MxN),
treats each column as an independant distribution and returns an 1xN entropy vector.
"""
function entropy(p)
    if any(size(p) .== 1)
        return -sum(p.*log2.(p))
    else
        return -sum(p.*log2.(p), dims = 1)
    end
end


"""
    get_y(x, algorithm = "nn")

Returns the corresponding context vector `y` to the data vector `x`.
Several types of context are possible, and are selected by the `type` argument.
type = "nn" is the next neighbor to every point of `x`.
type = "an" stands for "adjacent neighbor", or the previous and next element of every element of x.
"""
function get_y(x, type = "nn")
    if type == "nn"
        y = x[2:end]
    elseif type == "an"
        y = Any[]
        for idx in collect(2:length(x)-1)
            push!(y, (x[idx-1], x[idx+1]))
        end
    end
    return y
end

function findin(el, array)
    return findall(x -> x == el, array)[1]
end

"""
    get_pxy(x,y)
Returns joint probability of x and y from normalized 2D histogram.
Automatically replaces 0 by 10^-16 to prevent Inf and NaN values during KL computation.
"""
function get_pxy(x,y)
    X = x[1:end-1]
    if length(X) != length(y)
        X = X[2:end]
    end
    count = 0
    x_elts, y_elts = sort(unique(X)), sort(unique(y))
    pxy = zeros(size(x_elts)[1], size(y_elts)[1])
    for (idx, el) in enumerate(X)
        xpos, ypos = findin(el, x_elts), findin(y[idx], y_elts)
        pxy[xpos, ypos] += 1
        count += 1
    end
    pxy[pxy .== 0.0] .= 10^-16 #replacing zeros

    return pxy./count
end

"""
    returns the conditional probability of y given x.
"""
function get_py_x(x, y)
    pyx = get_pxy(x, y)
    inv_px = 1 ./ get_px(x, y)
    cond_pyx = inv_px.*pyx
    return cond_pyx' #transposition in order to get the y dimension in dims = 1.
end

function get_py_x(pxy)
    inv_px = 1 ./ get_px(pxy)
    cond_pyx = inv_px.*pxy
    return cond_pyx' #transposition in order to get the y dimension in dims = 1.
end

"""
    Computes marginal probability of 'x'.
"""
function get_px(x, y)
    return sum(get_pxy(x, y), dims = 2)
end

function get_px(pxy)
    return sum(pxy, dims = 2)
end

"""
    Computes marginal probability of 'y'.
"""
function get_py(x, y)
    return sum(get_pxy(x, y), dims = 1)
end

function get_py(pxy)
    return sum(pxy, dims = 1)
end

function get_qt_x_init(px, algorithm)
    if algorithm == "IB"
        qt_x_init = rand(size(px)[1], size(px)[1])  # randomly initialize qt_x_init
        qt_x_init = qt_x_init.*(1 ./ sum(qt_x_init, dims = 1)) #normalize qt_x_init
        return qt_x_init
    elseif algorithm == "DIB"
        index = shuffle(collect(1:size(px)[1]))
        qt_x_init = zeros(size(px)[1], size(px)[1])
        for (i,v) in enumerate(index)
            qt_x_init[v,i] = 1
        end
        return replace_zeros(qt_x_init)
    else
        throw("Algorithm must be 'IB' or 'DIB', but $algorithm was given.")
    end
end

"""
    init_values(x, y)

Initializes and return clustering porbability matrices `qt_x`, `q` and `qy_t`.
'qt_x' refers to the probability of a category in x to be mapped to an element of t.
'q' is the overall probability of an element of t to occur. 'qy_t' is the probability that an element of t
corresponds to a given element of `y`.
"""
function init_values(x, y, algorithm = "IB")
    pxy = get_pxy(x, y)
    px = get_px(x, y)
    qt_x_init = get_qt_x_init(px, algorithm)
    qt_init = qt_x_init*px
    qy_t_init = pxy'*qt_x_init'.*(1 ./ qt_init)'
    return (qt_x_init, qt_init, qy_t_init)
end

function init_values(pxy, algorithm = "IB")
    px = get_px(pxy)
    qt_x_init = get_qt_x_init(px, algorithm)
    qt_init = qt_x_init*px
    qy_t_init = pxy'*qt_x_init'.*(1 ./ qt_init)'
    return (qt_x_init, qt_init, qy_t_init)
end

"""
    mapping(input)

Maps an input array of any type (string, float, any...) to an array of Int usable by the algorithm.
The structure of the array is preserved as the mapping is one-to-one.
"""
function mapping(input)
    categories = unique(input)
    mapped_series = Int64[]
    for value in input
        for (idx, ctg) in enumerate(categories)
            if ctg == value
                push!(mapped_series, idx)
            end
        end
    end
    mappingDict = Dict(unique(mapped_series) .=> categories)
    return mapped_series, mappingDict
end

"""
Interface for IB clustering computation.
Input data 'x' must be 1D. β controls the amount of clustering : the smaller the beta, the greater the clustering (and information loss).
Can be initalized directly with the data via IB(x), or if you already have the histograms/probability of co-occurences of 'x' and 'y' via IB(pxy).
To use the deterministic IB algorithm, set algorithm to "DIB" like IB(x, "DIB").
"""
mutable struct IB
    algorithm::String
    β::Real
    x
    px
    py
    pxy
    py_x
    colDict #dict from original value to mapped values. used for representation in print_results.
    qt_x
    qt
    qy_t

    IB(x, y, β = 100, algorithm::String = "IB") = new(algorithm, β, x, get_px(x,y), get_py(x,y), get_pxy(x,y), get_py_x(x,y), nothing, init_values(x, y, algorithm)...)
    IB(x::Array{Float64,1}, β = 100, algorithm::String = "IB") = new(algorithm, β, x, get_px(x, get_y(x)), get_py(x, get_y(x)), get_pxy(x, get_y(x)), get_py_x(x, get_y(x)), nothing, init_values(x,  get_y(x), algorithm)...)
    IB(x::Array{Any,1}, β = 100, algorithm::String = "IB") = new(algorithm, β, x, get_px(x, get_y(x)), get_py(x, get_y(x)), get_pxy(x, get_y(x)), get_py_x(x, get_y(x)), nothing, init_values(x,  get_y(x), algorithm)...)
    IB(x::Array{Int64,1}, β = 100, algorithm::String = "IB") = new(algorithm, β, x, get_px(x, get_y(x)), get_py(x, get_y(x)), get_pxy(x, get_y(x)), get_py_x(x, get_y(x)), nothing, init_values(x,  get_y(x), algorithm)...)
    IB(pxy::Array{Float64,2}, β = 100, algorithm::String = "IB") = new(algorithm, β, nothing, get_px(pxy), get_py(pxy), pxy, get_py_x(pxy), nothing, init_values(pxy, algorithm)...)
    IB(x::Array{String,1}, β = 100, algorithm::String = "IB") = new(algorithm, β, mapped_init_values(x, algorithm)...)
end


"""
    Returns all initialization values for the IB struct in cases where input is not a real valued array.
"""
function mapped_init_values(input, algorithm)
    mapped_x, mapping_dict = mapping(input)
    px = get_px(mapped_x, get_y(mapped_x))
    py = get_py(mapped_x, get_y(mapped_x))
    pxy = get_pxy(mapped_x, get_y(mapped_x))
    py_x = get_py_x(mapped_x, get_y(mapped_x))
    qt_x, qt, qy_t = init_values(pxy, algorithm)
    return (mapped_x, px, py, pxy, py_x, mapping_dict, qt_x, qt, qy_t)
end

"""
    Re-initialies values of IB model
"""
function init!(m::IB)
    m.qt_x, m.qt, m.qy_t = init_values(m.pxy, m.algorithm)
end

"""
    Peforms q(t|x) update step.
"""
function qt_x_step!(m::IB)
    kl_part = log.(m.qt)' .- m.β*KL(m.py_x, m.qy_t)
    if m.algorithm == "DIB"
        m.qt_x = zeros(size(m.qt_x))
        for x_index in 1:size(kl_part)[1]
            m.qt_x[argmax(kl_part[x_index,:]), x_index] = 1
        end
        m.qt_x = replace_zeros(m.qt_x)
    elseif m.algorithm == "IB"
        q = replace_zeros(exp.(kl_part))
        m.qt_x = q'.*(1 ./ sum(q', dims = 1)) #transpose and normalize q.
    else
        throw(ArgumentError("'algorithm' unknown must be either \"IB\" or \"DIB\""))
    end
end

"""
    Peforms q(t) update step and drops unused clusters.
"""
function qt_step!(m::IB, clust_thres; report = true)
    m.qt = m.qt_x*m.px
    dropped = m.qt .< clust_thres
    if any(dropped) && size(m.qt)[1] > 1
        if report
            println("$(length(dropped[dropped])) clusters dropped, $(length(dropped[.~dropped])) ramaining.")
        end
        m.qt = m.qt[.~dropped]
        m.qt_x = drop_rows(m.qt_x, dropped)
        m.qt_x = m.qt_x.*(1 ./ sum(m.qt_x, dims = 1)) #normalize cluster distribution again
        m.qt = m.qt_x*m.px
    end
end

"""
    Remove rows of 2D array 'value_arr'. Boolean array 'drop_arr' indicates which rows to drop.
"""
function drop_rows(value_arr, drop_arr)
    new_size = length(drop_arr[.~drop_arr])
    dropped_arr = zeros(new_size, size(value_arr)[2])
    count = 1
    for (idx, drop) in enumerate(drop_arr)
        if ~drop
            dropped_arr[count,:] = value_arr[idx,:]
            count += 1
        end
    end
    return dropped_arr
end


"""
    Performs q(y|t) update step.
"""
function qy_t_step!(m::IB)
    qt_y = m.qt_x*m.pxy
    m.qy_t = qt_y'.*(1 ./ m.qt')
end

"""
    replace_zeros(arr::Array{Float64,2}, cut = 10^-305)

Replaces zeros (values below 'cut') by 'cut' to avoid exponential term causing NaN at high beta values.
"""
function replace_zeros(arr::Array{Float64,2}, cut = 10^-305)
    arr[arr .< cut] .= cut
    return arr
end


"""
    Performs one iteration of the IB algorithm.
"""
function IB_step!(m::IB, clust_thres; report = true)
    qt_x_step!(m)
    qt_step!(m, clust_thres; report = report)
    qy_t_step!(m)
end

"""
    Checks if absolute or relative convergence criterion on L function is satisfied.
"""
function check_conv(L_prev, L_current, abs_thres, rel_thres)
    abs_change = abs(L_prev - L_current)
    rel_change = abs_change/L_prev
    if abs_change <= abs_thres || rel_change <= rel_thres
        return true
    else
        return false
    end
end

function calc_metrics(m::IB)
    ht = entropy(m.qt)
    hy_t = entropy(m.qy_t)*m.qt
    iyt = entropy(m.py) - hy_t[1] # [1] cause of dumb julia rules focing dot product of vector to be 2 dimensional.
    ht_x = entropy(m.qt_x)*m.px
    ixt = entropy(m.qt) - ht_x[1] # same as iyt.
    L = ht - ht_x[1] - m.β*iyt    # same
    return (ht, ixt, iyt, L)
end

"""
    IB_optimize(model::IB; abs_thres = 10^-3, rel_thres = 0.01, ct_thres = 5, n_conv = 1000)

Optimizes IB model of the data.
Keywords arguments are :
    - abs_thres : threshold for convergence in absoulte change of the optimization function (lagrangian) L.
    - rel_thres : threshold for convergence in relative change of L.
    - ct_thres : number of steps where L has to be < to abs_thres or rel_thres to declare convergence.
    - n_conv : limit of the number of iterations for the optimization algorithm.
    - scan (bool): when using 'DIB' algorithm, checks and merges clusters if it leads to reduction of cost function.
    - report : wether or not to say when and how convergence was reached.
"""
function IB_optimize!(model::IB; merge_thres = 5*10^-2, abs_thres = 10^-10, rel_thres = 0, ct_thres = 25, n_conv = 1000, scan = true, report = true)
    consec_conv = 0
    total_steps = 0
    L = calc_metrics(model)[4]
    while consec_conv <= ct_thres && total_steps <= n_conv
        IB_step!(model, merge_thres; report = report)
        if model.algorithm == "DIB" && scan
            brute_optimize!(model, false)
        end
        L_current = calc_metrics(model)[4]
        if check_conv(L, L_current, abs_thres, rel_thres)
            consec_conv += 1
        else
            consec_conv = 0
        end
        L = L_current
        total_steps += 1
        if total_steps == n_conv && report
            println("Convergence not reached after $n_conv iterations, stopping optimization.")
        elseif consec_conv == ct_thres && report
            println("Convergence reached after $total_steps iterations, stopping optimization")
        end
    end
end

"""
    scans through the qt probabilities and merges clusters where values are too low.
"""
function merge!(m::IB, clust_thres = 10^-5)
    dropped = m.qt .< clust_thres
    if any(dropped) && size(m.qt)[1] > 1
        m.qt = m.qt[.~dropped]
        m.qt_x = drop_rows(m.qt_x, dropped)
        m.qt_x = m.qt_x.*(1 ./ sum(m.qt_x, dims = 1)) #normalize cluster distribution again
        m.qt = m.qt_x*m.px
        m.qy_t = drop_rows(m.qy_t', dropped)'
    end
end

"""
    brute_optimize!(m::IB, findbest = true)

Performs bute force optimization of IB model, looking at all possible cluster combinations
and choosing the merges that reduce the cost function.
"""
function brute_optimize!(m::IB, findbest = true)
    if m.algorithm == "IB"
        println("IB agorithm must be DIB.")
        return
    end
    if size(m.qt_x)[1] == 1
        return
    end
    anybetter = false
    best = deepcopy(m)
    for t1 in range(1, stop = size(m.qt_x)[1]-1)
        for t2 in range(2, stop = size(m.qt_x)[1])
            if ~anybetter || findbest
                tmp = deepcopy(best)
                tmp.qt_x[t1, findall(x -> x == 1, m.qt_x[t2,:])] .= 1
                tmp.qt_x[t2,:] .= 0
                IB_step!(tmp, 0)
                # check if cost function L reduced relative to best so far
                if calc_metrics(tmp)[4] < calc_metrics(best)[4]
                    best = deepcopy(tmp)
                    mergedt1 = t1
                    mergedt2 = t2
                    anybetter = true
                    print("merged clusters $mergedt1 and $mergedt2 \n")
                end
            end
        end
    end
    merge!(best)
    if anybetter
        m.qt_x = best.qt_x
        m.qt = best.qt
        m.qy_t = best.qy_t
        return true
    else
        return false
    end
end

"""
    Scans the IB plane with various vlaues of beta to get the optimal curve in the IB plane.
"""
function get_IB_curve(m::IB, start = 0.1, stop = 400, step = 0.05; glob = false)
    beta = collect(start:step:stop)
    ixt, iyt = [], []
    for (idx, b) in enumerate(beta)
        m.β = b
        init!(m)
        if glob
            search_optima!(m)
        else
            IB_optimize!(m; report = false)
        end
        x, y = calc_metrics(m)[2:3]
        if idx == 1
            push!(ixt, x)
            push!(iyt, y)
        elseif y >= iyt[end]
            push!(ixt, x)
            push!(iyt, y)
        end
    end
    return ixt, iyt
end

"""
    search_optima!(m::IB, n_iter = 10000)

The original IB algorithm can converge to a local optima. To find the global optima,
search_optima! performs 'n_iter' optimization of input model 'm', randomly re-initializing it at every iteration.
The best optimization is kept and 'm' is modified in-place.
"""
function search_optima!(m::IB, n_iter = 10000)
    if m.algorithm == "DIB"
        println("DIB algorithm is already giving optima.")
        return
    end
    best = deepcopy(m)
    for i in 1:n_iter
        init!(m)
        IB_optimize!(m; report = false)
        if calc_metrics(m)[4] < calc_metrics(best)[4]
            best.qt_x = m.qt_x
            best.qt = m.qt
            best.qy_t = m.qy_t
        end
    end
    m.qt_x = best.qt_x
    m.qt = best.qt
    m.qy_t = best.qy_t
end

"""
    Helper function that puts all similar clusters next to each other for ease of readability.
"""
function group_equivalent!(df)
    for i = 1:size(df,2)-1
        for j = i+1:size(df,2)
            if df[!,i] == df[!,j]
                idxs = collect(1:size(df,2))
                idxs[j], idxs[i+1] = idxs[i+1], idxs[j]
                select!(df, idxs)
                break
            end
        end
    end
end

"""
    print_results(m::IB, disp_thres = 0.1)

Displays the clustering of data for an optimized input model 'm'.

Since the IB algorithm is not deterministic, some probabilities in the clustering matrix 'm.qt_x' are non-zero.
For ease of readability, every probability under 'disp_thres' is displayed as 0, everything else as 1.
The result is a 2D matrix where the rows represent the clusters and the column the initial categories.
The 1 tell to which cluster which category belongs.
If several 1 are present for a single category, it means that the optimization wasn't able to assign
a single cluster to this category for the chosen β value. You might try a diffent β, or use the DIB variant of the algorithm.
"""
function print_results(m::IB, disp_thres = 0.1)
    if ~isnothing(m.x)
        df = isnothing(m.colDict) ? DataFrame(m.qt_x .> disp_thres, Symbol.(sort(unique(m.x)))) : DataFrame(m.qt_x .> disp_thres, [Symbol.(m.colDict[i]) for i in sort(unique(m.x))])
        group_equivalent!(df)
        display(convert.(Int,df))
    else
        @warn "Initial data is probability distribution, categories will be displayed as x1, x2 ..., xn."
        df = DataFrame(m.qt_x .> disp_thres, :auto)
        group_equivalent!(df)
        display(convert.(Int,df))
    end
end

export IB, search_optima!, brute_optimize!, IB_optimize!, calc_metrics, get_IB_curve, get_y, print_results, mapping

end
