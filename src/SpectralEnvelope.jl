module SpectralEnvelope

using FFTW
using Polynomials
using LinearAlgebra
using Statistics

function detrend(x,y; deg = 1)
    fit = polyfit(x,y,deg)
    detrended = y - fit(x)
    return detrended
end

"""
smoothens the given series by mixen each points with it's neighboors.

    input:
    - data: the series you want to smooth
    - m : the numbers of points to be involved in the mixing process.

    output:
    - smoothened series, shorter of m points than the original.

"""

function smooth(data::Array{Float64,1}; m = 5)
    if m != 0
        smoothed_data = Float64[]
        weigths = Int[]
        for w in collect(-m:m)
            p = abs(2*w)+1
            if p >= 10
                append!(weigths,p)
            else
                append!(weigths,10)
            end
        end
        p_tot = sum([1/w for w in weigths])
        for i in 1:length(data)
            if i > m && i < length(data) - m
                smoothed_Iw = 0
                for (l,k) in enumerate(collect(-m:m))
                    smoothed_Iw += (1/weigths[l])*data[i+k]
                end
                append!(smoothed_data, smoothed_Iw/p_tot)
            end
        end
        return smoothed_data
    elseif m == 0
        return data
    end
end


"""
One-hot encodes the time-series. For k categories,
each category gets a k-1 unit vector, except the last one which
is mapped to zeros(k - 1).

    Example:

    vectorize([1,2,3,1,2]) returns :
    [[1,0],[0,1],[0,0],[1,0]]

"""

function vectorize(data)
    categories = unique(data)
    deleteat!(categories, length(categories))
    sorted_data = zeros(length(categories),length(data))
    for (t,d) in enumerate(data)
        for (i,j) in enumerate(categories)
            if d == j
                sorted_data[i,t] = 1
            end
        end
    end
    return sorted_data
end


"""
Splits the time-series into overlapping equal-lengthed chunks of given size.

    input :
    - x : the time-series
    - window : the size of the equal-lengthed chunks
    - step : size of the step by which the window is sliding along the time-series
    the smaller the value of size, the bigger the overlap between the different chunks.

    output :
    -  overlapping equal-lengthed windowed time-series
"""

function partitioning(x,window::Int,step::Int)
    return [x[i:i+window] for i in 1:step:length(x) if  i + window <= length(x)]
end

function partitioning(x::Array{Float64,2},window::Int,step::Int)
    return [x[:,i:i+window-1] for i in 1:step:length(x[1,:]) if  i + window <= length(x[1,:])]
end


"""
returns the power spectrum of a given time serie,

    input :
    - time serie
    - window for averaging (default = length(time serie)/10)
    - step for the window's sliding

    output :
    - value of the power spectrum

You will have to provide the frequencies yourself
"""

function power_spectrum(x::Array{Float64,1}, window::Int, step::Int)
    ts = partitioning(x,window,step)
    ps = [real(fft(i .- mean(i)).*conj(fft(i .- mean(i)))) for i in ts]
    pxx = mean(ps)
    return [pxx[i] for i in 1:div(window,2)]
end


"""
    periodogram_matrix(ts::Array{Float64,2}, demean = false; m=2)

Computes the k x k periodogram matrix (k number of categories).
If 'average' true, segments data and computes the average periodogram matrix over all segments.
m is the smoothing parameter.
"""
function periodogram_matrix(ts::Array{Float64,2}, average = false; m=2)
    if average == false
        periodo = zeros(length(ts[:,1]), length(ts[:,1]), length(ts[1,:]))
        DFT = zeros(Complex{Float64}, length(ts[:,1]), length(ts[1,:]))
        for i in 1:length(ts[:,1])
            DFT[i,:] = fft(ts[i,:] .- mean(ts[i,:]))
        end
        for i in 1:length(ts[:,1])
            for j in 1:length(ts[:,1])
                tmp = smooth(real((1/length(ts[1,:]))*DFT[i,:].*conj(DFT[j,:])); m = m)
                for w in 1:length(ts[1,:])-(2*m+1)
                    periodo[i,j,w] = tmp[w]
                end
            end
        end
        return periodo
    elseif average == true
        window = Int(floor(length(ts[1,:])/2))
        step = Int(floor(length(ts[1,:])/4))
        dims = length(ts[:,1])
        periodo = zeros(length(ts[:,1]), length(ts[:,1]), window)
        seg_ts = partitioning(ts, window, step)
        DFT = []
        for el in seg_ts
            tmp = zeros(Complex, dims, window)
            for i in 1:length(el[:,1])
                tmp[i,:] = fft(el[i,:] .- mean(el[i,:]) )
            end
            push!(DFT, tmp)
        end
        for el in DFT
            for i in 1:dims
                for j in 1:dims
                    tmp = smooth(real((1/window)*el[i,:].*conj(el[j,:])); m = m)
                    for w in 1:window-(2*m+1)
                        periodo[i,j,w] += tmp[w]
                    end
                end
            end
        end
        return periodo./length(seg_ts)
    end
end


"""
Computes the covariance-variance matrix of a given time-series.
"""
function varcov(ts::Array{Float64,2})
    cov_matrix = zeros(length(ts[:,1]), length(ts[:,1]))
    for i in 1:length(ts[:,1])
        for j in 1:length(ts[:,1])
            cov_matrix[i,j] = sum((1/(length(ts[1,:])-1))*(ts[i,:].-mean(ts[i,:])).*(ts[j,:].-mean(ts[j,:])))
        end
    end
    return cov_matrix
end


"""
Computes the spectral envelope of the given time-series.

    input:

    - ts: the time series to analyse
    - m : the degree of smoothing wished.
          m corresponds to the number of neighbooring points that are mixed
          with given point to realize the smoothing.

    output:

    - Frequencies : list of points corresponding to the involved frequencies.
                    contained in [0,0.5]
    - amplitude : values of the spectral envelope for each given frequency point.
    - eigenvectors : the optimal scaling for the different categories, for each frequency point.
    - categories : the list of categories present in the data.

"""
function spectral_envelope(ts; m = 3)
    result = Float64[]
    vec = vectorize(ts)
    eigvec = zeros(length(vec[:,1]), trunc(Int,length(ts)/2))
    S = inv(sqrt(varcov(vec)))
    # alternative definition of S
    # S = inv(varcov(vec))
    p = periodogram_matrix(vec; m = m)
    f_len = trunc(Int,length(p[1,1,:])/2)
    for i in 1:f_len
        ev = findmax(real.(eigvals(S*p[:,:,i]*S/f_len)))
        # ev = findmax(real.(eigvals(S*p[:,:,i]/f_len)))
        append!(result, ev[1])
        b = S*eigvecs(S*p[:,:,i])[ev[2],:]
        b = b./sqrt(sum(b.^2))
        eigvec[:,i] = real.(b)
    end
    return collect(1:length(result))/length(p[1,1,:]), result[1:end], eigvec
end

"""
     get_mappings(data, freq; m = 3)

Returns the optimal mappings corresponding to the frequency 'freq'.
Prints the position and strengh of peak at 'f' for control purposes.
input :
    - data : input categorical time series
    - freq : the frequency where one wants the optimal mappings
    - m : smoothing parameters (see 'spectral_envelope')
returns :
    - mappings : the mapping at the peak associated with the desired frequency.

example:
    # we have some data, plotting the spectral envelope, we see a peak at 0.33
    m = get_mappings(data, 0.33)
"""
function get_mappings(data, freq; m = 3)
    categories = unique(data)
    f, se, eigvecs = spectral_envelope(data; m = m)
    delta = f[2] - f[1]
    window = Int(div(0.04*length(f),1))
    peak_pos = findmin([abs(freq - delta*i) for i in 1:length(f)])[2]
    true_pos = findmax(se[peak_pos-window:peak_pos+window])[2] + peak_pos - window - 1
    mappings = [string(categories[i]," : ",round(eigvecs[i,true_pos]; digits = 2)) for i in 1:length(eigvecs[:,1])]
    if length(categories) != length(eigvecs[:,1])
        push!(mappings, string(unique(data)[end]," : ", 0.0))
    end
    println("position of peak: ",round(f[true_pos],digits = 2)," strengh of peak: ",round(se[true_pos], digits = 2))
    return mappings
end

function findmax_in(xserie,yserie,xlim)
    range = []
    real_pos = []
    for (index,value) in enumerate(xserie)
        if value>= xlim[1] && value <= xlim[2]
            append!(range,yserie[index])
            append!(real_pos,index)
        end
    end
    max,pos = findmax(range)
    return max, xserie[real_pos[pos]], real_pos[pos]
end


export spectral_envelope, get_mappings, detrend, smooth, power_spectrum

end
