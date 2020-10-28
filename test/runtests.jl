using Test
using DelimitedFiles
using CategoricalTimeSeries

#testing the serial dependences functions by reproducing the results from C. Weiss's book "An Introduction to Discrete-Valued Time Series".
#I was not able to get my hands the full time-series, C. weiss says it is ~1300 elements long, which explains the slight numerical differences in values.
pewee = Int64.(readdlm("pewee.txt", ',')[1,:])
@test round(cramer_coefficient(pewee, 4)[1], digits = 2) == 0.46
@test round(cohen_coefficient(pewee, 4)[1], digits = 2) == 0.55
#this one was not in the book, but the plot correctly reproduces the results.
@test round(theils_u(pewee,[1])[1], digits = 2) == 0.55
@test round(H(pewee), digits = 2) == 1.49

#testing spectral envelope method by reproducing David's stoffer results.
DNA = readdlm("DNA_data.txt")
x,y,e = spectral_envelope(DNA; m=0)
@test round(spectral_envelope(DNA)[2][5]; digits = 3) == round(0.175; digits = 3)

#testing integerIB
#bach chorale data, taken from Nori Jacoby's webpage.
bach_histogram = [0.00001 392 14 582 896 262 106; 12 0.0001 15 5 674 17 182; 13 6 0.0001 44 12 32 3; 239 134 16 0.0001 380 13 173; 2099 25 46 93 0.00001 211 2; 96 206 25	95	172	0.00001 36; 499 1 4 9 24 20 0.00001]
pxy = bach_histogram./sum(bach_histogram)
model = IB(pxy, 0.95, "DIB")
brute_optimize!(model)
h, ixt, iyt, L = calc_metrics(model)
@test size(model.qt_x, 1) == 5
@test round(h, digits = 2) == 2.06
@test round(ixt, digits = 2) == 2.06
@test round(iyt, digits = 2) == 0.78
@test round(L, digits = 2) == 1.32

