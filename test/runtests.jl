using Test
using DelimitedFiles
using CategoricalTimeSeries

#testing the serial dependences functions by reproducing the results from C. Weiss's book "An Introduction to Discrete-Valued Time Series".
#I was not able to get my hands the full time-series, C. weiss says it is ~1300 elements long, which explains the slight numerical differences in values.
@testset "AssociationMeasurement" begin
  pewee_path = joinpath(dirname(dirname(pathof(CategoricalTimeSeries))), "test", "pewee.txt")
  pewee = Int64.(readdlm(pewee_path, ',')[1,:])
  @test round(cramer_coefficient(pewee, Int64(4))[1], digits = 2) == 0.46
  @test round(cohen_coefficient(pewee, Int64(4))[1], digits = 2) == 0.55
  #this one was not in the book, but the plot correctly reproduces the results.
  @test round(theils_u(pewee,[1])[1], digits = 2) == 0.55
  @test round(H(pewee), digits = 2) == 1.49
end

@testset "SpectralEnvelope" begin
  #testing spectral envelope method by reproducing David's stoffer results.
  DNA_path = joinpath(dirname(dirname(pathof(CategoricalTimeSeries))), "test", "DNA_data.txt")
  DNA = readdlm("DNA_data.txt", ',')[1,:]
  x,y,e = spectral_envelope(DNA; m=0)
  @test round(spectral_envelope(DNA)[2][5]; digits = 3) == round(0.003; digits = 3)
  mappings = get_mappings(DNA, 0.33)
  @test round(mappings["A"]; digits = 2) == round(0.54; digits = 2)
end

#testing integerIB
#bach chorale data, taken from Nori Jacoby's webpage.
@testset "integerIB" begin
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
  #toy time-series presenting one predictable pattern "b" -> "a"
  test_ts = ["a", "b", "a", "b", "a", "c", "d", "b", "a", "a", "d", "b", "a", "b", "a", "d", "c", "d", "b", "a", "d", "a", "c", "b", "a", "a", "b", "a", "c", "b", "a"]
  model = IB(test_ts, 500)
  IB_optimize!(model)
  @test (size(model.qt_x, 1) == 2 || size(model.qt_x, 1) == 3) #we expect the algorithm to cluster at least 2 labels together.
end

#testing motif recognition
@testset "MotifRecognition" begin
  confirmation_path = joinpath(dirname(dirname(pathof(CategoricalTimeSeries))), "test", "confirmation")
  data = readdlm("confirmation")
  pitch = mod.(data, 12)
  intervals = pitch[2:end] .- pitch[1:end-1]
  m = detect_motifs(intervals, 7, 1; iters = 700, tolerance = 0.7)
  @test m[1].shape == [-1.0, -2.0, 10.0, -10.0, 2.0, 3.0, 5.0]
  consensus_shape = m[1].shape
  similar_motifs = find_motifs(intervals, consensus_shape, 1)
  @test similar_motifs.instances[1] ==  [-1.0, -2.0, 10.0, -10.0, 2.0, 3.0, 5.0]
end

