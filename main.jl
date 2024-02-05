using StatsBase
using Random
using DataFrames
using CSV
using Dates

include("src/Instance.jl")
include("src/greedy.jl")
include("src/ga.jl")

function main()
    travels = CSV.read("data/travelMatrixosrm.csv", DataFrame) |> Matrix{Float64}
    n = 1000
    m = 100
    seed = 1
    maxW = 10
    maxR = 240 
    C = 1000 
    T = 480.0 
    S = 1.0 
    P = 2.0
    rng = MersenneTwister(seed)
    inst = Instance(travels, n, m, C, T; rng = rng, parkingTime = P, serviceTime = S, maxWeights = maxW, maxReadyTime = maxR)
    greedy_value = greedy2(inst; type = both)[2]
    time = 3
    genetic_value = geneticAlgorithm(inst; timeLimit = time, rng = rng, crossover_type = GREEDY, mutation_type = SWAP, selection_type = CUP, populationSize = 40)[end][2]
    greedy_value/genetic_value
end

main()
