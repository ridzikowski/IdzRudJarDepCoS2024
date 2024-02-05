struct Instance
    n::Int64 # number of orders
    m::Int64 # number of destination
    readyTimes::Vector{Int64}
    weights::Vector{Int64}
    destination::Vector{Int64}
    travelTimes::Matrix{Float64}
    capacity::Int64
    timeLimit::Float64
    parkingTime::Float64
    serviceTime::Float64  
end

@enum Move SWAP REVERSE BEST_INSERTION
@enum Selection ROULETTE CUP
@enum Crossover RANDOM ORDER GREEDY

function Instance(n::Int64, m::Int64, capacity::Int64, timeLimit::Float64; rng::AbstractRNG = Random.MersenneTwister(1), parkingTime::Float64 = 2.0, serviceTime::Float64 = 1.0, maxWeights::Int64 = 5, maxReadyTime::Int64 = 20, speed::Float64 = 1.0)::Instance
    cords = []
    for _ in range(1, m + 1)
        push!(cords, (rand(rng) * 100,rand(rng) *100 ))
    end
    travelTimes = zeros(Float64,m + 1,m +1)
    for i in range(1,m + 1)
        for j in range(1,m + 1)
            travelTimes[i,j] = sqrt((cords[j][1]-cords[i][1])^2 + (cords[j][2]-cords[i][2])^2) * speed
        end
    end
    w = rand(rng, 1:maxWeights, n)
    d = rand(rng, 1:m, n)
    r = rand(rng, 1:maxReadyTime, n)
    Instance(n, m, r, w, d, travelTimes, capacity, timeLimit, parkingTime, serviceTime)
end

function Instance(fullTraveles::Matrix{Float64}, n::Int64, m::Int64, capacity::Int64, timeLimit::Float64; rng::AbstractRNG = Random.MersenneTwister(1), parkingTime::Float64 = 2.0, serviceTime::Float64 = 1.0, maxWeights::Int64 = 5, maxReadyTime::Int64 = 20, speed::Float64 = 1.0)::Instance
    @assert m <= 344
    hub = sample(rng, 1:3)
    location = sample(rng, 4:344, m, replace = false)
    allLocation = Int64[hub]
    append!(allLocation, location)
    travelTimes = zeros(Float64,m + 1,m +1)
    for i in 1:m+1
        for j in 1:m+1
            travelTimes[i,j] = fullTraveles[allLocation[i],allLocation[j]]
        end
    end
    w = rand(rng, 1:maxWeights, n)
    d = rand(rng, 1:m, n)
    r = rand(rng, 1:maxReadyTime, n)
    Instance(n, m, r, w, d, travelTimes, capacity, timeLimit, parkingTime, serviceTime)
end

function evaluate(π::Vector{Int64}, instance::Instance)::Float64
    totalTime = 0.0
    currentCapacity = 0
    currentTime = 0.0
    maxR = 0
    last =  0
    for idx in π[1:end]
        if idx == 0
            currentTime += maxR + instance.travelTimes[last + 1, 1] 
            if (currentCapacity > instance.capacity || currentTime > instance.timeLimit) 
                return Inf
            end

            totalTime += currentTime 

            currentCapacity = 0     
            currentTime = 0        
            maxR = 0 
            last = 0
        else
            if instance.readyTimes[idx] > maxR 
                maxR = instance.readyTimes[idx] 
            end
            destination = instance.destination[idx] 
            currentTime += instance.travelTimes[last + 1, destination + 1] + instance.serviceTime  
            currentCapacity += instance.weights[idx]  
            if last != destination 
                currentTime += instance.parkingTime 
            end
            last = destination   
        end
    end
    totalTime
end

function getRandom(instance::Instance; rng::AbstractRNG = Random.MersenneTwister(1))::Vector{Int64}
    
    πₓ = Int64[0]
    orders = collect(1:instance.n)
    shuffle!(rng, orders)
    currentCapacity = 0
    currentTime = 0
    last = 0
    maxR = 0
    for order in orders
        w = instance.weights[order]
        destination = instance.destination[order]
        travel = instance.travelTimes[last + 1, destination + 1] 
        if maxR < instance.readyTimes[order]
            maxR = instance.readyTimes[order]
        end
        combeback = instance.travelTimes[destination + 1, 1]
        approximate = maxR + instance.travelTimes[last + 1, destination + 1] + combeback + instance.serviceTime
        if destination != last
            approximate += instance.parkingTime 
        end
        if ((currentCapacity + w <= instance.capacity) && (currentTime + approximate  <= instance.timeLimit))
            currentCapacity += w
            currentTime += instance.travelTimes[last + 1, destination + 1] + instance.serviceTime
            if destination != last
                currentTime += instance.parkingTime 
            end
            append!(πₓ, order)
        else
            append!(πₓ, 0)
            append!(πₓ, order)
            currentCapacity = w
            currentTime = instance.travelTimes[1, destination + 1] + instance.serviceTime + instance.parkingTime
            maxR = instance.readyTimes[order]
        end
        last = destination
    end
    append!(πₓ, 0)
end