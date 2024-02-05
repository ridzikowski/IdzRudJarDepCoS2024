mutable struct Individual
    chromosome::Vector{Int64}
    phenotype::Float64
    fittnes::Float64
end 

function populationInit(instance::Instance, rng::AbstractRNG, size_of_poplution::Int64)::Tuple{Vector{Individual}, Vector{Int64}, Float64}
    population = Individual[]
    fᵦ = Inf
    πᵦ = []
    for x in 1:size_of_poplution
        
        π = getRandom(instance; rng = rng)
        fit = evaluate(π, instance)
        individual = Individual(π, fit, 1/fit)
        push!(population, individual)
        if fit < fᵦ
            fᵦ = fit
            πᵦ = π
        end
    end
    population, πᵦ, fᵦ
end

function selection(type::Selection, population::Vector{Individual}, rng::AbstractRNG, fᵦ::Float64; tournamentSize::Int64 = 4)::Vector{Individual}
    parents = Individual[]
    number = length(population)
    for _ in range(1,number)
        if type == ROULETTE
            idx = sample(rng,1:number)
            while  rand(rng) > population[idx].phenotype/fᵦ
                idx = sample(rng,1:number)
            end
            push!(parents, population[idx])
        elseif type == CUP
            participant = sample(rng,1:number)
            for _ in range(2, tournamentSize)
                opponent = sample(rng,1:number)
                if population[opponent].phenotype < population[participant].phenotype
                    participant = opponent
                end
            end
            push!(parents, population[participant])
        end
    end
    parents
end

function crossover(instance::Instance, type::Crossover, population::Vector{Individual}, rng::AbstractRNG; probability::Float64=0.8)::Vector{Individual}
    children = Individual[]
    number = length(population)
    half = floor(number/2) |> Int64
    all_orders = Set(range(1,maximum(population[1].chromosome)))
    for idx in range(1,half)
        if rand(rng) > probability
            push!(children, population[idx])
            push!(children, population[idx + half])
            continue
        end
        cut1 = sample(rng,2:min(length(population[idx].chromosome)-1,length(population[idx + half].chromosome)-1))
        cut2 = sample(rng,2:min(length(population[idx].chromosome)-1,length(population[idx + half].chromosome)-1))
        if cut2 < cut1
            cut1, cut2 = cut2, cut1
        end
        child1 = copy(population[idx       ].chromosome[cut1:cut2])
        child2 = copy(population[idx + half].chromosome[cut1:cut2])

        for child in [child1, child2]
            if child[1] != 0
                insert!(child, 1, 0)
            end
            orders = setdiff(all_orders,Set(child))
            currentCapacity = 0
            currentTime = 0.0
            maxR = 0
            last = 0
            for id in child
                if id == 0
                    currentCapacity = 0
                    currentTime = 0.0
                    maxR = 0
                    last = 0
                else
                    destination = instance.destination[id]
                    currentCapacity += instance.weights[id]
                    currentTime += instance.travelTimes[last + 1, destination + 1] + instance.serviceTime
                    if destination != last
                        currentTime += instance.parkingTime 
                    end
                    if instance.readyTimes[id] > maxR
                        maxR = instance.readyTimes[id]
                    end
                    last = destination
                end

            end
            while !isempty(orders)
                order = 0
                if type == RANDOM
                    order = sample(rng, collect(orders))
                elseif type == ORDER
                    order = collect(orders)[1]
                elseif type == GREEDY
                    order = collect(orders)[1]
                    source = child[end]
                    if source != 0
                        source = instance.destination[source]
                    end
                    travel = instance.travelTimes[source + 1, instance.destination[order] + 1]
                    for id in  collect(orders)[2:end]
                        minTravel = instance.travelTimes[source + 1, instance.destination[id] + 1]
                        if minTravel < travel
                            order = id
                            travel = minTravel
                        end
                    end
                end
                
                delete!(orders, order)
                destination = instance.destination[order]
                if instance.readyTimes[order] > maxR
                        maxR = instance.readyTimes[order]
                end
                # print("last :", destination, " dest: ", destination)
                approximate = maxR + instance.travelTimes[last + 1, destination + 1] + instance.travelTimes[ destination + 1, 1] + instance.serviceTime
                if destination != last
                    approximate += instance.parkingTime 
                end
                if (currentCapacity + instance.weights[order] <= instance.capacity && currentTime + approximate  <= instance.timeLimit )
                    currentCapacity += instance.weights[order]
                    currentTime +=  instance.travelTimes[last + 1, destination + 1] + instance.serviceTime
                    if destination != last
                        currentTime += instance.parkingTime 
                    end
                    
                    append!(child, order)
                else
                    append!(child, 0)
                    append!(child, order)
                    currentCapacity = instance.weights[order]
                    currentTime = instance.travelTimes[1, destination + 1] + instance.serviceTime + instance.parkingTime 
                    maxR = instance.readyTimes[order]
                    last = destination
                end
                last = destination
            end
            append!(child, 0)
            fit = evaluate(child, instance)
            push!(children,Individual(child,fit,1/fit))
        end
    end
    children
end

function mutation(instance::Instance, type::Move, population::Vector{Individual}, rng::AbstractRNG; probability::Float64 = 0.05)::Vector{Individual}
    size = length(population)
    for idx in range(1,size)
        if rand(rng) < probability
            individual = copy(population[idx].chromosome)
            cut1 = sample(rng,2:length(individual)-1)
            cut2 = sample(rng,2:length(individual)-1)
            if cut2 < cut1
                cut1, cut2 = cut2, cut1
            end
            if type == REVERSE
                reverse!(individual, cut1, cut2)
            elseif type == SWAP
                individual[cut1], individual[cut2] = individual[cut2], individual[cut1]
            end
            fit = evaluate(individual, instance)
            if fit != Inf
                population[idx].chromosome = individual
                population[idx].phenotype = fit
                population[idx].fittnes = 1/fit
            end
        end
    end
    population
end

function update(population::Vector{Individual}, πᵦ::Vector{Int64}, fᵦ::Float64)::Tuple{Vector{Int64}, Float64}
    size = length(population)
    for id in range(1,size)
        if population[id].phenotype < fᵦ
            πᵦ = population[id].chromosome
            fᵦ = population[id].phenotype
        end
    end
    πᵦ, fᵦ
end

function geneticAlgorithm(instance::Instance; iterationLimit::Int64=-1, timeLimit::Int64=-1, rng::AbstractRNG = Random.MersenneTwister(1), populationSize::Int64 = 64, selection_type::Selection = ROULETTE, crossover_type::Crossover = RANDOM, mutation_type::Move = REVERSE)::Vector{Tuple{Vector{Int64}, Float64}}
    population, πᵦ, fᵦ = populationInit(instance, rng, populationSize)
    t = round(sqrt(populationSize)) |> Int64
    iteration = 0
    start = Dates.now()
    output = Tuple{Vector{Int64}, Float64}[]
    point1 = true
    point2 = true
    point3 = true
    quarter = floor(timeLimit/4) |> Int64
    half = floor(timeLimit/2) |> Int64
    three = floor(3*timeLimit/4) |> Int64
    while (iterationLimit > 0 && iteration < iterationLimit) || (timeLimit > 0 && Dates.now() - start <  Second(timeLimit))
        iteration +=1
        population = selection(selection_type, population, rng, fᵦ; tournamentSize = t)
        population = crossover(instance, crossover_type, population, rng)
        πᵦ, fᵦ = update(population, πᵦ, fᵦ)
        population = mutation(instance, mutation_type, population, rng)
        if  (timeLimit > 0 && Dates.now() - start >  Second(timeLimit))
            while length(output) < 3
                push!(output, (πᵦ, fᵦ))
            end
            push!(output, (πᵦ, fᵦ))
            return output
        end
        πᵦ, fᵦ = update(population, πᵦ, fᵦ)
        if (point1 &&  Dates.now() - start >  Second(quarter))
            push!(output, (πᵦ, fᵦ))
            point1 = false
        end
        if (point2 &&  Dates.now() - start >  Second(half))
            push!(output, (πᵦ, fᵦ))
            point2 = false
        end
        if (point3 &&  Dates.now() - start >  Second(three))
            push!(output, (πᵦ, fᵦ))
            point3 = false
        end
    end
    
    push!(output, (πᵦ, fᵦ))
    output
end