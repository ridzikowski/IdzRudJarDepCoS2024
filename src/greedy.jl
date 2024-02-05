@enum GreedySortType ready weight both none 

function greedy(instance::Instance; type::GreedySortType = none)::Tuple{Vector{Int64}, Float64}
    vehicles = instance.n 
    π = zeros(Int64, vehicles + 1) 
    orders = collect(1:instance.n)
    if type == ready
        r = deepcopy(instance.readyTimes)
        perm = sortperm(r)
        orders = reverse(orders[perm])
    elseif type == weight
        w = deepcopy(instance.weights)
        perm = sortperm(w)
        orders = reverse(orders[perm])
    elseif  type == both
        r = deepcopy(instance.readyTimes)
        w = deepcopy(instance.weights)
        s = r .* w
        perm = sortperm(s)
        orders = reverse(orders[perm])
    end
    b = 2
    insert!(π, 2, orders[1])
    for i in orders[2:end]
        insert!(π, 2, i)
        best = evaluate(π, instance) 
        b = 2
        for k in 2:length(π)-2
            π[k] , π[k+1] = π[k+1] , π[k]
            value = evaluate(π, instance)
            if value < best
                best = value
                b = k + 1 
            end
        end
        deleteat!(π, length(π)-1)
        insert!(π, b, i)
    end
    π, evaluate(π, instance)
end

function greedy2(instance::Instance; type::GreedySortType = none)::Tuple{Vector{Int64}, Float64}
    π = zeros(Int64, 2) 
    orders = collect(1:instance.n)
    if type == ready
        r = deepcopy(instance.readyTimes)
        perm = sortperm(r)
        orders = reverse(orders[perm])
    elseif type == weight
        w = deepcopy(instance.weights)
        perm = sortperm(w)
        orders = reverse(orders[perm])
    elseif  type == both
        r = deepcopy(instance.readyTimes)
        w = deepcopy(instance.weights)
        s = r .* w
        perm = sortperm(s)
        orders = reverse(orders[perm])
    end
    b = 2
    insert!(π, 2, orders[1])
    for i in orders[2:end]
        if π[end-1] != 0
            append!(π, 0)
        end
        insert!(π, 2, i)
        best = evaluate(π, instance) 
        b = 2
        for k in 2:length(π)-2
            π[k] , π[k+1] = π[k+1] , π[k]
            value = evaluate(π, instance)
            if value < best
                best = value
                b = k + 1 
            end
        end
        deleteat!(π, length(π)-1)
        insert!(π, b, i)
    end
    π, evaluate(π, instance)
end