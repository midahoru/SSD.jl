################################################## Models

function gen_caps(data, params, cap_levels)
    # Gets the max demand per time period
    max_dem_t = maximum(sum.(eachcol(data.a)))
    Q_j3 = zeros(Float64, data.J)
    for j in 1:data.J
        Q_j3[j] = 1.25*max_dem_t / (data.J * data.FLR)
    end
    return round.(Q_j3 .* cap_levels', digits=params.round_digits)
end

function gen_costs(data, params, cost_levels)
    f = (x,y) -> euclidean((x,y), (maximum(data.coords_bounds), maximum(data.coords_bounds))./2)
    F_j3 = zeros(Float64, data.J)
    for j in 1:data.J
        F_j3[j] = data.FCR*f(data.Jcoords[1,j], data.Jcoords[2,j])
    end
    return round.(F_j3 .* cost_levels', digits=params.round_digits)
end



################################################## Linear model with cuts
function ini_ρ_h(data)
    I = 1:data.I
    J = 1:data.J 
    T = 1:data.t 
    K = 1:data.k

    ini_ρ_h = collect(range(0.1,step=0.4,0.9))
    H = 1:length(ini_ρ_h)
    ρ_h = Array{Float64}(undef,data.J,data.t,length(ini_ρ_h))

    for j in J, t in T, h in H
        ρ_h[j, t, h]= ini_ρ_h[h]
    end

    return ρ_h
end

function calc_ub(ub, x, y, data)
    I = 1:data.I
    J = 1:data.J
    λ = data.a
    C = data.C
    F = data.F
    Q = data.Q
    cv = data.cv
    D = data.D    
    K = 1:data.k
    T = 1:data.t

    term_1 = sum(F[j,k]*y[j,k] for j in J for k in K) 
    term_2 = sum(C[i,j,t]*x[i,j,t] for i in I for j in J for t in T)
    term_3 = 0
    for t in T
        term_3_t = 0
        for j in J
            den_1 = 2*sum(Q[j,k]*y[j,k] for k in K)*(sum(Q[j,k]*y[j,k] for k in K)-sum(λ[i,t]*x[i,j,t] for i in I))
            den_2 = sum(Q[j,k]*y[j,k] for k in K)
            if den_1 != 0 && den_2 != 0
                sum_temp = ((1+cv^2*sum(y[j,k] for k in K))*sum(λ[i,t]*x[i,j,t] for i in I)^2/(den_1))+
                +(sum(λ[i,t]*x[i,j,t] for i in I)/den_2)
                term_3_t += sum_temp
            end        
        end
        term_3 += (D/sum(λ[i,t] for i in I))*term_3_t
    end
    new_ub = term_1 + term_2 + term_3
    return minimum([ub, new_ub])
end

function calc_new_ρ(xq, yq, data)
    I = 1:data.I
    J = 1:data.J 
    λ = data.a
    T = 1:data.t 
    K = 1:data.k
    Q = data.Q    
    T = 1:data.t

    ρ_new = Array{Float64}(undef,data.J,data.t)
    for j in J
        for t in T
            num = sum(λ[i,t].*xq[i,j,t] for i in I)
            den = sum(Q[j,k]*yq[j,k] for k in K)
            if den != 0
                ρ_new[j,t]=num/den
            else
                ρ_new[j,t]=0
            end
        end
    end
    return ρ_new
end






# Calculates the maximum dissimilarity between two figures
function calc_max_diss(c1, d1, c2, d2, params)
    return round(maximum(pairwise(euclidean, extreme_p_generator(c1, d1), extreme_p_generator(c2, d2))), digits = params.round_digits)
end


################################################## miscellaneous

function set_maximum_time(params, s)
    params.max_time = s
end

function elapsed(status)
    b = Dates.now()
    ms = Dates.value(b - status.initTime)
    s = ms / 1000
end

function total_elapsed(status)
    ms = Dates.value(status.endTime - status.initTime)
    s = ms / 1000
end

function optimal(status)
    return status.endStatus == :optimal
end

function get_status(status)
    return status.endStatus
end

function get_number_groups(status)
    return status.endGroupsNb
end

function set_data(data::Data, name, length, heigth)
    data.name = name
    data.len = length
    data.hei = hei
    return data
end

function set_params(params::Parameters, max_time, rounding_digits, seed)
    params.max_time = max_time
    params.round_digits = rounding_digits
    params.rng = MersenneTwister(seed)
    return params
end

##################################################


