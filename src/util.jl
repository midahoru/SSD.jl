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
    J = 1:data.J 
    T = 1:data.t

    ini_ρ_h = collect(range(0.1,step=0.1,0.9))
    H = 1:length(ini_ρ_h)
    ρ_h = Array{Float64}(undef,data.J,data.t,length(ini_ρ_h))

    for j in J, t in T, h in H
        ρ_h[j, t, h]= ini_ρ_h[h]
    end

    return ρ_h
end


function calc_ub(x, y, data)
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

    ρ = calc_new_ρ(x, y, data)
    R = calc_new_R(ρ, data)
    w, z = calc_new_w_z(x, y, R, ρ, data)

    Fterm = dot(F, y)
    Cterm = dot(C, x)

    Dt = [data.D/sum(λ[i, t] for i in I) for t in T]
    
    last_term = 0.5 * sum(Dt[t] * (R[j, t] + ρ[j, t] + sum(data.cv^2 * (w[j, k, t] - z[j, k, t]) for k in K)) for j in J for t in T)

    ub = Fterm + Cterm + last_term
    return ub
end

function calc_ub(lb, ub, ρq, Rq, wq, data)
    Dt = [data.D/sum(data.a[i, t] for i in 1:data.I) for t in 1:data.t]

    new_ub = lb + sum(Dt[t]*(((ρq[j, t]/(1-ρq[j,t]))-Rq[j,t]) + data.cv^2*((ρq[j, t]/(1-ρq[j,t]))-
    sum(wq[j,k,t] for k in 1:data.k))) for j in 1:data.J for t in 1:data.t)
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

    ρ_new = zeros(data.J, data.t)
    for j in J, t in T
        num = sum(λ[i,t].*xq[i,j,t] for i in I)
        den = sum(Q[j,k]*yq[j,k] for k in K)
        if den != 0
            ρ_new[j,t]=num/den
        end
    end
    return ρ_new
end

function calc_new_R(ρ, data)
    R = zeros(data.J, data.t)
    for j in 1:data.J, t in 1:data.t
        R[j, t] = ρ[j, t] / (1 - ρ[j, t])
    end
    return R
end

function calc_new_w_z(x, y, R, ρ, data)
    z = zeros(data.J, data.k, data.t)
    w = zeros(data.J, data.k, data.t)
    for j in 1:data.J, k in 1:data.k, t in 1:data.t
        z[j, k, t] = ρ[j, t] * y[j, k]
        w[j, k, t] = R[j, t] * y[j, k]
    end
    return w, z
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


