################################################## Instances

function instance_gen(nI, nJ, coords_bounds, λ_bounds, r_bounds, cv, D, k, t, FLR, FCR, params)

    Icoords = round.(rand(params.rng, Uniform(coords_bounds[1], coords_bounds[2]), 2, nI),
        digits=params.round_digits)
    Jcoords = round.(rand(params.rng2, Uniform(coords_bounds[1], coords_bounds[2]), 2, nJ),
        digits=params.round_digits)
    dist = round.(pairwise(euclidean,Icoords, Jcoords), digits=params.round_digits)
    λ = round.(rand(params.rng, Uniform(λ_bounds[1], λ_bounds[2]), nI, t),
        digits=params.round_digits)
    r = round.(rand(params.rng, Uniform(r_bounds[1], r_bounds[2]), nI, nJ, t),
    digits=params.round_digits)
    C = round.(r .* dist, digits=params.round_digits)   

    write_file(nI, nJ, coords_bounds, cv, D, k, t, FLR, FCR, λ, λ_bounds, r_bounds, Icoords, Jcoords, C)
end

function write_file(nI, nJ, coords_bounds, cv, D, k, t, FLR, FCR, λ, λ_bounds, r_bounds,
     Icoords, Jcoords, C)
    filename = "instances/I_$nI J_$nJ $coords_bounds cv_$cv D_$D k_$k t_$t FLR_$FLR FCR_$FCR lam_$λ_bounds r_$r_bounds.txt"
    open(filename, "w") do f
        write(f, "I $nI\n")
        write(f, "J $nJ\n")
        write(f, "bounds $coords_bounds\n")
        write(f, "cv $cv\n")
        write(f, "D $D\n")
        write(f, "k $k\n")
        write(f, "t $t\n")
        write(f, "FLR $FLR\n")
        write(f, "FCR $FCR\n")
        write(f, "λbounds $λ_bounds\n\n")        
        write(f, "λ \n")
        for i in 1:size(Icoords,2)
            lams = λ[i,:]
            for l in lams
                write(f, "$l ")
            end
            write(f,"\n")
        end
        write(f,"\n")
        write(f, "Icoords\n")
        for i in 1:nI
            coords = Icoords[:,i]
            x = coords[1]
            y = coords[2]
            write(f, "$x $y\n")
        end
        write(f,"\n")
        write(f, "Jcoords\n")
        for j in 1:nJ
            coords = Jcoords[:,j]
            x = coords[1]
            y = coords[2]
            write(f, "$x $y\n")
        end 
        write(f,"\n")     
        write(f, "C\n")
        for t in 1:t
            Ct = C[:,:,t]
            for i in 1:nI
                coords = Ct[i,:]
                for j in coords
                    write(f, "$j ")
                end
                write(f, "\n")
            end
        end 
        write(f, "EOF")
    end
    println("\ninstance file $filename saved successfully\n")
end

function read_file(filename, data)    
    open("instances/$filename","r") do f
        λ_ = false
        i_coords = false
        j_coords = false
        c_ = false
        ini_count = 1      
        t_count = 1        
        
        while !eof(f)         
            line = readline(f)
            sline = split(line)
            if line == "EOF"
                break
            elseif line != ""
                if λ_
                    row = map(t -> parse(Float64, t), sline)
                    data.λ[ini_count,:] = row
                    ini_count += 1
                elseif i_coords
                    col = map(t -> parse(Float64, t), sline)
                    data.Icoords[:, ini_count] = col
                    ini_count += 1
                elseif j_coords
                    col = map(t -> parse(Float64, t), sline)
                    data.Jcoords[:, ini_count] = col
                    ini_count += 1
                elseif c_
                    row = map(t -> parse(Float64, t), sline)
                    data.C[ini_count, :, t_count] = row
                    if ini_count == data.I
                        ini_count = 0
                        t_count+=1
                    end                
                    ini_count+=1                
                elseif sline[1] == "I"
                    data.I = parse(Int64, sline[2])
                elseif sline[1] == "J"
                    data.J = parse(Int64, sline[2])                
                elseif sline[1] == "bounds"
                    lb = parse(Float64, replace(sline[2],"(" => "", ","=>""))
                    ub = parse(Float64, replace(sline[3],")" => ""))
                    data.coords_bounds = Tuple([lb, ub])
                elseif sline[1] == "cv"
                    data.cv = parse(Float64, sline[2])
                elseif sline[1] == "D"
                    data.D = parse(Float64, sline[2])
                elseif sline[1] == "k"
                    data.k = parse(Int64, sline[2])
                elseif sline[1] == "t"
                    data.t = parse(Int64, sline[2])
                    data.λ = Array{Float64}(undef, data.I, data.t)
                    data.Icoords = Matrix{Float64}(undef, 2, data.I)
                    data.Jcoords = Matrix{Float64}(undef, 2, data.J)
                    data.C = Array{Float64}(undef, data.I, data.J, data.t)
                elseif sline[1] == "FLR"
                    data.FLR = parse(Float64, sline[2])
                elseif sline[1] == "FCR"
                    data.FCR = parse(Float64, sline[2])
                elseif sline[1] == "λ"
                    λ_ = true
                elseif sline[1] == "Icoords"
                    i_coords = true
                    ini_count = 1
                elseif sline[1] == "Jcoords"
                    j_coords = true
                    ini_count = 1
                elseif sline[1] == "C"
                    c_ = true
                    ini_count = 1
                end
            elseif line == ""
                λ_ = false
                i_coords = false
                j_coords = false
                c_ = false            
            end
        end
    end
    println("\ninstance file $filename parsed successfully\n")
end

    
################################################## Models

function gen_caps(data, params, cap_levels)
    # Gets the max demand per time period
    max_dem_t = maximum(sum.(eachcol(data.λ)))
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

function calc_ub(ub, x, y, data)
    I = 1:data.I
    J = 1:data.J
    λ = data.λ
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
                sum_temp = ((1+cv^2*sum(y[j,k] for k in K)*sum(λ[i,t]*x[i,j,t] for i in I))/(den_1))+
                +(sum(λ[i,t]*x[i,j,t] for i in I)/den_2)
                #println(sum_temp)
                term_3_t += sum_temp
            end        
        end
        term_3 += D*term_3_t/sum(λ[i,t] for i in I)
    end
    new_ub = term_1 + term_2 + term_3
    println(term_1)
    println(term_2)
    println(term_3)
    println(ub)
    return minimum([ub, new_ub])
end

function calc_new_ρ(xq, yq, data)
    I = 1:data.I
    J = 1:data.J 
    T = 1:data.t 
    K = 1:data.k

    ρ_new = Array{Float64}(undef,data.J,data.t)
    for j in J
        for t in T
            num = sum(data.λ[i,t].*xq[i,j,t] for i in I)
            den = sum(data.Q[j,k]*yq[j,k] for k in K)
            if den > 0
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


