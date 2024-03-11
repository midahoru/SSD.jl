################################################## Instances

function instance_gen(nI, nJ, coords_bounds, λ_bounds, r_bounds, cv, D, k, t, params)

    Icoords = round.(rand(params.rng, Uniform(coords_bounds[1], coords_bounds[2]), 2, nI),
        digits=params.round_digits)
    Jcoords = round.(rand(params.rng2, Uniform(coords_bounds[1], coords_bounds[2]), 2, nJ),
        digits=params.round_digits)
    dist = round.(pairwise(euclidean,Icoords, Jcoords), digits=params.round_digits)
    λ = round(rand(params.rng, Uniform(λ_bounds[1], λ_bounds[2])),
        digits=params.round_digits)
    r = round.(rand(params.rng, Uniform(r_bounds[1], r_bounds[2]), nI, nJ, t),
    digits=params.round_digits)
    C = round.(r .* dist, digits=params.round_digits)

    write_file(nI, nJ, coords_bounds, cv, D, k, t, λ, λ_bounds, r_bounds, Icoords, Jcoords, C)
end

function write_file(nI, nJ, coords_bounds, cv, D, k, t, λ, λ_bounds, r_bounds, Icoords, Jcoords, C)
    fname = "instances/I_$nI J_$nJ $coords_bounds cv_$cv D_$D k_$k t_$t lam_$λ_bounds r_$r_bounds.txt"
    open(fname, "w") do f
        write(f, "I $nI\n")
        write(f, "J $nJ\n")
        write(f, "bounds $coords_bounds\n")
        write(f, "cv $cv\n")
        write(f, "D $D\n")
        write(f, "k $k\n")
        write(f, "t $t\n")
        write(f, "λ $λ\n")
        write(f, "λ bounds $λ_bounds\n\n")
        write(f, "Icoords\n")
        for i in 1:size(Icoords,2)
            coords = Icoords[:,i]
            x = coords[1]
            y = coords[2]
            write(f, "$x $y\n")
        end
        write(f, "\nJcoords\n")
        for j in 1:size(Jcoords,2)
            coords = Jcoords[:,j]
            x = coords[1]
            y = coords[2]
            write(f, "$x $y\n")
        end      
        write(f, "\nC\n")
        for t in 1:size(C,3)
            Ct = C[:,:,t]
            for j in 1:size(Ct, 1)
                coords = Ct[j,:]
                for x in coords
                    write(f, "$x ")
                end
                write(f, "\n")
            end
            write(f, "\n")
        end 
        write(f, "end")
    end
end

function read_file(filename)    
    start_coords = false
    open("instances/"+filename) do f
        lines = readlines(f)
        for line in lines
            sline = split(line)
            if sline[1] == "EOF"
                break
            elseif start_coords
                p = parse(Int64, sline[1])
                x = parse(Float64, sline[2])
                y = parse(Float64, sline[3])
                data.D[1, p] = x
                data.D[2, p] = y
            elseif sline[1] == "NAME"
                data.name = sline[3]
            elseif sline[1] == "DIMENSION"
                data.nnodes = parse(Int64, sline[3])
                data.D = zeros(Float64, 2, data.nnodes)
            elseif sline[1] == "EDGE_WEIGHT_TYPE"
                if sline[3] == "EUC_2D"
                    params.wtype = :round
                elseif sline[3] == "CEIL_2D"
                    params.wtype = :ceil
                elseif sline[3] == "GEOM"
                    params.wtype = :geom
                end
            elseif sline[1] == "NODE_COORD_SECTION"
                start_coords = true
            end
        end
    end
    println("instance file $filename parsed successfully")

end

    
################################################## initialClustering




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


