################################################## Instances

function instance_gen(nI, nJ, coords_bounds, a_bounds, r_bounds, cv, D, k, t, FLR, FCR, params)

    Icoords = round.(rand(params.rng, Uniform(coords_bounds[1], coords_bounds[2]), 2, nI),
        digits=params.round_digits)
    Jcoords = round.(rand(params.rng2, Uniform(coords_bounds[1], coords_bounds[2]), 2, nJ),
        digits=params.round_digits)
    dist = round.(pairwise(euclidean,Icoords, Jcoords), digits=params.round_digits)
    #= λ = round.(rand(params.rng, Uniform(λ_bounds[1], λ_bounds[2]), nI, t),
        digits=params.round_digits) =#
    a  = round.(rand(params.rng, Poisson(rand(params.rng, Uniform(a_bounds[1], a_bounds[2]))), nI, t),
        digits=0)
    r = round.(rand(params.rng, Uniform(r_bounds[1], r_bounds[2]), nI, nJ, t),
    digits=params.round_digits)
    C = round.(r .* dist, digits=params.round_digits)   

    write_file(nI, nJ, coords_bounds, cv, D, k, t, FLR, FCR, a, a_bounds, r_bounds, Icoords, Jcoords, C)
end

function write_file(nI, nJ, coords_bounds, cv, D, k, t, FLR, FCR, a, a_bounds, r_bounds,
     Icoords, Jcoords, C)
    filename = "instances/I_$nI J_$nJ $coords_bounds cv_$cv D_$D k_$k t_$t FLR_$FLR FCR_$FCR a_$a_bounds r_$r_bounds.txt"
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
        write(f, "abounds $a_bounds\n\n")        
        write(f, "a \n")
        for i in 1:size(Icoords,2)
            lams = a[i,:]
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
    open("$filename","r") do f
        A_ = false
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
                if A_
                    row = map(t -> parse(Float64, t), sline)
                    data.a[ini_count,:] = row
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
                    data.a = Array{Float64}(undef, data.I, data.t)
                    data.Icoords = Matrix{Float64}(undef, 2, data.I)
                    data.Jcoords = Matrix{Float64}(undef, 2, data.J)
                    data.C = Array{Float64}(undef, data.I, data.J, data.t)
                elseif sline[1] == "FLR"
                    data.FLR = parse(Float64, sline[2])
                elseif sline[1] == "FCR"
                    data.FCR = parse(Float64, sline[2])
                elseif sline[1] == "a"
                    A_ = true
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
                A_ = false
                i_coords = false
                j_coords = false
                c_ = false            
            end
        end
    end
    println("\ninstance file $filename parsed successfully\n")
end