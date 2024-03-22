function cutting_plane(data, params, status, ρ_h)
    # Initialize the bounds
    lb, ub = 0, sum(data.F)+sum(data.C)+data.D*sum(data.a)

    xq, yq, lb = [], [], 0
    q = 0
    while (ub-lb)/lb >= params.ϵ
        println("Current LB= $lb ; UB = $ub")
        xq, yq, lb = model_cuts(data, params, status, ρ_h)
        println("Optimal LB= $lb")
        
        ub = calc_ub(ub, xq, yq, data)

        ρ_new = calc_new_ρ(xq, yq, data)
        println("ρ_new = $ρ_new")
        println(maximum(ρ_new))
        println(minimum(ρ_new))

        cat(ρ_h, ρ_new, dims=3)
        q+=1
        println("Iter $q: LB= $lb ; UB = $ub")
    end
    status.nIter = q

    return xq, yq, lb   
end
