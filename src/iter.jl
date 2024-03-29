function cutting_plane(data, params, status, ρ_h)
    # Initialize the bounds
    lb, ub = 0, sum(data.F)+sum(data.C)+data.D*sum(data.a)

    xq, yq = [], []
    q = 0

    while (ub-lb)/ub >= params.ϵ
        xq, yq, ρq, Rq, wq, lb = model_cuts(data, params, status, ρ_h)

        ub = calc_ub(lb, ub, ρq, Rq, wq, data)

        ρ_h = cat(ρ_h, ρq, dims=3)
        q+=1
        println("Iter $q: LB= $lb ; UB = $ub")
    end
    status.nIter = q

    return lb, ub, xq, yq
end
