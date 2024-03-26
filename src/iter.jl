function cutting_plane(data, params, status, ρ_h)
    # Initialize the bounds
    lb, ub = 0, sum(data.F)+sum(data.C)+data.D*sum(data.a)

    xq, yq = [], []
    q = 0

    while (ub-lb)/ub >= params.ϵ
        xq, yq, ρ_q, Rq, wq, lb = model_cuts(data, params, status, ρ_h)
        newub = calc_ub(xq, yq, data)
        ub = min(ub, newub)
        ρ_h = cat(ρ_h, ρ_q, dims=3)
        q+=1
        println("Iter $q: LB= $lb ; UB = $ub")
    end
    status.nIter = q

    return lb, ub, xq, yq
end
