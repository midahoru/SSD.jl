function cuts_priori(data, params, status, ρ_h)
    # Initialize the bounds
    lb, ub = 0, sum(data.F)+sum(data.C)+data.D*sum(data.a)

    xq, yq, zq, ρq, wq, Rq = [], [], [], [], [], []
    q = 0

    while (ub-lb)/ub >= params.ϵ
        xq, yq, zq, ρq, wq, Rq, lb = model_cuts_priori(data, params, status, ρ_h)
        newub = calc_ub(xq, yq, data) # calc_ub(lb, ρq, Rq, wq, data)
        ub = min(ub, newub)
        ρ_h = cat(ρ_h, ρq, dims=3)
        q+=1
        println("Iter $q: LB= $lb ; UB = $ub")
    end
    status.nIter = q

    return xq, yq, zq, ρq, wq, Rq, lb, ub
end
