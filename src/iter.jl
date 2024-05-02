function cuts_priori(data, params, status)    

    # Initialize the bounds
    lb, ub = 0, sum(data.F)+sum(data.C)+data.D*sum(data.a)

    xq, yq, zq, ρq, wq, Rq = [], [], [], [], [], []
    ρ_h = ini_ρ_h(data)
    q = 0

    while (ub-lb)/ub >= params.ϵ
        xq, yq, zq, ρq, wq, Rq, new_lb = model_cuts_priori(data, params, status, ρ_h)
        if (new_lb-ub)/ub > params.ϵ 
            break
        end
        if status.endStatus == :infeasible || status.endStatus == :none
            break
        end
        lb = new_lb
        newub = calc_ub(lb, ρq, Rq, wq, data)#calc_ub(xq, yq, data) # 
        ub = min(ub, newub)
        ρ_h = cat(ρ_h, ρq, dims=3)
        q+=1
        println("Iter $q: LB= $lb ; UB = $ub")
    end
    status.nIter = q

    # return xq, yq, zq, ρq, wq, Rq, lb, ub
    return lb, ub
end
