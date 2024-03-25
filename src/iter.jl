function cutting_plane(data, params, status, ρ_h)
    # Initialize the bounds
    lb, ub = 0, sum(data.F)+sum(data.C)+data.D*sum(data.a)

    xq, yq, lb = [], [], 0
    q = 0
    while (ub-lb)/ub >= params.ϵ
        println("Current LB= $lb ; UB = $ub")
        xq, yq, ρ_q, Rq, lb = model_cuts(data, params, status, ρ_h)
        println("Optimal LB= $lb")
        
        ub = calc_ub(ub, xq, yq, data)

        other_ub = lb + sum(data.D/sum(data.a[i, t] for i in 1:data.I)*((ρ_q[j, t]/(1-ρ_q[j,t]))-Rq[j,t]) for j in 1:data.J for t in 1:data.t)
        println("other_ub - ub = $(other_ub - ub)")
        diff = [(ρ_q[j, t] / (1 - ρ_q[j, t])) - Rq[j, t] for j in 1:data.J for t in 1:data.t]
        println("max diff = $(maximum(diff))")
        println("min diff = $(minimum(diff))")
        # @assert(sum(abs.(ρ_q - ρ_new)) < 1e-5, "should be the same")
        # println("ρ_new = $ρ_new")
        # println(maximum(ρ_new))
        # println(minimum(ρ_new))
        println("rho = $ρ_q")
        cat(ρ_h, ρ_q, dims=3)
        q+=1
        println("Iter $q: LB= $lb ; UB = $ub")
    end
    status.nIter = q

    return xq, yq, lb   
end
