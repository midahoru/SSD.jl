#= function model_cuts(data, params, cost_levels, cap_levels)
    #cost_levels = [0.60, 0.85, 1, 1.15, 1.35]
    #cap_levels = [0.5, 0.75, 1, 1.25, 1.5]
    I = 1:data.I
    J = 1:data.J
    λ = data.λ
    F = data.F
    Q = data.Q
    C = data.C
    cv = data.cv
    D = data.D    
    K = 1:data.k
    T = 1:data.t
    FLR = data.FLR
    FCR = data.FCR
    M = data.M

    F = gen_costs(data, params, cost_levels)
    Q = gen_caps(data, params, cap_levels)   
    
    maxtime = max(1, params.max_time - elapsed())
    m = Model(optimizer_with_attributes(Gurobi.Optimizer,
                                    "OutputFlag" => 0,
                                    "Threads" => 1,
                                    "MIPFocus" => 1,
                                    "TimeLimit" => maxtime + 1,
                                    "SolutionLimit" => 1
                                    )
                )
    @variable(m, x[I,J,T], Bin)
    @variable(m, y[J,K], Bin)
    
    @variable(m, 0 <= z[J,K,T] <= 1 )
    @variable(m, 0 <= ρ[J,T] <= 1)

    @variable(m, 0 <= w[J,K,T])
    @variable(m, 0 <= R[J,T])
    
    @objective(m, Min, sum(F[j,k]*y[j,k] for j in J for k in K) +
    sum(C[i,j,t]*x[i,j,t] for i in I for j in J for t in T) +
    0.5*sum(D*(R[j,t] + ρ[j,t] + sum(cv^2*(w[j,k,t]-z[j,k,t]) for k in K)) for j in J for t in T))

    # Capacity cannot be exceeded and steady state is always true
    for j in J, t in T
        @constraint(m, sum(a[i,t]*x[i,j,t] for i in I) <= sum(Q[j,k]*y[j,k] for k in K))
    end

    # All customer zones need to be assigned to exactly one facility
    for i in I, t in T
        @constraint(m, sum(x[i,j,t] for j in J) == 1)
    end

    # At most one capacity level can be selected per facility
    for j in J
        @constraint(m, sum(y[j,k] for k in K) <= 1)
    end

    # 14 - 16 - 17 - 19
    for j in J, t in T
        @constraint(m, sum(a[i,t]*x[i,j,t] for i in I) - sum(Q[j,k]*z[j,k,t] for k in K) == 0)
        @constraint(m, sum(z[j,k,t] for k in K) - ρ[j,t] == 0)
        @constraint(m, R[j,t] = ρ[j,t]/(1-ρ[j,t]))
        @constraint(m, sum(w[j,k,t] for k in K)-R[j,t] == 0)
    end
    # 15 - 18
    for j in J, t in T, k in K
        @constraint(m, z[j,k,t] - y[j,k] <= 0)
        @constraint(m, w[j,k,t] - M*y[j,k])
    end 

    #=function lazycb(cb)
        xvals = callback_value.(cb, x)
        zvals = callback_value.(cb, z)
        
        # if callback_node_status(cb, m) == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        #     return
        # end
        if callback_node_status(cb, m) == MOI.CALLBACK_NODE_STATUS_INTEGER
            let
                sol = maximum_weighted_clique(nnodes, unfeas, xvals, true) 
                if sum(xvals[i] for i in sol) > 1+1e-7
                    println("adding clique inequality of type I --------------------")
                    con = @build_constraint(sum(x[u] for u in sol) <= 1)
                    MOI.submit(m, MOI.LazyConstraint(cb), con)
                    return
                end
            end
        end
        
        for k in K
            weights = [xvals[i] + zvals[k] - 1 for i in I]
            # sol = maximum_weighted_clique(nnodes, E_D[k:end], weights, false) # exact
            if !isempty(E_D[1:k-1])
                sol = maximum_weighted_clique(nnodes, E_D[1:k-1], weights, false) # heuristic                        
                nsol = length(sol)
                if nsol > 0 && sum(xvals[i] for i in sol) + (nsol - 1) * zvals[k] > nsol + 1e-1
                    println("adding clique inequality of type II --------------------")
                    con = @build_constraint(sum(x[i] for i in sol) + (nsol - 1) * z[k] <= nsol)
                    MOI.submit(m, MOI.LazyConstraint(cb), con)
                end
            end
        end
    end

    MOI.set(m, MOI.LazyConstraintCallback(), lazycb) =#

    optimize!(m)
    status = termination_status(m)
    if status == MOI.OPTIMAL || status == MOI.SOLUTION_LIMIT
        xval = value.(x)
        yval = value.(y)
        optval = objective_value(m)
        return xval, yval, optval
    else return [], [], 0
    end
end =#