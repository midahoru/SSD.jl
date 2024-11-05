function solve_lazy_cuts(data, params, status, ρ_h, M)
    I = 1:data.I
    J = 1:data.J
    λ = data.a
    C = data.C
    F = data.F
    Q = data.Q
    cv = data.cv
    D = data.D    
    K = 1:data.k
    T = 1:data.t
    Dt = [D/sum(λ[i, t] for i in I) for t in T]

    # ρ_h = ini_ρ_h(data)  
    # M = calc_big_M(data, ρ_h)
    H = 1:size(ρ_h[:,:,:],3)
    
    maxtime = max(1, params.max_time - elapsed(status))
    m = Model(optimizer_with_attributes(Gurobi.Optimizer,
                                    "OutputFlag" => 0,
                                    "Threads" => 1,
                                    "MIPGap" => 1e-5,
                                    "MIPFocus" => 2,
                                    "TimeLimit" => maxtime + 1,
                                    "LazyConstraints" => 1,
                                    "Presolve" => 0,
                                    )
                                    )

    # @variable(m, x[I,J,T], Bin)
    @variable(m, 0 <= x[I,J,T] <= 1)
    @variable(m, y[J,K], Bin)
    # @variable(m, 0 <= y[J,K] <= 1)
    
    @variable(m, 0 <= z[J,K,T] <= 1)
    @variable(m, 0 <= ρ[J,T])

    @variable(m, 0 <= w[J,K,T])
    @variable(m, 0 <= R[J,T])

    @expression(m, of_term1, sum(F[j,k]*y[j,k] for j in J for k in K))
    @expression(m, of_term2, sum(C[i,j,t]*x[i,j,t] for i in I for j in J for t in T))
    @expression(m, of_term3, 0.5*sum(Dt[t]*(R[j,t] + ρ[j,t] + sum(cv^2*(w[j,k,t]-z[j,k,t]) for k in K)) for j in J for t in T))

    @objective(m, Min, of_term1 + of_term2 + of_term3)
    
    # Capacity cannot be exceeded and steady state has to be conserved
    @constraint(m, [j in J, t in T], sum(λ[i,t]*x[i,j,t] for i in I) - sum(Q[j,k]*y[j,k] for k in K) <= 0, set_string_name = false)
    
    # All customer zones need to be assigned to exactly one facility
    @constraint(m, [i in I, t in T], sum(x[i,j,t] for j in J) == 1, set_string_name = false)
    
    # At most one capacity level can be selected per facility
    @constraint(m, [j in J], sum(y[j,k] for k in K) <= 1, set_string_name = false)
    
    # 13 - 15 - 16 - 18
    @constraint(m, [j in J, t in T], sum(λ[i,t]*x[i,j,t] for i in I) - sum(Q[j,k]*z[j,k,t] for k in K) == 0, set_string_name = false)
    @constraint(m, [j in J, t in T], sum(z[j,k,t] for k in K) - ρ[j,t] == 0, set_string_name = false)
    @constraint(m, [j in J, t in T], sum(w[j,k,t] for k in K)-R[j,t] == 0, set_string_name = false)
    # 14 - 17
    @constraint(m, [j in J, t in T, k in K], z[j,k,t] - y[j,k] <= 0, set_string_name = false)
    
    # Upper bound ρ
    @constraint(m, [j in J, t in T], -ρ[j, t] >= -(1-10e-7), set_string_name = false)   

    @constraint(m, [j in J, t in T, k in K], w[j,k,t] - M[j,t]*y[j,k] <= 0, set_string_name = false)
    @constraint(m, [j in J, t in T, h in H], (1-ρ_h[j,t,h])^2 * R[j,t] - ρ[j,t] >= -ρ_h[j,t,h]^2, set_string_name = false)
    
    update_M = false
    function lazycb(cb)

        if callback_node_status(cb, m) != MOI.CALLBACK_NODE_STATUS_INTEGER
            return
        end

        ρvals = callback_value.(cb, ρ)
        Rvals = callback_value.(cb, R)
        
        update_ρ = false
        for j in J, t in T
            # Add if violation to the left of the right most tangent and if not added before
            if Rvals[j,t] < ρvals[j,t]/(1-ρvals[j,t]) # && !(ρvals[j,t] in ρ_h[j,t,:])
                # println("adding cut --------------------")
                update_ρ = true
                if ρvals[j,t] < maximum(ρ_h[j,t,:])
                    con = @build_constraint((1-ρvals[j,t])^2 * R[j,t] - ρ[j,t] >= -ρvals[j,t]^2)
                    MOI.submit(m, MOI.LazyConstraint(cb), con)
                else
                    update_M = true
                end
                
            end
        end
        status.nIter += 1

        # Update ρ_h
        if update_ρ
            ρ_h = cat(ρ_h, ρvals, dims=3)
        end
        
    end
    # write_to_file(m, "lazy_cuts_lp_relax_$(status.nIter).lp")
    MOI.set(m, MOI.LazyConstraintCallback(), lazycb)
    
    optimize!(m)
    end_stat = termination_status(m)
    println("$end_stat")
    if end_stat == MOI.INFEASIBLE
        status.endStatus = :infeasible
        return 10e8, 10e8, 10e8, 10e8, [], []
    elseif end_stat == MOI.OPTIMAL
        status.endStatus = :optimal
    elseif end_stat == MOI.TIME_LIMIT
        status.endStatus = :tlim
    elseif end_stat == MOI.ITERATION_LIMIT
        status.endStatus = :itlim
    elseif end_stat == MOI.SOLUTION_LIMIT
        status.endStatus = :sollim        
    end
        # xval = value.(x)
        # yval = value.(y)
        # zval = value.(z)
        # ρval = value.(ρ)
        # wval = value.(w)
        # Rval = value.(R)
        # optval = objective_value(m)
        # return xval, yval, zval, ρval, wval, Rval, optval
    tests_feas = is_sol_feas(data, value.(y).data, value.(x).data)

    return objective_value(m), value(of_term1), value(of_term2), value(of_term3), value.(y), value.(x), tests_feas, ρ_h, update_M
end

function model_lazy_cuts(data, params, status, n_outer_cuts)

    ρ_h = ini_ρ_h(data, n_outer_cuts)  
    M = calc_big_M(data, ρ_h)
    update_M = true
    of= loc_cost= al_cost= con_cost= y= x= tests_feas = 0
    while update_M
        of, loc_cost, al_cost, con_cost, y, x, tests_feas, ρ_h, update_M = solve_lazy_cuts(data, params, status, ρ_h, M)
        M = calc_big_M(data, ρ_h)
        println("we will update M = $update_M")
    end
    return of, loc_cost, al_cost, con_cost, y, x, tests_feas, [n_outer_cuts, size(ρ_h)[3]]
end