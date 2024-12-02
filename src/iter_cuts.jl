function cuts_priori(data, params, status, n_outer_cuts)   

    xq, yq, zq, ρq, wq, Rq = [], [], [], [], [], []
    ρ_h = ini_ρ_h(data, n_outer_cuts)
    M = calc_big_M(data, ρ_h)
    lb_ls = []
    ub_dict = Dict()
    of_term1_lb = of_term2_lb = of_term3_lb = 0    
    n_vars, n_cons, n_nodes = 0, [], []

    I = 1:data.I
    J = 1:data.J
    T = 1:data.t
    K = 1:data.k

    # Initialize the bounds
    lb, ub = floor(sum([minimum(data.C[i,:,:]) for i in I])), sum(data.F)+sum(data.C)+data.D*sum(data.a)

    timer = Timer(params.max_time - elapsed(status))
    while true

        isopen(timer) || break
        yield()

        m_iter = ini_model_iter_cuts(data, params, status, ρ_h, M)

        while (ub-lb)/ub >= params.ϵ
            new_lb, of_term1_lb, of_term2_lb, of_term3_lb, yq, xq, zq, ρq, wq, Rq, n_vars, n_cons_temp, n_nodes_temp = solve_model_iter_cuts(m_iter, params, status)

            push!(n_cons,n_cons_temp)
            push!(n_nodes,n_nodes_temp)

            if (new_lb-ub)/ub > params.ϵ 
                break
            end
            if status.endStatus == :infeasible || status.endStatus == :none
                break
            end
            lb = new_lb
            push!(lb_ls,lb)
            newub, of_term1_ub, of_term2_ub, of_term3_ub = calc_ub(xq, yq, data)
            ub_dict[newub] =[of_term1_ub, of_term2_ub, of_term3_ub]
            ub = min(ub, newub)

            # Add new constraints for the linear outer approximation
            update_ρ = false
            update_M = false
            for j in J, t in T
            #     if !(ρq[j,t] in ρ_h[j,t])
            # for I in findall(==(1), Rq.data .< (ρq.data ./ (1 .- ρq.data))) 
            #     j, t = I[1], I[2]
                if Rq[j,t] + params.ϵ < ρq[j,t]/(1-ρq[j,t])
                    update_ρ = true
                    if ρq[j,t] + params.ϵ < maximum(ρ_h[j,t,:])
                        @constraint(m_iter, (1-ρq[j,t])^2 * m_iter[:R][j,t] - m_iter[:ρ][j,t] >= -ρq[j,t]^2)
                    else 
                        update_M = true
                    end
                end

                # Avoid having Nan in the M calc. This actually happens
                # if ρq > 1 - 10e-18, but we decrease the max to limit the
                # max value of M
                if ρq[j,t] > 1 - 10e-4
                    ρq[j,t] = 1 - 10e-4
                end
            end

            # Update ρ_h
            if update_ρ
                ρ_h = cat(ρ_h, ρq, dims=3)
            end

            if update_M    
                M = calc_big_M(data, ρ_h)
                for j in J, t in T, k in K
                    set_normalized_coefficient(m_iter[:cM][j,t,k], m_iter[:y][j,k], -M[j,t])
                end
            end
            
            status.nIter+=1
        end
        tests_feas = is_sol_feas(data, yq.data, xq.data)

        return lb, ub, of_term1_lb, of_term2_lb, of_term3_lb, ub_dict[ub][1], ub_dict[ub][2], ub_dict[ub][3], yq, xq, lb_ls, tests_feas, n_vars, n_cons, n_nodes, [n_outer_cuts, size(ρ_h)[3]]
    end
    return lb, ub, of_term1_lb, of_term2_lb, of_term3_lb, ub_dict[ub][1], ub_dict[ub][2], ub_dict[ub][3], yq, xq, lb_ls, tests_feas, n_vars, n_cons, n_nodes, [n_outer_cuts, size(ρ_h)[3]]
    
end


function ini_model_iter_cuts(data, params, status, ρ_h, M)
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
    H = 1:size(ρ_h,3)

    Dt = [D/sum(λ[i, t] for i in I) for t in T]
    
    maxtime = max(1, params.max_time - elapsed(status))
    m = Model(optimizer_with_attributes(Gurobi.Optimizer,
                                    "OutputFlag" => 1,
                                    "Threads" => 1,
                                    "MIPGap" => 1e-5,
                                    # "MIPFocus" => 1,
                                    "TimeLimit" => maxtime + 1,
                                    "LazyConstraints" => 1,
                                    # "Presolve" => 0,                                    
                                    # "Cuts" => 0,
                                    )
                                    )

    @variable(m, 0 <= x[I,J,T] <= 1 )
    @variable(m, y[J,K], Bin)
    
    @variable(m, 0 <= z[J,K,T] <= 1 )
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
    @constraint(m, [j in J, t in T, h in H], (1-ρ_h[j,t,h])^2 * R[j,t] - ρ[j,t] >= -ρ_h[j,t,h]^2, set_string_name = false)
    # 14 - 17
    @constraint(m, [j in J, t in T, k in K], z[j,k,t] - y[j,k] <= 0, set_string_name = false)
    @constraint(m, [j in J, t in T], ρ[j,t] <= 1-10e-7, set_string_name = false)

    @constraint(m, cM[j in J, t in T, k in K], w[j,k,t] - M[j,t]*y[j,k] <= 0, set_string_name = false)
    return m
end

function solve_model_iter_cuts(m, params, status)

    maxtime = max(1, params.max_time - elapsed(status))
    set_attribute(m, "TimeLimit", maxtime)

    optimize!(m)

    n_vars = num_variables(m)
    n_cons = num_constraints(m; count_variable_in_set_constraints = true)
    n_nodes = MOI.get(m, MOI.NodeCount())
    
    end_stat = termination_status(m)

    if end_stat == MOI.INFEASIBLE
        status.endStatus = :infeasible
        return [], [], [], [], [], [], 0        
    elseif end_stat == MOI.OPTIMAL
        status.endStatus = :optimal
    elseif end_stat == MOI.TIME_LIMIT
        status.endStatus = :tlim
    elseif end_stat == MOI.ITERATION_LIMIT
        status.endStatus = :itlim
    elseif end_stat == MOI.SOLUTION_LIMIT
        status.endStatus = :sollim        
    end

    xval = value.(m[:x])
    yval = value.(m[:y])
    zval = value.(m[:z])
    ρval = value.(m[:ρ])
    wval = value.(m[:w])
    Rval = value.(m[:R])
    optval = objective_value(m)
    op_val_term1 = value(m[:of_term1])
    op_val_term2 = value(m[:of_term2])
    op_val_term3 = value(m[:of_term3])
    return optval, op_val_term1, op_val_term2, op_val_term3, yval, xval, zval, ρval, wval, Rval, n_vars, n_cons, n_nodes    
end