function model_cuts_priori(data, params, status, ρ_h)
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
    M = data.M
    H = 1:size(ρ_h,3)

    Dt = [D/sum(λ[i, t] for i in I) for t in T]
    
    maxtime = max(1, params.max_time - elapsed(status))
    m = Model(optimizer_with_attributes(Gurobi.Optimizer,
                                    "OutputFlag" => 0,
                                    "Threads" => 1,
                                    "MIPGap" => 1e-5,
                                    #"MIPFocus" => 1,
                                    "TimeLimit" => maxtime + 1,
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
    0.5*sum(Dt[t]*(R[j,t] + ρ[j,t] + sum(cv^2*(w[j,k,t]-z[j,k,t]) for k in K)) for j in J for t in T))

    # Capacity cannot be exceeded and steady state has to be conserved
    @constraint(m, [j in J, t in T], sum(λ[i,t]*x[i,j,t] for i in I) - sum(Q[j,k]*y[j,k] for k in K) <= 0, base_name = "c3")
    
    # All customer zones need to be assigned to exactly one facility
    @constraint(m, [i in I, t in T], sum(x[i,j,t] for j in J) == 1, base_name = "c4")
    
    # At most one capacity level can be selected per facility
    @constraint(m, [j in J], sum(y[j,k] for k in K) <= 1, base_name = "c5")
    
    # 13 - 15 - 16 - 18
    @constraint(m, [j in J, t in T], sum(λ[i,t]*x[i,j,t] for i in I) - sum(Q[j,k]*z[j,k,t] for k in K) == 0, base_name = "c13")
    @constraint(m, [j in J, t in T], sum(z[j,k,t] for k in K) - ρ[j,t] == 0, base_name = "c15")
    @constraint(m, [j in J, t in T], sum(w[j,k,t] for k in K)-R[j,t] == 0, base_name = "c18")
    @constraint(m, [j in J, t in T, h in H], (1-ρ_h[j,t,h])^2 * R[j,t] - ρ[j,t] >= -ρ_h[j,t,h]^2, base_name = "c16")
    # 14 - 17
    @constraint(m, [j in J, t in T, k in K], z[j,k,t] - y[j,k] <= 0, base_name = "c14")

    M2 = calc_big_M(data, ρ_h)
    #M2 = calc_big_M_2(data, ρ_h)
    #println("Calc M =$M2")
    #println("ρ_h = $ρ_h")
    @constraint(m, [j in J, t in T, k in K], w[j,k,t] - M2[j,t]*y[j,k] <= 0, base_name = "c17")
    
    
    
    #write_to_file(m, "debug_cuts.lp")
    optimize!(m)
    end_stat = termination_status(m)
    if end_stat == MOI.OPTIMAL || end_stat == MOI.SOLUTION_LIMIT
        status.endStatus = :optimal
        if end_stat == MOI.SOLUTION_LIMIT
            status.endStatus = :tlim
        end
        xval = value.(x)
        yval = value.(y)
        zval = value.(z)
        ρval = value.(ρ)
        wval = value.(w)
        Rval = value.(R)
        optval = objective_value(m)
        return xval, yval, zval, ρval, wval, Rval, optval
    else return [], [], [], [], [], [], objective_value(m)
    end
end

