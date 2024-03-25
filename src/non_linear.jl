#= using Alpine
using Ipopt
using Gurobi
using JuMP
using MathOptInterface =#

function minlp(data, params, status)
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
    
    #ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
    #gurobi = optimizer_with_attributes(Gurobi.Optimizer, "output_flag" => false)
    
    maxtime = max(1, params.max_time - elapsed(status))
    #= m = Model(optimizer_with_attributes(Alpine.Optimizer,
                                    "nlp_solver" => ipopt,
                                    "mip_solver" => gurobi,
                                    "minlp_solver" => ipopt,
                                    #"OutputFlag" => 0,
                                    #"Threads" => 1,
                                    #"MIPFocus" => 1,
                                    #"TimeLimit" => maxtime + 1,
                                    #"SolutionLimit" => 1
                                    )
                                    ) =#
    m = Model(optimizer_with_attributes(Gurobi.Optimizer,
                                    "OutputFlag" => 0,
                                    "Threads" => 1,
                                    #"MIPFocus" => 1,
                                    "TimeLimit" => maxtime + 1,
                                    "NonConvex" => 2,
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
    0.5*sum((D/sum(λ[i,t] for i in I))*(R[j,t] + ρ[j,t] + sum(cv^2*(w[j,k,t]-z[j,k,t]) for k in K)) for j in J for t in T))

    # Capacity cannot be exceeded and steady state has to be conserved
    for j in J, t in T
        @constraint(m, sum(λ[i,t]*x[i,j,t] for i in I) <= sum(Q[j,k]*y[j,k] for k in K))
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
        @constraint(m, sum(λ[i,t]*x[i,j,t] for i in I) - sum(Q[j,k]*z[j,k,t] for k in K) == 0)
        @constraint(m, sum(z[j,k,t] for k in K) - ρ[j,t] == 0)
        @constraint(m, R[j,t] - R[j,t]*ρ[j,t] - ρ[j,t] == 0)
        @constraint(m, sum(w[j,k,t] for k in K)-R[j,t] == 0)
    end
    # 15 - 18
    for j in J, t in T, k in K
        @constraint(m, z[j,k,t] - y[j,k] <= 0)
        @constraint(m, w[j,k,t] - M*y[j,k] <= 0)
    end 
    #write_to_file(m, "debug_nlp.lp")
    optimize!(m)
    end_stat = termination_status(m)
    if end_stat == MOI.OPTIMAL || end_stat == MOI.SOLUTION_LIMIT
        status.endStatus = :optimal
        if end_stat == MOI.SOLUTION_LIMIT
            status.endStatus = :tlim
        end
        xval = value.(x)
        yval = value.(y)
        optval = objective_value(m)
        return xval, yval, optval
    else return [], [], objective_value(m)
    end
end