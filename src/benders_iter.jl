function model_benders_iter(data, params, status, types=["Gral"])
    Fterm = 0
    Cterm = 0
    Congterm = 0
    timer = Timer(params.max_time - elapsed(status))
    
    while true

        isopen(timer) || break
        yield()

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

        ρ_h = ini_ρ_h(data)  
        M = calc_big_M(data, ρ_h)
        Dt = [D/sum(λ[i, t] for i in I) for t in T]

        # Initialize the bounds
        lb, ub = 0, sum(data.F)+sum(data.C)+data.D*sum(data.a)

        # println("LB = $lb, UB = $ub")

        # Termination parameters
        t_lim_mp = 10
        t_lim_sp = 5

        # Monitor the bounds...
        lb_iter = Dict()
        ub_iter = Dict()
        # ... and the iterations
        n_iter_lp = 0
        n_iter_nodes = 0

        # Param for RHS perturbation Sherali
        μ = 0
        # Multiplier Fischetti
        w_fc = 0

        # Convergence criteria for the LP relax
        τ = 10e-4
        # Number of bounds to consider for calculating the Convergence
        n_bounds = 5
        last_bounds = zeros(n_bounds)
        conv = false

        # Initialize arrays for decision variables results
        yvals = Array{Float64}(undef, data.J, data.k)
        yvals_opt = Array{Float64}(undef, data.J, data.k)        
        αvals = [0 for t in T]
        xvals = []
        ρvals = []
        Rvals = []
        wvals = []
        zvals = []
        yvals_lp = []

        # Initialize linear relaxation of the master problem
        mp = ini_lp_benders(data, params, status)

        # Gurobi environments
        GRB_ENV_primals = Gurobi.Env()

        # Initialize the primal subproblems
        primals = Dict()
        int_y_ind = [data.k for j in J]
        int_y = gen_y(data, int_y_ind)
        # Initialize all the subproblems (one per period)
        for t in T
            primals[t] = ini_benders_sp_primal(int_y, data, ρ_h, t, μ, params, status,GRB_ENV_primals)
        end 

        # Initial interior point
        first_it_MW = false
        first_it_PK = true

        # Solve the LP relaxation
        println("\n---- STARTING ROOT NODE ----\n")
        while !conv && lb + 10e-5 < ub
            status.nIter += 1
            n_iter_lp += 1

            end_stat_lp, of_lp, yvals_lp, αvals_lp = benders_mp(mp)
            lb = of_lp
            println("\n /////////// Iter $(status.nIter)")            
            println("The LMP is $end_stat_lp")
            println("LB = $lb")
            last_bounds = append!(last_bounds, of_lp)[2:end]
            conv = last_bounds[end]-last_bounds[begin] <= τ*last_bounds[end] ? true : false
            
            # println("First OF $(last_bounds[begin]) / last OF $of_lp")
            # println("Convergeance is $conv because ϵ=$((last_bounds[end]-last_bounds[begin])/last_bounds[end])")
            # println("Yvals = $(yvals_lp.data)")
            # println("αvals = $(αvals_lp.data)")
            
            mp, primals, ρ_h, cost_sp, update_bounds, xvals, ρvals, Rvals, wvals, zvals, _first_it_MW, _int_y  = separate_cuts(mp, yvals_lp, αvals_lp, primals, ρ_h, μ, [], data, status)
            
            println("α=$(sum(αvals_lp)) and cost_sp = $cost_sp")
            if sum(αvals_lp) >= cost_sp
                println("It should be done")
                break
            end
            
            ub = ub > dot(F, yvals_lp) + cost_sp ? dot(F, yvals_lp) + cost_sp : ub
            println("UB = $ub and ub_Current = $(dot(F, yvals_lp) + cost_sp)")
            # println("ρ_h = $ρ_h")
            # ub = minimum([ub, dot(F, yvals_lp) + cost_sp])

            lb_iter[status.nIter] = lb
            ub_iter[status.nIter] = ub
            # if 29 <= status.nIter <= 31
            #     write_to_file(mp, "lp_relax_$(status.nIter).lp")
            # end
        end

        # print("Yvals rel = $yvals_lp")
        

        # Change the y vars as binary
        for j in J,  k in K
            set_binary(mp[:y][j,k])
        end

        set_start_value.(mp[:y], ceil.(yvals_lp))

        
        if "SH" in types # =="SH"
            μ = 10e-6
        end

        if "FC" in types # =="FC"
            w_fc = 1
        end    

        ub = sum(data.F)+sum(data.C)+data.D*sum(data.a)

        println("\n---- STARTING BRANCHING ----")
        println(" $(status.nIter) iterations so far")
        println("Starting with bounds lb=$lb and ub=$ub\n")
        while (ub-lb)/ub >= 10e-4 && lb + 10e-5 < ub && status.nIter < 125  # (ub-lb)/ub >= params.ϵ

            if "PK" in types
                τ_PK = 0.5
                if !first_it_PK
                    int_y = τ_PK.*int_y .+ (1-τ_PK).*yvals
                else
                    int_y = gen_y(data, int_y_ind) #yvals
                    first_it_PK = false
                end         
                
                for t in T 
                    sp_stat, _sp_val, _π0,  π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_dual(int_y, data, ρ_h, t, αvals)

                    if sp_stat == MOI.OPTIMAL
                    
                        status.nOptCuts += 1
                        H = 1:size(ρ_h[:,t,:],2)
                        expr = zero(AffExpr)
                        expr += sum(-π1[j]*Q[j,k]*mp[:y][j,k] for j in J for k in K)
                        expr += sum(-π4[j,k]*mp[:y][j,k] for j in J for k in K)
                        expr += sum(-π6[j,k]*M[j,t]*mp[:y][j,k] for j in J for k in K)
                        
                        expr += sum(π2[i] for i in I)
                        expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
                        expr += sum(-π10[i,j] for i in I for j in J)
                        expr += sum(-π11[j,k] for j in J for k in K)
                        expr += sum(-π12[j] for j in J)
                        mp[Symbol("PK_$(t)_$(status.nIter)" )] = @constraint(mp, mp[:α][t] >= expr, base_name="PK_$(t)_$(status.nIter)")
                        println("iter= $(status.nIter) adding Papadakos optimality cut for t=$t --------------------")
                    end
                end
            end

            end_stat, mp_of, yvals, αvals = benders_mp(mp)            

            if end_stat == MOI.INFEASIBLE
                println("The problem is infeasible")
                break 
            end
            lb = mp_of
            println("LB= $lb    UB=$ub")

            # println("Solving sp with $yvals")
            mp, primals, ρ_h, cost_sp, update_bounds, xvals, ρvals, Rvals, wvals, zvals, first_it_MW, int_y = separate_cuts(mp, yvals, αvals, primals, ρ_h, μ, types, data, status, first_it_MW, int_y, w_fc)
            # println("$αvals")
            # println(cost_sp)
            if sum(αvals) >= cost_sp
                println("It should be done")
                break
            elseif ub >= dot(F, yvals) + cost_sp
            # # Update bounds   
            # # if update_bounds
            # if ub >= dot(F, yvals) + cost_sp #&& update_bounds
                println("UB updated to $ub")
                
                
                ub = dot(F, yvals) + cost_sp
                # set_normalized_rhs(mp[:mp_ub], ub) #cost_sp)
                Fterm = dot(F, yvals)
                Cterm = dot(C, xvals)  
                Congterm = 0.5*sum(Dt[t] * (Rvals[j, t] + ρvals[j, t] + sum(cv^2 * (wvals[j, k, t] - zvals[j, k, t]) for k in K)) for j in J for t in T)                
                yvals_opt = yvals                
            end
            println("iter $(status.nIter): LB = $mp_of / UB = $(dot(F, yvals) + cost_sp) /  ub = $ub  / ϵ = $((ub-lb)/ub)\n")
            # ub = ub > dot(F, yvals) + cost_sp ? dot(F, yvals) + cost_sp : ub
            # ub = minimum([ub, dot(F, yvals) + cost_sp])
            
            # Fterm = dot(F, yvals)
            # Cterm = dot(C, xvals)  
            # Congterm = 0.5*sum(Dt[t] * (Rvals[j, t] + ρvals[j, t] + sum(cv^2 * (wvals[j, k, t] - zvals[j, k, t]) for k in K)) for j in J for t in T)
            # println("SP = $cost_sp while with vars Var cost = $(Cterm + Congterm)")

            status.nIter+=1
            n_iter_nodes += 1

            lb_iter[status.nIter] = lb
            ub_iter[status.nIter] = ub
        
            # return Fterm+Cterm+Congterm, Fterm+Cterm+Congterm, Fterm, Cterm, Congterm, yvals, xvals, status.nIter, status.nFeasCuts, status.nOptCuts, lb_iter, ub_iter, n_iter_lp, n_iter_nodes
            
            # println(yvals)
        end

        if "SH" in types                        
            M = calc_big_M(data, ρ_h)

            for t in T
                primal_sp = primals[t]
                # primal_sp = update_primal_comp(primal_sp, data, yval, M, t, 0, true)
                primal_sp = update_primal(primal_sp, data, yvals_opt, M, t, 0, ρ_h)
                
                _sp_stat, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, _sp_val, _sp_Cterm, _sp_Congterm, _π1, _π2, _π4, _π6, _π8, _π10, _π11, _π12 = benders_sp_primal(primal_sp, data, ρ_h, t)
                xvals[:,:,t] = sp_xval
                ρvals[:,t] = sp_ρval
                Rvals[:,t] = sp_Rval
                wvals[:,:,t] = sp_wval
                zvals[:,:,t] = sp_zval
            end
            Fterm = dot(F, yvals_opt)
            Cterm = dot(C, xvals)  
            Congterm = 0.5*sum(Dt[t] * (Rvals[j, t] + ρvals[j, t] + sum(cv^2 * (wvals[j, k, t] - zvals[j, k, t]) for k in K)) for j in J for t in T)        
        end

        return Fterm+Cterm+Congterm, Fterm+Cterm+Congterm, Fterm, Cterm, Congterm, yvals_opt, xvals, status.nIter, status.nFeasCuts, status.nOptCuts, lb_iter, ub_iter, n_iter_lp, n_iter_nodes
    end
    # Fterm = dot(F, yvals)
    # Cterm = dot(C, xvals)  
    # Congterm = 0.5*sum(Dt[t] * (Rvals[j, t] + ρvals[j, t] + sum(cv^2 * (wvals[j, k, t] - zvals[j, k, t]) for k in K)) for j in J for t in T)
    
    return Fterm+Cterm+Congterm, Fterm+Cterm+Congterm, Fterm, Cterm, Congterm, yvals_opt, xvals, status.nIter, status.nFeasCuts, status.nOptCuts, lb_iter, ub_iter, n_iter_lp, n_iter_nodes

end


function ini_mp_benders2(data, params, status, types) 
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
    # M = data.M

    Dt = [D/sum(λ[i, t] for i in I) for t in T]

    ρ_h = ini_ρ_h(data)  
    M = calc_big_M(data, ρ_h)
    
    maxtime = max(1, params.max_time - elapsed(status))

    # Param for RHS perturbation Sherali
    μ = 0
    # Multiplier Fischetti
    w_fc = 0

    # Master problem
    mp = Model(optimizer_with_attributes(Gurobi.Optimizer,
                                    "OutputFlag" => 0,
                                    "Threads" => 1,
                                    "MIPGap" => 1e-5,
                                    "MIPFocus" => 2,
                                    "TimeLimit" => maxtime + 1,
                                    "LazyConstraints" => 1,
                                    "Presolve" => 0,
                                    )
                                    )

    
    @variable(mp, y[J,K], Bin)
    @variable(mp, α[T] >= 0)

    @expression(mp, of_mp, sum(F[j,k]*y[j,k] for j in J for k in K))
    @expression(mp, α_mp, sum(α[t] for t in T))

    @objective(mp, Min, of_mp + α_mp)

    # @objective(mp, Min, sum(F[j,k]*y[j,k] for j in J for k in K) + sum(α[t] for t in T))

    # At most one capacity level can be selected per facility
    @constraint(mp, [j in J], sum(y[j,k] for k in K) <= 1, base_name = "c_mp_1")

    ## Valid inequalities
    # The total available capacity needs to be enough to handle the demand at all periods    
    ### @constraint(mp,  sum(Q[j,k]*y[j,k] for j in J for k in K) >= maximum(sum(λ,dims=1)))
    # Better by time
    @constraint(mp, [t in T], sum(Q[j,k]*y[j,k] for j in J for k in K) >= sum(λ[:,t]), base_name="c_mp_2")

    # Since we are allocating all the demand, the variable costs have to be greatter or equal than
    # the cost of allocating each node to the least expensive facility (closest one) plus
    # the congestion that each demand node would produce by itself if it was the onlye one allocated
    # to the facility with greatter capacity
    max_q = maximum(Q)
    min_con = minimum(λ, dims=1)
    # min_con = []
    # for t in T
    #     min_con_t = 0
    #     for i in I
    #         min_con_t += ((1+cv^2)/2) * (λ[i,t]^2/(max_q*(max_q-λ[i,t]))) + (λ[i,t]/max_q)
    #     end
    #     append!(min_con, min_con_t)
    # end
    @constraint(mp, [t in T], α[t] >= floor(sum([minimum(C[i,:,t]) for i in I]))  + 
    floor(Dt[t]*(((1+cv^2)/2) * (min_con[t]^2/(max_q*(max_q-min_con[t]))) + (min_con[t]/max_q))), base_name="c_mp_3")

    # for t in T
    #     println("For t=$t the lb is $(floor(sum([minimum(C[i,:,t]) for i in I]))) + $(floor(Dt[t]*min_con[t]))")
    # end

    ρ = Array{Float64}(undef, data.J, data.t)

    int_y_ind = [data.k for j in J]
    int_y = gen_y(data, int_y_ind)

    ### Warmstart (WS) ####
    if "WS1" in types         
        of, of_term1, of_term2, of_term3, y_ind_ws, best_x = heur_local_search_first_best(data, params, status, 50, 3)
        y_ws = gen_y(data, y_ind_ws)  
        set_start_value.(y, y_ws)
    elseif "WS2" in types
        all_sols=Dict()
        of, of_term1, of_term2, of_term3, y_ind_ws, best_x =  heur_nelder_mead_y(data, params, status, ini_y(data), all_sols, 30, 3)
        y_ws = gen_y(data, y_ind_ws)  
        set_start_value.(y, y_ws)
    end

    return mp
end

function benders_mp2(m)
    optimize!(m)
    end_stat = termination_status(m)
    yval = value.(m[:y])
    αval = value.(m[:α])
        
    if end_stat == MOI.DUAL_INFEASIBLE || end_stat == MOI.INFEASIBLE
        return end_stat, 10e8, yval, αval
    elseif end_stat == MOI.OPTIMAL
        return end_stat, objective_value(m), yval, αval
    else
        return end_stat, 10e8 , [], []
    end       
end

function ini_lp_benders2(data, params, status)
    I = 1:data.I
    J = 1:data.J
    λ = data.a
    C = data.C
    F = data.F
    Q = data.Q
    K = 1:data.k
    T = 1:data.t
    
    maxtime = max(1, params.max_time - elapsed(status))

    # Linear relaxation of the master problem
    lp = Model(optimizer_with_attributes(Gurobi.Optimizer,
                                    "OutputFlag" => 0,
                                    "Threads" => 1,
                                    "MIPGap" => 1e-5,
                                    "TimeLimit" => maxtime + 1,
                                    "Presolve" => 0,
                                    )
                                    )
    # lp = Model(optimizer_with_attributes(CPLEX.Optimizer,
    #                                 "CPXPARAM_ScreenOutput" => 0,
    #                                 "CPXPARAM_Threads" => 1,
    #                                 "CPXPARAM_MIP_Tolerances_MIPGap" => 1e-5,
    #                                 "CPXPARAM_TimeLimit" => maxtime + 1,
    #                                 "CPXPARAM_Preprocessing_Presolve" => 0,
    #                                 )
    #                                 )
                                    

    
    @variable(lp, 0 <= y[J,K] <= 1)
    @variable(lp, α[T] >= 0)

    @expression(lp, of_lp, sum(F[j,k]*y[j,k] for j in J for k in K))
    @expression(lp, α_lp, sum(α[t] for t in T))

    @objective(lp, Min, of_lp + α_lp)

    # At most one capacity level can be selected per facility
    @constraint(lp, [j in J], sum(y[j,k] for k in K) <= 1, base_name = "c_lp_1")

    ## Valid inequalities
    # The capacity must be enough to cover all the demand in each period
    @constraint(lp, [t in T], sum(Q[j,k]*y[j,k] for j in J for k in K) >= sum(λ[:,t]), base_name="c_lp_2")

    # Since we are allocating all the demand, the variable costs have to be greatter than or equal to
    # the cost of allocating each node to the least expensive facility (closest one)
    @constraint(lp, [t in T], α[t] >= floor(sum([minimum(C[i,:,t]) for i in I])), base_name="c_lp_3")

    # @constraint(lp, mp_α_ub, sum(α[t] for t in T) <= sum(data.F)+sum(data.C)+data.D*sum(data.a))
    # @constraint(lp, mp_ub, of_lp + α_lp <= sum(data.F)+sum(data.C)+data.D*sum(data.a))

    return lp
end

function separate_cuts2(m, yvals, αvals, primals, ρ_h, μ, types, data, status, first_it_MW=false, int_y=[], w_fc=0)
    xvals = Array{Float64}(undef, data.I, data.J, data.t)
    ρvals = Array{Float64}(undef, data.J, data.t)    
    Rvals = Array{Float64}(undef, data.J, data.t)    
    wvals = Array{Float64}(undef, data.J, data.k, data.t)
    zvals = Array{Float64}(undef, data.J, data.k, data.t)
    
    I = 1:data.I
    J = 1:data.J
    Q = data.Q
    K = 1:data.k
    T = 1:data.t
    
    ρ_k = Array{Float64}(undef,data.J,data.t)                    
    M = calc_big_M(data, ρ_h) 

    cost_sp = 0
    update_bounds = true

    for t in T
        ##### Classical Benders ######
        # Just one model
        primal_sp = primals[t]
        primal_sp = update_primal(primal_sp, data, yvals, M, t, μ, ρ_h)
        
        sp_stat, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, sp_val, Cterm, Congterm, π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_primal(primal_sp, data, ρ_h, t)
        cost_sp += sp_val

        xvals[:,:,t] = sp_xval
        ρvals[:,t] = sp_ρval
        Rvals[:,t] = sp_Rval
        wvals[:,:,t] = sp_wval
        zvals[:,:,t] = sp_zval
        

        H = 1:size(ρ_h[:,t,:],2)
        expr=zero(AffExpr)
        expr += sum(-π1[j]*Q[j,k]*m[:y][j,k] for j in J for k in K)
        expr += sum(-π4[j,k]*m[:y][j,k] for j in J for k in K)
        expr += sum(-π6[j,k]*M[j,t]*m[:y][j,k] for j in J for k in K)
        
        expr += sum(π2[i] for i in I)
        expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
        expr += sum(-π10[i,j] for i in I for j in J)
        expr += sum(-π11[j,k] for j in J for k in K)
        expr += sum(-π12[j] for j in J)
        # expr = sum(-π1[j]*Q[j,k]*m[:y][j,k] for j in J for k in K)+
        # sum(-π4[j,k]*m[:y][j,k] for j in J for k in K)+
        # sum(-π6[j,k]*M[j,t]*m[:y][j,k] for j in J for k in K)+
        # sum(π2[i] for i in I)+
        # sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)+
        # sum(-π10[i,j] for i in I for j in J)+
        # sum(-π11[j,k] for j in J for k in K)+
        # sum(-π12[j] for j in J)

        # Add a feasibility cut
        if sp_stat == MOI.DUAL_INFEASIBLE || sp_stat == MOI.INFEASIBLE
            status.nFeasCuts += 1 
            m[Symbol("feas_$(t)_$(status.nIter)")] =  @constraint(m, 0 >= expr, base_name="feas_$(t)_$(status.nIter)") 
            println("iter= $(status.nIter) adding Benders feasibility cut for t=$t ######################")   
            update_bounds = false                            
        
        elseif sp_stat == MOI.OPTIMAL
            # Add an optimality cut
            if  sp_val  > αvals[t] + 1e-4  
                # All but Fischetti
                if !("FC" in types)
                    status.nOptCuts += 1
                    m[Symbol("opt_$(t)_$(status.nIter)")] = @constraint(m, m[:α][t] >= expr, base_name="opt_$(t)_$(status.nIter)")
                    # println("Cut\n")
                    # println(expr)
                    println("iter= $(status.nIter) adding Benders optimality cut for t=$t --------------------") 
                end

                ### Magnanti-Wong ###
                if "MW" in types # =="MW"
                    
                    # ### Magnanti-Wong ###  
                    if first_it_MW
                        int_y = 0.5.*int_y .+ 0.5.*yvals
                    else
                        int_y = yvals
                        first_it_MW = true
                    end            
                    sp_stat, _sp_val, _π0,  π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_dual(int_y, data, ρ_h, t, αvals) #, types)

                    expr=zero(AffExpr)
                    expr += sum(-π1[j]*Q[j,k]*m[:y][j,k] for j in J for k in K)
                    expr += sum(-π4[j,k]*m[:y][j,k] for j in J for k in K)
                    expr += sum(-π6[j,k]*M[j,t]*m[:y][j,k] for j in J for k in K)
                    
                    expr += sum(π2[i] for i in I)
                    expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
                    expr += sum(-π10[i,j] for i in I for j in J)
                    expr += sum(-π11[j,k] for j in J for k in K)
                    expr += sum(-π12[j] for j in J)

                    status.nOptCuts += 1
                    @constraint(m, m[:α][t] >= expr) 
                    println("iter= $(status.nIter) adding Magnanti-Wong optimality cut for t=$t --------------------")
                end                                        
                
                ### Fischetti ###
                if "FC" in types
                    sp_stat, _sp_val, π0, π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_dual(yvals, data, ρ_h, t, αvals, w_fc, 1)

                    expr=zero(AffExpr)
                    expr += sum(-π1[j]*Q[j,k]*m[:y][j,k] for j in J for k in K)
                    expr += sum(-π4[j,k]*m[:y][j,k] for j in J for k in K)
                    expr += sum(-π6[j,k]*M[j,t]*m[:y][j,k] for j in J for k in K)
                    
                    expr += sum(π2[i] for i in I)
                    expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
                    expr += sum(-π10[i,j] for i in I for j in J)
                    expr += sum(-π11[j,k] for j in J for k in K)
                    expr += sum(-π12[j] for j in J)

                    status.nOptCuts += 1
                    # cons = @build_constraint(π0* m[:α][t] >=  expr)
                    @constraint(m, π0*m[:α][t] >= expr) 
                    println("iter= $(status.nIter) adding Fischetti optimality cut for t=$t --------------------")
                end
            else
                println("Sol is ok for $t")
            end
        end
        ρ_k[:,t] = sp_ρval
        primals[t] = primal_sp
    end
    ρ_h = cat(ρ_h, ρ_k, dims=3)
    return m, primals, ρ_h, cost_sp, update_bounds, xvals, ρvals, Rvals, wvals, zvals, first_it_MW, int_y
end