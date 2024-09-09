function model_benders_iter(data, params, status, types=["Gral"])

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

        # Param for RHS perturbation Sherali
        μ = 0
        # Multiplier Fischetti
        w_fc = 0

        # Master problem
        mp = ini_mp_benders(data, params, status, types)
        # println(mp)

        # Monitor the bounds
        lb_iter = Dict()
        ub_iter = Dict()
        
        # Initial values
        int_y_ind = [data.k for j in J]
        int_y = gen_y(data, int_y_ind)

        # Initialize all the subproblems (one per period)
        if "SH" in types # =="SH"
            μ = 10e-6
        end

        if "FC" in types # =="FC"
            w_fc = 1
        end

        # Primal subproblems
        primals = Dict()

        for t in T
            primals[t] = ini_benders_sp_primal(int_y, data, ρ_h, t, μ)
        end    

        yvals = Array{Float64}(undef, data.J, data.k)
        xvals = Array{Float64}(undef, data.I, data.J, data.t)
        ρvals = Array{Float64}(undef, data.J, data.t)    
        Rvals = Array{Float64}(undef, data.J, data.t)    
        wvals = Array{Float64}(undef, data.J, data.k, data.t)
        zvals = Array{Float64}(undef, data.J, data.k, data.t)

        # yval_temp = Array{Float64}(undef, data.J, data.k)

        while (ub-lb)/ub >= 10e-1 # params.ϵ

            ### Papadakos ###
            if "PK" in types   
                if status.nIter > 0          
                    int_y = 0.5.*int_y .+ 0.5.*yvals
                end 

                for t in T
                    primal_sp = primals[t]
                    # primal_sp = update_primal_comp(primal_sp, data, int_y, M, t, μ)
                    primal_sp = update_primal(primal_sp, data, int_y, M, t, μ, ρ_h)
                    
                    # sp_stat, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, sp_val, Cterm, Congterm,  π9 = benders_sp_primal_comp(primal_sp)
                    sp_stat, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, sp_val, Cterm, Congterm, π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_primal(primal_sp, data, ρ_h, t)
                    # Add cuts for the linear approximation
                    # primal_sp = add_cuts(primal_sp, sp_ρval, ρ_h, data, t, μ)                   
                    
                    status.nOptCuts += 1
                    # cons = @build_constraint(α[t] >= sp_val + sum(π9[j,k]*(y[j,k]-int_y[j,k]) for j in J for k in K))
                    if sp_stat == MOI.OPTIMAL
                        H = 1:size(ρ_h[:,t,:],2)
                        expr=zero(AffExpr)
                        expr += sum(-π1[j]*Q[j,k]*mp[:y][j,k] for j in J for k in K)
                        expr += sum(-π4[j,k]*mp[:y][j,k] for j in J for k in K)
                        expr += sum(-π6[j,k]*M[j,t]*mp[:y][j,k] for j in J for k in K)
                        
                        expr += sum(π2[i] for i in I)
                        expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
                        expr += sum(-π10[i,j] for i in I for j in J)
                        expr += sum(-π11[j,k] for j in J for k in K)
                        expr += sum(-π12[j] for j in J)
                        # cons = @build_constraint(mp[:α][t] >= expr)
                        # MOI.submit(mp, MOI.LazyConstraint(cb), cons) 
                        mp[Symbol("PK_$(t)_$(status.nIter)" )] = @constraint(mp, mp[:α][t] >= expr, base_name="PK_$(t)_$(status.nIter)")
                        println("iter= $(status.nIter) adding Papadakos optimality cut for t=$t --------------------")
                    end
                    
                    primals[t] = primal_sp
                end
            end
            # println(mp)
            end_stat, mp_of, yvals, αvals = benders_mp(mp)
            

            if end_stat == MOI.INFEASIBLE
                println("The problem is infeasible")
                break 
            end

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
                expr += sum(-π1[j]*Q[j,k]*mp[:y][j,k] for j in J for k in K)
                expr += sum(-π4[j,k]*mp[:y][j,k] for j in J for k in K)
                expr += sum(-π6[j,k]*M[j,t]*mp[:y][j,k] for j in J for k in K)
                
                expr += sum(π2[i] for i in I)
                expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
                expr += sum(-π10[i,j] for i in I for j in J)
                expr += sum(-π11[j,k] for j in J for k in K)
                expr += sum(-π12[j] for j in J)

                # Add a feasibility cut
                if sp_stat == MOI.DUAL_INFEASIBLE || sp_stat == MOI.INFEASIBLE
                    status.nFeasCuts += 1 
                    mp[Symbol("feas_$(t)_$(status.nIter)" )] =  @constraint(mp, 0 >= expr, base_name="feas_$(t)_$(status.nIter)") 
                    println("iter= $(status.nIter) adding Benders feasibility cut for t=$t ######################")   
                    update_bounds = false 
                                
                
                elseif sp_stat == MOI.OPTIMAL
                    # Add an optimality cut
                    if  sp_val  > αvals[t] - 1e-5

                        # All but Fischetti
                        if !("FC" in types)
                            status.nOptCuts += 1
                            # cons = @build_constraint(mp[:α][t] >= expr)
                            mp[Symbol("opt_$(t)_$(status.nIter)" )] = @constraint(mp, mp[:α][t] >= expr, base_name="opt_$(t)_$(status.nIter)")
                            # MOI.submit(mp, MOI.LazyConstraint(cb), cons) 
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
                            sp_stat, sp_valMW, _π0,  π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_dual(int_y, data, ρ_h, t, yvals, αvals, sp_val)

                            expr=zero(AffExpr)
                            expr += sum(-π1[j]*Q[j,k]*mp[:y][j,k] for j in J for k in K)
                            expr += sum(-π4[j,k]*mp[:y][j,k] for j in J for k in K)
                            expr += sum(-π6[j,k]*M[j,t]*mp[:y][j,k] for j in J for k in K)
                            
                            expr += sum(π2[i] for i in I)
                            expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
                            expr += sum(-π10[i,j] for i in I for j in J)
                            expr += sum(-π11[j,k] for j in J for k in K)
                            expr += sum(-π12[j] for j in J)

                            status.nOptCuts += 1
                            @constraint(mp, mp[:α][t] >= expr) 
                            println("iter= $(status.nIter) adding Magnanti-Wong optimality cut for t=$t --------------------")
                        end                                        
                        
                        # ### Papadakos ###
                        # if "PK" in types # =="PK"                
                        #     int_y = 0.5.*int_y .+ 0.5.*yvals          
                        #     # sp_stat,sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, sp_valPK, π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_primal(int_y, data, ρ_h, t)
                        #     primal_sp = update_primal(primal_sp, data, int_y, M, t, μ)

                        #     sp_stat, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, sp_val,  π9 = benders_sp_primal(primal_sp)
                        #     # Add cuts for the linear approximation
                        #     primal_sp = add_cuts(primal_sp, sp_ρval, ρ_h, data, t, μ)                   
                            
                        #     status.nOptCuts += 1
                        #     @constraint(mp, mp[:α][t] >= sp_val + sum(π9[j,k]*(mp[:y][j,k]-int_y[j,k]) for j in J for k in K))
                        #     println("iter= $(status.nIter) adding Papadakos optimality cut for t=$t --------------------") 
                        # end
                        
                        ### Fischetti ###
                        if "FC" in types # == "FC"
                            sp_stat, sp_valMW, π0, π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_dual(yvals, data, ρ_h, t, yvals, αvals, sp_val, types, w_fc, 1)

                            expr=zero(AffExpr)
                            expr += sum(-π1[j]*Q[j,k]*mp[:y][j,k] for j in J for k in K)
                            expr += sum(-π4[j,k]*mp[:y][j,k] for j in J for k in K)
                            expr += sum(-π6[j,k]*M[j,t]*mp[:y][j,k] for j in J for k in K)
                            
                            expr += sum(π2[i] for i in I)
                            expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
                            expr += sum(-π10[i,j] for i in I for j in J)
                            expr += sum(-π11[j,k] for j in J for k in K)
                            expr += sum(-π12[j] for j in J)

                            status.nOptCuts += 1
                            cons = @build_constraint(π0*α[t] >=  expr)
                            @constraint(mp, π0*mp[:α][t] >= expr) 
                            println("iter= $(status.nIter) adding Fischetti optimality cut for t=$t LLLLLLLLLLLLLAAAAAAAAAAAAAAAAAAA")
                        end  
                    else
                        println("Done")
                        break            
                    end
                end
                ρ_k[:,t] = sp_ρval
                primals[t] = primal_sp
            end        
            ρ_h = cat(ρ_h, ρ_k, dims=3)
            
            # Update bounds   
            if update_bounds     
                lb = mp_of
                # println("MP of = $mp_of while with Vars = $(dot(F, yvals)+sum(αvals[t] for t in T))")
                ub = ub > dot(F, yvals) + cost_sp ? dot(F, yvals) + cost_sp : ub
                ub = minimum([ub, dot(F, yvals) + cost_sp])
                Cterm = dot(C, xvals)  
                Congterm = 0.5*sum(Dt[t] * (Rvals[j, t] + ρvals[j, t] + sum(cv^2 * (wvals[j, k, t] - zvals[j, k, t]) for k in K)) for j in J for t in T)
                # println("SP = $cost_sp while with vars Var cost = $(Cterm + Congterm)")
                println("iter $(status.nIter): LB = $mp_of / UB = $(dot(F, yvals) + Cterm + Congterm) /  ub = $ub ")
            end

            status.nIter+=1

            lb_iter[status.nIter] = lb
            ub_iter[status.nIter] = ub
            # println(yvals)
        end


        # println(yvals)
        Fterm = dot(F, yvals)
        Cterm = dot(C, xvals)  
        Congterm = 0.5*sum(Dt[t] * (Rvals[j, t] + ρvals[j, t] + sum(cv^2 * (wvals[j, k, t] - zvals[j, k, t]) for k in K)) for j in J for t in T)

        return Fterm+Cterm+Congterm, Fterm+Cterm+Congterm, Fterm, Cterm, Congterm, yvals, xvals, status.nIter, status.nFeasCuts, status.nOptCuts, lb_iter, ub_iter

    end
    Fterm = dot(F, yvals)
    Cterm = dot(C, xvals)  
    Congterm = 0.5*sum(Dt[t] * (Rvals[j, t] + ρvals[j, t] + sum(cv^2 * (wvals[j, k, t] - zvals[j, k, t]) for k in K)) for j in J for t in T)

    return Fterm+Cterm+Congterm, Fterm+Cterm+Congterm, Fterm, Cterm, Congterm, yvals, xvals, status.nIter, status.nFeasCuts, status.nOptCuts, lb_iter, ub_iter


end


function ini_mp_benders(data, params, status, types) 
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

function benders_mp(m)

    optimize!(m)

    end_stat = termination_status(m)
    yval = value.(m[:y])
    αval = value.(m[:α])
        
    if end_stat == MOI.DUAL_INFEASIBLE || end_stat == MOI.INFEASIBLE
        return end_stat, -1, yval, αval
    elseif end_stat == MOI.OPTIMAL
        return end_stat, objective_value(m), yval, αval
    else
        return end_stat, -1 , [], []
    end       
end

