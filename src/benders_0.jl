function model_benders(data, params, status, types=["Gral"]) # ini_model_benders(data, params, status, types=["Gral"])
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

    Dt = [D/sum(λ[:, t]) for t in T]

    ρ_h = ini_ρ_h(data)  
    M = calc_big_M(data, ρ_h)
    
    maxtime = max(1, params.max_time - elapsed(status))

    # Monitor the bounds
    lb_iter = Dict()
    ub_iter = Dict()

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
    @constraint(mp, [t in T], sum(Q[j,k]*y[j,k] for j in J for k in K) >= sum(λ[:,t]))

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
    @constraint(mp, [t in T], α[t] >= floor(sum([minimum(C[i,:,t]) for i in I])))

    # for t in T
    #     println("For t=$t the lb is $(floor(sum([minimum(C[i,:,t]) for i in I]))) + $(Dt[t]*min_con[t])")
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

    # Just for one model
    primals = Dict()

    # Initialize all the subproblems (one per period)
    if "SH" in types # =="SH"
        μ = 10e-6
    end

    if "FC" in types # =="FC"
        w_fc = 1
    end

    for t in T
        # primals[t] = ini_benders_sp_primal_comp(int_y, data, ρ_h, t, μ)
        primals[t] = ini_benders_sp_primal(int_y, data, ρ_h, t, μ, params, status)
    end   
    

    ### Papadakos 1st iteration ####
    if "PK" in types # =="PK"
        M = calc_big_M(data, ρ_h)
        ρ_k = Array{Float64}(undef,data.J,data.t)
        
        for t in T
            primal_sp = primals[t]
            sp_stat, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, sp_val, Cterm, Congterm, π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_primal(primal_sp, data, ρ_h, t)           
            # _sp_stat, _sp_xval, sp_ρval, _sp_Rval, _sp_wval, _sp_zval, sp_val, Cterm, Congterm,  π9 = benders_sp_primal_comp(primal_sp)

            # Add cuts for the linear approximation
            # primal_sp = add_cuts(primal_sp, sp_ρval, ρ_h, data, t, μ)
                
            println("iter= $(status.nIter) adding Papadakos initial cut for t=$t @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")      
            status.nOptCuts += 1
            # @constraint(mp, α[t] >= sp_val + sum(π9[j,k]*(y[j,k]-int_y[j,k]) for j in J for k in K))
            H = 1:size(ρ_h[:,t,:],2)
            expr=zero(AffExpr)
            # expr += sp_val + sum(π9[j,k]*(y[j,k]-yvals[j,k]) for j in J for k in K)
            expr += sum(-π1[j]*Q[j,k]*y[j,k] for j in J for k in K)
            expr += sum(-π4[j,k]*y[j,k] for j in J for k in K)
            expr += sum(-π6[j,k]*M[j,t]*y[j,k] for j in J for k in K)
            
            expr += sum(π2[i] for i in I)
            expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
            expr += sum(-π10[i,j] for i in I for j in J)
            expr += sum(-π11[j,k] for j in J for k in K)
            expr += sum(-π12[j] for j in J)
            @constraint(mp, α[t] >= expr)
            ρ_k[:,t] = sp_ρval        
        end
        ρ_h = cat(ρ_h, ρ_k, dims=3)
    end

    first_it_MW = false

    function lazycb(cb)
        if callback_node_status(cb, mp) != MOI.CALLBACK_NODE_STATUS_INTEGER
            return
        end        

        yvals = callback_value.(cb, y)
        αvals = callback_value.(cb, α)

        status.nIter += 1        

        ρ_k = Array{Float64}(undef,data.J,data.t)
                    
        M = calc_big_M(data, ρ_h)  

        xval = Array{Float64}(undef, data.I, data.J, data.t)
        ρval = Array{Float64}(undef, data.J, data.t)    
        Rval = Array{Float64}(undef, data.J, data.t)    
        wval = Array{Float64}(undef, data.J, data.k, data.t)
        zval = Array{Float64}(undef, data.J, data.k, data.t)

        add_bounds = true
        
        for t in T   

            ##### Classical Benders ######
            # Just one model
            primal_sp = primals[t]
            # primal_sp  = update_primal_comp(primal_sp, data, yvals, M, t, μ)
            primal_sp = update_primal(primal_sp, data, yvals, M, t, μ, ρ_h)
            
            sp_stat, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, sp_val, Cterm, Congterm, π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_primal(primal_sp, data, ρ_h, t)
            # println("SP of val: $sp_val o $(Cterm+Congterm)")
            # sp_stat, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, sp_val, Cterm, Congterm,  π9 = benders_sp_primal_comp(primal_sp)
            # Add cuts for the linear approximation
            # primal_sp = add_cuts(primal_sp, sp_ρval, ρ_h, data, t, μ)

            xval[:,:,t] = sp_xval
            ρval[:,t] = sp_ρval
            Rval[:,t] = sp_Rval
            wval[:,:,t] = sp_wval
            zval[:,:,t] = sp_zval
            
            H = 1:size(ρ_h[:,t,:],2)
            # println("Start iter $(status.nIter), H  = $H ")
            expr=zero(AffExpr)
            # expr += sp_val + sum(π9[j,k]*(y[j,k]-yvals[j,k]) for j in J for k in K)
            expr += sum(-π1[j]*Q[j,k]*y[j,k] for j in J for k in K)
            expr += sum(-π4[j,k]*y[j,k] for j in J for k in K)
            expr += sum(-π6[j,k]*M[j,t]*y[j,k] for j in J for k in K)
            
            expr += sum(π2[i] for i in I)
            expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
            expr += sum(-π10[i,j] for i in I for j in J)
            expr += sum(-π11[j,k] for j in J for k in K)
            expr += sum(-π12[j] for j in J)
            # Add a feasibility cut
            # println("The sp of is $sp_val")
            # println("The sp is $(sp_stat)")
            # println("y : $(yvals.data)")
            # println("x : $(sp_xval.data)")
            if sp_stat == MOI.DUAL_INFEASIBLE || sp_stat == MOI.INFEASIBLE # Must include both cases
                status.nFeasCuts += 1
                # println("Expresion is $expr")
                feas_cuts = @build_constraint(0 >= expr)  
                MOI.submit(mp, MOI.LazyConstraint(cb), feas_cuts)
                println("iter= $(status.nIter) adding Benders feasibility cut for t=$t ######################")   
                add_bounds = false             
            
            elseif sp_stat == MOI.OPTIMAL

                # Add an optimality cut
                if  sp_val  > αvals[t] + 10e-5
                    
                    # All but Fischetti
                    if !("FC" in types) # !="FC"
                        status.nOptCuts += 1
                        # cons = @build_constraint(α[t] >= sp_val + sum(π9[j,k]*(y[j,k]-yvals[j,k]) for j in J for k in K))
                        cons = @build_constraint(α[t] >= expr)
                        MOI.submit(mp, MOI.LazyConstraint(cb), cons) 
                        println("iter= $(status.nIter) adding Benders optimality cut for t=$t --------------------") 
                    end
                    
                    ### Magnanti-Wong ###
                    if "MW" in types # =="MW"
                        
                        if first_it_MW
                            int_y = 0.5.*int_y .+ 0.5.*yvals
                        else
                            int_y = yvals
                            first_it_MW = true
                        end             
                        sp_stat, sp_valMW, _π0,  π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_dual(int_y, data, ρ_h, t, yvals, αvals, sp_val)

                        expr=zero(AffExpr)
                        expr += sum(-π1[j]*Q[j,k]*y[j,k] for j in J for k in K)
                        expr += sum(-π4[j,k]*y[j,k] for j in J for k in K)
                        expr += sum(-π6[j,k]*M[j,t]*y[j,k] for j in J for k in K)
                        
                        expr += sum(π2[i] for i in I)
                        expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
                        expr += sum(-π10[i,j] for i in I for j in J)
                        expr += sum(-π11[j,k] for j in J for k in K)
                        expr += sum(-π12[j] for j in J)

                        status.nOptCuts += 1
                        cons = @build_constraint(α[t] >=  expr)
                        MOI.submit(mp, MOI.LazyConstraint(cb), cons) 
                        println("iter= $(status.nIter) adding Magnanti-Wong optimality cut for t=$t --------------------")
                    end                                        
                    
                    ### Papadakos ###
                    if "PK" in types # =="PK"                
                        int_y = 0.5.*int_y .+ 0.5.*yvals   
                        # primal_sp = update_primal_comp(primal_sp, data, int_y, M, t, μ)
                        primal_sp = update_primal(primal_sp, data, int_y, M, t, μ, ρ_h)

                        # sp_stat, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, sp_val, Cterm, Congterm,  π9 = benders_sp_primal_comp(primal_sp)
                        sp_stat, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, sp_val, Cterm, Congterm, π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_primal(primal_sp, data, ρ_h, t)
                        # Add cuts for the linear approximation
                        # primal_sp = add_cuts(primal_sp, sp_ρval, ρ_h, data, t, μ)                   
                        
                        status.nOptCuts += 1
                        # cons = @build_constraint(α[t] >= sp_val + sum(π9[j,k]*(y[j,k]-int_y[j,k]) for j in J for k in K))
                        expr=zero(AffExpr)
                        expr += sum(-π1[j]*Q[j,k]*y[j,k] for j in J for k in K)
                        expr += sum(-π4[j,k]*y[j,k] for j in J for k in K)
                        expr += sum(-π6[j,k]*M[j,t]*y[j,k] for j in J for k in K)
                        
                        expr += sum(π2[i] for i in I)
                        expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
                        expr += sum(-π10[i,j] for i in I for j in J)
                        expr += sum(-π11[j,k] for j in J for k in K)
                        expr += sum(-π12[j] for j in J)
                        cons = @build_constraint(α[t] >= expr)
                        MOI.submit(mp, MOI.LazyConstraint(cb), cons) 
                        println("iter= $(status.nIter) adding Papadakos optimality cut for t=$t --------------------") 
                    end
                    
                    ### Fischetti ###
                    if "FC" in types # == "FC"
                        sp_stat, sp_valMW, π0, π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_dual(yvals, data, ρ_h, t, yvals, αvals, sp_val, types, w_fc, 1)

                        expr=zero(AffExpr)
                        expr += sum(-π1[j]*Q[j,k]*y[j,k] for j in J for k in K)
                        expr += sum(-π4[j,k]*y[j,k] for j in J for k in K)
                        expr += sum(-π6[j,k]*M[j,t]*y[j,k] for j in J for k in K)
                        
                        expr += sum(π2[i] for i in I)
                        expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
                        expr += sum(-π10[i,j] for i in I for j in J)
                        expr += sum(-π11[j,k] for j in J for k in K)
                        expr += sum(-π12[j] for j in J)

                        status.nOptCuts += 1
                        cons = @build_constraint(π0*α[t] >=  expr) #sp_valMW +
                        MOI.submit(mp, MOI.LazyConstraint(cb), cons) 
                        println("iter= $(status.nIter) adding Fischetti optimality cut for t=$t LLLLLLLLLLLLLAAAAAAAAAAAAAAAAAAA")
                    end
                else
                    add_bounds = false
                end
            end
            ρ_k[:,t] = sp_ρval

            # # Add the new constraints for the linear outter approximation
            # delete.(primal_sp, primal_sp[:c8_sp])            
            # unregister(primal_sp, :c8_sp)          
            # @constraint(primal_sp, c8_sp[j in J, h in H], (1-ρ_h_t[j,h])^2 * primal_sp[:R][j] - primal_sp[:ρ][j] >= -ρ_h_t[j,h]^2, base_name = "c8_sp")
            primals[t] = primal_sp

            ### Papadakos ###
            # if "PK" in types # =="PK"                
            #     int_y = 0.5.*int_y .+ 0.5.*yvals   
            #     # primal_sp = update_primal_comp(primal_sp, data, int_y, M, t, μ)
            #     primal_sp = update_primal(primal_sp, data, int_y, M, t, μ, ρ_h)

            #     # sp_stat, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, sp_val, Cterm, Congterm,  π9 = benders_sp_primal_comp(primal_sp)
            #     sp_stat, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, sp_val, Cterm, Congterm, π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_primal(primal_sp, data, ρ_h, t)
            #     # Add cuts for the linear approximation
            #     # primal_sp = add_cuts(primal_sp, sp_ρval, ρ_h, data, t, μ)                   
                
            #     status.nOptCuts += 1
            #     # cons = @build_constraint(α[t] >= sp_val + sum(π9[j,k]*(y[j,k]-int_y[j,k]) for j in J for k in K))
            #     expr=zero(AffExpr)
            #     expr += sum(-π1[j]*Q[j,k]*y[j,k] for j in J for k in K)
            #     expr += sum(-π4[j,k]*y[j,k] for j in J for k in K)
            #     expr += sum(-π6[j,k]*M[j,t]*y[j,k] for j in J for k in K)
                
            #     expr += sum(π2[i] for i in I)
            #     expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
            #     expr += sum(-π10[i,j] for i in I for j in J)
            #     expr += sum(-π11[j,k] for j in J for k in K)
            #     expr += sum(-π12[j] for j in J)
            #     cons = @build_constraint(α[t] >= expr)
            #     MOI.submit(mp, MOI.LazyConstraint(cb), cons) 
            #     println("iter= $(status.nIter) adding Papadakos optimality cut for t=$t --------------------") 
            # end
            
        end  
        # if add_bounds
        #     Fterm_i = dot(F, yvals)
        #     Cterm_i = dot(C, xval)   
        #     Congterm_i = 0.5 * sum(Dt[t] * (Rval[j, t] + ρval[j, t] + sum(cv^2 * (wval[j, k, t] - zval[j, k, t]) for k in K)) for j in J for t in T)
        #     upper_bound = Fterm_i + Cterm_i + Congterm_i
        #     ub_iter[status.nIter] = upper_bound

        #     lower_bound = callback_value.(cb, of_mp) + callback_value.(cb, α_mp)
        #     println("UB : $(upper_bound)")
        #     println("LB : $(lower_bound)")
        #     if length(ub_iter)==0
        #         lb_iter[status.nIter] = lower_bound
        #     elseif minimum(ub_iter)[2] >= lower_bound
        #         lb_iter[status.nIter] = lower_bound
        #     else
        #         lb_iter[status.nIter] = lb_iter[status.nIter-1]
        #     end
        #     # push!(lb_iter,lower_bound)
        # else
        #     println("Unfeas")
        #     lb_iter[status.nIter] = lb_iter[status.nIter-1]
        #     ub_iter[status.nIter] = ub_iter[status.nIter-1]
        # end
        
        # push!(ub_iter, upper_bound)      
        ρ_h = cat(ρ_h, ρ_k, dims=3)
        # println("End iter $(status.nIter), ρ_h = $ρ_h")
        # status.nIter += 1           
    end

    MOI.set(mp, MOI.LazyConstraintCallback(), lazycb)

    optimize!(mp)
    end_stat = termination_status(mp)
    if end_stat == MOI.INFEASIBLE
        status.endStatus = :infeasible
        return 0, 0

    elseif end_stat == MOI.OPTIMAL
        status.endStatus = :optimal
    elseif end_stat == MOI.TIME_LIMIT
        status.endStatus = :tlim
    elseif end_stat == MOI.ITERATION_LIMIT
        status.endStatus = :itlim
    elseif end_stat == MOI.SOLUTION_LIMIT
        status.endStatus = :sollim
    end
    
    n_nodes =  MOI.get(mp, MOI.NodeCount())

    yval = value.(mp[:y])
    xval = Array{Float64}(undef, data.I, data.J, data.t)
    ρval = Array{Float64}(undef, data.J, data.t)    
    Rval = Array{Float64}(undef, data.J, data.t)    
    wval = Array{Float64}(undef, data.J, data.k, data.t)
    zval = Array{Float64}(undef, data.J, data.k, data.t)
                
    M = calc_big_M(data, ρ_h)

    for t in T
        primal_sp = primals[t]
        # primal_sp = update_primal_comp(primal_sp, data, yval, M, t, 0, true)
        primal_sp = update_primal(primal_sp, data, yval, M, t, 0, ρ_h)
        
        _1, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, _2, _Ct, _Cont, _π1, _π2, _π4, _π6, _π8, _π10, _π11, _π12 = benders_sp_primal(primal_sp, data, ρ_h, t)
        # _1, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, _2, _Ct, _Cont, _π9 = benders_sp_primal_comp(primal_sp)
        # println("solve primal benders ok")
        xval[:,:,t] = sp_xval
        ρval[:,t] = sp_ρval
        Rval[:,t] = sp_Rval
        wval[:,:,t] = sp_wval
        zval[:,:,t] = sp_zval
    end

    Fterm = dot(F, yval)
    Cterm = dot(C, xval) #sum(C[i,j,t]*xval[i,j,t] for i in I for j in J for t in T)  
    Congterm = 0.5 * sum(Dt[t] * (Rval[j, t] + ρval[j, t] + sum(cv^2 * (wval[j, k, t] - zval[j, k, t]) for k in K)) for j in J for t in T)
    
    return objective_value(mp), Fterm+Cterm+Congterm, Fterm, Cterm, Congterm, yval, xval, n_nodes, status.nFeasCuts, status.nOptCuts, lb_iter, ub_iter
end

# function model_benders(data, params, status, types=["Gral"])
    
#     I = 1:data.I
#     J = 1:data.J
#     λ = data.a
#     C = data.C
#     F = data.F
#     cv = data.cv
#     D = data.D    
#     K = 1:data.k
#     T = 1:data.t

#     mp, ρ_h, primals, M = ini_model_benders(data, params, status, types)

#     optimize!(mp)   
        
#     end_stat = termination_status(mp)
#     if end_stat == MOI.INFEASIBLE
#         status.endStatus = :infeasible
#         return 0, 0

#     elseif end_stat == MOI.OPTIMAL
#         status.endStatus = :optimal
#     elseif end_stat == MOI.TIME_LIMIT
#         status.endStatus = :tlim
#     elseif end_stat == MOI.ITERATION_LIMIT
#         status.endStatus = :itlim
#     elseif end_stat == MOI.SOLUTION_LIMIT
#         status.endStatus = :sollim
#     end

#     n_nodes =  MOI.get(mp, MOI.NodeCount())

#     yval = value.(mp[:y])
    

#     xval = Array{Float64}(undef, data.I, data.J, data.t)
#     ρval = Array{Float64}(undef, data.J, data.t)    
#     Rval = Array{Float64}(undef, data.J, data.t)    
#     wval = Array{Float64}(undef, data.J, data.k, data.t)
#     zval = Array{Float64}(undef, data.J, data.k, data.t)
                
#     # M = calc_big_M(data, ρ_h)

#     for t in T
#         primal_sp = primals[t]
#         primal_sp = update_primal_comp(primal_sp, data, yval, M, t, 0, true)
        
#         # _1, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, _2, _π1, _π2, _π4, _π6, _π8, _π10, _π11, _π12 = benders_sp_primal(yval, data, ρ_h, t)
#         # _1, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, _2, _π9 = benders_sp_primal(yval, data, ρ_h, t)
#         _1, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, _2, _π9 = benders_sp_primal(primal_sp)
#         println("solve primal benders ok")
#         xval[:,:,t] = sp_xval
#         ρval[:,t] = sp_ρval
#         Rval[:,t] = sp_Rval
#         wval[:,:,t] = sp_wval
#         zval[:,:,t] = sp_zval

#         # xval[:,:,t] = value.(primal_sp[:x])
#         # ρval[:,t] = value.(primal_sp[:ρ])
#         # Rval[:,t] = value.(primal_sp[:R])
#         # wval[:,:,t] = value.(primal_sp[:w])
#         # zval[:,:,t] = value.(primal_sp[:z])
#     end

    
#     Dt = [D/sum(λ[i, t] for i in I) for t in T]

#     Fterm = dot(F, yval)
#     Cterm = dot(C, xval)  
#     Congterm = 0.5 * sum(Dt[t] * (Rval[j, t] + ρval[j, t] + sum(cv^2 * (wval[j, k, t] - zval[j, k, t]) for k in K)) for j in J for t in T)
    
#     return objective_value(mp), Fterm+Cterm+Congterm, Fterm, Cterm, Congterm, yval, xval, n_nodes, status.nFeasCuts, status.nOptCuts 
# end

function update_primal_comp(primal_model, data, yvals, M, t, μ, all=false)
    # primal_model = primals[t]
    for j in 1:data.J
        if all
            # set_normalized_rhs(primal_model[:c1_sp][j], -sum(data.Q[j,k]*yvals[j,k] for k in 1:data.k) - μ*sum(data.Q[j,k] for k in 1:data.k))
            set_normalized_rhs(primal_model[:c1_sp][j],  - μ*sum(data.Q[j,k] for k in 1:data.k))
        end
        for k in 1:data.k
            # set_normalized_rhs(primal_model[:c4_sp][j,k], -yvals[j,k] - μ)
            # set_normalized_rhs(primal_model[:c6_sp][j,k], -M[j,t]*yvals[j,k] - μ*M[j,t])
            if all
                set_normalized_rhs(primal_model[:c4_sp][j,k],  - μ)
                set_normalized_rhs(primal_model[:c6_sp][j,k], - μ*M[j,t])
            end

            if !(all) && μ != 0
                set_normalized_rhs(primal_model[:c6_sp][j,k], - μ*M[j,t])
            end
            set_normalized_coefficient(primal_model[:c6_sp][j,k], primal_model[:y][j,k], M[j,t])
            set_normalized_rhs(primal_model[:c9_sp][j,k], yvals[j,k])
        end
    end
    # primals[t] = primal_model
    return primal_model
end

# function update_primal(primal_model, data, yvals, M, t, μ, ρ_h, all=false)
#     for j in 1:data.J
#         if all
#             set_normalized_rhs(primal_model[:c1_sp][j], -sum(data.Q[j,k]*yvals[j,k] for k in 1:data.k) - μ*sum(data.Q[j,k] for k in 1:data.k))
#         end
#         for k in 1:data.k            
#             if all
#                 set_normalized_rhs(primal_model[:c4_sp][j,k], -yvals[j,k] - μ)
#                 set_normalized_rhs(primal_model[:c6_sp][j,k], -M[j,t]*yvals[j,k] - μ*M[j,t])
#             end

#             if !(all) && μ != 0
#                 set_normalized_rhs(primal_model[:c6_sp][j,k], -M[j,t]*yvals[j,k] - μ*M[j,t])
#             end
#         end
#     end

#     if size(ρ_h[:,t,:],2) > 9
#         primal_model[Symbol("c8_$(size(ρ_h[:,t,:],2))")] = @constraint(primal_model, [j in 1:data.J], (1-ρ_h[j,t,end])^2 * primal_model[:R][j] - primal_model[:ρ][j] >= -(1+μ)*ρ_h[j,t,end]^2)
#     end
#     return primal_model
# end
function update_primal(primal_model, data, yvals, M, t, μ, ρ_h)
    for j in 1:data.J
        set_normalized_rhs(primal_model[:c1_sp][j], -sum(data.Q[j,k]*yvals[j,k] for k in 1:data.k) - μ*sum(data.Q[j,k] for k in 1:data.k))
        for k in 1:data.k            
            set_normalized_rhs(primal_model[:c4_sp][j,k], -yvals[j,k] - μ)
            set_normalized_rhs(primal_model[:c6_sp][j,k], -M[j,t]*yvals[j,k] - μ*M[j,t])
        end
    end

    last_h = size(ρ_h[:,t,:],2)
    if last_h > 9        
        # println("Adding the extra 8 cuts")
        primal_model[Symbol("c8_$(last_h)")] = @constraint(primal_model, [j in 1:data.J], (1-ρ_h[j,t,last_h])^2 * primal_model[:R][j] - primal_model[:ρ][j] >= -(1+μ)*ρ_h[j,t,last_h]^2)
    end
    return primal_model
end

function add_cuts(primal_model, sp_ρval, ρ_h, data, t, μ)
    for j in 1:data.J
        if !(sp_ρval[j] in ρ_h[j,t,:])
            @constraint(primal_model, (1-sp_ρval[j])^2 * primal_model[:R][j] - primal_model[:ρ][j] >= -(1+μ)*sp_ρval[j]^2)
        end
    end
    return primal_model
end


# function benders_sp_primal(y, data, ρ_h, t)#y_fixed
#     I = 1:data.I
#     J = 1:data.J
#     λ = data.a
#     C = data.C
#     Q = data.Q
#     cv = data.cv
#     D = data.D    
#     K = 1:data.k

    
#     # H = 1:size(ρ_h[:,t,:],3)
#     H = 1:size(ρ_h[:,t,:],2) 
    
#     M = calc_big_M(data, ρ_h)

#     Dt = D/sum(λ[i, t] for i in I)

#     # Subproblem
#     m = Model(optimizer_with_attributes(Gurobi.Optimizer,
#                                     "OutputFlag" => 0,
#                                     "Threads" => 1,
#                                     "MIPGap" => 1e-5,
#                                     # "TimeLimit" => maxtime + 1,                                    
#                                     "LazyConstraints" => 1,
#                                     "Presolve" => 0,
#                                     )
#                                     )
                                    
#     # @variable(m, 0 <= y[J,K] <= 1)
#     # @variable(m, 0 <= x[I,J] <= 1)
#     # @variable(m, 0 <= z[J,K] <= 1)
#     # @variable(m, 0 <= ρ[J] <= 1)
#     @variable(m, 0 <= x[I,J])
#     @variable(m, 0 <= z[J,K])
#     @variable(m, 0 <= ρ[J])

#     @variable(m, 0 <= w[J,K])
#     @variable(m, 0 <= R[J])
    
#     @objective(m, Min, sum(C[i,j,t]*x[i,j] for i in I for j in J) +
#     0.5*sum(Dt*(R[j] + ρ[j] + sum(cv^2*(w[j,k]-z[j,k]) for k in K)) for j in J))

#     # Capacity cannot be exceeded and steady state has to be conserved
#     c1_sp = @constraint(m, [j in J], -sum(λ[i,t]*x[i,j] for i in I) >= -sum(Q[j,k]*y[j,k] for k in K), base_name = "c1_sp")
    
#     # All customer zones need to be assigned to exactly one facility
#     # c2_1_sp = @constraint(m, [i in I], sum(x[i,j] for j in J) >= 1, base_name = "c2_1_sp")
#     # c2_2_sp = @constraint(m, [i in I], -sum(x[i,j] for j in J) >= -1, base_name = "c2_2_sp")
    
#     c2_sp = @constraint(m, [i in I], sum(x[i,j] for j in J) == 1, base_name = "c2_2_sp")
    
#     @constraint(m, [j in J], sum(λ[i,t]*x[i,j] for i in I) - sum(Q[j,k]*z[j,k] for k in K) == 0, base_name = "c3_sp")

#     c4_sp = @constraint(m, [j in J, k in K], -z[j,k] >= -y[j,k], base_name = "c4_sp")

#     @constraint(m, [j in J], sum(z[j,k] for k in K) - ρ[j] == 0, base_name = "c5_sp")

#     c6_sp = @constraint(m, [j in J, k in K], -w[j,k] >= -M[j,t]*y[j,k], base_name = "c6_sp")

#     @constraint(m, [j in J], sum(w[j,k] for k in K) - R[j] == 0, base_name = "c7_sp")

#     c8_sp = @constraint(m, [j in J, h in H], (1-ρ_h[j,t,h])^2 * R[j] - ρ[j] >= -ρ_h[j,t,h]^2, base_name = "c8_sp") #; ρ_h[j,t,h] > 0

#     # c9_sp = @constraint(m, [j in J, k in K], y[j,k] == y_fixed[j,k], base_name = "c9_sp")
    
#     c10_sp = @constraint(m, [i in I, j in J], -x[i,j] >= -1, base_name = "c10_sp")
#     c11_sp = @constraint(m, [j in J, k in K], -z[j,k] >= -1, base_name = "c11_sp")
#     c12_sp = @constraint(m, [j in J], -ρ[j] >= -(1-10e-5), base_name = "c12_sp")

#     optimize!(m)
#     end_stat = termination_status(m)
#     ρval = value.(ρ)
#     xval = value.(x)
#     Rval = value.(R)
#     wval = value.(w)
#     zval = value.(z)
        
#     if end_stat == MOI.DUAL_INFEASIBLE || end_stat == MOI.INFEASIBLE
#         # Feasibility cut
#         return end_stat, xval, ρval, Rval, wval, zval, -1, dual.(c1_sp), dual.(c2_sp),dual.(c4_sp), dual.(c6_sp), dual.(c8_sp), dual.(c10_sp), dual.(c11_sp), dual.(c12_sp)
#     elseif end_stat == MOI.OPTIMAL
#         return end_stat, xval, ρval, Rval, wval, zval, objective_value(m), dual.(c1_sp), dual.(c2_sp), dual.(c4_sp), dual.(c6_sp), dual.(c8_sp), dual.(c10_sp), dual.(c11_sp), dual.(c12_sp)
#     else
#         println("other results $end_stat")
#         return end_stat, [], [], [], [], [], 0, 0, 0, 0, 0, 0, 0, 0, 0        
#     end       
# end


function ini_benders_sp_primal_comp(y_fixed, data, ρ_h, t, μ)
    I = 1:data.I
    J = 1:data.J
    λ = data.a
    C = data.C
    Q = data.Q
    cv = data.cv
    D = data.D    
    K = 1:data.k

    H = 1:size(ρ_h[:,t,:],2) 
    M = calc_big_M(data, ρ_h)
    Dt = D/sum(λ[i, t] for i in I)

    # Subproblem
    m = Model(optimizer_with_attributes(Gurobi.Optimizer,
                                    "OutputFlag" => 0,
                                    "Threads" => 1,
                                    # "MIPGap" => 1e-5,
                                    # "TimeLimit" => maxtime + 1,                                    
                                    # "LazyConstraints" => 1,
                                    "Presolve" => 0,
                                    )
                                    )
                                    
    @variable(m, 0 <= y[J,K] <= 1)
    @variable(m, 0 <= x[I,J])
    @variable(m, 0 <= z[J,K])
    @variable(m, 0 <= ρ[J])

    @variable(m, 0 <= w[J,K])
    @variable(m, 0 <= R[J])

    @expression(m, of_sp_allo, sum(C[i,j,t]*x[i,j] for i in I for j in J))
    @expression(m, of_sp_cong, 0.5*sum(Dt*(R[j] + ρ[j] + sum(cv^2*(w[j,k]-z[j,k]) for k in K)) for j in J))

    @objective(m, Min, of_sp_allo + of_sp_cong)
    
    # @objective(m, Min, sum(C[i,j,t]*x[i,j] for i in I for j in J) +
    # 0.5*sum(Dt*(R[j] + ρ[j] + sum(cv^2*(w[j,k]-z[j,k]) for k in K)) for j in J))

    # Capacity cannot be exceeded and steady state has to be conserved
    # @constraint(m, c1_sp[j in J], -sum(λ[i,t]*x[i,j] for i in I) >= -sum(Q[j,k]*y[j,k] for k in K), base_name = "c1_sp")
    @constraint(m, c1_sp[j in J], -sum(λ[i,t]*x[i,j] for i in I) >= -sum(Q[j,k]*y[j,k] for k in K) - μ*sum(Q[j,k] for k in K), base_name = "c1_sp")
    
    # All customer zones need to be assigned to exactly one facility    
    @constraint(m, c2_sp[i in I], sum(x[i,j] for j in J) == 1, base_name = "c2_2_sp")
    # @constraint(m, [i in I], sum(x[i,j] for j in J) >= 1 + μ, base_name = "c2_1_sp")
    # @constraint(m, [i in I], -sum(x[i,j] for j in J) >= -1 - μ, base_name = "c2_2_sp")
    
    @constraint(m, [j in J], sum(λ[i,t]*x[i,j] for i in I) - sum(Q[j,k]*z[j,k] for k in K) == 0, base_name = "c3_sp")

    @constraint(m, c4_sp[j in J, k in K], -z[j,k] >= - y[j,k] - μ  , base_name = "c4_sp") #- μ

    @constraint(m, [j in J], sum(z[j,k] for k in K) - ρ[j] == 0, base_name = "c5_sp")
     
    @constraint(m, c6_sp[j in J, k in K], -w[j,k] >= -M[j,t]*y[j,k] - μ*M[j,t] , base_name = "c6_sp") #- M[j,t]*μ

    @constraint(m, [j in J], sum(w[j,k] for k in K) - R[j] == 0, base_name = "c7_sp")

    @constraint(m, c8_sp[j in J, h in H], (1-ρ_h[j,t,h])^2 * R[j] - ρ[j] >= -(1+μ)ρ_h[j,t,h]^2, base_name = "c8_sp")

    @constraint(m, c9_sp[j in J, k in K], y[j,k] == y_fixed[j,k])
    
    @constraint(m, c10_sp[i in I, j in J], -x[i,j] >= -1 - μ , base_name = "c10_sp")
    @constraint(m, c11_sp[j in J, k in K], -z[j,k] >= -1 - μ , base_name = "c11_sp")
    @constraint(m, c12_sp[j in J], -ρ[j] >= -(1-10e-5) - μ , base_name = "c12_sp")

    return m
end

function benders_sp_primal_comp(m)

    optimize!(m)
    end_stat = termination_status(m)
    Cterm = value(m[:of_sp_allo])
    Congterm = value(m[:of_sp_cong])
    ρval = value.(m[:ρ])
    xval = value.(m[:x])
    Rval = value.(m[:R])
    wval = value.(m[:w])
    zval = value.(m[:z])
    π9val = dual.(m[:c9_sp])
        
    if end_stat == MOI.DUAL_INFEASIBLE || end_stat == MOI.INFEASIBLE
        # Feasibility cut
        # return end_stat, xval, ρval, Rval, wval, zval, -1, dual.(m[:c1_sp]), dual.(m[:c2_sp]),dual.(m[:c4_sp]), dual.(m[:c6_sp]), dual.(m[:c8_sp]), dual.(m[:c10_sp]), dual.(m[:c11_sp]), dual.(m[:c12_sp])
        return end_stat, xval, ρval, Rval, wval, zval, -1, -1, -1,  π9val
    elseif end_stat == MOI.OPTIMAL
        return end_stat, xval, ρval, Rval, wval, zval, objective_value(m), Cterm, Congterm, π9val
        # return end_stat, xval, ρval, Rval, wval, zval, objective_value(m), dual.(m[:c1_sp]), dual.(m[:c2_sp]), dual.(m[:c4_sp]), dual.(m[:c6_sp]), dual.(m[:c8_sp]), dual.(m[:c10_sp]), dual.(m[:c11_sp]), dual.(m[:c12_sp])
    else
        # return end_stat, [], [], [], [], [], 0, 0, 0, 0, 0, 0, 0, 0, 0
        return end_stat, [], [], [], [], [], -1, -1, -1, 0
    end       
end

function ini_benders_sp_primal(y_fixed, data, ρ_h, t, μ, params, status)
    I = 1:data.I
    J = 1:data.J
    λ = data.a
    C = data.C
    Q = data.Q
    cv = data.cv
    D = data.D    
    K = 1:data.k

    H = 1:size(ρ_h[:,t,:],2) 
    M = calc_big_M(data, ρ_h)
    Dt = D/sum(λ[i, t] for i in I)

    maxtime = max(1, params.max_time - elapsed(status))

    # Subproblem
    m = Model(optimizer_with_attributes(Gurobi.Optimizer,
                                    "OutputFlag" => 0,
                                    "Threads" => 1,
                                    "MIPGap" => 1e-5,
                                    "TimeLimit" => maxtime + 1,                                    
                                    # "LazyConstraints" => 1,
                                    "Presolve" => 0,
                                    )
                                    )
                                    
    @variable(m, 0 <= x[I,J])
    @variable(m, 0 <= z[J,K])
    @variable(m, 0 <= ρ[J])

    @variable(m, 0 <= w[J,K])
    @variable(m, 0 <= R[J])

    @expression(m, of_sp_allo, sum(C[i,j,t]*x[i,j] for i in I for j in J))
    @expression(m, of_sp_cong, 0.5*sum(Dt*(R[j] + ρ[j] + sum(cv^2*(w[j,k]-z[j,k]) for k in K)) for j in J))

    @objective(m, Min, of_sp_allo + of_sp_cong)

    # Capacity cannot be exceeded and steady state has to be conserved
    @constraint(m, c1_sp[j in J], -sum(λ[i,t]*x[i,j] for i in I) >= -sum(Q[j,k]*y_fixed[j,k] for k in K) - μ*sum(Q[j,k] for k in K), base_name = "c1_sp")
    
    # All customer zones need to be assigned to exactly one facility    
    @constraint(m, c2_sp[i in I], sum(x[i,j] for j in J) == 1, base_name = "c2_2_sp")
    
    @constraint(m, [j in J], sum(λ[i,t]*x[i,j] for i in I) - sum(Q[j,k]*z[j,k] for k in K) == 0, base_name = "c3_sp")

    @constraint(m, c4_sp[j in J, k in K], -z[j,k] >= - y_fixed[j,k] - μ  , base_name = "c4_sp")

    @constraint(m, [j in J], sum(z[j,k] for k in K) - ρ[j] == 0, base_name = "c5_sp")
     
    @constraint(m, c6_sp[j in J, k in K], -w[j,k] >= -M[j,t]*y_fixed[j,k] - μ*M[j,t] , base_name = "c6_sp")

    @constraint(m, [j in J], sum(w[j,k] for k in K) - R[j] == 0, base_name = "c7_sp")

    @constraint(m, c8_sp[j in J, h in H], (1-ρ_h[j,t,h])^2 * R[j] - ρ[j] >= -(1+μ)ρ_h[j,t,h]^2, base_name = "c8_sp")
    
    @constraint(m, c10_sp[i in I, j in J], -x[i,j] >= -1 - μ , base_name = "c10_sp")
    @constraint(m, c11_sp[j in J, k in K], -z[j,k] >= -1 - μ , base_name = "c11_sp")
    @constraint(m, c12_sp[j in J], -ρ[j] >= -(1-10e-5) - μ , base_name = "c12_sp")

    return m
end

function benders_sp_primal(m, data, ρ_h, t)

    optimize!(m)
    end_stat = termination_status(m)
    Cterm = value(m[:of_sp_allo])
    Congterm = value(m[:of_sp_cong])
    ρval = value.(m[:ρ])
    xval = value.(m[:x])
    Rval = value.(m[:R])
    wval = value.(m[:w])
    zval = value.(m[:z])
    π1val = dual.(m[:c1_sp])
    π2val = dual.(m[:c2_sp])
    π4val = dual.(m[:c4_sp])
    π6val = dual.(m[:c6_sp])
    # π8val = dual.(m[:c8_sp])
    π10val = dual.(m[:c10_sp])
    π11val = dual.(m[:c11_sp])
    π12val = dual.(m[:c12_sp])

    π8val = Array{Float64}(undef, data.J, size(ρ_h[:,t,:],2))
    π8val[:,1:9] = dual.(m[:c8_sp])
    
    if size(ρ_h[:,t,:],2) > 9
        for h in 10:size(ρ_h[:,t,:],2)
            π8val[:,h] = dual.(m[Symbol("c8_$h")])
        end
    end
        
    if end_stat == MOI.DUAL_INFEASIBLE || end_stat == MOI.INFEASIBLE
        # Feasibility cut
        return end_stat, xval, ρval, Rval, wval, zval, 10e8, 10e8, 10e8,  π1val, π2val, π4val, π6val, π8val, π10val, π11val, π12val
    elseif end_stat == MOI.OPTIMAL
        return end_stat, xval, ρval, Rval, wval, zval, objective_value(m), Cterm, Congterm,  π1val, π2val, π4val, π6val, π8val, π10val, π11val, π12val
    else
        return end_stat, [], [], [], [], [], 10e8, 10e8, 10e8, [], [], [], [], [], [], [], []
    end       
end


function benders_sp_dual(y, data, ρ_h, t, y_mp, α_mp, sp1_of, types=["Gral"], w_fc = 0, w0 = 1)
    I = 1:data.I
    J = 1:data.J
    λ = data.a
    C = data.C
    Q = data.Q
    cv = data.cv
    D = data.D    
    K = 1:data.k
    
    H = 1:size(ρ_h[:,t,:],2)    
    M = calc_big_M(data, ρ_h)   
    Dt = D/sum(λ[i, t] for i in I)

    # Subproblem
    m = Model(optimizer_with_attributes(Gurobi.Optimizer,
                                    "OutputFlag" => 0,
                                    "Threads" => 1,
                                    "MIPGap" => 1e-5,                                   
                                    "LazyConstraints" => 1,
                                    "Presolve" => 0,
                                    )
                                    )
    
    @variable(m, 0 <= π0)                                    
    @variable(m, 0 <= π1[J])
    @variable(m, π2[I])
    @variable(m, π3[J])
    @variable(m, 0 <= π4[J,K])
    @variable(m, π5[J])
    @variable(m, 0 <= π6[J,K])
    @variable(m, π7[J])
    @variable(m, 0 <= π8[J,H])
    @variable(m, 0 <= π10[I,J])
    @variable(m, 0 <= π11[J,K])
    @variable(m, 0 <= π12[J])
    
    # Original  ###############

    # @objective(m, Max, sum(-π1[j]*Q[j,k]*y[j,k] for j in J for k in K)+
    #     sum(π2[i] for i in I)+
    #     sum(-π4[j,k]*y[j,k] for j in J for k in K)+
    #     sum(-π6[j,k]*M[j,t]*y[j,k] for j in J for k in K)+
    #     sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)+
    #     sum(-π10[i,j] for i in I for j in J)+
    #     sum(-π11[j,k] for j in J for k in K)+
    #     sum(-π12[j] for j in J))

    # c1_sp_dual = @constraint(m, [i in I, j in J], - π1[j]*λ[i,t] + π2[i] + π3[j]*λ[i,t] - π10[i,j] <= C[i,j,t])
    # c2_sp_dual = @constraint(m, [j in J, h in H], -π7[j] + π8[j,h]*(1-ρ_h[j,t,h])^2 <= 0.5*Dt)
    # c3_sp_dual = @constraint(m, [j in J, h in H], -π5[j] - π8[j,h] - π12[j] <= 0.5*Dt)
    # c4_sp_dual = @constraint(m, [j in J, k in K], -π6[j,k] + π7[j] <= 0.5*Dt*cv^2)
    # c5_sp_dual = @constraint(m, [j in J, k in K], -π3[j]*Q[j,k] - π4[j,k] + π5[j] - π11[j,k] <= -0.5*Dt*cv^2)

    
    # @constraint(m, sum(-π1[j]*Q[j,k]*y_mp[j,k] for j in J for k in K)+
    # sum(π2[i] for i in I)+
    # sum(-π4[j,k]*y_mp[j,k] for j in J for k in K)+
    # sum(-π6[j,k]*M[j,t]*y_mp[j,k] for j in J for k in K)+
    # sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)+
    # sum(-π10[i,j] for i in I for j in J)+
    # sum(-π11[j,k] for j in J for k in K)+
    # sum(-π12[j] for j in J) == sp1_of)


    # General ##################

    w1  = [w_fc for j in J]
    w2 = [w_fc for i in I]
    w4 = [w_fc for j in J, k in K]
    w6 = [w_fc for j in J, k in K]
    w8 = [w_fc for j in J, h in H]
    w10 = [w_fc for i in I, j in J]
    w11 = [w_fc for j in J, k in K]
    w12 = [w_fc for j in J]

    @objective(m, Max, sum(-π1[j]*Q[j,k]*y[j,k] for j in J for k in K)+
    sum(π2[i] for i in I)+
    sum(-π4[j,k]*y[j,k] for j in J for k in K)+
    sum(-π6[j,k]*M[j,t]*y[j,k] for j in J for k in K)+
    sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)+
    sum(-π10[i,j] for i in I for j in J)+
    sum(-π11[j,k] for j in J for k in K)+
    sum(-π12[j] for j in J) - π0*α_mp[t])

    @constraint(m, [i in I, j in J], - π1[j]*λ[i,t] + π2[i] + π3[j]*λ[i,t] - π10[i,j] <= π0*C[i,j,t])
    @constraint(m, [j in J, h in H], -π7[j] + π8[j,h]*(1-ρ_h[j,t,h])^2 <= π0*0.5*Dt)
    @constraint(m, [j in J, h in H], -π5[j] - π8[j,h] - π12[j] <= π0*0.5*Dt)
    @constraint(m, [j in J, k in K], -π6[j,k] + π7[j] <= π0*0.5*Dt*cv^2)#[j,k])
    @constraint(m, [j in J, k in K], -π3[j]*Q[j,k] - π4[j,k] + π5[j] - π11[j,k] <= -π0*0.5*Dt*cv^2)#[j,k])   

    @constraint(m, sum(w1[j]*π1[j] + w12[j]*π12[j] for j in J) + sum(w2[i]*π2[i] for i in I) + 
    sum(w4[j,k]*π4[j,k] + w6[j,k]*π6[j,k] + w11[j,k]*π11[j,k] for j in J for k in K) +
    sum(w8[j, h]*π8[j,h] for j in J for h in H) + sum(w10[i,j]*π10[i,j] for i in I for j in J) + w0*π0 == 1)
    

    if "MW" in types # =="MW"
        @constraint(m, sum(-π1[j]*Q[j,k]*y_mp[j,k] for j in J for k in K)+
        sum(π2[i] for i in I)+
        sum(-π4[j,k]*y_mp[j,k] for j in J for k in K)+
        sum(-π6[j,k]*M[j,t]*y_mp[j,k] for j in J for k in K)+
        sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)+
        sum(-π10[i,j] for i in I for j in J)+
        sum(-π11[j,k] for j in J for k in K)+
        sum(-π12[j] for j in J) == sp1_of)
    end
    

    optimize!(m)
    end_stat = termination_status(m)
    if end_stat == MOI.INFEASIBLE
        # Feasibility cut
        return end_stat, -1, value.(π0), value.(π1), value.(π2), value.(π4), value.(π6), value.(π8), value.(π10), value.(π11), value.(π12) #value.(π0)
    elseif end_stat == MOI.OPTIMAL
        return end_stat, objective_value(m), value.(π0), value.(π1), value.(π2), value.(π4), value.(π6), value.(π8), value.(π10), value.(π11), value.(π12)
    else
        println("other results $end_stat")
        return end_stat, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0   
    end       
end