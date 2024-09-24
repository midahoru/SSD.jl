function model_benders(data, params, status, types=["Gral"])
    # Gurobi environment
    GRB_ENV = Gurobi.Env()

    # Magnanti-Wang
    first_it_MW = false

    # Parameters for the different acceleration techniques
    # RHS perturbation Sherali    
    if "SH" in types
        μ = 10e-6
    else
        μ = 0
    end
    # Multiplier Fischetti
    if "FC" in types
        w_fc = 1
    else
        w_fc = 0
    end

    # Problem parameters
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

    # Initialize the Master Problem with the y var (location) relaxed
    mp = ini_benders_mp(data, params, status, GRB_ENV)

    # Initialize all the subproblems (one per period)
    # Includes the restricted sp and the relaxations
    primals = Dict()
    primals_k = Dict()
    primals_j = Dict()
    primals_k_j = Dict()
    
    int_y_ind = [data.k for j in J]
    int_y = gen_y(data, int_y_ind)
    for t in T
        primals[t] = ini_benders_sp_primal(int_y, data, ρ_h, t, μ, params, status, GRB_ENV)
        primals_k[t] = ini_benders_sp_primal(int_y, data, ρ_h, t, μ, params, status, GRB_ENV, ["K"])
        primals_j[t] = ini_benders_sp_primal(int_y, data, ρ_h, t, μ, params, status, GRB_ENV, ["J"])
        primals_k_j[t] = ini_benders_sp_primal(int_y, data, ρ_h, t, μ, params, status, GRB_ENV, ["K","J"])
    end 


    ##############################################
    ######### Warmstart and initial cuts #########

    all_sols=Dict()
    y_NM_vert = ini_y(data)
    of_NM, _of_term1_NM, _of_term2_NM, _of_term3_NM, y_ind_NM, _best_x_NM, y_NM_vert = heur_nelder_mead(primals, μ, data, params, status, y_NM_vert, GRB_ENV, all_sols, min(data.I,50), 2)
    int_y = gen_y(data, y_ind_NM)
    
    # Add the cut for the resulting vertex
    # TODO: Validate if using the final simplex is a good idea
    unique_index_NM = indexin(unique(y_NM_vert), y_NM_vert)
    for y_NM in y_NM_vert[unique_index_NM]
        println("NM: Adding NM cut for $y_NM")
        for t in T
            primal_sp = primals[t]
            primal_sp = update_sp_primal(primal_sp, data, gen_y(data, y_NM), M, t, μ, ρ_h)
            
            _sp_stat, _sp_xval, _sp_ρval, _sp_Rval, _sp_wval, sp_zval, _sp_val, _Allocterm, _Congterm, π1, π2, π4_vec, π6_vec, π8, π10, π11, π12 = solve_benders_sp_primal(primal_sp, data, ρ_h, t)
            π4 = π4_vec[1]
            π6 = π6_vec[1]

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
            
            mp[Symbol("opt_NM_$(t)_$(y_NM)")] = @constraint(mp, mp[:α][t] >= expr)
            status.nOptCuts += 1 
        end
    end
    # Set the initial bounds
    lb, ub = 0, of_NM #floor(sum([minimum(C[i,:,:]) for i in I]))

    ######### End warmstart #########
    #################################



    ##################################
    ######### LP relaxations #########

    # Convergence criteria
    tol = 10e-5 # To compare lb against ub
    τ = 10e-4 # Gap     
    n_bounds = 5 # Number of bounds to consider for calculating the Convergence   

    # Keep track of the bounds
    lb_iter = Dict()
    ub_iter = Dict()
    # ... and the iterations
    n_iter_lp = 0
    
    # Solve iteratively for different relaxations
    println("\n---- SOLVING RELAXATIONS ----\n")
    for (agg, prim) in zip([["K"], ["J"], ["K", "J"], []], [primals_k, primals_j, primals_k_j, primals])
        mp, lb_temp, ub_temp, lb_iter, ub_iter, n_iter_lp, yint, ub_yint = benders_iter(mp, prim, ρ_h, data, status, ub, lb_iter, ub_iter, n_iter_lp, tol, agg, μ, first_it_MW, int_y, w_fc, τ, n_bounds)
        println("With the agg $agg the bounds are lb=$lb_temp and ub=$ub_temp")
        # Update UB in case of having an integer solution
        if sum(isinteger.(yint)) == length(yint) && ub_yint < ub
            prinlnt("UB updated from $ub to $ub_yint")
            ub = ub_yint
        end
        if lb < lb_temp
            lb = lb_temp
        end
    end

    ######### End LP relaxations #########
    ######################################
    


    ######################################
    ######### Restricted problem #########

    # Set y vars as binary
    set_binary.(mp[:y])

    println("Using $y_ind_NM as warm-start with ub=$ub\n")
    set_start_value.(mp[:y], int_y)

    # The last feasible y
    yvals_opt = int_y #Array{Float64}(undef, data.J, data.k)
    stop_adding_cuts = false

    # Activate the callback function
    function lazycb(cb)  
        
        yvals = callback_value.(cb, mp[:y])
        αvals = callback_value.(cb, mp[:α])
        
        # println("Y =$(yvals)")
        # for t in T
        #     println("Alpha $t=$(αvals[t])")
        # end

        if callback_node_status(cb, mp) != MOI.CALLBACK_NODE_STATUS_INTEGER
            return
        end    

        of_mp_node = callback_value(cb, mp[:of_mp])
        α_mp_node = callback_value(cb, mp[:α_mp])
        
        lb_node = of_mp_node + α_mp_node
        println("LB node=$lb_node iter = $(status.nIter)")
        println("UB node=$ub iter = $(status.nIter)\n")
        # println("LB node=$lb_node iter = $(status.nIter)")
        # if lb_node + tol >= ub #|| stop_adding_cuts
        #     println("Return")
        #     # println("LB node=$lb_node iter = $(status.nIter)")
        #     # println("UB node=$ub iter = $(status.nIter)\n")
        #     return
        # end

        status.nIter += 1        

        #### TODO: TEST ####
        # ρ_k = Array{Float64}(undef,data.J,data.t)                    
        # M = calc_big_M(data, ρ_h)  

        # xval = Array{Float64}(undef, data.I, data.J, data.t)
        # ρval = Array{Float64}(undef, data.J, data.t)    
        # Rval = Array{Float64}(undef, data.J, data.t)    
        # wval = Array{Float64}(undef, data.J, data.k, data.t)
        # zval = Array{Float64}(undef, data.J, data.k, data.t)

        # cost_sp = 0
        
        # for t in T   
        #     ##### Classical Benders ######
        #     primal_sp = primals[t]
        #     primal_sp = update_sp_primal(primal_sp, data, yvals, M, t, μ, ρ_h)
            
        #     sp_stat, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, sp_val, Allocterm, Congterm, π1, π2, π4_vec, π6_vec, π8, π10, π11, π12 = solve_benders_sp_primal(primal_sp, data, ρ_h, t)
        #     # println("Cost sp $t=$sp_val")
        #     # println("Alpha $t=$(αvals[t])")
        #     cost_sp += sp_val
        #     π4 = π4_vec[1]
        #     π6 = π6_vec[1]

        #     xval[:,:,t] = sp_xval
        #     ρval[:,t] = sp_ρval
        #     Rval[:,t] = sp_Rval
        #     wval[:,:,t] = sp_wval
        #     zval[:,:,t] = sp_zval
            
        #     H = 1:size(ρ_h[:,t,:],2)
        #     expr=zero(AffExpr)
        #     expr += sum(-π1[j]*Q[j,k]*mp[:y][j,k] for j in J for k in K)
        #     expr += sum(-π4[j,k]*mp[:y][j,k] for j in J for k in K)
        #     expr += sum(-π6[j,k]*M[j,t]*mp[:y][j,k] for j in J for k in K)
        #     expr += sum(π2[i] for i in I)
        #     expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
        #     expr += sum(-π10[i,j] for i in I for j in J)
        #     expr += sum(-π11[j,k] for j in J for k in K)
        #     expr += sum(-π12[j] for j in J)
            
        #     if sp_stat == MOI.DUAL_INFEASIBLE || sp_stat == MOI.INFEASIBLE # Must include both cases
        #         status.nFeasCuts += 1
        #         feas_cut = @build_constraint(0 >= expr)  
        #         MOI.submit(mp, MOI.LazyConstraint(cb), feas_cut)
        #         println("iter= $(status.nIter) adding Benders feasibility cut for t=$t ######################")
            
        #     elseif sp_stat == MOI.OPTIMAL
        #         # println("++++++   sp val = $sp_val and α_t = $(αvals[t])")
        #         # Add an optimality cut
        #         if  sp_val + tol> αvals[t] #+ tol 
                    
        #             # All but Fischetti
        #             if !("FC" in types)
        #                 status.nOptCuts += 1
        #                 opt_cut_gen = @build_constraint(mp[:α][t] >= expr)
        #                 MOI.submit(mp, MOI.LazyConstraint(cb), opt_cut_gen) 
        #                 println("iter= $(status.nIter) adding Benders optimality cut for t=$t --------------------") 
        #             end
                    
        #             ### Magnanti-Wong ###
        #             if "MW" in types # =="MW"
        #                 τ_MG = 0.7
        #                 if first_it_MW
        #                     int_y = τ_MG.*int_y .+ (1-τ_MG).*yvals
        #                 else
        #                     int_y = yvals
        #                     first_it_MW = true
        #                 end
                        
        #                 sp_stat, sp_valMW, _π0,  π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_dual(int_y, data, ρ_h, t, αvals, GRB_ENV)
        #                 # println("π0 = $_π0")

        #                 expr=zero(AffExpr)
        #                 expr += sum(-π1[j]*Q[j,k]*y[j,k] for j in J for k in K)
        #                 expr += sum(-π4[j,k]*y[j,k] for j in J for k in K)
        #                 expr += sum(-π6[j,k]*M[j,t]*y[j,k] for j in J for k in K)
                        
        #                 expr += sum(π2[i] for i in I)
        #                 expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
        #                 expr += sum(-π10[i,j] for i in I for j in J)
        #                 expr += sum(-π11[j,k] for j in J for k in K)
        #                 expr += sum(-π12[j] for j in J)

        #                 status.nOptCuts += 1
        #                 cons = @build_constraint(α[t] >=  expr)
        #                 MOI.submit(mp, MOI.LazyConstraint(cb), cons) 
        #                 println("iter= $(status.nIter) adding Magnanti-Wong optimality cut for t=$t --------------------")
        #             end
                    
        #             ### Fischetti ###
        #             if "FC" in types # == "FC"
        #                 sp_stat, sp_valMW, π0, π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_dual(yvals, data, ρ_h, t, αvals, GRB_ENV, w_fc, 1)

        #                 expr=zero(AffExpr)
        #                 expr += sum(-π1[j]*Q[j,k]*y[j,k] for j in J for k in K)
        #                 expr += sum(-π4[j,k]*y[j,k] for j in J for k in K)
        #                 expr += sum(-π6[j,k]*M[j,t]*y[j,k] for j in J for k in K)
                        
        #                 expr += sum(π2[i] for i in I)
        #                 expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
        #                 expr += sum(-π10[i,j] for i in I for j in J)
        #                 expr += sum(-π11[j,k] for j in J for k in K)
        #                 expr += sum(-π12[j] for j in J)

        #                 status.nOptCuts += 1
        #                 cons = @build_constraint(π0*α[t] >=  expr) #sp_valMW +expr
        #                 MOI.submit(mp, MOI.LazyConstraint(cb), cons) 
        #                 println("iter= $(status.nIter) adding Fischetti optimality cut for t=$t LLLLLLLLLLLLLAAAAAAAAAAAAAAAAAAA")
        #             end
        #         else
        #             println("----- Result is opt for $t")
        #         end
        #     end
        #     ρ_k[:,t] = sp_ρval
        #     # primals[t] = primal_sp
            
        #     ### Papadakos ###
        #     if "PK" in types
        #         τ_PK = 0.5            
        #         int_y = τ_PK.*int_y .+ (1-τ_PK).*yvals   
                
        #         _sp_stat, _sp_val, _π0,  π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_dual(int_y, data, ρ_h, t, αvals, GRB_ENV)
        #         # Add cuts for the linear approximation
        #         # primal_sp = add_cuts(primal_sp, sp_ρval, ρ_h, data, t, μ)                   
                
        #         status.nOptCuts += 1
        #         # cons = @build_constraint(α[t] >= sp_val + sum(π9[j,k]*(y[j,k]-int_y[j,k]) for j in J for k in K))
        #         expr=zero(AffExpr)
        #         expr += sum(-π1[j]*Q[j,k]*mp[:y][j,k] for j in J for k in K)
        #         expr += sum(-π4[j,k]*mp[:y][j,k] for j in J for k in K)
        #         expr += sum(-π6[j,k]*M[j,t]*mp[:y][j,k] for j in J for k in K)
                
        #         expr += sum(π2[i] for i in I)
        #         expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
        #         expr += sum(-π10[i,j] for i in I for j in J)
        #         expr += sum(-π11[j,k] for j in J for k in K)
        #         expr += sum(-π12[j] for j in J)
        #         cons = @build_constraint(mp[:α][t] >= expr)
        #         MOI.submit(mp, MOI.LazyConstraint(cb), cons) 
        #         println("iter= $(status.nIter) adding Papadakos optimality cut for t=$t --------------------") 
        #     end

        # end
        #### TODO: END TEST #### 

        mp, primals, ρ_h, cost_sp, xvals, ρvals, Rvals, wvals, zvals, first_it_MW, int_y = separate_cuts(mp, yvals, αvals, primals, ρ_h, data, status, types, tol, true, cb, [], μ, first_it_MW, int_y, w_fc)

        # println("LB node=$lb_node iter = $(status.nIter)")
        # println("Y node = $(yvals.data)")
        println("Sum alpha = $(sum(αvals))")
        println("Subproblem cost =$cost_sp")        
        println("Y cost callback=$(callback_value(cb, mp[:of_mp]))")
        println("Y cost calc=$(dot(F, yvals))")
        println("Node feasible cost =$(callback_value(cb, mp[:of_mp]) + cost_sp)")
        
        ub_temp = of_mp_node + cost_sp
        if ub > ub_temp 
            println("UB updated to $ub_temp with y = $(round.(yvals.data))")
            ub = ub_temp               
            yvals_opt = yvals                
        end
        println("UB=$ub\n")
        # lb_iter[status.nIter] = lb_node
        # ub_iter[status.nIter] = ub
        # println("LB= $lb_node    UB=$ub")
        # println("Y node = $(yvals.data)\n")
        # ρ_h = cat(ρ_h, ρ_k, dims=3)

        # # Test Convergence
        # if sum(αvals)  > cost_sp + tol #abs(sum(αvals)-cost_sp) <= 10e-4
        #     # stop_adding_cuts = true
        #     yvals_opt = yvals
        #     println("Y = $yvals_opt ")
        # end
        # else        
        #     # Update bounds
        #     ub_temp = dot(F, yvals) + cost_sp #callback_value.(cb, of_mp) + cost_sp
        #     if ub > ub_temp 
        #         println("UB updated to $ub_temp with y = $(yvals.data)")
        #         ub = ub_temp               
        #         # yvals_opt = yvals                
        #     end
        #     lb_iter[status.nIter] = lb_node
        #     ub_iter[status.nIter] = ub
        #     # if status.nIter > 1
        #     #     lb_iter[status.nIter] = maximum([lb, lb_iter[status.nIter-1]])
        #     # else
        #     #     lb_iter[status.nIter] = lb
        #     # end
        #     println("LB= $lb_node    UB=$ub\n")                
        #     ρ_h = cat(ρ_h, ρ_k, dims=3)
        # end
        # yvals_opt = yvals         
    end
    set_attribute(mp, MOI.LazyConstraintCallback(), lazycb)

    println("\n---- STARTING BRANCHING ----")
    end_stat_mp, _of_mp, yval_mp, _αval_mp = solve_benders_mp(mp)
    println("Y end = $(yval_mp.data)")
    println("OF end = $_of_mp")
    println("End status is $end_stat_mp")
    if end_stat_mp == MOI.INFEASIBLE
        status.endStatus = :infeasible
        return 0, 0
    elseif end_stat_mp == MOI.OPTIMAL
        status.endStatus = :optimal
    elseif end_stat_mp == MOI.TIME_LIMIT
        status.endStatus = :tlim
    end
    
    n_nodes =  MOI.get(mp, MOI.NodeCount())
    yval= 0
    try
        yval = yvals_opt.data
    catch
        yval = yvals_opt
    end
    
    println("Y opt = $(yval)")
    # yval = yval.data
    xval = Array{Float64}(undef, data.I, data.J, data.t)
    ρval = Array{Float64}(undef, data.J, data.t)    
    Rval = Array{Float64}(undef, data.J, data.t)    
    wval = Array{Float64}(undef, data.J, data.k, data.t)
    zval = Array{Float64}(undef, data.J, data.k, data.t)
                
    M = calc_big_M(data, ρ_h)
    # println("solving last step with y= $(yvals_opt.data)")
    for t in T
        primal_sp = primals[t]
        primal_sp = update_sp_primal(primal_sp, data, yval, M, t, 0, ρ_h)
        
        _1, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, _2, _Ct, _Cont, _π1, _π2, _π4, _π6, _π8, _π10, _π11, _π12 = solve_benders_sp_primal(primal_sp, data, ρ_h, t)
        xval[:,:,t] = sp_xval
        ρval[:,t] = sp_ρval
        Rval[:,t] = sp_Rval
        wval[:,:,t] = sp_wval
        zval[:,:,t] = sp_zval
    end

    Fterm = dot(F, yval)
    Allocterm = dot(C, xval) #sum(C[i,j,t]*xval[i,j,t] for i in I for j in J for t in T)  
    Congterm = 0.5 * sum(Dt[t] * (Rval[j, t] + ρval[j, t] + sum(cv^2 * (wval[j, k, t] - zval[j, k, t]) for k in K)) for j in J for t in T)
    
    return objective_value(mp), Fterm+Allocterm+Congterm, Fterm, Allocterm, Congterm, yval, xval, n_nodes, status.nFeasCuts, status.nOptCuts, lb_iter, ub_iter
    ######### End Restricted problem #############
    ##############################################

end


function ini_benders_mp(data, params, status, GRB_ENV)
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
    mp = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV),
                                    "OutputFlag" => 0,
                                    "Threads" => 1,
                                    "MIPGap" => 1e-5,
                                    "TimeLimit" => maxtime + 1,
                                    # "LazyConstraints" => 1,
                                    "Presolve" => 0,
                                    )
                                    )                              

    
    @variable(mp, 0 <= y[J,K] <= 1)
    @variable(mp, α[T] >= 0)

    @expression(mp, of_mp, sum(F[j,k]*y[j,k] for j in J for k in K))
    @expression(mp, α_mp, sum(α[t] for t in T))

    @objective(mp, Min, of_mp + α_mp)

    # At most one capacity level can be selected per facility
    @constraint(mp, [j in J], sum(y[j,k] for k in K) <= 1, base_name = "c_mp_1")

    ## Valid inequalities
    # The capacity must be enough to cover all the demand in each period
    @constraint(mp, [t in T], sum(Q[j,k]*y[j,k] for j in J for k in K) >= sum(λ[:,t]), base_name="c_mp_2")

    # Since we are allocating all the demand, the variable costs have to be greatter than or equal to
    # the cost of allocating each node to the least expensive facility (closest one)
    @constraint(mp, [t in T], α[t] >= floor(sum([minimum(C[i,:,t]) for i in I])), base_name="c_mp_3")

    return mp
end

function solve_benders_mp(m)
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

function ini_benders_sp_primal(y_fixed, data, ρ_h, t, μ, params, status, GRB_ENV, agg=[])
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

    m = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV),
                                    "OutputFlag" => 0,
                                    "Threads" => 1,
                                    "MIPGap" => 1e-5,
                                    "TimeLimit" => maxtime + 1,
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
    @constraint(m, c2_sp[i in I], sum(x[i,j] for j in J) == 1, base_name = "c2_sp")
    
    @constraint(m, [j in J], sum(λ[i,t]*x[i,j] for i in I) - sum(Q[j,k]*z[j,k] for k in K) == 0, base_name = "c3_sp")

    @constraint(m, [j in J], sum(z[j,k] for k in K) - ρ[j] == 0, base_name = "c5_sp")

    # No relaxation
    if length(agg)==0
        @constraint(m, c4_sp[j in J, k in K], -z[j,k] >= - y_fixed[j,k] - μ  , base_name = "c4_sp")
        @constraint(m, c6_sp[j in J, k in K], -w[j,k] >= -M[j,t]*y_fixed[j,k] - μ*M[j,t] , base_name = "c6_sp")
    else
        # Aggregate on K
        if "K" in agg
            @constraint(m, c4k_sp[k in K], sum(-z[j,k] for j in J) >= sum(-y_fixed[j,k] - μ for j in J), base_name = "c4k_sp")
            @constraint(m, c6k_sp[k in K], sum(-w[j,k] for j in J) >= sum(-M[j,t]*y_fixed[j,k] - μ*M[j,t] for j in J) , base_name = "c6k_sp")
        end
        # Aggregate on J
        if "J" in agg
            @constraint(m, c4j_sp[j in J], sum(-z[j,k] for k in K) >= sum(-y_fixed[j,k] - μ for k in K), base_name = "c4j_sp")
            @constraint(m, c6j_sp[j in J], sum(-w[j,k] for k in K) >= sum(-M[j,t]*y_fixed[j,k] - μ*M[j,t] for k in K) , base_name = "c6j_sp")
        end
    end

    @constraint(m, [j in J], sum(w[j,k] for k in K) - R[j] == 0, base_name = "c7_sp")

    @constraint(m, c8_sp[j in J, h in H], (1-ρ_h[j,t,h])^2 * R[j] - ρ[j] >= -(1+μ)ρ_h[j,t,h]^2, base_name = "c8_sp")
    
    @constraint(m, c10_sp[i in I, j in J], -x[i,j] >= -1 - μ , base_name = "c10_sp")
    @constraint(m, c11_sp[j in J, k in K], -z[j,k] >= -1 - μ , base_name = "c11_sp")
    @constraint(m, c12_sp[j in J], -ρ[j] >= -(1-10e-5) - μ , base_name = "c12_sp")

    return m
end

function solve_benders_sp_primal(m, data, ρ_h, t, agg=[])

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

    # Constraints 4 and 6
    π4val = 0
    π6val = 0
    π4kval = 0
    π6kval = 0
    π4jval = 0
    π6jval = 0
    # No relaxation
    if length(agg)==0
        π4val = dual.(m[:c4_sp])
        π6val = dual.(m[:c6_sp])
    else
        # Aggregate on K
        if "K" in agg
            π4kval = dual.(m[:c4k_sp])
            π6kval = dual.(m[:c6k_sp])
        end
        # Aggregate on J
        if "J" in agg            
            π4jval = dual.(m[:c4j_sp])
            π6jval = dual.(m[:c6j_sp])
        end
    end
    
    π10val = dual.(m[:c10_sp])
    π11val = dual.(m[:c11_sp])
    π12val = dual.(m[:c12_sp])

    # Constraint 8
    π8val = Array{Float64}(undef, data.J, size(ρ_h[:,t,:],2))
    π8val[:,1:9] = dual.(m[:c8_sp])
    
    if size(ρ_h[:,t,:],2) > 9
        for h in 10:size(ρ_h[:,t,:],2)
            π8val[:,h] = dual.(m[Symbol("c8_$h")])
        end
    end
        
    if end_stat == MOI.DUAL_INFEASIBLE || end_stat == MOI.INFEASIBLE
        # Feasibility cut
        return end_stat, xval, ρval, Rval, wval, zval, 10e8, 10e8, 10e8,  π1val, π2val, [π4val, π4kval, π4jval], [π6val, π6kval, π6jval], π8val, π10val, π11val, π12val
    elseif end_stat == MOI.OPTIMAL
        # Optimality cut
        return end_stat, xval, ρval, Rval, wval, zval, objective_value(m), Cterm, Congterm,  π1val, π2val, [π4val, π4kval, π4jval], [π6val, π6kval, π6jval], π8val, π10val, π11val, π12val
    else
        # Other option (?)
        return end_stat, [], [], [], [], [], 10e8, 10e8, 10e8, [], [], [], [], [], [], [], []
    end       
end

function update_sp_primal(primal_model, data, yvals, M, t, μ, ρ_h, agg=[])
    for j in 1:data.J
        set_normalized_rhs(primal_model[:c1_sp][j], -sum(data.Q[j,k]*yvals[j,k] for k in 1:data.k) - μ*sum(data.Q[j,k] for k in 1:data.k))
        if length(agg)==0
            for k in 1:data.k            
                set_normalized_rhs(primal_model[:c4_sp][j,k], -yvals[j,k] - μ)
                set_normalized_rhs(primal_model[:c6_sp][j,k], -M[j,t]*yvals[j,k] - μ*M[j,t])
            end
        end
    end

    if "K" in agg
        for k in 1:data.k
            set_normalized_rhs(primal_model[:c4k_sp][k], sum(-yvals[j,k] - μ for j in 1:data.J))
            set_normalized_rhs(primal_model[:c6k_sp][k], sum(-M[j,t]*yvals[j,k] - μ*M[j,t] for j in 1:data.J))
        end
    end

    if "J" in agg
        for j in 1:data.J
            set_normalized_rhs(primal_model[:c4j_sp][j], sum(-yvals[j,k] - μ for k in 1:data.k))
            set_normalized_rhs(primal_model[:c6j_sp][j], sum(-M[j,t]*yvals[j,k] - μ*M[j,t] for k in 1:data.k))
        end
    end

    last_h = size(ρ_h[:,t,:],2)
    if last_h > 9        
        # println("Adding the extra 8 cuts")
        primal_model[Symbol("c8_$(last_h)")] = @constraint(primal_model, [j in 1:data.J], (1-ρ_h[j,t,last_h])^2 * primal_model[:R][j] - primal_model[:ρ][j] >= -(1+μ)*ρ_h[j,t,last_h]^2)
    end
    return primal_model
end

function benders_sp_dual(y, data, ρ_h, t, α_mp, GRB_ENV, w_fc = 0, w0 = 1)
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
    
    m = Model(optimizer_with_attributes(()->Gurobi.Optimizer(GRB_ENV),
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

    @objective(m, Max, sum(π1[j]*sum(-Q[j,k]*y[j,k] for k in K) for j in J)+
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
    @constraint(m, [j in J, k in K], -π6[j,k] + π7[j] <= π0*0.5*Dt*cv^2)
    @constraint(m, [j in J, k in K], -π3[j]*Q[j,k] - π4[j,k] + π5[j] - π11[j,k] <= -π0*0.5*Dt*cv^2)   

    @constraint(m, sum(w1[j]*π1[j] + w12[j]*π12[j] for j in J) + sum(w2[i]*π2[i] for i in I) + 
    sum(w4[j,k]*π4[j,k] + w6[j,k]*π6[j,k] + w11[j,k]*π11[j,k] for j in J for k in K) +
    sum(w8[j, h]*π8[j,h] for j in J for h in H) + sum(w10[i,j]*π10[i,j] for i in I for j in J) + w0*π0 == 1)
    
    # As in Papadakos, they use the independent method
    # if "MW" in types # =="MW"
    #     @constraint(m, sum(-π1[j]*Q[j,k]*y_mp[j,k] for j in J for k in K)+
    #     sum(π2[i] for i in I)+
    #     sum(-π4[j,k]*y_mp[j,k] for j in J for k in K)+
    #     sum(-π6[j,k]*M[j,t]*y_mp[j,k] for j in J for k in K)+
    #     sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)+
    #     sum(-π10[i,j] for i in I for j in J)+
    #     sum(-π11[j,k] for j in J for k in K)+
    #     sum(-π12[j] for j in J) == sp1_of)
    # end   
    # TODO: split initialization and solving, if is possible to update.
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

function separate_cuts(m, yvals, αvals, primals, ρ_h, data, status, types, tol, lazy=false, cb=[], agg=[], μ=0, first_it_MW=false, int_y=[], w_fc=0)
    I = 1:data.I
    J = 1:data.J
    Q = data.Q
    K = 1:data.k
    T = 1:data.t

    ρ_k = Array{Float64}(undef,data.J,data.t)                    
    M = calc_big_M(data, ρ_h) 
    
    xvals = Array{Float64}(undef, data.I, data.J, data.t)
    ρvals = Array{Float64}(undef, data.J, data.t)    
    Rvals = Array{Float64}(undef, data.J, data.t)    
    wvals = Array{Float64}(undef, data.J, data.k, data.t)
    zvals = Array{Float64}(undef, data.J, data.k, data.t)

    cost_sp = 0

    for t in T
        ##### Classical Benders ######
        primal_sp = primals[t]
        primal_sp = update_sp_primal(primal_sp, data, yvals, M, t, μ, ρ_h, agg)
        
        sp_stat, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, sp_val, Cterm, Congterm, π1, π2, π4_vec, π6_vec, π8, π10, π11, π12 = solve_benders_sp_primal(primal_sp, data, ρ_h, t, agg)
        cost_sp += sp_val

        π4 = π4_vec[1]
        π6 = π6_vec[1]
        π4k = π4_vec[2]
        π6k = π6_vec[2]
        π4j = π4_vec[3]
        π6j = π6_vec[3]

        xvals[:,:,t] = sp_xval
        ρvals[:,t] = sp_ρval
        Rvals[:,t] = sp_Rval
        wvals[:,:,t] = sp_wval
        zvals[:,:,t] = sp_zval        

        H = 1:size(ρ_h[:,t,:],2)
        expr=zero(AffExpr)
        expr += sum(-π1[j]*Q[j,k]*m[:y][j,k] for j in J for k in K)
        # expr += sum(-π4[j,k]*m[:y][j,k] for j in J for k in K)
        # expr += sum(-π6[j,k]*M[j,t]*m[:y][j,k] for j in J for k in K)

        if length(agg) == 0
            expr += sum(-π4[j,k]*m[:y][j,k] for j in J for k in K)
            expr += sum(-π6[j,k]*M[j,t]*m[:y][j,k] for j in J for k in K)
        else
            if "K" in agg
                for k in K
                    expr += sum(-π4k[k]*m[:y][j,k] for j in J) + sum(-π6k[k]*M[j,t]*m[:y][j,k] for j in J)
                end
            end
            if "J" in agg
                for j in J
                    expr += sum(-π4j[j]*m[:y][j,k] for k in K) + sum(-π6j[j]*M[j,t]*m[:y][j,k] for k in K)
                end
            end
        end       
        expr += sum(π2[i] for i in I)
        expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
        expr += sum(-π10[i,j] for i in I for j in J)
        expr += sum(-π11[j,k] for j in J for k in K)
        expr += sum(-π12[j] for j in J)

        # Add a feasibility cut
        if sp_stat == MOI.DUAL_INFEASIBLE || sp_stat == MOI.INFEASIBLE # Must include both cases
            status.nFeasCuts += 1
            if lazy
                feas_cut = @build_constraint(0 >= expr)  
                MOI.submit(m, MOI.LazyConstraint(cb), feas_cut)
            else
                feas_cut = @constraint(m, 0 >= expr) 
                m[Symbol("feas_$(t)_$(status.nIter)")] = feas_cut
            end
            println("iter= $(status.nIter) adding Benders feasibility cut for t=$t ######################")
        
        elseif sp_stat == MOI.OPTIMAL
            # Add an optimality cut
            if  sp_val + tol > αvals[t]   
                # All but Fischetti
                if !("FC" in types)
                    status.nOptCuts += 1
                    if lazy
                        opt_cut_gen = @build_constraint(m[:α][t] >= expr)
                        MOI.submit(m, MOI.LazyConstraint(cb), opt_cut_gen)
                    else
                        opt_cut_gen = @constraint(m, m[:α][t] >= expr)
                        m[Symbol("opt_$(t)_$(status.nIter)")] = opt_cut_gen
                    end
                    println("iter= $(status.nIter) adding Benders optimality cut for t=$t --------------------") 
                end

                ### Magnanti-Wong ###
                if "MW" in types
                    
                    # ### Magnanti-Wong ###  
                    if first_it_MW
                        int_y = 0.5.*int_y .+ 0.5.*yvals
                    else
                        int_y = yvals
                        first_it_MW = true
                    end            
                    sp_stat, _sp_val, _π0,  π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_dual(int_y, data, ρ_h, t, αvals, GRB_ENV) #, types)

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
                println("Solution is opt for $t")
            end
        end
        ρ_k[:,t] = sp_ρval
        # primals[t] = primal_sp

        ### Papadakos ###
        if "PK" in types
            τ_PK = 0.5            
            int_y = τ_PK.*int_y .+ (1-τ_PK).*yvals   
            
            _sp_stat, _sp_val, _π0,  π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_dual(int_y, data, ρ_h, t, αvals, GRB_ENV)
            
            status.nOptCuts += 1
            expr=zero(AffExpr)
            expr += sum(-π1[j]*Q[j,k]*m[:y][j,k] for j in J for k in K)
            expr += sum(-π4[j,k]*m[:y][j,k] for j in J for k in K)
            expr += sum(-π6[j,k]*M[j,t]*m[:y][j,k] for j in J for k in K)
            
            expr += sum(π2[i] for i in I)
            expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
            expr += sum(-π10[i,j] for i in I for j in J)
            expr += sum(-π11[j,k] for j in J for k in K)
            expr += sum(-π12[j] for j in J)

            cons_PK = @build_constraint(m[:α][t] >= expr)
            MOI.submit(m, MOI.LazyConstraint(cb), cons_PK) 
            println("iter= $(status.nIter) adding Papadakos optimality cut for t=$t --------------------") 
        end
    end
    # ρ_h = cat(ρ_h, ρ_k, dims=3)
    return m, primals, ρ_h, cost_sp, xvals, ρvals, Rvals, wvals, zvals, first_it_MW, int_y
end

function benders_iter(m, primals, ρ_h, data, status, ub, lb_iter, ub_iter,  n_iter_lp, tol, agg=[], μ=0, first_it_MW=false, int_y=[], w_fc=0, τ=10e-5, n_bounds=5)
    # Convergence test
    last_bounds = zeros(n_bounds)
    conv = false
    of_lp = 0
    ub_lp = ub
    yint = 0.2
    ub_yint = 10e8
    # Solve the problem using Benders iteratively
    while !conv
        status.nIter += 1
        n_iter_lp += 1

        end_stat_lp, of_lp, yvals_lp, αvals_lp = solve_benders_mp(m)
        if of_lp + tol >= ub_lp
            break
        end
        if end_stat_lp == MOI.INFEASIBLE
            println("The problem with aggregation $agg is infeasible ****")
            break
        end
        last_bounds = append!(last_bounds, of_lp)[2:end]
        conv = last_bounds[end]-last_bounds[begin] <= τ*last_bounds[end] ? true : false
        
        m, primals, ρ_h, cost_sp, _xvals, _ρvals, _Rvals, _wvals, _zvals, first_it_MW, int_y  = separate_cuts(m, yvals_lp, αvals_lp, primals, ρ_h, data, status, [], tol, false, [], agg, μ, first_it_MW, int_y, w_fc)
        
        if sum(αvals_lp) + tol >= cost_sp
            break
        end

        ub_lp_temp = dot(data.F, yvals_lp) + cost_sp
        if sum(isinteger.(yvals_lp)) == length(yvals_lp) && ub_lp_temp < ub_yint
            yint = yvals_lp
            prinlnt("UB updated from $ub to $ub_temp")
            ub_yint = ub_temp
        end
        
        ub_lp = ub_lp_temp < ub_lp ? ub_lp_temp : ub_lp

        lb_iter[status.nIter] = of_lp
        ub_iter[status.nIter] = ub_lp
    end
    return m, of_lp, ub_lp, lb_iter, ub_iter, n_iter_lp, yint, ub_yint
end