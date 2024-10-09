function model_benders(data, params, status, types=["Gral"])
    solver = "gurobi" #cplex or gurobi
    # Solver environment
    if solver == "gurobi"
        Solver_ENV = Gurobi.Env()
    elseif solver == "cplex" 
        Solver_ENV = CPLEX.Env()
    end

    # Parameters for the different acceleration techniques
    # Magnanti-Wang & Papadakos
    first_it_MW = true
    interior_y = gen_y(data, [data.k for j in 1:data.J])
    τ_MW = 0.5
    τ_PK = 0.5
    
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
    w0 = 1

    # Keep track of the bounds
    lb_iter = Dict()
    ub_iter = Dict()

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
    n_outter_cuts = 100
    ρ_h = ini_ρ_h(data, n_outter_cuts)  
    M = calc_big_M(data, ρ_h)
    Dt = [D/sum(λ[i, t] for i in I) for t in T]

    # Initialize the Master Problem with the y var (location) relaxed
    mp = ini_benders_mp(data, params, status, solver, Solver_ENV)

    # Initialize all the subproblems (one per period)
    # Includes the restricted sp and the relaxations
    primals = Dict()
    primals_k = Dict()
    primals_j = Dict()
    primals_k_j = Dict()

    duals = Dict()

    # Find the min cap level for wich the total cap of all facilities
    # is greatter than the max demand of all the periods
    # cap_0 = findfirst(>=(maximum(sum(λ, dims=1))), sum(Q, dims=1))[2]
    # println("Min cap $cap_0")
    
    int_y_ind = [data.k for j in J] # repeat([cap_0], data.J)  #
    int_y = gen_y(data, int_y_ind)
    for t in T
        primals[t] = ini_benders_sp_primal(int_y, data, ρ_h, t, μ, params, status, solver, Solver_ENV)
        primals_k[t] = ini_benders_sp_primal(int_y, data, ρ_h, t, μ, params, status, solver, Solver_ENV, ["K"])
        primals_j[t] = ini_benders_sp_primal(int_y, data, ρ_h, t, μ, params, status, solver, Solver_ENV, ["J"])
        primals_k_j[t] = ini_benders_sp_primal(int_y, data, ρ_h, t, μ, params, status, solver, Solver_ENV, ["K","J"])
        
        duals[t] = ini_benders_sp_dual(int_y, [0 for t in T], data, ρ_h, t, params, status, solver, Solver_ENV, w_fc, w0)
    end 


    ##############################################
    ######### Warmstart and initial cuts #########

    # all_sols=Dict()
    # y_NM_vert = ini_y(data)
    # of_NM, _of_term1_NM, _of_term2_NM, _of_term3_NM, y_ind_NM, _best_x_NM, y_NM_vert = heur_nelder_mead(primals, μ, n_outter_cuts, data, params, status, solver, y_NM_vert, Solver_ENV, all_sols, 10,2) #min(data.I,50), 2)
    # int_y = gen_y(data, y_ind_NM)
    
    # # Add the cut for the resulting vertex
    # # TODO: Validate if using the final simplex is a good idea
    # # TODO: Validate the addition of cuts: they might be cutting the sol
    # unique_index_NM = indexin(unique(y_NM_vert), y_NM_vert)
    # for y_NM in y_NM_vert[unique_index_NM]
    #     Allocost, Congcost, xval_sp, ρval_sp, Rval_sp, wval_sp, zval_sp, all_sp_stat, all_sp_feas, all_sp_vals, all_sp_duals = solve_benders_sp_primal(primals, data, params, status, solver, gen_y(data, y_NM), ρ_h, n_outter_cuts, μ)
    #     for t in T
    #         duals_sp = all_sp_duals[t]
    #         π1 = duals_sp["π1"]
    #         π2 = duals_sp["π2"]
    #         π4_vec = duals_sp["π4_vec"]
    #         π6_vec = duals_sp["π6_vec"]
    #         π8 = duals_sp["π8"]
    #         π10 = duals_sp["π10"]
    #         π11 = duals_sp["π11"]
    #         π12 = duals_sp["π12"] 
    
    #         π4 = π4_vec[1]
    #         π6 = π6_vec[1]

    #         H = 1:size(ρ_h[:,t,:],2)
    #         expr =  AffExpr(0)
    #         add_to_expression!(expr, sum(-π1[j]*sum(Q[j,k]*mp[:y][j,k] for k in K) for j in J))
    #         add_to_expression!(expr, sum(-π4[j,k]*mp[:y][j,k] for j in J for k in K))
    #         add_to_expression!(expr, sum(-π6[j,k]*M[j,t]*mp[:y][j,k] for j in J for k in K))
    #         add_to_expression!(expr, sum(π2[i] for i in I))
    #         add_to_expression!(expr, sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H))
    #         add_to_expression!(expr, sum(-π10[i,j] for i in I for j in J))
    #         add_to_expression!(expr, sum(-π11[j,k] for j in J for k in K))
    #         add_to_expression!(expr, sum(-π12[j] for j in J))
            
    #         mp[Symbol("opt_NM_$(t)_$(y_NM)")] = @constraint(mp, mp[:α][t] >= expr)
    #         status.nOptCuts += 1 
    #     end
    # end
    ######### End warmstart #########
    #################################



    ##################################
    ######### LP relaxations #########

    # Set the initial bounds
    lb, ub = floor(sum([minimum(C[i,:,:]) for i in I])), 10e5 #of_NM #
    lb_iter[status.nIter] = lb
    ub_iter[status.nIter] = ub

    # Convergence criteria
    tol = 10e-5 # To compare lb against ub
    τ = 10e-4 # Gap     
    n_bounds = 5 # Number of bounds to consider for calculating the Convergence 

    set_start_value.(mp[:y], int_y)
    
    # Solve iteratively for different relaxations
    # println("\n---- SOLVING RELAXATIONS ----\n")

    # Number of iterations in the relaxations
    n_iter_lp = 0
    # Iters per relaxation
    relax_iters = []
    for (agg, prim) in zip([["K"], ["J"], ["K", "J"], []], [primals_k, primals_j, primals_k_j, primals])        
    # for (agg, prim) in zip([[]], [primals])
        println("Starting relaxation $agg")
        lb_temp, ub_temp, lb_iter, ub_iter, n_iter_lp, iter_rel, yint, ub_yint, cuts_rel = benders_iter(mp, prim, ρ_h, data, params, status, solver, ub, lb_iter, ub_iter, n_iter_lp, tol, Solver_ENV, n_outter_cuts, agg, μ, first_it_MW, int_y, w_fc, τ, n_bounds)
        push!(relax_iters, iter_rel)
        println("With the agg $agg the bounds are lb=$lb_temp and ub=$ub_temp\n")
        # Update UB in case of having an integer solution
        if sum(isinteger.(yint)) == length(yint) && ub_yint < ub
            prinlnt("UB updated from $ub to $ub_yint")
            ub = ub_yint
        end
        if lb < lb_temp
            lb = lb_temp
        end
        # If one relaxation creates cuts, then stop
        # using the others relaxed models
        if cuts_rel > 0
            break
        end
    end

    ######### End LP relaxations #########
    ######################################
    


    ######################################
    ######### Restricted problem #########

    #println("Using $y_ind_NM as warm-start with ub=$ub\n")
    # set_start_value.(mp[:y], int_y)

    # The last feasible y
    yvals_opt = int_y #Array{Float64}(undef, data.J, data.k)
    stop_adding_cuts = false

    tested_nodes_y = []

    # Activate the callback function
    function lazyCB(cb)

        # if status.nIter > 30
        #     return
        # end 

        if callback_node_status(cb, mp) != MOI.CALLBACK_NODE_STATUS_INTEGER || stop_adding_cuts
            return
        end  

        yvals = round.(callback_value.(cb, mp[:y]))
        αvals = callback_value.(cb, mp[:α])

        status.nIter += 1
        
        of_mp_node = callback_value(cb, mp[:of_mp])
        α_mp_node = callback_value(cb, mp[:α_mp])
        
        lb_node = of_mp_node + α_mp_node
        println("LB node=$lb_node iter = $(status.nIter)")
        println("UB node=$ub iter = $(status.nIter)")

        Allocost, Congcost, _xval_sp, ρval_sp, _Rval_sp, _wval_sp, _zval_sp, all_sp_stat, all_sp_feas, all_sp_vals, all_sp_duals = solve_benders_sp_primal(primals, data, params, status, solver, yvals, ρ_h, n_outter_cuts, μ)

        ub_temp = of_mp_node + Allocost + Congcost

        # Update bounds
        if !(false in all_sp_feas)
            if ub_temp < ub || abs(ub-ub_temp)/ub < tol
                println("UB updated to $ub_temp with y = $(round.(yvals.data))")
                ub = ub_temp
                # println("UPDATED y from $(yvals_opt) to $(yvals.data)")
                yvals_opt = yvals           
            end
        end 

        if lb_node < ub
            lb = lb_node
        end       

        cuts_types = "all"
        if (ub-lb_node)/lb_node < -10e-7
            cuts_types = "feas"
        end

        all_sp_dual_vals = []
        ##### TODO Duals
        if "MW" in types || "PK" in types || "FC" in types

            # For MW, use the y values for the first iteration,
            # then calculate a new point between feasible solutions
            # If y is not feasible, then use the same last feasible solution
            if "MW" in types 
                if !first_it_MW && cuts_types == "all"
                    interior_y = τ_MW.*interior_y .+ (1-τ_MW).*yvals
                elseif first_it_MW
                    interior_y = yvals
                    first_it_MW = false
                end
            end

            # For PK, if y is not feasible for all the periods
            # then use the last interior point, otherwise, updateit
            if "PK" in types 
                if cuts_types == "all"          
                    interior_y = τ_PK.*interior_y .+ (1-τ_PK).*yvals
                end
            end

            if "FC" in types
                interior_y = yvals
            end
                      
            all_sp_dual_vals = solve_benders_sp_dual(duals, data, params, status, solver, interior_y, αvals, ρ_h, n_outter_cuts)
        end
        #####

        # interior_y,
        _cuts_gen = separate_cuts(mp, yvals, αvals, ρ_h, data, status, types, tol, Solver_ENV, true, cuts_types, all_sp_stat, all_sp_vals, all_sp_duals, all_sp_dual_vals, cb, μ, [], first_it_MW, interior_y, w_fc)
        # Limit the linear approximation to 100 cuts
        # TODO avoid adding the same cut twice
        if size(ρ_h[:,1,:],2) < 1
            ρ_h = cat(ρ_h, reduce(hcat,ρval_sp), dims=3)
        end

        lb_iter[status.nIter] = lb
        ub_iter[status.nIter] = ub
        
        println("\n")
        # y_ind_pr = []
        # for r in 1:size(yvals.data)[1]
        #     cap = findfirst(==(1), yvals.data[r,:])
        #     if cap === nothing
        #         cap = 0
        #     end
        #     push!(y_ind_pr, cap)
        # end
        # println("Y node = $(y_ind_pr)")
        # for t in T
        #     println("Alpha $t=$(αvals[t])")
        # end
        # println("Sum alpha = $(sum(αvals))")
        # println("Subproblem cost =$cost_sp")        
        # println("Y cost =$(callback_value(cb, mp[:of_mp]))")
        # println("Node feasible cost =$(dot(F, yvals) + cost_sp)\n")  
            
    end
    MOI.set(mp, MOI.LazyConstraintCallback(), lazyCB)

    # Set y vars as binary
    set_binary.(mp[:y])

    println("\n---- STARTING BRANCHING ----")
    end_stat_mp, _of_mp, yval_mp, _αval_mp = solve_benders_mp(mp)
    if end_stat_mp == MOI.INFEASIBLE
        status.endStatus = :infeasible
        return 0, 0
    elseif end_stat_mp == MOI.OPTIMAL
        status.endStatus = :optimal
    elseif end_stat_mp == MOI.TIME_LIMIT
        status.endStatus = :tlim
    end
    
    n_nodes =  MOI.get(mp, MOI.NodeCount())

    # println("Primal status = $(primal_status(mp))")
    if primal_status(mp) == MOI.NO_SOLUTION
        yval_mp = int_y        
    end
    
    yval=0   
    try
        yval = round.(yvals_opt.data)
    catch
        yval = round.(yvals_opt)
    end
    
    println("Y end = $(yval_mp)")
    println("OF end = $_of_mp")
    println("End status is $end_stat_mp")
    println("Y OPT = $(yval)")
    
    xval = Array{Float64}(undef, data.I, data.J, data.t)
    ρval = Array{Float64}(undef, data.J, data.t)    
    Rval = Array{Float64}(undef, data.J, data.t)    
    wval = Array{Float64}(undef, data.J, data.k, data.t)
    zval = Array{Float64}(undef, data.J, data.k, data.t)
                
    # M = calc_big_M(data, ρ_h)
    _Allocost, _Congcost, xval_sp, ρval_sp, Rval_sp, wval_sp, zval_sp, all_sp_stat, all_sp_feas, all_sp_vals, all_sp_duals = solve_benders_sp_primal(primals, data, params, status, solver, yval, ρ_h, n_outter_cuts, μ)
    for t in T
        xval[:,:,t] = xval_sp[t]
        ρval[:,t] = ρval_sp[t]
        Rval[:,t] = Rval_sp[t]
        wval[:,:,t] = wval_sp[t]
        zval[:,:,t] = zval_sp[t]
    end

    Fterm = sum(F[j,k]*yval[j,k] for j in J for k in K) #dot(F, yval)
    Allocterm = sum(C[i,j,t]*xval[i,j,t] for i in I for j in J for t in T) #dot(C, xval) #sum(C[i,j,t]*xval[i,j,t] for i in I for j in J for t in T)  
    Congterm = 0.5 * sum(Dt[t] * (Rval[j,t] + ρval[j,t] + sum(cv^2 * (wval[j,k,t] - zval[j,k,t]) for k in K)) for j in J for t in T)
    # println(n_iter_lp)

    test_cap, test_cap_sel, test_alloc = is_sol_feas(data, yval, xval)    
    
    return objective_value(mp), Fterm+Allocterm+Congterm, Fterm, Allocterm, Congterm, yval, xval, n_nodes, status.nFeasCuts, status.nOptCuts, lb_iter, ub_iter, test_cap, test_cap_sel, test_alloc, relax_iters
    ######### End Restricted problem #############
    ##############################################

end


function ini_benders_mp(data, params, status, solver, Solver_ENV)
    I = 1:data.I
    J = 1:data.J
    λ = data.a
    C = data.C
    F = data.F
    Q = data.Q
    K = 1:data.k
    T = 1:data.t
    
    maxtime = max(1, params.max_time - elapsed(status))

    if solver == "gurobi"
        mp = Model(optimizer_with_attributes(()->Gurobi.Optimizer(Solver_ENV),
                                    "OutputFlag" => 0,
                                    "Threads" => 1,
                                    "MIPGap" => 1e-5,
                                    "MIPFocus" => 3, # 0: balance feas sol & prove opt, 1:found feas sol, 2: prove opt, 3: bound
                                    "TimeLimit" => maxtime + 1,
                                    "LazyConstraints" => 1,
                                    "Presolve" => 0,                                    
                                    "Cuts" => 0,
                                    # "FeasibilityTol" => 1e-2
                                    )
                                    )
    elseif solver == "cplex" 
        mp = Model(optimizer_with_attributes(()->CPLEX.Optimizer(Solver_ENV),
                                "CPXPARAM_ScreenOutput" => 1,
                                "CPXPARAM_Threads" => 1,
                                "CPXPARAM_MIP_Tolerances_MIPGap" => 1e-5,
                                "CPXPARAM_TimeLimit" => maxtime + 1,
                                "CPXPARAM_Preprocessing_Presolve" => 0,                                
                                "CPXPARAM_MIP_Limits_CutsFactor" => 0,
                                )
                                ) 
        
        set_attribute(mp, "CPX_PARAM_EPINT", 10e-6)
        set_attribute(mp, "CPX_PARAM_PREIND", 0)
    end

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

function ini_benders_sp_primal(y_fixed, data, ρ_h, t, μ, params, status, solver, Solver_ENV, agg=[])
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

    if solver == "gurobi"
        m = Model(optimizer_with_attributes(()->Gurobi.Optimizer(Solver_ENV),
        "OutputFlag" => 0,
        "Threads" => 1,
        "MIPGap" => 1e-5,
        "TimeLimit" => maxtime + 1,
        "Presolve" => 0,
        "Cuts" => 0,
        # "OptimalityTol" => 1e-9
        )
        )

    elseif solver == "cplex" 
        m = Model(optimizer_with_attributes(()->CPLEX.Optimizer(Solver_ENV),
                                            "CPXPARAM_ScreenOutput" => 0,
                                            "CPXPARAM_Threads" => 1,
                                            "CPXPARAM_MIP_Tolerances_MIPGap" => 1e-5,
                                            "CPXPARAM_TimeLimit" => maxtime + 1,  
                                            "CPXPARAM_Preprocessing_Presolve" => 0,
                                            "CPXPARAM_MIP_Limits_CutsFactor" => 0,
                                            )
                                            )
        set_optimizer_attribute(m, CPLEX.PassNames(), true)
        set_attribute(m, "CPX_PARAM_EPINT", 10e-6)
        set_attribute(m, "CPX_PARAM_PREIND", 0)
    end                             
                                    
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

function solve_benders_sp_primal(primals, data, params, status, solver, yvals, ρ_h, n_outter_cuts, μ, agg=[])
    T = 1:data.t

    xval_sp = []
    ρval_sp = []
    Rval_sp = []
    wval_sp = []
    zval_sp = []
    all_sp_stat = []
    all_sp_feas = []
    all_sp_vals = []
    all_sp_duals = Dict()
    Allocost, Congcost = 0, 0

    for t in T
        update_sp_primal(primals[t], data, params, status, solver, yvals, t, n_outter_cuts, μ, ρ_h, agg)            
        
        optimize!(primals[t])
        end_stat = termination_status(primals[t])

        xval = value.(primals[t][:x])
        ρval = value.(primals[t][:ρ])        
        Rval = value.(primals[t][:R])
        wval = value.(primals[t][:w])
        zval = value.(primals[t][:z])

        π1val = dual.(primals[t][:c1_sp])
        π2val = dual.(primals[t][:c2_sp])
        # Constraints 4 and 6
        π4val = 0
        π6val = 0
        π4kval = 0
        π6kval = 0
        π4jval = 0
        π6jval = 0
        # No relaxation
        if length(agg)==0
            π4val = dual.(primals[t][:c4_sp])
            π6val = dual.(primals[t][:c6_sp])
        else
            # Aggregate on K
            if "K" in agg
                π4kval = dual.(primals[t][:c4k_sp])
                π6kval = dual.(primals[t][:c6k_sp])
            end
            # Aggregate on J
            if "J" in agg            
                π4jval = dual.(primals[t][:c4j_sp])
                π6jval = dual.(primals[t][:c6j_sp])
            end
        end    
        π10val = dual.(primals[t][:c10_sp])
        π11val = dual.(primals[t][:c11_sp])
        π12val = dual.(primals[t][:c12_sp])
        # Constraint 8
        π8val = Array{Float64}(undef, data.J, size(ρ_h[:,t,:],2))
        π8val[:,1:n_outter_cuts] = dual.(primals[t][:c8_sp])
        if size(ρ_h[:,t,:],2) > n_outter_cuts
            for h in (n_outter_cuts+1):size(ρ_h[:,t,:],2)
                π8val[:,h] = dual.(primals[t][Symbol("c8_$h")])
            end
        end

        # Store the results
        push!(xval_sp,xval)
        push!(ρval_sp,ρval)
        push!(Rval_sp,Rval)
        push!(wval_sp,wval)
        push!(zval_sp,zval)

        all_sp_duals[t] = Dict()
        all_sp_duals[t]["π1"] = π1val
        all_sp_duals[t]["π2"] = π2val
        all_sp_duals[t]["π4_vec"] = [π4val, π4kval, π4jval]
        all_sp_duals[t]["π6_vec"] = [π6val, π6kval, π6jval]
        all_sp_duals[t]["π8"] = π8val
        all_sp_duals[t]["π10"] = π10val
        all_sp_duals[t]["π11"] = π11val
        all_sp_duals[t]["π12"] = π12val

        push!(all_sp_stat, end_stat)

        if end_stat == MOI.DUAL_INFEASIBLE || end_stat == MOI.INFEASIBLE
            push!(all_sp_feas, false)
            Allocterm, Congterm = 10e8, 10e8
        elseif end_stat == MOI.OPTIMAL
            push!(all_sp_feas, true)
            Allocterm, Congterm = value(primals[t][:of_sp_allo]),  value(primals[t][:of_sp_cong])
        end
        Allocost += Allocterm
        Congcost += Congterm
        push!(all_sp_vals, Allocterm + Congterm)
    end
    return Allocost, Congcost, xval_sp, ρval_sp, Rval_sp, wval_sp, zval_sp, all_sp_stat, all_sp_feas, all_sp_vals, all_sp_duals 
end

function update_sp_primal(primal_model, data, params, status, solver, yvals, t, n_outter_cuts, μ, ρ_h, agg=[])
    maxtime = max(1, params.max_time - elapsed(status))
    if solver == "gurobi"
        set_attribute(primal_model, "TimeLimit", maxtime)
    elseif solver == "cplex"
        set_attribute(primal_model, "CPXPARAM_TimeLimit", maxtime)
    end

    M = calc_big_M(data, ρ_h)
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
    if last_h > n_outter_cuts        
        # Adding the extra c8 cuts
        primal_model[Symbol("c8_$(last_h)")] = @constraint(primal_model, [j in 1:data.J], (1-ρ_h[j,t,last_h])^2 * primal_model[:R][j] - primal_model[:ρ][j] >= -(1+μ)*ρ_h[j,t,last_h]^2)
    end
    return #primal_model
end

# function benders_sp_dual(y, data, ρ_h, t, α_mp, Solver_ENV, w_fc = 0, w0 = 1)
#     I = 1:data.I
#     J = 1:data.J
#     λ = data.a
#     C = data.C
#     Q = data.Q
#     cv = data.cv
#     D = data.D    
#     K = 1:data.k
    
#     H = 1:size(ρ_h[:,t,:],2)    
#     M = calc_big_M(data, ρ_h)   
#     Dt = D/sum(λ[i, t] for i in I)

#     # TODO for the moment, until I fix the dual initialization issue
#     solver="gurobi"

#     if solver == "gurobi"
#         m = Model(optimizer_with_attributes(()->Gurobi.Optimizer(Solver_ENV),
#         "OutputFlag" => 0,
#         "Threads" => 1,
#         "MIPGap" => 1e-5,
#         #"TimeLimit" => maxtime + 1,
#         "Presolve" => 0,
#         #"Cuts" => 0,
#         # "OptimalityTol" => 1e-9
#         )
#         )
#     elseif solver == "cplex" 
#         m = Model(optimizer_with_attributes(()->CPLEX.Optimizer(Solver_ENV),
#                                             "CPXPARAM_ScreenOutput" => 0,
#                                             "CPXPARAM_Threads" => 1,
#                                             "CPXPARAM_MIP_Tolerances_MIPGap" => 1e-5,
#                                             "CPXPARAM_TimeLimit" => maxtime + 1,  
#                                             "CPXPARAM_Preprocessing_Presolve" => 0,
#                                             "CPXPARAM_MIP_Limits_CutsFactor" => 0,
#                                             )
#                                             )
#         set_optimizer_attribute(m, CPLEX.PassNames(), true)
#         set_attribute(m, "CPX_PARAM_EPINT", 10e-6)
#         set_attribute(m, "CPX_PARAM_PREIND", 0)
#     end                  
    
#     @variable(m, 0 <= π0)                                    
#     @variable(m, 0 <= π1[J])
#     @variable(m, π2[I])
#     @variable(m, π3[J])
#     @variable(m, 0 <= π4[J,K])
#     @variable(m, π5[J])
#     @variable(m, 0 <= π6[J,K])
#     @variable(m, π7[J])
#     @variable(m, 0 <= π8[J,H])
#     @variable(m, 0 <= π10[I,J])
#     @variable(m, 0 <= π11[J,K])
#     @variable(m, 0 <= π12[J])
    
#     # Original  ###############

#     # @objective(m, Max, sum(-π1[j]*Q[j,k]*y[j,k] for j in J for k in K)+
#     #     sum(π2[i] for i in I)+
#     #     sum(-π4[j,k]*y[j,k] for j in J for k in K)+
#     #     sum(-π6[j,k]*M[j,t]*y[j,k] for j in J for k in K)+
#     #     sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)+
#     #     sum(-π10[i,j] for i in I for j in J)+
#     #     sum(-π11[j,k] for j in J for k in K)+
#     #     sum(-π12[j] for j in J))

#     # c1_sp_dual = @constraint(m, [i in I, j in J], - π1[j]*λ[i,t] + π2[i] + π3[j]*λ[i,t] - π10[i,j] <= C[i,j,t])
#     # c2_sp_dual = @constraint(m, [j in J, h in H], -π7[j] + π8[j,h]*(1-ρ_h[j,t,h])^2 <= 0.5*Dt)
#     # c3_sp_dual = @constraint(m, [j in J, h in H], -π5[j] - π8[j,h] - π12[j] <= 0.5*Dt)
#     # c4_sp_dual = @constraint(m, [j in J, k in K], -π6[j,k] + π7[j] <= 0.5*Dt*cv^2)
#     # c5_sp_dual = @constraint(m, [j in J, k in K], -π3[j]*Q[j,k] - π4[j,k] + π5[j] - π11[j,k] <= -0.5*Dt*cv^2)

    
#     # @constraint(m, sum(-π1[j]*Q[j,k]*y_mp[j,k] for j in J for k in K)+
#     # sum(π2[i] for i in I)+
#     # sum(-π4[j,k]*y_mp[j,k] for j in J for k in K)+
#     # sum(-π6[j,k]*M[j,t]*y_mp[j,k] for j in J for k in K)+
#     # sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)+
#     # sum(-π10[i,j] for i in I for j in J)+
#     # sum(-π11[j,k] for j in J for k in K)+
#     # sum(-π12[j] for j in J) == sp1_of)


#     # General ##################

#     w1  = [w_fc for j in J]
#     w2 = [w_fc for i in I]
#     w4 = [w_fc for j in J, k in K]
#     w6 = [w_fc for j in J, k in K]
#     w8 = [w_fc for j in J, h in H]
#     w10 = [w_fc for i in I, j in J]
#     w11 = [w_fc for j in J, k in K]
#     w12 = [w_fc for j in J]

#     @objective(m, Max, sum(π1[j]*sum(-Q[j,k]*y[j,k] for k in K) for j in J)+
#     sum(π2[i] for i in I)+
#     sum(-π4[j,k]*y[j,k] for j in J for k in K)+
#     sum(-π6[j,k]*M[j,t]*y[j,k] for j in J for k in K)+
#     sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)+
#     sum(-π10[i,j] for i in I for j in J)+
#     sum(-π11[j,k] for j in J for k in K)+
#     sum(-π12[j] for j in J) - π0*α_mp[t])

#     @constraint(m, [i in I, j in J], - π1[j]*λ[i,t] + π2[i] + π3[j]*λ[i,t] - π10[i,j] <= π0*C[i,j,t])
#     @constraint(m, [j in J, h in H], -π7[j] + π8[j,h]*(1-ρ_h[j,t,h])^2 <= π0*0.5*Dt)
#     @constraint(m, [j in J, h in H], -π5[j] - π8[j,h] - π12[j] <= π0*0.5*Dt)
#     @constraint(m, [j in J, k in K], -π6[j,k] + π7[j] <= π0*0.5*Dt*cv^2)
#     @constraint(m, [j in J, k in K], -π3[j]*Q[j,k] - π4[j,k] + π5[j] - π11[j,k] <= -π0*0.5*Dt*cv^2)   

#     @constraint(m, sum(w1[j]*π1[j] + w12[j]*π12[j] for j in J) + sum(w2[i]*π2[i] for i in I) + 
#     sum(w4[j,k]*π4[j,k] + w6[j,k]*π6[j,k] + w11[j,k]*π11[j,k] for j in J for k in K) +
#     sum(w8[j, h]*π8[j,h] for j in J for h in H) + sum(w10[i,j]*π10[i,j] for i in I for j in J) + w0*π0 == 1)
    
#     # As in Papadakos, they use the independent method
#     # if "MW" in types # =="MW"
#     #     @constraint(m, sum(-π1[j]*Q[j,k]*y_mp[j,k] for j in J for k in K)+
#     #     sum(π2[i] for i in I)+
#     #     sum(-π4[j,k]*y_mp[j,k] for j in J for k in K)+
#     #     sum(-π6[j,k]*M[j,t]*y_mp[j,k] for j in J for k in K)+
#     #     sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)+
#     #     sum(-π10[i,j] for i in I for j in J)+
#     #     sum(-π11[j,k] for j in J for k in K)+
#     #     sum(-π12[j] for j in J) == sp1_of)
#     # end   
#     # TODO: split initialization and solving, if is possible to update.
#     optimize!(m)
#     end_stat = termination_status(m)
#     if end_stat == MOI.INFEASIBLE
#         # Feasibility cut
#         return end_stat, -1, value.(π0), value.(π1), value.(π2), value.(π4), value.(π6), value.(π8), value.(π10), value.(π11), value.(π12) #value.(π0)
#     elseif end_stat == MOI.OPTIMAL
#         return end_stat, objective_value(m), value.(π0), value.(π1), value.(π2), value.(π4), value.(π6), value.(π8), value.(π10), value.(π11), value.(π12)
#     else
#         # println("other results $end_stat")
#         return end_stat, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0   
#     end       
# end

function ini_benders_sp_dual(yvals, αvals, data, ρ_h, t, params, status, solver, Solver_ENV, w_fc = 0, w0 = 1)
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

    if solver == "gurobi"
        m = Model(optimizer_with_attributes(()->Gurobi.Optimizer(Solver_ENV),
        "OutputFlag" => 0,
        "Threads" => 1,
        "MIPGap" => 1e-5,
        "TimeLimit" => maxtime + 1,
        "Presolve" => 0,
        #"Cuts" => 0,
        # "OptimalityTol" => 1e-9
        )
        )
    elseif solver == "cplex" 
        m = Model(optimizer_with_attributes(()->CPLEX.Optimizer(Solver_ENV),
                                            "CPXPARAM_ScreenOutput" => 0,
                                            "CPXPARAM_Threads" => 1,
                                            "CPXPARAM_MIP_Tolerances_MIPGap" => 1e-5,
                                            "CPXPARAM_TimeLimit" => maxtime + 1,  
                                            "CPXPARAM_Preprocessing_Presolve" => 0,
                                            "CPXPARAM_MIP_Limits_CutsFactor" => 0,
                                            )
                                            )
        set_optimizer_attribute(m, CPLEX.PassNames(), true)
        set_attribute(m, "CPX_PARAM_EPINT", 10e-6)
        set_attribute(m, "CPX_PARAM_PREIND", 0)
    end                  
    
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

    @objective(m, Max, sum(π1[j]*sum(-Q[j,k]*yvals[j,k] for k in K) for j in J)+
    sum(π2[i] for i in I)+
    sum(-π4[j,k]*yvals[j,k] for j in J for k in K)+
    sum(-π6[j,k]*M[j,t]*yvals[j,k] for j in J for k in K)+
    sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)+
    sum(-π10[i,j] for i in I for j in J)+
    sum(-π11[j,k] for j in J for k in K)+
    sum(-π12[j] for j in J) - π0*αvals[t])

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
    return m 
end

function solve_benders_sp_dual(duals, data, params, status, solver, yvals, αvals, ρ_h, n_outter_cuts, agg=[])    
    T = 1:data.t    
    all_sp_duals = Dict()

    for t in T
        update_sp_dual(duals[t], data, params, status, solver, yvals, αvals, t, n_outter_cuts, ρ_h, agg)            
        
        optimize!(duals[t])

        π0val = value.(duals[t][:π0])
        π1val = value.(duals[t][:π1])
        π2val = value.(duals[t][:π2])
        π3val = value.(duals[t][:π3])
        π5val = value.(duals[t][:π5])
        π7val = value.(duals[t][:π7])
        # Constraints 4 and 6
        π4val = 0
        π6val = 0
        π4kval = 0
        π6kval = 0
        π4jval = 0
        π6jval = 0
        # No relaxation
        if length(agg)==0
            π4val = value.(duals[t][:π4])
            π6val = value.(duals[t][:π6])
        else
            # Aggregate on K
            if "K" in agg
                π4kval = value.(duals[t][:π4k])
                π6kval = value.(duals[t][:π6k])
            end
            # Aggregate on J
            if "J" in agg            
                π4jval = value.(duals[t][:π4j])
                π6jval = value.(duals[t][:π6j])
            end
        end    
        π10val = value.(duals[t][:π10])
        π11val = value.(duals[t][:π11])
        π12val = value.(duals[t][:π12])
        # Constraint 8
        # π8val = Array{Float64}(undef, data.J, size(ρ_h[:,t,:],2))
        π8val = value.(duals[t][:π8])

        # Store the results
        all_sp_duals[t] = Dict()
        all_sp_duals[t]["π0"] = π0val
        all_sp_duals[t]["π1"] = π1val
        all_sp_duals[t]["π2"] = π2val
        all_sp_duals[t]["π3"] = π3val
        all_sp_duals[t]["π5"] = π5val
        all_sp_duals[t]["π7"] = π7val
        all_sp_duals[t]["π4_vec"] = [π4val, π4kval, π4jval]
        all_sp_duals[t]["π6_vec"] = [π6val, π6kval, π6jval]
        all_sp_duals[t]["π8"] = π8val
        all_sp_duals[t]["π10"] = π10val
        all_sp_duals[t]["π11"] = π11val
        all_sp_duals[t]["π12"] = π12val
    end
    return all_sp_duals
end

function update_sp_dual(dual, data, params, status, solver, yvals, αvals, t, n_outter_cuts, ρ_h, agg=[])
    maxtime = max(1, params.max_time - elapsed(status))
    if solver == "gurobi"
        set_attribute(dual, "TimeLimit", maxtime)
    elseif solver == "cplex"
        set_attribute(dual, "CPXPARAM_TimeLimit", maxtime)
    end

    M = calc_big_M(data, ρ_h)
    set_objective_coefficient(dual, dual[:π0], -αvals[t])
    for j in 1:data.J
        set_objective_coefficient(dual, dual[:π1][j], -sum(data.Q[j,k]*yvals[j,k] for k in 1:data.k))

        for k in 1:data.k            
            set_objective_coefficient(dual, dual[:π4][j,k], -yvals[j,k])
            set_objective_coefficient(dual, dual[:π6][j,k], -M[j,t]*yvals[j,k])
        end
    end
end

function separate_cuts(m, yvals, αvals, ρ_h, data, status, types, tol, Solver_ENV, lazy=false, cuts_types="feas", all_sp_stat=[], all_sp_vals=[], all_sp_duals=[], all_sp_dual_vals=[], cb=[], μ=0, agg=[], first_it_MW=false, interior_y=[], w_fc=0)
    I = 1:data.I
    J = 1:data.J
    Q = data.Q
    K = 1:data.k
    T = 1:data.t

    # ρ_k = Array{Float64}(undef,data.J,data.t)                    
    M = calc_big_M(data, ρ_h) 

    cuts_gen = false
    
    # xvals = Array{Float64}(undef, data.I, data.J, data.t)
    # ρvals = Array{Float64}(undef, data.J, data.t)    
    # Rvals = Array{Float64}(undef, data.J, data.t)    
    # wvals = Array{Float64}(undef, data.J, data.k, data.t)
    # zvals = Array{Float64}(undef, data.J, data.k, data.t)

    # cost_sp = 0
    # all_feas = true
    # opt_cuts = []

    for t in T
        ##### Classical Benders ######
        sp_stat = all_sp_stat[t]
        sp_val = all_sp_vals[t]
        duals_sp = all_sp_duals[t]

        π1 = duals_sp["π1"]
        π2 = duals_sp["π2"]
        π4_vec = duals_sp["π4_vec"]
        π6_vec = duals_sp["π6_vec"]
        π8 = duals_sp["π8"]
        π10 = duals_sp["π10"]
        π11 = duals_sp["π11"]
        π12 = duals_sp["π12"] 

        π4 = π4_vec[1]
        π6 = π6_vec[1]
        π4k = π4_vec[2]
        π6k = π6_vec[2]
        π4j = π4_vec[3]
        π6j = π6_vec[3]

        # xvals[:,:,t] = sp_xval
        # ρvals[:,t] = sp_ρval
        # Rvals[:,t] = sp_Rval
        # wvals[:,:,t] = sp_wval
        # zvals[:,:,t] = sp_zval        

        H = 1:size(ρ_h[:,t,:],2)
        expr = AffExpr(0) 
        add_to_expression!(expr, sum(π1[j]*sum(-Q[j,k]*m[:y][j,k] for k in K) for j in J))

        if length(agg) == 0
            add_to_expression!(expr, sum(-π4[j,k]*m[:y][j,k] for j in J for k in K))
            add_to_expression!(expr, sum(-π6[j,k]*M[j,t]*m[:y][j,k] for j in J for k in K))
        else
            if "K" in agg
                for k in K
                    add_to_expression!(expr, sum(-π4k[k]*m[:y][j,k] for j in J))
                    add_to_expression!(expr, sum(-π6k[k]*M[j,t]*m[:y][j,k] for j in J))
                end
            end
            if "J" in agg
                for j in J
                    add_to_expression!(expr, sum(-π4j[j]*m[:y][j,k] for k in K))
                    add_to_expression!(expr, sum(-π6j[j]*M[j,t]*m[:y][j,k] for k in K))
                end
            end
        end
        add_to_expression!(expr, sum(π2[i] for i in I))
        add_to_expression!(expr, sum(-π8[j,h]*(ρ_h[j,t,h]^2) for j in J for h in H))
        add_to_expression!(expr, sum(-π10[i,j] for i in I for j in J))
        add_to_expression!(expr, sum(-π11[j,k] for j in J for k in K))
        add_to_expression!(expr, sum(-π12[j] for j in J))
        
        # Add a feasibility cut
        if sp_stat == MOI.INFEASIBLE || sp_stat == MOI.DUAL_INFEASIBLE
            status.nFeasCuts += 1
            if lazy
                feas_cut = @build_constraint(0 >= expr)  
                MOI.submit(m, MOI.LazyConstraint(cb), feas_cut)
            else
                feas_cut = @constraint(m, 0 >= expr) 
                m[Symbol("feas_$(t)_$(status.nIter)")] = feas_cut
            end
            #println("iter= $(status.nIter) adding Benders feasibility cut for t=$t ######################")
        
        elseif sp_stat == MOI.OPTIMAL && cuts_types == "all" 
            # println("SP cost for $t = $sp_val")
            # println("Alpha cost for $t= $(αvals[t])")
            # Add an optimality cut
            if  sp_val > αvals[t] + 10e-2 #tol
            # if  (sp_val-αvals[t])/αvals[t] > 10e-2
                cuts_gen = true
                
                # All but Fischetti
                if !("FC" in types)
                # if "Gral" in types
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

                ### Magnanti-Wong or Fischetti ###
                if "MW" in types || "FC" in types
                    dual_sp_vars = all_sp_dual_vals[t]

                    π0 = dual_sp_vars["π0"]
                    π1 = dual_sp_vars["π1"]
                    π2 = dual_sp_vars["π2"]
                    π4_vec = dual_sp_vars["π4_vec"]
                    π6_vec = dual_sp_vars["π6_vec"]
                    π8 = dual_sp_vars["π8"]
                    π10 = dual_sp_vars["π10"]
                    π11 = dual_sp_vars["π11"]
                    π12 = dual_sp_vars["π12"] 

                    π4 = π4_vec[1]
                    π6 = π6_vec[1]
                    π4k = π4_vec[2]
                    π6k = π6_vec[2]
                    π4j = π4_vec[3]
                    π6j = π6_vec[3]

                    expr_dual = AffExpr(0) 
                    add_to_expression!(expr_dual, sum(π1[j]*sum(-Q[j,k]*m[:y][j,k] for k in K) for j in J))

                    if length(agg) == 0
                        add_to_expression!(expr_dual, sum(-π4[j,k]*m[:y][j,k] for j in J for k in K))
                        add_to_expression!(expr_dual, sum(-π6[j,k]*M[j,t]*m[:y][j,k] for j in J for k in K))
                    # TODO ver cómo incluir las agregaciones en el dual
                    else
                        if "K" in agg
                            for k in K
                                add_to_expression!(expr_dual, sum(-π4k[k]*m[:y][j,k] for j in J))
                                add_to_expression!(expr_dual, sum(-π6k[k]*M[j,t]*m[:y][j,k] for j in J))
                            end
                        end
                        if "J" in agg
                            for j in J
                                add_to_expression!(expr_dual, sum(-π4j[j]*m[:y][j,k] for k in K))
                                add_to_expression!(expr_dual, sum(-π6j[j]*M[j,t]*m[:y][j,k] for k in K))
                            end
                        end
                    end
                    add_to_expression!(expr_dual, sum(π2[i] for i in I))
                    add_to_expression!(expr_dual, sum(-π8[j,h]*(ρ_h[j,t,h]^2) for j in J for h in H))
                    add_to_expression!(expr_dual, sum(-π10[i,j] for i in I for j in J))
                    add_to_expression!(expr_dual, sum(-π11[j,k] for j in J for k in K))
                    add_to_expression!(expr_dual, sum(-π12[j] for j in J))

                    ### Magnanti-Wong ###
                    if "MW" in types                    
                        # τ_MW = 0.5
                        # if !first_it_MW
                        #     interior_y = τ_MW.*interior_y .+ (1-τ_MW).*yvals
                        # else
                        #     interior_y = yvals
                        #     first_it_MW = true
                        # end            
                        # _sp_stat, _sp_val, _π0,  π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_dual(interior_y, data, ρ_h, t, αvals, Solver_ENV)                    

                        if lazy
                            cons_MW = @build_constraint(m[:α][t] >= expr_dual)
                            MOI.submit(m, MOI.LazyConstraint(cb), cons_MW)

                            status.nOptCuts += 1
                            #println("iter= $(status.nIter) adding Magnanti-Wong optimality cut for t=$t --------------------")
                        end
                    end                                        
                    
                    ### Fischetti ###
                    if "FC" in types
                        # _sp_stat, _sp_val, π0, π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_dual(yvals, data, ρ_h, t, αvals, Solver_ENV, w_fc, 1)

                        # expr=zero(AffExpr)
                        # expr += sum(-π1[j]*Q[j,k]*m[:y][j,k] for j in J for k in K)
                        # expr += sum(-π4[j,k]*m[:y][j,k] for j in J for k in K)
                        # expr += sum(-π6[j,k]*M[j,t]*m[:y][j,k] for j in J for k in K)
                        
                        # expr += sum(π2[i] for i in I)
                        # expr += sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)
                        # expr += sum(-π10[i,j] for i in I for j in J)
                        # expr += sum(-π11[j,k] for j in J for k in K)
                        # expr += sum(-π12[j] for j in J)
                        
                        cons_FC = @build_constraint(π0*m[:α][t] >=  expr_dual)
                        MOI.submit(m, MOI.LazyConstraint(cb), cons_FC)
                        status.nOptCuts += 1
                        println("iter= $(status.nIter) adding Fischetti optimality cut for t=$t --------------------")
                    end
                end                
            else
                # println("Solution is opt for $t")
            end
        end
        # ρ_k[:,t] = sp_ρval

        ### Papadakos ###
        if "PK" in types
            cuts_gen = true
            # τ_PK = 0.5
            # if cuts_types == "all"          
            #     interior_y = τ_PK.*interior_y .+ (1-τ_PK).*yvals
            # end
            
            
            # _sp_stat, _sp_val, _π0,  π1, π2, π4, π6, π8, π10, π11, π12 = benders_sp_dual(interior_y, data, ρ_h, t, αvals, Solver_ENV)

            dual_sp_vars = all_sp_dual_vals[t]

            π0 = dual_sp_vars["π0"]
            π1 = dual_sp_vars["π1"]
            π2 = dual_sp_vars["π2"]
            π4_vec = dual_sp_vars["π4_vec"]
            π6_vec = dual_sp_vars["π6_vec"]
            π8 = dual_sp_vars["π8"]
            π10 = dual_sp_vars["π10"]
            π11 = dual_sp_vars["π11"]
            π12 = dual_sp_vars["π12"] 

            π4 = π4_vec[1]
            π6 = π6_vec[1]
            π4k = π4_vec[2]
            π6k = π6_vec[2]
            π4j = π4_vec[3]
            π6j = π6_vec[3]
            
            
            expr_dual = AffExpr(0) 
            add_to_expression!(expr_dual, sum(π1[j]*sum(-Q[j,k]*m[:y][j,k] for k in K) for j in J))

            if length(agg) == 0
                add_to_expression!(expr_dual, sum(-π4[j,k]*m[:y][j,k] for j in J for k in K))
                add_to_expression!(expr_dual, sum(-π6[j,k]*M[j,t]*m[:y][j,k] for j in J for k in K))
            # TODO ver cómo incluir las agregaciones en el dual
            else
                if "K" in agg
                    for k in K
                        add_to_expression!(expr_dual, sum(-π4k[k]*m[:y][j,k] for j in J))
                        add_to_expression!(expr_dual, sum(-π6k[k]*M[j,t]*m[:y][j,k] for j in J))
                    end
                end
                if "J" in agg
                    for j in J
                        add_to_expression!(expr_dual, sum(-π4j[j]*m[:y][j,k] for k in K))
                        add_to_expression!(expr_dual, sum(-π6j[j]*M[j,t]*m[:y][j,k] for k in K))
                    end
                end
            end
            add_to_expression!(expr_dual, sum(π2[i] for i in I))
            add_to_expression!(expr_dual, sum(-π8[j,h]*(ρ_h[j,t,h]^2) for j in J for h in H))
            add_to_expression!(expr_dual, sum(-π10[i,j] for i in I for j in J))
            add_to_expression!(expr_dual, sum(-π11[j,k] for j in J for k in K))
            add_to_expression!(expr_dual, sum(-π12[j] for j in J))

            cons_PK = @build_constraint(m[:α][t] >= expr_dual)
            MOI.submit(m, MOI.LazyConstraint(cb), cons_PK)
            status.nOptCuts += 1
            # println("iter= $(status.nIter) adding Papadakos optimality cut for t=$t --------------------") 
        end
    end
    # return interior_y, cuts_gen
    return cuts_gen
end

function benders_iter(m, prims, ρ_h, data, params, status, solver, ub, lb_iter, ub_iter, n_iter_lp, tol, Solver_ENV, n_outter_cuts, agg=[], μ=0, first_it_MW=false, int_y=[], w_fc=0, τ=10e-5, n_bounds=5)
    # Convergence test
    last_bounds = zeros(n_bounds)
    conv = false
    lb_lp = 0
    ub_lp = ub
    yint = 0.2
    ub_yint = 10e8
    iter_rel = 0
    cuts_rel = 0
    
    # Solve the problem using Benders iteratively
    while !conv
        status.nIter += 1
        n_iter_lp += 1
        iter_rel += 1

        end_stat_lp, of_lp, yvals_lp, αvals_lp = solve_benders_mp(m)
        if of_lp + tol >= ub_lp
            break
        end
        lb_lp = of_lp
        if end_stat_lp == MOI.INFEASIBLE
            println("The problem with aggregation $agg is infeasible ****")
            break
        end
        last_bounds = append!(last_bounds, lb_lp)[2:end]
        conv = last_bounds[end]-last_bounds[begin] <= τ*last_bounds[end] ? true : false

        Allocost, Congcost, _xval_sp, _ρval_sp, _Rval_sp, _wval_sp, _zval_sp, all_sp_stat, all_sp_feas, all_sp_vals, all_sp_duals = solve_benders_sp_primal(prims, data, params, status, solver, yvals_lp, ρ_h, n_outter_cuts, μ, agg)
        
        # _interior_y, 
        cuts_gen = separate_cuts(m, yvals_lp, αvals_lp, ρ_h, data, status, [], tol, Solver_ENV, false, "all", all_sp_stat, all_sp_vals, all_sp_duals, [], [], μ, agg)

        if cuts_gen
            cuts_rel += 1
        end
        
        ub_lp_temp = dot(data.F, yvals_lp) + Allocost + Congcost
        # Update the feasible, integer, solution
        if sum(isinteger.(yvals_lp)) == length(yvals_lp) && ub_lp_temp < ub_yint
            yint = yvals_lp
            prinlnt("UB updated from $ub to $ub_temp")
            ub_yint = ub_temp
        end
        
        # Update the UB
        ub_lp = ub_lp_temp < ub_lp ? ub_lp_temp : ub_lp
        println("New UB =$ub_lp")

        lb_iter[status.nIter] = lb_lp
        ub_iter[status.nIter] = ub_lp
        # println("Iter = $(status.nIter)")

        if sum(αvals_lp) + tol >= Allocost + Congcost
            break
        end
    end
    return lb_lp, ub_lp, lb_iter, ub_iter, n_iter_lp, iter_rel, yint, ub_yint, cuts_rel
end