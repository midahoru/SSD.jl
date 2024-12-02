function model_benders_M(data, params, status, relaxations=["Clust"], types=["B"], n_outer_cuts=32, ρ_h=[], M=[])
    solver = "gurobi" #cplex or gurobi
    # Solver environment
    if solver == "gurobi"
        Solver_ENV = Gurobi.Env()
    elseif solver == "cplex" 
        Solver_ENV = CPLEX.Env()
    end

    relaxations=["LP"]

    # Parameters for the different acceleration techniques
    # Magnanti-Wang & Papadakos
    first_it_MW = true
    interior_y = gen_y(data, repeat([data.k], data.J))
    τ_MW = 0.5
    τ_PK = 0.5
    
    # RHS perturbation Sherali    
    if "SH" in types
        μ =  0 #10e-6
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
    # and of the number of vars and constraints
    n_vars, n_cons = 0, []

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
    ρ_h = ini_ρ_h(data, n_outer_cuts)  
    M = calc_big_M(data, ρ_h)
    Dt = [D/sum(λ[i, t] for i in I) for t in T]

    # Initialize the Master Problem with the y var (location) relaxed
    mp = ini_benders_mp(data, params, status, solver, Solver_ENV)
    push!(n_cons, num_constraints(mp; count_variable_in_set_constraints = true))

    # Initialize all the subproblems (one per period)
    # Includes the restricted sp and the relaxations
    primals = Dict()
    duals = Dict()    
    dualsMW = Dict()

    # Find the min cap level for wich the total cap of all facilities
    # is greatter than the max demand of all the periods
    # cap_0 = findfirst(>=(maximum(sum(λ, dims=1))), sum(Q, dims=1))[2]
    # println("Min cap $cap_0")
    
    int_y_ind = repeat([0], data.J) #[0 for j in J]  #[data.k for j in J] # repeat([cap_0], data.J)  #
    int_y = gen_y(data, int_y_ind)
    for t in T
        if "SH" in types
            primals[t] = ini_benders_sp_primal(int_y, data, ρ_h, M, t, μ, params, status, solver, Solver_ENV)
        elseif "B" in types || "PK" in types || "FC" in types
            duals[t] = ini_benders_sp_dual(int_y, data, ρ_h, M, t, params, status, solver, Solver_ENV, types, w_fc, w0)
        elseif "MW" in types
            duals[t] = ini_benders_sp_dual(int_y, data, ρ_h, M, t, params, status, solver, Solver_ENV, ["B"], w_fc, w0)
            dualsMW[t] = ini_benders_sp_dual(int_y, data, ρ_h, M, t, params, status, solver, Solver_ENV, types, w_fc, w0)
        end
    end
    


    ##############################################
    ######### Warmstart and initial cuts #########

    # all_sols=Dict()
    # y_NM_vert = ini_y(data)
    # of_NM, _of_term1_NM, _of_term2_NM, _of_term3_NM, y_ind_NM, _best_x_NM, y_NM_vert = heur_nelder_mead(primals, μ, n_outer_cuts, data, params, status, solver, y_NM_vert, Solver_ENV, all_sols, 10,2) #min(data.I,50), 2)
    # int_y = gen_y(data, y_ind_NM)
    
    # # Add the cut for the resulting vertex
    # # TODO: Validate if using the final simplex is a good idea
    # # TODO: Validate the addition of cuts: they might be cutting the sol
    # unique_index_NM = indexin(unique(y_NM_vert), y_NM_vert)
    # for y_NM in y_NM_vert[unique_index_NM]
    #     Allocost, Congcost, xval_sp, ρval_sp, Rval_sp, wval_sp, zval_sp, all_sp_stat, all_sp_feas, all_sp_vals, all_sp_duals = solve_benders_sp_primal(primals, data, params, status, solver, gen_y(data, y_NM), ρ_h, n_outer_cuts, μ)
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
    relax_max_iters = 1e4

    set_start_value.(mp[:y], int_y)
    
    # Solve iteratively for different relaxations
    println("\n---- SOLVING RELAXATIONS ----\n")

    # Iters per relaxation
    relax_iters = []
    # Cuts per relaxation
    relax_cuts = []
    # lb per relaxation
    relax_lb = []

    # Define the relaxations to use
    rel_to_run = [] # [["K"], ["J"], ["K", "J"], []]
    #[primals_k, primals_j, primals_k_j, primals]
    max_relax = []
    if "K" in relaxations
        push!(rel_to_run, ["K"]) # = [types]
        push!(max_relax, 1)
        # sp_to_test = [primals_k]
    elseif "J" in relaxations
        push!(rel_to_run, ["J"]) # = [types]
        push!(max_relax, 1)
        # sp_to_test = [primals_j]
    elseif "KJ" in relaxations
        push!(rel_to_run, ["KJ"]) # = [["K", "J"]]
        push!(max_relax, 1)
        # sp_to_test = [primals_k_j]
    elseif "LP" in relaxations
        push!(rel_to_run, []) # = [[]]
        push!(max_relax, 1)
        # sp_to_test = [primals]
    elseif "Clust" in relaxations
        push!(max_relax, 5)
        n_clus = max(2,Int64(floor(data.J/10))) # Reduce in 1 order of magnitud
        res_clus = kmeans(data.Jcoords, n_clus)
        @assert nclusters(res_clus) == n_clus "The number of clusters is not the specified"
        # Clusters of facilities
        set_R = []
        for clus in 1:nclusters(res_clus)
            set_Rl=[]
            fac = findall(==(clus), assignments(res_clus))
            for f in fac
                for k in K
                    push!(set_Rl, (f,k))
                end
            end
            push!(set_R, set_Rl)
        end
        push!(rel_to_run, set_R)
    end

    n_iter_rel_clus = 1
    keep_relax = true
    yval_rel_sol = []
    αval_rel_sol = []
    if length(rel_to_run) > 150 
        for (agg, max_rel) in zip(rel_to_run, max_relax)
        # for agg in rel_to_run
            while n_iter_rel_clus <= max_rel && keep_relax
                # println("Starting relaxation $agg")

                prim_rel = Dict()
                for t in T
                    prim_rel[t] = ini_benders_sp_primal(int_y, data, ρ_h, M, t, 0, params, status, solver, Solver_ENV, agg)
                end
                                                                                                                                                       
                lb_temp, ub_temp, lb_iter, ub_iter, iter_rel, cuts_rel, yval_rel, αval_rel, xval_rel, ρval_rel, Rval_rel, wval_rel, zval_rel, ρ_h = benders_iter(mp, prim_rel, ["Rel"], ρ_h, M, data, params, status, solver, lb, ub, lb_iter, ub_iter, tol, n_outer_cuts, agg, 0, int_y, w_fc, τ, n_bounds, relax_max_iters)
                push!(relax_iters, iter_rel)
                push!(relax_cuts, cuts_rel)        
                push!(relax_lb, lb_temp)
                println("With this agg the bounds are lb=$lb_temp and ub=$ub_temp after $cuts_rel cuts\n")
                # Update UB in case of having an integer solution
                # TODO arreglar para considerar soluciones muy cercanas a un entero
                # if sum(isinteger.(yint)) == length(yint) && ub_yint < ub
                #     prinlnt("UB updated from $ub to $ub_yint")
                #     ub = ub_yint
                # end
                if lb < lb_temp
                    lb = lb_temp
                    yval_rel_sol = yval_rel
                    αval_rel_sol = αval_rel
                end
                
                if cuts_rel <= 5*data.t
                    keep_relax = false
                    if cuts_rel == 0
                        break
                    end
                end

                # Finds the pairs (j,k) violating the original constraints 
                if "Clust" in relaxations
                    calc_viol!(data, yval_rel, wval_rel, zval_rel, M, agg)
                end
                # If one relaxation creates cuts, then stop
                # using the others relaxed models
                # if cuts_rel > 0 
                #     break
                # end
                n_iter_rel_clus += 1
            end
        end
    end

    push!(n_cons, num_constraints(mp; count_variable_in_set_constraints = true))

    
    
    # Just to compare the results of the relaxations
    if ["K"] == types || ["J"] == types || ["KJ"] == types || ["LP"] == types || ["Clust"] == types
        return  relax_iters, relax_cuts, relax_lb
    end

    ######### End LP relaxations #########
    ######################################
    


    ######################################
    ######### Restricted problem #########

    # Update lower bound of the α from the relaxations
    # set_normalized_rhs([mp[:α_lb][t] for t in T], [αval_rel_sol[t] for t in T])
    # set_start_value.(mp[:y], ceil.(yval_rel_sol))   


    #println("Using $y_ind_NM as warm-start with ub=$ub\n")
    # set_start_value.(mp[:y], int_y)

    # First iteration in PK
    if "PK" in types
        sp_stats_PK, sp_feas_PK, sp_ofs_PK, sp_duals_PK = solve_benders_sp_dual(duals, data, params, status, solver, interior_y, repeat([0], data.t), ρ_h, M, n_outer_cuts, types)        
        for t in T
            expr = get_expr(mp, sp_duals_PK[t], data, ρ_h, M, t, types)
            if sp_stats_PK[t] == MOI.OPTIMAL
                status.nOptCuts += 1                
                opt_cut = @constraint(mp, mp[:α][t] >= expr) 
                mp[Symbol("PK_opt_cut_$(t)_$(status.nIter)")] = opt_cut
            end
        end
    end

    # The last feasible y
    yvals_opt = int_y #Array{Float64}(undef, data.J, data.k)

    update_ρ = false
    update_M = false

    # Activate the callback function
    function lazyCB(cb)

        if callback_node_status(cb, mp) != MOI.CALLBACK_NODE_STATUS_INTEGER # || update_M
            return
        end  

        yvals = round.(callback_value.(cb, mp[:y]))
        αvals = callback_value.(cb, mp[:α])

        status.nIter += 1
        
        of_mp_node = callback_value(cb, mp[:of_mp])
        α_mp_node = callback_value(cb, mp[:α_mp])
        
        lb_node = of_mp_node + α_mp_node
            
        ## Solve SPs
        if "SH" in types
            sp_vars, sp_stats, sp_feas, sp_ofs, sp_duals = solve_benders_sp_primal(primals, data, params, status, solver, yvals, ρ_h, M, n_outer_cuts, μ, update_ρ)
        elseif "B" in types || "PK" in types || "FC" in types
            sp_stats, sp_feas, sp_ofs, sp_duals = solve_benders_sp_dual(duals, data, params, status, solver, yvals, αvals, ρ_h, M, n_outer_cuts, types, update_ρ, w_fc, w0)
        elseif "MW" in types
            # This is the dual that provides feas cuts
            sp_stats, sp_feas, sp_ofs, sp_duals = solve_benders_sp_dual(duals, data, params, status, solver, yvals, αvals, ρ_h, M, n_outer_cuts, ["B"], update_ρ, w_fc, w0)
        end

        # Error getting the duals
        if -1 in sp_ofs
            return
        end
        
        ub_temp = of_mp_node + sum(sp_ofs)
        if "FC" in types
            println("summing")
            println(ub_temp)
            # println(sum([αvals[t] for t in T]))
            ub_temp += sum(αvals) #sp_duals[t]["π0"]*
            println(ub_temp)
        end        

        # Update bounds
        if !(false in sp_feas)
            if ub_temp < ub 
                ub = ub_temp
                yvals_opt = yvals                
                # println("LB node=$lb_node iter = $(status.nIter)")
                println("UB node=$ub iter = $(status.nIter)\n")
            end
        end 

        if lb_node < ub
            lb = lb_node
        end

        cuts_types = "all"
        # If the lb is greatter than the ub, then only add feasibility cuts
        if lb_node > ub + tol # || "MW" in types
            cuts_types = "feas"
            # If the MP sol is feasible, then return
            if !(false in sp_feas)
                return
            end
        end
        
        # println("Feas $sp_stats")
        # println("Big M = $M")
        # println("rhho = $ρ_h")
        # _cuts_gen, _cuts_rel, cuts_sep_other = separate_cuts(mp, αvals, ρ_h, M, data, status, types, true, cuts_types, all_sp_stat, all_sp_feas, all_sp_vals, all_sp_duals, all_sp_dual_vals, all_sp_dual_of, all_sp_dual_feas, cb, agg)
        _cuts_gen, _cuts_rel, cuts_sep_other = separate_cuts(mp, αvals, ρ_h, M, data, status, types, true, cuts_types, sp_stats, sp_feas, sp_ofs, sp_duals, cb)
        
        # For MW, use the y values for the first iteration,
        # then calculate a new point between feasible solutions
        # If y is not feasible, then use the same last feasible solution
        if "MW" in types && !(false in sp_feas)
            if !first_it_MW
                interior_y = τ_MW.*interior_y .+ (1-τ_MW).*yvals
            elseif first_it_MW
                if !(false in sp_feas)
                    interior_y = yvals
                else
                    interior_y = gen_y(data, repeat([data.k], data.J))
                end
                first_it_MW = false
            end
            sp_stats_MW, sp_feas_MW, sp_ofs_MW, sp_duals_MW = solve_benders_sp_dual(dualsMW, data, params, status, solver, interior_y, αvals, ρ_h, M, n_outer_cuts, types, update_ρ, w_fc, w0, yvals, sp_ofs)
            _cuts_gen, _cuts_rel, cuts_sep_other = separate_cuts(mp, αvals, ρ_h, M, data, status, ["B"], true, cuts_types, sp_stats_MW, sp_feas_MW, sp_ofs_MW, sp_duals_MW, cb)
                        
            # if !(false in sp_feas_MW)
            #     for t in T
            #         if sp_stats_MW[t] == MOI.OPTIMAL && sp_ofs_MW[t] > αvals[t] + 10e-2 #tol
            #             status.nOptCuts += 1
            #             expr = get_expr(mp, sp_duals_MW[t], data, ρ_h, M, t, types)
            #             opt_cut = @build_constraint(mp[:α][t] >= expr)
            #             MOI.submit(mp, MOI.LazyConstraint(cb), opt_cut)
            #         end
            #     end
            # end
        end

        # For PK, if y is not feasible for all the periods
        # then use the last interior point, otherwise, update it
        if "PK" in types
            if !(false in sp_feas) && !first_it_MW
                interior_y =  τ_PK.*interior_y .+ (1-τ_PK).*yvals #yvals #
            elseif first_it_MW
                if !(false in sp_feas)
                    interior_y = yvals
                else
                    interior_y = gen_y(data, repeat([data.k], data.J))
                end
                first_it_MW = false
            end
            sp_stats_PK, sp_feas_PK, sp_ofs_PK, sp_duals_PK = solve_benders_sp_dual(duals, data, params, status, solver, interior_y, αvals, ρ_h, M, n_outer_cuts, types, false, w_fc, w0) #update_ρ
            _cuts_gen, _cuts_rel, cuts_sep_other = separate_cuts(mp, αvals, ρ_h, M, data, status, ["B"], true, cuts_types, sp_stats_PK, sp_feas_PK, sp_ofs_PK, sp_duals_PK, cb)
            # if !(false in sp_feas)
            #     interior_y = yvals
            # end
            # for t in T
            #     if sp_stats_PK[t] == MOI.OPTIMAL
            #         status.nOptCuts += 1
            #         expr = get_expr(mp, sp_duals_PK[t], data, ρ_h, M, t, types)
            #         opt_cut = @build_constraint(mp[:α][t] >= expr)
            #         MOI.submit(mp, MOI.LazyConstraint(cb), opt_cut)
            #     end
            # end
        end

        ρvals = Array{Float64}(undef, data.J, data.t)    
        Rvals = Array{Float64}(undef, data.J, data.t)
        for t in T
            if "B" in types || "PK" in types || "FC" in types
                ρvals[:,t] = -dual.(duals[t][:cρ_sp_d])
                Rvals[:,t] = -dual.(duals[t][:cR_sp_d])
            elseif "MW" in types
                ρvals[:,t] = -dual.(duals[t][:cρ_sp_d]) #dualsMW
                Rvals[:,t] = -dual.(duals[t][:cR_sp_d]) 
            elseif "SH" in types
                ρvals[:,t] = sp_vars[t]["ρ"].data
                Rvals[:,t] = sp_vars[t]["R"].data
            end
        end
        
        update_ρ = false

        if  !(false in sp_feas)  && 2 ==3
            for j in J, t in T
                # Update if violation
                if Rvals[j,t] + tol < ρvals[j,t]/(1-ρvals[j,t])
                    update_ρ = true
            # for I in findall(==(1), Rvals .+ tol .< (ρvals ./ (1 .- ρvals)))
                # j, t = I[1], I[2]
                    if ρvals[j,t] > maximum(ρ_h[j,t,:])
                        update_M = true
                        println("We will update M @@@@@@@@@@@@@@@@@@@@@@@@@@@@")
                    end
                    # Avoid having Nan in the M calc. This actually happens
                    # if ρvals > 1 - 10e-18, but we decrease the max to limit the
                    # max value of M
                    if ρvals[j,t] > 1 - 10e-4
                        ρvals[j,t] = 1 - 10e-4
                    end
                end
            end        
            
            if update_ρ
                ρ_h = cat(ρ_h, ρvals, dims=3)
                # if update_M    
                #     M = calc_big_M(data, ρ_h)
                # end
            end
        end

        lb_iter[status.nIter] = lb
        ub_iter[status.nIter] = ub

    end
    MOI.set(mp, MOI.LazyConstraintCallback(), lazyCB)

    # Set y vars as binary
    set_binary.(mp[:y])

    println("\n---- STARTING BRANCHING ----")
    n_vars = num_variables(mp)
    end_stat_mp, _of_mp, yval_mp, αval_mp = solve_benders_mp(mp)
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
    
    # println("Y end = $(yval_mp.data)")
    # println("OF end = $_of_mp")
    # println("End status is $end_stat_mp")
    # println("Y OPT = $(yval)")
    
    xval = Array{Float64}(undef, data.I, data.J, data.t)
    ρval = Array{Float64}(undef, data.J, data.t)    
    Rval = Array{Float64}(undef, data.J, data.t)    
    wval = Array{Float64}(undef, data.J, data.k, data.t)
    zval = Array{Float64}(undef, data.J, data.k, data.t)

    if "SH" in types
        sp_vars, _sp_stats, _sp_feas, _sp_ofs, _sp_duals = solve_benders_sp_primal(primals, data, params, status, solver, yval, ρ_h, M, n_outer_cuts, 0, update_ρ)
        for t in T
            vars = sp_vars[t]
            xval[:,:,t] = vars["x"]
            ρval[:,t] = vars["ρ"]
            Rval[:,t] = vars["R"]
            wval[:,:,t] = vars["w"]
            zval[:,:,t] = vars["z"]
        end
    else
        if "B" in types || "MW" in types || "PK" in types
            sp_stats, sp_feas, sp_ofs, sp_duals = solve_benders_sp_dual(duals, data, params, status, solver, yval, αval_mp, ρ_h, M, n_outer_cuts, ["B"], update_ρ, w_fc, w0)
        elseif "FC" in types
            println("y=$yval")
            sp_stats, sp_feas, sp_ofs, sp_duals = solve_benders_sp_dual(duals, data, params, status, solver, yval, repeat([0], data.t), ρ_h, M, n_outer_cuts, types, update_ρ, 0, w0)
            println("Feas $sp_feas")
            # println("Feas $sp_stats")
            # println("Duals $(sp_duals[1])")
            # println("OF = $(sum(sp_ofs)+sum(αval_mp))")
        elseif "MW1" in  types
            sp_stats, sp_feas, sp_ofs, sp_duals = solve_benders_sp_dual(dualsMW, data, params, status, solver, yval, αval_mp, ρ_h, M, n_outer_cuts, ["B"], update_ρ, w_fc, w0)
        end
            
        for t in T
            xval[:,:,t] = -dual.(duals[t][:cx_sp_d])
            ρval[:,t] = -dual.(duals[t][:cρ_sp_d])
            Rval[:,t] = -dual.(duals[t][:cR_sp_d])
            wval[:,:,t] = -dual.(duals[t][:cw_sp_d])
            zval[:,:,t] = -dual.(duals[t][:cz_sp_d])
        end
    end    

    Fterm = sum(F[j,k]*yval[j,k] for j in J for k in K) #dot(F, yval)
    Allocterm = sum(C[i,j,t]*xval[i,j,t] for i in I for j in J for t in T) #dot(C, xval) #sum(C[i,j,t]*xval[i,j,t] for i in I for j in J for t in T)  
    Congterm = 0.5 * sum(Dt[t] * (Rval[j,t] + ρval[j,t] + sum(cv^2 * (wval[j,k,t] - zval[j,k,t]) for k in K)) for j in J for t in T)

    tests_feas = is_sol_feas(data, yval, xval)
    # push!(n_cons, num_constraints(mp; count_variable_in_set_constraints = false))

    return Fterm+Allocterm+Congterm, Fterm, Allocterm, Congterm, yval, xval, status, lb_iter, ub_iter, tests_feas, relax_iters, relax_cuts, n_vars, n_cons, n_nodes, update_M, ρ_h
    ######### End Restricted problem #############
    ##############################################

end

function model_benders(data, params, status, relaxations=["Clust"], types=["B"], n_outer_cuts=32)

    ρ_h = ini_ρ_h(data, n_outer_cuts)  
    M = calc_big_M(data, ρ_h)
    update_M = true
    Fterm= Allocterm= Congterm= yval= xval= lb_iter= ub_iter= tests_feas= relax_iters= relax_cuts= n_vars= n_cons= n_nodes =0
    while update_M
        Fterm+Allocterm+Congterm, Fterm, Allocterm, Congterm, yval, xval, status, lb_iter, ub_iter, tests_feas, relax_iters, relax_cuts, n_vars, n_cons, n_nodes, update_M, ρ_h = model_benders_M(data, params, status, relaxations, types, n_outer_cuts, ρ_h, M)
        if update_M
            M = calc_big_M(data, ρ_h)
        end
        # println("we will update M = $update_M")
    end
    return Fterm+Allocterm+Congterm, Fterm, Allocterm, Congterm, yval, xval, status, lb_iter, ub_iter, tests_feas, relax_iters, relax_cuts, n_vars, n_cons, n_nodes, [n_outer_cuts, size(ρ_h)[3]]
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
    @constraint(mp, [j in J], sum(y[j,k] for k in K) <= 1, set_string_name = false)

    ## Valid inequalities
    # The capacity must be enough to cover all the demand in each period
    @constraint(mp, [t in T], sum(Q[j,k]*y[j,k] for j in J for k in K) >= sum(λ[:,t]), set_string_name = false)

    # Since we are allocating all the demand, the variable costs have to be greatter than or equal to
    # the cost of allocating each node to the least expensive facility (closest one)
    @constraint(mp, α_lb[t in T], α[t] >= floor(sum([minimum(C[i,:,t]) for i in I])), set_string_name = false)

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

function ini_benders_sp_primal(y_fixed, data, ρ_h, M, t, μ, params, status, solver, Solver_ENV, agg=[])
    I = 1:data.I
    J = 1:data.J
    λ = data.a
    C = data.C
    Q = data.Q
    cv = data.cv
    D = data.D    
    K = 1:data.k

    H = 1:size(ρ_h[:,t,:],2)
    Dt = D/sum(λ[:,t])

    maxtime = max(1, params.max_time - elapsed(status))

    if solver == "gurobi"
        m = Model(optimizer_with_attributes(()->Gurobi.Optimizer(Solver_ENV),
        "OutputFlag" => 0,
        "Threads" => 1,
        "MIPGap" => 1e-5,
        "TimeLimit" => maxtime + 1,
        "Presolve" => 0,
        "Cuts" => 0,
        "LogToConsole" => 0,
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
    @variable(m, 0 <= R[J])
    @variable(m, 0 <= ρ[J])
    @variable(m, 0 <= w[J,K])
    @variable(m, 0 <= z[J,K])

    @expression(m, of_sp_allo, sum(C[i,j,t]*x[i,j] for i in I for j in J))
    @expression(m, of_sp_cong, 0.5*Dt*sum((R[j] + ρ[j] + sum(cv^2*(w[j,k]-z[j,k]) for k in K)) for j in J))

    @objective(m, Min, of_sp_allo + of_sp_cong)

    # Capacity cannot be exceeded and steady state has to be conserved
    # Redundant
    # @constraint(m, c1_sp[j in J], -sum(λ[i,t]*x[i,j] for i in I) >= -sum(Q[j,k]*y_fixed[j,k] for k in K) - μ*sum(Q[j,k] for k in K), set_string_name = false)
    
    # All customer zones need to be assigned to exactly one facility    
    @constraint(m, c2_sp[i in I], sum(x[i,j] for j in J) == 1, set_string_name = false)
    
    # c3: Relation x - z
    @constraint(m, c3_sp[j in J], sum(λ[i,t]*x[i,j] for i in I) - sum(Q[j,k]*z[j,k] for k in K) == 0, set_string_name = false)
    # @constraint(m, c3_sp[j in J], sum(λ[i,t]*x[i,j] for i in I) - ρ[j]*sum(Q[j,k]*y_fixed[j,k] for k in K) == 0, set_string_name = false)

    # c4: Linear outer approximation of R
    @constraint(m, c4_sp[j in J, h in H], (1-ρ_h[j,t,h])^2*R[j] - ρ[j] >= -(1+μ)ρ_h[j,t,h]^2, set_string_name = false)
    
    
    # c5: Relation z - ρ
    # c6: Relation w - R
    # No relaxation or "LP"
    if length(agg)==0
        # c5
        @constraint(m, c5_sp[j in J, k in K], z[j,k] - ρ[j]*y_fixed[j,k] == 0, set_string_name = false)
        # c6
        @constraint(m, c6_sp[j in J, k in K], w[j,k] - R[j]*y_fixed[j,k] == 0, set_string_name = false)
    else
        # Aggregate on J
        if "J" in agg
            # c5
            @constraint(m, c5_sp[k in K], sum(z[j,k] for j in J) - sum(ρ[j]*y_fixed[j,k] for j in J) == 0, set_string_name = false)
            # c6
            @constraint(m, c6_sp[k in K], sum(w[j,k] for j in J) - sum(R[j]*y_fixed[j,k] for j in J) == 0, set_string_name = false)
        # Aggregate on K
        elseif "K" in agg
            # c5
            @constraint(m, c5_sp[j in J], sum(z[j,k] for k in K) - sum(ρ[j]*y_fixed[j,k] for k in K) == 0, set_string_name = false)
            # c6
            @constraint(m, c6_sp[j in J], sum(w[j,k] for k in K) - sum(R[j]*y_fixed[j,k] for k in K) == 0, set_string_name = false)
        # Aggregate on K and j
        elseif "KJ" in agg
            # c5
            @constraint(m, c5_sp, sum(z[j,k] for k in K for j in J) - sum(ρ[j]*y_fixed[j,k] for k in K for j in J) == 0, set_string_name = false)
            # c6
            @constraint(m, c6_sp, sum(w[j,k] for k in K for j in J) - sum(R[j]*y_fixed[j,k] for k in K for j in J) == 0, set_string_name = false)
        # Assume is a vector of vectors, with the cluster results
        else
            for n_clus in eachindex(agg)
                set_Rl = agg[n_clus]
                # c5
                m[Symbol("c5clus_sp_$n_clus")] = @constraint(m, sum(z[j,k] for (j,k) in set_Rl) - sum(ρ[j]*y_fixed[j,k] for (j,k) in set_Rl) == 0, set_string_name = false)
                # c6
                m[Symbol("c6clus_sp_$n_clus")] = @constraint(m, sum(w[j,k] for (j,k) in set_Rl) - sum(R[j]*y_fixed[j,k] for (j,k) in set_Rl) == 0, set_string_name = false)
            end
        end
    end
   
    # Upper bounds
    @constraint(m, c7_sp[i in I, j in J], -x[i,j] >= -1 - μ, set_string_name = false)
    # @constraint(m, c8_sp[j in J, k in K], -z[j,k] >= -1 - μ, set_string_name = false)
    @constraint(m, c9_sp[j in J], -ρ[j] >= -(1-10e-7) - μ, set_string_name = false)

    return m
end

function solve_benders_sp_primal(primals, data, params, status, solver, y_fixed, ρ_h, M, n_outer_cuts, μ, update_ρ=false, agg=[])
    T = 1:data.t
    J = 1:data.J
    sp_stats = []
    sp_feas = []
    sp_ofs = []
    sp_duals = Dict()
    sp_vars = Dict()

    for t in T
        update_sp_primal(primals[t], data, params, status, solver, y_fixed, t, n_outer_cuts, μ, ρ_h, M, update_ρ, agg)            
        
        optimize!(primals[t])
        end_stat = termination_status(primals[t])

        sp_duals[t] = Dict()
        sp_vars[t] = Dict()
        
        try # if has_values(primals[t])#primal_status(primals[t]) == MOI.FEASIBLE_POINT
            xval = value.(primals[t][:x])
            ρval = value.(primals[t][:ρ])        
            Rval = value.(primals[t][:R])
            wval = value.(primals[t][:w])
            zval = value.(primals[t][:z])

            # π1val = dual.(primals[t][:c1_sp])
            π2val = dual.(primals[t][:c2_sp])
            
            # c4 - π4
            n_tan = size(ρ_h[:,t,:],2)
            π4val = Array{Float64}(undef, data.J, n_tan)
            π4val[:,1:n_outer_cuts] = dual.(primals[t][:c4_sp])
            if  n_tan > n_outer_cuts
                for j in J, h in (n_outer_cuts+1):n_tan
                    π4val[j,h] = dual.(primals[t][Symbol("c4_$(j)_$(t)_$(h)")])
                end
            end

            π7val = dual.(primals[t][:c7_sp])
            # π8val = dual.(primals[t][:c8_sp])
            π9val = dual.(primals[t][:c9_sp])            

            # Store the results
            # sp_duals[t]["π1"] = π1val
            sp_duals[t]["π2"] = π2val
            sp_duals[t]["π4"] = π4val
            sp_duals[t]["π7"] = π7val
            # sp_duals[t]["π8"] = π8val
            sp_duals[t]["π9"] = π9val

            sp_vars[t]["x"] = xval
            sp_vars[t]["R"] = Rval
            sp_vars[t]["ρ"] = ρval
            sp_vars[t]["w"] = wval
            sp_vars[t]["z"] = zval

            push!(sp_stats, end_stat)

            if end_stat == MOI.DUAL_INFEASIBLE || end_stat == MOI.INFEASIBLE
                push!(sp_feas, false)
                Allocterm, Congterm = 10e8, 10e8
            elseif end_stat == MOI.OPTIMAL
                push!(sp_feas, true)
                Allocterm, Congterm = value(primals[t][:of_sp_allo]),  value(primals[t][:of_sp_cong])
            end
            push!(sp_ofs, Allocterm + Congterm)
        catch e #else
            push!(sp_feas, false)
            push!(sp_ofs, -1)
            push!(sp_stats, end_stat)
        end
    end
    return sp_vars, sp_stats, sp_feas, sp_ofs, sp_duals
end

function update_sp_primal(primal, data, params, status, solver, y_fixed, t, n_outer_cuts, μ, ρ_h, M, update_ρ=false, agg=[])
    maxtime = max(1, params.max_time - elapsed(status))
    if solver == "gurobi"
        set_attribute(primal, "TimeLimit", maxtime)
    elseif solver == "cplex"
        set_attribute(primal, "CPXPARAM_TimeLimit", maxtime)
    end

    J = 1:data.J  
    K = 1:data.k
    Q = data.Q

    # c1
    # set_normalized_rhs([primal[:c1_sp][j] for j in J], [-sum(Q[j,k]*y_fixed[j,k] for k in K) - μ*sum(Q[j,k] for k in K) for j in J])

    # c3
    # set_normalized_coefficient([primal[:c3_sp][j] for j in J], [primal[:ρ][j] for j in J], [-sum(Q[j,k]*y_fixed[j,k] for k in K) for j in J])

    # c4
    if update_ρ=="asdsad"
        n_tan = size(ρ_h[:,t,:],2)
        # Adding the extra c4 cuts
        for j in J
            primal[Symbol("c4_$(j)_$(t)_$(n_tan)")] = @constraint(primal, (1-ρ_h[j,t,n_tan])^2 * primal[:R][j] - primal[:ρ][j] >= -(1+μ)*ρ_h[j,t,n_tan]^2)
        end
    end

    # TODO modificar las otras agregaciones
    if length(agg)==0
        # for k in K
        #     for j in J
        #         set_normalized_coefficient(primal[:c5_sp][j,k], primal[:ρ][j], -y_fixed[j,k])
        #         set_normalized_coefficient(primal[:c6_sp][j,k], primal[:R][j], -y_fixed[j,k])
        #     end
        # end
        set_normalized_coefficient([primal[:c5_sp][j,k] for k in K for j in J], [primal[:ρ][j] for k in K for j in J], [-y_fixed[j,k] for k in K for j in J])
        set_normalized_coefficient([primal[:c6_sp][j,k] for k in K for j in J], [primal[:R][j] for k in K for j in J], [-y_fixed[j,k] for k in K for j in J])
    elseif "J" in agg
        for k in K
            set_normalized_rhs(primal[:c5_sp][k], sum(-y_fixed[j,k] - μ for j in J))
            set_normalized_rhs(primal[:c6_sp][k], sum(-M[j,t]*y_fixed[j,k] - μ*M[j,t] for j in J))
        end
    elseif "K" in agg
        for j in J
            set_normalized_rhs(primal[:c5_sp][j], sum(-y_fixed[j,k] - μ for k in K))
            set_normalized_rhs(primal[:c6_sp][j], sum(-M[j,t]*y_fixed[j,k] - μ*M[j,t] for k in K))
        end
    elseif "KJ" in agg
        set_normalized_rhs(primal[:c5_sp], sum(-y_fixed[j,k] - μ for k in K for j in J))
        set_normalized_rhs(primal[:c6_sp], sum(-M[j,t]*y_fixed[j,k] - μ*M[j,t] for k in K for j in J))
    else
        for n_clus in eachindex(agg)
            set_Rl = agg[n_clus]
            # set_normalized_rhs(primal[Symbol("c5clus_sp_$n_clus")], sum(-y_fixed[j,k] - μ for (j,k) in set_Rl))
            # set_normalized_rhs(primal[Symbol("c6clus_sp_$n_clus")], sum(-M[j,t]*y_fixed[j,k] - μ*M[j,t] for (j,k) in set_Rl))
            for (j,k1) in set_Rl
                m[Symbol("c5clus_sp_$n_clus")] = @constraint(m, sum(z[j,k] for (j,k) in set_Rl) - sum(ρ[j]*y_fixed[j,k] for (j,k) in set_Rl) == 0, set_string_name = false)
                set_normalized_coefficient(primal[Symbol("c5clus_sp_$n_clus")], primal[:ρ][j], sum(-y_fixed[j,k] for (j1,k) in set_Rl))
                set_normalized_coefficient(primal[Symbol("c5clus_sp_$n_clus")], primal[:R][j], sum(-y_fixed[j,k] for (j1,k) in set_Rl))
            end
        end
    end
end

function ini_benders_sp_dual(y_fixed, data, ρ_h, M, t, params, status, solver, Solver_ENV, types, w_fc=0, w0=1)
    I = 1:data.I
    J = 1:data.J
    λ = data.a
    C = data.C
    Q = data.Q
    cv = data.cv
    D = data.D    
    K = 1:data.k
    
    H = 1:size(ρ_h[:,t,:],2)   
    Dt = D/sum(λ[:,t])

    maxtime = max(1, params.max_time - elapsed(status))

    if solver == "gurobi"
        m = Model(optimizer_with_attributes(()->Gurobi.Optimizer(Solver_ENV),
        "OutputFlag" => 0,
        "Threads" => 1,
        "MIPGap" => 1e-5,
        "TimeLimit" => maxtime + 1,
        "Presolve" => 0,
        "Cuts" => 0,
        "LogToConsole" => 0,
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
    if "MW" in types || "PK" in types
        @expression(m, dual_of, sum(π1[j]*sum(-Q[j,k]*y_fixed[j,k] for k in K) for j in J))
        add_to_expression!(dual_of, sum(π2[i] for i in I))
        add_to_expression!(dual_of, sum(-π4[j,k]*y_fixed[j,k] for j in J for k in K))
        add_to_expression!(dual_of, sum(-π6[j,k]*M[j,t]*y_fixed[j,k] for j in J for k in K))
        add_to_expression!(dual_of, sum(-π8[j,h]*(ρ_h[j,t,h]^2) for j in J for h in H))
        add_to_expression!(dual_of, sum(-π10[i,j] for i in I for j in J))
        add_to_expression!(dual_of, sum(-π11[j,k] for j in J for k in K))
        add_to_expression!(dual_of, sum(-π12[j] for j in J)) #(1-10e-7)*

        @objective(m, Max, dual_of)

        @constraint(m, cx_sp_d[i in I, j in J], -π1[j]*λ[i,t] + π2[i] + π3[j]*λ[i,t] - π10[i,j]  <= C[i,j,t], set_string_name = false)
        @constraint(m, cR_sp_d[j in J], -π7[j] + sum(π8[j,h]*(1-ρ_h[j,t,h])^2 for h in H)  <= 0.5*Dt, set_string_name = false)
        @constraint(m, cρ_sp_d[j in J], -π5[j] + sum(-π8[j,h] for h in H) - π12[j]  <= 0.5*Dt, set_string_name = false)
        @constraint(m, cw_sp_d[j in J, k in K], -π6[j,k] + π7[j] <= 0.5*Dt*cv^2, set_string_name = false)
        @constraint(m, cz_sp_d[j in J, k in K], -π3[j]*Q[j,k] - π4[j,k] + π5[j] - π11[j,k] <= -0.5*Dt*cv^2, set_string_name = false)

        
        if "MW" in types
            @constraint(m, consMW, sum(-π1[j]*sum(Q[j,k]*y_fixed[j,k] for k in K) for j in J)+
            sum(π2[i] for i in I)+
            sum(-π4[j,k]*y_fixed[j,k] for j in J for k in K)+
            sum(-π6[j,k]*M[j,t]*y_fixed[j,k] for j in J for k in K)+
            sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)+
            sum(-π10[i,j] for i in I for j in J)+
            sum(-π11[j,k] for j in J for k in K)+
            sum(-π12[j] for j in J) == 1, set_string_name = false)
        end
        return m

    else
        # General ##################        
        w1  = w_fc # [w_fc for j in J]
        w2 = w_fc # [w_fc for i in I]
        w4 = w_fc # [w_fc for j in J, k in K]
        w6 = w_fc # [w_fc for j in J, k in K]
        w8 = w_fc # [w_fc for j in J, h in H]
        w10 = w_fc # [w_fc for i in I, j in J]
        w11 = w_fc # [w_fc for j in J, k in K]
        w12 = w_fc # [w_fc for j in J]

        @variable(m, 0 <= π0)

        @expression(m, dual_of, sum(π1[j]*sum(-Q[j,k]*y_fixed[j,k] for k in K) for j in J))
        add_to_expression!(dual_of, sum(π2[i] for i in I))
        add_to_expression!(dual_of, sum(-π4[j,k]*y_fixed[j,k] for j in J for k in K))
        add_to_expression!(dual_of, sum(-π6[j,k]*M[j,t]*y_fixed[j,k] for j in J for k in K))
        add_to_expression!(dual_of, sum(-π8[j,h]*(ρ_h[j,t,h]^2) for j in J for h in H))
        add_to_expression!(dual_of, sum(-π10[i,j] for i in I for j in J))
        add_to_expression!(dual_of, sum(-π11[j,k] for j in J for k in K))
        add_to_expression!(dual_of, sum(-π12[j] for j in J))
        add_to_expression!(dual_of, -π0)

        @objective(m, Max, dual_of)

        # @objective(m, Max, sum(π1[j]*sum(-Q[j,k]*yvals[j,k] for k in K) for j in J)+
        # sum(π2[i] for i in I)+
        # sum(-π4[j,k]*yvals[j,k] for j in J for k in K)+
        # sum(-π6[j,k]*M[j,t]*yvals[j,k] for j in J for k in K)+
        # sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)+
        # sum(-π10[i,j] for i in I for j in J)+
        # sum(-π11[j,k] for j in J for k in K)+
        # sum(-π12[j] for j in J) - π0*αvals[t])

        @constraint(m, cx_sp_d[i in I, j in J], -π1[j]*λ[i,t] + π2[i] + π3[j]*λ[i,t] - π10[i,j]  <= π0*C[i,j,t], set_string_name = false)
        @constraint(m, cR_sp_d[j in J], -π7[j] + sum(π8[j,h]*(1-ρ_h[j,t,h])^2 for h in H)  <= π0*0.5*Dt, set_string_name = false)
        @constraint(m, cρ_sp_d[j in J], -π5[j] + sum(-π8[j,h] for h in H) - π12[j]  <= π0*0.5*Dt, set_string_name = false)
        @constraint(m, cw_sp_d[j in J, k in K], -π6[j,k] + π7[j] <= π0*0.5*Dt*cv^2, set_string_name = false)
        @constraint(m, cz_sp_d[j in J, k in K], -π3[j]*Q[j,k] - π4[j,k] + π5[j] - π11[j,k] <= -π0*0.5*Dt*cv^2, set_string_name = false)

        # @constraint(m, c_norm_d, sum(w1[j]*π1[j] + w12[j]*π12[j] for j in J) + sum(w2[i]*π2[i] for i in I) + 
        # sum(w4[j,k]*π4[j,k] + w6[j,k]*π6[j,k] + w11[j,k]*π11[j,k] for j in J for k in K) +
        # sum(w8[j, h]*π8[j,h] for j in J for h in H) + sum(w10[i,j]*π10[i,j] for i in I for j in J) + w0*π0 == 1, set_string_name = false)
        @constraint(m, c_norm_d, sum(w1*π1[j] + w12*π12[j] for j in J) + sum(w2*π2[i] for i in I) + 
        sum(w4*π4[j,k] + w6*π6[j,k] + w11*π11[j,k] for j in J for k in K) +
        sum(w8*π8[j,h] for j in J for h in H) + sum(w10*π10[i,j] for i in I for j in J) + w0*π0 == 1, set_string_name = false)
        
        # if "MW" in types
        #     @constraint(m, consMW, sum(-π1[j]*sum(Q[j,k]*y_fixed[j,k] for k in K) for j in J)+
        #     sum(π2[i] for i in I)+
        #     sum(-π4[j,k]*y_fixed[j,k] for j in J for k in K)+
        #     sum(-π6[j,k]*M[j,t]*y_fixed[j,k] for j in J for k in K)+
        #     sum(-π8[j,h]*ρ_h[j,t,h]^2 for j in J for h in H)+
        #     sum(-π10[i,j] for i in I for j in J)+
        #     sum(-π11[j,k] for j in J for k in K)+
        #     sum(-π12[j] for j in J) == 1, set_string_name = false)
        # end
    end
    return m 
end

function solve_benders_sp_dual(duals, data, params, status, solver, y_fixed, αvals, ρ_h, M, n_outer_cuts, types, update_ρ=false, w_fc=0, w0=1, y_ref=[], sp_costs=[])
    T = 1:data.t
    J = 1:data.J
    sp_stats = []    
    sp_feas = []
    sp_ofs = []
    sp_duals = Dict()
    # println("In sp solve, update ρ = $update_ρ")

    for t in T
        update_sp_dual(duals[t], data, params, status, solver, y_fixed, αvals, t, n_outer_cuts, ρ_h, M, types, update_ρ, w_fc, w0, y_ref, sp_costs)
        
        optimize!(duals[t])
        end_stat = termination_status(duals[t])
        # println("Status dual $end_stat")

        try        
            π0val = 0
            if "FC" in types
                π0val = value.(duals[t][:π0])
            end
            π1val = value.(duals[t][:π1])
            π2val = value.(duals[t][:π2])
            π3val = value.(duals[t][:π3])
            π5val = value.(duals[t][:π5])
            π7val = value.(duals[t][:π7])
            π4val = value.(duals[t][:π4])
            π6val = value.(duals[t][:π6])       
            π10val = value.(duals[t][:π10])
            π11val = value.(duals[t][:π11])
            π12val = value.(duals[t][:π12])
            # Constraint 8
            n_tan = size(ρ_h[:,t,:],2)
            π8val = Array{Float64}(undef, data.J, n_tan)
            π8val[:,1:n_outer_cuts] = value.(duals[t][:π8])
            if  n_tan > n_outer_cuts
                for j in J, h in (n_outer_cuts+1):n_tan
                    π8val[j,h] = value.(duals[t][Symbol("π8_$(j)_$(t)_$(h)")])
                end
            end

            # Store the results
            sp_duals[t] = Dict()
            sp_duals[t]["π0"] = π0val
            sp_duals[t]["π1"] = π1val
            sp_duals[t]["π2"] = π2val
            sp_duals[t]["π3"] = π3val
            sp_duals[t]["π5"] = π5val
            sp_duals[t]["π7"] = π7val
            sp_duals[t]["π4"] = π4val        
            sp_duals[t]["π6"] = π6val
            sp_duals[t]["π8"] = π8val
            sp_duals[t]["π10"] = π10val
            sp_duals[t]["π11"] = π11val
            sp_duals[t]["π12"] = π12val

            push!(sp_stats, end_stat)

            if end_stat == MOI.DUAL_INFEASIBLE || end_stat == MOI.INFEASIBLE
                push!(sp_feas, false)
                push!(sp_ofs, 2*10e8)            
            elseif end_stat == MOI.OPTIMAL
                push!(sp_feas, true)  
                push!(sp_ofs, value(duals[t][:dual_of]))                      
            end
        catch e
            # @warn "Invalid input: $e"
            push!(sp_feas, false)
            push!(sp_ofs, -1)
            push!(sp_stats, end_stat)
        end       
    end
    return sp_stats, sp_feas, sp_ofs, sp_duals 
end

function update_sp_dual(dual, data, params, status, solver, y_fixed, αvals, t, n_outer_cuts, ρ_h, M, types, update_ρ=false, w_fc=0, w0=1, y_ref=[], sp_costs=[])
    maxtime = max(1, params.max_time - elapsed(status))
    if solver == "gurobi"
        set_attribute(dual, "TimeLimit", maxtime)
    elseif solver == "cplex"
        set_attribute(dual, "CPXPARAM_TimeLimit", maxtime)
    end

    I = 1:data.I
    J = 1:data.J 
    λ = data.a
    K = 1:data.k
    D = data.D
    Q = data.Q       
    Dt = D/sum(λ[:,t])

    if "FC" in types
        set_objective_coefficient(dual, dual[:π0], -αvals[t])
    end

    set_objective_coefficient(dual, [dual[:π1][j] for j in J], [sum(-Q[j,k]*y_fixed[j,k] for k in K) for j in J])
    set_objective_coefficient(dual, [dual[:π4][j,k] for j in J for k in K], [-y_fixed[j,k] for j in J for k in K])
    set_objective_coefficient(dual, [dual[:π6][j,k] for j in J for k in K], [-M[j,t]*y_fixed[j,k] for j in J for k in K])

    if "MW" in types        
        for j in J
            set_normalized_coefficient(dual[:consMW], dual[:π1][j], -sum(data.Q[j,k]*y_ref[j,k] for k in K))
            for k in K
                set_normalized_coefficient(dual[:consMW], dual[:π4][j,k], -y_ref[j,k])
                set_normalized_coefficient(dual[:consMW], dual[:π6][j,k], -M[j,t]*y_ref[j,k])
            end
        end
        set_normalized_rhs(dual[:consMW], sp_costs[t])
    end

    if update_ρ
        #OF
        of_expr = objective_function(dual)
        n_tan = size(ρ_h[:,t,:],2)
        # Adding the extra π8 vars
        for j in J 
            dual[Symbol("π8_$(j)_$(t)_$(n_tan)")] = @variable(dual, lower_bound = 0, set_string_name = false)
            add_to_expression!(of_expr, -dual[Symbol("π8_$(j)_$(t)_$(n_tan)")]*ρ_h[j,t,n_tan])
        end
        @objective(dual, Max, of_expr)

        # R and ρ constraints
        delete.(dual, dual[:cR_sp_d])
        delete.(dual, dual[:cρ_sp_d])
        unregister(dual, :cR_sp_d)
        unregister(dual, :cρ_sp_d)

        n_tan = size(ρ_h[:,t,:],2)
        @constraint(dual, cR_sp_d[j in J], -dual[:π7][j] + sum(dual[:π8][j,h]*(1-ρ_h[j,t,h])^2 for h in 1:n_outer_cuts) + sum(dual[Symbol("π8_$(j)_$(t)_$(h)")]*(1-ρ_h[j,t,h])^2 for h in (n_outer_cuts+1):n_tan)  <= 0.5*Dt, set_string_name = false)
        @constraint(dual, cρ_sp_d[j in J], -dual[:π5][j] + sum(-dual[:π8][j,h] for h in 1:n_outer_cuts) +sum(-dual[Symbol("π8_$(j)_$(t)_$(h)")] for h in (n_outer_cuts+1):n_tan) - dual[:π12][j]  <= 0.5*Dt, set_string_name = false)

        if "FC" in types
            # Normalization constraint
            w1  = w_fc # [w_fc for j in J]
            w2 = w_fc # [w_fc for i in I]
            w4 = w_fc # [w_fc for j in J, k in K]
            w6 = w_fc # [w_fc for j in J, k in K]
            w8 = w_fc # [w_fc for j in J, h in H]
            w10 = w_fc # [w_fc for i in I, j in J]
            w11 = w_fc # [w_fc for j in J, k in K]
            w12 = w_fc # [w_fc for j in J]

            delete(dual, dual[:c_norm_d])
            unregister(dual, :c_norm_d)

            @constraint(dual, c_norm_d, sum(w1*dual[:π1][j] + w12*dual[:π12][j] for j in J) + sum(w2*dual[:π2][i] for i in I) + 
            sum(w4*dual[:π4][j,k] + w6*dual[:π6][j,k] + w11*dual[:π11][j,k] for j in J for k in K) +
            sum(w8*dual[:π8][j,h] for j in J for h in 1:n_outer_cuts) + sum(w8*dual[Symbol("π8_$(j)_$(t)_$(h)")] for j in J for h in (n_outer_cuts+1):n_tan) + 
            sum(w10*dual[:π10][i,j] for i in I for j in J) + w0*dual[:π0] == 1, set_string_name = false)
        end
    end
    return nothing
end

function separate_cuts(m, αvals, ρ_h, M, data, status, types, lazy=false, cuts_types="feas", sp_stats=[], sp_feas=[], sp_ofs=[], sp_duals=[], cb=[], agg=[])
    I = 1:data.I
    J = 1:data.J
    Q = data.Q
    K = 1:data.k
    T = 1:data.t

    cuts_gen = false
    cuts_sep = 0
    cuts_sep_other = 0
    
    for t in T
        sp_stat = sp_stats[t]
        sp_of = sp_ofs[t]
        sp_du = sp_duals[t]

        expr = get_expr(m, sp_du, data, ρ_h, M, t, types, agg)
        
        # println("Expr = $expr")
        if !("FC" in types)
            # Add a feasibility cut
            if (sp_stat==MOI.INFEASIBLE || sp_stat==MOI.DUAL_INFEASIBLE || sp_stat==MOI.INFEASIBLE_OR_UNBOUNDED)
                status.nFeasCuts += 1
                cuts_sep += 1
                if lazy
                    feas_cut = @build_constraint(0 >= expr)  
                    MOI.submit(m, MOI.LazyConstraint(cb), feas_cut)
                else
                    feas_cut = @constraint(m, 0 >= expr) 
                    m[Symbol("feas_$(t)_$(status.nIter)")] = feas_cut
                end
                println("iter= $(status.nIter) adding feasibility cut for t=$t ######################")
            
            elseif sp_stat == MOI.OPTIMAL && !("MW" in types) && !("PK" in types) && cuts_types == "all"
                # Add an optimality cut
                if  sp_of > αvals[t] + 10e-2 #tol
                    cuts_gen = true
                    
                    ##### Classical Benders ######
                    status.nOptCuts += 1
                    cuts_sep += 1
                    if lazy
                        opt_cut_gen = @build_constraint(m[:α][t] >= expr)
                        MOI.submit(m, MOI.LazyConstraint(cb), opt_cut_gen)
                    else
                        opt_cut_gen = @constraint(m, m[:α][t] >= expr)
                        m[Symbol("opt_Gen_$(t)_$(status.nIter)")] = opt_cut_gen
                    end
                    println("iter= $(status.nIter) adding optimality cut for t=$t --------------------")
                else
                    # println("Solution is opt for $t")
                end
            end
        
        ### Fischetti ###
        elseif "FC" in types
            #expr = get_expr(m, sp_du, data, ρ_h, M, t, types, agg) # get_expr(m, dual_sp_vars, data, ρ_h, M, t, types, agg)
            if sp_stat != MOI.OPTIMAL #(sp_stat==MOI.INFEASIBLE || sp_stat==MOI.DUAL_INFEASIBLE || sp_stat==MOI.INFEASIBLE_OR_UNBOUNDED)
                π0 = sp_du["π0"]
                println("π0= $π0")
                # if sp_stat == MOI.OPTIMAL
                #     π0 = sp_du["π0"]
                # else
                #     π0 = 0
                # end
                if π0 != 0
                    status.nOptCuts += 1
                else
                    status.nFeasCuts += 1
                end
                cuts_sep += 1
                cuts_sep_other += 1
                if lazy
                    opt_cut = @build_constraint(π0*m[:α][t] >= expr)
                    MOI.submit(m, MOI.LazyConstraint(cb), opt_cut)
                else
                    opt_cut = @constraint(m, π0*m[:α][t] >= expr)
                    m[Symbol("cut_FC_$(t)_$(status.nIter)")] = opt_cut
                end
            end
            
        end
    end
    
    return cuts_gen, cuts_sep, cuts_sep_other
end

function get_expr(mp, duals, data, ρ_h, M, t, types, agg=[])
    I = 1:data.I
    J = 1:data.J
    Q = data.Q
    K = 1:data.k
    T = 1:data.t
    H = 1:size(ρ_h[:,t,:],2)

    # π1 = duals["π1"]
    π2 = duals["π2"]
    π4 = duals["π4"]
    π7 = duals["π7"]
    # π8 = duals["π8"]
    π9 = duals["π9"]

    expr = AffExpr(0) 
    # add_to_expression!(expr, sum(π1[j]*sum(-Q[j,k]*mp[:y][j,k] for k in K) for j in J))    
    add_to_expression!(expr, sum(π2[i] for i in I))
    add_to_expression!(expr, sum(-π4[j,h]*(ρ_h[j,t,h]^2) for j in J for h in H))
    add_to_expression!(expr, sum(-π7[i,j] for i in I for j in J))
    # add_to_expression!(expr, sum(-π8[j,k] for j in J for k in K))
    add_to_expression!(expr, sum(-π9[j]*(1-10e-7) for j in J))

    return expr
end

function calc_viol!(data, yvals, wvals, zvals, M, agg)
    T = 1:data.t

    for t in T
        zvalst = zvals[:,:,t]
        wvalst = wvals[:,:,t]
        for Rl in agg
            new_Rl4 = []
            new_Rl6 = []
            for (j,k) in Rl
                if zvalst[j,k] > yvals[j,k]
                    push!(new_Rl4, (j,k))
                end
                if wvalst[j,k] > M[j,t]*yvals[j,k]
                    push!(new_Rl6, (j,k))
                end
            end
            if length(new_Rl4) > 0 && !(new_Rl4 in agg)
                push!(agg, new_Rl4)
            end
            if length(new_Rl6) > 0 && !(new_Rl6 in agg)
                push!(agg, new_Rl6)
            end            
        end
    end  
end