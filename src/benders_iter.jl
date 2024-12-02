function model_benders_iter(data, params, status, relaxations=["Clust"], types=["B"], n_outer_cuts=32)
    solver = "gurobi" #cplex or gurobi
    # Solver environment
    if solver == "gurobi"
        Solver_ENV = Gurobi.Env()
    elseif solver == "cplex" 
        Solver_ENV = CPLEX.Env()
    end

    # # Parameters for the different acceleration techniques
    # # Magnanti-Wang & Papadakos
    # first_it_MW = true
    # interior_y = gen_y(data, repeat([data.k], data.J))
    # τ_MW = 0.5
    # τ_PK = 0.5
    
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
    # and of the number of vars and constraints
    n_vars, n_cons, n_nodes = 0, [], 0

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
    # n_outer_cuts = 32
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
    sp_m = Dict()

    # Find the min cap level for wich the total cap of all facilities
    # is greatter than the max demand of all the periods
    # cap_0 = findfirst(>=(maximum(sum(λ, dims=1))), sum(Q, dims=1))[2]
    # println("Min cap $cap_0")
    
    int_y_ind = repeat([0], data.J) #[0 for j in J]  #[data.k for j in J] # repeat([cap_0], data.J)  #
    int_y = gen_y(data, int_y_ind)   


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
    n_bounds_rel = 5 # Number of bounds to consider for calculating the Convergence
    max_iters_rel = 1000

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
    rel_to_run = []
    max_relax = []
    if "K" in relaxations
        push!(rel_to_run, ["K"]) # = [types]
        push!(max_relax, 1)
        # sp_to_test = [primals_k]
    elseif "J" in relaxations
        push!(rel_to_run, ["J"]) 
        push!(max_relax, 1)
    elseif "KJ" in relaxations
        push!(rel_to_run, ["KJ"]) 
        push!(max_relax, 1)
    elseif "LP" in relaxations
        push!(rel_to_run, [])
        push!(max_relax, 1)
    elseif "Clust" in relaxations
        push!(max_relax, 5)
        n_clus = max(2,Int64(floor(data.J/10)))
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
    if length(rel_to_run) > 0 
        for (agg, max_rel) in zip(rel_to_run, max_relax)
            while n_iter_rel_clus <= max_rel && keep_relax

                prim_rel = Dict()
                for t in T
                    prim_rel[t] = ini_benders_sp_primal(int_y, data, ρ_h, M, t, 0, params, status, solver, Solver_ENV, agg)
                end
                                                                                                                                                       
                lb_temp, ub_temp, lb_iter, ub_iter, iter_rel, cuts_rel, yval_rel, _αval_rel, xval_rel, ρval_rel, Rval_rel, wval_rel, zval_rel, ρ_h = benders_iter(mp, prim_rel, ["Rel"], ρ_h, M, data, params, status, solver, lb, ub, lb_iter, ub_iter, tol, n_outer_cuts, agg, 0, int_y, w_fc, τ, n_bounds_rel, max_iters_rel)
                push!(relax_iters, iter_rel)
                push!(relax_cuts, cuts_rel)        
                push!(relax_lb, lb_temp)
                println("With this agg the bounds are lb=$lb_temp and ub=$ub_temp after $cuts_rel cuts\n")

                if lb < lb_temp
                    lb = lb_temp
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

                n_iter_rel_clus += 1
            end
        end
    end

    # Set y vars as binary
    set_binary.(mp[:y])

    # Initialize the subproblems
    for t in T
        if "SH" in types
            sp_m[t] = ini_benders_sp_primal(int_y, data, ρ_h, M, t, μ, params, status, solver, Solver_ENV)
        elseif "B" in types || "PK" in types || "FC" in types
            sp_m[t] = ini_benders_sp_dual(int_y, data, ρ_h, M, t, params, status, solver, Solver_ENV, types, w_fc, w0)
        elseif "MW" in types
            sp_m[t] = ini_benders_sp_dual(int_y, data, ρ_h, M, t, params, status, solver, Solver_ENV, ["B"], w_fc, w0)
            dualsMW[t] = ini_benders_sp_dual(int_y, data, ρ_h, M, t, params, status, solver, Solver_ENV, types, w_fc, w0)
        end
    end

    n_bounds = 100
    max_iters = 10e6
    lb_temp, ub_temp, lb_iter, ub_iter, iter_rel, cuts_rel, yval, _αval, xval, ρval, Rval, wval, zval, ρ_h = benders_iter(mp, sp_m, types, ρ_h, M, data, params, status, solver, lb, ub, lb_iter, ub_iter, tol, n_outer_cuts, [], μ, int_y, w_fc, τ, n_bounds, max_iters, dualsMW)     
    
    push!(n_cons, num_constraints(mp; count_variable_in_set_constraints = true))
    
    Fterm = sum(F[j,k]*yval[j,k] for j in J for k in K)
    Allocterm = sum(C[i,j,t]*xval[i,j,t] for i in I for j in J for t in T)
    Congterm = 0.5 * sum(Dt[t] * (Rval[j,t] + ρval[j,t] + sum(cv^2 * (wval[j,k,t] - zval[j,k,t]) for k in K)) for j in J for t in T)

    tests_feas = is_sol_feas(data, yval, xval)

    return Fterm+Allocterm+Congterm, Fterm, Allocterm, Congterm, yval, xval, status, lb_iter, ub_iter, tests_feas, relax_iters, relax_cuts, n_vars, n_cons, n_nodes, [n_outer_cuts, size(ρ_h)[3]]
end

function benders_iter(mp, sp_m, types, ρ_h, M, data, params, status, solver, lb, ub, lb_iter, ub_iter, tol, n_outer_cuts, agg=[], μ=0, int_y=[], w_fc=0, τ=10e-5, n_bounds=5, max_iter=10e5, dualsMW=[])

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

    # Parameters for the different acceleration techniques
    # Magnanti-Wang & Papadakos
    first_it_MW = true
    interior_y = gen_y(data, repeat([data.k], data.J))
    τ_MW = 0.5
    τ_PK = 0.5
    w0 = 1
    
    # Convergence test
    last_bounds = zeros(n_bounds)
    conv = false
    lb_mp, ub_mp = lb, ub
    # Variables
    yvals_mp, αvals_mp = [], []
    xval = Array{Float64}(undef, data.I, data.J, data.t)
    ρval = Array{Float64}(undef, data.J, data.t)    
    Rval = Array{Float64}(undef, data.J, data.t)    
    wval = Array{Float64}(undef, data.J, data.k, data.t)
    zval = Array{Float64}(undef, data.J, data.k, data.t)

    iter_rel = 0
    cuts_rel = 0
    update_ρ = false

    # Solver_ENV = Gurobi.Env() # TESTING
    # solver="gurobi"
    
    # Solve the problem using Benders iteratively
    while !conv && iter_rel < max_iter
        status.nIter += 1
        iter_rel += 1

        # For PK, if y is not feasible for all the periods
        # then use the last interior point, otherwise, update it
        if "PK" in types
            if !first_it_MW
                interior_y =  τ_PK.*interior_y .+ (1-τ_PK).*yvals_mp
            else
                first_it_MW = false
            end
            sp_stats_PK, sp_feas_PK, sp_ofs_PK, sp_duals_PK = solve_benders_sp_dual(sp_m, data, params, status, solver, interior_y, repeat([0], data.t), ρ_h, M, n_outer_cuts, types)        
            for t in T
                expr = get_expr(mp, sp_duals_PK[t], data, ρ_h, M, t, types)
                if sp_stats_PK[t] == MOI.OPTIMAL
                    status.nOptCuts += 1                
                    opt_cut = @constraint(mp, mp[:α][t] >= expr) 
                    mp[Symbol("PK_opt_cut_$(t)_$(status.nIter)")] = opt_cut
                end
            end
        end


        end_stat_mp, of_mp, yvals_mp, αvals_mp = solve_benders_mp(mp)
        if of_mp + tol >= ub_mp
            println("Stop since mp of is > tan ub")
            break
        end
        
        if end_stat_mp == MOI.INFEASIBLE
            println("The problem with aggregation $agg is infeasible ****")
            conv = true
            break
        end

        ## Solve SPs
        if "SH" in types || "Rel" in types
            sp_vars, sp_stats, sp_feas, sp_ofs, sp_duals = solve_benders_sp_primal(sp_m, data, params, status, solver, yvals_mp, ρ_h, M, n_outer_cuts, μ, update_ρ, agg)
        elseif "B" in types || "PK" in types || "FC" in types
            sp_stats, sp_feas, sp_ofs, sp_duals = solve_benders_sp_dual(sp_m, data, params, status, solver, yvals_mp, αvals_mp, ρ_h, M, n_outer_cuts, types, update_ρ, w_fc, w0)
        elseif "MW" in types
            # This is the dual that provides feas cuts
            sp_stats, sp_feas, sp_ofs, sp_duals = solve_benders_sp_dual(sp_m, data, params, status, solver, yvals_mp, αvals_mp, ρ_h, M, n_outer_cuts, ["B"], update_ρ, w_fc, w0)
        # elseif "Rel" in types
        #     sp_vars, sp_stats, sp_feas, sp_ofs, sp_duals = solve_benders_sp_primal(sp_m, data, params, status, solver, yvals_mp, ρ_h, M, n_outer_cuts, μ, false, agg)
        end

        # sp_vars, sp_stats, sp_feas, sp_ofs, sp_duals = solve_benders_sp_primal(prims, data, params, status, solver, yvals_mp, ρ_h, M, n_outer_cuts, μ, false, agg)              
        
        # println("αvals_mp = $(sum(αvals_mp))")
        # println("sum of sp =$(sum(sp_ofs))")
        if sum(αvals_mp) + tol >= sum(sp_ofs)
            conv = true
            break
        end

        lb_mp = of_mp
        last_bounds = append!(last_bounds, lb_mp)[2:end]
        conv = last_bounds[end]-last_bounds[begin] <= τ*last_bounds[end] ? true : false

        cuts_types = "all"
        # If the lb is greatter than the ub, then only add feasibility cuts
        # if lb_mp > ub + tol
        #     cuts_types = "feas"
        # end

        if "Rel" in types
            _cuts_gen, cuts_rel_iter, _cuts_sep_other = separate_cuts(mp, αvals_mp, ρ_h, M, data, status, ["B"], false, cuts_types, sp_stats, sp_feas, sp_ofs, sp_duals, [], agg)
            cuts_rel += cuts_rel_iter
        else
            _cuts_gen, _cuts_rel, cuts_sep_other = separate_cuts(mp, αvals_mp, ρ_h, M, data, status, types, false, cuts_types, sp_stats, sp_feas, sp_ofs, sp_duals)
        end

        # For MW, use the y values for the first iteration,
        # then calculate a new point between feasible solutions
        # If y is not feasible, then use the same last feasible solution
        if "MW" in types && !(false in sp_feas)
            if !first_it_MW
                interior_y = τ_MW.*interior_y .+ (1-τ_MW).*yvals_mp
            elseif first_it_MW
                if !(false in sp_feas)
                    interior_y = yvals_mp
                else
                    interior_y = gen_y(data, repeat([data.k], data.J))
                end
                first_it_MW = false
            end
            sp_stats_MW, sp_feas_MW, sp_ofs_MW, sp_duals_MW = solve_benders_sp_dual(dualsMW, data, params, status, solver, interior_y, αvals_mp, ρ_h, M, n_outer_cuts, types, update_ρ, w_fc, w0, yvals_mp, sp_ofs)
            _cuts_gen, _cuts_rel, cuts_sep_other = separate_cuts(mp, αvals_mp, ρ_h, M, data, status, ["B"], false, cuts_types, sp_stats_MW, sp_feas_MW, sp_ofs_MW, sp_duals_MW)
        end

        update_ρ = false
        update_M = false

        if !("Rel" in types)
            ρvals = Array{Float64}(undef, data.J, data.t)    
            Rvals = Array{Float64}(undef, data.J, data.t)
            for t in T
                if "B" in types || "PK" in types || "FC" in types
                    ρvals[:,t] = -dual.(sp_m[t][:cρ_sp_d])
                    Rvals[:,t] = -dual.(sp_m[t][:cR_sp_d])
                elseif "MW" in types
                    ρvals[:,t] = -dual.(sp_m[t][:cρ_sp_d]) #dualsMW
                    Rvals[:,t] = -dual.(sp_m[t][:cR_sp_d]) 
                elseif "SH" in types
                    ρvals[:,t] = sp_vars[t]["ρ"].data
                    Rvals[:,t] = sp_vars[t]["R"].data
                end
            end

            # Thigten the outer approximation only if the sp are all feasible

            if  !(false in sp_feas) # && 2 ==3
                # for j in J, t in T
                #     # Update if violation
                #     if Rvals[j,t] + tol < ρvals[j,t]/(1-ρvals[j,t])
                for I in findall(==(1), Rvals .< (ρvals ./ (1 .- ρvals)))
                    update_ρ = true
                    j, t = I[1], I[2]
                    if ρvals[j,t] > maximum(ρ_h[j,t,:]) && !update_M
                        update_M = true
                    end
                    # end
                    # Avoid having Nan in the M calc. This actually happens
                    # if ρvals > 1 - 10e-18, but we decrease the max to limit the
                    # max value of M
                    if ρvals[j,t] > 1 - 10e-4
                        ρvals[j,t] = 1 - 10e-4
                    end

                end        
                # println("We are uptading p =$update_ρ")
                if update_ρ
                    ρ_h = cat(ρ_h, ρvals, dims=3)
                    if update_M    
                        M = calc_big_M(data, ρ_h)
                    end
                end
            end
        end
        
        ub_mp_temp = dot(data.F, yvals_mp) + sum(sp_ofs)
        
        # Update the UB
        if ub_mp_temp < ub_mp
            ub_mp =  ub_mp_temp
            println("LB=$lb_mp iter = $(status.nIter)")      
            println("UB=$ub_mp iter = $(status.nIter)\n")
        end

        lb_iter[status.nIter] = lb_mp
        ub_iter[status.nIter] = ub_mp

        if "SH" in types || "Rel" in types
            for t in T
                vars = sp_vars[t]
                xval[:,:,t] = vars["x"]
                ρval[:,t] = vars["ρ"]
                Rval[:,t] = vars["R"]
                wval[:,:,t] = vars["w"]
                zval[:,:,t] = vars["z"]
            end
        else                
            for t in T
                xval[:,:,t] = -dual.(sp_m[t][:cx_sp_d])
                ρval[:,t] = -dual.(sp_m[t][:cρ_sp_d])
                Rval[:,t] = -dual.(sp_m[t][:cR_sp_d])
                wval[:,:,t] = -dual.(sp_m[t][:cw_sp_d])
                zval[:,:,t] = -dual.(sp_m[t][:cz_sp_d])
            end
        end  
    end
    
    return lb_mp, ub_mp, lb_iter, ub_iter, iter_rel, cuts_rel, yvals_mp, αvals_mp, xval, ρval, Rval, wval, zval, ρ_h #, yint, ub_yint
end