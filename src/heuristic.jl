function heur_nelder_mead(primals, μ,  data, params, status, y_ix, GRB_ENV, all_sols=Dict(), max_iter = 100, max_iter_no_impro = 10)
    eps = 1e-2                   # tolerance for stopping criterion
    α = 1  #0.1 #                      # reflection scale (>0)
    γ = 2  #1.5 #                      # expansion scale (>1)
    ρNM = 1/2 #0.1 #                     # contraction scale (>0  and <=0.5)
    σ = 1/2  # 0.9 #                   # shrinking scale (<1)

    timer = Timer(params.max_time - elapsed(status))
    while true

        isopen(timer) || break
        yield()

        # Initial data
        # y_ix = ini_y(data)
        ρ_h = ini_ρ_h(data)

        sol_vert = Dict()
        # all_sols = Dict()
        # all_sols[[0 for i in 1:data.J]] = 10e5*(sum(data.F)+sum(data.C)+data.D*sum(data.a))

        
        int_y_ind = [data.k for j in 1:data.J]
        int_y = gen_y(data, int_y_ind)
        # μ = 0
        # primals = Dict()
        # for t in 1:data.t
        #     primals[t] = ini_benders_sp_primal(int_y, data, ρ_h, t, μ, params, status, GRB_ENV)
        # end 

        index_to_look_for = 1:size(y_ix)[1]

        iter = 0
        iter_no_imp = 0

        while iter <= max_iter && iter_no_imp <= max_iter_no_impro
            iter += 1
            println()
            println("----------------Starting iter $iter --------------------------")       
            for ind in index_to_look_for
                y_ind = y_ix[ind]
                if y_ind in keys(all_sols)
                    println("----- Retrieving sol for $y_ind -----")
                    of_sp = all_sols[y_ind]
                else
                    println("----- Solving SP for $y_ind -----")
                    res_calc, Scost, Ccost, Congcost, ρ_k, x_k = calc_cost_sp(y_ind, data, ρ_h, primals, μ, true)
                    of_sp = Scost + Ccost + Congcost
                    all_sols[y_ind] = of_sp
                    # all_sols, y_ind, of_sp = search_y_same_sol(y_ind, x_k, data, all_sols)
                    if res_calc
                        ρ_h = cat(ρ_h, ρ_k, dims=3)
                    end
                end                           
                sol_vert[ind] = of_sp            
            end

            # println("Sol vert $sol_vert")

            # Indexes to keep and to compare against
            s_order = sort(collect(sol_vert), by = tuple -> last(tuple), rev=false)
            keys_sorted_sols_vert = [s[1] for s in s_order]     
            id_keep = collect(keys_sorted_sols_vert)[1:end-1] 
            id_not = collect(keys_sorted_sols_vert)[end]
            index_to_look_for = [id_not]
            # println("y_ix = $y_ix")
            # println("id_keep $id_keep") 
            
            # Calculates the centroid
            y_cent = calc_centr(all_sols, data, y_ix, id_keep)
            # println("centroide $y_cent")
            # println("Centroid $y_cent")

            # Calculates the reflection
            y_ref = calc_ref(all_sols, data, y_cent, y_ix[id_not], α)
            # println("Reflection $y_ref")

            # Calculates the solution for the reflection
            if y_ref in keys(all_sols)
                of_sp_ref = all_sols[y_ref]
            else
                res_calc, Scost, Ccost, Congcost, ρ_k, x_k = calc_cost_sp(y_ref, data, ρ_h, primals, μ, true)
                of_sp_ref = Scost + Ccost + Congcost                
                all_sols[y_ref] = of_sp_ref
                if res_calc
                    ρ_h = cat(ρ_h, ρ_k, dims=3)
                end
                # all_sols, y_ref, of_sp_ref = search_y_same_sol(y_ref, x_k, data, all_sols)
            end 

            # Cases
            ## Reflection: If reflection is better than the 2nd worst, but not the best:
            # if minimum(values(sol_vert)) <= of_sp_ref &&  of_sp_ref < sort(collect(values(sol_vert)))[end-1]
            if sol_vert[keys_sorted_sols_vert[1]] <= of_sp_ref < sol_vert[keys_sorted_sols_vert[end-1]]
                println("####### Reflection #######")
                # y_ref_new = calc_new_y(all_sols, y_ref, y_ix[id_not])
                y_ix[id_not] = y_ref
                iter_no_imp = 0

            ## Expansion: If the reflection point is the best one
            elseif of_sp_ref < sol_vert[keys_sorted_sols_vert[1]]#minimum(values(sol_vert))
                # Calculates the expanded point
                y_exp = calc_exp(all_sols, data, y_cent, y_ref, γ)
                
                # Calculates the solution for the expanded point
                if y_exp in keys(all_sols)
                    of_sp_exp = all_sols[y_exp]
                else
                    res_calc, Scost, Ccost, Congcost, ρ_k, x_k = calc_cost_sp(y_exp, data, ρ_h, primals, μ, true)
                    of_sp_exp = Scost + Ccost + Congcost
                    all_sols[y_exp] = of_sp_exp
                    if res_calc
                        ρ_h = cat(ρ_h, ρ_k, dims=3)
                    end
                    # all_sols, y_exp, of_sp_exp = search_y_same_sol(y_exp, x_k, data, all_sols)
                end
                
                # Select the new vertex comparing expanded and reflected
                if of_sp_exp < of_sp_ref
                    println("####### Expansion #######")
                    # y_ix[id_not] = calc_new_y(all_sols, y_exp, y_ix[id_not])
                    y_ix[id_not] = y_exp
                    iter_no_imp = 0
                else
                    println("####### Reflection 2 #######")
                    # y_ix[id_not] = calc_new_y(all_sols, y_ref, y_ix[id_not])
                    y_ix[id_not] = y_ref
                    iter_no_imp = 0
                end

            ## Contraction (if we got here, then reflected is equal or worst than the 2nd worst): 
            else
                if of_sp_ref < sol_vert[id_not] # Reflection is best than worst vertex
                    y_reference = y_ref
                else #if of_sp_ref >= maximum(values(sol_vert)) # Reflection is worst than worst vertex           
                    y_reference = y_ix[id_not]
                end

                # Calculates contracted point
                y_con = calc_con(all_sols, data, y_cent, y_reference, ρNM)

                # Calculates the solution for the contracted point
                if y_con in keys(all_sols)
                    of_sp_con = all_sols[y_con]
                else
                    res_calc, Scost, Ccost, Congcost, ρ_k, x_k = calc_cost_sp(y_con, data, ρ_h, primals, μ, true)
                    of_sp_con = Scost + Ccost + Congcost
                    all_sols[y_con] = of_sp_con
                    if res_calc
                        ρ_h = cat(ρ_h, ρ_k, dims=3)
                    end                    
                    # all_sols, y_con, of_sp_con = search_y_same_sol(y_con, x_k, data, all_sols)
                end

                # Select the new vertex comparing contracted and reflected
                if of_sp_con < all_sols[y_reference] # Since we can be sure it has been already calculated
                    println("####### Contraction #######")
                    y_ix[id_not] = y_con
                    iter_no_imp = 0
                    # y_ix[id_not] = calc_new_y(all_sols, y_con, y_ix[id_not])                
                # Shrink
                else
                    index_to_look_for = []
                    ind_best = keys_sorted_sols_vert[1]
                    println("####### Shrink #######")                    
                    iter_no_imp += 1
                    for ind in 1:size(y_ix)[1]
                        if ind != ind_best
                            new_y = calc_shr(data, all_sols, y_ix[ind], y_ix[ind_best], σ)
                            if y_ix[ind] != new_y
                                y_ix[ind] = new_y
                            else
                                # max_rand = Int(floor(data.k/2))
                                max_rand = 1
                                rand_ad_k = rand(-max_rand:max_rand, data.J)
                                # y_replace_shrink = [k+r_k <= data.k ? k+r_k : k+r_k-1 <= data.k ? k+r_k-1 : k for (k,r_k) in zip(y_ix[ind_best],rand_ad_k)]
                                y_replace_shrink = [0 <= k+r_k <= data.k ? k+r_k : k for (k,r_k) in zip(y_ix[ind_best],rand_ad_k)]
                                println("Shrinked the same, change to $y_replace_shrink instead of $new_y")
                                y_ix[ind] = y_replace_shrink
                            end
                            push!(index_to_look_for, ind)                  
                        end          
                    end
                end
            end
            
        end
        # best_sol = sort(collect(all_sols), by = tuple -> last(tuple), rev=false)[1]
        best_argmin = argmin(all_sols)
        # res_calc, Scost, Ccost, Congcost, ρ_k, x_k = calc_cost_sp(best_sol[1], data, ρ_h, primals, μ, GRB_ENV)
        # return  best_sol[2], Scost, Ccost, Congcost, best_sol[1], x_k
        res_calc, Scost, Ccost, Congcost, ρ_k, x_k = calc_cost_sp(best_argmin, data, ρ_h, primals, μ, true)
        return  all_sols[best_argmin], Scost, Ccost, Congcost, best_argmin, x_k, y_ix
    
    end
    
    best_argmin = argmin(all_sols)
    res_calc, Scost, Ccost, Congcost, ρ_k, x_k = calc_cost_sp(best_argmin, data, ρ_h, primals, μ, true)
    return  all_sols[best_argmin], Scost, Ccost, Congcost, best_argmin, x_k
end

function heur_local_search(data, params, status, max_iter = 100, max_iter_no_impro = 10,
    op_up = 1, op_down = -1, α1 = 70, α2 = 10)

    timer = Timer(params.max_time - elapsed(status))

    while true
        isopen(timer) || break
        yield()

        
        # indexes = collect(1:data.J)
        ρ_h = ini_ρ_h(data)
        # GRB_ENV = Gurobi.Env()
        
        all_sols = Dict()

        best_y_ind = [data.k for j in 1:data.J]
        int_y = gen_y(data, best_y_ind)
        μ = 0
        primals = Dict()
        for t in 1:data.t
            primals[t] = ini_benders_sp_primal(int_y, data, ρ_h, t, μ, params, status)
        end 


        res_calc, Scost, Ccost, Congcost, ρ_k, best_x_k = calc_cost_sp(best_y_ind, data, ρ_h, primals, μ, true)
        best_cost = Scost + Ccost + Congcost
        all_sols[best_y_ind] = best_cost
        # all_sols, best_y_ind, best_cost = search_y_same_sol(best_y_ind, best_x_k, data, all_sols)
        if res_calc
            ρ_h = cat(ρ_h, ρ_k, dims=3)
        end

        test_y_ind = modify_y(best_y_ind, best_x_k, data, all_sols, params, status, op_down, op_up, α1, α2)
        println("Y modified from $best_y_ind to $test_y_ind")

        iter = 1
        iter_no_imp = 1

        while iter < max_iter  && iter_no_imp < max_iter_no_impro
            println("")
            println("Starting iter $iter")

            if !(test_y_ind in keys(all_sols))
                res_calc, Scost, Ccost, Congcost, ρ_k, x_k = calc_cost_sp(test_y_ind, data, ρ_h, primals, μ, true)
                new_cost = Scost + Ccost + Congcost
                all_sols[test_y_ind] = new_cost
                # all_sols, test_y_ind, new_cost = search_y_same_sol(test_y_ind, x_k, data, all_sols)
                ρ_h = cat(ρ_h, ρ_k, dims=3)

                if new_cost < best_cost
                    println("The min so far is $test_y_ind with OF $new_cost")
                    best_cost = new_cost
                    best_y_ind = test_y_ind
                    best_x_k = x_k
                    iter_no_imp = 0
                    test_y_ind = modify_y(best_y_ind, best_x_k, data, all_sols, params, status, op_down, op_up, α1, α2)
                    println("Y modified from $best_y_ind to $test_y_ind")
                else
                    iter_no_imp += 1
                    # test_y_ind = perturb_y(best_y_ind, data, all_sols, op_down, op_up)                    
                    println("new sol not better. Y perturbed from $test_y_ind")
                    test_y_ind = perturb_y(test_y_ind, data, params, status, all_sols, op_down, op_up)
                    println("to $test_y_ind") 
                end
            else
                iter_no_imp += 1
                # test_y_ind = perturb_y(best_y_ind, data, all_sols, op_down, op_up)
                println("no sol not in all sols. Y perturbed from $test_y_ind")
                test_y_ind = perturb_y(test_y_ind, data, params, status, all_sols, op_down, op_up)
                println(" to $test_y_ind")
            end                         
            iter += 1
        end
        status.nIter = iter
        return  best_cost, best_y_ind, best_x_k, argmin(all_sols)
        
    end 
    status.nIter = iter
    return  best_cost, best_y_ind, best_x_k, argmin(all_sols)   
    
end



function heur_local_search_first_best(data, params, status, max_iter = 100, max_iter_no_impro = 10, op_up = 1, op_down = -1)
    
    timer = Timer(params.max_time - elapsed(status))

    while true
        isopen(timer) || break
        yield()

        y_ind = repeat([data.k], data.J)
        ρ_h = ini_ρ_h(data)
        GRB_ENV = Gurobi.Env()

        # int_y_ind = [data.k for j in 1:data.J]
        int_y = gen_y(data, y_ind)
        μ = 0
        primals = Dict()
        for t in 1:data.t
            primals[t] = ini_benders_sp_primal(int_y, data, ρ_h, t, μ, params, status)
        end 

        iter = 0
        iter_no_imp = 0

        all_sols = Dict()
        all_sols_detailed = Dict()

        res_calc, Scost, Ccost, Congcost, ρ_k, x_ind = calc_cost_sp(y_ind, data, ρ_h, primals, μ, true)
        best_min = Scost + Ccost + Congcost
        all_sols[y_ind]=best_min
        all_sols_detailed[y_ind]  = [Scost, Ccost, Congcost, x_ind] 
        # all_sols, y_ind = search_y_same_sol(y_ind, x_k, data, all_sols)
        if res_calc
            ρ_h = cat(ρ_h, ρ_k, dims=3)
        end

        all_ops = [(j,op) for j in 1:data.J for op in [op_down, op_up]]

        improved = true

        while iter < max_iter && iter_no_imp < max_iter_no_impro && improved
            improved = false                      

            shuffle!(all_ops)
            y_test = deepcopy(y_ind)
            x_k = deepcopy(x_ind)
            ind_op = 1
            last_op = (-1,-1)
            
            while !improved && ind_op <= length(all_ops)
                println("Starting iter $iter")

                (j,op) = all_ops[ind_op]
                # Test feasibility
                println("Test operator $op on fac $j and base $y_test")
                if (j,op) != last_op && 0 <= y_test[j] + op <= data.k   
                    y_test[j] += op
                    println("Test $y_test")
                    if !(y_test in keys(all_sols))
                        res_calc, Scost, Ccost, Congcost, ρ_k, x_k = calc_cost_sp(y_test, data, ρ_h, primals, μ, true)
                        new_cost = Scost + Ccost + Congcost
                        all_sols[y_test] = new_cost
                        all_sols_detailed[y_test]  = [Scost, Ccost, Congcost, x_k]                       
                        # all_sols, y_test = search_y_same_sol(y_test, x_k, data, all_sols)
                        if res_calc
                            ρ_h = cat(ρ_h, ρ_k, dims=3)
                        end

                        if new_cost < best_min
                            best_min = new_cost
                            y_ind = y_test
                            x_ind = x_k
                            improved = true
                            iter_no_imp = 0
                            last_op = (j,-op)
                            
                        else
                            y_test[j] -= op 
                        end
                    else
                        y_test[j] -= op
                    end
                    # else
                    #     new_cost = all_sols[y_test]
                    # end
                    # if new_cost < best_min
                    #     best_min = new_cost
                    #     y_ind = y_test
                    #     improved = true
                    #     iter_no_imp = 0
                    # else
                    #     iter_no_imp += 1
                    # end 
                end

                if !improved
                    iter_no_imp += 1
                    max_rand = 2
                    rand_ad_k = rand(-max_rand:max_rand, data.J)
                    y_ind = [0 <= k+r_k <= data.k ? k+r_k : k for (k,r_k) in zip(y_ind,rand_ad_k)]
                end
                ind_op += 1
                iter += 1 
            end                           
            # iter += 1
        end
        status.nIter = iter
        best_argmin = argmin(all_sols)
        return  all_sols[best_argmin], all_sols_detailed[best_argmin][1], all_sols_detailed[best_argmin][2], all_sols_detailed[best_argmin][3], best_argmin, all_sols_detailed[best_argmin][4]
    end
    status.nIter = iter
    best_argmin = argmin(all_sols)
    return  all_sols[best_argmin], all_sols_detailed[best_argmin][1], all_sols_detailed[best_argmin][2], all_sols_detailed[best_argmin][3], best_argmin, all_sols_detailed[best_argmin][4] 
end