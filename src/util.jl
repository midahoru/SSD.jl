################################################## Models

function gen_caps(data, params, cap_levels)
    # Gets the max demand per time period
    max_dem_t = maximum(sum.(eachcol(data.a)))
    Q_j3 = zeros(Float64, data.J)
    for j in 1:data.J
        Q_j3[j] = 1.25*max_dem_t / (data.J * data.FLR)
    end
    return round.(Q_j3 .* cap_levels', digits=params.round_digits)
end

function gen_costs(data, params, cost_levels)
    delta_x  = maximum(maximum.([data.Jcoords[1,:], data.Icoords[1,:]]))-minimum(minimum.([data.Jcoords[1,:], data.Icoords[1,:]]))
    delta_y = maximum(maximum.([data.Jcoords[2,:], data.Icoords[2,:]]))-minimum(minimum.([data.Jcoords[2,:], data.Icoords[2,:]]))
    x_cen = delta_x/2 + minimum(minimum.([data.Jcoords[1,:], data.Icoords[1,:]]))
    y_cen = delta_y/2 + minimum(minimum.([data.Jcoords[2,:], data.Icoords[2,:]]))
    f = (x,y) -> euclidean((x,y), (x_cen, y_cen))
    # f = (x,y) -> euclidean((x,y), (50,50))
    F_j3 = zeros(Float64, data.J)
    for j in 1:data.J
        F_j3[j] = data.FCR*f(data.Jcoords[1,j], data.Jcoords[2,j])
    end
    return round.(F_j3 .* cost_levels', digits=params.round_digits)
end

function solve_ssd(data, params, status, solve_method)
    if solve_method == "nlp"
        of, of_term1, of_term2, of_term3, y, x = minlp(data, params, status)
        return round(of, digits=params.round_digits), round(of, digits=params.round_digits), round(of_term1, digits=params.round_digits), round(of_term2, digits=params.round_digits), round(of_term3, digits=params.round_digits), convert_y_to_print(y.data, data), x.data, status.endStatus

    elseif solve_method == "iter_cuts"
        lb, ub, of_term1_lb, of_term2_lb, of_term3_lb, of_term1_ub, of_term2_ub, of_term3_ub, y, x = cuts_priori(data, params, status)
        return round(lb, digits=params.round_digits), round(ub, digits=params.round_digits), round(of_term1_lb, digits=params.round_digits), round(of_term2_lb, digits=params.round_digits), round(of_term3_lb, digits=params.round_digits), round(of_term1_ub, digits=params.round_digits), round(of_term2_ub, digits=params.round_digits), round(of_term3_ub, digits=params.round_digits), convert_y_to_print(y.data, data), x.data, status.endStatus

    
    
    
    elseif solve_method == "lazy_cuts"
        of, of_term1, of_term2, of_term3, y, x, test_cap, test_cap_sel, test_alloc = model_lazy_cuts(data, params, status)
        return round(of, digits=params.round_digits), round(of, digits=params.round_digits), round(of_term1, digits=params.round_digits), round(of_term2, digits=params.round_digits), round(of_term3, digits=params.round_digits), convert_y_to_print(y.data, data), x.data, status.endStatus, test_cap, test_cap_sel, test_alloc, [] #convert_y_to_print(y.data, data), x.data

    elseif solve_method == "benders"
        of, of2, of_term1, of_term2, of_term3, y, x, nodes, f_cut, opt_cut, lb, ub, test_cap, test_cap_sel, test_alloc, relax_iters  = model_benders(data, params, status)
        return round(of, digits=params.round_digits), round(of2, digits=params.round_digits), round(of_term1, digits=params.round_digits), round(of_term2, digits=params.round_digits), round(of_term3, digits=params.round_digits), convert_y_to_print(y, data), x, status.endStatus, nodes, f_cut, opt_cut, lb, ub, test_cap, test_cap_sel, test_alloc, relax_iters
        # return of, of, of_term1, of_term2, of_term3, y, x, nodes, f_cut, opt_cut      

    elseif solve_method == "bendersMW"
        of, of2, of_term1, of_term2, of_term3, y, x, nodes, f_cut, opt_cut, lb, ub, test_cap, test_cap_sel, test_alloc, relax_iters  = model_benders(data, params, status, ["MW"])
        return round(of, digits=params.round_digits), round(of2, digits=params.round_digits), round(of_term1, digits=params.round_digits), round(of_term2, digits=params.round_digits), round(of_term3, digits=params.round_digits), convert_y_to_print(y, data), x, status.endStatus, nodes, f_cut, opt_cut, lb, ub, test_cap, test_cap_sel, test_alloc, relax_iters

    elseif solve_method == "bendersPK"
        of, of2, of_term1, of_term2, of_term3, y, x, nodes, f_cut, opt_cut, lb, ub, test_cap, test_cap_sel, test_alloc, relax_iters  = model_benders(data, params, status, ["PK"])
        return round(of, digits=params.round_digits), round(of2, digits=params.round_digits), round(of_term1, digits=params.round_digits), round(of_term2, digits=params.round_digits), round(of_term3, digits=params.round_digits), convert_y_to_print(y, data), x, status.endStatus, nodes, f_cut, opt_cut, lb, ub, test_cap, test_cap_sel, test_alloc, relax_iters

    
    
    
    
    
    elseif solve_method == "bendersSH"
        of, of2, of_term1, of_term2, of_term3, y, x, nodes, f_cut, opt_cut, lb, ub  = model_benders(data, params, status, ["SH"])
        return round(of, digits=params.round_digits), round(of2, digits=params.round_digits), round(of_term1, digits=params.round_digits), round(of_term2, digits=params.round_digits), round(of_term3, digits=params.round_digits), convert_y_to_print(y.data, data), x, status.endStatus, nodes, f_cut, opt_cut, lb, ub

    elseif solve_method == "bendersFC"
        of, of2, of_term1, of_term2, of_term3, y, x, nodes, f_cut, opt_cut, lb, ub  = model_benders(data, params, status, ["FC"])
        return round(of, digits=params.round_digits), round(of2, digits=params.round_digits), round(of_term1, digits=params.round_digits), round(of_term2, digits=params.round_digits), round(of_term3, digits=params.round_digits), convert_y_to_print(y.data, data), x, status.endStatus, nodes, f_cut, opt_cut, lb, ub

    elseif solve_method == "bendersIter"
        of, of2, of_term1, of_term2, of_term3, y, x, nodes, f_cut, opt_cut, lb, ub, iter_lp, iter_nodes  = model_benders_iter(data, params, status)
        return round(of, digits=params.round_digits), round(of2, digits=params.round_digits), round(of_term1, digits=params.round_digits), round(of_term2, digits=params.round_digits), round(of_term3, digits=params.round_digits), convert_y_to_print(y.data, data), x, status.endStatus, nodes, f_cut, opt_cut, lb, ub, iter_lp, iter_nodes

    elseif solve_method == "heur_NM"
        of, Scost, Ccost, Congcost, best_y, best_x = heur_nelder_mead(data, params, status, ini_y(data), Dict(), 1000, 50)    
        return round(of, digits=params.round_digits), round(of, digits=params.round_digits), round(Scost, digits=params.round_digits), round(Ccost, digits=params.round_digits), round(Congcost, digits=params.round_digits), best_y, best_x, status.endStatus

    elseif solve_method == "heur_LS"
        return heur_local_search(data, params, status, 1000, 50)
        
    elseif solve_method == "heur_FB"
        of, of_term1, of_term2, of_term3, best_y, best_x = heur_local_search_first_best(data, params, status, 1000, 300)
        return round(of, digits=params.round_digits), round(of, digits=params.round_digits), round(of_term1, digits=params.round_digits), round(of_term2, digits=params.round_digits), round(of_term3, digits=params.round_digits), best_y, best_x, status.endStatus      
    end
end
    




################################################## Linear model with cuts
function ini_ρ_h(data, n=9)
    J = 1:data.J 
    T = 1:data.t

    ini_ρ_h = collect(range(0.1,step=1/(n+1),0.9))
    H = 1:length(ini_ρ_h)
    ρ_h = Array{Float64}(undef,data.J,data.t,length(ini_ρ_h))

    for j in J, t in T, h in H
        ρ_h[j, t, h]= ini_ρ_h[h]
    end
    return ρ_h
end


function calc_ub(x, y, data)
    I = 1:data.I
    J = 1:data.J
    λ = data.a
    C = data.C
    F = data.F
    cv = data.cv
    D = data.D    
    K = 1:data.k
    T = 1:data.t

    ρ = calc_new_ρ(x, y, data)
    R = calc_new_R(ρ, data)
    w, z = calc_new_w_z(x, y, R, ρ, data)

    Fterm = dot(F, y)
    Cterm = dot(C, x)

    Dt = [D/sum(λ[i, t] for i in I) for t in T]
    
    last_term = 0.5 * sum(Dt[t] * (R[j, t] + ρ[j, t] + sum(cv^2 * (w[j, k, t] - z[j, k, t]) for k in K)) for j in J for t in T)

    ub = Fterm + Cterm + last_term
    return ub, Fterm, Cterm, last_term
end

function calc_ub(lb, ρq, Rq, wq, data)
    Dt = [data.D/sum(data.a[i, t] for i in 1:data.I) for t in 1:data.t]

    new_ub = lb + sum(Dt[t]*(((ρq[j, t]/(1-ρq[j,t]))-Rq[j,t]) + data.cv^2*((ρq[j, t]/(1-ρq[j,t]))-
    sum(wq[j,k,t] for k in 1:data.k))) for j in 1:data.J for t in 1:data.t)
    return new_ub
end

function calc_new_ρ(xq, yq, data)
    I = 1:data.I
    J = 1:data.J 
    λ = data.a
    T = 1:data.t 
    K = 1:data.k
    Q = data.Q    
    T = 1:data.t

    ρ_new = zeros(data.J, data.t)
    for j in J, t in T
        num = sum(λ[i,t]*xq[i,j,t] for i in I)
        den = sum(Q[j,k]*yq[j,k] for k in K)
        if den != 0
            ρ_new[j,t]=num/den
        else
            ρ_new[j,t]=num/rand(1)
        end
    end
    return ρ_new
end

function calc_new_R(ρ, data)
    R = zeros(data.J, data.t)
    for j in 1:data.J, t in 1:data.t
        R[j, t] = ρ[j, t] / (1 - ρ[j, t])
    end
    return R
end

function calc_new_w_z(x, y, R, ρ, data)
    z = zeros(data.J, data.k, data.t)
    w = zeros(data.J, data.k, data.t)
    for j in 1:data.J, k in 1:data.k, t in 1:data.t
        z[j, k, t] = ρ[j, t] * y[j, k]
        w[j, k, t] = R[j, t] * y[j, k]
    end
    return w, z
end

function calc_big_M(data, ρ_h)
    J = 1:data.J
    T = 1:data.t
    H = 1:size(ρ_h,3)

    M = zeros(data.J, data.t)
    for j in J, t in T
        max_in_1 = 0
        for h in H
            value_in_1 = 1/(1-ρ_h[j,t,h])^2 -ρ_h[j,t,h]^2/(1-ρ_h[j,t,h])^2
            if value_in_1 > max_in_1
                max_in_1 = value_in_1
            end
        end
        M[j,t]=max_in_1
    end
    return M
end

function calc_big_M_2(data, ρ_h)
    J = 1:data.J
    T = 1:data.t
    M = zeros(data.J, data.t)
    for j in J, t in T
        max_ρ = maximum(ρ_h[j,t,:])
        M[j,t]=1+ 1/(1-max_ρ)^2 -max_ρ^2/(1-max_ρ)^2
    end
    return M
end 

################################################## Nelder-Mead
# Cost of the Subproblem
function calc_cost_sp(y0, data, ρ_h, primals, μ, n_outter_cuts, gen_yb = false)
    ρ_k = Array{Float64}(undef,data.J,data.t)
    M = calc_big_M(data, ρ_h) 

    x = zeros(data.I, data.J, data.t)

    y = gen_yb == true ? gen_y(data, y0) : y0

    # Check total capacity of 1st stage variables
    if dot(data.Q, y) <= maximum(sum(data.a,dims=1))
        non_opt_val = 10e5*(sum(data.F)+sum(data.C)+data.D*sum(data.a))
        return false, non_opt_val, non_opt_val, non_opt_val, ρ_k, x    
    else
        # Location cost
        Loccost = dot(data.F, y)
        # Assignation cost
        Allocost = 0
        # Congestion cost
        Congcost = 0

        Allocost, Congcost, xval_sp, ρval_sp, Rval_sp, wval_sp, zval_sp, all_sp_stat, all_sp_feas, all_sp_vals, all_sp_duals = solve_benders_sp_primal(primals, data, y, ρ_h, n_outter_cuts, μ)
        
        for t in 1:data.t

            # sp_stat, sp_xval, sp_ρval, sp_Rval, sp_wval, sp_zval, Alloterm, Congterm, π1, π2, π4, π6, π8, π10, π11, π12 = solve_benders_sp_primal(primals[t], data, ρ_h, t, n_outter_cuts)
            sp_stat = all_sp_stat[t]
            
            if sp_stat == MOI.OPTIMAL
                ρ_k[:,t] = ρval_sp[t]  
                x[:,:,t] = xval_sp[t]       
            end
            # # Add costs
            # Allocost += Alloterm
            # Congcost += Congterm
        end
        return true, Loccost, Allocost, Congcost, ρ_k, x
    end
end

function search_y_same_sol(y_ind, x_k, data, all_sols)
    y = gen_y(data, y_ind)
    # Sol is feasible?
    y_feas = sum(data.Q[j,k]*y[j,k] for j in 1:data.J for k in 1:data.k) >= maximum(sum(data.a,dims=1))

    # Set-up cost
    of_sp = sum(data.F .* y)

    dem_j = [[sum(data.a[i,t]*x_k[i,j,t] for i in 1:data.I)  for t in 1:data.t] for j in 1:data.J]
    max_dems = []
    for j in 1:data.J
        append!(max_dems, maximum(dem_j[j]))
    end

    # y_best = deepcopy(y_ind)
    min_cap_y_ind = deepcopy(y_ind)

    for ind in eachindex(y_ind)
        keep = true
        cap_down = 1
        # Search equivalent solutions with less cap
        while keep && y_ind[ind] - cap_down >= 0

            cap_ref = y_ind[ind] - cap_down == 0 ? 0 : data.Q[ind, y_ind[ind]-cap_down]
            
            if cap_ref >= max_dems[ind]
                new_y_ind = deepcopy(y_ind)
                new_y_ind[ind] -= cap_down
                if !(new_y_ind in keys(all_sols))
                    # println("@@@@@@@@@@@@   New sol found : $new_y_ind and is equivalent to $y_ind")
                    new_y = gen_y(data, new_y_ind)
                    new_cost = all_sols[y_ind]-of_sp+sum(data.F .* new_y)
                    all_sols[new_y_ind] = new_cost
                    # if new_cost < all_sols[y_best]
                    #     y_best = new_y_ind
                    # end
                end                
                cap_down += 1
            else
                min_cap_y_ind[ind] = y_ind[ind] - cap_down + 1
                keep = false
            end
        end
        # Search equivalent solutions with greatter cap
        # if y_feas && y_ind[ind] > 0 && data.Q[ind, y_ind[ind]] > max_dems[ind] && y_ind[ind] < data.k
        #     keep_up = true
        #     cap_up = 1
        #     while keep_up && y_ind[ind] + cap_up <= data.k
        #         new_y_ind = deepcopy(y_ind)
        #         new_y_ind[ind] += cap_up
        #         if !(new_y_ind in keys(all_sols))
        #             println("@@@@@@@@@@@@####   New sol found : $new_y_ind and is equivalent to $y_ind")
        #             new_y = gen_y(data, new_y_ind)
        #             new_cost = all_sols[y_ind]-of_sp+sum(data.F .* new_y)
        #             all_sols[new_y_ind] = new_cost
        #             # if new_cost < all_sols[y_best]
        #             #     println(" !!!!!!!!!!!!!!!!!!! A new sol with more cap is better $new_y_ind !!!!!!!!!!!!!!!")
        #             #     y_best = new_y_ind                        
        #             # end
        #         end                
        #         cap_up += 1
        #     end
        # end
    end
    min_cap_y = gen_y(data, min_cap_y_ind)
    all_sols[min_cap_y_ind] = all_sols[y_ind]-of_sp+sum(data.F .* min_cap_y)
    println("Best equivalent y is $min_cap_y_ind")
    return all_sols, min_cap_y_ind, all_sols[min_cap_y_ind] #y_best
end

function modify_y(y_ind, x, data, all_sols, params, status, op_down, op_up, α1 = 60, α2 = 10)
    # Select index based on workload            
    dem_j = [[sum(data.a[i,t]*x[i,j,t] for i in 1:data.I)  for t in 1:data.t] for j in 1:data.J]
    max_dems = []
    for j in 1:data.J
        append!(max_dems, maximum(dem_j[j]))
    end

    y = gen_y(data, y_ind)
    cap_y = sum(data.Q.*y, dims=2)    
    cap_rem = vec((cap_y .- max_dems)*100 ./cap_y)
    replace!(cap_rem, NaN=>0)

    ind_down = [ind for ind in eachindex(cap_rem) if cap_rem[ind] > α1]
    ind_up = [ind for ind in eachindex(cap_rem) if cap_rem[ind] <= α2]

    # all_ops = [(j,op_down) for j in ind_down]
    all_ops = [(j,0) for j in ind_down]
    append!(all_ops, [(j,op_up) for j in ind_up]) 

    y_ret = deepcopy(y_ind)

    if size(all_ops,1) > 0
        for (j,op) in all_ops
            y_ret[j] = op == 0 ? 0 : 0 <= y_ret[j] + op <= data.k ? y_ret[j] + op : y_ret[j]
        end
        mat_y_ret = gen_y(data, y_ret)
        if sum(data.Q .* mat_y_ret) >= maximum(sum(data.a,dims=1))
            return y_ret
        end
    end

    # y_ret = perturb_y(y_ind, data, all_sols, op_down, op_up)
    y_ini = ini_y(data)
    y_ini[end]=y_ind
    of_ret, y_ret = heur_nelder_mead_y(data, params, status, y_ini, all_sols, 50, 10)
    println("Vector perturbed using NM from $y_ind to $y_ret")
    
    
    return y_ret
end

function perturb_y(y_ind, data, params, status, all_sols, op_down, op_up)
    # all_ops = [(j,o) for j in 1:data.J for o in [op_down, op_up]]
    # shuffle!(all_ops)

    # for (j,op) in all_ops
    #     y_ret = deepcopy(y_ind)
    #     y_ret[j] = 0 <= y_ret[j] + op <= data.k ? y_ret[j] + op : y_ret[j]
    #     if !(y_ret in keys(all_sols))
    #         return y_ret
    #     end
    # end

    # y_ret = deepcopy(y_ind)
    # for (j,op) in all_ops
    #     y_ret[j] = 0 <= y_ret[j] + op <= data.k ? y_ret[j] + op : y_ret[j]
    #     if !(y_ret in keys(all_sols))
    #         return y_ret
    #     end
    # end
    y_test = deepcopy(y_ind)

    while minimum(y_test) < data.k
        for ind in findall(<(data.k),y_test)
            y_test[ind] = y_test[ind] + 1 <= data.k ? y_test[ind] + 1 : y_test[ind]

            if !(y_test in keys(all_sols))
                y_ini = ini_y(data)
                y_ini[end]=y_test
                of_ret, _Scost, _Ccost, _Congcost, y_ret, x_ret = heur_nelder_mead(data, params, status, y_ini, all_sols, 20, 5)
                if !(y_ret in keys(all_sols))
                    println("Vector perturbed using NM from $y_test to $y_ret")
                    return y_ret
                end
            end
        end
    end
    

    println("Random perturbation")
    iter = 0
    while iter < 10
        y_ret = rand(0:data.k, data.J)
        if !(y_ret in keys(all_sols))
            return y_ret
        end
        iter += 1
    end
    return rand(0:data.k, data.J) 
end


# Initial vertices of the simplex
function ini_y(data)
    vecs = []
    for i in 1:data.J
        y_ind = Int.(floor(data.k/2) .+ zeros(data.J))
        y_temp = gen_y(data, y_ind)
        # If the initial solution is infeasible or just to tight
        # increase all caps by one
        if sum(data.Q .* y_temp) <= maximum(sum(data.a,dims=1))
            y_ind = y_ind .+ 1
        end
        # if i > 1
        #     y_ind[i] = Int(data.k)
        # end

        y_ind[i] = data.k-1
        push!(vecs,y_ind)
    end
    return push!(vecs,[data.k for i in 1:data.J])
    
end

# Creates an array for the y variable from a vector of
# capacity index
function gen_y(data, y_ind)
    y = zeros(data.J, data.k)
    for j_ind in 1:size(y_ind)[1]
        k_ind = y_ind[j_ind]
        if k_ind > 0
            y[j_ind, k_ind] = 1
        end
    end
    return y
end

function round_vertex(y_raw, data)
    # return Int.(round.([min(max(0,i), data.k) for i in y_raw]))
    return Int.(round.([min(max(0,i), data.k) for i in y_raw]))
end

function calc_centr(all_sols, data, y_ix, id_keep)
    # return Int.(round.([sum(i[j] for i in y_ix[id_keep]) for j in 1:length(y_ix[id_keep])] ./ length(y_ix[id_keep])))
    # return round.([sum(i[j] for i in y_ix[id_keep]) for j in 1:length(y_ix[id_keep])] ./ length(y_ix[id_keep]))
    y_cent_raw = [sum(i[j] for i in y_ix[id_keep]) for j in 1:length(y_ix[id_keep])] ./ length(y_ix[id_keep])
    return calc_new_y(all_sols, y_cent_raw, y_ix, data, "general")
    # return y_cent_raw
end

function calc_ref(all_sols, data, y_cent, y_not, α)
    y_ref_raw =  y_cent + α .*(y_cent - y_not)
    # y_ref = Int.(round.([min(max(0,i), data.k) for i in y_ref_raw]))
    return calc_new_y(all_sols, y_ref_raw, y_not, data, "general")
end

function calc_exp(all_sols, data, y_cent, y_ref, γ)
    y_exp_raw = y_cent + γ .*(y_ref - y_cent)
    # y_exp = Int.(round.([min(max(0,i), data.k) for i in y_exp_raw]))
    return calc_new_y(all_sols, y_exp_raw, y_ref, data, "general")
end

function calc_con(all_sols, data, y_cent, y_reference, ρNM)
    y_con_raw = y_cent + ρNM .*(y_reference - y_cent)
    # y_con = Int.(round.([min(max(0,i), data.k) for i in y_con_raw]))
    return calc_new_y(all_sols, y_con_raw, y_reference, data, "general")
end

function calc_shr(data, all_sols, y_old, y_best, σ)
    y_shr_raw = y_best + σ .*(y_old - y_best)
    # y_shr = Int.(round.([min(max(0,i), data.k) for i in y_shr_raw]))
   # println("Shrinked $y_shr from $y_old using $y_old and the best is $y_best")
    return calc_new_y(all_sols, y_shr_raw, y_old, data, "shrink")
end

function calc_new_y(all_sols, y_new_raw, y_old, data, type="general")
    y_new_raw = round_vertex(y_new_raw, data)
    return y_new_raw
    y_new = deepcopy(y_new_raw)
    previous_vert = collect(keys(all_sols))

    if y_new in previous_vert# y_ix
        # if type == "general"
        max_var =  max.(0, y_old .- y_new)
        keep = true
        while keep
            # If no more increase is possible, start decreasing
            # if (sum(max_var) == 0 || type == "shrink") && sum(y_new) >= 1            
                # rand_ind = rand(1:length(y_new))
                # if y_new[rand_ind] > 0
                #     y_new[rand_ind] -= 1
                # end
            if sum(max_var) == 0  && sum(y_new) >= 1
                rand_ind = rand(findall(>(0), y_new))
                y_new[rand_ind] -= 1

            # Increase if possible    
            elseif sum(max_var) > 0 #&& type == "general"
                var_ind = findall(>(0), max_var)
                rand_ind = rand(var_ind)
                max_var[rand_ind] -= 1
                y_new[rand_ind] += 1

            elseif sum(max_var) == 0  && sum(y_new) == 1
                keep = false
            end

            if !(y_new in previous_vert) || sum(y_new) == 0
                return y_new   
            end
        end
    else
        return y_new
    end   
end


function heur_sp(y, data, ρ_h, t, GRB_ENV)
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
    m = Model(() -> Gurobi.Optimizer(GRB_ENV))
    # Set the Gurobi parameters
    set_optimizer_attribute(m, "OutputFlag", 0)
    set_optimizer_attribute(m, "Threads", 1)
    set_optimizer_attribute(m, "MIPGap", 1e-5)
    set_optimizer_attribute(m, "LazyConstraints", 1)
    set_optimizer_attribute(m, "Presolve", 0)
                                    
    @variable(m, 0 <= x[I,J])
    @variable(m, 0 <= z[J,K])
    @variable(m, 0 <= ρ[J])

    @variable(m, 0 <= w[J,K])
    @variable(m, 0 <= R[J])

    @expression(m, of_term2, sum(C[i,j,t]*x[i,j] for i in I for j in J))
    @expression(m, of_term3, 0.5*sum(Dt*(R[j] + ρ[j] + sum(cv^2*(w[j,k]-z[j,k]) for k in K)) for j in J))

    @objective(m, Min, of_term2 + of_term3)
    
    # @objective(m, Min, sum(C[i,j,t]*x[i,j] for i in I for j in J) +
    # 0.5*sum(Dt*(R[j] + ρ[j] + sum(cv^2*(w[j,k]-z[j,k]) for k in K)) for j in J))

    # Capacity cannot be exceeded and steady state has to be conserved
    c1_sp = @constraint(m, [j in J], -sum(λ[i,t]*x[i,j] for i in I) >= -sum(Q[j,k]*y[j,k] for k in K), base_name = "c1_sp")
    
    # All customer zones need to be assigned to exactly one facility
    c2_1_sp = @constraint(m, [i in I], sum(x[i,j] for j in J) >= 1, base_name = "c2_1_sp")
    c2_2_sp = @constraint(m, [i in I], -sum(x[i,j] for j in J) >= -1, base_name = "c2_2_sp")
    
    @constraint(m, [j in J], sum(λ[i,t]*x[i,j] for i in I) - sum(Q[j,k]*z[j,k] for k in K) == 0, base_name = "c3_sp")

    c4_sp = @constraint(m, [j in J, k in K], -z[j,k] >= -y[j,k], base_name = "c4_sp")

    @constraint(m, [j in J], sum(z[j,k] for k in K) - ρ[j] == 0, base_name = "c5_sp")

    c6_sp = @constraint(m, [j in J, k in K], -w[j,k] >= -M[j,t]*y[j,k], base_name = "c6_sp")
    # c6_sp = @constraint(m, [j in J, k in K], w[j,k] <= M*y[j,k], base_name = "c6_sp")

    @constraint(m, [j in J], sum(w[j,k] for k in K) - R[j] == 0, base_name = "c7_sp")

    c8_sp = @constraint(m, [j in J, h in H], (1-ρ_h[j,t,h])^2 * R[j] - ρ[j] >= -ρ_h[j,t,h]^2, base_name = "c8_sp")

    # c9_sp = @constraint(m, [j in J, k in K], y[j,k] == y_fixed[j,k], base_name = "c9_sp")
    
    c10_sp = @constraint(m, [i in I, j in J], -x[i,j] >= -1, base_name = "c10_sp")
    c11_sp = @constraint(m, [j in J, k in K], -z[j,k] >= -1, base_name = "c11_sp")
    c12_sp = @constraint(m, [j in J], -ρ[j] >= -(1-10e-5), base_name = "c12_sp")

    optimize!(m)
    non_opt_val = 10e5*(sum(data.F)+sum(data.C)+data.D*sum(data.a))
    end_stat = termination_status(m)
    if end_stat == MOI.DUAL_INFEASIBLE || end_stat == MOI.INFEASIBLE
        # Feasibility cut
        return end_stat, non_opt_val, non_opt_val, non_opt_val, [], []
    elseif end_stat == MOI.OPTIMAL        
        ρval = value.(ρ)
        xval = value.(x)
        return end_stat, objective_value(m), value(of_term2), value(of_term3), ρval, xval
    else
        println("other results $end_stat")
        return end_stat, non_opt_val, non_opt_val, non_opt_val, [], []       
    end       
end

################################################## General
function is_sol_feas(data, y, x)
    test_cap = true
    test_cap_sel = true
    test_alloc = true
    for t in 1:data.t
        test_cap = sum(round.(x[:,:,t]'data.a[:,t], digits=3) .> sum(y.*data.Q,dims=2)) > 0 ? false : true
        for i in 1:data.I
            test_alloc = sum(x[i,j,t] for j in 1:data.J) == 1 ? true : false
        end
    end
    for j in 1:data.J
        test_cap_sel = sum(y[j,k] for k in 1:data.k) == 1 ? true : sum(y[j,k] for k in 1:data.k) == 0 ? true : false
    end
    return test_cap, test_cap_sel, test_alloc
end


################################################## miscellaneous

function convert_y_to_print(y, data)
    y_cap = []
    for j in 1:data.J
        # Use in case of continuous y
        # append!(y_cap,0)
        if sum(y[j,:]) > 0.1
            append!(y_cap,findfirst(>=(0.9), y[j,:]))
        else
            append!(y_cap,0)
        end
    end
    return y_cap
end

function set_maximum_time(params, s)
    params.max_time = s
end

function elapsed(status)
    b = Dates.now()
    ms = Dates.value(b - status.initTime)
    s = ms / 1000
end

function total_elapsed(status)
    ms = Dates.value(status.endTime - status.initTime)
    s = ms / 1000
end

function optimal(status)
    return status.endStatus == :optimal
end

function get_status(status)
    return status.endStatus
end

function get_number_groups(status)
    return status.endGroupsNb
end

function set_data(data::Data, name, length, heigth)
    data.name = name
    data.len = length
    data.hei = hei
    return data
end

function set_params(params::Parameters, max_time, rounding_digits, seed)
    params.max_time = max_time
    params.round_digits = rounding_digits
    params.rng = MersenneTwister(seed)
    return params
end

##################################################


