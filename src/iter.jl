# function cuts_priori(data, params, status)    

#     # Initialize the bounds
#     lb, ub = 0, sum(data.F)+sum(data.C)+data.D*sum(data.a)

#     xq, yq, zq, ρq, wq, Rq = [], [], [], [], [], []
#     ρ_h = ini_ρ_h(data)
#     q = 0
#     ub_dict = Dict()
#     of_term1_lb = 0
#     of_term2_lb = 0
#     of_term3_lb = 0

#     while (ub-lb)/ub >= params.ϵ
#         new_lb, of_term1_lb, of_term2_lb, of_term3_lb, yq, xq, zq, ρq, wq, Rq = model_iter_cuts(data, params, status, ρ_h)
#         if (new_lb-ub)/ub > params.ϵ 
#             break
#         end
#         if status.endStatus == :infeasible || status.endStatus == :none
#             break
#         end
#         lb = new_lb
#         newub, of_term1_ub, of_term2_ub, of_term3_ub = calc_ub(xq, yq, data) #calc_ub(lb, ρq, Rq, wq, data)#
#         ub_dict[newub] =[of_term1_ub, of_term2_ub, of_term3_ub]
#         ub = min(ub, newub)
#         ρ_h = cat(ρ_h, ρq, dims=3)
#         q+=1
#         # println("Iter $q: LB= $lb ; UB = $ub")
#     end
#     status.nIter = q

#     # return xq, yq, zq, ρq, wq, Rq, lb, ub
#     return lb, ub, of_term1_lb, of_term2_lb, of_term3_lb, ub_dict[ub][1], ub_dict[ub][2], ub_dict[ub][3], yq, xq
# end

function cuts_priori(data, params, status)

    # Initialize the bounds
    lb, ub = 0, sum(data.F)+sum(data.C)+data.D*sum(data.a)

    xq, yq, zq, ρq, wq, Rq = [], [], [], [], [], []
    ρ_h = ini_ρ_h(data)
    q = 0
    ub_dict = Dict()
    of_term1_lb = 0
    of_term2_lb = 0
    of_term3_lb = 0

    J = 1:data.J
    T = 1:data.t
    K = 1:data.k

    m_iter = ini_model_iter_cuts(data, params, status, ρ_h)

    while (ub-lb)/ub >= params.ϵ
        new_lb, of_term1_lb, of_term2_lb, of_term3_lb, yq, xq, zq, ρq, wq, Rq = solve_model_iter_cuts(m_iter, status)

        if (new_lb-ub)/ub > params.ϵ 
            break
        end
        if status.endStatus == :infeasible || status.endStatus == :none
            break
        end
        lb = new_lb
        newub, of_term1_ub, of_term2_ub, of_term3_ub = calc_ub(xq, yq, data)
        ub_dict[newub] =[of_term1_ub, of_term2_ub, of_term3_ub]
        ub = min(ub, newub)

        # Add the new constraints for the linear outter approximation
        @constraint(m_iter, [j in J, t in T], (1-ρq[j,t])^2 * m_iter[:R][j,t] - m_iter[:ρ][j,t] >= -ρq[j,t]^2)
        
        # Modify the constraint with the big M
        ρ_h = cat(ρ_h, ρq, dims=3)
        M = calc_big_M(data, ρ_h)
        for j in J, t in T, k in K
            set_normalized_coefficient(m_iter[:cM][j,t,k], m_iter[:y][j,k], -M[j,t])
        end
        q+=1
    end
    status.nIter = q

    return lb, ub, of_term1_lb, of_term2_lb, of_term3_lb, ub_dict[ub][1], ub_dict[ub][2], ub_dict[ub][3], yq, xq
end
