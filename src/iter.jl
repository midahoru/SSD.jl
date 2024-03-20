function cutting_plane(data, params, status, ρ_h, ϵ)
    # Initialize the bounds
    lb, ub = 0, sum(data.F)+sum(data.C)+data.D*sum(data.λ)

    xq, yq, lb = [], [], 0
    q = 0
    while (ub-lb)/lb >= ϵ || q < 20
        xq, yq, lb = model_cuts(data, params, status, ρ_h)
        ub = calc_ub(ub, xq, yq, data)

        ρ_new = calc_new_ρ(xq, yq, data)

        cat(ρ_h, ρ_new, dims=3)
        q+=1
    end

    return xq, yq, lb    
end
