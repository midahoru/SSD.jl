mutable struct Parameters
    max_time::Float64    
    ϵ::Float64
    round_digits::Int64
    rng::MersenneTwister
end

mutable struct Data
    name::String
    Icoords::Matrix{Float64}
    Jcoords::Matrix{Float64}    
    dist::Matrix{Float}
    λ::Vector{Float64}
    C::Matrix{Float64}
    cv::Float64
    D::Int64    
    k::Int64
    t::Int64
end

mutable struct FeasibleSolution
    coords::Matrix{Float64}
    caps::Vector{Float64}
end


mutable struct SolverStatus
    initTime::DateTime
    endTime::DateTime
    ok::Bool # if false, optimization has been aborted due to time limits
    endStatus::Symbol
end

default_params() = Parameters(30*60, 10^-5, 4, MersenneTwister(0))
init_solver_status() = SolverStatus(Dates.now(), Dates.now(), true, :none)