mutable struct Parameters
    max_time::Float64    # In seconds
    ϵ::Float64
    round_digits::Int64
    rng::MersenneTwister
    rng2::MersenneTwister
end

mutable struct Data
    I::Int64
    J::Int64
    Icoords::Matrix{Float64}
    Jcoords::Matrix{Float64}
    dist::Matrix{Float64}
    a::Matrix{Float64}
    F::Matrix{Float64}
    Q::Matrix{Float64}
    C::Array{Float64}
    cv::Float64
    D::Int64    
    k::Int64
    t::Int64
    FLR::Float64
    FCR::Float64
    M::Int64
end

mutable struct FeasibleSolution
    coords::Matrix{Float64}
    caps::Vector{Float64}
end


mutable struct SolverStatus
    initTime::DateTime
    endTime::DateTime
    endStatus::Symbol
    nIter::Int64
    nFeasCuts::Int64
    nOptCuts::Int64
end

default_params() = Parameters(3*60*60, 10^-5, 4, MersenneTwister(0), MersenneTwister(15))
default_data() = Data(0, 0, Array{Float64,}(undef,0,0), Array{Float64}(undef,0,0), #(0,1),
Array{Float64}(undef,0, 0), Array{Float64}(undef,0, 0),
Array{Float64}(undef,0, 0), Array{Float64}(undef,0, 0), Array{Float64}(undef,0,0,0),
1, 1, 1, 1, 1, 1, 1000)
init_solver_status() = SolverStatus(Dates.now(), Dates.now(), :none, 0, 0, 0)