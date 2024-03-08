module SSD

    using Random 
    using Dates
    using Gurobi
    using JuMP
    using MathOptInterface
    using Distributions
    using Distances

    include("typedefs.jl")
    include("iter.jl")
    include("util.jl")

end # module SSD
