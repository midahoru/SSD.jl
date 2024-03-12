# SSD
Service system design with service time variability and time varying demand

## Basic usage

### Instance creation
```julia
# Import the package
using SSD

# Specify the attributes defining the instances
coords_bounds=(0,100) # Bounds for the coords of the nodes
r_bounds=(0.3, 0.4) # Bounds of r
k=5 # Number of capacity levels
t=3 # Time periods
λ_bounds=(80,120) # Bounds for λ   

nI=250 # Number of customers zones
nJ=25 # Number of potential facilities locations

cvs = [0.5, 1, 1.5] # Options for cv
Ds = [1, 10, 25, 50, 100] # Options for D
for cv in cvs
    for D in Ds
        SSD.instance_gen(nI, nJ, coords_bounds, λ_bounds, r_bounds, cv, D, k, t, params0)
    end
end
```

### Solving the exact non-linear problem
```julia
# Import the package
using SSD

# Create the data container
data = default_data()


# Read the instance to solve
filename = "I_50 J_5 (0, 100) cv_0.5 D_1 k_5 t_3 lam_(80, 120) r_(0.3, 0.4).txt"
SSD.read_instance()
