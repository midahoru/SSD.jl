# SSD
Service system design with service time variability and time varying demand

## Basic usage

### Instance creation
```julia
# Import the package
using SSD

# Specify the attributes defining the instances
coords_bounds=(0,100)
r_bounds=(0.3, 0.4)
k=5
t=3
λ_bounds=(80,120)

nI=250
nJ=25

for cv in [0.5, 1, 1.5]
    for D in [1, 10, 25, 50, 100]
        SSD.instance_gen(nI, nJ, coords_bounds, λ_bounds, r_bounds, cv, D, k, t, params0)
    end
end
