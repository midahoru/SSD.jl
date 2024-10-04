# SSD
Service system design with service time variability and time varying demand

## Basic usage

### Instance creation
```julia
# Import the package
using SSD

# Create default parameters
params = SSD.default_params()

# Specify the attributes defining the instances
coords_bounds=(0,100) # Bounds for the coords of the nodes
r_bounds=(0.3, 0.4) # Bounds of r
k=5 # Number of capacity levels
t=3 # Time periods
a_bounds=(80,120) # Bounds for λ   

nodes_sizes = [(50,5),(100,10), (150,15), (200,20), (250,25)] #(customers, facilities)
cvs = [0.5, 1, 1.5] # Options for cv
Dt = [1, 10, 25, 50, 100] # Options for D
FLR = [0.4, 0.6, 0.8] # Facility load ratio
FCR = [2,4,8,10] # Facility capacity ratio
for (nI, nJ) in nodes_sizes
    for cv in cvs
        for D in Dt
            for flr in FLR
                for fcr in FCR
                    SSD.instance_gen(nI, nJ, coords_bounds, a_bounds, r_bounds, cv, D, k, t, flr, fcr, params)
                end
            end
        end
    end
end
```

### Solve the exact problem with a non-linear constraint
```julia
# Import the package
using SSD
using CPUTime

# Create default parameters
params = SSD.default_params()
# Create the data container
data = SSD.default_data()




# Read the instance to solve
directory = "instances/Own/"
filename= "1.txt"
SSD.read_file(directory*filename, data)

# Define the different cost levels
cost_levels = [0.60, 0.85, 1, 1.15, 1.35]
data.F = SSD.gen_costs(data, params, cost_levels)
# Define the different capacity levels
cap_levels = [0.5, 0.75, 1, 1.25, 1.5]
data.Q = SSD.gen_caps(data, params, cap_levels)

status = SSD.init_solver_status()

solve_method = "bendersPK" # or "lazy_cuts, benders, bendersMW

el = @CPUelapsed res = SSD.solve_ssd(data, params, status, solve_method)
if solve_method == "iter_cuts"
    println("$solve_method: lb=$(res[1]), ub=$(res[2]),
        of_term1_lb=$(res[3]), of_term2_lb=$(res[4]), of_term3_lb=$(res[5]),
        of_term1_ub=$(res[6]), of_term2_ub=$(res[7]), of_term3_ub=$(res[8]),
        y=$(res[9]), CPU_t=$el")
elseif solve_method in ["benders", "bendersWS1", "bendersWS2", "bendersMW", "bendersPK", "bendersSH", "bendersFC", "bendersIter","bendersIterPK"]
    println("$solve_method: of=$(res[2]), of_term1=$(res[3]), of_term2=$(res[4]), of_term3=$(res[5]), y=$(res[6]), CPU_t=$el, nodes=$(res[9]),
        feas=$(res[10]), opt=$(res[11])")
else
    println("$solve_method: of=$(res[2]), of_term1=$(res[3]), of_term2=$(res[4]), of_term3=$(res[5]), y=$(res[6]), CPU_t=$el")
end
```


### Solve the problem with the linear approximation
```julia
# Import the package
using SSD

# Create default parameters
params = SSD.default_params()
# Create the data container
data = SSD.default_data()


# Read the instance to solve
filename = "instances/I_50 J_5 (0, 100) cv_0.5 D_1 k_5 t_3 FLR_0.4 FCR_2 a_(80, 120) r_(0.3, 0.4).txt"
SSD.read_file(filename, data)

# Define the different cost levels
cost_levels = [0.60, 0.85, 1, 1.15, 1.35]
data.F = SSD.gen_costs(data, params, cost_levels)
# Define the different capacity levels
cap_levels = [0.5, 0.75, 1, 1.25, 1.5]
data.Q = SSD.gen_caps(data, params, cap_levels)

# Create the set ρ_h with values between 0 and 1
# without the extreme points of the domain     
ρ_h = SSD.ini_ρ_h(data)

status = SSD.init_solver_status()

lb, ub, x, y = SSD.cutting_plane(data, params, status, ρ_h, ϵ)

# Gets the number of iterations
println("N. iter = $(status.nIter)")








```
