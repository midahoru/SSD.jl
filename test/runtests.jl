using SSD

# Create default parameters
params = SSD.default_params()
# Create the data container
data = SSD.default_data()


# Read the instance to solve
filename = "instances/I_250 J_25 (0, 100) cv_0.5 D_1 k_5 t_3 FLR_0.4 FCR_2 a_(80, 120) r_(0.3, 0.4).txt"
# filename = "I_50 J_5 (0, 100) cv_0.5 D_1 k_5 t_3 FLR_0.4 FCR_2 a_(80, 120) r_(0.3, 0.4).txt"
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

x, y, cost = SSD.cutting_plane(data, params, status, ρ_h)
