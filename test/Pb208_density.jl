using DelimitedFiles, Interpolations, Plots
include("../dirac.jl")

# test data generated from Hartree.f
# format: x S(x) V(x) R(x) A(x)
test_data = readdlm("test/Pb208FldsFSUGarnet.csv") 
xs = test_data[:, 1]
Ss = test_data[:, 2]
Vs = test_data[:, 3]
Rs = test_data[:, 4]
As = test_data[:, 5]

S_interp = linear_interpolation(xs, Ss)
V_interp = linear_interpolation(xs, Vs)
R_interp = linear_interpolation(xs, Rs)
A_interp = linear_interpolation(xs, As)

κ = -1
p = true
N = 82
r_max = maximum(xs)
E_min = 880
E_max = 939
divs = 400

(ρ_s, ρ_v) = calculateNucleonDensity(N, p, S_interp, V_interp, R_interp, A_interp, r_max, divs, E_min, E_max)

rs = range(0, r_max, length=divs+1)

plot(rs, ρ_s, label="ρₛ(r)")
plot!(rs, ρ_v, label="ρᵥ(r)")
xlabel!("r (fm)")
