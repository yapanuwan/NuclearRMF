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

p = true
N = 82
r_max = maximum(xs)
E_min = 880
E_max = 939

κs, Es = findAllOrbitals(p, S_interp, V_interp, R_interp, A_interp, r_max, E_min, E_max)
occ = fillNucleons(N, κs, Es)

scatter(κs, Es, label=(p ? "proton" : "neutron") * " spectrum")
annotate!(κs .+ 0.3, Es, occ)
xlabel!("κ")
ylabel!("E (MeV)")