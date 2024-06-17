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

Es = collect(840:0.5:940)
boundaryVals = [boundaryValue(-1, M_n, E, S_interp, V_interp, maximum(xs))^2 for E in Es]

plot(Es, boundaryVals, yscale=:log10, label="g(r_max)^2")
xlabel!("E (MeV)")
