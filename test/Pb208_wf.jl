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
r_max = maximum(xs)
E_min = 880
E_max = 939

approxE = findEs(κ, p, S_interp, V_interp, R_interp, A_interp, r_max, E_min, E_max) |> minimum
groundE = refineEs(κ, p, S_interp, V_interp, R_interp, A_interp, r_max, [approxE])[1]

println("ground state E = $groundE")

divs = 50
wf = solveWf(κ, p, groundE, S_interp, V_interp, R_interp, A_interp, r_max, divs)
rs = range(0, r_max, length=divs+1)
gs = wf[1, :]
fs = wf[2, :]

plot(rs, gs, label="g(r)")
plot!(rs, fs, label="f(r)")
xlabel!("r (fm)")
