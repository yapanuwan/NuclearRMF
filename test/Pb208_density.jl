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

N_p = 82
N_n = 126
r_max = maximum(xs)
E_min = 880
E_max = 939
divs = 400

rs = range(0, r_max, length=divs+1)

(ρ_sp, ρ_vp) = calculateNucleonDensity(N_p, true, S_interp, V_interp, R_interp, A_interp, r_max, divs, E_min, E_max)

p_sp = plot(rs, ρ_sp, xlabel="r (fm)", label="ρₛₚ(r) calculated")
p_vp = plot(rs, ρ_vp, xlabel="r (fm)", label="ρᵥₚ(r) calculated")

(ρ_sn, ρ_vn) = calculateNucleonDensity(N_n, false, S_interp, V_interp, R_interp, A_interp, r_max, divs, E_min, E_max)

p_sn = plot(rs, ρ_sn, xlabel="r (fm)", label="ρₛₙ(r) calculated")
p_vn = plot(rs, ρ_vn, xlabel="r (fm)", label="ρᵥₙ(r) calculated")

# benchmark data generated from Hartree.f
# format: x  Rhos(n)  Rhov(n)  Rhot(n)  Rhos(p)  Rhov(p)  Rhot(p)
bench_data = readdlm("test/Pb208DensFSUGarnet.csv")
xs = bench_data[:, 1]
rho_sn = bench_data[:, 2]
rho_vn = bench_data[:, 3]
rho_sp = bench_data[:, 5]
rho_vp = bench_data[:, 6]

plot!(p_sp, xs, rho_sp, label="ρₛₚ(r) benchmark")
plot!(p_vp, xs, rho_vp, label="ρᵥₚ(r) benchmark")
plot!(p_sn, xs, rho_sn, label="ρₛₙ(r) benchmark")
plot!(p_vn, xs, rho_vn, label="ρᵥₙ(r) benchmark")

plot(p_sp, p_vp, p_sn, p_vn, layout=(2, 2))
