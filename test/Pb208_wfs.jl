using DelimitedFiles, Interpolations, Plots
include("../mesons.jl")

# test data generated from Hartree.f
# format: x  Rhos(n)  Rhov(n)  Rhot(n)  Rhos(p)  Rhov(p)  Rhot(p)
test_data = readdlm("test/Pb208DensFSUGarnet.csv")
xs = test_data[:, 1]
rho_sn = test_data[:, 2]
rho_vn = test_data[:, 3]
rho_sp = test_data[:, 5]
rho_vp = test_data[:, 6]

ρ_sn = linear_interpolation(xs, rho_sn)
ρ_vn = linear_interpolation(xs, rho_vn)
ρ_sp = linear_interpolation(xs, rho_sp)
ρ_vp = linear_interpolation(xs, rho_vp)

r_max = maximum(xs)
divs = 400

wfs = solveWfs(ρ_sp, ρ_vp, ρ_sn, ρ_vn, r_max, divs)

rs = range(0, r_max, length=divs+1)

plot(rs, transpose(wfs))
xlabel!("r (fm)")
