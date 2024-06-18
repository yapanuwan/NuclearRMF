using DifferentialEquations, Roots

ħc = 197.327 # ħc in MeVfm
M_n = 939.5654133 # Neutron mass in MeV/c2
M_p = 938.2720813 # Proton mass in MeV/c2

"The spherical Dirac equation that returns du=[dg, df] in-place where
(g, f) are the reduced radial components evaluated at r,
κ is the generalized angular momentum,
M is the mass in MeV/c2,
E in the energy in MeV,
S(r) & V(r) are functions corresponding to scalar and vector potentials in MeV,
r is the radius in fm.
Reference: P. Giuliani, K. Godbey, E. Bonilla, F. Viens, and J. Piekarewicz, Frontiers in Physics 10, (2023)."
function dirac!(du, (g, f), (κ, M, E, S, V), r)
    du[1] = -(κ/r) * g + (E + M - S(r) - V(r)) * f / ħc
    du[2] =  (κ/r) * f - (E - M + S(r) - V(r)) * g / ħc
end

"Solve the Dirac equation and return g(r=r_max) for given scalar and vector potentials where
κ is the generalized angular momentum,
M is the mass in MeV/c2,
E in the energy in MeV,
S(r) & V(r) are functions corresponding to scalar and vector potentials in MeV,
r_max is the outer boundary in fm,
r_min (=r_max/1000) is inside boundary in fm which cannot be 0 due to the centrifugal term."
function boundaryValue(κ, M, E, S, V, r_max, r_min=r_max/1000)
    prob = ODEProblem(dirac!, [0, 1], (r_min, r_max))
    sol = solve(prob, RK4(), p=(κ, M, E, S, V))
    return sol(r_max)[1]
end

"Find all bound energies between E_min (=0) and E_max (=M) where
κ is the generalized angular momentum,
M is the mass in MeV/c2,
S(r) & V(r) are functions corresponding to scalar and vector potentials in MeV,
r_max is the outer boundary in fm,
r_min (=r_max/1000) is inside boundary in fm which cannot be 0 due to the centrifugal term."
function findEs(κ, M, S, V, r_max, r_min=r_max/1000, E_min=0, E_max=M)
    f(E) = boundaryValue(κ, M, E, S, V, r_max, r_min)
    return find_zeros(f, (E_min, E_max))
end
