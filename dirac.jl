using DifferentialEquations, Roots

const ħc = 197.327 # ħc in MeVfm
const M_n = 939.5654133 # Neutron mass in MeV/c2
const M_p = 938.2720813 # Proton mass in MeV/c2

const r_reg = 1E-6 # regulator for the centrifugal term in fm

"The spherical Dirac equation that returns du=[dg, df] in-place where
    u=[g, f] are the reduced radial components evaluated at r,
    κ is the generalized angular momentum,
    p is true for proton and false for neutron,
    E in the energy in MeV,
    Φ0, W0, B0, A0 are the mean-field potentials (couplings included) in MeV as functions of r in fm,
    r is the radius in fm.
Reference: P. Giuliani, K. Godbey, E. Bonilla, F. Viens, and J. Piekarewicz, Frontiers in Physics 10, (2023)"
function dirac!(du, u, (κ, p, E, Φ0, W0, B0, A0), r)
    M = p ? M_p : M_n
    common1 = E - W0(r) - (p - 0.5) * B0(r) - p * A0(r)
    common2 = M - Φ0(r)
    (g, f) = u
    du[1] = -(κ/(r + r_reg)) * g + (common1 + common2) * f / ħc
    du[2] =  (κ/(r + r_reg)) * f - (common1 - common2) * g / ħc
end

"Solve the Dirac equation and return g(r=r_max) for given scalar and vector potentials where
    r_max is the outer boundary in fm,
    the other parameters are the same from dirac!(...)."
function boundaryValue(κ, p, E, Φ0, W0, B0, A0, r_max)
    prob = ODEProblem(dirac!, [0, 1], (0, r_max))
    sol = solve(prob, RK4(), p=(κ, p, E, Φ0, W0, B0, A0))
    return sol(r_max)[1]
end

"Find all bound energies between E_min (=0) and E_max (=mass) where
    the other parameters are the same from dirac!(...)."
function findEs(κ, p, Φ0, W0, B0, A0, r_max, E_min=0, E_max=(p ? M_p : M_n))
    f(E) = boundaryValue(κ, p, E, Φ0, W0, B0, A0, r_max)
    return find_zeros(f, (E_min, E_max))
end
