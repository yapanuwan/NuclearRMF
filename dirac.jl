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

"Solve the Dirac equation and return the wave function u(r)=[g(r), f(r)] where
    divs is the number of mesh divisions if the solution should be discretized as a 2×(1+divs) matrix (keep divs=0 to obtain an interpolating function),
    refine determines whether to switch to high-precision mode and optimize the energy beforehand (assuming a bound state),
    the other parameters are the same from dirac!(...)."
function solveWf(κ, p, E, Φ0, W0, B0, A0, r_max, divs::Int=0, refine=true)
    if refine
        dtype = BigFloat
        algo = Feagin12()

        f(E_) = boundaryValue(κ, p, E_, Φ0, W0, B0, A0, r_max; dtype=dtype, algo=algo)
        E = find_zero(f, convert(dtype, E))
    else
        dtype = Float64
        algo = RK()
    end

    prob = ODEProblem(dirac!, convert.(dtype, [0, 1]), (0, r_max))
    sol = solve(prob, algo, p=(κ, p, E, Φ0, W0, B0, A0), saveat=(divs == 0 ? [] : r_max/divs))
    return divs == 0 ? sol : hcat(sol.u...)
end

"Solve the Dirac equation and return g(r=r_max) where
    r_max is the outer boundary in fm,
    the other parameters are the same from dirac!(...)."
function boundaryValue(κ, p, E, Φ0, W0, B0, A0, r_max; dtype=Float64, algo=RK4())
    prob = ODEProblem(dirac!, convert.(dtype, [0, 1]), (0, r_max))
    sol = solve(prob, algo, p=(κ, p, E, Φ0, W0, B0, A0), saveat=[r_max], save_idxs=[1])
    return sol[1, 1]
end

"Find all bound energies between E_min (=0) and E_max (=mass) where
    the other parameters are the same from dirac!(...)."
function findEs(κ, p, Φ0, W0, B0, A0, r_max, E_min=0, E_max=(p ? M_p : M_n))
    f(E) = boundaryValue(κ, p, E, Φ0, W0, B0, A0, r_max)
    return find_zeros(f, (E_min, E_max))
end

"Find all orbitals and return two lists containing κ values and corresponding energies for a single species where
    the other parameters are defined above"
function findAllOrbitals(p, Φ0, W0, B0, A0, r_max, E_min=0, E_max=(p ? M_p : M_n))
    κs = Int[]
    Es = Float64[]
    # start from κ=-1 and go both up and down
    for direction in [-1, 1]
        for κ in direction * (1:100) # cutoff is 100
            new_Es = findEs(κ, p, Φ0, W0, B0, A0, r_max, E_min, E_max)
            if isempty(new_Es); break; end
            append!(Es, new_Es)
            append!(κs, fill(κ, length(new_Es)))
        end
    end
    return (κs, Es)
end

"For a given list of κ values with corresponding energies, attempt to fill N lowest lying orbitals and return occupancy numbers"
function fillNucleons(N::Int, κs, Es)
    sort_i = sortperm(Es)

    occ = zeros(Int, length(κs))

    for i in sort_i
        if N ≤ 0; break; end;
        max_occ = 2 * j_κ(κs[i]) + 1
        occ[i] = min(max_occ, N)
        N -= occ[i]
    end

    N == 0 || @warn "All orbitals could not be filled"
    return occ
end

"Total angular momentum j for a given κ value"
j_κ(κ::Int) = abs(κ) - 1/2

"Orbital angular momentum l for a given κ value"
l_κ(κ::Int) = abs(κ) - (κ < 0) # since true = 1 and false = 0

"Calculate scalar and vector densities of a nucleon species on [0,r_max] divided into (divs+1) points and returns them as vectors (ρ_s, ρ_v) where
    the other parameters are defined above"
function calculateNucleonDensity(N, p, Φ0, W0, B0, A0, r_max, divs, E_min=0, E_max=(p ? M_p : M_n))
    κs, Es = findAllOrbitals(p, Φ0, W0, B0, A0, r_max, E_min, E_max)
    occs = fillNucleons(N, κs, Es)
    
    r2s = (collect ∘ range)(0, r_max, length=divs+1).^2 |> transpose

    ρ_s = zeros(divs + 1)
    ρ_v = zeros(divs + 1)
    
    for (κ, E, occ) in zip(κs, Es, occs)
        wf = solveWf(κ, p, E, Φ0, W0, B0, A0, r_max, divs) # TODO: Needs to be normalized
        wf2 = wf .* wf
        ρ = (occ / (4 * pi)) * (wf2 ./ r2s)  # 2j+1 factor is accounted in the occupancy number

        ρ_s += ρ[1, :] - ρ[2, :] # g^2 - f^2
        ρ_v += ρ[1, :] + ρ[2, :] # g^2 + f^2
    end

    return (ρ_s, ρ_v)
end