using DifferentialEquations

const ħc = 197.327 # ħc in MeVfm

# Values taken from B. G. Todd and J. Piekarewicz, Relativistic Mean-Field Study of Neutron-Rich Nuclei, Phys. Rev. C 67, 044317 (2003)
const m_ρ = 763
const m_ω = 782.5
const m_s = 508.194
const g2_s = 104.3871
const g2_v = 165.5854
const g2_ρ = 79.6000
const κ = 3.8599
const λ = -0.01591
const ζ = 0.00
const Λv = 0 # value unknown

const r_reg = 1E-6 # regulator for the centrifugal term in fm

"Coupled radial meson equations that returns ddu=[ddΦ0, ddW0, ddB0, ddA0] in-place where
    du=[dΦ0, dW0, dB0, dA0] are the first derivatives of meson fields evaluated at r,
    u=[Φ0, W0, B0, A0] are the meson fields evaluated at r,
    κ is the generalized angular momentum,
    ρ_sp, ρ_vp are the scalar and vector proton densities as funtions of r in fm,
    ρ_sn, ρ_vn are the scalar and vector neutron densities as funtions of r in fm,
    Φ0, W0, B0, A0 are the mean-field potentials (couplings included) in MeV as functions of r in fm,
    r is the radius in fm.
Reference: P. Giuliani, K. Godbey, E. Bonilla, F. Viens, and J. Piekarewicz, Frontiers in Physics 10, (2023)"
function mesons!(ddu, du, u, (ρ_sp, ρ_vp, ρ_sn, ρ_vn), r)
    (dΦ0, dW0, dB0, dA0) = du
    (Φ0, W0, B0, A0) = u

    ddΦ0 = -(2/(r + r_reg)) * dΦ0 + m_s^2 * Φ0 / ħc + g2_s * ((κ/2) * Φ0^2 + (λ/6) * Φ0^3 - (ρ_sp(r) + ρ_sn(r))) / ħc
    ddW0 = -(2/(r + r_reg)) * dW0 + m_ω^2 * W0 / ħc + g2_v * ((ζ/6) * W0^3 + 2 * Λv * B0^2 * W0 - (ρ_vp(r) + ρ_vn(r))) / ħc
    ddB0 = -(2/(r + r_reg)) * dB0 + m_ρ^2 * B0 / ħc + g2_ρ * (2 * Λv * W0^2 * B0 - (ρ_vp(r) - ρ_vn(r)) / 2) / ħc
    ddA0 = -(2/(r + r_reg)) * dA0 - ρ_vp(r) / ħc # e goes here?

    ddu[1] = ddΦ0
    ddu[2] = ddW0
    ddu[3] = ddB0 
    ddu[4] = ddA0
end

"Solve meson equations and return the wave functions u(r)=[Φ0(r), W0(r), B0(r), A0(r)] where
    divs is the number of mesh divisions so solution would be returned as a 4×(1+divs) matrix,
    the other parameters are the same from mesons!(...)."
function solveWfs(ρ_sp, ρ_vp, ρ_sn, ρ_vn, r_max, divs)
    prob = SecondOrderODEProblem(mesons!, BigFloat.([0, 0, 0, 0]), BigFloat.([0, 0, 0, 0]), (r_max, 0))
    sol = solve(prob, Feagin14(), p=(ρ_sp, ρ_vp, ρ_sn, ρ_vn), saveat=r_max/divs)
    wfs = hcat(sol.u...)

    return wfs
end
