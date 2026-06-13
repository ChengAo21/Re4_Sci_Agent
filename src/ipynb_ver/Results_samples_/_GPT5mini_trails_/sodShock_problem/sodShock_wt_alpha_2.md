### 🛌🏻 Re4gent with Multi-Modal Review

#### Configs 🏷️:
MaxCount: 2, MultiModalReview: False, AlternativeNum: 2, MaxComments: 5, ConsModule_ON: True, CUT_OUTEXT: 2000

### prob_todo

The Euler equations for compressible flow in 1-D are given by:
\begin{cases}
\frac{\partial \rho}{\partial t} + \frac{\partial (\rho u)}{\partial x} = 0, \\
\frac{\partial (\rho u)}{\partial t} + \frac{\partial (\rho u^2 + p)}{\partial x} = 0, \\
\frac{\partial (\rho E)}{\partial t} + \frac{\partial (\rho E u + p u)}{\partial x} = 0,
\end{cases}
in the interval \( x \in [0,1] \), where \( E = \frac{1}{2} u^2 + \frac{p}{(\gamma - 1)\rho} \), \( \gamma = 1.4 \).
At \( t = 0 \), the initial condition is
(\rho, u, p) = 
\begin{cases}
(1.0, 0.0, 1.0), & 0 < x \leq 0.5, \\
(0.125, 0.0, 0.1), & 0.5 \leq x < 1,
\end{cases}

Implement a stable and efficient method to solve this problem.
Strictly use visc instead of nu throughout the code to avoid confusion with u.
Ensure all constants are passed explicitly to functions or defined globally to prevent NameError.
Plot the density, velocity, and pressure at t=0.1,0.2,0.3 in three figures, respectively.
Just save figs do not use plt.show() in the code.

[HINTS]:
Print diagnostic information that helps assess the physical correctness of the solution, such as shock position and the absence of spurious oscillations.
Global conservation outputs are intended for qualitative reference only and should not be treated as strict acceptance criteria.
Do not use 'np.trapz' as it is removed in NumPy 2.0, implement the integration manually.


### expanded_prob
We must solve the 1-D compressible Euler equations (hyperbolic system of conservation laws) on x in [0,1] with Riemann initial data (the classical Sod shock tube). The unknowns are density rho, momentum rho*u, and total energy rho*E. The fluxes are nonlinear; the solution develops discontinuities (shock, contact, rarefaction) that must be captured without spurious oscillations and while preserving physically admissible states (rho>0, p>0). Numerical implementation must be conservative, stable under a CFL constraint, and efficient for repeated runs at final times t = 0.1, 0.2, 0.3. Practical programming constraints: use the name visc (not nu) if an artificial/dissipative viscosity parameter is needed; pass all constants explicitly (gamma, visc, CFL, domain limits, Nx, final times, etc.) to functions to avoid NameError; do not call np.trapz — compute integrals manually; save figures to files and do not call plt.show().

Primary mathematical and numerical challenges:
- Nonlinearity and hyperbolicity: characteristic speeds depend on local state; must choose flux solvers that capture shocks and contacts robustly and satisfy entropy condition.
- Discontinuities in initial data: require shock-capturing schemes (finite-volume) and non-oscillatory reconstruction to avoid Gibbs phenomena near shocks and contacts.
- Positivity preservation: numerical method must avoid producing negative density or pressure; high-order reconstructions and linearized solvers (e.g., Roe) may violate positivity unless corrective steps are taken.
- Accuracy vs robustness tradeoff: low-order diffusive schemes are robust but smear features; high-order schemes (MUSCL, WENO) preserve sharpness but need limiters/positivity fixers and are more expensive.
- Time integration and stability: explicit time stepping constrained by CFL condition based on maximum signal speed; SSPRK methods help maintain nonlinear stability for TVD-limited reconstructions.
- Conservation diagnostics: compute global mass, momentum, and energy via discrete integrals (manual sum of cell average * dx) to monitor conservation; careful when using boundary conditions or artificial viscosity that modify global invariants.
- Implementation hygiene: ensure all physical and numerical constants are passed explicitly; name artificial viscosity variable visc; implement boundary conditions (transmissive/outflow is typical for Sod) consistently with finite-volume stencil operations; save plots without showing; provide diagnostic outputs (shock location from density gradient peak, min/max checks, total variation) to assess solution quality.

The accepted solution will be a conservative finite-volume shock-capturing solver with either a robust approximate Riemann flux (HLLC or Rusanov/HLL) combined with a non-oscillatory reconstruction (MUSCL with limiter or WENO) and an SSP time integrator. Important diagnostics: shock/contact positions, absence of spurious oscillations, global integral checks (mass, momentum, energy) computed by manual discrete sums, and minima of density/pressure for physical admissibility.

### solution_plans
	Current Stage [A/2]
solu_name='Plan A: Finite-Volume MUSCL-Hancock + HLLC Riemann Solver + SSPRK3' content='Governing idea:\nUse a conservative finite-volume scheme with second-order spatial accuracy via a MUSCL-Hancock predictor-corrector reconstruction and the HLLC approximate Riemann solver for fluxes. Time integration uses a 3-stage SSPRK (strong-stability-preserving Runge–Kutta) to keep nonlinear stability with limiters. Transmissive (zero-gradient) boundary conditions are applied at domain ends. Diagnostics are computed from cell averages using manual discrete integrals.\n\nAlgorithmic steps (numbered):\n1) Define constants and inputs explicitly: gamma, visc (if using optional artificial viscosity), CFL, x_left=0.0, x_right=1.0, Nx, final_times=[0.1,0.2,0.3], output filenames. Build uniform cell centers and cell width dx = (x_right-x_left)/Nx.\n2) Initialize conservative variables U = [rho, rho*u, rho*E] from piecewise initial conditions at cell centers. Convert primitives->conservatives with E = 0.5*u^2 + p/((gamma-1)*rho).\n3) Implement boundary ghost cells (2 layers) with transmissive (extrapolate) boundary conditions to support the reconstruction stencils.\n4) Reconstruction (MUSCL-Hancock): for each primitive variable (rho, u, p) compute limited slopes using a robust limiter (minmod or van Albada). Reconstruct left/right states at each cell interface and use Hancock predictor to evolve these reconstructed states by half time-step dt/2 using local flux gradients to get time-centered states.\n5) Compute interface fluxes using the HLLC Riemann solver applied to the predicted left/right primitive states. Ensure conversion between primitive and conservative variables is robust and includes safeguards to keep rho and p positive; if positivity is violated, fallback to more diffusive flux (HLL/Rusanov) or apply a local positivity fix.\n6) Update cell conservative variables by integrating flux differences: U^{n+1} = U^n - (dt/dx) * (F_{i+1/2} - F_{i-1/2}). Use SSPRK3 time-stepping: perform three stage updates with recomputed fluxes and boundary fills each substage. Choose dt via CFL: dt = CFL * dx / max_i(|u_i| + c_i) where c_i = sqrt(gamma*p_i/rho_i).\n7) Optional: add small linear artificial viscosity term (if using visc): e.g., Q_i = visc * max_lambda * (U_{i+1}-2U_i+U_{i-1}) to damp grid-scale oscillations. Always call this parameter visc explicitly and test sensitivity. Implement visc conservatively (as flux-like diffusive term) so global conservation is preserved except for intended dissipation.\n8) At each desired output time (0.1,0.2,0.3) or when crossing them, compute diagnostics: shock position (index of maximum density gradient or location of largest |d rho/dx|), total mass/momentum/energy via manual sums Mass = sum(rho_i * dx) over cells, Momentum = sum((rho*u)_i * dx), Energy = sum((rho*E)_i * dx). Compute min(p) and min(rho), and total variation TV(rho)=sum(|rho_{i+1}-rho_i|). Print these diagnostic values for qualitative assessment.\n9) Save plots of density, velocity (u = (rho*u)/rho), and pressure (p via equation of state) for each requested time into separate figure files. Do not call plt.show(). Ensure consistent axis limits and labels, and overlay or annotate shock position if desired.\n\nStability / accuracy / complexity / efficiency limitations:\n- Stability is limited by the CFL condition; for MUSCL-Hancock + SSPRK3 use CFL ≲ 0.5 (typical safe choice 0.2–0.4 depending on limiter and HLLC implementation). Larger CFL may produce instability or loss of TVD properties.\n- Accuracy: scheme is second-order in smooth regions but reduces to first-order at discontinuities due to limiting; contact discontinuity resolution depends on HLLC quality and limiter choice.\n- Efficiency: computational cost per time step includes reconstruction and HLLC evaluation per interface (O(N) per stage). SSPRK3 triples the flux computations relative to single-step methods, but yields better nonlinear stability.\n- Robustness: HLLC captures contacts better than HLL/Rusanov but requires careful handling of intermediate states; Roe linearization is more accurate but may produce nonphysical states (negative p/rho) without entropy fix or positivity limiter. If required, fallback to more diffusive Rusanov/HLL for troubled cells.\n- Positivity: to avoid negative density/pressure, include protective checks and fallbacks. Artificial viscosity (visc) may help but should be tuned; inappropriately large visc will oversmooth the solution.\n- Conservation: scheme is conservative by finite-volume formulation; adding non-conservative source-like corrections (improper visc implementation) can break conservation—implement visc as conservative flux differences.\n- Implementation hygiene: make sure all constants are passed explicitly to functions; use manual discrete sums for diagnostics (avoid np.trapz).'

	Current Stage [B/2]
solu_name='Plan B: High-order WENO5 Reconstruction + HLLC Flux + SSPRK3 (sharp, higher-order option)' content='Governing idea:\nUse a high-order finite-volume method with 5th-order WENO reconstruction for non-oscillatory high-resolution spatial interpolation to interfaces, combined with the HLLC Riemann solver for accurate shock/contact capture. Time-advance with SSPRK3 for nonlinear stability. This yields higher accuracy in the smooth regions (e.g., inside rarefaction) and much sharper contacts than second-order schemes, at the cost of extra computational work and complexity (stencil of 5 points).\n\nAlgorithmic steps (numbered):\n1) Set domain and constants explicitly: gamma, visc (optional small dissipation), CFL, Nx, x_left, x_right, final_times=[0.1,0.2,0.3]. Compute cell centers and dx.\n2) Initialize conservative variables U from primitive initial data; verify positive rho and p; compute primitive state arrays for reconstruction.\n3) Implement ghost cells (at least 3 on each side for WENO5) with transmissive/outflow extrapolation to supply stencil values.\n4) WENO5 reconstruction: for each primitive (or characteristic) field, perform WENO5 to compute left/right interface states at every interface. For improved robustness and reduced spurious oscillations across strong shocks, implement reconstruction in characteristic variables (i.e., transform primitive vectors to characteristic space using local Roe-average or local flux Jacobian eigenvectors), perform WENO limiting there, then transform back to primitive space.\n5) At each interface use the HLLC Riemann solver with the reconstructed left/right states to compute numerical flux. Ensure positivity by checking reconstructed states; if a negative pressure or density would occur, fall back to more robust reconstruction (component-wise WENO) or use HLL/Rusanov flux for that interface.\n6) Time integration: SSPRK3 with dt computed by CFL using the maximum signal speed (|u|+c). Update U conservatively at each substage and reapply boundary ghost cell fills before reconstruction.\n7) Optionally include an explicit dissipative filter or small visc-based diffusive flux Q that is proportional to the second difference of U scaled by local max_wave_speed to damp unresolved oscillations; provide visc as an explicit parameter and document its effect in diagnostics.\n8) Diagnostics and post-processing: at each requested time compute manual discrete integrals for Mass, Momentum, Energy: sum(U_i * dx) across cells (no np.trapz). Locate shock/contact via maxima of |d rho/dx| or by tracking density jump position. Compute minima of rho and p and total variation TV(rho). Print these values as run diagnostics.\n9) Output: for each requested time, compute primitive fields (rho, u, p) from U and save three figures (one each for density, velocity, pressure) to files. Use clear filenames that include time and method. Do not call plt.show().\n\nStability / accuracy / complexity / efficiency limitations:\n- Stability: explicit SSPRK + WENO requires small dt controlled by CFL; for WENO5 typical stable CFL is ≲ 0.4 (often chosen 0.2–0.35 for safety). WENO in characteristic space improves robustness but increases cost.\n- Accuracy: WENO5 attains fifth-order accuracy in smooth regions and avoids Gibbs near discontinuities, preserving sharp features like contacts better than second-order methods. However, near steep gradients it reduces locally to lower order.\n- Complexity: WENO5 uses wider stencils and characteristic decomposition (if used), increasing arithmetic intensity and memory accesses; cost per time step is roughly an order of magnitude higher than a simple first-order method and several times higher than MUSCL-Hancock.\n- Efficiency: use vectorized loops and avoid Python-level loops where possible; if performance is critical consider numba/C or compiled kernels for WENO weights and flux computations.\n- Robustness: high-order reconstructions can produce negative densities/pressures in extreme cases—implement positivity fixers or fallback to component-wise or more diffusive reconstructions at troubled cells. HLLC is accurate but requires robust state handling.\n- Conservation: finite-volume update is conservative; any added filter or visc must be introduced as conservative flux differences to preserve integrated invariants (apart from intended dissipative losses).\n- Implementation hygiene: always pass gamma, visc, CFL, and other constants explicitly to functions; compute diagnostics via manual sums (Mass = sum(rho_i*dx), etc.) and implement shock-location detection consistently. Ensure figure saving uses explicit file names and omit plt.show().'



### technical_spec
	Current Stage [A/3]
Solver architecture:
- Finite-volume explicit solver with cell-centered conservative variables U (3 x Ntot including ghosts).
- MUSCL second-order spatial reconstruction in primitive variables with a three-argument minmod limiter.
- HLLC approximate Riemann solver at interfaces for robust shock/contact resolution; optional viscous dissipation added as a conservative interface flux scaled by 'visc'.
- SSPRK3 time integrator for strong-stability-preserving time advancement.
- Preallocated arrays for primitives, slopes, interface states, fluxes and RHS reused at each step to avoid allocations inside the time loop.
- Diagnostics: shock position (max density gradient), global integrals (mass, momentum, energy) computed by manual discrete sums, minima of rho/p, and total variation of density.
- Outputs: three PNG files (density_times.png, velocity_times.png, pressure_times.png) with curves at t=0.1,0.2,0.3.

Data flow:
U (conservative) -> cons_to_prim -> primitive arrays -> limited slopes -> interface left/right primitives -> HLLC fluxes (+ viscous diffusive flux) -> RHS = -dF/dx -> SSPRK3 updates U -> repeat until requested output times.


	Current Stage [B/3]
This solver implements a finite-volume MUSCL (second-order) solver for the 1D Euler equations (Sod shock tube) with HLLC Riemann solver, optional interface viscous dissipation (visc), and SSPRK3 time integration. Data flow: conservative U -> cons_to_prim -> primitive arrays -> limited slopes (minmod3) -> left/right interface primitives (qL,qR) -> HLLC fluxes -> optional diffusive flux -> RHS = -dF/dx -> SSPRK3 update of U. Arrays (U, prim, slope, qL, qR, flux_if, rhs) are preallocated once and reused; boundary conditions are transmissive (zero-gradient) by copying nearest interior cells into ghosts. Diagnostics include shock position (max |drho|), global integrals (mass/momentum/energy via manual sums), minima of rho/p, total variation of density, and count of pressure repairs. Snapshots at requested times are saved and plotted into three PNG files (density_times.png, velocity_times.png, pressure_times.png). Key robustness fixes: correct interface indexing for divergence, assertion for interface shape, guard against zero dt, vectorized conservative pressure-floor repair, and clear output-saving logic using tolerances and saved flags.

	Current Stage [C/3]
- Modular finite-volume MUSCL solver for 1D Euler (Sod tube) with HLLC Riemann flux and SSPRK3 time integration.
- Data flow: conservative U (3xNtot) -> cons_to_prim -> primitive arrays -> limited slopes (minmod3) -> interface reconstructions qL/qR -> HLLC fluxes -> optional viscous diffusive flux -> RHS = -div(flux)/dx -> SSPRK3 update on U.
- Arrays preallocated once: U, prim, slope, qL, qR, flux_if, rhs; no reallocation inside time loop.
- Boundary conditions: transmissive zero-gradient by replicating nearest interior cell into all ghost slots; BCs enforced before RHS and after each RK substage to avoid ghost contamination.
- Conservation/corrections: discrete divergence corrected with ng offset; density-floor repair rescales momentum and energy conservatively; pressure-floor repair adjusts total energy conservatively.
- Diagnostics: mass/momentum/energy integrals, min(rho)/min(p), total variation of density, improved shock locator (max |drho| and right-most positive density jump), counts of density/pressure repairs.
- Outputs: saves three PNG figures (density_times.png, velocity_times.png, pressure_times.png) for requested snapshot times; prints quantitative diagnostics for each snapshot.



### python_codes

#### Script block1:

```python
#!/usr/bin/env python3
"""
MUSCL finite-volume solver for 1D Euler equations (Sod shock tube)
- Second-order MUSCL reconstruction (primitive variables) with minmod limiter
- HLLC Riemann solver at interfaces
- SSPRK3 time integration
- Optional interface viscous dissipation controlled by 'visc'

Outputs: three saved figures showing density, velocity, pressure at t=0.1,0.2,0.3
Diagnostics printed at each saved time.
"""

import numpy as np
import matplotlib.pyplot as plt

# --------------------------- Utility converters ---------------------------

def prim_to_cons(rho, u, p, gamma):
    """Convert primitive vars to conservative vector U = [rho, rho*u, rho*E].
    Inputs are arrays of same shape.
    """
    mom = rho * u
    E_spec = 0.5 * u * u + p / ((gamma - 1.0) * rho)
    Econs = rho * E_spec
    U = np.vstack((rho, mom, Econs))
    return U


def cons_to_prim(U, gamma):
    """Convert conservative U (3 x N) to primitive rho,u,p arrays.
    Returns rho,u,p as arrays.
    """
    rho = U[0].copy()
    mom = U[1].copy()
    Econs = U[2].copy()
    # avoid division by zero
    u = np.zeros_like(rho)
    mask = rho > 0.0
    u[mask] = mom[mask] / rho[mask]
    # p = (gamma-1)*(rho*E - 0.5*(rho*u^2))
    kinetic = 0.5 * mom * u
    p = (gamma - 1.0) * (Econs - kinetic)
    return rho, u, p


def flux_from_prim(rho, u, p, gamma):
    """Compute physical flux vector F = [rho*u, rho*u^2 + p, (rho*E + p)*u]."""
    mom = rho * u
    E_spec = 0.5 * u * u + p / ((gamma - 1.0) * rho)
    Econs = rho * E_spec
    F0 = mom
    F1 = mom * u + p
    F2 = (Econs + p) * u
    return np.vstack((F0, F1, F2))

# ----------------------------- Limiters ----------------------------------

def minmod3(a, b, c):
    """Elementwise three-argument minmod limiter."""
    result = np.zeros_like(a)
    s = np.sign(a) + np.sign(b) + np.sign(c)
    mask_pos = (s == 3)
    mask_neg = (s == -3)
    result[mask_pos] = np.minimum(np.minimum(np.abs(a[mask_pos]), np.abs(b[mask_pos])), np.abs(c[mask_pos]))
    result[mask_neg] = -np.minimum(np.minimum(np.abs(a[mask_neg]), np.abs(b[mask_neg])), np.abs(c[mask_neg]))
    return result

# ----------------------------- HLLC Solver -------------------------------

def hllc_flux(primL, primR, gamma):
    """Compute HLLC flux for left/right primitive states at one interface.
    primL, primR are tuples (rho,u,p) scalars or 1D arrays shape (Ninterface,)
    Returns flux array shape (3, ...)
    """
    rhoL, uL, pL = primL
    rhoR, uR, pR = primR
    # sound speeds
    cL = np.sqrt(np.maximum(gamma * pL / rhoL, 0.0))
    cR = np.sqrt(np.maximum(gamma * pR / rhoR, 0.0))
    # conservative states
    UL = prim_to_cons(rhoL, uL, pL, gamma)
    UR = prim_to_cons(rhoR, uR, pR, gamma)
    # fluxes
    FL = flux_from_prim(rhoL, uL, pL, gamma)
    FR = flux_from_prim(rhoR, uR, pR, gamma)
    # wave speed estimates
    SL = np.minimum(uL - cL, uR - cR)
    SR = np.maximum(uL + cL, uR + cR)

    # avoid divide by zero: if SL >= SR fallback to HLL
    use_hll = False
    mask_bad = SL >= SR
    if np.any(mask_bad):
        use_hll = True

    # middle wave speed S_M
    denom = (rhoL * (SL - uL) - rhoR * (SR - uR))
    # safe compute S_M (where denom not zero)
    S_M = np.zeros_like(denom)
    safe = denom != 0.0
    S_M[safe] = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR))[safe] / denom[safe]
    # build star states
    # E_spec per unit mass
    EL = UL[2] / rhoL
    ER = UR[2] / rhoR
    rhoSL = rhoL * (SL - uL) / (SL - S_M)
    rhoSR = rhoR * (SR - uR) / (SR - S_M)
    UstarL = np.zeros_like(UL)
    UstarR = np.zeros_like(UR)
    UstarL[0] = rhoSL
    UstarL[1] = rhoSL * S_M
    UstarL[2] = rhoSL * (EL + (S_M - uL) * (S_M + pL / (rhoL * (SL - uL))))
    UstarR[0] = rhoSR
    UstarR[1] = rhoSR * S_M
    UstarR[2] = rhoSR * (ER + (S_M - uR) * (S_M + pR / (rhoR * (SR - uR))))

    # allocate flux
    F = np.zeros_like(FL)

    # Regions
    # SL >= 0 -> FL
    mask_FL = SL >= 0.0
    if np.any(mask_FL):
        F[:, mask_FL] = FL[:, mask_FL]
    # SL < 0 <= S_M -> F = FL + SL*(U*_L - UL)
    mask_Lstar = (SL < 0.0) & (S_M >= 0.0)
    if np.any(mask_Lstar):
        F[:, mask_Lstar] = FL[:, mask_Lstar] + (SL[mask_Lstar] * (UstarL[:, mask_Lstar] - UL[:, mask_Lstar]))
    # S_M < 0 <= SR -> F = FR + SR*(U*_R - UR)
    mask_Rstar = (S_M < 0.0) & (SR >= 0.0)
    if np.any(mask_Rstar):
        F[:, mask_Rstar] = FR[:, mask_Rstar] + (SR[mask_Rstar] * (UstarR[:, mask_Rstar] - UR[:, mask_Rstar]))
    # SR < 0 -> FR
    mask_FR = SR < 0.0
    if np.any(mask_FR):
        F[:, mask_FR] = FR[:, mask_FR]

    # For any bad cells where SL>=SR fallback to HLL (robust but diffusive)
    if use_hll:
        # HLL flux: (SR*FL - SL*FR + SR*SL*(UR-UL)) / (SR - SL)
        denom_hll = (SR - SL)
        # avoid division by zero
        denom_hll[denom_hll == 0.0] = 1e-16
        F_hll = (SR * FL - SL * FR + SR * SL * (UR - UL)) / denom_hll
        F[:, mask_bad] = F_hll[:, mask_bad]

    return F

# ----------------------------- RHS builder -------------------------------

def compute_RHS(U, dx, gamma, visc,
                prim_arr, slope, qL, qR, flux_interface, rhs_out):
    """Compute RHS = -dF/dx for interior cells. Uses preallocated arrays passed in.
    Arrays shapes:
      U: (3, Ntot)
      prim_arr: (3, Ntot)
      slope: (3, Ntot)
      qL, qR: (3, Ninterfaces)   (interfaces are between cells)
      flux_interface: (3, Ninterfaces)
      rhs_out: (3, Ncells)
    """
    Ntot = U.shape[1]
    # Apply transmissive (zero-gradient) BCs by copying nearest interior cells into ghosts
    ng = 2
    # left ghosts
    U[:, 0:ng] = U[:, ng:2*ng][:, ::-1]  # simple mirror (extrapolate); keep zero-gradientish
    # right ghosts
    U[:, -ng:] = U[:, -2*ng:-ng][:, ::-1]

    # primitives
    rho, u, p = cons_to_prim(U, gamma)
    prim_arr[0, :] = rho
    prim_arr[1, :] = u
    prim_arr[2, :] = p

    # limited slopes (minmod3) on primitive variables
    theta = 1.5
    # compute differences
    dqL = prim_arr[:, 1:-1] - prim_arr[:, 0:-2]
    dqR = prim_arr[:, 2:] - prim_arr[:, 1:-1]
    dqC = 0.5 * (prim_arr[:, 2:] - prim_arr[:, 0:-2])
    # apply minmod3
    slope[:, 1:-1] = minmod3(theta * dqL, dqC, theta * dqR)
    # set slopes at boundaries to zero (ghosts)
    slope[:, 0] = slope[:, 1]
    slope[:, -1] = slope[:, -2]

    # reconstruct at interfaces: between cell i and i+1 we take left from i and right from i+1
    # number of interfaces = Ntot - 1
    Nint = Ntot - 1
    # qL at interface i+1/2 comes from cell i: q_i + 0.5*slope_i
    qL[:, :] = prim_arr[:, :-1] + 0.5 * slope[:, :-1]
    # qR at interface i+1/2 comes from cell i+1: q_{i+1} - 0.5*slope_{i+1}
    qR[:, :] = prim_arr[:, 1:] - 0.5 * slope[:, 1:]

    # Ensure positive density and pressure in reconstructed states (simple floor)
    eps = 1e-12
    qL[0, :] = np.maximum(qL[0, :], eps)
    qR[0, :] = np.maximum(qR[0, :], eps)
    qL[2, :] = np.maximum(qL[2, :], eps)
    qR[2, :] = np.maximum(qR[2, :], eps)

    # Riemann fluxes at interfaces (HLLC)
    flux_interface[:, :] = hllc_flux((qL[0, :], qL[1, :], qL[2, :]),
                                     (qR[0, :], qR[1, :], qR[2, :]), gamma)

    # Add small diffusive (viscous) flux to stabilize if visc > 0
    if visc > 0.0:
        # compute local max wave speed at interfaces
        cL = np.sqrt(np.maximum(gamma * qL[2, :] / qL[0, :], 0.0))
        cR = np.sqrt(np.maximum(gamma * qR[2, :] / qR[0, :], 0.0))
        smax = np.maximum(np.abs(qL[1, :]) + cL, np.abs(qR[1, :]) + cR)
        # diffusive flux proportional to state jump
        flux_interface += (-0.5 * visc) * (smax[np.newaxis, :] * (prim_to_cons(qR[0, :], qR[1, :], qR[2, :], gamma) -
                                                                    prim_to_cons(qL[0, :], qL[1, :], qL[2, :], gamma)))

    # compute RHS for interior cells (exclude ghost cells)
    # flux_interface has shape (3, Ntot-1); RHS for cell i is -(F_{i+1/2} - F_{i-1/2})/dx
    rhs_out[:, :] = 0.0
    rhs_out[:, 1:-1] = - (flux_interface[:, 2:] - flux_interface[:, 1:-1]) / dx
    # Note: rhs_out for ghost cells will be ignored by integrator
    return rhs_out

# ----------------------------- Main solver --------------------------------

def solve_sod(nx=400, x_left=0.0, x_right=1.0,
              gamma=1.4, visc=0.0, CFL=0.45,
              final_times=(0.1, 0.2, 0.3)):
    # grid
    Nx = nx
    ng = 2  # ghost cells
    Ntot = Nx + 2 * ng
    x = np.linspace(x_left + 0.5 * (x_right - x_left) / Nx,
                    x_right - 0.5 * (x_right - x_left) / Nx,
                    Nx)
    dx = (x_right - x_left) / Nx

    # allocate arrays (reuse to avoid allocations inside time loop)
    U = np.zeros((3, Ntot))
    prim = np.zeros((3, Ntot))
    slope = np.zeros((3, Ntot))
    qL = np.zeros((3, Ntot - 1))
    qR = np.zeros((3, Ntot - 1))
    flux_if = np.zeros((3, Ntot - 1))
    rhs = np.zeros((3, Ntot))

    # initialize cell-centered interior primitive states (Sod tube)
    # interior indices
    i0 = ng
    i1 = ng + Nx
    # set primitives for interior
    # left state for x <= 0.5
    xcell = np.linspace(x_left + dx / 2.0, x_right - dx / 2.0, Nx)
    rho0 = np.where(xcell <= 0.5, 1.0, 0.125)
    u0 = np.zeros_like(xcell)
    p0 = np.where(xcell <= 0.5, 1.0, 0.1)
    U[:, i0:i1] = prim_to_cons(rho0, u0, p0, gamma)

    # fill ghost cells with nearest interior (transmissive)
    U[:, :ng] = U[:, ng:2*ng][:, ::-1]
    U[:, -ng:] = U[:, -2*ng:-ng][:, ::-1]

    # time stepping
    t = 0.0
    Tfinal_all = sorted(final_times)
    next_output_index = 0
    saved = {t_out: None for t_out in final_times}

    # diagnostics storage for plotting
    save_rho = {t_out: None for t_out in final_times}
    save_u = {t_out: None for t_out in final_times}
    save_p = {t_out: None for t_out in final_times}

    # main time loop until last final time
    Tfinal = max(final_times)

    # Precompute max iterations guard
    max_steps = 500000
    step = 0
    while t < Tfinal and step < max_steps:
        # compute primitive for dt estimation
        rho, u, p = cons_to_prim(U, gamma)
        c = np.sqrt(np.maximum(gamma * p / rho, 0.0))
        max_wave = np.max(np.abs(u) + c)
        if max_wave <= 0.0:
            dt = 1e-6
        else:
            dt = CFL * dx / max_wave
        # avoid stepping past next output
        next_output_time = Tfinal_all[next_output_index] if next_output_index < len(Tfinal_all) else Tfinal
        if t + dt > next_output_time:
            dt = next_output_time - t

        # SSPRK3 stages
        Un = U.copy()

        # stage 1
        compute_RHS(U, dx, gamma, visc, prim, slope, qL, qR, flux_if, rhs)
        U1 = Un + dt * rhs

        # stage 2
        compute_RHS(U1, dx, gamma, visc, prim, slope, qL, qR, flux_if, rhs)
        U2 = 0.75 * Un + 0.25 * (U1 + dt * rhs)

        # stage 3
        compute_RHS(U2, dx, gamma, visc, prim, slope, qL, qR, flux_if, rhs)
        U = (1.0 / 3.0) * Un + (2.0 / 3.0) * (U2 + dt * rhs)

        t += dt
        step += 1

        # after update ensure physical admissibility (small floor)
        rho_new = U[0, :]
        mom_new = U[1, :]
        E_new = U[2, :]
        eps = 1e-12
        neg_rho = rho_new < eps
        if np.any(neg_rho):
            rho_new[neg_rho] = eps
            U[0, :] = rho_new
            # adjust momentum and energy accordingly if needed
        # compute pressure and floor
        _, u_tmp, p_tmp = cons_to_prim(U, gamma)
        neg_p = p_tmp < 1e-12
        if np.any(neg_p):
            # if negative pressure occurs, slightly increase it by resetting internal energy
            # conservative fix: set pressure floor uniformly where needed
            p_tmp[neg_p] = 1e-12
            # reconstruct total energy from new pressure
            idx = np.where(neg_p)[0]
            for j in idx:
                rhoj = U[0, j]
                muj = U[1, j]
                uj = muj / rhoj
                E_spec = 0.5 * uj * uj + p_tmp[j] / ((gamma - 1.0) * rhoj)
                U[2, j] = rhoj * E_spec

        # check if we reached the next required output time
        while next_output_index < len(Tfinal_all) and abs(t - Tfinal_all[next_output_index]) < 1e-10 or (
                next_output_index < len(Tfinal_all) and t >= Tfinal_all[next_output_index] - 1e-14):
            t_out = Tfinal_all[next_output_index]
            # compute primitives for interior cells
            rho_i, u_i, p_i = cons_to_prim(U[:, ng:-ng], gamma)
            save_rho[t_out] = rho_i.copy()
            save_u[t_out] = u_i.copy()
            save_p[t_out] = p_i.copy()

            # diagnostics
            # shock position estimated by max density gradient
            drho = np.abs(np.diff(rho_i))
            max_idx = np.argmax(drho)
            xcell_centers = x
            # pick position between cells
            x_shock = 0.5 * (xcell_centers[max_idx] + xcell_centers[max_idx + 1])
            # global integrals (manual sums)
            mass = np.sum(rho_i) * dx
            momentum = np.sum((rho_i * u_i)) * dx
            energy = np.sum(((U[2, ng:-ng]))) * dx  # U[2] interior already rho*E
            min_rho = np.min(rho_i)
            min_p = np.min(p_i)
            TV_rho = np.sum(np.abs(np.diff(rho_i)))

            print(f"\nSnapshot at t = {t_out:.3f}")
            print(f"  Shock position (est.): x = {x_shock:.5f}")
            print(f"  Mass = {mass:.6f}, Momentum = {momentum:.6f}, Energy = {energy:.6f}")
            print(f"  min(rho) = {min_rho:.6e}, min(p) = {min_p:.6e}, TV(rho) = {TV_rho:.6e}")

            next_output_index += 1
            if next_output_index >= len(Tfinal_all):
                break
            next_output_time = Tfinal_all[next_output_index]

    if step >= max_steps:
        print("Warning: reached max steps limit")

    # plotting: overlay three times into three separate figures
    times = list(Tfinal_all)

    # prepare x for interior
    x_interior = x

    # density
    plt.figure(figsize=(8, 4))
    for t_out in times:
        rho_plot = save_rho[t_out]
        plt.plot(x_interior, rho_plot, label=f't={t_out:.2f}')
    plt.xlabel('x')
    plt.ylabel('Density')
    plt.title('Density at different times (MUSCL-Hancock-like + HLLC)')
    plt.legend()
    plt.grid(True)
    plt.savefig('density_times.png', dpi=200)
    plt.close()

    # velocity
    plt.figure(figsize=(8, 4))
    for t_out in times:
        u_plot = save_u[t_out]
        plt.plot(x_interior, u_plot, label=f't={t_out:.2f}')
    plt.xlabel('x')
    plt.ylabel('Velocity')
    plt.title('Velocity at different times (MUSCL-like + HLLC)')
    plt.legend()
    plt.grid(True)
    plt.savefig('velocity_times.png', dpi=200)
    plt.close()

    # pressure
    plt.figure(figsize=(8, 4))
    for t_out in times:
        p_plot = save_p[t_out]
        plt.plot(x_interior, p_plot, label=f't={t_out:.2f}')
    plt.xlabel('x')
    plt.ylabel('Pressure')
    plt.title('Pressure at different times (MUSCL-like + HLLC)')
    plt.legend()
    plt.grid(True)
    plt.savefig('pressure_times.png', dpi=200)
    plt.close()

    print('\nSaved figures: density_times.png, velocity_times.png, pressure_times.png')
    return save_rho, save_u, save_p


if __name__ == '__main__':
    # Explicit constants passed to solver
    gamma_val = 1.4
    visc_val = 0.0  # small positive value (e.g., 0.1) increases dissipation
    CFL_val = 0.45
    Nx = 400
    final_times_tuple = (0.1, 0.2, 0.3)

    solve_sod(nx=Nx, x_left=0.0, x_right=1.0,
              gamma=gamma_val, visc=visc_val, CFL=CFL_val,
              final_times=final_times_tuple)
```



#### Script block2:

```python
#!/usr/bin/env python3
"""
MUSCL finite-volume solver for 1D Euler equations (Sod shock tube)
- Second-order MUSCL reconstruction (primitive variables) with minmod limiter
- HLLC Riemann solver at interfaces
- SSPRK3 time integration
- Optional interface viscous dissipation controlled by 'visc'

Outputs: three saved figures showing density, velocity, pressure at t=0.1,0.2,0.3
"""

import numpy as np
import matplotlib.pyplot as plt

# --------------------------- Utility converters ---------------------------

def prim_to_cons(rho, u, p, gamma):
    """Convert primitive vars to conservative vector U = [rho, rho*u, rho*E]."""
    mom = rho * u
    E_spec = 0.5 * u * u + p / ((gamma - 1.0) * rho)
    Econs = rho * E_spec
    U = np.vstack((rho, mom, Econs))
    return U


def cons_to_prim(U, gamma):
    """Convert conservative U (3 x N) to primitive rho,u,p arrays."""
    rho = U[0].copy()
    mom = U[1].copy()
    Econs = U[2].copy()
    # avoid division by zero
    u = np.zeros_like(rho)
    mask = rho > 0.0
    u[mask] = mom[mask] / rho[mask]
    kinetic = 0.5 * mom * u
    p = (gamma - 1.0) * (Econs - kinetic)
    return rho, u, p


def flux_from_prim(rho, u, p, gamma):
    """Compute physical flux vector F = [rho*u, rho*u^2 + p, (rho*E + p)*u]."""
    mom = rho * u
    E_spec = 0.5 * u * u + p / ((gamma - 1.0) * rho)
    Econs = rho * E_spec
    F0 = mom
    F1 = mom * u + p
    F2 = (Econs + p) * u
    return np.vstack((F0, F1, F2))

# ----------------------------- Limiters ----------------------------------

def minmod3(a, b, c):
    """Elementwise three-argument minmod limiter."""
    result = np.zeros_like(a)
    s = np.sign(a) + np.sign(b) + np.sign(c)
    mask_pos = (s == 3)
    mask_neg = (s == -3)
    if np.any(mask_pos):
        result[mask_pos] = np.minimum(np.minimum(np.abs(a[mask_pos]), np.abs(b[mask_pos])), np.abs(c[mask_pos]))
    if np.any(mask_neg):
        result[mask_neg] = -np.minimum(np.minimum(np.abs(a[mask_neg]), np.abs(b[mask_neg])), np.abs(c[mask_neg]))
    return result

# ----------------------------- HLLC Solver -------------------------------

def hllc_flux(primL, primR, gamma):
    """Compute HLLC flux for left/right primitive states at interfaces.
    primL, primR are tuples (rho,u,p) where each element can be arrays.
    Returns flux array shape (3, Ninterfaces).
    """
    rhoL, uL, pL = primL
    rhoR, uR, pR = primR
    # sound speeds
    cL = np.sqrt(np.maximum(gamma * pL / rhoL, 0.0))
    cR = np.sqrt(np.maximum(gamma * pR / rhoR, 0.0))
    # conservative states
    UL = prim_to_cons(rhoL, uL, pL, gamma)
    UR = prim_to_cons(rhoR, uR, pR, gamma)
    # fluxes
    FL = flux_from_prim(rhoL, uL, pL, gamma)
    FR = flux_from_prim(rhoR, uR, pR, gamma)
    # wave speed estimates
    SL = np.minimum(uL - cL, uR - cR)
    SR = np.maximum(uL + cL, uR + cR)

    # handle bad ordering
    mask_bad = SL >= SR

    # middle wave speed S_M
    denom = (rhoL * (SL - uL) - rhoR * (SR - uR))
    S_M = np.zeros_like(denom)
    safe = denom != 0.0
    if np.any(safe):
        S_M[safe] = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR))[safe] / denom[safe]

    # E_spec per unit mass
    EL = UL[2] / rhoL
    ER = UR[2] / rhoR
    # avoid division by zero for star densities
    rhoSL = np.where((SL - S_M) != 0.0, rhoL * (SL - uL) / (SL - S_M), 0.0)
    rhoSR = np.where((SR - S_M) != 0.0, rhoR * (SR - uR) / (SR - S_M), 0.0)
    UstarL = np.zeros_like(UL)
    UstarR = np.zeros_like(UR)
    UstarL[0] = rhoSL
    UstarL[1] = rhoSL * S_M
    UstarL[2] = rhoSL * (EL + (S_M - uL) * (S_M + pL / (rhoL * (SL - uL))))
    UstarR[0] = rhoSR
    UstarR[1] = rhoSR * S_M
    UstarR[2] = rhoSR * (ER + (S_M - uR) * (S_M + pR / (rhoR * (SR - uR))))

    # allocate flux
    F = np.zeros_like(FL)

    # Regions
    mask_FL = SL >= 0.0
    if np.any(mask_FL):
        F[:, mask_FL] = FL[:, mask_FL]
    mask_Lstar = (SL < 0.0) & (S_M >= 0.0)
    if np.any(mask_Lstar):
        F[:, mask_Lstar] = FL[:, mask_Lstar] + (SL[mask_Lstar] * (UstarL[:, mask_Lstar] - UL[:, mask_Lstar]))
    mask_Rstar = (S_M < 0.0) & (SR >= 0.0)
    if np.any(mask_Rstar):
        F[:, mask_Rstar] = FR[:, mask_Rstar] + (SR[mask_Rstar] * (UstarR[:, mask_Rstar] - UR[:, mask_Rstar]))
    mask_FR = SR < 0.0
    if np.any(mask_FR):
        F[:, mask_FR] = FR[:, mask_FR]

    # fallback to HLL where SL >= SR
    if np.any(mask_bad):
        denom_hll = (SR - SL)
        denom_hll[denom_hll == 0.0] = 1e-16
        F_hll = (SR * FL - SL * FR + SR * SL * (UR - UL)) / denom_hll
        F[:, mask_bad] = F_hll[:, mask_bad]

    return F

# ----------------------------- RHS builder -------------------------------

def apply_transmissive_BC(U, ng):
    """Apply transmissive (zero-gradient) BCs by copying nearest interior cells into ghosts."""
    # left ghosts <- first interior cells
    U[:, :ng] = U[:, ng:2 * ng].copy()
    # right ghosts <- last interior cells
    U[:, -ng:] = U[:, -2 * ng:-ng].copy()


def compute_RHS(U, dx, gamma, visc,
                prim_arr, slope, qL, qR, flux_interface, rhs_out):
    """Compute RHS = -dF/dx for interior cells. Uses preallocated arrays passed in.
    Arrays shapes:
      U: (3, Ntot)
      prim_arr: (3, Ntot)
      slope: (3, Ntot)
      qL, qR: (3, Ninterfaces)   (interfaces are between cells)
      flux_interface: (3, Ninterfaces)
      rhs_out: (3, Ntot)
    """
    Ntot = U.shape[1]
    ng = 2

    # boundary conditions: transmissive (zero-gradient)
    apply_transmissive_BC(U, ng)

    # primitives
    rho, u, p = cons_to_prim(U, gamma)
    prim_arr[0, :] = rho
    prim_arr[1, :] = u
    prim_arr[2, :] = p

    # limited slopes (minmod3) on primitive variables
    theta = 1.5
    dqL = prim_arr[:, 1:-1] - prim_arr[:, 0:-2]
    dqR = prim_arr[:, 2:] - prim_arr[:, 1:-1]
    dqC = 0.5 * (prim_arr[:, 2:] - prim_arr[:, 0:-2])
    slope[:, 1:-1] = minmod3(theta * dqL, dqC, theta * dqR)
    # set slopes at ghosts equal to nearest interior slope (simple)
    slope[:, 0] = slope[:, 1]
    slope[:, -1] = slope[:, -2]

    # reconstruct at interfaces: between cell i and i+1 we take left from i and right from i+1
    # number of interfaces = Ntot - 1
    Nint = Ntot - 1
    assert qL.shape[1] == Nint and qR.shape[1] == Nint and flux_interface.shape[1] == Nint, \
        "Interface arrays must have length Ntot-1"

    # qL at interface i+1/2 comes from cell i: q_i + 0.5*slope_i
    qL[:, :] = prim_arr[:, :-1] + 0.5 * slope[:, :-1]
    # qR at interface i+1/2 comes from cell i+1: q_{i+1} - 0.5*slope_{i+1}
    qR[:, :] = prim_arr[:, 1:] - 0.5 * slope[:, 1:]

    # Ensure positive density and pressure in reconstructed states (simple floor)
    eps = 1e-12
    qL[0, :] = np.maximum(qL[0, :], eps)
    qR[0, :] = np.maximum(qR[0, :], eps)
    qL[2, :] = np.maximum(qL[2, :], eps)
    qR[2, :] = np.maximum(qR[2, :], eps)

    # Riemann fluxes at interfaces (HLLC)
    flux_interface[:, :] = hllc_flux((qL[0, :], qL[1, :], qL[2, :]),
                                     (qR[0, :], qR[1, :], qR[2, :]), gamma)

    # Add small diffusive (viscous) flux to stabilize if visc > 0
    if visc > 0.0:
        cL = np.sqrt(np.maximum(gamma * qL[2, :] / qL[0, :], 0.0))
        cR = np.sqrt(np.maximum(gamma * qR[2, :] / qR[0, :], 0.0))
        smax = np.maximum(np.abs(qL[1, :]) + cL, np.abs(qR[1, :]) + cR)
        flux_interface += (-0.5 * visc) * (smax[np.newaxis, :] * (prim_to_cons(qR[0, :], qR[1, :], qR[2, :], gamma) -
                                                                    prim_to_cons(qL[0, :], qL[1, :], qL[2, :], gamma)))

    # compute RHS for interior cells (exclude ghost cells)
    rhs_out[:, :] = 0.0
    # Correct divergence indexing: flux_interface[j] corresponds to F_{j+1/2}; cell i uses F_{i+1/2} - F_{i-1/2}
    # therefore rhs[:, i] = - (flux_interface[:, i] - flux_interface[:, i-1]) / dx. For interior cells i=1..Ntot-2 -> use indices 1:-1
    rhs_out[:, 1:-1] = - (flux_interface[:, 1:] - flux_interface[:, :-1]) / dx
    return rhs_out

# ----------------------------- Main solver --------------------------------

def solve_sod(nx=400, x_left=0.0, x_right=1.0,
              gamma=1.4, visc=0.0, CFL=0.45,
              final_times=(0.1, 0.2, 0.3)):
    # grid
    Nx = nx
    ng = 2  # ghost cells
    Ntot = Nx + 2 * ng
    dx = (x_right - x_left) / Nx
    xcell = np.linspace(x_left + 0.5 * dx, x_right - 0.5 * dx, Nx)

    # allocate arrays (reuse to avoid allocations inside time loop)
    U = np.zeros((3, Ntot))
    prim = np.zeros((3, Ntot))
    slope = np.zeros((3, Ntot))
    qL = np.zeros((3, Ntot - 1))
    qR = np.zeros((3, Ntot - 1))
    flux_if = np.zeros((3, Ntot - 1))
    rhs = np.zeros((3, Ntot))

    # initialize cell-centered interior primitive states (Sod tube)
    i0 = ng
    i1 = ng + Nx
    rho0 = np.where(xcell <= 0.5, 1.0, 0.125)
    u0 = np.zeros_like(xcell)
    p0 = np.where(xcell <= 0.5, 1.0, 0.1)
    U[:, i0:i1] = prim_to_cons(rho0, u0, p0, gamma)

    # fill ghost cells with nearest interior (transmissive)
    apply_transmissive_BC(U, ng)

    # quick unit-check for BC behavior (simple test)
    assert np.allclose(U[:, :ng], U[:, ng:2 * ng]), "Transmissive BC left ghosts incorrect"
    assert np.allclose(U[:, -ng:], U[:, -2 * ng:-ng]), "Transmissive BC right ghosts incorrect"

    # time stepping
    t = 0.0
    Tfinal_all = sorted(list(final_times))
    saved_flags = {t_out: False for t_out in Tfinal_all}
    save_rho = {}
    save_u = {}
    save_p = {}

    Tfinal = max(final_times)
    max_steps = 500000
    step = 0
    tol = 1e-12
    min_dt = 1e-14

    # cumulative diagnostic for pressure repairs
    total_pressure_repairs = 0

    while t < Tfinal - 1e-14 and step < max_steps:
        # compute primitive for dt estimation
        rho, u, p = cons_to_prim(U, gamma)
        c = np.sqrt(np.maximum(gamma * p / rho, 0.0))
        max_wave = np.max(np.abs(u) + c)
        if max_wave <= 0.0:
            dt = 1e-6
        else:
            dt = CFL * dx / max_wave

        # if stepping would overshoot the next unsaved snapshot, reduce dt
        unsaved = [tt for tt in Tfinal_all if not saved_flags[tt] and tt >= t]
        next_unsaved = min(unsaved) if unsaved else None
        if next_unsaved is not None and t + dt > next_unsaved:
            dt = next_unsaved - t

        # Guard against zero (or negative) dt due to rounding
        if dt <= tol:
            # If we're at or extremely close to a snapshot time, advance time to that snapshot and save without RK update
            if next_unsaved is not None and next_unsaved >= t - tol:
                t = next_unsaved
                # save snapshot from current U (no update)
                rho_i, u_i, p_i = cons_to_prim(U[:, ng:-ng], gamma)
                save_rho[t] = rho_i.copy()
                save_u[t] = u_i.copy()
                save_p[t] = p_i.copy()
                # diagnostics
                drho = np.abs(np.diff(rho_i))
                max_idx = np.argmax(drho)
                x_shock = 0.5 * (xcell[max_idx] + xcell[max_idx + 1])
                mass = np.sum(rho_i) * dx
                momentum = np.sum(rho_i * u_i) * dx
                energy = np.sum(U[2, ng:-ng]) * dx
                min_rho = np.min(rho_i)
                min_p = np.min(p_i)
                TV_rho = np.sum(np.abs(np.diff(rho_i)))
                print(f"\nSnapshot at t = {t:.6f} (dt ~ 0, saved without RK)")
                print(f"  Shock position (est.): x = {x_shock:.6f}")
                print(f"  Mass = {mass:.6f}, Momentum = {momentum:.6f}, Energy = {energy:.6f}")
                print(f"  min(rho) = {min_rho:.6e}, min(p) = {min_p:.6e}, TV(rho) = {TV_rho:.6e}")
                saved_flags[t] = True
                # continue to next loop iteration
                continue
            else:
                # fallback small positive dt to make some progress
                dt = max(min_dt, tol)

        # SSPRK3 stages
        Un = U.copy()

        # stage 1
        compute_RHS(U, dx, gamma, visc, prim, slope, qL, qR, flux_if, rhs)
        U1 = Un + dt * rhs

        # stage 2
        compute_RHS(U1, dx, gamma, visc, prim, slope, qL, qR, flux_if, rhs)
        U2 = 0.75 * Un + 0.25 * (U1 + dt * rhs)

        # stage 3
        compute_RHS(U2, dx, gamma, visc, prim, slope, qL, qR, flux_if, rhs)
        U = (1.0 / 3.0) * Un + (2.0 / 3.0) * (U2 + dt * rhs)

        t += dt
        step += 1

        # ensure physical admissibility (vectorized, conservative energy repair)
        rho_new = U[0, :]
        mom_new = U[1, :]
        E_new = U[2, :]
        eps = 1e-14
        neg_rho = rho_new < eps
        if np.any(neg_rho):
            rho_new[neg_rho] = eps
            U[0, :] = rho_new
            # momentum and energy left as-is; this is a strong fix

        # vectorized pressure check and conservative repair
        rho_tmp = U[0, :]
        u_tmp = np.zeros_like(rho_tmp)
        mask = rho_tmp > 0.0
        u_tmp[mask] = U[1, mask] / rho_tmp[mask]
        E_spec = U[2, :] / rho_tmp - 0.5 * u_tmp * u_tmp
        p_tmp = (gamma - 1.0) * rho_tmp * E_spec
        neg_p = p_tmp < 1e-14
        if np.any(neg_p):
            # conservative repair: enforce minimum specific internal energy corresponding to pressure floor
            p_floor = 1e-14
            e_floor = p_floor / ((gamma - 1.0) * rho_tmp[neg_p])
            E_spec_new = E_spec.copy()
            E_spec_new[neg_p] = np.maximum(E_spec_new[neg_p], e_floor)
            U[2, :] = rho_tmp * (0.5 * u_tmp * u_tmp + E_spec_new)
            nfix = np.sum(neg_p)
            total_pressure_repairs += int(nfix)
            print(f"  Pressure repairs at t={t:.6f}: {int(nfix)} cells (cumulative {total_pressure_repairs})")

        # after advance check for snapshots to save (use tolerance-based condition)
        for t_out in Tfinal_all:
            if (not saved_flags[t_out]) and (t >= t_out - tol):
                rho_i, u_i, p_i = cons_to_prim(U[:, ng:-ng], gamma)
                save_rho[t_out] = rho_i.copy()
                save_u[t_out] = u_i.copy()
                save_p[t_out] = p_i.copy()
                # diagnostics
                drho = np.abs(np.diff(rho_i))
                max_idx = np.argmax(drho)
                x_shock = 0.5 * (xcell[max_idx] + xcell[max_idx + 1])
                mass = np.sum(rho_i) * dx
                momentum = np.sum((rho_i * u_i)) * dx
                energy = np.sum(U[2, ng:-ng]) * dx
                min_rho = np.min(rho_i)
                min_p = np.min(p_i)
                TV_rho = np.sum(np.abs(np.diff(rho_i)))
                print(f"\nSnapshot at t = {t_out:.6f}")
                print(f"  Shock position (est.): x = {x_shock:.6f}")
                print(f"  Mass = {mass:.6f}, Momentum = {momentum:.6f}, Energy = {energy:.6f}")
                print(f"  min(rho) = {min_rho:.6e}, min(p) = {min_p:.6e}, TV(rho) = {TV_rho:.6e}")
                saved_flags[t_out] = True

    if step >= max_steps:
        print("Warning: reached max steps limit")

    # ensure all requested snapshots saved
    for t_out in Tfinal_all:
        if not saved_flags[t_out]:
            print(f"Warning: requested snapshot at t={t_out} was not saved")

    # plotting: overlay three times into three separate figures
    times = Tfinal_all

    # density
    plt.figure(figsize=(8, 4))
    for t_out in times:
        rho_plot = save_rho.get(t_out)
        if rho_plot is not None:
            plt.plot(xcell, rho_plot, label=f't={t_out:.2f}')
    plt.xlabel('x')
    plt.ylabel('Density')
    plt.title('Density at different times (MUSCL + HLLC)')
    plt.legend()
    plt.grid(True)
    plt.savefig('density_times.png', dpi=200)
    plt.close()

    # velocity
    plt.figure(figsize=(8, 4))
    for t_out in times:
        u_plot = save_u.get(t_out)
        if u_plot is not None:
            plt.plot(xcell, u_plot, label=f't={t_out:.2f}')
    plt.xlabel('x')
    plt.ylabel('Velocity')
    plt.title('Velocity at different times (MUSCL + HLLC)')
    plt.legend()
    plt.grid(True)
    plt.savefig('velocity_times.png', dpi=200)
    plt.close()

    # pressure
    plt.figure(figsize=(8, 4))
    for t_out in times:
        p_plot = save_p.get(t_out)
        if p_plot is not None:
            plt.plot(xcell, p_plot, label=f't={t_out:.2f}')
    plt.xlabel('x')
    plt.ylabel('Pressure')
    plt.title('Pressure at different times (MUSCL + HLLC)')
    plt.legend()
    plt.grid(True)
    plt.savefig('pressure_times.png', dpi=200)
    plt.close()

    print('\nSaved figures: density_times.png, velocity_times.png, pressure_times.png')
    print(f'Total pressure repairs performed: {total_pressure_repairs}')
    return save_rho, save_u, save_p


if __name__ == '__main__':
    # Explicit constants passed to solver
    gamma_val = 1.4
    visc_val = 0.0  # small positive value (e.g., 0.1) increases dissipation
    CFL_val = 0.45
    Nx = 400
    final_times_tuple = (0.1, 0.2, 0.3)

    solve_sod(nx=Nx, x_left=0.0, x_right=1.0,
              gamma=gamma_val, visc=visc_val, CFL=CFL_val,
              final_times=final_times_tuple)
```


#### Script block3:

```python
#!/usr/bin/env python3
"""
Revised MUSCL finite-volume solver for 1D Euler equations (Sod shock tube)
- Second-order MUSCL reconstruction (primitive variables) with minmod limiter
- HLLC Riemann solver at interfaces
- SSPRK3 time integration
- Optional interface viscous dissipation controlled by 'visc'

Key fixes from review:
- Corrected discrete divergence indexing (uses ghost count 'ng')
- Enforce transmissive BCs after each RK substage (prevents ghost corruption)
- apply_transmissive_BC now replicates nearest interior cell into all ghost slots
- Conservative density-floor repair: momentum and energy rescaled consistently
- Improved shock-location reporting: both max |drho| and right-most large positive jump

Outputs: three saved figures showing density, velocity, pressure at t=0.1,0.2,0.3
"""

import numpy as np
import matplotlib.pyplot as plt

# --------------------------- Utility converters ---------------------------

def prim_to_cons(rho, u, p, gamma):
    """Convert primitive vars to conservative vector U = [rho, rho*u, rho*E]."""
    mom = rho * u
    E_spec = 0.5 * u * u + p / ((gamma - 1.0) * rho)
    Econs = rho * E_spec
    U = np.vstack((rho, mom, Econs))
    return U


def cons_to_prim(U, gamma):
    """Convert conservative U (3 x N) to primitive rho,u,p arrays."""
    rho = U[0].copy()
    mom = U[1].copy()
    Econs = U[2].copy()
    u = np.zeros_like(rho)
    mask = rho > 0.0
    u[mask] = mom[mask] / rho[mask]
    kinetic = 0.5 * mom * u
    p = (gamma - 1.0) * (Econs - kinetic)
    return rho, u, p


def flux_from_prim(rho, u, p, gamma):
    """Compute physical flux vector F = [rho*u, rho*u^2 + p, (rho*E + p)*u]."""
    mom = rho * u
    E_spec = 0.5 * u * u + p / ((gamma - 1.0) * rho)
    Econs = rho * E_spec
    F0 = mom
    F1 = mom * u + p
    F2 = (Econs + p) * u
    return np.vstack((F0, F1, F2))


# ----------------------------- Limiters ----------------------------------

def minmod3(a, b, c):
    """Elementwise three-argument minmod limiter."""
    result = np.zeros_like(a)
    s = np.sign(a) + np.sign(b) + np.sign(c)
    mask_pos = (s == 3)
    mask_neg = (s == -3)
    if np.any(mask_pos):
        result[mask_pos] = np.minimum(np.minimum(np.abs(a[mask_pos]), np.abs(b[mask_pos])), np.abs(c[mask_pos]))
    if np.any(mask_neg):
        result[mask_neg] = -np.minimum(np.minimum(np.abs(a[mask_neg]), np.abs(b[mask_neg])), np.abs(c[mask_neg]))
    return result


# ----------------------------- HLLC Solver -------------------------------

def hllc_flux(primL, primR, gamma):
    """Compute HLLC flux for left/right primitive states at interfaces.
    primL, primR are tuples (rho,u,p) where each element can be arrays.
    Returns flux array shape (3, Ninterfaces).
    """
    rhoL, uL, pL = primL
    rhoR, uR, pR = primR
    # sound speeds
    cL = np.sqrt(np.maximum(gamma * pL / rhoL, 0.0))
    cR = np.sqrt(np.maximum(gamma * pR / rhoR, 0.0))
    # conservative states
    UL = prim_to_cons(rhoL, uL, pL, gamma)
    UR = prim_to_cons(rhoR, uR, pR, gamma)
    # fluxes
    FL = flux_from_prim(rhoL, uL, pL, gamma)
    FR = flux_from_prim(rhoR, uR, pR, gamma)
    # wave speed estimates
    SL = np.minimum(uL - cL, uR - cR)
    SR = np.maximum(uL + cL, uR + cR)

    # handle bad ordering
    mask_bad = SL >= SR

    # middle wave speed S_M
    denom = (rhoL * (SL - uL) - rhoR * (SR - uR))
    S_M = np.zeros_like(denom)
    safe = denom != 0.0
    if np.any(safe):
        S_M[safe] = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR))[safe] / denom[safe]

    # E_spec per unit mass
    EL = UL[2] / rhoL
    ER = UR[2] / rhoR
    # avoid division by zero for star densities
    rhoSL = np.where((SL - S_M) != 0.0, rhoL * (SL - uL) / (SL - S_M), 0.0)
    rhoSR = np.where((SR - S_M) != 0.0, rhoR * (SR - uR) / (SR - S_M), 0.0)
    UstarL = np.zeros_like(UL)
    UstarR = np.zeros_like(UR)
    UstarL[0] = rhoSL
    UstarL[1] = rhoSL * S_M
    UstarL[2] = rhoSL * (EL + (S_M - uL) * (S_M + pL / (rhoL * (SL - uL))))
    UstarR[0] = rhoSR
    UstarR[1] = rhoSR * S_M
    UstarR[2] = rhoSR * (ER + (S_M - uR) * (S_M + pR / (rhoR * (SR - uR))))

    # allocate flux
    F = np.zeros_like(FL)

    # Regions
    mask_FL = SL >= 0.0
    if np.any(mask_FL):
        F[:, mask_FL] = FL[:, mask_FL]
    mask_Lstar = (SL < 0.0) & (S_M >= 0.0)
    if np.any(mask_Lstar):
        F[:, mask_Lstar] = FL[:, mask_Lstar] + (SL[mask_Lstar] * (UstarL[:, mask_Lstar] - UL[:, mask_Lstar]))
    mask_Rstar = (S_M < 0.0) & (SR >= 0.0)
    if np.any(mask_Rstar):
        F[:, mask_Rstar] = FR[:, mask_Rstar] + (SR[mask_Rstar] * (UstarR[:, mask_Rstar] - UR[:, mask_Rstar]))
    mask_FR = SR < 0.0
    if np.any(mask_FR):
        F[:, mask_FR] = FR[:, mask_FR]

    # fallback to HLL where SL >= SR
    if np.any(mask_bad):
        denom_hll = (SR - SL)
        denom_hll[denom_hll == 0.0] = 1e-16
        F_hll = (SR * FL - SL * FR + SR * SL * (UR - UL)) / denom_hll
        F[:, mask_bad] = F_hll[:, mask_bad]

    return F


# ----------------------------- RHS builder -------------------------------

def apply_transmissive_BC(U, ng):
    """Apply transmissive (zero-gradient) BCs by copying the nearest interior cell
    into all ghost slots (vectorized).
    """
    # left ghosts <- nearest interior (index ng)
    U[:, :ng] = U[:, ng:ng+1]
    # right ghosts <- nearest interior (index -ng-1)
    U[:, -ng:] = U[:, -ng-1:-ng]


def compute_RHS(U, dx, gamma, visc,
                prim_arr, slope, qL, qR, flux_interface, rhs_out):
    """Compute RHS = -dF/dx for interior (physical) cells. Uses preallocated arrays.
    Arrays shapes:
      U: (3, Ntot)
      prim_arr: (3, Ntot)
      slope: (3, Ntot)
      qL, qR: (3, Ninterfaces)   (interfaces are between cells)
      flux_interface: (3, Ninterfaces)
      rhs_out: (3, Ntot)
    """
    Ntot = U.shape[1]
    ng = 2

    # boundary conditions: transmissive (zero-gradient)
    apply_transmissive_BC(U, ng)

    # primitives
    rho, u, p = cons_to_prim(U, gamma)
    prim_arr[0, :] = rho
    prim_arr[1, :] = u
    prim_arr[2, :] = p

    # limited slopes (minmod3) on primitive variables
    theta = 1.5
    dqL = prim_arr[:, 1:-1] - prim_arr[:, 0:-2]
    dqR = prim_arr[:, 2:] - prim_arr[:, 1:-1]
    dqC = 0.5 * (prim_arr[:, 2:] - prim_arr[:, 0:-2])
    slope[:, 1:-1] = minmod3(theta * dqL, dqC, theta * dqR)
    # set slopes at ghosts equal to nearest interior slope
    slope[:, 0] = slope[:, 1]
    slope[:, -1] = slope[:, -2]

    # reconstruct at interfaces: between cell i and i+1 we take left from i and right from i+1
    # number of interfaces = Ntot - 1
    Nint = Ntot - 1
    assert qL.shape[1] == Nint and qR.shape[1] == Nint and flux_interface.shape[1] == Nint, \
        "Interface arrays must have length Ntot-1"

    # qL at interface i+1/2 comes from cell i: q_i + 0.5*slope_i
    qL[:, :] = prim_arr[:, :-1] + 0.5 * slope[:, :-1]
    # qR at interface i+1/2 comes from cell i+1: q_{i+1} - 0.5*slope_{i+1}
    qR[:, :] = prim_arr[:, 1:] - 0.5 * slope[:, 1:]

    # Ensure positive density and pressure in reconstructed states (simple floor)
    eps = 1e-12
    qL[0, :] = np.maximum(qL[0, :], eps)
    qR[0, :] = np.maximum(qR[0, :], eps)
    qL[2, :] = np.maximum(qL[2, :], eps)
    qR[2, :] = np.maximum(qR[2, :], eps)

    # Riemann fluxes at interfaces (HLLC)
    flux_interface[:, :] = hllc_flux((qL[0, :], qL[1, :], qL[2, :]),
                                     (qR[0, :], qR[1, :], qR[2, :]), gamma)

    # Add small diffusive (viscous) flux to stabilize if visc > 0
    if visc > 0.0:
        cL = np.sqrt(np.maximum(gamma * qL[2, :] / qL[0, :], 0.0))
        cR = np.sqrt(np.maximum(gamma * qR[2, :] / qR[0, :], 0.0))
        smax = np.maximum(np.abs(qL[1, :]) + cL, np.abs(qR[1, :]) + cR)
        flux_interface += (-0.5 * visc) * (smax[np.newaxis, :] * (prim_to_cons(qR[0, :], qR[1, :], qR[2, :], gamma) -
                                                                    prim_to_cons(qL[0, :], qL[1, :], qL[2, :], gamma)))

    # compute RHS for physical interior cells (exclude ghost cells)
    rhs_out[:, :] = 0.0
    # Correct divergence indexing using ng offset
    Nphys = Ntot - 2 * ng
    # flux_interface[:, j] = F_{j+1/2}; for cell i the divergence is (F_{i+1/2} - F_{i-1/2})
    # for physical cells i in [ng .. Ntot-ng-1] the corresponding flux_interface indices run
    # i -> flux_interface index i, and i-1 -> flux_interface index i-1. Thus:
    rhs_out[:, ng:-ng] = - (flux_interface[:, ng:ng + Nphys] - flux_interface[:, ng - 1:ng - 1 + Nphys]) / dx
    return rhs_out


# ----------------------------- Main solver --------------------------------

def solve_sod(nx=400, x_left=0.0, x_right=1.0,
              gamma=1.4, visc=0.0, CFL=0.45,
              final_times=(0.1, 0.2, 0.3)):
    # grid
    Nx = nx
    ng = 2  # ghost cells (kept consistent)
    Ntot = Nx + 2 * ng
    dx = (x_right - x_left) / Nx
    xcell = np.linspace(x_left + 0.5 * dx, x_right - 0.5 * dx, Nx)

    # allocate arrays (reuse to avoid allocations inside time loop)
    U = np.zeros((3, Ntot))
    prim = np.zeros((3, Ntot))
    slope = np.zeros((3, Ntot))
    qL = np.zeros((3, Ntot - 1))
    qR = np.zeros((3, Ntot - 1))
    flux_if = np.zeros((3, Ntot - 1))
    rhs = np.zeros((3, Ntot))

    # initialize cell-centered interior primitive states (Sod tube)
    i0 = ng
    i1 = ng + Nx
    rho0 = np.where(xcell <= 0.5, 1.0, 0.125)
    u0 = np.zeros_like(xcell)
    p0 = np.where(xcell <= 0.5, 1.0, 0.1)
    U[:, i0:i1] = prim_to_cons(rho0, u0, p0, gamma)

    # fill ghost cells with nearest interior (transmissive)
    apply_transmissive_BC(U, ng)

    # quick unit-check for BC behavior
    assert np.allclose(U[:, :ng], U[:, ng:ng+1].repeat(ng, axis=1)), "Transmissive BC left ghosts incorrect"
    assert np.allclose(U[:, -ng:], U[:, -ng-1:-ng].repeat(ng, axis=1)), "Transmissive BC right ghosts incorrect"

    # time stepping
    t = 0.0
    Tfinal_all = sorted(list(final_times))
    saved_flags = {t_out: False for t_out in Tfinal_all}
    save_rho = {}
    save_u = {}
    save_p = {}

    Tfinal = max(final_times)
    max_steps = 500000
    step = 0
    tol = 1e-12
    min_dt = 1e-14

    # cumulative diagnostic for pressure/density repairs
    total_pressure_repairs = 0
    total_density_repairs = 0

    # helper: improved shock locator and reporting
    def report_snapshot(t_report, U_in):
        rho_i, u_i, p_i = cons_to_prim(U_in[:, ng:-ng], gamma)
        # primary shock locator: max |drho|
        drho = np.diff(rho_i)
        absdrho = np.abs(drho)
        max_idx = int(np.argmax(absdrho))
        # secondary: right-most large positive jump (shock tends to be positive density jump)
        thresh = 0.10 * (np.max(absdrho) if np.max(absdrho) > 0 else 1.0)
        candidates = np.where(drho > thresh)[0]
        if candidates.size > 0:
            rightmost_pos_idx = int(candidates[-1])
        else:
            rightmost_pos_idx = max_idx
        x_shock_maxabs = 0.5 * (xcell[max_idx] + xcell[max_idx + 1])
        x_shock_rightpos = 0.5 * (xcell[rightmost_pos_idx] + xcell[rightmost_pos_idx + 1])
        mass = np.sum(rho_i) * dx
        momentum = np.sum(rho_i * u_i) * dx
        energy = np.sum(U_in[2, ng:-ng]) * dx
        min_rho = np.min(rho_i)
        min_p = np.min(p_i)
        TV_rho = np.sum(np.abs(np.diff(rho_i)))
        print(f"\nSnapshot at t = {t_report:.6f}")
        print(f"  Shock pos (max|drho|): x = {x_shock_maxabs:.6f}, right-most +jump: x = {x_shock_rightpos:.6f}")
        print(f"  Mass = {mass:.6f}, Momentum = {momentum:.6f}, Energy = {energy:.6f}")
        print(f"  min(rho) = {min_rho:.6e}, min(p) = {min_p:.6e}, TV(rho) = {TV_rho:.6e}")
        return rho_i, u_i, p_i

    # main time loop
    while t < Tfinal - 1e-14 and step < max_steps:
        # compute primitive for dt estimation
        rho, u, p = cons_to_prim(U, gamma)
        c = np.sqrt(np.maximum(gamma * p / rho, 0.0))
        max_wave = np.max(np.abs(u) + c)
        if max_wave <= 0.0:
            dt = 1e-6
        else:
            dt = CFL * dx / max_wave

        # if stepping would overshoot the next unsaved snapshot, reduce dt
        unsaved = [tt for tt in Tfinal_all if not saved_flags[tt] and tt >= t]
        next_unsaved = min(unsaved) if unsaved else None
        if next_unsaved is not None and t + dt > next_unsaved:
            dt = next_unsaved - t

        # Guard against zero (or negative) dt due to rounding
        if dt <= tol:
            # If we're at or extremely close to a snapshot time, advance time to that snapshot and save without RK update
            if next_unsaved is not None and next_unsaved >= t - tol:
                t = next_unsaved
                rho_i, u_i, p_i = report_snapshot(t, U)
                save_rho[t] = rho_i.copy()
                save_u[t] = u_i.copy()
                save_p[t] = p_i.copy()
                saved_flags[t] = True
                continue
            else:
                dt = max(min_dt, tol)

        # SSPRK3 stages (ensure BCs enforced after each stage)
        Un = U.copy()

        # stage 1
        compute_RHS(U, dx, gamma, visc, prim, slope, qL, qR, flux_if, rhs)
        U1 = Un + dt * rhs
        apply_transmissive_BC(U1, ng)

        # stage 2
        compute_RHS(U1, dx, gamma, visc, prim, slope, qL, qR, flux_if, rhs)
        U2 = 0.75 * Un + 0.25 * (U1 + dt * rhs)
        apply_transmissive_BC(U2, ng)

        # stage 3
        compute_RHS(U2, dx, gamma, visc, prim, slope, qL, qR, flux_if, rhs)
        U = (1.0 / 3.0) * Un + (2.0 / 3.0) * (U2 + dt * rhs)
        # enforce BCs on updated solution prior to diagnostics
        apply_transmissive_BC(U, ng)

        t += dt
        step += 1

        # ensure physical admissibility (vectorized, conservative energy & density repair)
        rho_new = U[0, :].copy()
        mom_new = U[1, :].copy()
        E_new = U[2, :].copy()

        # density floor: conservative scaling of momentum and total energy
        rho_floor = 1e-8
        small_mask = rho_new < rho_floor
        if np.any(small_mask):
            rho_before = U[0, :].copy()
            mom_before = U[1, :].copy()
            E_before = U[2, :].copy()
            # compute old velocity and specific internal energy safely
            u_old = np.zeros_like(rho_before)
            mask_pos = rho_before > 0.0
            u_old[mask_pos] = mom_before[mask_pos] / rho_before[mask_pos]
            E_spec_old = np.zeros_like(rho_before)
            E_spec_old[mask_pos] = E_before[mask_pos] / rho_before[mask_pos] - 0.5 * u_old[mask_pos] * u_old[mask_pos]
            # for zero-rho entries assign minimal E_spec
            if np.any(~mask_pos):
                E_spec_old[~mask_pos] = 1e-12
                u_old[~mask_pos] = 0.0
            # apply floors and rescale conservatively
            rho_new[small_mask] = rho_floor
            mom_new[small_mask] = u_old[small_mask] * rho_new[small_mask]
            E_new[small_mask] = rho_new[small_mask] * (0.5 * u_old[small_mask] * u_old[small_mask] + E_spec_old[small_mask])
            U[0, :] = rho_new
            U[1, :] = mom_new
            U[2, :] = E_new
            nfix = int(np.sum(small_mask))
            total_density_repairs += nfix
            print(f"  Density repairs at t={t:.6f}: {nfix} cells (cumulative {total_density_repairs})")

        # pressure floor: conservative repair already implemented (adjust total energy)
        rho_tmp = U[0, :]
        u_tmp = np.zeros_like(rho_tmp)
        mask = rho_tmp > 0.0
        u_tmp[mask] = U[1, mask] / rho_tmp[mask]
        E_spec = U[2, :] / rho_tmp - 0.5 * u_tmp * u_tmp
        p_tmp = (gamma - 1.0) * rho_tmp * E_spec
        p_floor = 1e-10
        neg_p = p_tmp < p_floor
        if np.any(neg_p):
            p_floor_local = p_floor
            e_floor = p_floor_local / ((gamma - 1.0) * rho_tmp[neg_p])
            E_spec_new = E_spec.copy()
            E_spec_new[neg_p] = np.maximum(E_spec_new[neg_p], e_floor)
            U[2, :] = rho_tmp * (0.5 * u_tmp * u_tmp + E_spec_new)
            nfix = int(np.sum(neg_p))
            total_pressure_repairs += nfix
            print(f"  Pressure repairs at t={t:.6f}: {nfix} cells (cumulative {total_pressure_repairs})")

        # after advance check for snapshots to save (use tolerance-based condition)
        for t_out in Tfinal_all:
            if (not saved_flags[t_out]) and (t >= t_out - tol):
                rho_i, u_i, p_i = report_snapshot(t_out, U)
                save_rho[t_out] = rho_i.copy()
                save_u[t_out] = u_i.copy()
                save_p[t_out] = p_i.copy()
                saved_flags[t_out] = True

    if step >= max_steps:
        print("Warning: reached max steps limit")

    # ensure all requested snapshots saved
    for t_out in Tfinal_all:
        if not saved_flags[t_out]:
            print(f"Warning: requested snapshot at t={t_out} was not saved")

    # plotting: overlay three times into three separate figures
    times = Tfinal_all

    # density
    plt.figure(figsize=(8, 4))
    for t_out in times:
        rho_plot = save_rho.get(t_out)
        if rho_plot is not None:
            plt.plot(xcell, rho_plot, label=f't={t_out:.2f}')
    plt.xlabel('x')
    plt.ylabel('Density')
    plt.title('Density at different times (MUSCL + HLLC)')
    plt.legend()
    plt.grid(True)
    plt.savefig('density_times.png', dpi=200)
    plt.close()

    # velocity
    plt.figure(figsize=(8, 4))
    for t_out in times:
        u_plot = save_u.get(t_out)
        if u_plot is not None:
            plt.plot(xcell, u_plot, label=f't={t_out:.2f}')
    plt.xlabel('x')
    plt.ylabel('Velocity')
    plt.title('Velocity at different times (MUSCL + HLLC)')
    plt.legend()
    plt.grid(True)
    plt.savefig('velocity_times.png', dpi=200)
    plt.close()

    # pressure
    plt.figure(figsize=(8, 4))
    for t_out in times:
        p_plot = save_p.get(t_out)
        if p_plot is not None:
            plt.plot(xcell, p_plot, label=f't={t_out:.2f}')
    plt.xlabel('x')
    plt.ylabel('Pressure')
    plt.title('Pressure at different times (MUSCL + HLLC)')
    plt.legend()
    plt.grid(True)
    plt.savefig('pressure_times.png', dpi=200)
    plt.close()

    print('\nSaved figures: density_times.png, velocity_times.png, pressure_times.png')
    print(f'Total pressure repairs performed: {total_pressure_repairs}')
    print(f'Total density repairs performed: {total_density_repairs}')
    return save_rho, save_u, save_p


if __name__ == '__main__':
    # Explicit constants passed to solver
    gamma_val = 1.4
    visc_val = 0.0  # small positive value (e.g., 0.1) increases dissipation
    CFL_val = 0.45
    Nx = 400
    final_times_tuple = (0.1, 0.2, 0.3)

    solve_sod(nx=Nx, x_left=0.0, x_right=1.0,
              gamma=gamma_val, visc=visc_val, CFL=CFL_val,
              final_times=final_times_tuple)
```




### runtime_outputs

#### Output block1


Runtime Error:
Traceback (most recent call last):
  File "C:\Users\admin\AppData\Local\Temp\ipykernel_45176\2663542217.py", line 138, in execute_code_tool
    exec(code, global_scope, global_scope)
  File "<string>", line 433, in <module>
  File "<string>", line 298, in solve_sod
  File "<string>", line 220, in compute_RHS
ValueError: could not broadcast input array from shape (3,401) into shape (3,402)



#### Output block2


Snapshot at t = 0.100000
  Shock position (est.): x = 0.675000
  Mass = 0.562500, Momentum = 0.090000, Energy = 1.375000
  min(rho) = 1.250000e-01, min(p) = 1.000000e-01, TV(rho) = 8.851488e-01

Snapshot at t = 0.200000
  Shock position (est.): x = 0.852500
  Mass = 0.562500, Momentum = 0.180000, Energy = 1.375000
  min(rho) = 1.250000e-01, min(p) = 1.000000e-01, TV(rho) = 8.816819e-01

Snapshot at t = 0.300000
  Shock position (est.): x = 0.780000
  Mass = 0.558780, Momentum = 0.263674, Energy = 1.358815
  min(rho) = 2.594841e-01, min(p) = 2.934428e-01, TV(rho) = 7.456568e-01

Saved figures: density_times.png, velocity_times.png, pressure_times.png
Total pressure repairs performed: 0



#### Output block3


Snapshot at t = 0.100000
  Shock pos (max|drho|): x = 0.675000, right-most +jump: x = 0.675000
  Mass = 0.562500, Momentum = 0.090000, Energy = 1.375000
  min(rho) = 1.250000e-01, min(p) = 1.000000e-01, TV(rho) = 8.851488e-01

Snapshot at t = 0.200000
  Shock pos (max|drho|): x = 0.852500, right-most +jump: x = 0.852500
  Mass = 0.562500, Momentum = 0.180000, Energy = 1.375000
  min(rho) = 1.250000e-01, min(p) = 1.000000e-01, TV(rho) = 8.816819e-01

Snapshot at t = 0.300000
  Shock pos (max|drho|): x = 0.780000, right-most +jump: x = 0.780000
  Mass = 0.558796, Momentum = 0.263669, Energy = 1.358867
  min(rho) = 2.604887e-01, min(p) = 2.951108e-01, TV(rho) = 7.446552e-01

Saved figures: density_times.png, velocity_times.png, pressure_times.png
Total pressure repairs performed: 0
Total density repairs performed: 0




### review_decision
accept

### review_comments
	Current Stage [A/3]
1. MAJOR

Category: runtime
Issue: Shape mismatch / broadcasting error when forming RHS: rhs_out[:,1:-1] expects (3,Ntot-2) but the code uses flux_interface[:,2:] - flux_interface[:,1:-1] which produces (3,Ntot-3). This causes ValueError: could not broadcast input array from shape (3,401) into shape (3,402).

Recommendation: Fix the interface indexing for the divergence. Replace the RHS line with the consistent indexing: rhs_out[:, 1:-1] = - (flux_interface[:, 1:] - flux_interface[:, :-1]) / dx. In general keep the convention that flux_interface[j] is F_{j+1/2}, so cell i uses flux_interface[i] and flux_interface[i-1]. Add an assertion or shape-check (e.g., assert flux_interface.shape[1] == U.shape[1]-1) to catch similar mismatches early.


2. MAJOR

Category: runtime
Issue: Possible zero time-step and/or infinite-loop when dt is reduced to exactly zero to land on an output time (dt = next_output_time - t may be zero), which may lead to no progress and spinning or saving logic confusion.

Recommendation: Guard against dt==0 before performing the SSPRK stages. If dt <= 0 (within machine tol), either (a) advance t directly to next_output_time and skip the RK update (compute snapshot using current U), or (b) set dt to a small positive minimum (e.g., dt = max(min_dt, computed_dt)) and proceed. Also restructure the output-saving condition to check 'if t >= t_out - tol and not yet saved' instead of relying on exact equality.


3. MAJOR

Category: structure
Issue: Ghost-cell boundary implementation uses reversed slices (e.g. U[:, 0:ng] = U[:, ng:2*ng][:, ::-1]) which mirrors interior values rather than performing transmissive/outflow (zero-gradient) boundary conditions. This is inconsistent with the stated transmissive BC intent and may introduce spurious reflections.

Recommendation: For transmissive (zero-gradient) BCs simply copy the nearest interior cells into ghosts: U[:, :ng] = U[:, ng:ng+ng] and U[:, -ng:] = U[:, -ng- ng:-ng] (or U[:, :ng] = U[:, ng:2*ng].copy() and U[:, -ng:] = U[:, -2*ng:-ng].copy()). If you actually want mirrored (reflecting) BCs, document that and apply an appropriate sign flip for momentum. Make the BC type configurable and add a small unit test that checks BC behavior.


4. MINOR

Category: accuracy
Issue: Ad-hoc pressure-floor fix and energy correction breaks strict conservation and may hide robustness problems. The code resets negative pressures by changing total energy in-place through a loop which can be expensive and non-conservative.

Recommendation: Use a conservative positivity-preserving fix or limiter: (a) enforce positivity on reconstructed states before the Riemann solve (you already floor rho and p there), (b) if negative pressure appears after update, try reducing the time-step (subcycling) or applying a local conservative fix that rescales internal energy while preserving total energy adjustments carefully, (c) vectorize the energy repair (avoid Python loops) if you must keep it, and print a diagnostic count of fixes so you can assess algorithmic issues rather than masking them silently.


5. MINOR

Category: structure
Issue: Output saving/loop condition logic is brittle and uses a complex boolean that can be confusing: while next_output_index < len(Tfinal_all) and abs(t - Tfinal_all[next_output_index]) < 1e-10 or ( ... ). This can lead to missed or duplicated saves and makes the control flow hard to follow.

Recommendation: Replace the complex boolean with a clearer pattern: after each successful advance, loop over remaining target times and if t >= t_out - tol and that t_out not yet saved, record the snapshot. Use a small tolerance (e.g., tol=1e-12) and maintain a boolean array/dict of saved flags. This is simpler, less error-prone, and easier to reason about.



	Current Stage [B/3]
1. MAJOR

Category: accuracy
Issue: Off-by-one in computing the discrete divergence (RHS) — rhs is filled for indices 1:-1 instead of the physical interior range (ng:-ng).

Recommendation: This is a correctness bug that breaks conservation and shifts which cells are updated. rhs_out should be filled for the actual physical cells (i = ng .. Ntot-ng-1). Concretely replace the line
  rhs_out[:, 1:-1] = - (flux_interface[:, 1:] - flux_interface[:, :-1]) / dx
with an expression that uses the ng offset, e.g.:
  rhs_out[:, ng:-ng] = - (flux_interface[:, ng:-ng+1] - flux_interface[:, ng-1:-ng]) / dx
(or equivalent indexing derived from your definition that flux_interface[:,j] = F_{j+1/2}). After this fix re-run and verify that global mass/momentum/energy are conserved to numerical tolerance and that shock/contact positions move monotonically with time.


2. MAJOR

Category: structure
Issue: Boundary-condition (ghost) handling and RK updates: ghost cells are overwritten only at the start of compute_RHS, but the SSPRK stages update the whole U array (including ghosts).

Recommendation: Either (a) enforce transmissive BCs after each RK substage (i.e. call apply_transmissive_BC right after computing U1 and U2, and before compute_RHS for the next stage), or (b) restrict updates to interior cells only (update U[:, ng:-ng] <- U[:, ng:-ng] + dt*rhs[:, ng:-ng]). Updating ghosts and then overwriting them only at the next RHS call can introduce transient nonphysical states that affect flux computations and conservation. Implementing one of these will improve robustness and conservation.


3. MINOR

Category: structure
Issue: apply_transmissive_BC copies a slice of interior cells into ghosts (U[:, :ng] = U[:, ng:2*ng]) rather than replicating the nearest interior cell into all ghost slots.

Recommendation: For a true zero-gradient (transmissive) BC it's more appropriate to copy the nearest interior cell value to all ghost cells, e.g.
  for i in range(ng):
      U[:, i] = U[:, ng]
      U[:, -1-i] = U[:, -ng-1]
or vectorized: U[:, :ng] = U[:, ng:ng+1] and U[:, -ng:] = U[:, -ng-1:-ng]. This avoids introducing spurious structure in ghost bands and is simpler/clearer.


4. MINOR

Category: accuracy
Issue: Shock-position detection via argmax(|drho|) can pick the contact or another strong gradient instead of the shock and lead to non-monotonic reported shock positions.

Recommendation: Improve the shock locator: search for the largest density jump on the expected right-most jump (e.g. choose the largest |drho| beyond a threshold x_min or in the region where density is larger than some value), or find the location of the largest positive density jump (shock usually produces a positive jump). Alternatively compute the shock speed from Rankine–Hugoniot posthoc. At minimum report both max |drho| index and the index of the right-most large jump so users can cross-check.


5. MINOR

Category: accuracy
Issue: Conservative repairs and numerical floors: the pressure/density floors and repairs are very small and applied globally after updates; some repairs are performed non-conservatively (e.g. setting rho floor without adjusting momentum/energy).

Recommendation: When enforcing rho floor, also rescale momentum and energy consistently (e.g. set mom = u_old * rho_new or scale momentum proportionally), or prefer using a positivity-preserving limiter at reconstruction/flux stage. For pressure floor repairs you already adjust total energy conservatively — apply the same careful conservative approach when fixing density. Also consider using a slightly larger floor (e.g. 1e-8 depending on machine precision and problem scales) to avoid frequent tiny repairs, and print a summary if repairs occur.



	Current Stage [C/3]
1. MINOR

Category: structure
Issue: apply_transmissive_BC uses slicing that relies on NumPy broadcasting to fill multiple ghost cells (U[:, -ng:] = U[:, -ng-1:-ng]).

Recommendation: Make the ghost filling explicit and robust for any ng (e.g. U[:, :ng] = np.repeat(U[:, ng:ng+1], ng, axis=1); U[:, -ng:] = np.repeat(U[:, -ng-1:-ng], ng, axis=1)). This is clearer and avoids surprises if ng is changed later.


2. MINOR

Category: accuracy
Issue: Shock locator (report_snapshot) can pick the wrong discontinuity (contact/rarefaction or spurious jump) and returns non-monotonic positions across snapshots.

Recommendation: Use a more robust shock-finding criterion (e.g. locate global maximum of |drho| restricted to a band around expected shock trajectory, or detect strong compression via negative velocity gradient / positive pressure jump). Also print the index and local density jump to help diagnosing which feature was detected.


3. MINOR

Category: accuracy
Issue: Small global conservation drift visible in runtime diagnostics (mass and energy change at later time).

Recommendation: This is likely from boundary treatment or numerical dissipation. To reduce drift: ensure fluxes are computed consistently at physical boundaries (you may want to compute and report boundary fluxes), double-check HLLC star-state energy formula, and consider increasing resolution or reducing visc/limiter strength for a conservation study.


4. MINOR

Category: accuracy
Issue: HLLC implementation contains compact algebra for star-region energy that is easy to get wrong; it should be double-checked numerically/analytically.

Recommendation: Add unit tests comparing the HLLC flux against known Riemann problems (including exact Sod at a single interface) and/or compare to a simpler Rusanov/HLL flux to validate energy behavior. If any discrepancies appear, rederive or adopt a well-tested HLLC implementation and include comments documenting the algebraic form used.





### iteration_count
3

### rev_image_description
	Current Stage [A/3]


	Current Stage [B/3]


	Current Stage [C/3]




