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
Plot the density, velocity, and pressure at t=0.4 in one figure.
Just save figs do not use plt.show() in the code.

[HINTS]:
Print necessary solving information to facilitate the reliability check of the solution.
Do not use 'np.trapz' as it is removed in NumPy 2.0, implement the integration manually.


### expanded_prob
We are asked to solve the 1-D compressible Euler equations (Sod shock tube) on x in [0,1] with an initial discontinuity at x=0.5. The system is hyperbolic and conservative: density rho, momentum rho*u, and total energy rho*E are the conserved variables. The solution for this Riemann problem contains a left-moving rarefaction, a contact discontinuity, and a right-moving shock. The core numerical goal is to compute the solution up to t = 0.4 and plot density, velocity, and pressure in a single figure (saving the figure to disk without plt.show()).

Primary mathematical and numerical challenges:
- Hyperbolicity and wave structure: correct capturing of three characteristic waves (two acoustic waves and a contact) demands a shock-capturing conservative scheme that respects the Rankine–Hugoniot conditions.
- Discontinuities: avoid nonphysical oscillations (Gibbs phenomena) near shocks and contacts while maintaining sharp profiles — requires nonlinear limiters or Riemann-based dissipation.
- Positivity: schemes must preserve positive density and pressure; some linear high-order reconstructions can produce negative internal energy.
- Stability: explicit time-stepping requires satisfying a CFL restriction based on maximum wave speed s_max = max(|u| + c) where c = sqrt(gamma*p/rho).
- Accuracy vs dissipation: first-order methods are robust but overly diffusive; higher-order reconstructions reduce dissipation but require limiters and careful time integration (e.g., SSP RK) to control spurious oscillations.
- Conservation: numerical method must be conservative to correctly capture shock speeds and global conserved quantities.
- Practical coding constraints: avoid naming conflicts (do not use variable name nu for viscosity; use visc), pass all physical / numerical constants explicitly to routines to avoid NameError, and avoid np.trapz (implement manual integration when needed). Also print run diagnostics (CFL, dt, iteration count, min/max of rho/p) to facilitate result verification.

Boundary conditions and implementation choices (recommended): use simple transmissive (zero-gradient) boundaries or ghost cells filled by constant extrapolation for this Riemann problem. Use a conservative finite-volume discretization with cell-averaged conserved variables. Compute fluxes at interfaces either using an approximate Riemann solver (Rusanov/HLL/HLLC/Roe) or a central scheme augmented with controlled artificial viscosity ('visc').

Outputs and verification:
- Produce and save a single figure with three subplots (density, velocity, pressure) at t=0.4.
- Print solver diagnostics at regular intervals: current time, dt, CFL, max wave speed, min/max of rho, p, total conserved integrals (computed by manual quadrature), and iteration count to aid reproducibility and reliability checks.

### solution_plans
	Current Stage [A/2]
solu_name='MUSCL-Hancock finite-volume with HLLC Riemann solver and slope limiter' content='Governing idea:\n- Use a conservative finite-volume discretization with second-order spatial accuracy via MUSCL (piecewise linear) reconstruction and a Hancock predictor step for time-centering; compute numerical interface fluxes with an HLLC approximate Riemann solver to capture contact/entropy features sharply. Advance in time with a second-order SSP Runge–Kutta or the MUSCL-Hancock two-step method. Enforce CFL using s_max = max(|u| + c). Pass all constants explicitly (gamma, CFL, t_end, N, visc if used for optional small artificial viscosity) to every function.\n\nAlgorithmic steps (numbered):\n1) Initialization:\n   - Set domain [0,1], number of cells N, cell centers x_i, cell width dx = 1/N.\n   - Initialize primitive variables (rho, u, p) according to the Riemann initial data; convert to conserved variables U = [rho, rho*u, rho*E] with E = 0.5*u^2 + p/((gamma-1)*rho).\n   - Define boundary treatment: ghost cells on left/right with transmissive (copy nearest interior cell) extrapolation.\n2) Time loop setup:\n   - Provide constants explicitly: gamma, CFL, t_end, max_iter.\n   - At each step compute primitive variables from conserved U, local sound speed c = sqrt(gamma*p/rho), and maximum signal speed s_max = max(|u| + c) over domain.\n   - Compute dt = CFL * dx / s_max (and ensure last step lands exactly on t_end if desired).\n3) MUSCL reconstruction (for each conserved variable or for primitive variables—prefer primitives to preserve positivity):\n   - Compute limited slopes sigma_i using a slope limiter (e.g., minmod, van Leer). Ensure limiter function is passed gamma/CFL etc as needed.\n   - Reconstruct left/right states at each interface: q_L = q_i + 0.5*sigma_i, q_R = q_{i+1} - 0.5*sigma_{i+1}.\n4) Hancock predictor (half-step):\n   - Evolve reconstructed states by half-step: q_L* = q_L - (dt/2/dx)*(F(q_{i+1}) - F(q_i)) and similarly for q_R* using fluxes computed from reconstructed states. (Use primitive -> conserved conversions as necessary.)\n5) Riemann solve at interfaces:\n   - For each interface, convert q_L*, q_R* to conserved form and compute flux using the HLLC approximate Riemann solver (requires estimates of wave speeds; use Einfeldt or Davis estimates). HLLC preserves contact and resolves intermediate states with low dissipation.\n6) Full-step update (conservative):\n   - U_i^{n+1} = U_i^{n} - (dt/dx)*(F_{i+1/2} - F_{i-1/2}). Optionally apply a small global or adaptive visc artificial diffusion term if positivity issues arise; if used, call it visc and pass visc as an explicit parameter.\n7) Boundary update and diagnostics:\n   - Update ghost cells, compute diagnostics: current time, dt, CFL, s_max, min(rho), min(p), total mass/momentum/energy via manual quadrature (sum(U*r_i*dx) etc), iteration count; print these.\n8) Continue until t >= t_end. After completion, compute primitive variables and save a figure with three subplots (rho, u, p) at t=0.4 using plt.savefig(filename) but do not call plt.show().\n\nStability / accuracy / complexity / efficiency limitations:\n- Stability: Explicit time-stepping requires dt <= CFL*dx/s_max; strong shocks reduce dt via large s_max. HLLC is stable for standard CFL (<=0.5 typically for second-order), but using too-aggressive limiters or no limiter may cause oscillations.\n- Accuracy: Second-order spatial plus second-order temporal accuracy; contact is well-resolved by HLLC but small smearing remains due to limiter choice. Using primitive-variable reconstruction is safer for positivity but must be limited properly near discontinuities.\n- Complexity: HLLC involves more algebra and per-interface calculations (wave-speed estimates) than simple Rusanov; cost per step ~ O(N) with a larger constant. Memory cost modest with ghost cells and temporary arrays for reconstructed states.\n- Efficiency: Reasonable for 1D; vectorize loops where possible. For production or multi-D, consider Riemann-solver optimizations. If speed is critical, a lower-cost approximate solver (Rusanov) reduces CPU but increases diffusion.\n- Additional caution: Ensure all functions take gamma, CFL, visc, and boundary choices as explicit arguments to avoid NameError; handle small negative pressures by floor limiting to a tiny positive value rather than letting the scheme crash.'

	Current Stage [B/2]
solu_name='Central finite-volume scheme with Jameson-style artificial visc (visc) and SSP-RK time integration' content='Governing idea:\n- Use a robust central finite-volume scheme (staggered or cell-centered central flux) that avoids the complexity of Riemann solvers by adding controlled, scale-selective artificial viscosity terms (second- and fourth-difference terms à la Jameson) to provide dissipation near shocks while preserving accuracy in smooth regions. Name the viscosity coefficient visc (explicit parameter). Use strong-stability-preserving Runge–Kutta (SSP RK3) for time advancement. Pass gamma, CFL, visc, N, and other constants explicitly to all functions.\n\nAlgorithmic steps (numbered):\n1) Initialization (same as Plan 1):\n   - Domain, N, dx, initial primitive and conserved variables, ghost cells for transmissive BCs.\n2) Compute central fluxes:\n   - At each interface compute a central (Rusanov-free) flux estimate: F_c = 0.5*(F(U_L) + F(U_R)). Do primitive reconstructions with a limiter (or use simple piecewise-constant for first-order central scheme).\n3) Add artificial visc terms (Jameson-style):\n   - Define second-difference dissipation D2_i = visc2 * (U_{i+1} - 2*U_i + U_{i-1}), and optionally a fourth-difference term D4_i = -visc4 * (U_{i+2} - 4*U_{i+1} + 6*U_i - 4*U_{i-1} + U_{i-2}). Use local sensors to turn on D4 selectively (e.g., based on pressure or density gradients). Use single parameter names visc2 and visc4 but ensure they are passed explicitly; shorter implementations may set visc4 = visc2*alpha.\n   - Combine interface residual: R_i = -(F_c_{i+1/2} - F_c_{i-1/2})/dx + (D2_i + D4_i)/dx (sign conventions consistent with conservative update). The visc terms damp high-frequency oscillations while minimally affecting smooth regions when properly tuned.\n4) Time integration with SSP-RK3:\n   - Use three sub-steps: U^(1) = U^n + dt*R(U^n), U^(2) = 3/4 U^n + 1/4(U^(1) + dt*R(U^(1))), U^(n+1) = 1/3 U^n + 2/3(U^(2) + dt*R(U^(2))). Ensure R recomputed at each stage.\n5) CFL and dt:\n   - Compute primitive variables, sound speeds, s_max, then dt = CFL*dx/s_max. Pass gamma and CFL explicitly to the dt routine.\n6) Positivity and limiters:\n   - After each stage, check rho and p, and if pressure < p_min or rho < rho_min, apply a small positivity fix (e.g., reset to floor or reduce dt). Print warnings if corrective action taken.\n7) Diagnostics and outputs:\n   - At intervals print time, dt, CFL, s_max, min/max rho and p, and global conserved integrals computed by manual summation (do not use np.trapz). Save final figure with rho, u, p at t=0.4 via plt.savefig and avoid plt.show().\n\nStability / accuracy / complexity / efficiency limitations:\n- Stability: The artificial visc approach is robust for strong shocks, but stability depends on proper tuning of visc parameters; overly small visc yields oscillations, overly large visc smears shocks. The scheme remains explicit and must satisfy CFL <= value determined empirically (often <=0.5 for SSP-RK3 with visc terms).\n- Accuracy: Central + visc approach can be made formally high-order in smooth regions (with high-order D4 tuned off there) but is more diffusive around contacts compared with Riemann-solver-based methods. Contacts are typically less sharply resolved than with HLLC/Roe.\n- Complexity: Implementation is simpler (no Riemann solver) but requires careful implementation of visc terms and sensor logic; proper ghost-cell handling for higher-order differences requires extra ghost cells (two on each side for D4).\n- Efficiency: Per-step arithmetic can be lower-than-HLLC if visc terms and central fluxes are vectorized; however, smaller dt may be required if visc terms cause constraints. Tuning visc reduces wasted CPU but requires experimentation.\n- Additional cautions: Use explicit parameter visc (and visc2/visc4 if separate) everywhere; compute all integrals by explicit summation (sum(conservative*dx) or more accurate composite quadrature coded manually). Print detailed run diagnostics to enable tuning of visc and CFL.'



### technical_spec
	Current Stage [A/3]
Second-order finite-volume MUSCL scheme with SSP RK2 time integration and HLLC Riemann solver. Data flow:
- Initialize conserved variables (U) with ghost cells.
- Time-stepping loop: apply transmissive BCs; compute primitive variables; determine dt from CFL and max(|u|+c); compute RHS via MUSCL reconstruction (minmod limiter) and HLLC fluxes; advance with SSP RK2. Diagnostics printed periodically include dt, max signal speed, min/max density and pressure, and integrated conserved quantities (mass, momentum, energy) computed by manual summation. Final primitives are plotted and saved to disk. All physical/numerical constants are passed explicitly to functions to avoid NameError.

	Current Stage [B/3]
Second-order finite-volume MUSCL solver with SSP RK2 and HLLC Riemann flux.

Architecture overview:
- Initialization: build initial primitive states with ghost cells and convert to conserved variables.
- Preallocation: allocate large reusable arrays (rhs, Fpad) once before time loop to avoid repeated memory allocations.
- Time integration: explicit SSP RK2 (Heun-like TVD). Each RK stage computes RHS via MUSCL reconstruction and HLLC fluxes.
- Reconstruction: piecewise-linear MUSCL with minmod limiter computed via explicit slicing (no np.roll) to avoid wrap-around across domain boundaries.
- Riemann solver: HLLC with early-case handling (SL>=0, SR<=0) and guarded computation of SM and star states; denominators clamped to avoid division by (near) zero.
- Boundary conditions: transmissive (zero-gradient) via explicit np.repeat copies into ghost cells.
- Diagnostics: print dt, max signal speed, min/max density and pressure, and integrated conserved quantities (mass/momentum/energy) via manual summation. Also initial/ending conservation checks.
- Output: saves labeled figure of density, velocity, and pressure at final time (no plt.show()).

Data flow:
U (conserved, shape (3,n_total)) -> cons_to_prim -> compute slopes -> reconstruct interface primitives -> hllc_flux -> assemble Fpad -> compute rhs -> RK update -> positivity fix -> loop.

All physical/numerical constants are passed explicitly to functions (gamma, dx, CFL, ng, etc.).

	Current Stage [C/3]
This module implements a 1-D Euler solver (Sod shock tube) using a second-order MUSCL finite-volume scheme, SSP RK2 time integration, and HLLC Riemann flux. Data flow:
- Initialize primitive states on N physical cells with ng ghost cells, convert to conserved variables U.
- Preallocate arrays (rhs, Fext and F_if) once in solver() before the time loop to avoid repeated allocations.
- Each RK stage: apply transmissive BCs, cons->prim, compute slopes (minmod) with explicit slicing, reconstruct interface primitives, compute HLLC interface fluxes (returns F_if), place into a padded flux array Fext, compute rhs = -(Fext[:,i+1]-Fext[:,i])/dx for each cell.
- HLLC includes standard SL/SR estimates, Rankine-Hugoniot SM, p_star consistency, guarded denominators, and diagnostics/assertions for shapes and non-finite values.
- Boundary conditions: transmissive via apply_bc which is used both at initialization and during time-stepping to ensure consistent ghost filling.
- Diagnostics: prints dt, smax, min/max primitive values, and integrated conserved quantities (mass/momentum/energy) computed manually (sum*dx). Conservation checks printed at the end.
- Output: saves a labeled figure with density, velocity, and pressure at t_final (no plt.show()).

Key fixes vs previous run:
- Fixed broadcasting error by making interface dimension explicit and using a padded flux array Fext of shape (3, n_total+1); assertions check F_if shape == n_total-1.
- Centralized ghost-cell handling: initial_sod now calls apply_bc to fill ghosts.
- Added defensive assertions and informative diagnostics when shapes or nonfinite values occur.
- Preallocated arrays outside the time loop (rhs, Fext, F_if) and never reallocated inside the time loop.

All physical/numerical constants are passed explicitly to functions. The optional 'visc' argument is kept but documented as unused.



### python_codes

#### Script block1:

```python
#!/usr/bin/env python3
"""
Sod shock tube (1D Euler) solver using MUSCL reconstruction + HLLC Riemann solver
with a second-order SSP RK2 time integrator.

- Uses primitive-variable slope limiting (minmod) for MUSCL reconstruction.
- HLLC Riemann solver resolves contact and shock well.
- Transmissive boundary conditions via ghost cells.
- Prints solver diagnostics and saves a figure (density, velocity, pressure) at t=0.4.

All physical / numerical constants are passed explicitly to functions.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

# ---------------------- Utility conversions and fluxes ----------------------

def prim_to_cons(rho, u, p, gamma):
    # conservative variables: [rho, rho*u, rho*E]
    E = 0.5 * u * u + p / ((gamma - 1.0) * rho)
    return np.stack((rho, rho * u, rho * E), axis=0)


def cons_to_prim(U, gamma):
    # U shape (3, n)
    rho = U[0].copy()
    mom = U[1].copy()
    E_tot = U[2].copy()
    # protect against tiny/negative density
    rho = np.maximum(rho, 1e-12)
    u = mom / rho
    e_internal = E_tot - 0.5 * u * u
    p = np.maximum((gamma - 1.0) * rho * e_internal, 1e-12)
    return rho, u, p


def flux_from_cons(U, gamma):
    # compute flux vector F(U) for Euler: [rho*u, rho*u^2 + p, (rho*E + p)*u]
    rho = U[0]
    mom = U[1]
    E_tot = U[2]
    u = mom / rho
    # internal energy and pressure
    e_internal = E_tot / rho - 0.5 * u * u
    p = (gamma - 1.0) * rho * e_internal
    F0 = mom
    F1 = mom * u + p
    F2 = (E_tot + p) * u
    return np.stack((F0, F1, F2), axis=0)

# ------------------------- Slope limiter (minmod) --------------------------

def minmod_limiter(a, b):
    # minmod for arrays a and b
    m = np.zeros_like(a)
    mask = (a * b) > 0.0
    m[mask] = np.where(np.abs(a[mask]) < np.abs(b[mask]), a[mask], b[mask])
    return m

# ------------------------------ HLLC solver --------------------------------

def hllc_flux( rhoL, uL, pL, rhoR, uR, pR, gamma ):
    # Inputs are arrays (or scalars) of left and right primitive states
    # Outputs flux array shape (3, n_interfaces)
    # Left conserved
    EL = pL / ((gamma - 1.0) * rhoL) + 0.5 * uL * uL
    ER = pR / ((gamma - 1.0) * rhoR) + 0.5 * uR * uR
    UL = np.stack((rhoL, rhoL * uL, rhoL * EL), axis=0)
    UR = np.stack((rhoR, rhoR * uR, rhoR * ER), axis=0)
    FL = flux_from_cons(UL, gamma)
    FR = flux_from_cons(UR, gamma)

    # sound speeds
    cL = np.sqrt(np.maximum(gamma * pL / rhoL, 0.0))
    cR = np.sqrt(np.maximum(gamma * pR / rhoR, 0.0))

    # Estimate wave speeds (Davis)
    SL = np.minimum(uL - cL, uR - cR)
    SR = np.maximum(uL + cL, uR + cR)

    # Prevent SL >= SR
    # If SL >= 0 -> return FL; if SR <= 0 -> return FR; else HLLC
    # Compute SM (middle wave speed)
    denom = (rhoL * (SL - uL) - rhoR * (SR - uR))
    # avoid division by zero
    small = 1e-12
    denom = np.where(np.abs(denom) < small, np.sign(denom) * small + small, denom)
    SM = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / denom

    # Compute star states
    # UL* = rhoL*(SL - uL)/(SL - SM) * [1, SM, E_L + (SM - uL)*(SM + pL/(rhoL*(SL - uL)))]
    coefL = rhoL * (SL - uL) / (SL - SM)
    coefR = rhoR * (SR - uR) / (SR - SM)

    E_starL = EL + (SM - uL) * (SM + pL / (rhoL * (SL - uL)))
    E_starR = ER + (SM - uR) * (SM + pR / (rhoR * (SR - uR)))

    UL_star = np.stack((coefL,
                        coefL * SM,
                        coefL * E_starL), axis=0)
    UR_star = np.stack((coefR,
                        coefR * SM,
                        coefR * E_starR), axis=0)

    # Now assemble flux
    # initialize with zeros
    F = np.zeros_like(FL)

    # Cases
    # SL >= 0 -> FL
    maskL = SL >= 0.0
    if np.any(maskL):
        F[:, maskL] = FL[:, maskL]

    # SR <= 0 -> FR
    maskR = SR <= 0.0
    if np.any(maskR):
        F[:, maskR] = FR[:, maskR]

    # SL < 0 < SR -> HLL or HLLC
    mask_middle = (~maskL) & (~maskR)
    if np.any(mask_middle):
        # subcases: SM >= 0 -> left star region else right star
        mask_SM_pos = (SM >= 0.0) & mask_middle
        mask_SM_neg = (SM < 0.0) & mask_middle
        if np.any(mask_SM_pos):
            idx = mask_SM_pos
            F[:, idx] = FL[:, idx] + SL[idx] * (UL_star[:, idx] - UL[:, idx])
        if np.any(mask_SM_neg):
            idx = mask_SM_neg
            F[:, idx] = FR[:, idx] + SR[idx] * (UR_star[:, idx] - UR[:, idx])

    return F

# ------------------------ RHS computation (finite-volume) -------------------

def compute_rhs(U, gamma, dx, ng):
    # U shape (3, n_total); returns RHS dU/dt same shape, ghost cells zeroed
    # Apply primitive recon, MUSCL (minmod), HLLC fluxes
    n_total = U.shape[1]

    # Convert to primitive
    rho, u, p = cons_to_prim(U, gamma)

    # Prepare slopes (minmod) for each primitive variable
    # Use simple two-point minmod between forward and backward differences
    dq_fwd = np.roll(rho, -1) - rho
    dq_bwd = rho - np.roll(rho, 1)
    slope_rho = minmod_limiter(dq_bwd, dq_fwd)

    dq_fwd = np.roll(u, -1) - u
    dq_bwd = u - np.roll(u, 1)
    slope_u = minmod_limiter(dq_bwd, dq_fwd)

    dq_fwd = np.roll(p, -1) - p
    dq_bwd = p - np.roll(p, 1)
    slope_p = minmod_limiter(dq_bwd, dq_fwd)

    # Reconstruct interface states: left state at i+1/2 is from cell i: qR = q_i + 0.5*slope_i
    # right state at i+1/2 is from cell i+1: qL = q_{i+1} - 0.5*slope_{i+1}
    # number of interfaces = n_total - 1
    n_if = n_total - 1
    rho_L = rho[:-1] + 0.5 * slope_rho[:-1]
    rho_R = rho[1:] - 0.5 * slope_rho[1:]
    u_L = u[:-1] + 0.5 * slope_u[:-1]
    u_R = u[1:] - 0.5 * slope_u[1:]
    p_L = p[:-1] + 0.5 * slope_p[:-1]
    p_R = p[1:] - 0.5 * slope_p[1:]

    # Ensure positive density/pressure for interface states
    rho_L = np.maximum(rho_L, 1e-12)
    rho_R = np.maximum(rho_R, 1e-12)
    p_L = np.maximum(p_L, 1e-12)
    p_R = np.maximum(p_R, 1e-12)

    # Compute HLLC flux at each interface
    F_if = hllc_flux(rho_L, u_L, p_L, rho_R, u_R, p_R, gamma)

    # Compute residual for interior cells only (physical cells indices ng..n_total-ng-1)
    rhs = np.zeros_like(U)
    # For cells i = 1..n_total-2 (interfaces from 0..n_total-2), residual for cell i is -(F_{i+1/2} - F_{i-1/2})/dx
    # We'll compute for all cells but ghost rows stay zero later
    # F_if index k corresponds to interface between cells k and k+1
    # So for cell i, left interface is i-1, right is i
    # Compute flux difference array: shape (3, n_total) with zeros at boundaries
    Fpad = np.zeros((3, n_total))
    Fpad[:, 1:-1] = F_if[:, 1:-1]  # assign interior interfaces; idx 1..n_total-2
    # Now residual
    rhs[:, :] = - (Fpad[:, 1:] - Fpad[:, :-1]) / dx

    # Zero-out ghost cell rhs
    rhs[:, :ng] = 0.0
    rhs[:, -ng:] = 0.0

    return rhs

# ------------------------------ Boundary BCs --------------------------------

def apply_bc(U, ng):
    # transmissive (zero-gradient) BC: copy edge interior values into ghost cells
    # U shape (3, n_total)
    # left ghosts
    U[:, :ng] = U[:, ng:ng+1]
    # right ghosts
    U[:, -ng:] = U[:, -ng-1:-ng]

# ----------------------------- Initial condition ----------------------------

def initial_sod(N, gamma, ng):
    # N = number of physical cells
    n_total = N + 2 * ng
    dx = 1.0 / N
    # cell centers for physical cells
    x = (np.arange(N) + 0.5) * dx
    # allocate arrays including ghosts
    rho = np.zeros(n_total)
    u = np.zeros(n_total)
    p = np.zeros(n_total)

    # fill physical cells
    # left state for x <= 0.5
    left_mask = x <= 0.5
    right_mask = x > 0.5
    rho_ph = np.where(left_mask, 1.0, 0.125)
    u_ph = np.zeros_like(rho_ph)
    p_ph = np.where(left_mask, 1.0, 0.1)

    # place into arrays with ghost offset
    rho[ng:ng+N] = rho_ph
    u[ng:ng+N] = u_ph
    p[ng:ng+N] = p_ph

    # fill ghost cells by simple copy of nearest physical cell
    rho[:ng] = rho[ng]
    rho[-ng:] = rho[ng+N-1]
    u[:ng] = u[ng]
    u[-ng:] = u[ng+N-1]
    p[:ng] = p[ng]
    p[-ng:] = p[ng+N-1]

    U = prim_to_cons(rho, u, p, gamma)
    return U, x, dx

# ---------------------------- Time stepping loop ----------------------------

def solver(N=400, gamma=1.4, CFL=0.5, t_end=0.4, visc=0.0, max_steps=20000, ng=2):
    # Set up
    U, x_phys, dx = initial_sod(N, gamma, ng)
    n_total = U.shape[1]

    t = 0.0
    step = 0

    # initial primitive to compute dt
    rho, u, p = cons_to_prim(U, gamma)
    c = np.sqrt(np.maximum(gamma * p / rho, 0.0))
    smax = np.max(np.abs(u) + c)
    dt = CFL * dx / (smax + 1e-12)

    print(f"Starting solver: N={N}, dx={dx:.5e}, CFL={CFL}, t_end={t_end}, gamma={gamma}")
    print(f"Initial dt={dt:.5e}, smax={smax:.5e}")

    # For diagnostics frequency
    diag_interval = max(1, int(0.02 / (CFL * dx / (smax + 1e-12))))

    while t < t_end and step < max_steps:
        # Enforce boundary conditions
        apply_bc(U, ng)

        # compute primitive & dt based on current state
        rho, u, p = cons_to_prim(U, gamma)
        c = np.sqrt(np.maximum(gamma * p / rho, 0.0))
        smax = np.max(np.abs(u) + c)
        dt = CFL * dx / (smax + 1e-12)
        # adjust dt to not overshoot t_end
        if t + dt > t_end:
            dt = t_end - t

        # Stage 1 (k1)
        rhs1 = compute_rhs(U, gamma, dx, ng)
        U1 = U + dt * rhs1

        # Apply BCs for intermediate state
        apply_bc(U1, ng)

        # Stage 2 (k2)
        rhs2 = compute_rhs(U1, gamma, dx, ng)
        U_new = 0.5 * (U + U1 + dt * rhs2)  # SSPRK2 (Heun-like TVD)

        # Positivity fix: ensure rho>0 and p>0 by adjusting energy if necessary
        rho_new, u_new, p_new = cons_to_prim(U_new, gamma)
        # If any negative, apply floors and correct energy conservatively
        rho_floor = 1e-12
        p_floor = 1e-12
        bad_rho = rho_new < rho_floor
        if np.any(bad_rho):
            rho_new[bad_rho] = rho_floor
            U_new[0] = rho_new
            U_new[1] = rho_new * u_new
            U_new[2] = rho_new * (0.5 * u_new * u_new + p_new / ((gamma - 1.0) * rho_new))

        bad_p = p_new < p_floor
        if np.any(bad_p):
            p_new[bad_p] = p_floor
            U_new[2] = rho_new * (0.5 * u_new * u_new + p_new / ((gamma - 1.0) * rho_new))

        # Advance
        U = U_new
        t += dt
        step += 1

        # Diagnostics printout
        if (step % max(1, diag_interval) == 0) or (t >= t_end) or (step < 10):
            # compute global conserved integrals (manual quadrature: sum over cells * dx)
            rho_phys = U[0, ng:-ng]
            mom_phys = U[1, ng:-ng]
            E_phys = U[2, ng:-ng]
            mass = np.sum(rho_phys) * dx
            momentum = np.sum(mom_phys) * dx
            energy = np.sum(E_phys) * dx
            rho_min = np.min(rho_phys)
            rho_max = np.max(rho_phys)
            # compute primitive for pressure min
            _, u_phys, p_phys = cons_to_prim(U[:, ng:-ng], gamma)
            p_min = np.min(p_phys)
            p_max = np.max(p_phys)

            print(f"step={step:6d}, t={t:.5f}, dt={dt:.3e}, smax={smax:.3e}, mass={mass:.6f}, E={energy:.6f}")
            print(f"   rho_min={rho_min:.6e}, rho_max={rho_max:.6e}, p_min={p_min:.6e}, p_max={p_max:.6e}")

    if step >= max_steps:
        print("Warning: reached maximum number of steps before t_end")

    # Final physical solution
    rho_final, u_final, p_final = cons_to_prim(U[:, ng:-ng], gamma)
    return x_phys, rho_final, u_final, p_final, step, t, dt

# ------------------------ Main execution and plotting -----------------------

def main():
    # Parameters (passed explicitly)
    gamma = 1.4
    CFL = 0.6
    t_end = 0.4
    N = 400  # number of physical cells
    visc = 0.0  # optional artificial viscosity parameter (unused in HLLC variant)

    x, rho, u, p, steps, t_final, dt_final = solver(N=N, gamma=gamma, CFL=CFL, t_end=t_end, visc=visc)

    print(f"Finished: steps={steps}, t_final={t_final:.6f}, last_dt={dt_final:.3e}")
    print(f"Final rho: min={rho.min():.6e}, max={rho.max():.6e}")
    print(f"Final p:   min={p.min():.6e}, max={p.max():.6e}")

    # Plot and save figure with density, velocity, pressure
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    axes[0].plot(x, rho, 'b-', lw=1.5)
    axes[0].set_ylabel('Density')
    axes[0].set_title(f'Sod shock tube at t={t_final:.3f} (N={N})')
    axes[0].grid(True)

    axes[1].plot(x, u, 'r-', lw=1.5)
    axes[1].set_ylabel('Velocity')
    axes[1].grid(True)

    axes[2].plot(x, p, 'g-', lw=1.5)
    axes[2].set_ylabel('Pressure')
    axes[2].set_xlabel('x')
    axes[2].grid(True)

    plt.tight_layout()
    fname = 'sod_musch_hllc.png'
    plt.savefig(fname, dpi=150)
    # do not call plt.show() as requested
    print(f"Saved figure to {fname}")


if __name__ == '__main__':
    main()
```



#### Script block2:

```python
#!/usr/bin/env python3
"""
Revised Sod shock tube solver (1D Euler) using MUSCL + HLLC with SSP RK2.
Improvements and fixes per review:
- Fixed conserved->primitive internal energy bug.
- Fixed interface flux padding indexing to avoid shape mismatches.
- Replaced np.roll-based slopes with explicit slicing (no wrap-around).
- Preallocated large arrays once outside time loop (rhs, Fpad) to avoid reallocations.
- Hardened HLLC with early-case handling, denominator clamping and basic diagnostics.
- apply_bc uses explicit np.repeat for ghost filling.
- Diagnostic prints of initial states and conservation checks at end.

Notes:
- Uses `visc` parameter but currently unused (kept for API compatibility).
- Saves figure (no plt.show()).
"""

import numpy as np
import matplotlib.pyplot as plt

# ---------------------- Utility conversions and fluxes ----------------------

def prim_to_cons(rho, u, p, gamma):
    # conservative variables: [rho, rho*u, rho*E]
    E = 0.5 * u * u + p / ((gamma - 1.0) * rho)
    return np.stack((rho, rho * u, rho * E), axis=0)


def cons_to_prim(U, gamma, rho_floor=1e-12, p_floor=1e-12):
    # U shape (3, n)
    rho = U[0].copy()
    mom = U[1].copy()
    E_tot = U[2].copy()
    # protect against tiny/negative density
    rho = np.maximum(rho, rho_floor)
    u = mom / rho
    # Correct internal energy per unit mass
    e_internal = E_tot / rho - 0.5 * u * u
    p = np.maximum((gamma - 1.0) * rho * e_internal, p_floor)
    return rho, u, p


def flux_from_cons(U, gamma):
    # compute flux vector F(U) for Euler: [rho*u, rho*u^2 + p, (rho*E + p)*u]
    rho = U[0]
    mom = U[1]
    E_tot = U[2]
    u = mom / rho
    e_internal = E_tot / rho - 0.5 * u * u
    p = (gamma - 1.0) * rho * e_internal
    F0 = mom
    F1 = mom * u + p
    F2 = (E_tot + p) * u
    return np.stack((F0, F1, F2), axis=0)


# ------------------------- Slope limiter (minmod) --------------------------

def minmod_limiter(a, b):
    # elementwise minmod for 1D arrays a and b
    m = np.zeros_like(a)
    mask = (a * b) > 0.0
    am = np.abs(a[mask])
    bm = np.abs(b[mask])
    m[mask] = np.where(am < bm, a[mask], b[mask])
    return m


# ------------------------------ HLLC solver --------------------------------

def hllc_flux(rhoL, uL, pL, rhoR, uR, pR, gamma):
    # Inputs arrays for left/right primitive states at interfaces
    # Outputs flux array shape (3, n_interfaces)
    small = 1e-12
    # left/right total energy per unit mass (E = e_internal + 0.5 u^2)
    EL = pL / ((gamma - 1.0) * rhoL) + 0.5 * uL * uL
    ER = pR / ((gamma - 1.0) * rhoR) + 0.5 * uR * uR
    UL = np.stack((rhoL, rhoL * uL, rhoL * EL), axis=0)
    UR = np.stack((rhoR, rhoR * uR, rhoR * ER), axis=0)
    FL = flux_from_cons(UL, gamma)
    FR = flux_from_cons(UR, gamma)

    # sound speeds
    cL = np.sqrt(np.maximum(gamma * pL / rhoL, 0.0))
    cR = np.sqrt(np.maximum(gamma * pR / rhoR, 0.0))

    # Estimate wave speeds (Davis)
    SL = np.minimum(uL - cL, uR - cR)
    SR = np.maximum(uL + cL, uR + cR)

    n_if = rhoL.size
    F = np.zeros((3, n_if))

    # Case SL >= 0 -> purely left flux
    maskL = SL >= 0.0
    if np.any(maskL):
        F[:, maskL] = FL[:, maskL]

    # Case SR <= 0 -> purely right flux
    maskR = SR <= 0.0
    if np.any(maskR):
        F[:, maskR] = FR[:, maskR]

    # Middle region: SL < 0 < SR
    mask_middle = (~maskL) & (~maskR)
    if np.any(mask_middle):
        idx = mask_middle
        # Compute SM only for middle
        denom = (rhoL * (SL - uL) - rhoR * (SR - uR))
        # robust clamping: where denom is too small, set to sign*small
        denom = np.where(np.abs(denom) < small, small * np.sign(denom + small), denom)
        SM = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / denom

        # Coefficients for star states, guard denominators SL-SM and SR-SM
        SL_minus_SM = SL - SM
        SR_minus_SM = SR - SM
        SL_minus_SM = np.where(np.abs(SL_minus_SM) < small, small * np.sign(SL_minus_SM + small), SL_minus_SM)
        SR_minus_SM = np.where(np.abs(SR_minus_SM) < small, small * np.sign(SR_minus_SM + small), SR_minus_SM)

        coefL = rhoL * (SL - uL) / SL_minus_SM
        coefR = rhoR * (SR - uR) / SR_minus_SM

        # Star-region energies per unit mass
        E_starL = EL + (SM - uL) * (SM + pL / (rhoL * (SL - uL) + small))
        E_starR = ER + (SM - uR) * (SM + pR / (rhoR * (SR - uR) + small))

        UL_star = np.stack((coefL, coefL * SM, coefL * E_starL), axis=0)
        UR_star = np.stack((coefR, coefR * SM, coefR * E_starR), axis=0)

        # Subcase where SM >= 0 -> left star; else right star
        mask_SM_pos = (SM >= 0.0) & mask_middle
        if np.any(mask_SM_pos):
            F[:, mask_SM_pos] = FL[:, mask_SM_pos] + SL[mask_SM_pos] * (
                UL_star[:, mask_SM_pos] - UL[:, mask_SM_pos]
            )
        mask_SM_neg = (SM < 0.0) & mask_middle
        if np.any(mask_SM_neg):
            F[:, mask_SM_neg] = FR[:, mask_SM_neg] + SR[mask_SM_neg] * (
                UR_star[:, mask_SM_neg] - UR[:, mask_SM_neg]
            )

        # Basic diagnostic check for NaNs/Infs in computed fluxes
        if not np.all(np.isfinite(F[:, idx])):
            # report a compact diagnostic
            nan_idx = np.where(~np.isfinite(F[:, idx].sum(axis=0)))[0]
            print("HLLC warning: nonfinite fluxes in middle region at local indices:", nan_idx)

    return F


# ------------------------ RHS computation (finite-volume) -------------------

def compute_rhs(U, gamma, dx, ng, rhs, Fpad):
    # U shape (3, n_total); rhs and Fpad are preallocated arrays to fill
    n_total = U.shape[1]

    # Convert to primitive
    rho, u, p = cons_to_prim(U, gamma)

    # Prepare slopes (minmod) for each primitive variable using explicit slices
    slope_rho = np.zeros_like(rho)
    slope_u = np.zeros_like(u)
    slope_p = np.zeros_like(p)

    # interior indices 1..n_total-2
    if n_total > 2:
        dq_plus = rho[2:] - rho[1:-1]
        dq_minus = rho[1:-1] - rho[:-2]
        slope_rho[1:-1] = minmod_limiter(dq_minus, dq_plus)

        dq_plus = u[2:] - u[1:-1]
        dq_minus = u[1:-1] - u[:-2]
        slope_u[1:-1] = minmod_limiter(dq_minus, dq_plus)

        dq_plus = p[2:] - p[1:-1]
        dq_minus = p[1:-1] - p[:-2]
        slope_p[1:-1] = minmod_limiter(dq_minus, dq_plus)

    # Reconstruct interface states between cells k and k+1 for k=0..n_total-2
    rho_L = rho[:-1] + 0.5 * slope_rho[:-1]
    rho_R = rho[1:] - 0.5 * slope_rho[1:]
    u_L = u[:-1] + 0.5 * slope_u[:-1]
    u_R = u[1:] - 0.5 * slope_u[1:]
    p_L = p[:-1] + 0.5 * slope_p[:-1]
    p_R = p[1:] - 0.5 * slope_p[1:]

    # Ensure positive density/pressure for interface states
    rho_L = np.maximum(rho_L, 1e-12)
    rho_R = np.maximum(rho_R, 1e-12)
    p_L = np.maximum(p_L, 1e-12)
    p_R = np.maximum(p_R, 1e-12)

    # Compute HLLC flux at each interface
    F_if = hllc_flux(rho_L, u_L, p_L, rho_R, u_R, p_R, gamma)

    # Place interface fluxes into preallocated Fpad: Fpad has shape (3, n_total)
    # Interface k corresponds to Fpad[:, k] for k=0..n_total-2
    Fpad.fill(0.0)
    Fpad[:, :-1] = F_if

    # Compute residuals: rhs = - (F_{i+1/2} - F_{i-1/2}) / dx
    # That is equivalent to differences of adjacent columns of Fpad
    # Fpad[:, j] is flux at interface between cells j and j+1
    rhs[:, :] = - (Fpad[:, 1:] - Fpad[:, :-1]) / dx

    # Zero-out ghost cell rhs
    rhs[:, :ng] = 0.0
    rhs[:, -ng:] = 0.0

    return rhs


# ------------------------------ Boundary BCs --------------------------------

def apply_bc(U, ng):
    # transmissive (zero-gradient) BC: copy edge interior values into ghost cells
    # use explicit repeat for clarity
    U[:, :ng] = np.repeat(U[:, ng:ng+1], ng, axis=1)
    U[:, -ng:] = np.repeat(U[:, -ng-1:-ng], ng, axis=1)


# ----------------------------- Initial condition ----------------------------

def initial_sod(N, gamma, ng):
    # N = number of physical cells
    n_total = N + 2 * ng
    dx = 1.0 / N
    # cell centers for physical cells
    x = (np.arange(N) + 0.5) * dx
    # allocate arrays including ghosts
    rho = np.zeros(n_total)
    u = np.zeros(n_total)
    p = np.zeros(n_total)

    # fill physical cells
    left_mask = x <= 0.5
    right_mask = x > 0.5
    rho_ph = np.where(left_mask, 1.0, 0.125)
    u_ph = np.zeros_like(rho_ph)
    p_ph = np.where(left_mask, 1.0, 0.1)

    # place into arrays with ghost offset
    rho[ng:ng+N] = rho_ph
    u[ng:ng+N] = u_ph
    p[ng:ng+N] = p_ph

    # fill ghost cells by simple copy of nearest physical cell
    rho[:ng] = rho[ng]
    rho[-ng:] = rho[ng+N-1]
    u[:ng] = u[ng]
    u[-ng:] = u[ng+N-1]
    p[:ng] = p[ng]
    p[-ng:] = p[ng+N-1]

    U = prim_to_cons(rho, u, p, gamma)
    return U, x, dx


# ---------------------------- Time stepping loop ----------------------------

def solver(N=400, gamma=1.4, CFL=0.5, t_end=0.4, visc=0.0, max_steps=20000, ng=2):
    # Set up
    U, x_phys, dx = initial_sod(N, gamma, ng)
    n_total = U.shape[1]

    # Preallocate arrays used every time step to avoid repeated allocation
    rhs = np.zeros_like(U)
    Fpad = np.zeros((3, n_total))

    # initial primitive to compute dt
    rho0, u0, p0 = cons_to_prim(U, gamma)
    c0 = np.sqrt(np.maximum(gamma * p0 / rho0, 0.0))
    smax = np.max(np.abs(u0) + c0)
    dt = CFL * dx / (smax + 1e-12)

    # diagnostics initial
    print(f"Starting solver: N={N}, dx={dx:.5e}, CFL={CFL}, t_end={t_end}, gamma={gamma}")
    print(f"Initial dt={dt:.5e}, smax={smax:.5e}")
    # verify initial left/right pressures exactly as prescribed
    print(f"Initial sample primitives (at few cells):")
    phys_slice = slice(ng, ng+min(5, N))
    print("  rho[0..] (first few phys):", rho0[phys_slice][:5])
    print("  p left sample:", p0[ng])
    print("  p right sample:", p0[-ng-1])

    # initial conserved integrals for mass conservation check
    mass0 = np.sum(U[0, ng:-ng]) * dx
    momentum0 = np.sum(U[1, ng:-ng]) * dx
    energy0 = np.sum(U[2, ng:-ng]) * dx

    # diagnostics interval based on dt
    diag_interval = max(1, int(0.02 / max(dt, 1e-12)))

    t = 0.0
    step = 0

    while t < t_end and step < max_steps:
        # Enforce boundary conditions
        apply_bc(U, ng)

        # compute primitive & dt based on current state
        rho, u, p = cons_to_prim(U, gamma)
        c = np.sqrt(np.maximum(gamma * p / rho, 0.0))
        smax = np.max(np.abs(u) + c)
        dt = CFL * dx / (smax + 1e-12)
        if t + dt > t_end:
            dt = t_end - t

        # Stage 1 (k1)
        rhs = compute_rhs(U, gamma, dx, ng, rhs, Fpad)
        U1 = U + dt * rhs

        # Apply BCs for intermediate state
        apply_bc(U1, ng)

        # Stage 2 (k2)
        rhs = compute_rhs(U1, gamma, dx, ng, rhs, Fpad)
        U_new = 0.5 * (U + U1 + dt * rhs)  # SSPRK2

        # Positivity fix
        rho_new, u_new, p_new = cons_to_prim(U_new, gamma)
        rho_floor = 1e-12
        p_floor = 1e-12
        bad_rho = rho_new < rho_floor
        if np.any(bad_rho):
            rho_new[bad_rho] = rho_floor
            U_new[0] = rho_new
            U_new[1] = rho_new * u_new
            U_new[2] = rho_new * (0.5 * u_new * u_new + p_new / ((gamma - 1.0) * rho_new))

        bad_p = p_new < p_floor
        if np.any(bad_p):
            p_new[bad_p] = p_floor
            U_new[2] = rho_new * (0.5 * u_new * u_new + p_new / ((gamma - 1.0) * rho_new))

        # Advance
        U = U_new
        t += dt
        step += 1

        # Diagnostics printout
        if (step % max(1, diag_interval) == 0) or (t >= t_end) or (step < 10):
            rho_phys = U[0, ng:-ng]
            mom_phys = U[1, ng:-ng]
            E_phys = U[2, ng:-ng]
            mass = np.sum(rho_phys) * dx
            momentum = np.sum(mom_phys) * dx
            energy = np.sum(E_phys) * dx
            rho_min = np.min(rho_phys)
            rho_max = np.max(rho_phys)
            _, u_phys, p_phys = cons_to_prim(U[:, ng:-ng], gamma)
            p_min = np.min(p_phys)
            p_max = np.max(p_phys)

            print(f"step={step:6d}, t={t:.5f}, dt={dt:.3e}, smax={smax:.3e}, mass={mass:.6f}, E={energy:.6f}")
            print(f"   rho_min={rho_min:.6e}, rho_max={rho_max:.6e}, p_min={p_min:.6e}, p_max={p_max:.6e}")

    if step >= max_steps:
        print("Warning: reached maximum number of steps before t_end")

    # Final physical solution
    rho_final, u_final, p_final = cons_to_prim(U[:, ng:-ng], gamma)

    # conservation checks
    mass_final = np.sum(rho_final) * dx
    momentum_final = np.sum((rho_final * u_final)) * dx
    energy_final = np.sum(U[2, ng:-ng]) * dx
    print("Conservation check:")
    print(f"  mass: initial={mass0:.8e}, final={mass_final:.8e}, diff={mass_final-mass0:.3e}")
    print(f"  momentum: initial={momentum0:.8e}, final={momentum_final:.8e}, diff={momentum_final-momentum0:.3e}")
    print(f"  energy: initial={energy0:.8e}, final={energy_final:.8e}, diff={energy_final-energy0:.3e}")

    # sanity checks
    if mass_final - mass0 > 1e-6:
        print("Warning: mass change exceeds tolerance")
    if np.min(p_final) < 1e-12 or np.min(rho_final) < 1e-12:
        print("Warning: small pressure/density encountered in final solution")

    return x_phys, rho_final, u_final, p_final, step, t, dt


# ------------------------ Main execution and plotting -----------------------

def main():
    # Parameters (passed explicitly)
    gamma = 1.4
    CFL = 0.6
    t_end = 0.4
    N = 400  # number of physical cells
    visc = 0.0  # kept for API compatibility (unused)

    x, rho, u, p, steps, t_final, dt_final = solver(N=N, gamma=gamma, CFL=CFL, t_end=t_end, visc=visc)

    print(f"Finished: steps={steps}, t_final={t_final:.6f}, last_dt={dt_final:.3e}")
    print(f"Final rho: min={rho.min():.6e}, max={rho.max():.6e}")
    print(f"Final p:   min={p.min():.6e}, max={p.max():.6e}")

    # Plot and save figure with density, velocity, pressure
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    axes[0].plot(x, rho, 'b-', lw=1.5)
    axes[0].set_ylabel('Density')
    axes[0].set_title(f'Sod shock tube at t={t_final:.3f} (N={N})')
    axes[0].grid(True)

    axes[1].plot(x, u, 'r-', lw=1.5)
    axes[1].set_ylabel('Velocity')
    axes[1].grid(True)

    axes[2].plot(x, p, 'g-', lw=1.5)
    axes[2].set_ylabel('Pressure')
    axes[2].set_xlabel('x')
    axes[2].grid(True)

    plt.tight_layout()
    fname = 'sod_musch_hllc_revised.png'
    plt.savefig(fname, dpi=150)
    print(f"Saved figure to {fname}")


if __name__ == '__main__':
    main()
```



#### Script block3:

```python
#!/usr/bin/env python3
"""
Sod shock tube (1D Euler) using MUSCL (minmod), SSP RK2 and HLLC.
Notes:
- All constants are passed explicitly to functions.
- Ghosts are filled by apply_bc centrally (initialization and each step).
- Preallocate rhs, Fext, F_if once in solver to avoid allocations inside the time loop.
- The 'visc' parameter is currently unused (kept for API compatibility).
"""

import numpy as np
import matplotlib.pyplot as plt

# ---------------------- Primitive/Conservative conversions -----------------

def prim_to_cons(rho, u, p, gamma):
    E = 0.5 * u * u + p / ((gamma - 1.0) * rho)
    return np.stack((rho, rho * u, rho * E), axis=0)


def cons_to_prim(U, gamma, rho_floor=1e-12, p_floor=1e-12):
    rho = U[0].copy()
    mom = U[1].copy()
    E_tot = U[2].copy()
    rho = np.maximum(rho, rho_floor)
    u = mom / rho
    e_internal = E_tot / rho - 0.5 * u * u
    p = np.maximum((gamma - 1.0) * rho * e_internal, p_floor)
    return rho, u, p


def flux_from_cons(U, gamma):
    rho = U[0]
    mom = U[1]
    E_tot = U[2]
    u = mom / rho
    e_internal = E_tot / rho - 0.5 * u * u
    p = (gamma - 1.0) * rho * e_internal
    F0 = mom
    F1 = mom * u + p
    F2 = (E_tot + p) * u
    return np.stack((F0, F1, F2), axis=0)

# ------------------------- Slope limiter (minmod) --------------------------

def minmod_limiter(a, b):
    m = np.zeros_like(a)
    mask = (a * b) > 0.0
    am = np.abs(a[mask])
    bm = np.abs(b[mask])
    m[mask] = np.where(am < bm, a[mask], b[mask])
    return m

# ------------------------------ HLLC solver --------------------------------

def hllc_flux(rhoL, uL, pL, rhoR, uR, pR, gamma):
    """Compute HLLC flux for arrays of interface left/right primitive states.
    Returns F_if shape (3, n_if).
    """
    small = 1e-12
    # total energy per unit mass
    EL = pL / ((gamma - 1.0) * rhoL) + 0.5 * uL * uL
    ER = pR / ((gamma - 1.0) * rhoR) + 0.5 * uR * uR
    UL = np.stack((rhoL, rhoL * uL, rhoL * EL), axis=0)
    UR = np.stack((rhoR, rhoR * uR, rhoR * ER), axis=0)
    FL = flux_from_cons(UL, gamma)
    FR = flux_from_cons(UR, gamma)

    # sound speeds
    cL = np.sqrt(np.maximum(gamma * pL / rhoL, 0.0))
    cR = np.sqrt(np.maximum(gamma * pR / rhoR, 0.0))

    # wave speed estimates (Davis)
    SL = np.minimum(uL - cL, uR - cR)
    SR = np.maximum(uL + cL, uR + cR)

    # number of interfaces
    n_if = rhoL.size
    F = np.zeros((3, n_if))

    # quick masks
    maskL = SL >= 0.0
    maskR = SR <= 0.0
    mask_middle = (~maskL) & (~maskR)

    # left-only and right-only
    if np.any(maskL):
        F[:, maskL] = FL[:, maskL]
    if np.any(maskR):
        F[:, maskR] = FR[:, maskR]

    if np.any(mask_middle):
        idx = mask_middle
        # robust SM via Rankine-Hugoniot
        denom = (rhoL * (SL - uL) - rhoR * (SR - uR))
        denom = np.where(np.abs(denom) < small, small * np.sign(denom + small), denom)
        SM = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / denom

        # star region pressures (consistent with momentum jump)
        p_starL = pL + rhoL * (SL - uL) * (SM - uL)
        p_starR = pR + rhoR * (SR - uR) * (SM - uR)
        # average p_star for robustness
        p_star = 0.5 * (p_starL + p_starR)

        # prevent nonphysical p_star
        p_star = np.maximum(p_star, small)

        # compute star densities
        SL_minus_SM = SL - SM
        SR_minus_SM = SR - SM
        SL_minus_SM = np.where(np.abs(SL_minus_SM) < small, small * np.sign(SL_minus_SM + small), SL_minus_SM)
        SR_minus_SM = np.where(np.abs(SR_minus_SM) < small, small * np.sign(SR_minus_SM + small), SR_minus_SM)

        rho_starL = rhoL * (SL - uL) / SL_minus_SM
        rho_starR = rhoR * (SR - uR) / SR_minus_SM

        # star-region conserved states (UL*, UR*)
        UstarL = np.zeros((3, idx.sum()))
        UstarR = np.zeros((3, idx.sum()))

        # energy in star region via total energy consistency
        # E_star = rho_star * (E + (SM - u) * (SM + p/(rho*(S - u))))
        E_starL = EL[idx] + (SM[idx] - uL[idx]) * (SM[idx] + pL[idx] / (rhoL[idx] * (SL[idx] - uL[idx]) + small))
        E_starR = ER[idx] + (SM[idx] - uR[idx]) * (SM[idx] + pR[idx] / (rhoR[idx] * (SR[idx] - uR[idx]) + small))

        UstarL[0, :] = rho_starL[idx]
        UstarL[1, :] = rho_starL[idx] * SM[idx]
        UstarL[2, :] = rho_starL[idx] * E_starL

        UstarR[0, :] = rho_starR[idx]
        UstarR[1, :] = rho_starR[idx] * SM[idx]
        UstarR[2, :] = rho_starR[idx] * E_starR

        # For middle region decide left or right star based on sign of SM
        mask_SM_pos = (SM >= 0.0) & mask_middle
        if np.any(mask_SM_pos):
            F[:, mask_SM_pos] = FL[:, mask_SM_pos] + SL[mask_SM_pos] * (
                UstarL[:, mask_SM_pos[idx]] - UL[:, mask_SM_pos]
            )
        mask_SM_neg = (SM < 0.0) & mask_middle
        if np.any(mask_SM_neg):
            F[:, mask_SM_neg] = FR[:, mask_SM_neg] + SR[mask_SM_neg] * (
                UstarR[:, mask_SM_neg[idx]] - UR[:, mask_SM_neg]
            )

        # check for nonfinite fluxes and provide diagnostics
        if not np.all(np.isfinite(F[:, idx])):
            bad = ~np.isfinite(F[:, idx]).all(axis=0)
            bad_indices = np.where(idx)[0][bad]
            print("HLLC fatal: nonfinite fluxes at interfaces:", bad_indices)
            for ib in bad_indices:
                iL = ib
                iR = ib + 1
                print(" Interface", ib)
                print("  rhoL,uL,pL =", rhoL[ib], uL[ib], pL[ib])
                print("  rhoR,uR,pR =", rhoR[ib], uR[ib], pR[ib])
            raise RuntimeError("Nonfinite HLLC flux encountered; see diagnostics above.")

    return F

# ------------------------------ Boundary BCs --------------------------------

def apply_bc(U, ng):
    # transmissive (zero-gradient) BC: copy nearest interior values into ghost cells
    U[:, :ng] = np.repeat(U[:, ng:ng+1], ng, axis=1)
    U[:, -ng:] = np.repeat(U[:, -ng-1:-ng], ng, axis=1)

# ----------------------------- Initial condition ----------------------------

def initial_sod(N, gamma, ng):
    # create physical cell centers
    n_total = N + 2 * ng
    dx = 1.0 / N
    x = (np.arange(N) + 0.5) * dx

    # allocate primitive arrays including ghosts
    rho = np.zeros(n_total)
    u = np.zeros(n_total)
    p = np.zeros(n_total)

    # fill physical cells
    left_mask = x <= 0.5
    rho_ph = np.where(left_mask, 1.0, 0.125)
    u_ph = np.zeros_like(rho_ph)
    p_ph = np.where(left_mask, 1.0, 0.1)

    rho[ng:ng+N] = rho_ph
    u[ng:ng+N] = u_ph
    p[ng:ng+N] = p_ph

    # convert to conserved and let apply_bc fill ghost cells consistently
    U = prim_to_cons(rho, u, p, gamma)
    apply_bc(U, ng)
    return U, x, dx

# ------------------------ RHS computation (finite-volume) -------------------

def compute_rhs(U, gamma, dx, ng, rhs, Fext, F_if):
    n_total = U.shape[1]
    n_if_expected = n_total - 1

    rho, u, p = cons_to_prim(U, gamma)

    # compute slopes with explicit slicing
    slope_rho = np.zeros_like(rho)
    slope_u = np.zeros_like(u)
    slope_p = np.zeros_like(p)

    if n_total > 2:
        dq_plus = rho[2:] - rho[1:-1]
        dq_minus = rho[1:-1] - rho[:-2]
        slope_rho[1:-1] = minmod_limiter(dq_minus, dq_plus)

        dq_plus = u[2:] - u[1:-1]
        dq_minus = u[1:-1] - u[:-2]
        slope_u[1:-1] = minmod_limiter(dq_minus, dq_plus)

        dq_plus = p[2:] - p[1:-1]
        dq_minus = p[1:-1] - p[:-2]
        slope_p[1:-1] = minmod_limiter(dq_minus, dq_plus)

    # reconstruct interface primitives for interfaces k=0..n_total-2
    rho_L = rho[:-1] + 0.5 * slope_rho[:-1]
    rho_R = rho[1:] - 0.5 * slope_rho[1:]
    u_L = u[:-1] + 0.5 * slope_u[:-1]
    u_R = u[1:] - 0.5 * slope_u[1:]
    p_L = p[:-1] + 0.5 * slope_p[:-1]
    p_R = p[1:] - 0.5 * slope_p[1:]

    # protect positivity
    rho_L = np.maximum(rho_L, 1e-12)
    rho_R = np.maximum(rho_R, 1e-12)
    p_L = np.maximum(p_L, 1e-12)
    p_R = np.maximum(p_R, 1e-12)

    # compute interface fluxes
    F_if_local = hllc_flux(rho_L, u_L, p_L, rho_R, u_R, p_R, gamma)

    # defensive checks
    assert F_if_local.ndim == 2 and F_if_local.shape[0] == 3, "HLLC returned unexpected shape"
    if F_if_local.shape[1] != n_if_expected:
        raise ValueError(f"HLLC returned {F_if_local.shape[1]} interfaces, expected {n_if_expected}")

    # copy into provided F_if buffer then into extended flux array Fext
    F_if[:, :] = F_if_local

    # Fext has shape (3, n_total+1). Place interface fluxes at indices 1..n_total-1
    Fext.fill(0.0)
    Fext[:, 1:-1] = F_if
    # pad endpoints for transmissive behaviour with nearest interface flux
    Fext[:, 0] = F_if[:, 0]
    Fext[:, -1] = F_if[:, -1]

    # compute rhs = - (F_{i+1/2} - F_{i-1/2}) / dx for each cell i
    rhs[:, :] = - (Fext[:, 1:] - Fext[:, :-1]) / dx

    # zero-out ghost-cell rhs
    rhs[:, :ng] = 0.0
    rhs[:, -ng:] = 0.0

    # check for nonfinite rhs
    if not np.all(np.isfinite(rhs)):
        bad = ~np.isfinite(rhs).all(axis=0)
        bad_idx = np.where(bad)[0]
        print("Nonfinite rhs at cell indices:", bad_idx)
        raise RuntimeError("Nonfinite RHS encountered; aborting.")

    return rhs

# ----------------------------- Time stepping loop ---------------------------

def solver(N=400, gamma=1.4, CFL=0.5, t_end=0.4, visc=0.0, max_steps=20000, ng=2):
    """Main solver. 'visc' parameter is currently unused (kept for API compatibility).
    Returns x (physical centers), rho, u, p (physical arrays), and final timing info.
    """
    U, x_phys, dx = initial_sod(N, gamma, ng)
    n_total = U.shape[1]
    n_if = n_total - 1

    # preallocate reusable arrays once
    rhs = np.zeros_like(U)
    Fext = np.zeros((3, n_total + 1))
    F_if = np.zeros((3, n_if))

    # initial dt based on primitive state
    rho0, u0, p0 = cons_to_prim(U, gamma)
    c0 = np.sqrt(np.maximum(gamma * p0 / rho0, 0.0))
    smax = np.max(np.abs(u0) + c0)
    dt = CFL * dx / (smax + 1e-12)

    print(f"Starting solver: N={N}, dx={dx:.5e}, CFL={CFL}, t_end={t_end}, gamma={gamma}")
    print(f"Initial dt={dt:.5e}, smax={smax:.5e}")
    phys_slice = slice(ng, ng+min(5, N))
    print("  rho[0..] (first few phys):", rho0[phys_slice][:5])
    print("  p left sample:", p0[ng], " p right sample:", p0[ng+N-1])

    mass0 = np.sum(U[0, ng:-ng]) * dx
    momentum0 = np.sum(U[1, ng:-ng]) * dx
    energy0 = np.sum(U[2, ng:-ng]) * dx

    diag_interval = max(1, int(0.02 / max(dt, 1e-12)))

    t = 0.0
    step = 0

    while t < t_end and step < max_steps:
        apply_bc(U, ng)

        rho, u, p = cons_to_prim(U, gamma)
        c = np.sqrt(np.maximum(gamma * p / rho, 0.0))
        smax = np.max(np.abs(u) + c)
        dt = CFL * dx / (smax + 1e-12)
        if t + dt > t_end:
            dt = t_end - t

        # Stage 1
        compute_rhs(U, gamma, dx, ng, rhs, Fext, F_if)
        U1 = U + dt * rhs
        apply_bc(U1, ng)

        # Stage 2
        compute_rhs(U1, gamma, dx, ng, rhs, Fext, F_if)
        U_new = 0.5 * (U + U1 + dt * rhs)

        # positivity fix (small floors)
        rho_new, u_new, p_new = cons_to_prim(U_new, gamma)
        rho_floor = 1e-12
        p_floor = 1e-12
        if np.any(rho_new < rho_floor):
            bad = rho_new < rho_floor
            rho_new[bad] = rho_floor
            U_new[0] = rho_new
            U_new[1] = rho_new * u_new
            U_new[2] = rho_new * (0.5 * u_new * u_new + p_new / ((gamma - 1.0) * rho_new))
        if np.any(p_new < p_floor):
            bad = p_new < p_floor
            p_new[bad] = p_floor
            U_new[2] = rho_new * (0.5 * u_new * u_new + p_new / ((gamma - 1.0) * rho_new))

        U = U_new
        t += dt
        step += 1

        if (step % max(1, diag_interval) == 0) or (t >= t_end) or (step < 10):
            rho_phys = U[0, ng:-ng]
            mom_phys = U[1, ng:-ng]
            E_phys = U[2, ng:-ng]
            mass = np.sum(rho_phys) * dx
            momentum = np.sum(mom_phys) * dx
            energy = np.sum(E_phys) * dx
            rho_min = np.min(rho_phys)
            rho_max = np.max(rho_phys)
            _, u_phys, p_phys = cons_to_prim(U[:, ng:-ng], gamma)
            p_min = np.min(p_phys)
            p_max = np.max(p_phys)

            print(f"step={step:6d}, t={t:.5f}, dt={dt:.3e}, smax={smax:.3e}, mass={mass:.6f}, E={energy:.6f}")
            print(f"   rho_min={rho_min:.6e}, rho_max={rho_max:.6e}, p_min={p_min:.6e}, p_max={p_max:.6e}")

    if step >= max_steps:
        print("Warning: reached maximum number of steps before t_end")

    rho_final, u_final, p_final = cons_to_prim(U[:, ng:-ng], gamma)

    mass_final = np.sum(rho_final) * dx
    momentum_final = np.sum((rho_final * u_final)) * dx
    energy_final = np.sum(U[2, ng:-ng]) * dx
    print("Conservation check:")
    print(f"  mass: initial={mass0:.8e}, final={mass_final:.8e}, diff={mass_final-mass0:.3e}")
    print(f"  momentum: initial={momentum0:.8e}, final={momentum_final:.8e}, diff={momentum_final-momentum0:.3e}")
    print(f"  energy: initial={energy0:.8e}, final={energy_final:.8e}, diff={energy_final-energy0:.3e}")

    if mass_final - mass0 > 1e-6:
        print("Warning: mass change exceeds tolerance")
    if np.min(p_final) < 1e-12 or np.min(rho_final) < 1e-12:
        print("Warning: small pressure/density encountered in final solution")

    return x_phys, rho_final, u_final, p_final, step, t, dt

# ------------------------ Main execution and plotting -----------------------

def main():
    gamma = 1.4
    CFL = 0.6
    t_end = 0.4
    N = 400
    visc = 0.0  # currently unused; kept for API compatibility

    x, rho, u, p, steps, t_final, dt_final = solver(N=N, gamma=gamma, CFL=CFL, t_end=t_end, visc=visc)

    print(f"Finished: steps={steps}, t_final={t_final:.6f}, last_dt={dt_final:.3e}")
    print(f"Final rho: min={rho.min():.6e}, max={rho.max():.6e}")
    print(f"Final p:   min={p.min():.6e}, max={p.max():.6e}")

    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    axes[0].plot(x, rho, 'b-', lw=1.5)
    axes[0].set_ylabel('Density')
    axes[0].set_title(f'Sod shock tube at t={t_final:.3f} (N={N})')
    axes[0].grid(True)

    axes[1].plot(x, u, 'r-', lw=1.5)
    axes[1].set_ylabel('Velocity')
    axes[1].grid(True)

    axes[2].plot(x, p, 'g-', lw=1.5)
    axes[2].set_ylabel('Pressure')
    axes[2].set_xlabel('x')
    axes[2].grid(True)

    plt.tight_layout()
    fname = 'sod_musch_hllc_revised.png'
    plt.savefig(fname, dpi=150)
    print(f"Saved figure to {fname}")


if __name__ == '__main__':
    main()
```




### runtime_outputs

#### Output block1

Starting solver: N=400, dx=2.50000e-03, CFL=0.6, t_end=0.4, gamma=1.4
Initial dt=1.26773e-03, smax=1.18322e+00

Runtime Error:
Traceback (most recent call last):
  File "C:\Users\admin\AppData\Local\Temp\ipykernel_50164\2663542217.py", line 138, in execute_code_tool
    exec(code, global_scope, global_scope)
  File "<string>", line 380, in <module>
  File "<string>", line 350, in main
  File "<string>", line 282, in solver
  File "<string>", line 189, in compute_rhs
ValueError: could not broadcast input array from shape (3,401) into shape (3,402)



#### Output block2

Starting solver: N=400, dx=2.50000e-03, CFL=0.6, t_end=0.4, gamma=1.4
Initial dt=1.26773e-03, smax=1.18322e+00
Initial sample primitives (at few cells):
  rho[0..] (first few phys): [1. 1. 1. 1. 1.]
  p left sample: 1.0
  p right sample: 0.1

Runtime Error:
Traceback (most recent call last):
  File "C:\Users\admin\AppData\Local\Temp\ipykernel_50164\2663542217.py", line 138, in execute_code_tool
    exec(code, global_scope, global_scope)
  File "<string>", line 421, in <module>
  File "<string>", line 392, in main
  File "<string>", line 310, in solver
  File "<string>", line 204, in compute_rhs
ValueError: could not broadcast input array from shape (3,403) into shape (3,404)



#### Output block3

<string>:17: RuntimeWarning: invalid value encountered in divide
Starting solver: N=400, dx=2.50000e-03, CFL=0.6, t_end=0.4, gamma=1.4
Initial dt=1.26773e-03, smax=1.18322e+00
  rho[0..] (first few phys): [1. 1. 1. 1. 1.]
  p left sample: 1.0  p right sample: 0.1
step=     1, t=0.00127, dt=1.268e-03, smax=1.183e+00, mass=0.562500, E=1.375000
   rho_min=1.250000e-01, rho_max=1.000000e+00, p_min=1.000000e-01, p_max=1.000000e+00
step=     2, t=0.00209, dt=8.242e-04, smax=1.820e+00, mass=0.562500, E=1.375000
   rho_min=1.250000e-01, rho_max=1.000000e+00, p_min=1.000000e-01, p_max=1.000000e+00
step=     3, t=0.00288, dt=7.861e-04, smax=1.908e+00, mass=0.562500, E=1.375000
   rho_min=1.250000e-01, rho_max=1.000000e+00, p_min=1.000000e-01, p_max=1.000000e+00
step=     4, t=0.00365, dt=7.735e-04, smax=1.939e+00, mass=0.562500, E=1.375000
   rho_min=1.250000e-01, rho_max=1.000000e+00, p_min=1.000000e-01, p_max=1.000000e+00
step=     5, t=0.00439, dt=7.415e-04, smax=2.023e+00, mass=0.562500, E=1

... [Truncated: Content too long (8659 chars). Hidden middle part.] ...

2.611324e-01, rho_max=9.999991e-01, p_min=2.960670e-01, p_max=9.999987e-01
step=   555, t=0.38022, dt=6.792e-04, smax=2.209e+00, mass=0.538890, E=1.271014
   rho_min=2.611339e-01, rho_max=9.999870e-01, p_min=2.960661e-01, p_max=9.999819e-01
step=   570, t=0.39041, dt=6.792e-04, smax=2.209e+00, mass=0.536367, E=1.259865
   rho_min=2.611336e-01, rho_max=9.998734e-01, p_min=2.960665e-01, p_max=9.998228e-01
step=   585, t=0.40000, dt=7.969e-05, smax=2.209e+00, mass=0.533996, E=1.249385
   rho_min=2.611320e-01, rho_max=9.992365e-01, p_min=2.960663e-01, p_max=9.989314e-01
Conservation check:
  mass: initial=5.62500000e-01, final=5.33995850e-01, diff=-2.850e-02
  momentum: initial=0.00000000e+00, final=3.10528471e-01, diff=3.105e-01
  energy: initial=1.37500000e+00, final=1.24938498e+00, diff=-1.256e-01
Finished: steps=585, t_final=0.400000, last_dt=7.969e-05
Final rho: min=2.611320e-01, max=9.992365e-01
Final p:   min=2.960663e-01, max=9.989314e-01
Saved figure to sod_musch_hllc_revised.png




### review_decision
revise

### review_comments
	Current Stage [A/3]
1. MAJOR

Category: runtime
Issue: Runtime ValueError when assembling flux array: attempted to broadcast array of shape (3,401) into shape (3,402) in compute_rhs during assignment to Fpad.

Recommendation: Fix the indexing when constructing the padded interface-flux array. A robust and simple approach is to place the interface fluxes into the first n_total-1 columns of Fpad and then compute the cell residuals from differences of adjacent columns. Replace the current lines

    Fpad = np.zeros((3, n_total))
    Fpad[:, 1:-1] = F_if[:, 1:-1]
    rhs[:, :] = - (Fpad[:, 1:] - Fpad[:, :-1]) / dx

with, for example:

    Fpad = np.zeros((3, n_total))
    Fpad[:, :-1] = F_if    # F_if has shape (3, n_total-1)
    rhs[:, :] = - (Fpad[:, 1:] - Fpad[:, :-1]) / dx

This aligns interface index k with Fpad[:, k] (interface between cells k and k+1) and avoids off-by-one size mismatches.


2. MAJOR

Category: accuracy
Issue: Incorrect conversion from conserved to primitive variables in cons_to_prim: internal energy is computed as E_tot - 0.5*u*u instead of (E_tot / rho) - 0.5*u*u. This yields erroneous pressures (and thus wrong characteristic speeds and fluxes).

Recommendation: Correct the internal energy computation to divide the total energy by density before subtracting kinetic energy per unit mass. Replace the e_internal line with:

    e_internal = E_tot / rho - 0.5 * u * u

and keep the pressure floor as you have it. Verify by printing a few primitive values after initialization and before the time loop, and check that initial left/right pressures are exactly the prescribed values (1.0 and 0.1).


3. MINOR

Category: structure
Issue: Slope computation uses np.roll, which implicitly creates periodic coupling (wrap-around) between the last and first cell. Although ghost cells are assigned via apply_bc, using np.roll can be fragile and confusing and may produce incorrect slopes at domain edges.

Recommendation: Compute slopes using explicit slicing to avoid wrap-around and to make the treatment of ghost cells explicit. For example, for interior i = 1..n_total-2:

    dq_plus = rho[2:] - rho[1:-1]
    dq_minus = rho[1:-1] - rho[:-2]
    slope_rho[1:-1] = minmod_limiter(dq_minus, dq_plus)

and set slope_rho[:1] = 0 and slope_rho[-1:] = 0 (or mirror). Do the same for u and p. This prevents accidental coupling across domain boundaries and is clearer to readers.


4. MINOR

Category: accuracy
Issue: Potential numerical robustness issues in HLLC solver: (a) SL/SR/SM estimation and denominator safeguards are present but nonstandard; (b) star-state energy formula and coefficient handling can be fragile when SL≈SM or SR≈SM; (c) the code uses small heuristics which can introduce bias.

Recommendation: After fixing cons_to_prim, add unit tests / diagnostic prints for SL, SR, SM ranges and check for NaNs/Infs. Consider using a well-tested HLLC implementation pattern (explicitly handle SL>=0 and SR<=0 early, and when SL<0<SR compute HLLC state with careful algebra). If you keep the current formulas, at least clamp denominators symmetrically (e.g. denom += small * np.sign(denom) or set denom[...] = small where |denom|<small) and document the choice. Also add assertions or np.nan checks to catch instabilities early.


5. MINOR

Category: structure
Issue: Various implementation and clarity issues that affect maintainability (unused visc parameter, diagnostic interval computation, and broadcasting reliance in apply_bc).

Recommendation: A few targeted improvements:
- Either remove the unused visc parameter or use it (and document it).
- Compute diag_interval in a clearer way (e.g. diag_interval = max(1, int(0.02 / dt))) after dt is known rather than mixing smax inside the expression).
- For apply_bc, it's fine to rely on numpy broadcasting, but consider using explicit repeats or np.tile when copying a single column into ng ghosts for readability, e.g. U[:, :ng] = np.repeat(U[:, ng:ng+1], ng, axis=1).
- Add simple automated checks at the end: mass conservation tolerance, min pressure/density not below floor, and optional comparison to a reference/analytic Sod solution for plausibility.



	Current Stage [B/3]
1. MAJOR

Category: runtime
Issue: Broadcasting error in compute_rhs: ValueError: could not broadcast input array from shape (3,403) into shape (3,404).

Recommendation: The root cause is a shape mismatch between the number of interface fluxes (n_if = n_total - 1) returned by hllc_flux and the slice of Fpad you assign them into. Make the interface dimension explicit and allocate/assign consistently. Options:
- Simpler fix: compute n_if = n_total - 1 and replace the current Fpad handling with an explicit interface-flux array F_if_padded of length n_total+1 (or length n_total with clear mapping). Example robust approach:
    n_if = n_total - 1
    # padded flux array with one extra entry at each end
    Fext = np.zeros((3, n_if + 1))    # shape (3, n_total)
    Fext[:, 1:-0] = F_if              # place F_if into the interior interfaces (choose consistent offsets)
    rhs[:, :] = - (Fext[:, 1:] - Fext[:, :-1]) / dx
  Or even clearer, avoid Fpad and compute rhs using F_if directly for interior cells:
    # For physical cells i = 0..n_total-1, interface fluxes are F_if[0..n_if-1]
    rhs[:, ng:-ng] = - (F_if[:, ng:ng+N] - F_if[:, ng-1:ng+N-1]) / dx
  Also add shape assertions before assignment to raise an informative error if shapes disagree (e.g. assert F_if.shape[1] == n_total - 1).


2. MINOR

Category: structure
Issue: Ghost-cell handling is inconsistent: initial_sod fills ghosts manually while apply_bc provides a reusable transmissive BC function. This risks subtle indexing inconsistencies.

Recommendation: Centralize ghost-cell filling. Call apply_bc(U, ng) immediately after creating U in initial_sod (or remove the manual ghost filling in initial_sod). This ensures a single correct implementation of boundary conditions and reduces the chance of off-by-one errors between initialization and time stepping.


3. MINOR

Category: accuracy
Issue: HLLC star-state energy (E_starL/E_starR) and SM computation use nonstandard/fragile algebra (denominator clamping, ad-hoc terms). This can produce inaccurate star energies or nonphysical states.

Recommendation: Verify and simplify the HLLC implementation against a canonical reference (e.g. Toro 'Riemann Solvers and Numerical Methods for Fluid Dynamics' or Batten et al.). Use the standard robust formulae:
- compute SL, SR (e.g. Davis or Einfeldt estimates), compute SM via the Rankine–Hugoniot relation,
- compute p_star from momentum jump, then construct UL*, UR* consistently.
Add unit tests comparing HLLC flux against e.g. a reference HLLC implementation on simple states (including strong shocks and near-vacuum) and check positivity (rho>0, p>0) of star states.


4. MINOR

Category: structure
Issue: Parameter 'visc' is passed around but unused and current code comments say it's unused. Leaving unused API arguments can confuse readers and maintainers.

Recommendation: Either implement the intended viscous/artificial viscosity behavior or remove the parameter from the public API. If kept for future extension, document clearly (docstring) that it is currently unused and default is 0.0.


5. MINOR

Category: runtime
Issue: Limited defensive programming and diagnostics for shape/NaN issues; current NaN diagnostic prints local indices but does not raise or show shapes, making debugging harder.

Recommendation: Add explicit shape assertions and richer diagnostics before assignments. For example:
    assert F_if.ndim == 2
    assert F_if.shape[0] == 3
    assert F_if.shape[1] == n_total - 1, f'expected {n_total-1} interfaces, got {F_if.shape[1]}'
When nonfinite values are detected, print the primitive states around the failing interface (rhoL/uL/pL/rhoR/uR/pR) and raise an exception so execution halts with a clear message during development. This avoids silent corruption and helps identify the exact location of the bug.



	Current Stage [C/3]
1. MAJOR

Category: runtime
Issue: Runtime warning 'invalid value encountered in divide' observed at startup; non-finite divisions may occur during cons->prim or HLLC computations.

Recommendation: Pinpoint the source of the divide-by-zero by enabling NumPy error raising (e.g. np.seterr(all='raise')) while debugging or print shapes and min/max of denominators. Enforce safe floors before any division (e.g. apply rho_floor to densities prior to dividing momentum by rho, and protect terms like rho*(S-u) before using them as denominators). Replace tiny 1e-12 floor with a somewhat larger, problem-appropriate floor (e.g. 1e-8) during debugging to reveal indexing issues. Once fixed, switch back to conservative floors and remove numpy warnings.


2. MAJOR

Category: accuracy
Issue: Significant change in global conserved quantities reported (mass dropped from 0.5625 to 0.533996, momentum moved from 0 to ~0.31, energy decreased by ~0.1256). This is large and unexpected for a closed shock-tube calculation and indicates either boundary-treatment issues or a bug in flux assembly/time integration.

Recommendation: Decide whether the domain is meant to be closed or open. If closed, use reflecting BCs (mirror ghost cells) or ensure zero net flux at outer boundaries. If transmissive/open is intended, verify that the loss of mass/momentum/energy equals the time integral of the computed boundary fluxes (i.e. check conservation by comparing change in integrals to boundary flux contributions). Instrument the code to print the boundary fluxes (F_if at leftmost and rightmost real interfaces) each step and integrate them in time to verify the source of the change. If discrepancies remain, inspect compute_rhs and the mapping between F_if and Fext — the current Fext indexing and endpoint assignment is brittle and easy to get off-by-one.


3. MAJOR

Category: structure
Issue: HLLC implementation uses complicated boolean masking and mixed boolean-indexing (e.g. mask_SM_pos[idx] and mask_SM_neg[idx]) that is hard to follow and likely to cause subtle broadcasting/indexing mistakes.

Recommendation: Refactor the HLLC code to use explicit index arrays from np.where and operate on compact local arrays for the 'middle' interfaces. For example: idx = np.where(mask_middle)[0]; then build local arrays using idx to compute SM, p_star, rho_star, UstarL, UstarR and assign F[:, idx] directly. Avoid double-masking expressions like mask_SM_pos[idx] inside indexing for mixed-size arrays; they are error-prone. Add unit tests for HLLC using hand-crafted left/right states (including sonic and nearly-vacuum cases) to exercise each branch (SL>0, SR<0, middle SM>0/SM<0) and check that fluxes satisfy Rankine-Hugoniot consistency and are finite.


4. MAJOR

Category: structure
Issue: Flux assembly into Fext is unintuitive and can easily introduce off-by-one mapping errors between interfaces and cells (Fext size = n_total+1, F_if size = n_total-1, endpoints are set by copying interior interface values). This complicates reasoning about flux differences and boundary fluxes.

Recommendation: Simplify the flux bookkeeping. Use a clear and documented mapping: let cells be indexed 0..n_total-1, interfaces 0..n_total (n_total+1 values) with Fext[j] ≡ flux at x_{j-1/2}; then populate Fext[1:-1] with F_if[0: n_total-1], and compute endpoints Fext[0] and Fext[-1] by calling the Riemann solver with (ghost,leftmost) and (rightmost,ghost) reconstructed states respectively (i.e. compute boundary interface fluxes explicitly from reconstructed ghost values). This makes the boundary fluxes explicit so you can validate the net conservation by summing boundary contributions. Add assertions and small tests to check that sum(rhs) * dx + boundary_flux_integral == 0 (within roundoff) for a simple advected constant test.


5. MINOR

Category: structure
Issue: Code robustness and clarity issues: (a) apply_bc uses np.repeat with awkward slices which is cryptic, (b) 'visc' is unused (commented as kept for API), (c) tiny floors (1e-12) may mask bugs or produce underflow; (d) diagnostic frequency computed from dt can be extreme.

Recommendation: a) Replace apply_bc with clearer code (e.g. U[:, :ng] = U[:, ng:ng+1].repeat(ng, axis=1) or a simple loop) and add comments. b) Either implement visc (if intended) or remove it from the API to avoid confusion. c) Use slightly larger floors during development (1e-8..1e-10), and expose floors as function arguments to tune. d) Make diag_interval a fixed small integer or clamp it to a sensible range, and consider printing aggregated diagnostics (mass/momentum/energy and boundary fluxes) less frequently. Also add unit-tests (mass advection, isentropic vortex if extended) and a small regression test comparing to an exact Sod solver to validate accuracy.





### iteration_count
3

### rev_image_description
	Current Stage [A/3]


	Current Stage [B/3]


	Current Stage [C/3]




