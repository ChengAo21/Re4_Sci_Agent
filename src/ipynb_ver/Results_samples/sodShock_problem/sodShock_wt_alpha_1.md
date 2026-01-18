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
We must solve the one-dimensional compressible Euler equations (conservation of mass, momentum, and energy) on x in [0,1] with a Riemann initial condition (left state ρ=1.0,u=0.0,p=1.0 for x<=0.5, right state ρ=0.125,u=0.0,p=0.1 for x>0.5) and γ=1.4. The goal is to produce stable and efficient time-evolving solutions and plot primitive fields (density, velocity, pressure) at t=0.1,0.2,0.3. Practical coding constraints require that all physical/numerical constants are explicitly passed to functions (or defined globally) to prevent NameError, the viscosity parameter should be named visc (never nu), and numerical integration in diagnostics must avoid np.trapz (implement trapezoidal integration manually). Additionally produce diagnostic outputs: shock position, checks for spurious oscillations, and global-conservation summaries.

Primary mathematical and numerical challenges:
- Hyperbolic system with discontinuous initial data: self-steepening waves, shock formation, contact discontinuity and rarefaction fan require shock-capturing schemes that maintain monotonicity without excessive smearing.
- Nonlinear Riemann physics: exact Riemann solution is available for benchmarking but expensive; approximate Riemann solvers (Roe, HLL, HLLC) have tradeoffs in robustness and accuracy (e.g., contact resolution vs. positivity).
- Conservation form and conservative update: spatial discretization must be conservative to preserve integral quantities; time integration must respect CFL stability condition.
- Avoiding spurious oscillations: high-order reconstructions need limiters (MUSCL, WENO) or monotonicity-preserving strategies; positivity-preserving corrections may be needed to avoid negative density/pressure.
- Stability/time-stepping: explicit schemes require timestep control via CFL based on maximum wave speed; high-order SSP RK is preferable for stability with non-linear limiters.
- Boundary conditions: choose appropriate open/transmissive boundaries consistent with characteristic directions to avoid reflections.
- Diagnostics and verification: robust diagnostics to locate shock/contact positions, measure total mass/momentum/energy, detect non-physical values, and quantify oscillations.
- Performance and complexity: tradeoff between accuracy and computational cost; higher-order methods (WENO5 + RK3) provide superior resolution of small-scale features but cost more per cell than second-order MUSCL + HLLC.
- Numerical integration for diagnostics: implement manual trapezoidal rule to compute integrals, not np.trapz.

Given these challenges we present two alternative, well-structured solution strategies that meet the coding constraints and target the required plots and diagnostics.

### solution_plans
	Current Stage [A/2]
solu_name='Finite-Volume MUSCL-Hancock with HLLC Riemann Solver and SSP-RK2' content='Governing idea:\nUse a second-order conservative finite-volume method with linear MUSCL reconstruction (slope limiter), the pointwise HLLC approximate Riemann solver for robust shock/contact capturing, and a Hancock predictor-corrector type update or a two-stage SSP-RK2 time integrator. Enforce positivity-preserving fixes and CFL timestep control. Compute diagnostics (shock position, conservative totals via manual trapezoidal rule, detection of negative/NaN values).\n\nAlgorithmic steps (numbered):\n1) Domain and data structures\n   - Discretize x ∈ [0,1] into N uniform cells; store cell-centered conserved vector U = [ρ, ρu, ρE] and derived primitives W = [ρ,u,p].\n   - Pass constants explicitly: gamma, CFL, final times list [0.1,0.2,0.3], visc (if any artificial viscosity term is used), N, boundary condition type.\n\n2) Initial conditions\n   - Initialize primitives from the Riemann initial data; convert to conserved variables U.\n\n3) Reconstruction (MUSCL)\n   - For each primitive variable (or better, for characteristic-projected primitive variables), compute limited slopes using a robust limiter (minmod or van Leer). Use explicit functions that accept constants and arrays.\n   - Reconstruct left/right states at each cell interface: W_L, W_R (second-order in space).\n\n4) HLLC Flux computation\n   - For each interface compute flux F_{i+1/2} = HLLC(W_L, W_R, gamma). Implementation notes: compute wave estimates (S_L, S_R) and contact speed S_*; ensure pressure/density positivity. All constants passed explicitly.\n\n5) Time integration (SSP-RK2 variant)\n   - Compute timestep dt = CFL * min(Δx / max(|u|+c)) where c = sqrt(gamma*p/ρ).\n   - Stage 1 (predictor): U^* = U^n + dt * RHS(U^n) with RHS from finite-volume flux differences: RHS_i = -(F_{i+1/2}-F_{i-1/2})/Δx.\n   - Reconstruct primitives from U^* and apply positivity fix if necessary.\n   - Stage 2 (corrector): U^{n+1} = 0.5*(U^n + U^* + dt * RHS(U^*)). This yields second-order temporal accuracy and TVD properties.\n   - Alternatively, implement the MUSCL-Hancock (predictor-corrector) with half-step evolution of reconstructed states using local flux Jacobians.\n\n6) Boundary conditions\n   - Use transmissive (zero-gradient/outflow) boundary conditions consistent with outward characteristics; implement ghost cells and fill them each timestep.\n\n7) Positivity and limiting\n   - After updates, check for ρ>ρ_min and p>p_min (small positive thresholds); if violated, apply a conservative floor or revert to first-order update locally (i.e., throttle slope to zero) to restore positivity.\n\n8) Diagnostics and shock detection\n   - At requested output times (0.1,0.2,0.3) compute primitives and produce arrays for plotting.\n   - Shock/contact position: compute the location of maximum |∂ρ/∂x| (approximate gradient via central differences) and report x_shock. Track its time evolution for sanity checks.\n   - Global conserved quantities: compute integrals of conserved fields using manual trapezoidal rule implementation (do not call np.trapz). Print mass, momentum, energy at each saved time and their relative change from initial values.\n   - Oscillation detection: compute local extrema count in density and/or a measure of total variation TV = Σ |ρ_{i+1}-ρ_i| and report if TV increases spuriously compared to fine-grid baseline.\n\n9) Plotting and output\n   - Produce three figures (density, velocity, pressure) with curves at t=0.1,0.2,0.3 overlaid or separate per requirement. Save files to disk; do not call plt.show().\n\nStability / accuracy / complexity / efficiency limitations:\n- Accuracy: nominally second-order in space and time. Contact resolution depends on the HLLC solver and limiter; some smearing remains compared with higher-order schemes.\n- Stability: CFL-limited explicit scheme; CFL must be <= ~0.5 (tune depending on limiter and RK variant). Positivity-preserving fixes may locally reduce accuracy.\n- Robustness: HLLC is robust for shocks and preserves contact waves better than HLL; Roe may be sharper but needs entropy fixes and positivity control.\n- Complexity: per time-step cost is O(N) with a modest constant (reconstruction + Riemann solves). Memory O(N).\n- Efficiency tradeoffs: characteristic-projection reconstruction improves accuracy but increases computation; using primitives for slope limiting is cheaper but may be less accurate across strong shocks.\n- Edge cases: for very coarse grids, oscillations or negative pressures may occur and require fallback to first-order updates locally.\n\nStop.'

	Current Stage [B/2]
solu_name='High-order WENO5 with Local Lax-Friedrichs Flux Splitting and SSP-RK3 (+ optional small visc artificial dissipation named visc)' content='Governing idea:\nUse a high-order (fifth-order) finite-difference/finite-volume hybrid WENO5 reconstruction to obtain interface states for flux evaluation, use local Lax-Friedrichs (Rusanov) flux splitting for robustness and simplicity (or HLLC if contact sharpness is critical), and integrate in time with a strong-stability-preserving third-order Runge-Kutta (SSP-RK3) for improved temporal accuracy. Optionally include a small explicit artificial visc (labelled visc in code) term for targeted dissipation if spurious oscillations are detected. Ensure positivity-preserving limiter and explicit passing of constants.\n\nAlgorithmic steps (numbered):\n1) Domain and variables\n   - Use uniform grid with N cells; store conserved U in cell centers. Pass gamma, CFL, visc (artificial viscosity coefficient, possibly zero), final times list, N, boundary condition type to all helper functions.\n\n2) Initialization\n   - Set primitive variables from Riemann initial condition and convert to conserved variables.\n\n3) WENO5 reconstruction\n   - For each component of flux-contributing quantities (commonly primitive variables or characteristic-projected conserved), reconstruct left and right states at each interface using WENO5 stencil and nonlinear weights. Implement smoothness indicators and epsilon to avoid division by zero; implement in explicit functions with constants passed.\n   - Optionally perform characteristic projection: compute local Jacobian eigenvectors to reconstruct in characteristic space to reduce spurious oscillations across shocks.\n\n4) Flux evaluation with LF splitting (or HLLC)\n   - Compute physical flux F(U) from reconstructed states. Use local Lax-Friedrichs flux: F_{i+1/2} = 0.5*(F(U_L)+F(U_R)) - 0.5*alpha*(U_R-U_L) where alpha = max(|u|+c) over the two states.\n   - Alternatively, compute HLLC flux if contact resolution must be prioritized.\n\n5) Optional explicit artificial visc term\n   - Add a low-order viscous flux term in conserved variables: F_visc_i = visc * (U_{i+1} - 2U_i + U_{i-1})/Δx; visc must be passed explicitly and tuned small to damp Gibbs oscillations. This term should be conservative and discrete consistent with finite-volume update.\n\n6) SSP-RK3 time stepping\n   - Compute dt via CFL based on maximum wave speed across domain. Use SSP-RK3 stages:\n     a) U^{(1)} = U^n + dt * RHS(U^n)\n     b) U^{(2)} = 3/4*U^n + 1/4*(U^{(1)} + dt*RHS(U^{(1)}))\n     c) U^{n+1} = 1/3*U^n + 2/3*(U^{(2)} + dt*RHS(U^{(2)}))\n   - After each stage enforce positivity-preserving correction or revert to local first-order update if needed.\n\n7) Boundary conditions\n   - Use ghost cells and populate via transmissive/outflow conditions or characteristic-based extrapolation.\n\n8) Diagnostics and shock localization\n   - At desired times, compute primitives and save for plotting.\n   - Locate sharp features by finding interface index where |ρ_{i+1}-ρ_i| is maximal; interpolate to estimate shock position. Print shock position and shock speed approximation (Δx/Δt of shock movement).\n   - Compute global integrals via manual trapezoidal rule for mass/momentum/energy and print relative changes.\n   - Detect spurious oscillations using local extrema counts or by comparing high-frequency content (e.g., discrete total variation) to a coarse baseline.\n\n9) Plotting and output\n   - Save three separate figure files for density, velocity, and pressure at t=0.1,0.2,0.3; do not call plt.show().\n\nStability / accuracy / complexity / efficiency limitations:\n- Accuracy: WENO5+RK3 is nominally 5th order in smooth regions and resolves contact/rarefaction sharply; shocks are captured with oscillation control but are still smeared over a few cells.\n- Stability: explicit method bound by CFL; WENO+RK3 typically requires CFL <= ~0.6 depending on implementation; artificial visc (visc) may permit slightly larger CFL but must be used cautiously.\n- Robustness: WENO is more expensive and complex (characteristic projection recommended for best behavior); positivity-preserving/entropy corrections are required to avoid negative density/pressure in under-resolved regions.\n- Complexity: per time-step cost is higher (wider stencils and nonlinear weights) — still O(N) but with larger constants and more memory accesses; parallelization is straightforward but careful halo handling needed.\n- Efficiency tradeoffs: superior accuracy for a given N may permit coarser grids overall, but absolute runtime higher; implement vectorized operations and minimize Python-level loops for performance (pass constants explicitly to compiled-like functions).\n- Use of visc: artificial visc can hide poor resolution; tune visc small and test for convergence. Ensure visc is always called visc in code (never nu) and passed explicitly to all functions that use it.\n- Diagnostics: manual trapezoidal integration must be used to compute global quantities; ensure ghost cells are excluded from integrals.\n\nStop.'



### technical_spec
	Current Stage [A/3]
Finite-volume MUSCL-Hancock (second-order) solver with HLLC approximate Riemann solver and SSP-RK2 time integration.

Architecture:
- Grid and state arrays include NG ghost cells; primitives and conserved are converted via dedicated functions.
- Reconstruction: piecewise-linear MUSCL with minmod limiter computed once per RHS evaluation (vectorized).
- Fluxes: HLLC solver vectorized over interfaces to compute numerical fluxes.
- Time integration: two-stage Heun/SSP-RK2 with positivity floor enforcement between stages.
- Diagnostics: shock position located by max density gradient; global conserved integrals computed with manual trapezoidal integration; total variation reported.
- Outputs: primitives at requested times stored and plotted; three PNG files saved (density, velocity, pressure).

Data flow:
Initial primitives -> conserved U with ghosts -> time-stepping loop:
  apply BCs -> compute slopes -> reconstruct interface states -> HLLC fluxes -> RHS -> RK stages -> positivity fix -> store snapshots -> diagnostics -> continue

All numeric constants (GAMMA, CFL, visc, N_CELL, NG) are set at the top and passed into functions where required.

	Current Stage [B/3]
This script implements a vectorized 1D finite-volume MUSCL-Hancock solver for the Euler equations using a correct conserved-variable definition and a properly implemented HLLC Riemann solver. Data flow: initialize primitives -> convert to conserved (energy density) -> apply transmissive BCs -> compute limited slopes (MUSCL) -> reconstruct interface states -> compute HLLC fluxes at interfaces -> form RHS via explicit flux differences over interior cells -> advance in time with SSP-RK2 (Heun) and positivity floors -> store snapshots and diagnostics. All matrix constructions (grid, index ranges, etc.) are done once outside inner loops where possible. Diagnostics include shock position (from interior gradient), total variation, and relative changes in global conserved integrals computed via manual trapezoidal integration.

	Current Stage [C/3]
Solver overview:
- Finite-volume 1D Euler solver (MUSCL-Hancock) with HLLC Riemann solver and SSP-RK2 time integration.
- Data flow: build grid once -> initialize primitives -> convert to conserved -> apply BCs -> time loop: compute RHS via MUSCL reconstruction + HLLC fluxes -> RK2 update with positivity floors -> snapshot outputs.
- Robustness additions to address instability: immediate clipping of rho/p to safe floors, protected denominators (EPS), and HLL fallback for any invalid HLLC result. Time-step steering avoids dt==0 stalls by explicitly capturing snapshots when remaining time is numerically negligible.
- Diagnostics: shock position, total variation, and manual trapezoidal global conserved integrals are printed at each output. A small single-interface flux test compares HLL and HLLC as a sanity unit test.



### python_codes

#### Script block1:

```python
#!/usr/bin/env python3
"""
MUSCL-Hancock finite-volume solver for 1D Euler equations with HLLC Riemann solver.
Solves Sod shock tube on x in [0,1] with gamma=1.4 and outputs density, velocity,
pressure at t=0.1, 0.2, 0.3. Saves three figures (density, velocity, pressure).
All important constants are defined at top and passed to functions where appropriate.
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import os

# --------------------------- Global constants ---------------------------
GAMMA = 1.4   # specific heat ratio
CFL = 0.5     # CFL number
N_CELL = 400  # number of computational cells (not counting ghost cells)
NG = 2        # number of ghost cells on each side
visc = 0.0    # optional artificial viscosity coefficient (kept as 'visc')
DX = None     # set later after domain build
X0, X1 = 0.0, 1.0  # domain
OUTPUT_TIMES = [0.1, 0.2, 0.3]
RHO_FLOOR = 1e-8
P_FLOOR = 1e-8
MAX_STEPS = 200000

# --------------------------- Utility functions ---------------------------

def manual_trapz(x, y):
    """Manual trapezoidal integration of y over x arrays (both 1D, same length).
    Avoids np.trapz.
    """
    # assumes x is sorted and equally spaced usually; general implementation
    n = len(x)
    if n < 2:
        return 0.0
    s = 0.0
    for i in range(n - 1):
        s += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i])
    return s


def minmod(a, b):
    """Minmod limiter for two arguments (elementwise)"""
    out = np.zeros_like(a)
    mask = (a * b) > 0.0
    out[mask] = np.where(np.abs(a[mask]) < np.abs(b[mask]), a[mask], b[mask])
    return out

# --------------------------- Primitive/Conserved conversions ---------------------------

def prim_to_cons(rho, u, p, gamma):
    E = p / (gamma - 1.0) + 0.5 * rho * u * u
    return np.vstack((rho, rho * u, rho * E))


def cons_to_prim(U, gamma):
    rho = U[0]
    mom = U[1]
    E = U[2]
    # avoid division by zero
    u = np.where(rho > 0.0, mom / rho, 0.0)
    p = (gamma - 1.0) * (E - 0.5 * rho * u * u)
    return rho, u, p

# --------------------------- Boundary conditions ---------------------------

def apply_transmissive_bc(U, ng):
    """Simple zero-gradient (transmissive) BC using ghost cells."""
    # U shape (3, Ntot)
    U[:, :ng] = U[:, ng:ng + 1]  # left ghost cells equal to first interior
    U[:, -ng:] = U[:, -ng - 1:-ng]  # right ghost cells equal to last interior

# --------------------------- Reconstruction (MUSCL) ---------------------------

def compute_slopes(prim, dx):
    """Compute limited slopes for each primitive variable. prim shape (3, Ntot)
    Returns slopes same shape.
    """
    slopes = np.zeros_like(prim)
    # central differences as candidates
    dqL = prim[:, 1:-1] - prim[:, :-2]
    dqR = prim[:, 2:] - prim[:, 1:-1]
    # minmod limiter
    limited = minmod(dqL, dqR)
    slopes[:, 1:-1] = limited
    # set slopes in ghost region to zero
    slopes[:, 0] = slopes[:, 1]
    slopes[:, -1] = slopes[:, -2]
    return slopes

# --------------------------- HLLC Riemann solver ---------------------------

def hllc_flux(UL, UR, gamma):
    """Compute HLLC flux for left/right primitive states.
    UL, UR: primitive arrays [rho, u, p] each of shape (3, M)
    Returns flux array shape (3, M)
    """
    rhoL, uL, pL = UL
    rhoR, uR, pR = UR

    # Conserved
    EL = pL / (gamma - 1.0) + 0.5 * rhoL * uL * uL
    ER = pR / (gamma - 1.0) + 0.5 * rhoR * uR * uR
    ULc = np.vstack((rhoL, rhoL * uL, rhoL * EL))
    URc = np.vstack((rhoR, rhoR * uR, rhoR * ER))

    # Fluxes
    FL = np.vstack((rhoL * uL,
                    rhoL * uL * uL + pL,
                    (rhoL * EL + pL) * uL))
    FR = np.vstack((rhoR * uR,
                    rhoR * uR * uR + pR,
                    (rhoR * ER + pR) * uR))

    # sound speeds
    cL = np.sqrt(np.maximum(gamma * pL / rhoL, 0.0))
    cR = np.sqrt(np.maximum(gamma * pR / rhoR, 0.0))

    # wave speed estimates (Davis/Roe-type estimate)
    SL = np.minimum(uL - cL, uR - cR)
    SR = np.maximum(uL + cL, uR + cR)

    # avoid zero denominator
    denom = (rhoL * (SL - uL) - rhoR * (SR - uR))
    # compute S_star safely
    S_star = np.zeros_like(SL)
    mask = denom != 0.0
    S_star[mask] = (pR[mask] - pL[mask] + rhoL[mask] * uL[mask] * (SL[mask] - uL[mask]) -
                    rhoR[mask] * uR[mask] * (SR[mask] - uR[mask])) / denom[mask]
    # initialize flux
    F = np.zeros_like(FL)

    # Case 1: SL >= 0 -> flux = FL
    mask1 = SL >= 0.0
    if np.any(mask1):
        F[:, mask1] = FL[:, mask1]

    # Case 4: SR <= 0 -> flux = FR
    mask4 = SR <= 0.0
    if np.any(mask4):
        F[:, mask4] = FR[:, mask4]

    # Case 2: SL < 0 <= S_star
    mask2 = (~mask1) & (~mask4) & (S_star >= 0.0)
    if np.any(mask2):
        # compute U_star_L
        rhoSL = rhoL[mask2] * (SL[mask2] - uL[mask2]) / (SL[mask2] - S_star[mask2])
        uSL = S_star[mask2]
        # energy star (Toro's formula)
        ESL = EL[mask2] + (S_star[mask2] - uL[mask2]) * (
            S_star[mask2] + pL[mask2] / (rhoL[mask2] * (SL[mask2] - uL[mask2])))
        UstL = np.vstack((rhoSL,
                          rhoSL * uSL,
                          rhoSL * ESL))
        F[:, mask2] = FL[:, mask2] + SL[mask2] * (UstL - ULc[:, mask2])

    # Case 3: S_star < 0 <= SR
    mask3 = (~mask1) & (~mask4) & (S_star < 0.0)
    if np.any(mask3):
        rhoSR = rhoR[mask3] * (SR[mask3] - uR[mask3]) / (SR[mask3] - S_star[mask3])
        uSR = S_star[mask3]
        ESR = ER[mask3] + (S_star[mask3] - uR[mask3]) * (
            S_star[mask3] + pR[mask3] / (rhoR[mask3] * (SR[mask3] - uR[mask3])))
        UstR = np.vstack((rhoSR,
                          rhoSR * uSR,
                          rhoSR * ESR))
        F[:, mask3] = FR[:, mask3] + SR[mask3] * (UstR - URc[:, mask3])

    return F

# --------------------------- RHS evaluator ---------------------------

def compute_rhs(U, dx, gamma, ng, visc_coeff):
    """Compute right-hand side (time derivative) for conserved variables.
    U shape (3, Ntot)
    returns dUdt shape (3, Ntot_interior) (excluding ghosts)
    """
    # Apply BCs first
    apply_transmissive_bc(U, ng)

    # primitives
    rho, u, p = cons_to_prim(U, gamma)
    prim = np.vstack((rho, u, p))

    # slopes and reconstructed states at interfaces
    slopes = compute_slopes(prim, dx)
    Ntot = prim.shape[1]
    # interfaces count = Ntot-1
    # left state at i+1/2 is from cell i: prim[:, i] + 0.5*slopes[:, i]
    # right state at i+1/2 from cell i+1: prim[:, i+1] - 0.5*slopes[:, i+1]
    UL = prim[:, :-1] + 0.5 * slopes[:, :-1]
    UR = prim[:, 1:] - 0.5 * slopes[:, 1:]

    # compute fluxes at interfaces via HLLC (vectorized over interfaces)
    F_int = hllc_flux(UL, UR, gamma)

    # compute conservative RHS: -(F_{i+1/2} - F_{i-1/2}) / dx for interior cells
    dU = np.zeros_like(U)
    # flux differences across interfaces -> assign to cells 0..Ntot-1 using interfaces 0..Ntot-2
    # interior cells i use F_int[:, i+1] and F_int[:, i]
    # so dU[:, i] = -(F_int[:, i+1] - F_int[:, i]) / dx
    dU[:, :-1] += - (F_int - np.roll(F_int, 1, axis=1))[:, :-1] / dx

    # The above uses roll which includes wrap-around; but we only trust interior cells; ensure
    # dU for ghost cells will be ignored downstream

    # optional artificial viscosity term (second-derivative) in conserved variables
    if visc_coeff > 0.0:
        # second difference for each conserved component
        lap = np.zeros_like(U)
        lap[:, 1:-1] = U[:, 2:] - 2.0 * U[:, 1:-1] + U[:, :-2]
        dU[:, :] += visc_coeff * lap / (dx * dx)

    # Return dU for interior cells only (exclude ghost cells)
    return dU

# --------------------------- Initial condition ---------------------------

def initial_condition(xc, gamma):
    rho0 = np.where(xc <= 0.5, 1.0, 0.125)
    u0 = np.zeros_like(xc)
    p0 = np.where(xc <= 0.5, 1.0, 0.1)
    U0 = prim_to_cons(rho0, u0, p0, gamma)
    return U0

# --------------------------- Shock locator and diagnostics ---------------------------

def shock_position(xc, rho, ng):
    # compute gradient magnitude excluding ghosts
    grad = np.abs(np.diff(rho))
    # indices relative to full array correspond to interfaces
    imax = np.argmax(grad[ng - 1: -ng]) + (ng - 1)
    # approximate shock position between cells imax and imax+1
    xpos = 0.5 * (xc[imax] + xc[imax + 1])
    return xpos, grad[imax]

# --------------------------- Main solver ---------------------------

def run_solver(N_cell, ng, gamma, cfl, dx_user, output_times, visc_coeff):
    # build grid including ghost cells
    Ntot = N_cell + 2 * ng
    # cell centers
    dx = (X1 - X0) / N_cell
    x_centers = np.linspace(X0 + 0.5 * dx, X1 - 0.5 * dx, N_cell)
    # expand to include ghost centers: we create extended array of length Ntot
    x_ext = np.zeros(Ntot)
    x_ext[ng:ng + N_cell] = x_centers
    # fill ghost centers by extrapolation
    for i in range(ng):
        x_ext[i] = x_ext[ng] - (ng - i) * dx
        x_ext[-1 - i] = x_ext[ng + N_cell - 1] + (i + 1) * dx

    # initial condition
    U = np.zeros((3, Ntot))
    U[:, ng:ng + N_cell] = initial_condition(x_centers, gamma)
    # apply BCs
    apply_transmissive_bc(U, ng)

    # store initial global conserved integrals
    rho0, u0, p0 = cons_to_prim(U, gamma)
    mass0 = manual_trapz(x_ext[ng:ng + N_cell], rho0[ng:ng + N_cell])
    mom0 = manual_trapz(x_ext[ng:ng + N_cell], (rho0 * u0)[ng:ng + N_cell])
    energy0 = manual_trapz(x_ext[ng:ng + N_cell], (U[2, ng:ng + N_cell]))

    t = 0.0
    step = 0
    out_index = 0
    outputs = []  # list of (t, x_centers, rho, u, p)
    shock_positions = []

    # max allowed steps guard
    while (out_index < len(output_times)) and (step < MAX_STEPS):
        # compute primitives for CFL
        rho, u, p = cons_to_prim(U, gamma)
        c = np.sqrt(np.maximum(gamma * p / rho, 0.0))
        max_speed = np.nanmax(np.abs(u) + c)
        if max_speed <= 0:
            dt = 1e-6
        else:
            dt = cfl * dx / max_speed

        # don't step past next output time
        t_next_output = output_times[out_index]
        if t + dt > t_next_output:
            dt = t_next_output - t

        # Stage 1: forward Euler
        rhs = compute_rhs(U, dx, gamma, ng, visc_coeff)
        U1 = U + dt * rhs
        # enforce positivity floor in primitives after stage
        rho1, u1, p1 = cons_to_prim(U1, gamma)
        p1 = np.maximum(p1, P_FLOOR)
        rho1 = np.maximum(rho1, RHO_FLOOR)
        U1[:, :] = prim_to_cons(rho1, u1, p1, gamma)

        # Stage 2: Heun / explicit trapezoidal
        rhs1 = compute_rhs(U1, dx, gamma, ng, visc_coeff)
        U_new = U + 0.5 * dt * (rhs + rhs1)

        # positivity fix
        rho_new, u_new, p_new = cons_to_prim(U_new, gamma)
        p_new = np.maximum(p_new, P_FLOOR)
        rho_new = np.maximum(rho_new, RHO_FLOOR)
        U[:, :] = prim_to_cons(rho_new, u_new, p_new, gamma)

        t += dt
        step += 1

        # capture outputs at specified times
        if abs(t - t_next_output) < 1e-12 or t >= t_next_output:
            # extract interior primitives
            rho_i, u_i, p_i = cons_to_prim(U, gamma)
            ri = rho_i[ng:ng + N_cell].copy()
            ui = u_i[ng:ng + N_cell].copy()
            pi = p_i[ng:ng + N_cell].copy()
            outputs.append((t, x_centers.copy(), ri, ui, pi))

            # diagnostics: shock position
            xpos, gradmax = shock_position(x_ext, rho_i, ng)
            shock_positions.append((t, xpos, gradmax))

            # global conserved integrals
            mass = manual_trapz(x_ext[ng:ng + N_cell], ri)
            momentum = manual_trapz(x_ext[ng:ng + N_cell], ri * ui)
            energy = manual_trapz(x_ext[ng:ng + N_cell], U[2, ng:ng + N_cell])

            tv = np.sum(np.abs(np.diff(ri)))

            print(f"[Output t={t:.4f}] step={step} dt={dt:.4e}")
            print(f"  Shock pos ~ {xpos:.6f} (grad max {gradmax:.4e})")
            print(f"  Mass: {mass:.6f} (rel Δ {(mass - mass0)/mass0:.3e})")
            print(f"  Momentum: {momentum:.6f} (rel Δ {(momentum - mom0)/max(1e-12, abs(mom0)):.3e})")
            print(f"  Energy: {energy:.6f} (rel Δ {(energy - energy0)/energy0:.3e})")
            print(f"  TV(rho): {tv:.6e}")
            print(f"  rho range: [{ri.min():.4e}, {ri.max():.4e}], p range: [{pi.min():.4e}, {pi.max():.4e}]\n")

            out_index += 1

        # safeguard for time floating error
        if t >= output_times[-1] - 1e-12:
            break

    return outputs, shock_positions

# --------------------------- Plotting ---------------------------

def save_plots(outputs, out_dir='plots'):
    os.makedirs(out_dir, exist_ok=True)
    times = [o[0] for o in outputs]
    # density figure
    plt.figure(figsize=(8, 4))
    for t, x, rho, u, p in outputs:
        plt.plot(x, rho, label=f't={t:.2f}')
    plt.xlabel('x')
    plt.ylabel('Density')
    plt.title('Density at output times')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    fname = os.path.join(out_dir, 'density.png')
    plt.savefig(fname)
    plt.close()

    # velocity figure
    plt.figure(figsize=(8, 4))
    for t, x, rho, u, p in outputs:
        plt.plot(x, u, label=f't={t:.2f}')
    plt.xlabel('x')
    plt.ylabel('Velocity')
    plt.title('Velocity at output times')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    fname = os.path.join(out_dir, 'velocity.png')
    plt.savefig(fname)
    plt.close()

    # pressure figure
    plt.figure(figsize=(8, 4))
    for t, x, rho, u, p in outputs:
        plt.plot(x, p, label=f't={t:.2f}')
    plt.xlabel('x')
    plt.ylabel('Pressure')
    plt.title('Pressure at output times')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    fname = os.path.join(out_dir, 'pressure.png')
    plt.savefig(fname)
    plt.close()

# --------------------------- Run and save outputs ---------------------------

if __name__ == '__main__':
    print('Running 1D Euler MUSCL-Hancock (MUSCL + HLLC, SSP RK2) solver...')
    outputs, shock_positions = run_solver(N_CELL, NG, GAMMA, CFL, DX, OUTPUT_TIMES, visc)
    print('Finished time integration. Saving plots...')
    save_plots(outputs)
    print('Plots saved in ./plots: density.png, velocity.png, pressure.png')
    print('\nShock positions history:')
    for t, xpos, g in shock_positions:
        print(f'  t={t:.4f} -> x_shock~{xpos:.6f} (grad {g:.3e})')

    # Basic spurious oscillation check: count negative / NaN values in final snapshot
    if outputs:
        t, x, rho, u, p = outputs[-1]
        bad = np.sum((rho <= 0) | (p <= 0) | np.isnan(rho) | np.isnan(p))
        if bad == 0:
            print('\nNo negative densities or pressures detected in final snapshot.')
        else:
            print(f'\nWarning: {bad} cells with non-physical values in final snapshot.')

```


#### Script block2:

#!/usr/bin/env python3
"""
Corrected MUSCL-Hancock finite-volume solver for 1D Euler equations with HLLC Riemann solver.
Sod shock tube on x in [0,1], gamma=1.4. Saves three figures (density, velocity, pressure)
at times 0.1, 0.2, 0.3 in './plots'.

Key fixes from previous version:
- Consistent conserved variable definition: U = [rho, rho*u, E_total] where E_total = p/(gamma-1) + 0.5*rho*u^2 (energy density).
- Correct, vectorized HLLC implementation using energy density and correct flux forms.
- Flux divergence computed with explicit indexing (no np.roll), updating only interior cells to avoid wraparound and shape mismatches.
- compute_rhs returns dU with same shape as U; non-zero entries only in interior cells.
- Transmissive BCs and shock locator operate robustly on interior-only data.
- Manual trapezoidal integration implemented (no np.trapz usage).

The code is modular and passes all essential constants explicitly where applicable.
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# --------------------------- Global constants ---------------------------
GAMMA = 1.4   # specific heat ratio
CFL = 0.5     # CFL number
N_CELL = 400  # number of computational cells (interior)
NG = 2        # number of ghost cells on each side
visc = 0.0    # artificial viscosity coefficient (kept name 'visc')
X0, X1 = 0.0, 1.0  # domain
OUTPUT_TIMES = [0.1, 0.2, 0.3]
RHO_FLOOR = 1e-8
P_FLOOR = 1e-8
MAX_STEPS = 200000

# --------------------------- Utility functions ---------------------------

def manual_trapz(x, y):
    """Manual trapezoidal integration of y over x arrays (1D, same length)."""
    n = len(x)
    if n < 2:
        return 0.0
    s = 0.0
    for i in range(n - 1):
        s += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i])
    return s


def minmod(a, b):
    """Elementwise minmod limiter for arrays a and b."""
    out = np.zeros_like(a)
    mask = (a * b) > 0.0
    out[mask] = np.where(np.abs(a[mask]) < np.abs(b[mask]), a[mask], b[mask])
    return out

# --------------------------- Primitive/Conserved conversions ---------------------------

def prim_to_cons(rho, u, p, gamma):
    """Convert primitive (rho,u,p) to conserved U = [rho, rho*u, E_total]
    where E_total is energy density = p/(gamma-1) + 0.5 * rho * u^2.
    Each input can be 1D arrays; returns U shape (3, N).
    """
    E_total = p / (gamma - 1.0) + 0.5 * rho * u * u
    return np.vstack((rho, rho * u, E_total))


def cons_to_prim(U, gamma):
    """Convert conserved U (shape (3,N)) to primitives (rho,u,p).
    Enforces safe division for u and returns arrays of shape (N,).
    """
    rho = U[0].copy()
    mom = U[1].copy()
    E_total = U[2].copy()
    # avoid division by zero
    u = np.zeros_like(rho)
    nonzero = rho > 0.0
    u[nonzero] = mom[nonzero] / rho[nonzero]
    # pressure from energy density
    p = (gamma - 1.0) * (E_total - 0.5 * rho * u * u)
    return rho, u, p

# --------------------------- Boundary conditions ---------------------------

def apply_transmissive_bc(U, ng):
    """Zero-gradient (transmissive) BC using ghost cells. U shape (3, Ntot)."""
    # Left ghosts replicate first interior cell
    U[:, :ng] = np.repeat(U[:, ng:ng + 1], ng, axis=1)
    # Right ghosts replicate last interior cell
    U[:, -ng:] = np.repeat(U[:, -(ng + 1):-ng], ng, axis=1)

# --------------------------- Reconstruction (MUSCL) ---------------------------

def compute_slopes(prim):
    """Compute limited slopes for each primitive variable using minmod.
    prim shape (3, Ntot); returns slopes same shape. Limited at interior cells.
    """
    slopes = np.zeros_like(prim)
    # central differences candidates for interior points 1..Ntot-2
    dqL = prim[:, 1:-1] - prim[:, :-2]
    dqR = prim[:, 2:] - prim[:, 1:-1]
    limited = minmod(dqL, dqR)
    slopes[:, 1:-1] = limited
    # for boundary cells, copy nearest interior slope (keeps array well-behaved)
    slopes[:, 0] = slopes[:, 1]
    slopes[:, -1] = slopes[:, -2]
    return slopes

# --------------------------- HLLC Riemann solver ---------------------------

def hllc_flux(UL, UR, gamma):
    """Compute HLLC flux for left/right primitive states UL, UR.
    UL, UR: primitive arrays [rho, u, p] each of shape (3, M)
    Returns flux array shape (3, M).

    The conserved variables are U = [rho, rho*u, E_total] where E_total is energy density.
    Flux: F = [rho*u, rho*u^2 + p, (E_total + p)*u]
    """
    rhoL, uL, pL = UL
    rhoR, uR, pR = UR

    # Energies (energy density)
    EL = pL / (gamma - 1.0) + 0.5 * rhoL * uL * uL
    ER = pR / (gamma - 1.0) + 0.5 * rhoR * uR * uR

    # Conserved vectors (per volume)
    ULc = np.vstack((rhoL, rhoL * uL, EL))
    URc = np.vstack((rhoR, rhoR * uR, ER))

    # Fluxes for left and right states
    FL = np.vstack((rhoL * uL,
                    rhoL * uL * uL + pL,
                    (EL + pL) * uL))
    FR = np.vstack((rhoR * uR,
                    rhoR * uR * uR + pR,
                    (ER + pR) * uR))

    # sound speeds
    cL = np.sqrt(np.maximum(gamma * pL / rhoL, 0.0))
    cR = np.sqrt(np.maximum(gamma * pR / rhoR, 0.0))

    # wave speed estimates
    SL = np.minimum(uL - cL, uR - cR)
    SR = np.maximum(uL + cL, uR + cR)

    # S_star (contact wave) standard formula
    denom = (rhoL * (SL - uL) - rhoR * (SR - uR))
    S_star = np.zeros_like(SL)
    valid = denom != 0.0
    if np.any(valid):
        S_star[valid] = (
            pR[valid] - pL[valid] + rhoL[valid] * uL[valid] * (SL[valid] - uL[valid])
            - rhoR[valid] * uR[valid] * (SR[valid] - uR[valid])
        ) / denom[valid]

    M = UL.shape[1]
    F = np.zeros_like(FL)

    # Regions:
    # If SL >= 0 -> use left flux
    maskL = SL >= 0.0
    if np.any(maskL):
        F[:, maskL] = FL[:, maskL]

    # If SR <= 0 -> use right flux
    maskR = SR <= 0.0
    if np.any(maskR):
        F[:, maskR] = FR[:, maskR]

    # Middle regions: SL < 0 < SR
    middle = (~maskL) & (~maskR)
    if np.any(middle):
        # split middle into left-star (SL < 0 <= S_star) and right-star (S_star < 0 <= SR)
        left_star = middle & (S_star >= 0.0)
        right_star = middle & (S_star < 0.0)

        if np.any(left_star):
            idx = left_star
            # rho in star left
            rhoSL = rhoL[idx] * (SL[idx] - uL[idx]) / (SL[idx] - S_star[idx])
            uSL = S_star[idx]
            # specific energy e = E/rho
            eL = EL[idx] / rhoL[idx]
            # compute specific energy in star
            denom_inner = rhoL[idx] * (SL[idx] - uL[idx])
            # avoid division by zero - small mask handled by numerator/denom previously
            eSL = eL + (S_star[idx] - uL[idx]) * (
                S_star[idx] + pL[idx] / np.maximum(denom_inner, 1e-16)
            )
            ESL = rhoSL * eSL
            UstL = np.vstack((rhoSL, rhoSL * uSL, ESL))
            F[:, idx] = FL[:, idx] + SL[idx] * (UstL - ULc[:, idx])

        if np.any(right_star):
            idx = right_star
            rhoSR = rhoR[idx] * (SR[idx] - uR[idx]) / (SR[idx] - S_star[idx])
            uSR = S_star[idx]
            eR = ER[idx] / rhoR[idx]
            denom_inner = rhoR[idx] * (SR[idx] - uR[idx])
            eSR = eR + (S_star[idx] - uR[idx]) * (
                S_star[idx] + pR[idx] / np.maximum(denom_inner, 1e-16)
            )
            ESR = rhoSR * eSR
            UstR = np.vstack((rhoSR, rhoSR * uSR, ESR))
            F[:, idx] = FR[:, idx] + SR[idx] * (UstR - URc[:, idx])

    return F

# --------------------------- RHS evaluator ---------------------------

def compute_rhs(U, dx, gamma, ng, visc_coeff):
    """Compute RHS dU/dt with same shape as U. Only interior cells (ng:ng+N_cell)
    will have non-zero updates. Uses transmissive BCs internally.
    """
    # Apply BCs first
    apply_transmissive_bc(U, ng)

    # primitives
    rho, u, p = cons_to_prim(U, gamma)
    prim = np.vstack((rho, u, p))

    # slopes and reconstructed states at interfaces
    slopes = compute_slopes(prim)
    Ntot = prim.shape[1]

    # left and right states at interfaces i+1/2 (interfaces indices 0..Ntot-2)
    UL = prim[:, :-1] + 0.5 * slopes[:, :-1]
    UR = prim[:, 1:] - 0.5 * slopes[:, 1:]

    # fluxes at interfaces via HLLC
    F_int = hllc_flux(UL, UR, gamma)  # shape (3, Ntot-1)

    # initialize dU
    dU = np.zeros_like(U)

    # compute flux divergence only for interior computational cells i in [ng, ng+N_cell-1]
    # For cell i, right interface is index i, left interface is i-1 in F_int indexing
    # so dU[:, i] = - (F_int[:, i] - F_int[:, i-1]) / dx
    i0 = ng
    i1 = Ntot - ng  # exclusive upper bound
    # ensure indices valid: F_int has interfaces 0..Ntot-2 so i in [1..Ntot-2]
    # we only update cells i0..i1-1 inclusive
    dU[:, i0:i1] = - (F_int[:, i0:i1] - F_int[:, i0 - 1:i1 - 1]) / dx

    # optional artificial viscosity term (second-derivative) in conserved variables
    if visc_coeff > 0.0:
        lap = np.zeros_like(U)
        lap[:, 1:-1] = U[:, 2:] - 2.0 * U[:, 1:-1] + U[:, :-2]
        dU[:, :] += visc_coeff * lap / (dx * dx)

    return dU

# --------------------------- Initial condition ---------------------------

def initial_condition(xc, gamma):
    rho0 = np.where(xc <= 0.5, 1.0, 0.125)
    u0 = np.zeros_like(xc)
    p0 = np.where(xc <= 0.5, 1.0, 0.1)
    U0 = prim_to_cons(rho0, u0, p0, gamma)
    return U0

# --------------------------- Shock locator and diagnostics ---------------------------

def shock_position(xc_centers, rho, ng, N_cell):
    """Locate shock position from max interior gradient of rho.
    rho is full array with ghosts; xc_centers are full centers including ghosts.
    Returns approximate x position and the maximum gradient value.
    """
    rho_int = rho[ng:ng + N_cell]
    x_int = xc_centers[ng:ng + N_cell]
    grad = np.abs(np.diff(rho_int))
    if len(grad) == 0:
        return 0.0, 0.0
    imax = np.argmax(grad)
    # shock between cell imax and imax+1 -> in full indexing it's ng+imax and ng+imax+1
    xpos = 0.5 * (x_int[imax] + x_int[imax + 1])
    return float(xpos), float(grad[imax])

# --------------------------- Main solver ---------------------------

def run_solver(N_cell, ng, gamma, cfl, output_times, visc_coeff):
    """Run time integration and return stored outputs and shock positions.
    All arrays include ghost cells; outputs store interior-only primitive arrays.
    """
    Ntot = N_cell + 2 * ng
    dx = (X1 - X0) / N_cell
    # interior cell centers
    x_centers = np.linspace(X0 + 0.5 * dx, X1 - 0.5 * dx, N_cell)
    # extended centers including ghosts
    x_ext = np.zeros(Ntot)
    x_ext[ng:ng + N_cell] = x_centers
    # fill ghost centers by extrapolation
    for i in range(ng):
        x_ext[i] = x_ext[ng] - (ng - i) * dx
        x_ext[-1 - i] = x_ext[ng + N_cell - 1] + (i + 1) * dx

    # initial condition and U with ghosts
    U = np.zeros((3, Ntot))
    U[:, ng:ng + N_cell] = initial_condition(x_centers, gamma)
    apply_transmissive_bc(U, ng)

    # store initial conserved integrals for diagnostics
    rho0, u0, p0 = cons_to_prim(U, gamma)
    mass0 = manual_trapz(x_ext[ng:ng + N_cell], rho0[ng:ng + N_cell])
    mom0 = manual_trapz(x_ext[ng:ng + N_cell], (rho0 * u0)[ng:ng + N_cell])
    energy0 = manual_trapz(x_ext[ng:ng + N_cell], (U[2, ng:ng + N_cell]))

    t = 0.0
    step = 0
    out_index = 0
    outputs = []  # list of (t, x_centers, rho, u, p)
    shock_positions = []

    # time loop
    while (out_index < len(output_times)) and (step < MAX_STEPS):
        # primitives for CFL
        rho, u, p = cons_to_prim(U, gamma)
        c = np.sqrt(np.maximum(gamma * p / np.maximum(rho, 1e-16), 0.0))
        max_speed = np.nanmax(np.abs(u) + c)
        if max_speed <= 0:
            dt = 1e-6
        else:
            dt = cfl * dx / max_speed

        # do not step past next output time
        t_next_output = output_times[out_index]
        if t + dt > t_next_output:
            dt = t_next_output - t
            if dt <= 1e-14:
                # if extremely small due to roundoff, force capture below without advancing
                dt = 0.0

        # Stage 1
        rhs = compute_rhs(U, dx, gamma, ng, visc_coeff)
        U1 = U + dt * rhs
        # enforce positivity floors in primitives after stage 1
        rho1, u1, p1 = cons_to_prim(U1, gamma)
        p1 = np.maximum(p1, P_FLOOR)
        rho1 = np.maximum(rho1, RHO_FLOOR)
        U1[:, :] = prim_to_cons(rho1, u1, p1, gamma)

        # Stage 2
        rhs1 = compute_rhs(U1, dx, gamma, ng, visc_coeff)
        U_new = U + 0.5 * dt * (rhs + rhs1)

        # positivity fix
        rho_new, u_new, p_new = cons_to_prim(U_new, gamma)
        p_new = np.maximum(p_new, P_FLOOR)
        rho_new = np.maximum(rho_new, RHO_FLOOR)
        U[:, :] = prim_to_cons(rho_new, u_new, p_new, gamma)

        t += dt
        step += 1

        # capture outputs at specified times (or when dt==0 due to rounding)
        if (abs(t - t_next_output) < 1e-12) or (t >= t_next_output) or (dt == 0.0):
            rho_i, u_i, p_i = cons_to_prim(U, gamma)
            ri = rho_i[ng:ng + N_cell].copy()
            ui = u_i[ng:ng + N_cell].copy()
            pi = p_i[ng:ng + N_cell].copy()
            outputs.append((t, x_centers.copy(), ri, ui, pi))

            # diagnostics
            xpos, gradmax = shock_position(x_ext, rho_i, ng, N_cell)
            shock_positions.append((t, xpos, gradmax))

            mass = manual_trapz(x_ext[ng:ng + N_cell], ri)
            momentum = manual_trapz(x_ext[ng:ng + N_cell], ri * ui)
            energy = manual_trapz(x_ext[ng:ng + N_cell], U[2, ng:ng + N_cell])

            tv = np.sum(np.abs(np.diff(ri)))

            print(f"[Output t={t:.4f}] step={step} dt={dt:.4e}")
            print(f"  Shock pos ~ {xpos:.6f} (grad max {gradmax:.4e})")
            print(f"  Mass: {mass:.6f} (rel Δ {(mass - mass0)/mass0:.3e})")
            print(f"  Momentum: {momentum:.6f} (rel Δ {(momentum - mom0)/max(1e-12, abs(mom0)):.3e})")
            print(f"  Energy: {energy:.6f} (rel Δ {(energy - energy0)/energy0:.3e})")
            print(f"  TV(rho): {tv:.6e}")
            print(f"  rho range: [{ri.min():.4e}, {ri.max():.4e}], p range: [{pi.min():.4e}, {pi.max():.4e}]\n")

            out_index += 1

        # break if we've reached final output time (guard for floating point)
        if t >= output_times[-1] - 1e-12:
            break

    return outputs, shock_positions

# --------------------------- Plotting ---------------------------

def save_plots(outputs, out_dir='plots'):
    os.makedirs(out_dir, exist_ok=True)

    # density
    plt.figure(figsize=(8, 4))
    for t, x, rho, u, p in outputs:
        plt.plot(x, rho, label=f't={t:.2f}')
    plt.xlabel('x')
    plt.ylabel('Density')
    plt.title('Density at output times')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'density.png'))
    plt.close()

    # velocity
    plt.figure(figsize=(8, 4))
    for t, x, rho, u, p in outputs:
        plt.plot(x, u, label=f't={t:.2f}')
    plt.xlabel('x')
    plt.ylabel('Velocity')
    plt.title('Velocity at output times')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'velocity.png'))
    plt.close()

    # pressure
    plt.figure(figsize=(8, 4))
    for t, x, rho, u, p in outputs:
        plt.plot(x, p, label=f't={t:.2f}')
    plt.xlabel('x')
    plt.ylabel('Pressure')
    plt.title('Pressure at output times')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'pressure.png'))
    plt.close()

# --------------------------- Run and save outputs ---------------------------

if __name__ == '__main__':
    print('Running 1D Euler MUSCL-Hancock (MUSCL + HLLC, SSP RK2) solver...')
    outputs, shock_positions = run_solver(N_CELL, NG, GAMMA, CFL, OUTPUT_TIMES, visc)
    print('Finished time integration. Saving plots...')
    save_plots(outputs)
    print('Plots saved in ./plots: density.png, velocity.png, pressure.png')
    print('\nShock positions history:')
    for t, xpos, g in shock_positions:
        print(f'  t={t:.4f} -> x_shock~{xpos:.6f} (grad {g:.3e})')

    # Basic spurious oscillation check: detect non-physical values in final snapshot
    if outputs:
        t, x, rho, u, p = outputs[-1]
        bad = np.sum((rho <= 0) | (p <= 0) | np.isnan(rho) | np.isnan(p))
        if bad == 0:
            print('\nNo negative densities or pressures detected in final snapshot.')
        else:
            print(f'\nWarning: {bad} cells with non-physical values in final snapshot.')



#### Script block3:

```python
#!/usr/bin/env python3
"""
Robust MUSCL-Hancock finite-volume 1D Euler solver with HLLC (with HLL fallback).
Sod shock tube on x in [0,1], saves density/velocity/pressure plots at output times.
Key robustness measures added per review:
 - Clip rho and p to floors immediately after primitive recovery
 - Protect denominators with EPS and use vectorized masks
 - If HLLC produces invalid values, fallback to HLL flux
 - Avoid dt==0 stall by capturing snapshots when remaining time is numerically tiny
 - Manual trapezoidal integration (no np.trapz)
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# --------------------------- Global constants ---------------------------
GAMMA = 1.4
CFL = 0.5
N_CELL = 400
NG = 2
visc = 0.0  # artificial viscosity coefficient
X0, X1 = 0.0, 1.0
OUTPUT_TIMES = [0.1, 0.2, 0.3]
RHO_FLOOR = 1e-8
P_FLOOR = 1e-8
EPS = 1e-12  # small number to protect denom/sqrt
MAX_STEPS = 200000
MIN_DT = 1e-14  # minimum meaningful dt for stepping

# --------------------------- Utility functions ---------------------------

def manual_trapz(x, y):
    """Manual trapezoidal integration over 1D arrays."""
    n = len(x)
    if n < 2:
        return 0.0
    s = 0.0
    for i in range(n - 1):
        s += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i])
    return s


def minmod(a, b):
    """Elementwise minmod limiter."""
    out = np.zeros_like(a)
    mask = (a * b) > 0.0
    out[mask] = np.where(np.abs(a[mask]) < np.abs(b[mask]), a[mask], b[mask])
    return out

# --------------------------- Primitive/Conserved conversions ---------------------------

def prim_to_cons(rho, u, p, gamma):
    """Convert primitives to conserved (energy density form).
    Inputs are 1D arrays of same length; returns U shape (3,N).
    """
    E_total = p / (gamma - 1.0) + 0.5 * rho * u * u
    return np.vstack((rho, rho * u, E_total))


def cons_to_prim(U, gamma, rho_floor=RHO_FLOOR, p_floor=P_FLOOR, log_warnings=False):
    """Convert conserved U=(3,N) to primitives (rho,u,p). Clips rho/p to floors.
    Returns (rho,u,p) arrays.
    """
    rho = U[0].copy()
    mom = U[1].copy()
    E_total = U[2].copy()

    # safe velocity division
    u = np.zeros_like(rho)
    nz = rho > 0.0
    u[nz] = mom[nz] / rho[nz]

    # pressure from energy density
    p = (gamma - 1.0) * (E_total - 0.5 * rho * u * u)

    # detect non-physical values and clip
    clipped = False
    if np.any(rho < rho_floor) or np.any(p < p_floor):
        clipped = True
    rho = np.maximum(rho, rho_floor)
    p = np.maximum(p, p_floor)

    if log_warnings and clipped:
        print("[cons_to_prim] Applied floors: clipped density/pressure to floors to maintain positivity.")

    return rho, u, p

# --------------------------- Boundary conditions ---------------------------

def apply_transmissive_bc(U, ng):
    """Apply zero-gradient transmissive BC by copying nearest interior cell into ghosts.
    U has shape (3, Ntot).
    """
    # left ghosts
    U[:, :ng] = np.repeat(U[:, ng:ng + 1], ng, axis=1)
    # right ghosts
    U[:, -ng:] = np.repeat(U[:, -(ng + 1):-ng], ng, axis=1)

# --------------------------- Reconstruction (MUSCL) ---------------------------

def compute_slopes(prim):
    """Compute limited slopes (minmod) for primitives prim shape (3,Ntot).
    Returns slopes same shape.
    """
    slopes = np.zeros_like(prim)
    dqL = prim[:, 1:-1] - prim[:, :-2]
    dqR = prim[:, 2:] - prim[:, 1:-1]
    limited = minmod(dqL, dqR)
    slopes[:, 1:-1] = limited
    slopes[:, 0] = slopes[:, 1]
    slopes[:, -1] = slopes[:, -2]
    return slopes

# --------------------------- Riemann solvers: HLL & HLLC (with robust protections) ---------------------------

def hll_flux(ULc, URc, FL, FR, SL, SR):
    """Compute HLL flux robustly. Arrays shape (3, M) for vectors and (M,) for speeds."""
    denom = SR - SL
    denom_safe = denom.copy()
    denom_safe[np.abs(denom_safe) < EPS] = EPS
    # broadcast speeds to shape (3, M)
    Fhll = (SR * FL - SL * FR + SL * SR * (URc - ULc)) / denom_safe
    return Fhll


def hllc_flux(UL, UR, gamma):
    """Vectorized HLLC flux with protections and HLL fallback.
    UL, UR: primitive arrays (rho,u,p) each shape (3,M)
    Returns flux array shape (3,M).
    """
    rhoL, uL, pL = UL
    rhoR, uR, pR = UR

    # clip inputs to floors to avoid negative/zero densities or pressures
    rhoL = np.maximum(rhoL, RHO_FLOOR)
    rhoR = np.maximum(rhoR, RHO_FLOOR)
    pL = np.maximum(pL, P_FLOOR)
    pR = np.maximum(pR, P_FLOOR)

    # energies (per volume)
    EL = pL / (gamma - 1.0) + 0.5 * rhoL * uL * uL
    ER = pR / (gamma - 1.0) + 0.5 * rhoR * uR * uR

    ULc = np.vstack((rhoL, rhoL * uL, EL))
    URc = np.vstack((rhoR, rhoR * uR, ER))

    FL = np.vstack((rhoL * uL,
                    rhoL * uL * uL + pL,
                    (EL + pL) * uL))
    FR = np.vstack((rhoR * uR,
                    rhoR * uR * uR + pR,
                    (ER + pR) * uR))

    cL = np.sqrt(np.maximum(gamma * pL / rhoL, 0.0))
    cR = np.sqrt(np.maximum(gamma * pR / rhoR, 0.0))

    # wave speed estimates
    SL = np.minimum(uL - cL, uR - cR)
    SR = np.maximum(uL + cL, uR + cR)

    # ensure arrays
    M = UL.shape[1]
    F = np.zeros_like(FL)

    # masks
    maskL = SL >= 0.0
    maskR = SR <= 0.0
    middle = (~maskL) & (~maskR)

    # safe denom for S_star
    denom = rhoL * (SL - uL) - rhoR * (SR - uR)
    denom_safe = np.where(np.abs(denom) < EPS, np.sign(denom) * EPS + EPS, denom)

    S_star = np.zeros_like(SL)
    # compute S_star only where middle is True
    if np.any(middle):
        idx = middle
        numer = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR))
        S_star[idx] = numer[idx] / denom_safe[idx]

    # Regions where SL >= 0
    if np.any(maskL):
        F[:, maskL] = FL[:, maskL]
    # Regions where SR <= 0
    if np.any(maskR):
        F[:, maskR] = FR[:, maskR]

    # Middle region SL < 0 < SR -> compute star states
    if np.any(middle):
        # left-star region where S_star >= 0
        left_star = middle & (S_star >= 0.0)
        if np.any(left_star):
            idx = left_star
            denom_L = (SL[idx] - S_star[idx])
            denom_L_safe = np.where(np.abs(denom_L) < EPS, EPS, denom_L)
            rhoSL = rhoL[idx] * (SL[idx] - uL[idx]) / denom_L_safe
            uSL = S_star[idx]
            # specific energy e = E/rho
            eL = EL[idx] / rhoL[idx]
            # protected inner denom
            inner_denom = rhoL[idx] * (SL[idx] - uL[idx])
            inner_denom_safe = np.where(np.abs(inner_denom) < EPS, EPS, inner_denom)
            eSL = eL + (S_star[idx] - uL[idx]) * (S_star[idx] + pL[idx] / inner_denom_safe)
            ESL = rhoSL * eSL
            UstL = np.vstack((rhoSL, rhoSL * uSL, ESL))
            F[:, idx] = FL[:, idx] + SL[idx] * (UstL - ULc[:, idx])

        # right-star region where S_star < 0
        right_star = middle & (S_star < 0.0)
        if np.any(right_star):
            idx = right_star
            denom_R = (SR[idx] - S_star[idx])
            denom_R_safe = np.where(np.abs(denom_R) < EPS, EPS, denom_R)
            rhoSR = rhoR[idx] * (SR[idx] - uR[idx]) / denom_R_safe
            uSR = S_star[idx]
            eR = ER[idx] / rhoR[idx]
            inner_denom = rhoR[idx] * (SR[idx] - uR[idx])
            inner_denom_safe = np.where(np.abs(inner_denom) < EPS, EPS, inner_denom)
            eSR = eR + (S_star[idx] - uR[idx]) * (S_star[idx] + pR[idx] / inner_denom_safe)
            ESR = rhoSR * eSR
            UstR = np.vstack((rhoSR, rhoSR * uSR, ESR))
            F[:, idx] = FR[:, idx] + SR[idx] * (UstR - URc[:, idx])

    # If any invalid numbers, compute robust HLL flux for those interfaces
    invalid = ~np.isfinite(F).all(axis=0)
    if np.any(invalid):
        # compute HLL flux for invalid interfaces
        Fhll = hll_flux(ULc, URc, FL, FR, SL, SR)
        F[:, invalid] = Fhll[:, invalid]
        print(f"[hllc_flux] Detected invalid HLLC flux at {np.sum(invalid)} interfaces; used HLL fallback.")

    return F

# --------------------------- RHS evaluator ---------------------------

def compute_rhs(U, dx, gamma, ng, visc_coeff):
    """Compute RHS dU/dt, returns array shape (3,Ntot). Only interior cells will be non-zero.
    Applies transmissive BCs internally and performs basic sanity checks.
    """
    # Apply BCs
    apply_transmissive_bc(U, ng)

    # primitives with clipping to floors
    rho, u, p = cons_to_prim(U, gamma)
    # ensure sound-speed safe
    rho = np.maximum(rho, RHO_FLOOR)
    p = np.maximum(p, P_FLOOR)

    prim = np.vstack((rho, u, p))

    # slopes and reconstructed states at interfaces
    slopes = compute_slopes(prim)
    Ntot = prim.shape[1]

    # left/right states at interfaces (i+1/2: 0..Ntot-2)
    UL = prim[:, :-1] + 0.5 * slopes[:, :-1]
    UR = prim[:, 1:] - 0.5 * slopes[:, 1:]

    # compute fluxes at interfaces robustly
    F_int = hllc_flux(UL, UR, gamma)  # shape (3, Ntot-1)

    # initialize dU
    dU = np.zeros_like(U)

    # cell indices: 0..Ntot-1 ; interfaces 0..Ntot-2 (between cells i and i+1)
    # interior cells i in [ng .. ng+N_cell-1]
    i0 = ng
    i1 = Ntot - ng  # exclusive
    # compute divergence: dU[:, i] = -(F(i+1/2) - F(i-1/2))/dx
    # F_int indexing: interface k is between cells k and k+1
    # so for cell i, left interface idx = i-1, right interface idx = i
    dU[:, i0:i1] = - (F_int[:, i0:i1] - F_int[:, i0 - 1:i1 - 1]) / dx

    # optional artificial viscosity term (Laplacian in conserved variables)
    if visc_coeff > 0.0:
        lap = np.zeros_like(U)
        lap[:, 1:-1] = U[:, 2:] - 2.0 * U[:, 1:-1] + U[:, :-2]
        dU[:, :] += visc_coeff * lap / (dx * dx)

    # sanity checks
    assert dU.shape == U.shape
    if not np.all(np.isfinite(dU)):
        print('[compute_rhs] Warning: RHS contains non-finite values; replacing with zeros at those entries.')
        bad = ~np.isfinite(dU)
        dU[bad] = 0.0

    return dU

# --------------------------- Initial condition ---------------------------

def initial_condition(xc, gamma):
    rho0 = np.where(xc <= 0.5, 1.0, 0.125)
    u0 = np.zeros_like(xc)
    p0 = np.where(xc <= 0.5, 1.0, 0.1)
    return prim_to_cons(rho0, u0, p0, gamma)

# --------------------------- Shock locator and diagnostics ---------------------------

def shock_position(xc_centers, rho, ng, N_cell):
    rho_int = rho[ng:ng + N_cell]
    x_int = xc_centers[ng:ng + N_cell]
    grad = np.abs(np.diff(rho_int))
    if len(grad) == 0:
        return 0.0, 0.0
    imax = np.argmax(grad)
    xpos = 0.5 * (x_int[imax] + x_int[imax + 1])
    return float(xpos), float(grad[imax])

# --------------------------- Main solver ---------------------------

def run_solver(N_cell, ng, gamma, cfl, output_times, visc_coeff):
    Ntot = N_cell + 2 * ng
    dx = (X1 - X0) / N_cell
    x_centers = np.linspace(X0 + 0.5 * dx, X1 - 0.5 * dx, N_cell)

    # extended centers including ghosts
    x_ext = np.zeros(Ntot)
    x_ext[ng:ng + N_cell] = x_centers
    for i in range(ng):
        x_ext[i] = x_ext[ng] - (ng - i) * dx
        x_ext[-1 - i] = x_ext[ng + N_cell - 1] + (i + 1) * dx

    # U with ghosts
    U = np.zeros((3, Ntot))
    U[:, ng:ng + N_cell] = initial_condition(x_centers, gamma)
    apply_transmissive_bc(U, ng)

    # initial conserved integrals for diagnostics
    rho0, u0, p0 = cons_to_prim(U, gamma)
    mass0 = manual_trapz(x_ext[ng:ng + N_cell], rho0[ng:ng + N_cell])
    mom0 = manual_trapz(x_ext[ng:ng + N_cell], (rho0 * u0)[ng:ng + N_cell])
    energy0 = manual_trapz(x_ext[ng:ng + N_cell], U[2, ng:ng + N_cell])

    t = 0.0
    step = 0
    out_index = 0
    outputs = []
    shock_positions = []

    # time loop
    while (out_index < len(output_times)) and (step < MAX_STEPS):
        # primitives (clipped)
        rho, u, p = cons_to_prim(U, gamma)
        rho = np.maximum(rho, RHO_FLOOR)
        p = np.maximum(p, P_FLOOR)

        # sound speeds protected
        c = np.sqrt(np.maximum(gamma * p / np.maximum(rho, EPS), 0.0))
        max_speed = np.nanmax(np.abs(u) + c)
        if not np.isfinite(max_speed) or max_speed <= 0.0:
            dt = 1e-6
        else:
            dt = cfl * (dx / max_speed)

        # do not step past next output time
        t_next_output = output_times[out_index]
        if t + dt > t_next_output:
            remaining = t_next_output - t
            if remaining <= 1e-14:
                # Remaining time effectively zero: capture snapshot without advancing
                rho_i, u_i, p_i = cons_to_prim(U, gamma)
                ri = rho_i[ng:ng + N_cell].copy()
                ui = u_i[ng:ng + N_cell].copy()
                pi = p_i[ng:ng + N_cell].copy()
                outputs.append((t, x_centers.copy(), ri, ui, pi))
                xpos, gradmax = shock_position(x_ext, rho_i, ng, N_cell)
                shock_positions.append((t, xpos, gradmax))
                print(f"[Output t={t:.12f}] Captured snapshot without advancing due to tiny remaining time.")
                out_index += 1
                continue
            dt = remaining
            # ensure dt not ridiculously small
            if dt < MIN_DT:
                dt = max(dt, MIN_DT)

        # Stage 1
        rhs = compute_rhs(U, dx, gamma, ng, visc_coeff)
        U1 = U + dt * rhs
        # positivity after stage 1
        rho1, u1, p1 = cons_to_prim(U1, gamma)
        p1 = np.maximum(p1, P_FLOOR)
        rho1 = np.maximum(rho1, RHO_FLOOR)
        U1[:, :] = prim_to_cons(rho1, u1, p1, gamma)

        # Stage 2
        rhs1 = compute_rhs(U1, dx, gamma, ng, visc_coeff)
        U_new = U + 0.5 * dt * (rhs + rhs1)

        # positivity fix
        rho_new, u_new, p_new = cons_to_prim(U_new, gamma)
        p_new = np.maximum(p_new, P_FLOOR)
        rho_new = np.maximum(rho_new, RHO_FLOOR)
        U[:, :] = prim_to_cons(rho_new, u_new, p_new, gamma)

        t += dt
        step += 1

        # capture outputs when reached or exceeded target time
        if (t >= t_next_output - 1e-12):
            rho_i, u_i, p_i = cons_to_prim(U, gamma)
            ri = rho_i[ng:ng + N_cell].copy()
            ui = u_i[ng:ng + N_cell].copy()
            pi = p_i[ng:ng + N_cell].copy()
            outputs.append((t, x_centers.copy(), ri, ui, pi))

            xpos, gradmax = shock_position(x_ext, rho_i, ng, N_cell)
            shock_positions.append((t, xpos, gradmax))

            mass = manual_trapz(x_ext[ng:ng + N_cell], ri)
            momentum = manual_trapz(x_ext[ng:ng + N_cell], ri * ui)
            energy = manual_trapz(x_ext[ng:ng + N_cell], U[2, ng:ng + N_cell])
            tv = np.sum(np.abs(np.diff(ri)))

            print(f"[Output t={t:.4f}] step={step} dt={dt:.4e}")
            print(f"  Shock pos ~ {xpos:.6f} (grad max {gradmax:.4e})")
            print(f"  Mass: {mass:.6f} (rel Δ {(mass - mass0)/mass0:.3e})")
            print(f"  Momentum: {momentum:.6f} (rel Δ {(momentum - mom0)/max(1e-12, abs(mom0)):.3e})")
            print(f"  Energy: {energy:.6f} (rel Δ {(energy - energy0)/energy0:.3e})")
            print(f"  TV(rho): {tv:.6e}")
            print(f"  rho range: [{ri.min():.4e}, {ri.max():.4e}], p range: [{pi.min():.4e}, {pi.max():.4e}]\n")

            out_index += 1

        # safety break if we reached final output
        if t >= output_times[-1] - 1e-12:
            break

    return outputs, shock_positions

# --------------------------- Plotting ---------------------------

def save_plots(outputs, out_dir='plots'):
    os.makedirs(out_dir, exist_ok=True)

    plt.figure(figsize=(8, 4))
    for t, x, rho, u, p in outputs:
        plt.plot(x, rho, label=f't={t:.2f}')
    plt.xlabel('x')
    plt.ylabel('Density')
    plt.title('Density at output times')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'density.png'))
    plt.close()

    plt.figure(figsize=(8, 4))
    for t, x, rho, u, p in outputs:
        plt.plot(x, u, label=f't={t:.2f}')
    plt.xlabel('x')
    plt.ylabel('Velocity')
    plt.title('Velocity at output times')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'velocity.png'))
    plt.close()

    plt.figure(figsize=(8, 4))
    for t, x, rho, u, p in outputs:
        plt.plot(x, p, label=f't={t:.2f}')
    plt.xlabel('x')
    plt.ylabel('Pressure')
    plt.title('Pressure at output times')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'pressure.png'))
    plt.close()

# --------------------------- Small unit test for interface flux (sanity) ---------------------------

def test_interface_flux():
    # classical Sod left/right states
    UL = np.array([1.0, 0.0, 1.0])
    UR = np.array([0.125, 0.0, 0.1])
    UL_arr = UL.reshape((3, 1))
    UR_arr = UR.reshape((3, 1))
    F_hllc = hllc_flux(UL_arr, UR_arr, GAMMA)
    # HLL for comparison
    rhoL, uL, pL = UL
    rhoR, uR, pR = UR
    EL = pL / (GAMMA - 1.0) + 0.5 * rhoL * uL * uL
    ER = pR / (GAMMA - 1.0) + 0.5 * rhoR * uR * uR
    ULc = np.vstack((rhoL, rhoL * uL, EL))
    URc = np.vstack((rhoR, rhoR * uR, ER))
    FL = np.vstack((rhoL * uL, rhoL * uL * uL + pL, (EL + pL) * uL))
    FR = np.vstack((rhoR * uR, rhoR * uR * uR + pR, (ER + pR) * uR))
    cL = np.sqrt(max(GAMMA * pL / rhoL, 0.0))
    cR = np.sqrt(max(GAMMA * pR / rhoR, 0.0))
    SL = min(uL - cL, uR - cR)
    SR = max(uL + cL, uR + cR)
    F_hll = hll_flux(ULc, URc, FL, FR, np.array([SL]), np.array([SR]))
    print('[test_interface_flux] HLLC flux:', F_hllc.ravel())
    print('[test_interface_flux] HLL   flux:', F_hll.ravel())

# --------------------------- Run and save outputs ---------------------------

if __name__ == '__main__':
    print('Running 1D Euler MUSCL-Hancock solver (robust HLLC with HLL fallback)...')
    # small unit test
    test_interface_flux()

    outputs, shock_positions = run_solver(N_CELL, NG, GAMMA, CFL, OUTPUT_TIMES, visc)
    print('Finished time integration. Saving plots...')
    save_plots(outputs)
    print('Plots saved in ./plots: density.png, velocity.png, pressure.png')
    print('\nShock positions history:')
    for t, xpos, g in shock_positions:
        print(f'  t={t:.4f} -> x_shock~{xpos:.6f} (grad {g:.3e})')

    if outputs:
        t, x, rho, u, p = outputs[-1]
        bad = np.sum((rho <= 0) | (p <= 0) | np.isnan(rho) | np.isnan(p))
        if bad == 0:
            print('\nNo negative densities or pressures detected in final snapshot.')
        else:
            print(f'\nWarning: {bad} cells with non-physical values in final snapshot.')
```



### runtime_outputs

#### Output block1

Running 1D Euler MUSCL-Hancock (MUSCL + HLLC, SSP RK2) solver...

Runtime Error:
Traceback (most recent call last):
  File "C:\Users\admin\AppData\Local\Temp\ipykernel_45176\2663542217.py", line 138, in execute_code_tool
    exec(code, global_scope, global_scope)
  File "<string>", line 397, in <module>
  File "<string>", line 289, in run_solver
  File "<string>", line 203, in compute_rhs
ValueError: operands could not be broadcast together with shapes (3,403) (3,402) (3,403) 



#### Output block2

Running 1D Euler MUSCL-Hancock (MUSCL + HLLC, SSP RK2) solver...
<string>:50: RuntimeWarning: overflow encountered in multiply
<string>:130: RuntimeWarning: overflow encountered in multiply
<string>:133: RuntimeWarning: overflow encountered in multiply
<string>:184: RuntimeWarning: overflow encountered in multiply
<string>:189: RuntimeWarning: invalid value encountered in add
<string>:202: RuntimeWarning: overflow encountered in multiply
<string>:148: RuntimeWarning: invalid value encountered in divide
<string>:189: RuntimeWarning: overflow encountered in multiply
<string>:332: RuntimeWarning: invalid value encountered in multiply
<string>:341: RuntimeWarning: invalid value encountered in multiply
[Output t=0.0011] step=5 dt=0.0000e+00
  Shock pos ~ 0.500000 (grad max 1.9757e+06)
  Mass: 9376.500094 (rel Δ 1.671e+04)
  Momentum: 104441222108388819388786981344866420195613840325150893839330546456075255711770200333922754202850820096.000000 (rel Δ 1.044e+113)
  Energy: nan (rel Δ nan)
  TV

... [Truncated: Content too long (2165 chars). Hidden middle part.] ...

 Mass: 9376.500094 (rel Δ 1.671e+04)
  Momentum: 104441222108388819388786981344866420195613840325150893839330546456075255711770200333922754202850820096.000000 (rel Δ 1.044e+113)
  Energy: nan (rel Δ nan)
  TV(rho): 7.500754e+06
  rho range: [1.0000e-08, 1.9757e+06], p range: [nan, nan]

[Output t=0.0011] step=7 dt=0.0000e+00
  Shock pos ~ 0.500000 (grad max 1.9757e+06)
  Mass: 9376.500094 (rel Δ 1.671e+04)
  Momentum: 104441222108388819388786981344866420195613840325150893839330546456075255711770200333922754202850820096.000000 (rel Δ 1.044e+113)
  Energy: nan (rel Δ nan)
  TV(rho): 7.500754e+06
  rho range: [1.0000e-08, 1.9757e+06], p range: [nan, nan]

Finished time integration. Saving plots...
Plots saved in ./plots: density.png, velocity.png, pressure.png

Shock positions history:
  t=0.0011 -> x_shock~0.500000 (grad 1.976e+06)
  t=0.0011 -> x_shock~0.500000 (grad 1.976e+06)
  t=0.0011 -> x_shock~0.500000 (grad 1.976e+06)

Warning: 4 cells with non-physical values in final snapshot.



#### Output block3

Running 1D Euler MUSCL-Hancock solver (robust HLLC with HLL fallback)...
[test_interface_flux] HLLC flux: [0.43026035 0.49090909 1.16170294]
[test_interface_flux] HLL   flux: [0.51765698 0.55       1.33111795]
[Output t=0.1000] step=174 dt=2.9081e-04
  Shock pos ~ 0.677500 (grad max 5.1541e-02)
  Mass: 0.561094 (rel Δ -6.332e-15)
  Momentum: 0.090000 (rel Δ 9.000e+10)
  Energy: 1.371563 (rel Δ -5.990e-15)
  TV(rho): 8.777411e-01
  rho range: [1.2500e-01, 1.0000e+00], p range: [1.0000e-01, 1.0000e+00]

[Output t=0.2000] step=350 dt=2.1844e-04
  Shock pos ~ 0.852500 (grad max 5.2827e-02)
  Mass: 0.561094 (rel Δ -1.029e-14)
  Momentum: 0.180000 (rel Δ 1.800e+11)
  Energy: 1.371563 (rel Δ -1.036e-14)
  TV(rho): 8.770324e-01
  rho range: [1.2500e-01, 1.0000e+00], p range: [1.0000e-01, 1.0000e+00]

[Output t=0.3000] step=526 dt=2.6189e-04
  Shock pos ~ 0.780000 (grad max 1.9211e-02)
  Mass: 0.557190 (rel Δ -6.957e-03)
  Momentum: 0.263336 (rel Δ 2.633e+11)
  Energy: 1.354632 (rel Δ -1.234e-02)
  TV(rho): 7.404972e-01
  rho range: [2.6121e-01, 1.0000e+00], p range: [2.9621e-01, 1.0000e+00]

Finished time integration. Saving plots...
Plots saved in ./plots: density.png, velocity.png, pressure.png

Shock positions history:
  t=0.1000 -> x_shock~0.677500 (grad 5.154e-02)
  t=0.2000 -> x_shock~0.852500 (grad 5.283e-02)
  t=0.3000 -> x_shock~0.780000 (grad 1.921e-02)

No negative densities or pressures detected in final snapshot.




### review_decision
revise

### review_comments
	Current Stage [A/3]
1. MAJOR

Category: runtime
Issue: Broadcast/shapes mismatch when forming flux divergence in compute_rhs: the code does dU[:, :-1] += - (F_int - np.roll(F_int,1,axis=1))[:, :-1] / dx which attempts to add an array of shape (3,402) into slice of shape (3,403) and raises a ValueError.

Recommendation: Replace the roll-based approach with explicit, shape-consistent indexing and only update the interior cells. For example, if F_int has shape (3, Ntot-1) then compute flux differences as flux_diff = F_int[:,1:] - F_int[:,:-1] (shape (3, Ntot-2)) and update dU[:,1:-1] = -flux_diff/dx. More robust: perform updates only for cells ng:ng+N_cell with dU[:, ng:ng+N_cell] = -(F_int[:, ng+1:ng+1+N_cell] - F_int[:, ng:ng+N_cell]) / dx. Avoid np.roll to prevent wraparound artifacts.


2. MAJOR

Category: accuracy
Issue: Inconsistent and incorrect definitions of energy in primitive <-> conserved conversions. prim_to_cons computes E incorrectly and multiplies by rho twice; cons_to_prim assumes the third conserved component is energy density but prim_to_cons returns a mis-scaled quantity. This produces physically wrong conserved variables and fluxes.

Recommendation: Make conserved energy the energy density (total energy per unit volume): E_total = p/(gamma-1) + 0.5 * rho * u**2, and return U = [rho, rho*u, E_total]. In prim_to_cons: return np.vstack((rho, rho*u, p/(gamma-1.0) + 0.5 * rho * u * u)). cons_to_prim should then compute u = mom/rho and p = (gamma-1)*(E_total - 0.5 * rho * u**2). Fix all places that assumed a different scaling.


3. MAJOR

Category: accuracy
Issue: HLLC implementation uses inconsistent conserved quantities and flux formulas: ULc/URc and flux third components are built with wrong factors (e.g. sometimes multiplying energy by rho again). The star-state formulas use these inconsistent variables and will produce incorrect fluxes (and possibly NaNs/negatives).

Recommendation: Re-derive and re-implement HLLC using the consistent definition that U = [rho, rho*u, E_total] where E_total is energy density. Then set flux F = [rho*u, rho*u**2 + p, (E_total + p)*u]. Compute S_L and S_R as done, S_star by the standard formula, and construct U*_L and U*_R according to Toro (1999) or standard references. If HLLC is proving tricky, temporarily replace with HLL (simpler) to get a robust working solver, then reintroduce HLLC once conserved/flux definitions are corrected. Ensure all intermediate formulas use energy density (not multiplied by rho again).


4. MINOR

Category: structure
Issue: compute_rhs docstring and behavior inconsistent: docstring says returns dU for interior cells only but function returns full dU; the code also uses a roll-based flux difference that allows wrap-around (ghost->physical). This makes boundary handling unclear and error-prone.

Recommendation: Document and implement consistent semantics: compute_rhs should return a dU array with the same shape as U and only non-zero in interior cells, or explicitly return dU_interior matching U[:, ng:ng+N_cell]. Remove roll-based arithmetic and use explicit indexing over ng..ng+N_cell to compute flux differences. Update the docstring accordingly. Also remove or use the unused 'dx_user' argument in run_solver to avoid confusion.


5. MINOR

Category: structure
Issue: Additional robustness/quality improvements: positivity enforcement, CFL and diagnostics—some are present but there are improvements left (and some unused variables like DX). Also shock locator uses diff on full array which may mis-index ghosts.

Recommendation: 1) After fixing conserved/primitive conversion and fluxes, re-run positivity checks. 2) Remove unused DX parameter or use it; tighten dt floor logic (avoid dt being set to very small positive when next output time equals current t within floating tolerances). 3) Improve shock_position to operate on interior cell values (e.g. take grad = abs(np.diff(rho[ng:ng+N_cell])) and compute interface index relative to interior grid). 4) Add unit tests: mass/momentum/energy conservation diagnostics (relative change tolerances), and compare solution to reference Sod tube (exact or high-resolution) to validate shock/contact/rarefaction positions.



	Current Stage [B/3]
1. MAJOR

Category: runtime
Issue: Simulation becomes numerically unstable almost immediately: overflow and invalid-value warnings (multiply/ divide -> inf/NaN), huge non-physical totals (mass/momentum blow-up) and NaN pressure/energy in output.

Recommendation: Root cause appears to be numerical blow-up in the Riemann/HLLC star-state energy computations and subsequent propagation of NaNs. Immediate mitigations: (a) clamp rho and p to safe floors before computing sound speeds and square roots (use rho = max(rho, RHO_FLOOR), p = max(p, P_FLOOR)); (b) avoid divides by tiny numbers in HLLC by using numpy where/masks and an EPS (e.g. 1e-12) to protect denominators; (c) simplify temporarily to the HLL solver (more robust) to confirm other parts of the code are correct, then reintroduce HLLC once star-energy formulas are verified; (d) add unit tests for a single interface and compare HLL/HLLC fluxes against known references to debug energy formulas.


2. MAJOR

Category: structure
Issue: Time-stepping logic can set dt to zero when remaining time to output is below threshold, which creates a stall/near-infinite loop and repeated identical outputs (t does not advance).

Recommendation: Do not set dt to 0.0 to 'capture' outputs. Options: (a) if remaining_time < tol, set dt = remaining_time (but ensure it's >= machine_eps or a min_dt), (b) or directly set t = t_next_output and compute a single RHS-based update or skip the time step and capture snapshot, (c) always increment step/time when dt==0 to avoid infinite loops. Replace the current dt==0 branch with a safe advance of time or an explicit output-capture without looping.


3. MAJOR

Category: accuracy
Issue: HLLC star-state energy computation and use of denominators (denom_inner) is fragile/likely incorrect, which leads to overflows. The implementation computes eSL/eSR with expressions that can be ill-conditioned and inconsistent with conservative star energy formulas.

Recommendation: Revisit the HLLC star state algebra. Use the standard, numerically-stable vectorized formula for the conservative U* and flux F* (see Toro 'Riemann Solvers' or trusted references). If you keep the current approach, (a) compute p_star robustly using the S_star expression, (b) compute U*_L and U*_R using conservative relations with protected denominators, or (c) implement HLL first (F_HLL = (SR*FL - SL*FR + SL*SR*(UR - UL))/(SR - SL)) which avoids star energy formulas and is more robust for debugging. Add extensive unit tests comparing mass/momentum/energy fluxes across a single interface with analytical expectations to validate HLLC implementation.


4. MAJOR

Category: accuracy
Issue: Physical quantities used in CFL and primitive recovery are not protected from non-physical values early enough (sound speed uses p and rho that can be zero/negative), producing invalid dt and sqrt of negative numbers.

Recommendation: Enforce positivity floors immediately after primitive recovery and before any operation that computes sound speeds, wave speeds, or time-step: in cons_to_prim or immediately after calling it, apply rho = maximum(rho, RHO_FLOOR) and p = maximum(p, P_FLOOR). Compute c = sqrt(gamma * p / rho) after clipping. Log warnings when floors are triggered to help debug where/why negativity originates.


5. MINOR

Category: structure
Issue: Flux divergence/indexing is error-prone and hard to reason about; combined with ghosts and interface indexing this may introduce off-by-one errors. Also the transmit BC assignment and ghost-centers filling should be double-checked for correctness.

Recommendation: Refactor index bookkeeping for clarity: use explicit cell-index loop or clear variable names (e.g. iface indices from 0..Ntot-2 and cell indices 0..Ntot-1) and assert shapes. Add small-N unit tests (N_cell=4, ng=2) where fluxes and dU can be checked by hand. Verify apply_transmissive_bc by asserting ghosts equal nearest interior cells after calling. Add assertions on array shapes and finite-ness (np.isfinite) in compute_rhs to catch problems earlier.



	Current Stage [C/3]
1. MAJOR

Category: accuracy
Issue: Net non-physical momentum generation and measurable loss of conservation by final time (momentum grows to ≈0.26, mass/energy drift up to ~0.7–1.2%). This violates expected symmetry of the Sod problem (initial momentum is zero) and indicates the time-stepping/positivity fixes are breaking conservation.

Recommendation: Do not rebuild conserved variables from clipped primitives after each RK stage in a way that arbitrarily changes momentum/energy. Instead: (a) apply positivity/limiters in a conservative manner, e.g. enforce positivity on reconstructed interface states before flux evaluation (Zhang & Shu style) or apply a conservative correction that modifies only the internal energy component while preserving mass and momentum; or (b) if you must clip primitives after an update, recompute the conserved energy consistently from the original momentum (E = p/(γ-1) + 0.5 * rho * (mom_old/rho)^2) so momentum remains exactly the same as in U. Add stepwise checks (per-step integrated net flux through boundaries) to locate where conservation is lost.


2. MINOR

Category: accuracy
Issue: Shock position detection is unstable / non-monotonic (shock x: 0.6775 → 0.8525 → 0.7800). The current locator simply uses max(|Δρ|) which can pick the contact/rarefaction instead of the shock and is sensitive to numerical smoothing.

Recommendation: Use a more robust detection: compute gradient of pressure (or density) and require the sign and magnitude consistent with a shock (e.g. pressure jump rather than contact); apply light smoothing or thresholding before argmax; or track the shock position in time and prefer the largest gradient that is near the previous shock position (monotonicity check). Optionally compute shock location from a threshold crossing (ρ drops below some middle value) or estimate shock speed from Rankine–Hugoniot relations between left/right averaged states.


3. MINOR

Category: structure
Issue: Reported relative change for momentum uses division by max(1e-12, abs(mom0)), producing enormous and misleading relative Δ when the initial momentum is zero (mom0≈0).

Recommendation: When baseline is nearly zero, report absolute change instead of relative. Replace denominator logic with something like denom = max(abs(m0), 1e-8) and if abs(m0) < tol print absolute Δ and note baseline≈0. This avoids huge meaningless relative numbers and makes diagnostics interpretable.


4. MAJOR

Category: accuracy
Issue: Positivity fixes are applied naively (cons_to_prim -> clip -> prim_to_cons) inside run_solver for intermediate RK stages, which can alter momentum and energy in a non-conservative way and mask the true cause of non-physical values.

Recommendation: Move positivity enforcement earlier (on reconstructed interface states) and/or use conservative positivity routines. If you must enforce positivity after a stage, modify only the total energy to ensure p >= p_floor while preserving momentum: compute u = mom/rho, set internal_energy = max(internal_energy, p_floor/(γ-1)), then set total energy = rho*internal_energy + 0.5*rho*u^2. This preserves mass and momentum exactly. Also add assertions/diagnostics to detect which cells required positivity fixes and count them per time step.


5. MINOR

Category: structure
Issue: Diagnostics and debugging aids are insufficient to trace where conservation breaks (no per-step flux balance, no count of positivity interventions, and BC-induced artifacts are not checked).

Recommendation: Add per-step diagnostics when conservation drifts are detected: compute net flux across domain boundaries (sum of F at leftmost and rightmost interfaces) and print/report if net flux ≠ 0 beyond roundoff. Also report number of cells where positivity floors were applied each snapshot and where HLL fallback triggered. Temporarily enable higher verbosity to print these quantities until conservation bug is eliminated.





### iteration_count
3

### rev_image_description
	Current Stage [A/3]


	Current Stage [B/3]


	Current Stage [C/3]




