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
We must solve the one-dimensional compressible Euler equations (mass, momentum, energy) on x in [0,1] with the classic Sod shock-tube initial condition: left state (rho=1.0,u=0.0,p=1.0) for x<=0.5 and right state (rho=0.125,u=0.0,p=0.1) for x>0.5. The system is hyperbolic and produces shocks, contact discontinuities and rarefactions. The numerical solver must be stable (no spurious oscillations), efficient, and able to produce reliable diagnostics (shock position, absence of nonphysical oscillations) and crude conservation checks (mass, momentum, energy computed by manual quadrature). Implementation constraints: always use a variable named visc (never nu) for any viscous/artificial-viscosity coefficient; pass all physical/numerical constants explicitly to functions or define them clearly in a global constants block to avoid NameError; avoid using np.trapz (implement manual integration); save plots for density, velocity, pressure at t=0.1, 0.2, 0.3 (three separate figures) and do not call plt.show(); print diagnostics. Primary mathematical/numerical challenges: (1) capturing discontinuities sharply without Gibbs oscillations (requires nonlinear limiting / nonlinear reconstruction such as MUSCL or WENO); (2) preserving positivity of density and pressure; (3) maintaining stability under the CFL restriction / choosing an appropriate time integrator (SSP-RK family or explicit Euler with predictor-corrector); (4) computing robust inter-cell fluxes (Riemann solver selection: HLLC, Roe with entropy fix, or Lax-Friedrichs); (5) ensuring global conservation up to truncation error while reporting diagnostics via manual integration; (6) balancing accuracy vs computational cost (high-order reconstructions give better resolution but cost more and require stronger stabilization).

### solution_plans
	Current Stage [A/2]
solu_name='Finite-Volume Godunov-type solver: MUSCL + HLLC + SSP-RK3' content='Governing idea:\n- Use a conservative finite-volume discretization on a uniform cell-centered grid. Reconstruct primitive variables with a slope-limited linear (MUSCL) reconstruction to obtain left/right states at interfaces. Use a robust approximate Riemann solver (HLLC) to compute interface fluxes. Advance in time using a strong-stability-preserving third-order Runge–Kutta (SSP-RK3). Maintain positivity of density and pressure by applying small floors after each stage. Use a user-specified visc (artificial-viscosity coefficient) only if extra smoothing is needed (implemented as a second-order dissipative term or via a flux-blend with Lax-Friedrichs) and pass visc and all other constants explicitly.\n\nAlgorithmic steps (numbered):\n1) Grid and initialization\n   a) Choose N cells, cell-centers x_i and cell width dx = 1/N. Initialize conservative variables U_i = [rho, rho*u, rho*E] from given primitive initial condition.\n   b) Define constants object: gamma, CFL, visc, t_snapshots = [0.1,0.2,0.3], final time = max(t_snapshots). Ensure these constants are passed to all functions.\n\n2) Reconstruction (MUSCL)\n   a) Convert conservative U to primitives W = [rho,u,p].\n   b) Compute limited slopes for each primitive component using a chosen limiter (minmod or van Leer). Slope choice affects oscillations vs diffusion.\n   c) Reconstruct left and right states at each interface: W_{i+1/2}^L = W_i + 0.5*slope_i, W_{i+1/2}^R = W_{i+1} - 0.5*slope_{i+1}.\n\n3) Riemann solver (HLLC) at each interface\n   a) Given left/right primitive states, compute sound speeds and estimate wave speeds (S_L, S_R) using standard estimates (e.g., Einfeldt or Davis estimates). Compute middle-state flux via HLLC formula to preserve contact wave.\n   b) Optionally, if visc>0, augment the interface flux with a small viscous/diffusive flux: F_total = F_HLLC - visc * (U^R - U^L)/dx, where visc is small and explicitly passed.\n\n4) Boundary conditions\n   a) Use transmissive/outflow (zero-gradient) boundary conditions implemented by ghost cells: copy nearest interior state to ghost cells.\n\n5) Time integration (SSP-RK3)\n   a) Compute time step dt = CFL * dx / max(|u|+c) across domain, with dt also bounded by remaining time to next snapshot.\n   b) SSP-RK3 stages: U^{(1)} = U^n + dt*L(U^n); U^{(2)} = 3/4 U^n + 1/4 (U^{(1)} + dt*L(U^{(1)})); U^{n+1} = 1/3 U^n + 2/3 (U^{(2)} + dt*L(U^{(2)})). Here L(U) = - (F_{i+1/2} - F_{i-1/2})/dx.\n   c) After each substage, enforce physical floors: rho = max(rho, rho_min), p = max(p, p_min), and recompute conserved variables accordingly.\n\n6) Diagnostics (print/update)\n   a) At each saved snapshot time (0.1, 0.2, 0.3): compute primitives, save figure data (do not call plt.show()).\n   b) Global conserved quantities: compute mass, momentum, energy via manual trapezoidal integration: sum over cells of cell-averages times dx, or implement manual trapezoid: integral = sum_{i=0}^{N-1} 0.5*(Q_i+Q_{i+1})*dx where Q_{N} is boundary copy—practically, for cell-centered constant cell-average Q_i, total = sum_i Q_i*dx.\n   c) Shock position detection: compute discrete gradient of density (|rho_{i+1}-rho_i|/dx) and locate x where this gradient is maximal; print x_shock and its time evolution.\n   d) Print max/min of rho and p to monitor negativity and spurious oscillations.\n\n7) Output and plotting\n   a) For each snapshot, save three figures (density, velocity, pressure) individually to files with descriptive names. Do not call plt.show().\n\nStability / accuracy / complexity / efficiency limitations:\n- Stability: explicit RK + finite-volume requires CFL <= ~0.5 (SSP-RK3) computed from max(|u|+c); if visc is used as explicit diffusion it further constrains dt via visc*dt/dx^2 for parabolic stability if visc introduces second-order diffusion; avoid large visc.\n- Accuracy: MUSCL with slope limiters is second-order in smooth regions but reverts to first order at discontinuities. Choice of limiter trades sharpness (compressive limiters) vs robustness (minmod).\n- Shock capturing: HLLC resolves contacts better than HLL/Lax-Friedrichs but may require an entropy fix in rare Roe-like variants; HLLC is robust for Sod problem and widely used.\n- Complexity: per time step cost is O(N) for reconstruction + fluxes; memory O(N). SSP-RK3 requires three flux evaluations per step.\n- Efficiency: small memory footprint and straightforward; however high-order methods (WENO/DG) can yield sharper profiles for larger cost. Using visc as artificial viscosity increases diffusion and can smooth shocks excessively if set too large. Ensure all constants including visc are provided to every function to prevent NameError.\n'

	Current Stage [B/2]
solu_name='High-resolution 5th-order WENO + HLLC flux with SSP-RK3 and optional adaptive viscous smoothing' content='Governing idea:\n- Use a high-order, shock-capturing finite-difference/finite-volume hybrid: reconstruct left/right states at interfaces with 5th-order WENO to obtain non-oscillatory high-order accuracy in smooth regions and stable behavior at discontinuities. Use HLLC approximate Riemann solver for numerical fluxes. Integrate in time with SSP-RK3. Optionally include a carefully designed, coefficient-controlled artificial viscosity term (named visc everywhere) that is activated only near detected discontinuities to avoid smearing smooth regions. Pass gamma, CFL, visc and any thresholds explicitly to every function.\n\nAlgorithmic steps (numbered):\n1) Grid and initialization\n   a) Define N, dx, cell-centered or staggered formulation; initialize primitive variables from Sod initial condition and compute conservative variables.\n   b) Define constants: gamma, CFL, visc_base, visc_factor, discontinuity_threshold, t_snapshots.\n\n2) High-order reconstruction (WENO5)\n   a) For each primitive variable (rho,u,p) compute left/right interface states using classical finite-difference WENO-5 stencil and smoothness indicators.\n   b) Use local characteristic projection optionally for better nonlinear stability (project primitive vectors onto characteristic fields, perform WENO there, reconstruct back). This is more expensive but reduces spurious oscillations across strong waves.\n\n3) Riemann solver (HLLC)\n   a) Compute fluxes at each interface from reconstructed states. If characteristic WENO used, convert reconstructed states back to primitives before Riemann solve.\n\n4) Adaptive artificial viscosity (optional but recommended for very coarse grids)\n   a) Compute a shock sensor S_i (e.g., scaled absolute second derivative of pressure or normalized density jump). Where S_i > discontinuity_threshold, set local visc_i = visc_base * (S_i) * visc_factor; else visc_i = 0.\n   b) Add a dissipative flux term proportional to visc_i * (U^R - U^L)/dx at the interface. Use variable visc per-interface but always refer to it as visc.\n\n5) Time stepping (SSP-RK3)\n   a) Compute dt from CFL and global maximum wave speed. Use RK stages as in plan A.\n   b) After each stage enforce positivity floors on rho and p, recompute U.\n\n6) Boundary conditions\n   a) Use ghost cells consistent with WENO stencil width (two ghost cells on each side). Apply transmissive/outflow (copy interior to ghost cells) or physical reflective if required.\n\n7) Diagnostics\n   a) Manual integration for global mass/momentum/energy: compute total = sum_i U_i * dx (explicit sum without np.trapz). For cell-centered constant U_i, total = sum_i U_i*dx.\n   b) Detect shock/contact positions via steep gradient or zero crossing in Riemann invariants; find index of maximum |rho_{i+1}-rho_i| and map to x.\n   c) Compute and print norms that indicate spurious oscillations: max overshoot/undershoot of density and pressure compared to expected physical bounds, L1 difference to a reference solution if available.\n\n8) Output and plotting\n   a) At each snapshot time save three figures (density, velocity, pressure) separately, using filenames that include time. Do not call plt.show().\n\nStability / accuracy / complexity / efficiency limitations:\n- Stability: WENO5 with SSP-RK3 is stable under the usual hyperbolic CFL restriction (CFL~0.3-0.6 depending on implementation). Adding visc as an explicit diffusion term may require a parabolic stability constraint dt <= O(dx^2/visc) if visc is large or implemented as second-derivative diffusion, so keep visc small and activate it adaptively.\n- Accuracy: WENO5 is high-order (5th) in smooth regions and avoids Gibbs near discontinuities; however, near discontinuities accuracy falls to first order. Characteristic projection improves robustness but increases cost.\n- Complexity: per time step cost is higher than MUSCL: WENO requires wider stencils and more arithmetic (O(5N) operations with larger constant), and if characteristic decomposition used, additional matrix transforms per interface increase cost by ~O(N*states^2).\n- Efficiency: yields sharper resolution on a given grid than MUSCL but requires more CPU; adaptive visc helps avoid oscillations on coarse grids but requires careful tuning of thresholds; WENO implementations are more complex and more sensitive to boundary ghost-cell handling.\n- Implementation caution: ensure all constants (gamma, CFL, visc_base, visc_factor, discontinuity_threshold) are passed explicitly to reconstruction, flux, and time-stepping functions to avoid NameError; implement manual integration via explicit summation for conservation diagnostics; and always use the variable name visc (never nu) for any viscosity-related parameter.\n'



### technical_spec
	Current Stage [A/3]
Finite-volume MUSCL solver using cell-centered grid with ghost cells. Data flow: initialize conserved variables U -> apply transmissive BCs -> convert U to primitives W -> compute limited slopes (minmod) -> reconstruct left/right states at interfaces -> compute interface fluxes via HLLC (optionally add viscous diffusive flux using visc parameter) -> assemble RHS = -(F_{i+1/2}-F_{i-1/2})/dx -> advance in time with SSP-RK3 (three stages) enforcing positivity floors after each stage. Diagnostics: manual integration of conserved quantities by summation, shock position via max density gradient, min/max checks for rho and p. Outputs: three PNG figures (density, velocity, pressure) each showing curves at t=0.1,0.2,0.3. All physical and numerical constants are explicitly defined in CONSTS and passed into routines.

	Current Stage [B/3]
This solver implements a finite-volume MUSCL-HLLC Euler 1D Sod problem with SSP-RK3 time integration. Key design choices:
- Data layout: conserved U shape (3, N+2) with two ghost cells. Cell centers x_cell length N.
- Reconstruction: minmod slope limiter (cell-centered) used to build left/right states at N+1 interfaces.
- Riemann solver: HLLC implemented following Toro (1997) formulas. The function accepts all relevant constants and does not read globals.
- Time stepping: SSP-RK3 with positivity enforcement after each stage.
- Efficiency: workspace arrays (W, slopes, WL_if, WR_if, F_if) are preallocated once in solve_euler and reused each call to compute_RHS to avoid reallocation in the time loop.
- Diagnostics: conserved totals, boundary fluxes, max CFL diagnostics printed for debugging; improved shock detection combines pressure and density gradients.
- Tests: unit checks ensure identical left/right states produce flux equal to physical flux (no artifact).
- Outputs: three PNG figures (density, velocity, pressure) at t=0.1,0.2,0.3.
All constants are passed explicitly via CONSTS and into functions; no hidden global dependencies remain.

	Current Stage [C/3]
Finite-volume MUSCL-HLLC 1D Euler solver with SSP-RK3. Dataflow:
- initial_condition builds conserved U with 2 ghost cells and cell centers x_cell.
- Preallocated workspace arrays (W, slopes, WL_if, WR_if, F_if) are created once in solve_euler and reused.
- compute_RHS: applies BCs, converts U->primitive W, computes limited slopes, reconstructs WL/WR at interfaces, calls hllc_flux (robust denom regularization), computes interface fluxes and RHS. Emits flux-diff diagnostics when debug.
- Time integration: SSP-RK3 with positivity enforcement optionally applied after each stage; energy change due to enforcement reported when enabled.
- Diagnostics: unit tests (identical-state and Sod-like), runtime NaN/Inf checks, conservation residual per time step computed as totals_new - totals_old + dt*(F_right - F_left), shock locator uses max |d rho/dx|. Outputs saved as PNGs.



### python_codes

#### Script block1:

```python
#!/usr/bin/env python3
"""
MUSCL finite-volume solver with HLLC Riemann solver and SSP-RK3 time stepping
for the 1D Euler equations (Sod shock tube).
Saves three figures: density, velocity, pressure at t=0.1,0.2,0.3
All physical/numerical constants are defined explicitly in CONSTS.
"""
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

# -------------------- Global constants (passed/used explicitly) --------------------
CONSTS = {
    'gamma': 1.4,
    'CFL': 0.45,
    'visc': 0.0,            # artificial viscosity coefficient (use small >0 if needed)
    'rho_min': 1e-8,
    'p_min': 1e-10,
    'N': 400,               # number of cells
    'x0': 0.0,
    'x1': 1.0,
    't_snapshots': [0.1, 0.2, 0.3],
}

# -------------------- Utility conversions --------------------

def prim_to_cons(rho, u, p, gamma):
    E = 0.5 * u * u + p / ((gamma - 1.0) * rho)
    return np.array([rho, rho * u, rho * E])


def cons_to_prim(U, gamma):
    rho = U[0]
    u = U[1] / rho
    E = U[2] / rho
    p = (gamma - 1.0) * rho * (E - 0.5 * u * u)
    return rho, u, p

# -------------------- Minmod slope limiter --------------------

def minmod(a, b):
    res = np.zeros_like(a)
    mask = (a * b) > 0.0
    res[mask] = np.where(np.abs(a[mask]) < np.abs(b[mask]), a[mask], b[mask])
    return res


def compute_slopes(W):
    # W shape: (3, N+2) including ghost cells
    # slopes for interior cells (1..N) using centered differences with minmod limiter
    dwf = W[:, 1:-1] - W[:, :-2]    # forward-backward combined (W_{i+1} - W_{i-1})
    dfr = W[:, 2:] - W[:, 1:-1]     # forward (W_{i+1} - W_i)
    dlf = W[:, 1:-1] - W[:, :-2]    # backward (W_i - W_{i-1})
    # minmod of forward and backward slopes
    slope = minmod(dfr, dlf)
    # pad to full length N+2 with zeros at ghosts
    slopes = np.zeros_like(W)
    slopes[:, 1:-1] = slope
    return slopes

# -------------------- HLLC Riemann solver --------------------

def hllc_flux(WL, WR, gamma, visc, dx):
    # WL, WR: primitive states arrays [rho,u,p] for each interface (shape (3, M))
    # returns flux array shape (3, M)
    rhoL, uL, pL = WL[0], WL[1], WL[2]
    rhoR, uR, pR = WR[0], WR[1], WR[2]

    # ensure positivity for sound speeds
    rhoL = np.maximum(rhoL, CONSTS['rho_min'])
    rhoR = np.maximum(rhoR, CONSTS['rho_min'])
    pL = np.maximum(pL, CONSTS['p_min'])
    pR = np.maximum(pR, CONSTS['p_min'])

    cL = np.sqrt(gamma * pL / rhoL)
    cR = np.sqrt(gamma * pR / rhoR)

    # Conserved variables
    UL = np.vstack([rhoL, rhoL * uL, rhoL * (0.5 * uL * uL + pL / ((gamma - 1.0) * rhoL))])
    UR = np.vstack([rhoR, rhoR * uR, rhoR * (0.5 * uR * uR + pR / ((gamma - 1.0) * rhoR))])

    # Physical fluxes
    FL = np.vstack([rhoL * uL,
                    rhoL * uL * uL + pL,
                    (UL[2] + pL) * uL])
    FR = np.vstack([rhoR * uR,
                    rhoR * uR * uR + pR,
                    (UR[2] + pR) * uR])

    # Wave speed estimates (Davis/Einfeldt type)
    S_L = np.minimum(uL - cL, uR - cR)
    S_R = np.maximum(uL + cL, uR + cR)

    # Avoid pathological cases
    S_L = np.minimum(S_L, np.minimum(uL, uR) - 1e-6)
    S_R = np.maximum(S_R, np.maximum(uL, uR) + 1e-6)

    # Contact speed S_M from Rankine-Hugoniot
    num = pR - pL + rhoL * uL * (S_L - uL) - rhoR * uR * (S_R - uR)
    den = rhoL * (S_L - uL) - rhoR * (S_R - uR)
    # protect division
    den = np.where(np.abs(den) < 1e-12, 1e-12 * np.sign(den + 1e-12), den)
    S_M = num / den

    # Star region states
    # U*_L and U*_R
    rho_star_L = rhoL * (S_L - uL) / (S_L - S_M)
    rho_star_R = rhoR * (S_R - uR) / (S_R - S_M)

    E_L = UL[2] / rhoL
    E_R = UR[2] / rhoR

    # compute specific total energy in star regions
    E_star_L = E_L + (S_M - uL) * (S_M + pL / (rhoL * (S_L - uL)))
    E_star_R = E_R + (S_M - uR) * (S_M + pR / (rhoR * (S_R - uR)))

    U_star_L = np.vstack([rho_star_L,
                          rho_star_L * S_M,
                          rho_star_L * E_star_L])
    U_star_R = np.vstack([rho_star_R,
                          rho_star_R * S_M,
                          rho_star_R * E_star_R])

    # HLLC flux assembly
    # Initialize flux array
    M = FL.shape[1]
    F = np.zeros_like(FL)

    # Regions
    # S_R <= 0 -> use FR
    mask_FR = S_R <= 0
    F[:, mask_FR] = FR[:, mask_FR]

    # S_L >= 0 -> use FL
    mask_FL = S_L >= 0
    F[:, mask_FL] = FL[:, mask_FL]

    # S_L < 0 < S_M -> FL + S_L*(U*_L - UL)
    mask_L = (S_L < 0) & (S_M >= 0)
    F[:, mask_L] = FL[:, mask_L] + S_L[mask_L] * (U_star_L[:, mask_L] - UL[:, mask_L])

    # S_M < 0 < S_R -> FR + S_R*(U*_R - UR)
    mask_R = (S_M < 0) & (S_R > 0)
    F[:, mask_R] = FR[:, mask_R] + S_R[mask_R] * (U_star_R[:, mask_R] - UR[:, mask_R])

    # Add small artificial viscous/diffusive flux if requested
    if visc is not None and visc > 0.0:
        diff = (UR - UL) / dx
        F -= visc * diff

    return F

# -------------------- RHS (spatial operator) --------------------

def compute_RHS(U, dx, gamma, visc):
    # U shape: (3, N+2) including ghost cells; compute RHS for interior cells 1..N
    # Apply boundary conditions (transmissive/outflow)
    U[:, 0] = U[:, 1]
    U[:, -1] = U[:, -2]

    # Convert to primitive W arrays
    # W shape (3, N+2)
    W = np.zeros_like(U)
    rho = U[0, :]
    u = U[1, :] / rho
    E = U[2, :] / rho
    p = (gamma - 1.0) * rho * (E - 0.5 * u * u)
    W[0, :] = rho
    W[1, :] = u
    W[2, :] = p

    # compute slopes
    slopes = compute_slopes(W)

    # reconstruct interface states for interfaces i=0..N (between cells i and i+1)
    # left at interface i is cell i, right is cell i+1
    # W_L_interface[:, i] = W[:, i] + 0.5*slopes[:, i]
    WL = W + 0.5 * slopes
    WR = np.zeros_like(W)
    WR[:, :-1] = W[:, 1:] - 0.5 * slopes[:, 1:]
    # Now WL[:,0..N] corresponds left states for interfaces 0..N, WR same

    # Compute fluxes at interfaces
    # Pass arrays of primitives for interfaces (shape (3, N+1))
    WL_if = WL[:, :-1]
    WR_if = WR[:, :-1]
    F_if = hllc_flux(WL_if, WR_if, gamma, visc, dx)

    # Compute RHS for interior cells 1..N
    RHS = np.zeros_like(U)
    # For cell i (1..N), flux difference F_{i+1/2} - F_{i-1/2}
    # F_if indices: interface k between cell k and k+1, so cell i has interfaces i-1 and i
    RHS[:, 1:-1] = -(F_if[:, 1:] - F_if[:, :-1]) / dx
    return RHS

# -------------------- Initial condition --------------------

def initial_condition(N, x0, x1):
    dx = (x1 - x0) / N
    x_cell = (np.arange(N) + 0.5) * dx + x0
    # allocate U with ghost cells: size N+2
    U = np.zeros((3, N + 2))
    # set primitives
    rho = np.where(x_cell <= 0.5, 1.0, 0.125)
    u = np.zeros_like(rho)
    p = np.where(x_cell <= 0.5, 1.0, 0.1)
    # fill interior cells 1..N
    for i in range(N):
        U[:, i + 1] = prim_to_cons(rho[i], u[i], p[i], CONSTS['gamma'])
    # ghost cells set by BC routine later
    return U, x_cell, dx

# -------------------- Diagnostics and manual integration --------------------

def compute_conserved_totals(U, dx):
    # interior cells only
    Uint = U[:, 1:-1]
    mass = np.sum(Uint[0, :] * dx)
    momentum = np.sum(Uint[1, :] * dx)
    energy = np.sum(Uint[2, :] * dx)
    return mass, momentum, energy


def primitives_from_U(U, gamma):
    rho = U[0, 1:-1]
    u = U[1, 1:-1] / rho
    E = U[2, 1:-1] / rho
    p = (gamma - 1.0) * rho * (E - 0.5 * u * u)
    return rho, u, p

# -------------------- Main solver --------------------

def solve_euler(CONSTS):
    gamma = CONSTS['gamma']
    CFL = CONSTS['CFL']
    visc = CONSTS['visc']
    N = CONSTS['N']
    x0 = CONSTS['x0']
    x1 = CONSTS['x1']
    snapshots = list(CONSTS['t_snapshots'])
    t_final = max(snapshots)

    U, x_cell, dx = initial_condition(N, x0, x1)

    # pre-allocate arrays reused in loop for efficiency
    RHS = np.zeros_like(U)

    t = 0.0
    # results storage for snapshots
    results = {}

    # Initial global conserved totals
    mass0, mom0, E0 = compute_conserved_totals(U, dx)
    print(f"Initial totals -> mass: {mass0:.6f}, momentum: {mom0:.6e}, energy: {E0:.6f}")

    # time loop
    while t < t_final - 1e-12:
        # compute primitive states to get max wave speed
        rho, u, p = primitives_from_U(U, gamma)
        c = np.sqrt(np.maximum(1e-12, gamma * p / rho))
        max_speed = np.max(np.abs(u) + c)
        dt = CFL * dx / max_speed
        # limit dt to not overshoot next snapshot or final time
        next_snapshot = min([s for s in snapshots if s > t] + [t_final])
        if t + dt > next_snapshot:
            dt = next_snapshot - t
        if dt <= 0.0:
            break

        # SSP-RK3 stages
        U0 = U.copy()

        # Stage 1
        RHS = compute_RHS(U, dx, gamma, visc)
        U1 = U0 + dt * RHS
        # enforce floors
        enforce_physical(U1, CONSTS)

        # Stage 2
        RHS = compute_RHS(U1, dx, gamma, visc)
        U2 = 0.75 * U0 + 0.25 * (U1 + dt * RHS)
        enforce_physical(U2, CONSTS)

        # Stage 3
        RHS = compute_RHS(U2, dx, gamma, visc)
        U_new = (1.0 / 3.0) * U0 + (2.0 / 3.0) * (U2 + dt * RHS)
        enforce_physical(U_new, CONSTS)

        U = U_new
        t += dt

        # If we reached a snapshot time (within tolerance), store result and print diagnostics
        if any(abs(t - s) < 1e-10 or t >= s and s not in results for s in snapshots):
            # store for any snapshot times <= t not yet stored
            for s in [s for s in snapshots if s <= t and s not in results]:
                rho_s, u_s, p_s = primitives_from_U(U, gamma)
                results[s] = (rho_s.copy(), u_s.copy(), p_s.copy(), x_cell.copy())
                # diagnostics
                mass, mom, E = compute_conserved_totals(U, dx)
                # shock position via max density gradient
                grad = np.abs(np.diff(rho_s)) / dx
                j = np.argmax(grad)
                x_shock = 0.5 * (x_cell[j] + x_cell[j + 1])
                print(f"t={s:.3f}: shock_position ~ {x_shock:.5f}, rho_min={rho_s.min():.6e}, rho_max={rho_s.max():.6f}, p_min={p_s.min():.6e}, p_max={p_s.max():.6f}")
                print(f"    Totals -> mass: {mass:.6f} (init {mass0:.6f}), momentum: {mom:.6e}, energy: {E:.6f}")

    return results


def enforce_physical(U, CONSTS):
    # enforce floors on density and pressure, updating conserved variables accordingly
    gamma = CONSTS['gamma']
    rho_min = CONSTS['rho_min']
    p_min = CONSTS['p_min']
    # interior cells
    for i in range(1, U.shape[1] - 1):
        rho = U[0, i]
        mom = U[1, i]
        E = U[2, i]
        rho = max(rho, rho_min)
        u = mom / rho
        total_E = E / rho
        p = (gamma - 1.0) * rho * (total_E - 0.5 * u * u)
        if p < p_min:
            p = p_min
            # reconstruct energy from p
            total_E = p / ((gamma - 1.0) * rho) + 0.5 * u * u
        U[0, i] = rho
        U[1, i] = rho * u
        U[2, i] = rho * total_E

# -------------------- Plotting --------------------

def save_plots(results, CONSTS):
    times = sorted(results.keys())
    xcell = results[times[0]][3]

    # Density plot
    plt.figure(figsize=(8, 4))
    for t in times:
        rho, u, p, x = results[t]
        plt.plot(x, rho, label=f"t={t:.2f}")
    plt.xlabel('x')
    plt.ylabel('Density')
    plt.title('Density at snapshots')
    plt.legend()
    plt.grid(True)
    plt.savefig('density_snapshots.png', dpi=150)
    plt.close()

    # Velocity plot
    plt.figure(figsize=(8, 4))
    for t in times:
        rho, u, p, x = results[t]
        plt.plot(x, u, label=f"t={t:.2f}")
    plt.xlabel('x')
    plt.ylabel('Velocity')
    plt.title('Velocity at snapshots')
    plt.legend()
    plt.grid(True)
    plt.savefig('velocity_snapshots.png', dpi=150)
    plt.close()

    # Pressure plot
    plt.figure(figsize=(8, 4))
    for t in times:
        rho, u, p, x = results[t]
        plt.plot(x, p, label=f"t={t:.2f}")
    plt.xlabel('x')
    plt.ylabel('Pressure')
    plt.title('Pressure at snapshots')
    plt.legend()
    plt.grid(True)
    plt.savefig('pressure_snapshots.png', dpi=150)
    plt.close()

# -------------------- Run solver and save outputs --------------------

def main():
    results = solve_euler(CONSTS)
    save_plots(results, CONSTS)
    print('Saved figures: density_snapshots.png, velocity_snapshots.png, pressure_snapshots.png')

if __name__ == '__main__':
    main()
```


#### Script block2:

```python
#!/usr/bin/env python3
"""
MUSCL finite-volume solver with HLLC Riemann solver and SSP-RK3 time stepping
for the 1D Euler equations (Sod shock tube).
Saves three figures: density, velocity, pressure at t=0.1,0.2,0.3
All physical/numerical constants are defined explicitly in CONSTS and passed into routines.
"""
import numpy as np
import matplotlib.pyplot as plt

# -------------------- Global constants (passed/used explicitly) --------------------
CONSTS = {
    'gamma': 1.4,
    'CFL': 0.45,
    'visc': 0.0,            # artificial viscosity coefficient
    'rho_min': 1e-8,
    'p_min': 1e-10,
    'N': 400,               # number of cells
    'x0': 0.0,
    'x1': 1.0,
    't_snapshots': [0.1, 0.2, 0.3],
}

# -------------------- Utility conversions --------------------

def prim_to_cons(rho, u, p, gamma):
    E = 0.5 * u * u + p / ((gamma - 1.0) * rho)
    return np.array([rho, rho * u, rho * E])


def cons_to_prim(U, gamma):
    rho = U[0]
    u = U[1] / rho
    E = U[2] / rho
    p = (gamma - 1.0) * rho * (E - 0.5 * u * u)
    return rho, u, p

# -------------------- Minmod slope limiter --------------------

def minmod(a, b):
    res = np.zeros_like(a)
    mask = (a * b) > 0.0
    res[mask] = np.where(np.abs(a[mask]) < np.abs(b[mask]), a[mask], b[mask])
    return res


def compute_slopes(W, slopes):
    # W shape: (3, N+2); slopes is preallocated same shape
    dfr = W[:, 2:] - W[:, 1:-1]
    dlf = W[:, 1:-1] - W[:, :-2]
    slope = minmod(dfr, dlf)
    slopes[:, :] = 0.0
    slopes[:, 1:-1] = slope
    return slopes

# -------------------- Apply boundary conditions --------------------

def apply_bc(U):
    # transmissive/outflow BCs: copy adjacent interior cell into ghosts
    U[:, 0] = U[:, 1]
    U[:, -1] = U[:, -2]

# -------------------- HLLC Riemann solver (vectorized) --------------------

def hllc_flux(WL, WR, gamma, visc, rho_min, p_min, dx):
    # WL, WR: primitive states arrays [3, M] where rows are rho,u,p
    rhoL = np.maximum(WL[0, :], rho_min)
    uL = WL[1, :]
    pL = np.maximum(WL[2, :], p_min)

    rhoR = np.maximum(WR[0, :], rho_min)
    uR = WR[1, :]
    pR = np.maximum(WR[2, :], p_min)

    cL = np.sqrt(gamma * pL / rhoL)
    cR = np.sqrt(gamma * pR / rhoR)

    UL = np.vstack([rhoL, rhoL * uL, rhoL * (0.5 * uL * uL + pL / ((gamma - 1.0) * rhoL))])
    UR = np.vstack([rhoR, rhoR * uR, rhoR * (0.5 * uR * uR + pR / ((gamma - 1.0) * rhoR))])

    FL = np.vstack([rhoL * uL,
                    rhoL * uL * uL + pL,
                    (UL[2, :] + pL) * uL])
    FR = np.vstack([rhoR * uR,
                    rhoR * uR * uR + pR,
                    (UR[2, :] + pR) * uR])

    # Wave speed estimates
    S_L = np.minimum(uL - cL, uR - cR)
    S_R = np.maximum(uL + cL, uR + cR)

    # Safeguards
    S_L = np.minimum(S_L, np.minimum(uL, uR) - 1e-12)
    S_R = np.maximum(S_R, np.maximum(uL, uR) + 1e-12)

    num = pR - pL + rhoL * uL * (S_L - uL) - rhoR * uR * (S_R - uR)
    den = rhoL * (S_L - uL) - rhoR * (S_R - uR)
    small = 1e-12
    den_safe = np.where(np.abs(den) < small, np.sign(den) * small + small, den)
    S_M = num / den_safe

    # star densities
    rho_star_L = rhoL * (S_L - uL) / (S_L - S_M)
    rho_star_R = rhoR * (S_R - uR) / (S_R - S_M)

    # specific total energies
    E_L = UL[2, :] / rhoL
    E_R = UR[2, :] / rhoR

    E_star_L = E_L + (S_M - uL) * (S_M + pL / (rhoL * (S_L - uL)))
    E_star_R = E_R + (S_M - uR) * (S_M + pR / (rhoR * (S_R - uR)))

    U_star_L = np.vstack([rho_star_L,
                          rho_star_L * S_M,
                          rho_star_L * E_star_L])
    U_star_R = np.vstack([rho_star_R,
                          rho_star_R * S_M,
                          rho_star_R * E_star_R])

    M = WL.shape[1]
    F = np.zeros_like(FL)

    mask_FR = S_R <= 0.0
    mask_FL = S_L >= 0.0
    mask_L = (S_L < 0.0) & (S_M >= 0.0)
    mask_R = (S_M < 0.0) & (S_R > 0.0)

    if np.any(mask_FL):
        F[:, mask_FL] = FL[:, mask_FL]
    if np.any(mask_FR):
        F[:, mask_FR] = FR[:, mask_FR]
    if np.any(mask_L):
        F[:, mask_L] = FL[:, mask_L] + S_L[mask_L] * (U_star_L[:, mask_L] - UL[:, mask_L])
    if np.any(mask_R):
        F[:, mask_R] = FR[:, mask_R] + S_R[mask_R] * (U_star_R[:, mask_R] - UR[:, mask_R])

    # optional viscous diffusive flux
    if visc is not None and visc > 0.0:
        diff = (UR - UL) / dx
        F = F - visc * diff

    return F

# -------------------- RHS (spatial operator) --------------------

def compute_RHS(U, dx, gamma, visc, rho_min, p_min, W, slopes, WL_if, WR_if, F_if):
    # U shape: (3, N+2). This routine will NOT mutate U externally.
    Ur = U.copy()
    apply_bc(Ur)

    # primitives into provided workspace W
    rho = Ur[0, :]
    u = Ur[1, :] / rho
    E = Ur[2, :] / rho
    p = (gamma - 1.0) * rho * (E - 0.5 * u * u)
    W[0, :] = rho
    W[1, :] = u
    W[2, :] = p

    # slopes (in-place)
    compute_slopes(W, slopes)

    # reconstruct interface states into WL_if, WR_if (shape (3, N+1))
    WL_if[:, :] = W[:, :-1] + 0.5 * slopes[:, :-1]
    WR_if[:, :] = W[:, 1:] - 0.5 * slopes[:, 1:]

    # fluxes at interfaces
    F = hllc_flux(WL_if, WR_if, gamma, visc, rho_min, p_min, dx)
    F_if[:, :] = F

    # Compute RHS for interior cells 1..N
    RHS = np.zeros_like(U)
    RHS[:, 1:-1] = -(F_if[:, 1:] - F_if[:, :-1]) / dx
    return RHS, F_if[:, 0], F_if[:, -1]

# -------------------- Initial condition --------------------

def initial_condition(N, x0, x1, gamma):
    dx = (x1 - x0) / N
    x_cell = (np.arange(N) + 0.5) * dx + x0
    U = np.zeros((3, N + 2))
    rho = np.where(x_cell <= 0.5, 1.0, 0.125)
    u = np.zeros_like(rho)
    p = np.where(x_cell <= 0.5, 1.0, 0.1)
    for i in range(N):
        U[:, i + 1] = prim_to_cons(rho[i], u[i], p[i], gamma)
    return U, x_cell, dx

# -------------------- Diagnostics and manual integration --------------------

def compute_conserved_totals(U, dx):
    Uint = U[:, 1:-1]
    mass = np.sum(Uint[0, :] * dx)
    momentum = np.sum(Uint[1, :] * dx)
    energy = np.sum(Uint[2, :] * dx)
    return mass, momentum, energy


def primitives_from_U(U, gamma):
    rho = U[0, 1:-1]
    u = U[1, 1:-1] / rho
    E = U[2, 1:-1] / rho
    p = (gamma - 1.0) * rho * (E - 0.5 * u * u)
    return rho, u, p

# -------------------- Positivity enforcement --------------------

def enforce_physical(U, CONSTS):
    gamma = CONSTS['gamma']
    rho_min = CONSTS['rho_min']
    p_min = CONSTS['p_min']
    for i in range(1, U.shape[1] - 1):
        rho = max(U[0, i], rho_min)
        mom = U[1, i]
        E_tot = U[2, i]
        u = mom / rho
        total_E = E_tot / rho
        p = (gamma - 1.0) * rho * (total_E - 0.5 * u * u)
        if p < p_min:
            p = p_min
            total_E = p / ((gamma - 1.0) * rho) + 0.5 * u * u
        U[0, i] = rho
        U[1, i] = rho * u
        U[2, i] = rho * total_E

# -------------------- Unit tests for HLLC algebra --------------------

def unit_tests(CONSTS):
    gamma = CONSTS['gamma']
    rho_min = CONSTS['rho_min']
    p_min = CONSTS['p_min']
    visc = CONSTS['visc']
    dx = (CONSTS['x1'] - CONSTS['x0']) / CONSTS['N']

    # Test 1: identical left/right states -> flux should equal physical flux
    rho = np.array([1.0])
    u = np.array([0.3])
    p = np.array([0.7])
    WL = np.vstack([rho, u, p])
    WR = WL.copy()
    F = hllc_flux(WL, WR, gamma, visc, rho_min, p_min, dx)
    # physical flux
    c = np.sqrt(gamma * p / rho)
    UL = np.vstack([rho, rho * u, rho * (0.5 * u * u + p / ((gamma - 1.0) * rho))])
    FL = np.vstack([rho * u,
                    rho * u * u + p,
                    (UL[2, :] + p) * u])
    err = np.max(np.abs(F - FL))
    if err > 1e-12:
        print(f"Unit test WARNING: identical-state HLLC flux mismatch err={err}")
    else:
        print("Unit test passed: identical-state HLLC flux matches physical flux")

# -------------------- Main solver --------------------

def solve_euler(CONSTS):
    gamma = CONSTS['gamma']
    CFL = CONSTS['CFL']
    visc = CONSTS['visc']
    rho_min = CONSTS['rho_min']
    p_min = CONSTS['p_min']
    N = CONSTS['N']
    x0 = CONSTS['x0']
    x1 = CONSTS['x1']
    snapshots = sorted(list(CONSTS['t_snapshots']))
    t_final = max(snapshots)

    U, x_cell, dx = initial_condition(N, x0, x1, gamma)

    # pre-allocate workspace arrays reused in compute_RHS
    W = np.zeros_like(U)
    slopes = np.zeros_like(U)
    WL_if = np.zeros((3, N + 1))
    WR_if = np.zeros((3, N + 1))
    F_if = np.zeros((3, N + 1))

    t = 0.0
    prev_t = 0.0
    results = {}

    mass0, mom0, E0 = compute_conserved_totals(U, dx)
    print(f"Initial totals -> mass: {mass0:.6f}, momentum: {mom0:.6e}, energy: {E0:.6f}")

    # run unit tests before time stepping
    unit_tests(CONSTS)

    # time loop
    while t < t_final - 1e-12:
        rho, u, p = primitives_from_U(U, gamma)
        c = np.sqrt(np.maximum(1e-12, gamma * p / rho))
        max_speed = np.max(np.abs(u) + c)
        dt = CFL * dx / max_speed
        # limit dt to not overshoot final time
        if t + dt > t_final:
            dt = t_final - t
        if dt <= 0.0:
            break

        # SSP-RK3
        U0 = U.copy()

        # Stage 1
        RHS, F_left, F_right = compute_RHS(U, dx, gamma, visc, rho_min, p_min, W, slopes, WL_if, WR_if, F_if)
        U1 = U0 + dt * RHS
        enforce_physical(U1, CONSTS)

        # Stage 2
        RHS, F_left2, F_right2 = compute_RHS(U1, dx, gamma, visc, rho_min, p_min, W, slopes, WL_if, WR_if, F_if)
        U2 = 0.75 * U0 + 0.25 * (U1 + dt * RHS)
        enforce_physical(U2, CONSTS)

        # Stage 3
        RHS, F_left3, F_right3 = compute_RHS(U2, dx, gamma, visc, rho_min, p_min, W, slopes, WL_if, WR_if, F_if)
        U_new = (1.0 / 3.0) * U0 + (2.0 / 3.0) * (U2 + dt * RHS)
        enforce_physical(U_new, CONSTS)

        # track net boundary flux estimate (use last evaluation's fluxes for diagnostics)
        F_left_net = F_left3.copy()
        F_right_net = F_right3.copy()

        U = U_new
        prev_t = t
        t += dt

        # store snapshots crossed by this step
        for s in snapshots:
            if prev_t < s <= t and s not in results:
                rho_s, u_s, p_s = primitives_from_U(U, gamma)
                results[s] = (rho_s.copy(), u_s.copy(), p_s.copy(), x_cell.copy())
                # diagnostics
                mass, mom, E = compute_conserved_totals(U, dx)
                grad_rho = np.abs(np.diff(rho_s)) / dx
                grad_p = np.abs(np.diff(p_s)) / dx
                combined = grad_rho * grad_p
                j = np.argmax(combined)
                x_shock = 0.5 * (x_cell[j] + x_cell[j + 1])
                print(f"t={s:.3f}: shock_position ~ {x_shock:.5f}, rho_min={rho_s.min():.6e}, rho_max={rho_s.max():.6f}, p_min={p_s.min():.6e}, p_max={p_s.max():.6f}")
                print(f"    Totals -> mass: {mass:.6f} (init {mass0:.6f}), momentum: {mom:.6e}, energy: {E:.6f}")
                # boundary flux diagnostics and CFL
                print(f"    boundary_flux_left: {F_left_net.tolist()}, boundary_flux_right: {F_right_net.tolist()}")
                print(f"    max_wave_speed: {max_speed:.6f}, dt: {dt:.6e}")

    return results

# -------------------- Plotting --------------------

def save_plots(results, CONSTS):
    times = sorted(results.keys())
    xcell = results[times[0]][3]

    plt.figure(figsize=(8, 4))
    for t in times:
        rho, u, p, x = results[t]
        plt.plot(x, rho, label=f"t={t:.2f}")
    plt.xlabel('x')
    plt.ylabel('Density')
    plt.title('Density at snapshots')
    plt.legend()
    plt.grid(True)
    plt.savefig('density_snapshots.png', dpi=150)
    plt.close()

    plt.figure(figsize=(8, 4))
    for t in times:
        rho, u, p, x = results[t]
        plt.plot(x, u, label=f"t={t:.2f}")
    plt.xlabel('x')
    plt.ylabel('Velocity')
    plt.title('Velocity at snapshots')
    plt.legend()
    plt.grid(True)
    plt.savefig('velocity_snapshots.png', dpi=150)
    plt.close()

    plt.figure(figsize=(8, 4))
    for t in times:
        rho, u, p, x = results[t]
        plt.plot(x, p, label=f"t={t:.2f}")
    plt.xlabel('x')
    plt.ylabel('Pressure')
    plt.title('Pressure at snapshots')
    plt.legend()
    plt.grid(True)
    plt.savefig('pressure_snapshots.png', dpi=150)
    plt.close()

# -------------------- Run solver and save outputs --------------------

def main():
    results = solve_euler(CONSTS)
    save_plots(results, CONSTS)
    print('Saved figures: density_snapshots.png, velocity_snapshots.png, pressure_snapshots.png')

if __name__ == '__main__':
    main()
```}


#### Script block3:

```python
#!/usr/bin/env python3
"""
Robust MUSCL-HLLC finite-volume solver for 1D Euler (Sod shock tube).
- SSP-RK3 time stepping
- Minmod limiter, preallocated work arrays
- HLLC with robust denominator regularization and runtime NaN/Inf checks
- Diagnostics: conservation residuals, flux norms, shock locator (max |d rho/dx|)
- Saves three PNG figures (density, velocity, pressure) at requested snapshots

All constants are passed explicitly via CONSTS.
"""

import numpy as np
import matplotlib.pyplot as plt

# -------------------- Constants --------------------
CONSTS = {
    'gamma': 1.4,
    'CFL': 0.20,            # lowered CFL for extra robustness
    'visc': 0.0,            # artificial viscosity coefficient
    'rho_min': 1e-8,
    'p_min': 1e-10,
    'N': 400,
    'x0': 0.0,
    'x1': 1.0,
    't_snapshots': [0.1, 0.2, 0.3],
    'debug': True,          # print extended diagnostics
    'apply_enforce': True,  # toggle positivity enforcement
    'report_enforce': True, # report energy change due to positivity fixer
}

# -------------------- Utility conversions --------------------

def prim_to_cons(rho, u, p, gamma):
    E = 0.5 * u * u + p / ((gamma - 1.0) * rho)
    return np.array([rho, rho * u, rho * (E)])


def cons_to_prim(U, gamma):
    rho = U[0]
    u = U[1] / rho
    E = U[2] / rho
    p = (gamma - 1.0) * rho * (E - 0.5 * u * u)
    return rho, u, p

# -------------------- Minmod slope limiter --------------------

def minmod(a, b):
    res = np.zeros_like(a)
    mask = (a * b) > 0.0
    res[mask] = np.where(np.abs(a[mask]) < np.abs(b[mask]), a[mask], b[mask])
    return res


def compute_slopes(W, slopes):
    # W shape: (3, N+2); slopes is preallocated same shape
    dfr = W[:, 2:] - W[:, 1:-1]
    dlf = W[:, 1:-1] - W[:, :-2]
    slope = minmod(dfr, dlf)
    slopes[:, :] = 0.0
    slopes[:, 1:-1] = slope
    return slopes

# -------------------- Boundary conditions --------------------

def apply_bc(U, bc_type='transmissive'):
    # transmissive/outflow or reflective (wall) bc
    if bc_type == 'transmissive':
        U[:, 0] = U[:, 1]
        U[:, -1] = U[:, -2]
    elif bc_type == 'reflective':
        U[0, 0] = U[0, 1]
        U[1, 0] = -U[1, 1]
        U[2, 0] = U[2, 1]
        U[0, -1] = U[0, -2]
        U[1, -1] = -U[1, -2]
        U[2, -1] = U[2, -2]
    else:
        raise ValueError('Unknown bc_type')

# -------------------- HLLC Riemann solver (vectorized) --------------------

def hllc_flux(WL, WR, gamma, visc, rho_min, p_min, dx, debug=False):
    """
    WL, WR: primitive states arrays [3, M] rows are rho,u,p
    Returns flux F of shape (3, M) and raises on algebraic breakdowns.
    """
    rhoL = np.maximum(WL[0, :], rho_min)
    uL = WL[1, :]
    pL = np.maximum(WL[2, :], p_min)

    rhoR = np.maximum(WR[0, :], rho_min)
    uR = WR[1, :]
    pR = np.maximum(WR[2, :], p_min)

    cL = np.sqrt(gamma * pL / rhoL)
    cR = np.sqrt(gamma * pR / rhoR)

    UL = np.vstack([rhoL, rhoL * uL, rhoL * (0.5 * uL * uL + pL / ((gamma - 1.0) * rhoL))])
    UR = np.vstack([rhoR, rhoR * uR, rhoR * (0.5 * uR * uR + pR / ((gamma - 1.0) * rhoR))])

    FL = np.vstack([rhoL * uL,
                    rhoL * uL * uL + pL,
                    (UL[2, :] + pL) * uL])
    FR = np.vstack([rhoR * uR,
                    rhoR * uR * uR + pR,
                    (UR[2, :] + pR) * uR])

    # Wave speed estimates (HLLC standard estimates)
    S_L = np.minimum(uL - cL, uR - cR)
    S_R = np.maximum(uL + cL, uR + cR)

    # Safeguard small separation
    S_L = np.minimum(S_L, np.minimum(uL, uR) - 1e-12)
    S_R = np.maximum(S_R, np.maximum(uL, uR) + 1e-12)

    # Compute middle wave speed S_M robustly
    num = pR - pL + rhoL * uL * (S_L - uL) - rhoR * uR * (S_R - uR)
    den = rhoL * (S_L - uL) - rhoR * (S_R - uR)
    small = 1e-12
    sign_den = np.sign(den)
    sign_den = np.where(sign_den == 0.0, 1.0, sign_den)
    den_safe = np.where(np.abs(den) < small, sign_den * small, den)

    S_M = num / den_safe

    # Basic algebraic checks
    if not np.all(np.isfinite(S_M)):
        raise FloatingPointError('Non-finite S_M in HLLC: {}'.format(S_M))

    # star densities
    rho_star_L = rhoL * (S_L - uL) / (S_L - S_M)
    rho_star_R = rhoR * (S_R - uR) / (S_R - S_M)

    # energies per mass
    E_L = UL[2, :] / rhoL
    E_R = UR[2, :] / rhoR

    # Toro (1997) style star energy calculation (robustified)
    # E*_k = E_k + (S_M - u_k) * (S_M + p_k/(rho_k*(S_k - u_k)))
    denomL = (S_L - uL)
    denomR = (S_R - uR)
    denomL_safe = np.where(np.abs(denomL) < small, np.sign(denomL) * small, denomL)
    denomR_safe = np.where(np.abs(denomR) < small, np.sign(denomR) * small, denomR)

    E_star_L = E_L + (S_M - uL) * (S_M + pL / (rhoL * denomL_safe))
    E_star_R = E_R + (S_M - uR) * (S_M + pR / (rhoR * denomR_safe))

    U_star_L = np.vstack([rho_star_L,
                          rho_star_L * S_M,
                          rho_star_L * E_star_L])
    U_star_R = np.vstack([rho_star_R,
                          rho_star_R * S_M,
                          rho_star_R * E_star_R])

    M = WL.shape[1]
    F = np.zeros_like(FL)

    mask_FR = S_R <= 0.0
    mask_FL = S_L >= 0.0
    mask_L = (S_L < 0.0) & (S_M >= 0.0)
    mask_R = (S_M < 0.0) & (S_R > 0.0)

    if np.any(mask_FL):
        F[:, mask_FL] = FL[:, mask_FL]
    if np.any(mask_FR):
        F[:, mask_FR] = FR[:, mask_FR]
    if np.any(mask_L):
        F[:, mask_L] = FL[:, mask_L] + S_L[mask_L] * (U_star_L[:, mask_L] - UL[:, mask_L])
    if np.any(mask_R):
        F[:, mask_R] = FR[:, mask_R] + S_R[mask_R] * (U_star_R[:, mask_R] - UR[:, mask_R])

    # optional viscous diffusive flux
    if visc is not None and visc > 0.0:
        diff = (UR - UL) / dx
        F = F - visc * diff

    # final safety checks
    if not np.all(np.isfinite(F)):
        raise FloatingPointError('Non-finite fluxes in HLLC')
    if np.any(rho_star_L <= 0.0) or np.any(rho_star_R <= 0.0):
        raise ValueError('Non-positive star density in HLLC')

    return F

# -------------------- RHS (spatial operator) --------------------

def compute_RHS(U, dx, gamma, visc, rho_min, p_min, W, slopes, WL_if, WR_if, F_if, bc_type='transmissive', debug=False):
    # U shape: (3, N+2). This routine does not mutate U externally.
    Ur = U.copy()
    apply_bc(Ur, bc_type=bc_type)

    # primitives into provided workspace W
    rho = Ur[0, :]
    u = Ur[1, :] / rho
    E = Ur[2, :] / rho
    p = (gamma - 1.0) * rho * (E - 0.5 * u * u)
    W[0, :] = rho
    W[1, :] = u
    W[2, :] = p

    # slopes (in-place)
    compute_slopes(W, slopes)

    # reconstruct interface states into WL_if, WR_if (shape (3, N+1))
    WL_if[:, :] = W[:, :-1] + 0.5 * slopes[:, :-1]
    WR_if[:, :] = W[:, 1:] - 0.5 * slopes[:, 1:]

    # fluxes at interfaces
    F = hllc_flux(WL_if, WR_if, gamma, visc, rho_min, p_min, dx, debug=debug)
    F_if[:, :] = F

    # Diagnostics on flux differences
    diff = F_if[:, 1:] - F_if[:, :-1]
    diff_max = np.max(np.abs(diff))
    diff_l1 = np.sum(np.abs(diff))
    if debug:
        print(f"    flux_diff: max={diff_max:.6e}, L1={diff_l1:.6e}")

    # Compute RHS for interior cells 1..N
    RHS = np.zeros_like(U)
    RHS[:, 1:-1] = -(F_if[:, 1:] - F_if[:, :-1]) / dx

    # NaN/Inf check
    if not np.all(np.isfinite(RHS)):
        raise FloatingPointError('Non-finite RHS computed')

    return RHS, F_if[:, 0].copy(), F_if[:, -1].copy()

# -------------------- Initial condition --------------------

def initial_condition(N, x0, x1, gamma):
    dx = (x1 - x0) / N
    x_cell = (np.arange(N) + 0.5) * dx + x0
    U = np.zeros((3, N + 2))
    rho = np.where(x_cell <= 0.5, 1.0, 0.125)
    u = np.zeros_like(rho)
    p = np.where(x_cell <= 0.5, 1.0, 0.1)
    for i in range(N):
        U[:, i + 1] = prim_to_cons(rho[i], u[i], p[i], gamma)
    return U, x_cell, dx

# -------------------- Diagnostics and manual integration --------------------

def compute_conserved_totals(U, dx):
    Uint = U[:, 1:-1]
    mass = np.sum(Uint[0, :] * dx)
    momentum = np.sum(Uint[1, :] * dx)
    energy = np.sum(Uint[2, :] * dx)
    return mass, momentum, energy


def primitives_from_U(U, gamma):
    rho = U[0, 1:-1]
    u = U[1, 1:-1] / rho
    E = U[2, 1:-1] / rho
    p = (gamma - 1.0) * rho * (E - 0.5 * u * u)
    return rho, u, p

# -------------------- Positivity enforcement --------------------

def enforce_physical(U, CONSTS, report=False):
    gamma = CONSTS['gamma']
    rho_min = CONSTS['rho_min']
    p_min = CONSTS['p_min']
    energy_before = np.sum(U[2, 1:-1])
    for i in range(1, U.shape[1] - 1):
        rho = max(U[0, i], rho_min)
        mom = U[1, i]
        E_tot = U[2, i]
        u = mom / rho
        total_E = E_tot / rho
        p = (gamma - 1.0) * rho * (total_E - 0.5 * u * u)
        if p < p_min:
            # Non-conservative fix: raise pressure by adjusting total energy only
            p = p_min
            total_E = p / ((gamma - 1.0) * rho) + 0.5 * u * u
        U[0, i] = rho
        U[1, i] = rho * u
        U[2, i] = rho * total_E
    energy_after = np.sum(U[2, 1:-1])
    if report:
        return energy_after - energy_before
    return 0.0

# -------------------- Unit tests for HLLC algebra --------------------

def unit_tests(CONSTS):
    gamma = CONSTS['gamma']
    rho_min = CONSTS['rho_min']
    p_min = CONSTS['p_min']
    visc = CONSTS['visc']
    dx = (CONSTS['x1'] - CONSTS['x0']) / CONSTS['N']

    # Test 1: identical left/right states -> flux should equal physical flux
    rho = np.array([1.0])
    u = np.array([0.3])
    p = np.array([0.7])
    WL = np.vstack([rho, u, p])
    WR = WL.copy()
    F = hllc_flux(WL, WR, gamma, visc, rho_min, p_min, dx)
    UL = np.vstack([rho, rho * u, rho * (0.5 * u * u + p / ((gamma - 1.0) * rho))])
    FL = np.vstack([rho * u,
                    rho * u * u + p,
                    (UL[2, :] + p) * u])
    err = np.max(np.abs(F - FL))
    if err > 1e-12:
        print(f"Unit test WARNING: identical-state HLLC flux mismatch err={err}")
    else:
        print("Unit test passed: identical-state HLLC flux matches physical flux")

    # Test 2: Sod-like strong jump checks
    WL = np.vstack([np.array([1.0]), np.array([0.0]), np.array([1.0])])
    WR = np.vstack([np.array([0.125]), np.array([0.0]), np.array([0.1])])
    try:
        F = hllc_flux(WL, WR, gamma, visc, rho_min, p_min, dx)
    except Exception as e:
        raise AssertionError('HLLC failed on Sod-like states: ' + str(e))
    # check finiteness
    if not np.all(np.isfinite(F)):
        raise AssertionError('HLLC returned non-finite flux for Sod-like test')
    print('Unit test passed: HLLC Sod-like flux finite')

# -------------------- Main solver --------------------

def solve_euler(CONSTS):
    gamma = CONSTS['gamma']
    CFL = CONSTS['CFL']
    visc = CONSTS['visc']
    rho_min = CONSTS['rho_min']
    p_min = CONSTS['p_min']
    N = CONSTS['N']
    x0 = CONSTS['x0']
    x1 = CONSTS['x1']
    snapshots = sorted(list(CONSTS['t_snapshots']))
    t_final = max(snapshots)
    debug = CONSTS.get('debug', False)

    U, x_cell, dx = initial_condition(N, x0, x1, gamma)

    # pre-allocate workspace arrays reused in compute_RHS
    W = np.zeros_like(U)
    slopes = np.zeros_like(U)
    WL_if = np.zeros((3, N + 1))
    WR_if = np.zeros((3, N + 1))
    F_if = np.zeros((3, N + 1))

    t = 0.0
    prev_t = 0.0
    results = {}

    mass0, mom0, E0 = compute_conserved_totals(U, dx)
    print(f"Initial totals -> mass: {mass0:.6f}, momentum: {mom0:.6e}, energy: {E0:.6f}")

    # run unit tests before time stepping
    unit_tests(CONSTS)

    # time loop
    step = 0
    while t < t_final - 1e-12:
        rho, u, p = primitives_from_U(U, gamma)
        c = np.sqrt(np.maximum(1e-12, gamma * p / rho))
        max_speed = np.max(np.abs(u) + c)
        dt = CFL * dx / max_speed
        # limit dt to not overshoot final time
        if t + dt > t_final:
            dt = t_final - t
        if dt <= 0.0:
            break

        U0 = U.copy()
        totals_before = compute_conserved_totals(U0, dx)

        # Stage 1
        RHS, F_left, F_right = compute_RHS(U, dx, gamma, visc, rho_min, p_min, W, slopes, WL_if, WR_if, F_if, debug=debug)
        if not np.all(np.isfinite(RHS)):
            raise FloatingPointError('RHS contained non-finite values at Stage 1')
        U1 = U0 + dt * RHS
        if CONSTS.get('apply_enforce', True):
            de = enforce_physical(U1, CONSTS, report=CONSTS.get('report_enforce', False))
            if debug and CONSTS.get('report_enforce', False):
                print(f"    enforce_physical Stage1 energy delta: {de:.6e}")

        # Stage 2
        RHS, F_left2, F_right2 = compute_RHS(U1, dx, gamma, visc, rho_min, p_min, W, slopes, WL_if, WR_if, F_if, debug=debug)
        U2 = 0.75 * U0 + 0.25 * (U1 + dt * RHS)
        if CONSTS.get('apply_enforce', True):
            de = enforce_physical(U2, CONSTS, report=CONSTS.get('report_enforce', False))
            if debug and CONSTS.get('report_enforce', False):
                print(f"    enforce_physical Stage2 energy delta: {de:.6e}")

        # Stage 3
        RHS, F_left3, F_right3 = compute_RHS(U2, dx, gamma, visc, rho_min, p_min, W, slopes, WL_if, WR_if, F_if, debug=debug)
        U_new = (1.0 / 3.0) * U0 + (2.0 / 3.0) * (U2 + dt * RHS)
        if CONSTS.get('apply_enforce', True):
            de = enforce_physical(U_new, CONSTS, report=CONSTS.get('report_enforce', False))
            if debug and CONSTS.get('report_enforce', False):
                print(f"    enforce_physical Stage3 energy delta: {de:.6e}")

        # NaN/Inf and positivity checks after full RK step
        if not np.all(np.isfinite(U_new)):
            raise FloatingPointError('Non-finite solution after RK update')
        if np.any(U_new[0, 1:-1] <= 0.0):
            raise ValueError('Non-positive density after RK update')

        # conservation residual check using last-evaluated boundary fluxes
        totals_after = compute_conserved_totals(U_new, dx)
        F_left_vec = F_left3
        F_right_vec = F_right3
        # residual = totals_new - totals_old + dt*(F_right - F_left)
        residual = np.array(totals_after) - np.array(totals_before) + dt * (F_right_vec - F_left_vec)
        res_norm = np.linalg.norm(residual, ord=1)
        if debug:
            print(f"step {step:04d}, t={t:.6f} -> t+dt={t+dt:.6f}: conservation_residual L1={res_norm:.3e}, residual={residual}")

        # update
        U = U_new
        prev_t = t
        t += dt
        step += 1

        # store snapshots crossed by this step
        for s in snapshots:
            if prev_t < s <= t and s not in results:
                rho_s, u_s, p_s = primitives_from_U(U, gamma)
                results[s] = (rho_s.copy(), u_s.copy(), p_s.copy(), x_cell.copy())
                # improved shock locator: index of max |d rho/dx|
                grad_rho = np.abs(np.diff(rho_s)) / dx
                j = np.argmax(grad_rho)
                x_shock = 0.5 * (x_cell[j] + x_cell[j + 1])
                print(f"t={s:.3f}: shock_position ~ {x_shock:.5f}, rho_min={rho_s.min():.6e}, rho_max={rho_s.max():.6f}, p_min={p_s.min():.6e}, p_max={p_s.max():.6f}")
                mass, mom, E = totals_after
                print(f"    Totals -> mass: {mass:.6f} (init {mass0:.6f}), momentum: {mom:.6e}, energy: {E:.6f}")
                print(f"    boundary_flux_left: {F_left_vec.tolist()}, boundary_flux_right: {F_right_vec.tolist()}")
                print(f"    max_wave_speed: {max_speed:.6f}, dt: {dt:.6e}")

    return results

# -------------------- Plotting --------------------

def save_plots(results, CONSTS):
    times = sorted(results.keys())
    xcell = results[times[0]][3]

    plt.figure(figsize=(8, 4))
    for t in times:
        rho, u, p, x = results[t]
        plt.plot(x, rho, label=f"t={t:.2f}")
    plt.xlabel('x')
    plt.ylabel('Density')
    plt.title('Density at snapshots')
    plt.legend()
    plt.grid(True)
    plt.savefig('density_snapshots.png', dpi=150)
    plt.close()

    plt.figure(figsize=(8, 4))
    for t in times:
        rho, u, p, x = results[t]
        plt.plot(x, u, label=f"t={t:.2f}")
    plt.xlabel('x')
    plt.ylabel('Velocity')
    plt.title('Velocity at snapshots')
    plt.legend()
    plt.grid(True)
    plt.savefig('velocity_snapshots.png', dpi=150)
    plt.close()

    plt.figure(figsize=(8, 4))
    for t in times:
        rho, u, p, x = results[t]
        plt.plot(x, p, label=f"t={t:.2f}")
    plt.xlabel('x')
    plt.ylabel('Pressure')
    plt.title('Pressure at snapshots')
    plt.legend()
    plt.grid(True)
    plt.savefig('pressure_snapshots.png', dpi=150)
    plt.close()

# -------------------- Run solver and save outputs --------------------

def main():
    results = solve_euler(CONSTS)
    save_plots(results, CONSTS)
    print('Saved figures: density_snapshots.png, velocity_snapshots.png, pressure_snapshots.png')

if __name__ == '__main__':
    main()
```



### runtime_outputs

#### Output block1

Initial totals -> mass: 0.562500, momentum: 0.000000e+00, energy: 1.375000
t=0.100: shock_position ~ 0.67750, rho_min=1.250000e-01, rho_max=1.000000, p_min=1.000000e-01, p_max=1.000000
    Totals -> mass: 0.562500 (init 0.562500), momentum: 9.000000e-02, energy: 1.375000
t=0.200: shock_position ~ 0.85250, rho_min=1.250000e-01, rho_max=1.000000, p_min=1.000000e-01, p_max=1.000000
    Totals -> mass: 0.562500 (init 0.562500), momentum: 1.800000e-01, energy: 1.375000
t=0.300: shock_position ~ 0.78000, rho_min=2.612068e-01, rho_max=1.000000, p_min=2.961842e-01, p_max=1.000000
    Totals -> mass: 0.558770 (init 0.562500), momentum: 2.636478e-01, energy: 1.358833
Saved figures: density_snapshots.png, velocity_snapshots.png, pressure_snapshots.png



#### Output block2

Initial totals -> mass: 0.562500, momentum: 0.000000e+00, energy: 1.375000
Unit test passed: identical-state HLLC flux matches physical flux
t=0.100: shock_position ~ 0.67750, rho_min=1.250000e-01, rho_max=1.000000, p_min=1.000000e-01, p_max=1.000000
    Totals -> mass: 0.562500 (init 0.562500), momentum: 9.005972e-02, energy: 1.375000
    boundary_flux_left: [0.0, 1.0, 0.0], boundary_flux_right: [0.0, 0.09999999999999998, 0.0]
    max_wave_speed: 2.192934, dt: 5.130113e-04
t=0.200: shock_position ~ 0.85250, rho_min=1.250000e-01, rho_max=1.000000, p_min=1.000000e-01, p_max=1.000000
    Totals -> mass: 0.562500 (init 0.562500), momentum: 1.801237e-01, energy: 1.375000
    boundary_flux_left: [0.0, 1.0, 0.0], boundary_flux_right: [0.0, 0.09999999999999998, 0.0]
    max_wave_speed: 2.191835, dt: 5.132685e-04
t=0.300: shock_position ~ 0.18250, rho_min=2.612246e-01, rho_max=1.000000, p_min=2.962161e-01, p_max=1.000000
    Totals -> mass: 0.558770 (init 0.562500), momentum: 2.636478e-01, energy: 1.358834
    boundary_flux_left: [0.0, 1.0, 0.0], boundary_flux_right: [0.247748223116389, 0.5311872843181369, 1.0947096706542383]
    max_wave_speed: 2.208417, dt: 3.551398e-04
Saved figures: density_snapshots.png, velocity_snapshots.png, pressure_snapshots.png



#### Output block3

Initial totals -> mass: 0.562500, momentum: 0.000000e+00, energy: 1.375000
Unit test passed: identical-state HLLC flux matches physical flux
Unit test passed: HLLC Sod-like flux finite
    flux_diff: max=1.161703e+00, L1=4.083927e+00
    enforce_physical Stage1 energy delta: 0.000000e+00
    flux_diff: max=1.077356e+00, L1=4.158104e+00
    enforce_physical Stage2 energy delta: 0.000000e+00
    flux_diff: max=1.111874e+00, L1=4.097315e+00
    enforce_physical Stage3 energy delta: 0.000000e+00
step 0000, t=0.000000 -> t+dt=0.000423: conservation_residual L1=2.220e-16, residual=[0.00000000e+00 0.00000000e+00 2.22044605e-16]
    flux_diff: max=1.064202e+00, L1=4.116669e+00
    enforce_physical Stage1 energy delta: 0.000000e+00
    flux_diff: max=9.744506e-01, L1=4.111631e+00
    enforce_physical Stage2 energy delta: 0.000000e+00
    flux_diff: max=1.017976e+00, L1=4.108614e+00
    enforce_physical Stage3 energy delta: 0.000000e+00
step 0001, t=0.000423 -> t+dt=0.000764: conservation_residu

... [Truncated: Content too long (586960 chars). Hidden middle part.] ...

.299859: conservation_residual L1=2.315e-11, residual=[-4.92141775e-13  4.10824819e-12  1.85497136e-11]
    flux_diff: max=1.866733e-02, L1=2.239285e+00
    enforce_physical Stage1 energy delta: 0.000000e+00
    flux_diff: max=1.866152e-02, L1=2.239298e+00
    enforce_physical Stage2 energy delta: 0.000000e+00
    flux_diff: max=1.866435e-02, L1=2.239287e+00
    enforce_physical Stage3 energy delta: 0.000000e+00
step 1311, t=0.299859 -> t+dt=0.300000: conservation_residual L1=7.371e-12, residual=[1.06673771e-13 1.49585377e-12 5.76815337e-12]
t=0.300: shock_position ~ 0.78000, rho_min=2.612141e-01, rho_max=1.000000, p_min=2.961985e-01, p_max=1.000000
    Totals -> mass: 0.558770 (init 0.562500), momentum: 2.636477e-01, energy: 1.358833
    boundary_flux_left: [0.0, 1.0, 0.0], boundary_flux_right: [0.2477526472769368, 0.5311850769149118, 1.0947135813119193]
    max_wave_speed: 2.208440, dt: 1.408477e-04
Saved figures: density_snapshots.png, velocity_snapshots.png, pressure_snapshots.png




### review_decision
accept

### review_comments
	Current Stage [A/3]
1. MAJOR

Category: accuracy
Issue: Global momentum and (to a lesser extent) energy are not conserved: reported momentum grows to O(10^-1) despite symmetric initial data (should remain ≈0). Mass also drifts slightly by t=0.3.

Recommendation: This is a MAJOR correctness issue. First instrument conservation diagnostics: compute net flux through domain boundaries each time step and print it to see whether momentum is created at boundaries or inside the domain. Add an assertion/test with two identical neighbouring cells to check that the numerical flux returns equal and opposite interface fluxes (F_{i+1/2} should be antisymmetric when WL==WR). If the flux is asymmetric, carefully re-derive and fix the HLLC formulas (especially the star-region U* and E* expressions). Compare your HLLC flux on simple states with a trusted reference implementation or textbook formula.


2. MINOR

Category: structure
Issue: Functions mutate inputs and rely on hidden globals: compute_RHS mutates U (sets ghost cells) and hllc_flux reads CONSTS['rho_min']/['p_min'] internally instead of receiving all constants explicitly.

Recommendation: Make boundary-condition application an explicit, separate function (e.g. apply_bc(U, CONSTS)) that is called prior to reconstruction/RHS evaluation. Ensure compute_RHS does not have side-effects beyond returning the RHS (or clearly document and control them). Pass rho_min and p_min (and any other constants) into hllc_flux via its argument list instead of referencing the global CONSTS inside the function so NameError/hidden-dependency problems are avoided and unit testing becomes straightforward.


3. MAJOR

Category: accuracy
Issue: Potential errors in the HLLC star-state energy and density calculations (E_star, rho_star) — these are sensitive and a small algebraic bug can break conservation and produce spurious net momentum.

Recommendation: Verify the algebra for S_L, S_R, S_M and the U*_L/U*_R and F*_L/F*_R formulas against a standard reference (Toro or other proven HLLC implementations). Add unit tests: (1) identical left/right states => flux == physical flux (no numerical artifact), (2) symmetric Riemann problem (Sod) but mirrored in code to verify zero net momentum. Consider replacing with a tested HLLC implementation or temporarily use a robust fallback (HLL or local Lax-Friedrichs/Rusanov) to isolate whether the issue lies in HLLC algebra.


4. MINOR

Category: accuracy
Issue: Shock-position detection and diagnostics are brittle (shock location found by max density gradient without guarding against picking the contact/expansion edge) and the snapshot selection condition is fragile.

Recommendation: Improve shock detection by combining density and pressure gradients or by searching for the leading shock where pressure increases strongly. Replace the snapshot trigger condition with a clearer approach: after advancing a time-step, check for any snapshot times t_snap s.t. previous_time < s <= current_time and then interpolate in time if required. Also print boundary fluxes and max CFL diagnostics to help debugging.


5. MINOR

Category: structure
Issue: Small coding issues: hllc_flux takes dx but uses global CONSTS for minima; WR is initialized with zeros and a trailing column left as zeros (harmless but confusing); unused imports (math.sqrt). The code sometimes accesses globals where explicit arguments were requested.

Recommendation: Clean up by passing all constants in function signatures (gamma, visc, rho_min, p_min). Initialize WR consistently to avoid confusion. Remove unused imports. Add unit tests and short reproducible test cases to detect regressions early.



	Current Stage [B/3]
1. MAJOR

Category: accuracy
Issue: Incorrect handling of small/zero denominator in S_M computation in hllc_flux leading to invalid star states and nonphysical fluxes.

Recommendation: Fix the denominator regularization in hllc_flux. The current line

den_safe = np.where(np.abs(den) < small, np.sign(den) * small + small, den)

is incorrect: when den is negative and small this expression yields zero (and when den==0 it yields small). That can produce S_M = num/0 or otherwise corrupt S_M and the star states, which will produce nonphysical fluxes and global conservation errors. Replace with a robust regularization such as

den_safe = np.where(np.abs(den) < small, np.sign(den) * small, den)

(or handle den==0 by falling back to a safe algebraic expression for S_M, e.g. an average velocity), and add explicit checks for NaN/Inf in S_M, rho_star_*, and U_star_* and a unit test that fails loudly if any NaN/Inf appears.


2. MAJOR

Category: accuracy
Issue: Substantial global conservation violations (growing momentum, mass loss) and appearance of large boundary fluxes indicate nonphysical flux computation or non-conservative corrections.

Recommendation: Investigate and eliminate the growing/non-physical net momentum and the mass/energy drift. Specific actionable steps:
- Add a diagnostic that tracks the discrete net flux across domain boundaries at every time step and verifies that mass/momentum/energy change equals negative of integrated boundary fluxes (i.e. check global conservation residual = totals_new - totals_old + dt*(F_right - F_left) ).
- Run the solver with closed/reflection BCs (or periodic) for a symmetric test; with symmetric ICs total momentum must remain ~0; if it does not, the flux assembly or sign convention is wrong.
- Instrument the RHS assembly to print the maximum and L1 norm of (F_if[:,1:] - F_if[:,:-1]) to ensure no rogue outliers appear.
- Temporarily disable enforce_physical to see whether the conservation errors come from flux computation or from the positivity fixer (enforce_physical modifies energy non-conservatively). If enforce_physical is required, replace it with a positivity-preserving conservative limiter (or apply minimal corrections that keep mass/momentum/energy balanced and document the residuals).


3. MAJOR

Category: structure
Issue: Insufficient algebraic/unit tests for HLLC and no runtime NaN/Inf/assert checks to catch algebraic breakdowns early.

Recommendation: Add robustness checks and unit tests for HLLC algebra and solver stability: 
- Extend unit_tests to cover strong left/right jumps (Sod-like states) and check that HLLC returns finite, physically admissible star densities (rho_star_* > rho_min) and pressures. Fail the test on NaN/Inf or negative densities/pressures.
- After fixing den_safe, add runtime asserts in debug mode to ensure no NaNs/Infs in fluxes or U after each stage and that rho>0 and p>0 in all cells.
- Lower CFL temporarily (e.g. 0.2) to check for stability; if solution stabilizes, the HLLC numerics or reconstruction may be marginal and needs revisiting.

These tests will localize whether the issue is algebraic (HLLC), reconstruction/limiter, or time stepping.


4. MINOR

Category: accuracy
Issue: Shock position diagnostic is fragile and picked inconsistent features; this makes qualitative checks unreliable.

Recommendation: Improve shock/feature diagnostics. The current shock locator uses combined grad_rho * grad_p which can pick the contact or other spurious feature and produces inconsistent positions (e.g. jumps reported in your run). Use a more robust metric such as the index of maximum |d rho/dx| (or smoothed absolute gradient) for the shock, optionally exclude the contact region by thresholding based on pressure jump. Also report the gradient magnitude and optionally plot a diagnostic of the gradient field when debugging.


5. MINOR

Category: structure
Issue: Positivity enforcement is opaque and there are missing runtime checks and documentation; unit testing and comments could be improved.

Recommendation: Small implementation and clarity improvements:
- In hllc_flux, guard E_star formulas and add comments with references for the specific algebra used (Toro, equations referenced). That will help future reviewers trace the algebraic signs.
- When enforcing positivity, document that this step is non-conservative (it modifies energy) and print the total energy change introduced by enforce_physical at each stage to quantify the correction.
- Add a quick NaN/Inf check after compute_RHS and after RK stages to fail loudly instead of silently producing wrong results.
- Consider adding a tighter unit test comparing to a reference Sod solution (e.g. using an exact solver) at low N to validate shock/contact positions quantitatively.

These will make debugging and validation much easier.



	Current Stage [C/3]
1. MINOR

Category: runtime
Issue: Code executed cleanly and produced saved figures; unit tests and runtime checks passed.

Recommendation: Accept. The code runs without errors (unit tests pass) and produces the requested output figures and diagnostics. Numerical method choices (MUSCL + minmod, HLLC, SSP-RK3) are appropriate for the Sod problem and were implemented in a robust, vectorized manner. Continue with small, optional refinements below if further accuracy or diagnostics are desired.


2. MINOR

Category: accuracy
Issue: Small but measurable total-mass/energy drift over the run; likely due to non-conservative positivity enforcement and boundary treatment.

Recommendation: Document and (optionally) reduce the small global conservation drift observed by t=0.3. The reported mass changed from 0.562500 to 0.558770 (~0.66% relative). This is not a correctness failure for an outflow/transmissive run, but it is worth noting in results and attributing causes (non-conservative positivity fixes, transmissive BCs, cumulative roundoff/flux at boundaries). If tighter global conservation is needed: (a) avoid non-conservative fixes where possible, (b) use conservative positivity limiters, or (c) track and correct boundary fluxes more carefully.


3. MINOR

Category: structure
Issue: Positivity enforcement is non-conservative and can alter global energy/mass; currently handled but should be clearly documented.

Recommendation: Make the non-conservative positivity fixer explicit in the user-facing documentation and consider switching to a conservative positivity-preserving limiter for production runs. The current enforce_physical modifies total energy non-conservatively (commented as such), which is acceptable for robustness but should be flagged since it changes global totals. If reporting of the energy correction is desired, ensure CONSTS['report_enforce'] is set to True (it is in the submission).


4. MINOR

Category: accuracy
Issue: Potential for sharper contact resolution and reduced smearing with alternative limiters or characteristic projection. No correctness problem, only a suggestion to boost solution quality.

Recommendation: Optional algorithmic improvements for accuracy (not required to accept): (a) characteristic-wise limiting or a less dissipative limiter (e.g. van Albada, MC) will sharpen contacts and reduce smearing; (b) small numerical viscosity (CONSTS['visc']>0) or a slope sensor could help when extremely strong shocks or noisy initial states appear; (c) compute L1/L2 errors vs exact Sod solution to quantify accuracy at snapshots for regression testing.


5. MINOR

Category: structure
Issue: Minor clarity/documentation improvements to make diagnostics easier to interpret for other developers/users.

Recommendation: Small code/clarity suggestion: add a short docstring comment above enforce_physical and the conservation-residual computation explaining the sign convention used and how boundary fluxes are incorporated into the residual. This will help future readers interpret the printed residuals and the boundary_flux vectors that are reported.





### iteration_count
3

### rev_image_description
	Current Stage [A/3]


	Current Stage [B/3]


	Current Stage [C/3]




