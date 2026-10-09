### 🛌🏻 Re4gent with Multi-Modal Review

#### Configs 🏷️:
MaxCount: 2, MultiModalReview: False, AlternativeNum: 2, MaxComments: 5, ConsModule_ON: True, CUT_OUTEXT: 2000

### prob_todo

The governing PDEs for 2D unsteady incompressible NS equations with explicit forcing are given by:
\begin{cases}
\frac{\partial u }{\partial t} + u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} + \frac{\partial p}{\partial x} - \frac{1}{\mathrm{Re}} \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right) = f_x, (x,y) \in \Omega,\\
\frac{\partial v }{\partial t} + u \frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} + \frac{\partial p}{\partial y} - \frac{1}{\mathrm{Re}} \left( \frac{\partial^2 v}{\partial x^2} + \frac{\partial^2 v}{\partial y^2} \right) = f_y, (x,y) \in \Omega,\\
\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0, (x,y) \in \Omega,\\
\end{cases}

The Reynolds number Re = 100.
The domain is $\Omega = [0, 2] \times [0, 1]$.
The forcing term $f = (f_x, f_y)$ is given by:
\begin{equation}
f_x = 0, \quad f_y = -\sin(\pi x) \sin(\pi y) \sin(\pi t).
\end{equation}

The boundary conditions are:
[Top/Bottom Walls] ($y=0$ and $y=1$): No-slip condition
$(u, v) = (0, 0)$.
[Inlet] (Left, $x=0$):
$u(0, y, t) = \sin(\pi y) \left(\sin(\pi t) + \sin(3\pi t) + \sin(5\pi t) \right)$
$v(0, y, t) = 0$
[Outlet] (Right, $x=2$):
Homogeneous Neumann for velocity: $\frac{\partial u}{\partial x} = 0, \frac{\partial v}{\partial x} = 0$
Dirichlet for pressure: $p(2, y, t) = 0$
The initial condition is $u(x, y, 0) = v(x, y, 0) = 0$ everywhere.

Implement a stable and efficient method to solve this problem.
Implement reasonable acceleration strategies to reduce computational cost.
Define visc = 1.0 / Re at the GLOBAL scope (outside main), and strictly use visc instead of nu throughout the code to avoid confusion with u.
Ensure all constants are passed explicitly to functions or defined globally to prevent NameError.
Simulate until t=1.0. Plot contours of u, v, and p at the final step in one figure using 'RdBu_r' colormap.
Just save figs do not use plt.show() in the code.

[HINTS]:
Never use the keyword argument 'tol' in any SciPy solver. Always use 'atol'. 
Do not use 'np.trapz' as it is removed in NumPy 2.0.
Use indexing='ij' for meshgrid and order='F' for all ravel/reshape operations to ensure consistent shapes and correct coordinate alignment.
Use Chorin's projection method on a Staggered (MAC) Grid using the upwind scheme to ensure stability.
Focus on generating runnable code with reasonable physical results, do not over-optimize.
Print concise progress information ONLY every 10% of total steps.


### expanded_prob
We must compute the unsteady 2D incompressible Navier–Stokes flow on Omega = [0,2]x[0,1] up to t=1.0 with Re=100 (visc = 1/Re). Forcing is zero in x and fy = -sin(pi x) sin(pi y) sin(pi t). Boundary conditions: no-slip on y=0,1; time-dependent Dirichlet inlet u(0,y,t)=sin(pi y)*(sin(pi t)+sin(3pi t)+sin(5pi t)), v(0,y,t)=0; homogeneous Neumann at outlet for velocity (du/dx = dv/dx = 0) and Dirichlet p=0 at x=2. Initial velocity is zero. The goal is an efficient, stable, reproducible solver producing u,v,p contours at final time (use 'RdBu_r' and save figures only). Implementation requirements/hints to respect: use a staggered (MAC) grid, Chorin projection method, upwind convective discretization, indexing='ij', order='F' for ravel/reshape, define visc globally, pass constants explicitly to functions, use atol (not tol) for SciPy solvers, avoid np.trapz, print progress every 10% of time steps, and favor stability/efficiency choices.

Primary mathematical and numerical challenges:
- Incompressibility constraint: numerical enforcement (projection) without pressure-velocity decoupling (avoid checkerboarding), hence staggered (MAC) grid.
- Pressure Poisson equation: large, sparse elliptic solve at every time step; proper boundary conditions (Dirichlet at outlet, mixed elsewhere) and compatibility with staggered discretization.
- Nonlinear convective term: stability-limited when treated explicitly; choice of discretization (upwind vs central) trades stability for accuracy; consistent face-based discretization on MAC is required.
- Time stepping: balance between stability (CFL and viscous diffusion constraints) and efficiency (implicitness reduces stiffness but increases cost per step). Semi-implicit (explicit convective, implicit viscous) or fully explicit RK SSP choices are natural.
- Discrete boundary conditions on staggered grid: mapping Dirichlet velocity/inlet profile and homogeneous Neumann outlet (ghost points) for staggered unknowns must be handled carefully to keep divergence continuity.
- Linear solver performance: Poisson and Helmholtz systems dominate runtime; need robust solvers and reuse of factorizations or multigrid/preconditioning to accelerate.
- Accuracy vs robustness: high-order convective schemes or flux limiters increase accuracy but complicate, while first-order upwind is robust but dissipative.
- Implementation reproducibility: global visc constant; explicit passing of all constants to functions; use of indexing/order conventions for array shapes and raveling.

All of the above must be reflected in the implementation plan so the result is physically plausible, stable to t=1.0, and computationally efficient for a moderate grid resolution.

### solution_plans
	Current Stage [A/2]
solu_name='Semi-implicit second-order projection on MAC (Crank–Nicolson viscous + Adams–Bashforth convective)' content="Governing idea:\n- Use a staggered (MAC) grid to avoid pressure–velocity decoupling. Treat convective terms explicitly with a stable upwind-biased discretization and viscous terms implicitly with Crank–Nicolson (second-order) to relax the diffusion timestep restriction. Use an incremental Chorin projection: predict velocity (Helmholtz solves), solve Poisson for pressure increment, correct velocity and update pressure.\n\nAlgorithmic steps (numbered):\n1) Grid and global constants: define Nx, Ny, dt with CFL estimator, create MAC grid (u at vertical faces, v at horizontal faces, p at cell centers) using meshgrid(..., indexing='ij'). Define visc = 1.0/Re globally.\n2) Initialize u,v,p arrays (Fortran order order='F' for ravel/reshape) and apply initial Dirichlet BCs (u=v=0). Precompute geometric arrays (dx, dy) and face-centered areas/volumes.\n3) Assemble constant sparse matrices once (and store factorizations if available):\n   - Helmholtz operator for u and v: I - (dt/2)*visc*Laplace (discretized on respective staggered unknowns with appropriate BCs). Coeffs constant if dt fixed; reuse factorization (sparse LU or AMG preconditioner). Use SciPy sparse formats; when using iterative solvers, pass atol and preconditioner.\n   - Pressure Poisson matrix: discrete Laplacian on cell centers with Dirichlet p at outlet (x=2) and Neumann/compatibility elsewhere consistent with staggered divergence; assemble once and factor or provide AMG preconditioner.\n4) Time stepping loop n=0..Nt-1 until t_final:\n   a) Compute explicit convective terms at faces using a monotone upwind flux (1st-order upwind) constructed from face velocities; for second-order accuracy in smooth regions use a limited reconstruction (optional) but default to upwind for stability.\n   b) Form RHS for predictor: for each velocity component use second-order Adams–Bashforth for convective terms (needs one previous convective field) and Crank–Nicolson for viscous part; assemble RHS_f = u^n + dt*(3/2 conv^n - 1/2 conv^{n-1}) + dt*(1/2 visc Lap u^n) (algebraic arrangement to form Helmholtz RHS).\n   c) Solve Helmholtz linear systems for intermediate face-centered velocities u* and v* using precomputed factorization or an iterative solver with recycled preconditioner. Use atol not tol in SciPy solvers.\n   d) Enforce boundary velocities on u* and v* consistent with physical BCs (inlet Dirichlet, walls no-slip via ghost points, outlet Neumann implemented in operator and RHS).\n   e) Compute divergence of u* at cell centers and solve pressure Poisson for the pressure increment phi: Laplace(phi) = (1/dt) div(u*). Use precomputed Poisson factorization or AMG preconditioner; impose p=0 at outlet cells explicitly and compatible Neumann elsewhere.\n   f) Correct velocities: u^{n+1} = u* - dt * grad_x(phi) (face-centered gradient), v^{n+1} = v* - dt * grad_y(phi). Update pressure: p^{n+1} = p^{n} + phi (incremental form).\n   g) Update convective histories for AB2 and apply physical velocity BCs exactly at faces (Dirichlet inlet/profile, no-slip walls). Compute diagnostics (CFL, max divergence). Print progress only when current step crosses 10%,20%,...,100% of total steps.\n5) After time loop, postprocess: interpolate to cell centers if required for plotting (consistent Fortran-order array shapes) and save figure with 3 contours (u,v,p) on subplots using 'RdBu_r' colormap. Do not call plt.show().\n\nStability / accuracy / complexity / efficiency limitations:\n- Stability: explicit convective AB2 or AB schemes still constrained by CFL: dt must satisfy CFL_u = max(|u|) dt/dx < C (C≈0.5 for AB2/upwind); implicit viscous Crank–Nicolson removes restrictive diffusive dt ~ dx^2/(2 visc) but convective limit remains. If CFL violated, reduce dt or use more implicit convective treatment.\n- Accuracy: scheme is formally second-order in time (CN + AB2) and second-order in space for viscous terms; with 1st-order upwind convective discretization convective accuracy drops to first order. Using a limited higher-order upwind (MUSCL) can restore higher accuracy at modest complexity.\n- Linear solvers cost: solving two Helmholtz systems + one Poisson per time step dominates CPU. Complexity per step depends on solver: direct sparse LU O(N^(3/2)) memory/time for 2D; iterative with good AMG preconditioner can be near O(N). Precomputing factorizations (if dt constant) drastically reduces per-step cost but increases memory.\n- Implementation complexity: assembling correct staggered operators and BC-consistent Poisson operator is nontrivial and must be carefully coded; ghost-cell handling must be consistent with indexing='ij' and order='F'.\n- Efficiency tips: reuse sparse factorizations or AMG preconditioners, vectorize array operations, use order='F' for ravel/reshape to match matrix storage; monitor and adapt dt via CFL to avoid unnecessary tiny steps."

	Current Stage [B/2]
solu_name='Explicit SSP-RK3 projection on MAC with algebraic multigrid Poisson (fully explicit convection + projection)' content="Governing idea:\n- Use a strong-stability-preserving explicit Runge–Kutta (SSP-RK3) for the convective and viscous terms (both treated explicitly) combined with Chorin projection at each RK substage. Use first-order upwind for convection to ensure monotonic behavior. Use an efficient AMG solver (pyamg or SciPy's multigrid if available) for the pressure Poisson to accelerate the elliptic solve. This approach simplifies implementation and avoids Helmholtz solves but requires smaller dt due to diffusion and convection CFL limits.\n\nAlgorithmic steps (numbered):\n1) Grid and constants: define MAC grid (indexing='ij'). Set visc = 1.0/Re globally. Choose dx, dy and compute dt from conservative CFL estimator dt = CFL * min(dx/maxU, dy/maxV, 0.5*dx*dx/visc) with a small safety factor; initialize velocities to zero.\n2) Preassemble pressure Poisson matrix on cell centers once and create an AMG preconditioner/factorization object if available; ensure Dirichlet p=0 at outlet and discrete Neumann elsewhere are encoded to keep solvable (compatibility enforced by projection RHS integral zero or adjust boundary ghost values).\n3) For each global time step until t_final do:\n   a) For each SSP-RK3 substage: compute convective fluxes at faces with first-order upwind, compute viscous Laplacian at faces with central second-order stencils, form explicit RHS for face velocities including forcing fx,fy.\n   b) Advance temporary face velocities to substage value explicitly per SSP-RK3 coefficients.\n   c) After the substage explicit update, enforce face Dirichlet BCs (inlet profile at x=0 changing in time, no-slip at y=0,1) and apply outlet Neumann via ghost values.\n   d) Compute divergence of the stage velocity field at cell centers and solve Poisson: Laplace(phi) = (1/dt) div(u_stage). Use AMG-preconditioned CG/BiCGSTAB with atol set (do not use tol). Apply Dirichlet p=0 at outlet in linear system.\n   e) Project: correct face velocities u_new = u_stage - dt * grad(phi) and accumulate stage solution per SSP-RK3 formula.\n   f) After full RK step, update p by adding phi (or set p to phi depending on projection variant) and enforce global BCs.\n   g) Print progress only at 10% increments.\n4) Postprocessing: map staggered fields for plotting (use Fortran-order arrangements), create a figure with three subplots showing contours of u, v, p with colormap 'RdBu_r' and save to disk (no plt.show()).\n\nStability / accuracy / complexity / efficiency limitations:\n- Stability: fully explicit treatment enforces dt <= min(C1*dx/|u|max, C2*dy/|v|max, C3*dx^2/visc). For Re=100 and fine grids the diffusion constraint can be severe (dt ∝ dx^2/visc), potentially requiring many time steps. SSP-RK3 keeps nonlinear stability for upwind convection but cannot avoid diffusion timestep restriction.\n- Accuracy: SSP-RK3 gives third-order strong stability for convective-dominated problems in smooth regions; however using first-order upwind for convection reduces spatial accuracy. Optionally replace upwind with a higher-order flux-limited scheme to improve accuracy at cost of complexity.\n- Linear solves: only a single Poisson solve per substage or per RK step is required (you may opt to project only once per full RK step to reduce cost at some accuracy penalty). AMG greatly reduces Poisson cost (near-linear scaling), but availability of pyamg or a good preconditioner is crucial.\n- Complexity: per-step cost is lower than semi-implicit Helmholtz+Poisson combined when using many iterative cycles; however, the smaller dt may offset that advantage by increasing number of steps. For moderate grid sizes and Re=100, this approach can be competitive and simpler to implement.\n- Implementation notes: projection after each RK substage increases Poisson solves (costly); to save cost, project only at RK stages' end but monitor divergence. Be careful with Neumann/Dirichlet BC compatibility on the Poisson solve; explicitly enforce mass-conservation compatibility to avoid drift.\n- Efficiency tips: reuse AMG hierarchy between solves, vectorize arrays, use order='F' for reshape/ravel, pass atol to SciPy solvers, and adjust print frequency to 10% step intervals as required."



### technical_spec
	Current Stage [A/3]
This script implements a 2D incompressible Navier–Stokes solver on a staggered MAC grid using a projection method with an explicit SSP-RK3 time integrator and first-order upwind treatment for convective terms (implemented via cell-centered upwind). The Poisson pressure solve is assembled once as a sparse matrix (mixed Dirichlet at outlet, Neumann elsewhere) and factorized for repeated solves. The code respects required indexing/order conventions (indexing='ij', order='F'), uses a global visc constant, enforces boundary conditions explicitly, prints progress every 10% of steps, and saves contours of u, v, and p at final time.

Main modules/functions:
- create_grid: builds MAC grid coordinates and spacings.
- assemble_poisson: builds the sparse Poisson matrix on cell centers with Dirichlet at outlet and Neumann elsewhere.
- apply_velocity_bcs: enforces Dirichlet/Neumann velocity BCs on staggered arrays.
- compute_cell_center_velocities: interpolates face velocities to cell centers.
- compute_convective_terms: constructs first-order upwind convective terms by computing derivatives on cell centers then averaging back to faces.
- laplacian_u/v: computes viscous Laplacian on faces with ghost treatments matching BCs.
- divergence: computes divergence at cell centers from face velocities.
- project: solves Poisson (using pre-factorized matrix) and corrects velocities, enforces pressure BC.
- time integrator: SSP-RK3 with projection at every substage for robustness.

Quantitative outputs: progress printouts at 10% intervals with time, step, max |u|, |v|, max divergence; final saved figure 'ns_mac_results.png' containing contour plots of u, v, p.

	Current Stage [B/3]
This script implements a 2D incompressible Navier–Stokes solver on a staggered MAC grid using an explicit SSP-RK3 integrator and a projection (Chorin) step. Key design points:
- Grid: MAC staggering (u on vertical faces, v on horizontal faces, p at cell centers). Arrays are created Fortran-contiguous to match Fortran-order indexing used for the Poisson assembly/solve.
- Poisson: assembled once on cell centers with mixed BC (p=0 at outlet x=2), assembled as the standard negative-definite Laplacian L (off-diagonals +1/dx^2, diag = -sum). The solver factorizes A once (LU) for repeated solves. The projection uses rhs = div/dt consistent with sign convention.
- Time-stepping: explicit SSP-RK3; convective terms are computed with a first-order upwind discretization implemented with array operations (vectorized) on cell-centered velocities; viscous terms are explicit and computed by vectorized second differences on faces.
- Boundary conditions: apply_velocity_bcs enforces inlet Dirichlet for u, zero v on walls/inlet, no-slip approximated by setting u adjacent to horizontal walls to zero, and outlet Neumann for u/v by mirroring the adjacent interior values. BC enforcement is applied consistently before any operator calls.
- Stability/safety: after each substage the code checks for overflows/NaNs and will abort if values exceed a safety threshold (1e6) with informative error.
- Consistency: all ravel/reshape operations use order='F'; arrays are Fortran-contiguous; asserts verify shapes and memory order when needed.
- Outputs: Progress prints every 10% of total steps with concise diagnostics; final contour figure of u, v, p saved to ns_mac_results.png (no plt.show()).

	Current Stage [C/3]
This script implements a 2D incompressible Navier–Stokes solver on a MAC staggered grid using an explicit SSP-RK3 integrator and a Chorin projection step. Key fixes and architecture:

- Global viscous coefficient: visc is defined at module scope and used everywhere.
- Grid: MAC staggering with Fortran-contiguous arrays and indexing='ij'. Mesh creation done once.
- Boundary masks: explicit boolean masks are built (u_dirichlet_mask, v_dirichlet_mask, p_dirichlet_mask) to identify which face velocities are Dirichlet (prescribed) and which cell-centered pressures are Dirichlet (outlet). These masks are used to assemble a masked Poisson matrix that corresponds to div( mask * grad ). This is the crucial fix: the Poisson operator reflects which faces the projection is allowed to change.
- Poisson assembly: assemble_poisson_masked constructs the sparse matrix once before the time loop. It sets exact identity rows for Dirichlet pressure nodes (outlet). Off-diagonals between adjacent cell-centers are included only if the face between them is correctable (mask True). Neumann boundaries are handled by not including coupling across faces that are non-correctable, which is consistent with the masked gradient operator. The resulting matrix is factorized once (LU) for repeated solves.
- Projection: project computes the divergence of the tentative velocity (full divergence) and solves A_masked * phi = div_full/dt. Because A_masked corresponds to divergence of the masked gradient, the projection only changes velocities at faces allowed by the masks. During the velocity correction step, only faces marked as correctable are updated; Dirichlet faces remain equal to their prescribed values (we still re-apply BCs after projection for consistency). This implements the review recommendation: contributions from prescribed faces are not treated as unknowns in the masked operator.
- Time stepping: SSP-RK3 with explicit convective (first-order upwind) and viscous terms; BCs are applied before operators. After each RK substage we project using the masked Poisson operator and only update correctable faces.
- Safety checks: overflow/NaN checks remain. Progress prints occur every 10% of total steps (concise). After each accepted time step the code prints max |div| to monitor divergence control.
- Output: final contours of u, v, p (cell-centered u_c, v_c, p) are plotted with 'RdBu_r' colormap and saved to ns_mac_results.png (no plt.show()).

Data flow summary: grid -> masks -> assemble masked Poisson -> time loop: apply BCs -> convective/viscous update -> tentative velocities -> compute full divergence -> solve masked Poisson for phi -> compute masked gradient corrections -> update only correctable faces -> accept stage -> proceed.

This approach fixes the main issues raised in the review by ensuring the Poisson operator and RHS are consistent with which face velocities are truly unknown for the projection, and by enforcing pressure Dirichlet rows as identity for robustness.



### python_codes

#### Script block1:

#!/usr/bin/env python
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time

# Global viscous coefficient
Re = 100.0
visc = 1.0 / Re

# ----------------------- Utility / Grid functions -----------------------
def create_grid(Nx, Ny, Lx=2.0, Ly=1.0):
    dx = Lx / float(Nx)
    dy = Ly / float(Ny)
    # u at vertical faces: x = 0..Nx (Nx+1), y centers at j+0.5
    xu = np.linspace(0.0, Lx, Nx+1)
    yu = (np.arange(Ny) + 0.5) * dy
    # v at horizontal faces: x centers at i+0.5, y = 0..Ny (Ny+1)
    xv = (np.arange(Nx) + 0.5) * dx
    yv = np.linspace(0.0, Ly, Ny+1)
    # p at cell centers
    xp = (np.arange(Nx) + 0.5) * dx
    yp = (np.arange(Ny) + 0.5) * dy
    Xp, Yp = np.meshgrid(xp, yp, indexing='ij')
    Xu, Yu = np.meshgrid(xu, yu, indexing='ij')
    Xv, Yv = np.meshgrid(xv, yv, indexing='ij')
    return dx, dy, Xu, Yu, Xv, Yv, Xp, Yp

# ----------------------- BC and field helpers -----------------------
def inlet_u_profile(y, t):
    return np.sin(np.pi * y) * (np.sin(np.pi * t) + np.sin(3*np.pi*t) + np.sin(5*np.pi*t))

def apply_velocity_bcs(u, v, Xu, Yu, Xv, Yv, t):
    # u: shape (Nx+1, Ny); v: shape (Nx, Ny+1)
    # Inlet (x=0): u at i=0 Dirichlet; v at i=0 column of v faces corresponds to x ~ dx/2, but problem states v(0,y,t)=0; enforce v leftmost ghost/face near inlet to zero
    # No-slip top/bottom: v at j=0 and j=Ny = 0; ensure u has antisymmetric ghost for walls handled in Laplacian calculation; but we will explicitly set v boundaries and enforce u at walls via ghost logic in laplacian
    # Outlet (x=2): homogeneous Neumann for u and v: implemented implicitly via ghost mirroring in laplacian and gradient calculations
    # Apply inlet u
    # Xu has x positions for u faces; inlet corresponds to Xu==0 which is i=0
    # Build inlet u profile for all y positions in Yu[:,]
    u_in = inlet_u_profile(Yu[0,:], t)  # Yu[0,:] corresponds to x=0 column y coordinates
    u[0, :] = u_in
    # Enforce v Dirichlet at inlet (v at x~dx/2 leftmost index i=0)
    v[0, :] = 0.0
    # Enforce v at top/bottom walls to zero (j=0 and j=Ny)
    v[:, 0] = 0.0
    v[:, -1] = 0.0
    # For u top/bottom walls, enforce no-slip in the following manner:
    # u does not sit on the wall; horizontal walls require ghost symmetry; we'll not directly set u at j boundary but can ensure correct values via apply after projection
    return

# ----------------------- Discrete operators -----------------------
def compute_cell_center_velocities(u, v):
    # interpolate face velocities to cell centers
    # u_c shape (Nx,Ny) = average of u[i] and u[i+1]
    u_c = 0.5 * (u[0:-1, :] + u[1:, :])
    # v_c shape (Nx,Ny) = average of v[:, j] and v[:, j+1]
    v_c = 0.5 * (v[:, 0:-1] + v[:, 1:])
    return u_c, v_c


def compute_convective_terms(u, v, dx, dy):
    # Use cell-centered upwind approach: compute derivatives on cell-centered velocities and map back to faces by averaging.
    Nx_p = u.shape[0] - 1
    Ny_p = u.shape[1]
    # cell-centered velocities
    u_c, v_c = compute_cell_center_velocities(u, v)
    # Prepare derivative arrays at cell centers
    du_dx = np.zeros_like(u_c)
    du_dy = np.zeros_like(u_c)
    dv_dx = np.zeros_like(v_c)
    dv_dy = np.zeros_like(v_c)
    # upwind for du/dx
    for i in range(Nx_p):
        # For interior i use neighbors
        for j in range(Ny_p):
            uc = u_c[i, j]
            # du/dx
            if i == 0:
                # forward difference
                du_dx[i, j] = (u_c[i+1, j] - u_c[i, j]) / dx if Nx_p > 1 else 0.0
            elif i == Nx_p-1:
                # backward difference
                du_dx[i, j] = (u_c[i, j] - u_c[i-1, j]) / dx
            else:
                if uc > 0:
                    du_dx[i, j] = (u_c[i, j] - u_c[i-1, j]) / dx
                else:
                    du_dx[i, j] = (u_c[i+1, j] - u_c[i, j]) / dx
            # du/dy
            if j == 0:
                du_dy[i, j] = (u_c[i, j+1] - u_c[i, j]) / dy
            elif j == Ny_p-1:
                du_dy[i, j] = (u_c[i, j] - u_c[i, j-1]) / dy
            else:
                vc = v_c[i, j]
                if vc > 0:
                    du_dy[i, j] = (u_c[i, j] - u_c[i, j-1]) / dy
                else:
                    du_dy[i, j] = (u_c[i, j+1] - u_c[i, j]) / dy
            # dv/dx
            if i == 0:
                dv_dx[i, j] = (v_c[i+1, j] - v_c[i, j]) / dx if Nx_p > 1 else 0.0
            elif i == Nx_p-1:
                dv_dx[i, j] = (v_c[i, j] - v_c[i-1, j]) / dx
            else:
                uc2 = u_c[i, j]
                if uc2 > 0:
                    dv_dx[i, j] = (v_c[i, j] - v_c[i-1, j]) / dx
                else:
                    dv_dx[i, j] = (v_c[i+1, j] - v_c[i, j]) / dx
            # dv/dy
            if j == 0:
                dv_dy[i, j] = (v_c[i, j+1] - v_c[i, j]) / dy
            elif j == Ny_p-1:
                dv_dy[i, j] = (v_c[i, j] - v_c[i, j-1]) / dy
            else:
                if v_c[i, j] > 0:
                    dv_dy[i, j] = (v_c[i, j] - v_c[i, j-1]) / dy
                else:
                    dv_dy[i, j] = (v_c[i, j+1] - v_c[i, j]) / dy
    # convective terms at cell centers
    conv_u_c = u_c * du_dx + v_c * du_dy
    conv_v_c = u_c * dv_dx + v_c * dv_dy
    # map convective terms back to faces by averaging neighboring cell-centered convective values
    # conv_u_face shape (Nx+1, Ny) -> average conv_u_c left and right cells where available
    conv_u = np.zeros_like(u)
    # interior faces i=1..Nx-1
    conv_u[1:-1, :] = 0.5 * (conv_u_c[0:-1, :] + conv_u_c[1:, :])
    # leftmost face (i=0): use conv_u_c[0]
    conv_u[0, :] = conv_u_c[0, :]
    # rightmost face (i=Nx): use conv_u_c[-1]
    conv_u[-1, :] = conv_u_c[-1, :]
    # conv_v_face shape (Nx, Ny+1) average vertically
    conv_v = np.zeros_like(v)
    conv_v[:, 1:-1] = 0.5 * (conv_v_c[:, 0:-1] + conv_v_c[:, 1:])
    conv_v[:, 0] = conv_v_c[:, 0]
    conv_v[:, -1] = conv_v_c[:, -1]
    return conv_u, conv_v


def laplacian_u(u, dx, dy):
    # compute Laplacian at u faces with BC ghosting logic
    Nx_u, Ny = u.shape
    Lu = np.zeros_like(u)
    # x-second derivative
    for i in range(Nx_u):
        if i == 0:
            # left boundary Dirichlet (inlet) -> skip interior update (we will not use Laplacian for i=0); set to zero
            continue
        elif i == Nx_u - 1:
            # rightmost face Neumann mirrored ghost u[Nx+1] = u[Nx-1]
            # d2/dx2 = (u[i-1] - 2u[i] + u_ghost)/dx^2 = (u[i-1] - 2u[i] + u[i-1]) / dx^2
            d2x = 2.0 * (u[i-1, :] - u[i, :]) / dx**2
        else:
            d2x = (u[i+1, :] - 2.0 * u[i, :] + u[i-1, :]) / dx**2
        Lu[i, :] += d2x
    # y-second derivative, enforce no-slip walls via antisymmetric ghost: u_ghost = -u_neighbor
    for j in range(Ny):
        if j == 0:
            # ghost j-1 = -u[:,0]
            d2y = (u[:, 1] - 2.0 * u[:, 0] + (-u[:, 0])) / dy**2
        elif j == Ny - 1:
            d2y = ((-u[:, -1]) - 2.0 * u[:, -1] + u[:, -2]) / dy**2
        else:
            d2y = (u[:, j+1] - 2.0 * u[:, j] + u[:, j-1]) / dy**2
        Lu[:, j] += d2y
    return Lu


def laplacian_v(v, dx, dy):
    Nx, Ny_v = v.shape
    Lv = np.zeros_like(v)
    # y-second derivative
    for j in range(Ny_v):
        if j == 0:
            # bottom wall v=0 Dirichlet -> no update (we won't update v at j=0)
            continue
        elif j == Ny_v - 1:
            # top wall v=0 Dirichlet -> skip
            continue
        else:
            d2y = (v[:, j+1] - 2.0 * v[:, j] + v[:, j-1]) / dy**2
        Lv[:, j] += d2y
    # x-second derivative, enforce inlet/outlet BCs: inlet v Dirichlet zero at i approx 0 -> treat i=0 as Dirichlet
    for i in range(Nx):
        if i == 0:
            continue
        elif i == Nx - 1:
            # rightmost: Neumann mirrored ghost v_ghost = v[Nx-2]
            d2x = 2.0 * (v[i-1, :] - v[i, :]) / dx**2
        else:
            d2x = (v[i+1, :] - 2.0 * v[i, :] + v[i-1, :]) / dx**2
        Lv[i, :] += d2x
    return Lv

# ----------------------- Pressure Poisson assembly -----------------------

def assemble_poisson(Nx, Ny, dx, dy):
    # Build A matrix for -Laplace (positive definite) on cell centers with mixed BCs
    N = Nx * Ny
    data = []
    rows = []
    cols = []
    def idx(i, j):
        return i + j * Nx  # Fortran order mapping (order='F')
    for j in range(Ny):
        for i in range(Nx):
            row = idx(i, j)
            diag = 0.0
            # left neighbor (i-1)
            if i - 1 >= 0:
                rows.append(row); cols.append(idx(i-1, j)); data.append(-1.0 / dx**2); diag += 1.0 / dx**2
            else:
                # inlet boundary (Neumann): skip (ghost = self) -> no contribution
                pass
            # right neighbor (i+1)
            if i + 1 <= Nx - 1:
                rows.append(row); cols.append(idx(i+1, j)); data.append(-1.0 / dx**2); diag += 1.0 / dx**2
            else:
                # outlet boundary at x=2: Dirichlet p = 0
                # contribute diag term (1/dx^2) and RHS will get p_b*1/dx^2 but p_b=0 -> only diag increment
                diag += 1.0 / dx**2
            # down neighbor (j-1)
            if j - 1 >= 0:
                rows.append(row); cols.append(idx(i, j-1)); data.append(-1.0 / dy**2); diag += 1.0 / dy**2
            else:
                # bottom wall Neumann: skip
                pass
            # up neighbor (j+1)
            if j + 1 <= Ny - 1:
                rows.append(row); cols.append(idx(i, j+1)); data.append(-1.0 / dy**2); diag += 1.0 / dy**2
            else:
                # top wall Neumann: skip
                pass
            # diagonal
            rows.append(row); cols.append(row); data.append(diag)
    A = sp.csr_matrix((data, (rows, cols)), shape=(N, N))
    return A

# ----------------------- Divergence and projection -----------------------

def divergence(u, v, dx, dy):
    # u shape (Nx+1, Ny), v shape (Nx, Ny+1)
    div = (u[1:, :] - u[0:-1, :]) / dx + (v[:, 1:] - v[:, 0:-1]) / dy
    return div


def project(u, v, dx, dy, dt, poisson_solver, Nx, Ny):
    # compute divergence and solve Poisson for phi, then correct velocities
    div = divergence(u, v, dx, dy)
    rhs = (div / dt).ravel(order='F')
    # Solve A * phi = rhs, with A assembled as positive definite earlier. Note A was assembled as -Laplace positive diag.
    phi_vec = poisson_solver(rhs)
    phi = phi_vec.reshape((Nx, Ny), order='F')
    # Correct u: for u faces i=1..Nx-1 use phi[i] - phi[i-1]
    # u shape (Nx+1, Ny)
    # interior faces
    gradx = np.zeros_like(u)
    # for i=1..Nx-1
    gradx[1:-1, :] = (phi[1:, :] - phi[0:-1, :]) / dx
    # leftmost face i=0 is inlet Dirichlet and not corrected
    # rightmost face i=Nx uses phi_outlet (p at x=2) which is Dirichlet 0
    gradx[-1, :] = (0.0 - phi[-1, :]) / dx
    u -= dt * gradx
    # Correct v: for v faces j=1..Ny-1
    grady = np.zeros_like(v)
    grady[:, 1:-1] = (phi[:, 1:] - phi[:, 0:-1]) / dy
    # bottom and top faces are Dirichlet v=0 and not corrected
    v -= dt * grady
    # Update pressure by adding phi (incremental form)
    p = phi.copy()
    return u, v, p

# ----------------------- Main solver -----------------------

def main():
    # Grid and time parameters
    Nx = 80  # number of pressure cells in x
    Ny = 40  # number of pressure cells in y
    Lx = 2.0
    Ly = 1.0
    dx, dy, Xu, Yu, Xv, Yv, Xp, Yp = create_grid(Nx, Ny, Lx, Ly)
    # Time step selection via CFL estimate
    # characteristic velocity scale ~ 3 (sum of sines). Use safety factor
    maxU_est = 3.0
    CFL = 0.4
    dt_conv = CFL * min(dx / maxU_est, dy / maxU_est)
    # viscous constraint for explicit diffusion
    dt_diff = 0.25 * min(dx*dx, dy*dy) / visc
    dt = min(dt_conv, dt_diff)
    # restrict dt modestly
    dt = min(dt, 0.008)
    t_final = 1.0
    Nt = int(np.ceil(t_final / dt))
    dt = t_final / Nt  # adjust dt to match final time exactly
    print(f"Grid: Nx={Nx}, Ny={Ny}, dx={dx:.4e}, dy={dy:.4e}, dt={dt:.4e}, steps={Nt}")

    # Initialize fields
    u = np.zeros((Nx+1, Ny), order='F')
    v = np.zeros((Nx, Ny+1), order='F')
    p = np.zeros((Nx, Ny), order='F')

    # Preassemble Poisson matrix and factorize
    A = assemble_poisson(Nx, Ny, dx, dy)
    A_csc = A.tocsc()
    lu = spla.splu(A_csc)
    def poisson_solver(rhs):
        # solve using LU factorization; rhs expected 1D numpy array
        return lu.solve(rhs)

    # Time stepping: SSP-RK3 with projection at each substage
    t = 0.0
    # For printing progress every 10%
    next_print = 0.1
    print_interval = max(1, Nt // 10)
    start_time = time.time()

    # initial BC enforcement
    apply_velocity_bcs(u, v, Xu, Yu, Xv, Yv, t)

    for step in range(1, Nt+1):
        t_stage = t
        # SSP-RK3 coefficients
        u0 = u.copy(order='F')
        v0 = v.copy(order='F')
        p_accum = np.zeros_like(p)
        # Three stages
        u_stage = u.copy(order='F')
        v_stage = v.copy(order='F')
        for stage in range(3):
            if stage == 0:
                a = 1.0
                b = 0.0
                dt_stage = dt
                t_sub = t + 0.0 * dt
            elif stage == 1:
                a = 3.0 / 4.0
                b = 1.0 / 4.0
                dt_stage = dt
                t_sub = t + dt
            else:
                a = 1.0 / 3.0
                b = 2.0 / 3.0
                dt_stage = dt
                t_sub = t + dt
            # compute convective terms and viscous Laplacian
            conv_u, conv_v = compute_convective_terms(u_stage, v_stage, dx, dy)
            Lu = laplacian_u(u_stage, dx, dy)
            Lv = laplacian_v(v_stage, dx, dy)
            # external forcing
            fx = np.zeros_like(u)
            # fy = -sin(pi x) sin(pi y) sin(pi t) at v-face positions: we need fy at v faces (x centers xv, y faces yv)
            # Xv, Yv are positions for v faces
            fy = -np.sin(np.pi * Xv) * np.sin(np.pi * Yv) * np.sin(np.pi * t_sub)
            # explicit update for velocities (forward Euler sub-step within RK stage)
            u_tent = u_stage + dt_stage * ( - conv_u + visc * Lu + fx )
            v_tent = v_stage + dt_stage * ( - conv_v + visc * Lv + fy )
            # apply BCs before projection
            apply_velocity_bcs(u_tent, v_tent, Xu, Yu, Xv, Yv, t_sub)
            # Project to divergence-free field
            u_proj, v_proj, phi = project(u_tent.copy(order='F'), v_tent.copy(order='F'), dx, dy, dt_stage, poisson_solver, Nx, Ny)
            # after projection enforce velocity BCs again (Dirichlet nodes)
            apply_velocity_bcs(u_proj, v_proj, Xu, Yu, Xv, Yv, t_sub)
            # SSP combination
            if stage == 0:
                u_stage = u_proj.copy(order='F')
                v_stage = v_proj.copy(order='F')
                p_accum = phi.copy(order='F')
            elif stage == 1:
                u_stage = a * u0 + b * u_proj
                v_stage = a * v0 + b * v_proj
                p_accum = a * p_accum + b * phi
            else:
                u_stage = a * u0 + b * u_proj
                v_stage = a * v0 + b * v_proj
                p_accum = a * p_accum + b * phi
        # End of RK stages; accept u_stage/v_stage as new fields
        u[:] = u_stage
        v[:] = v_stage
        # update pressure field: incremental add of p_accum (we choose to set p = p_accum for convenience)
        p[:] = p_accum
        t += dt
        # Diagnostics and progress
        if step % print_interval == 0 or step == Nt:
            div_max = np.max(np.abs(divergence(u, v, dx, dy)))
            print(f"Step {step}/{Nt}, t={t:.4f}, max|u|={np.max(np.abs(u)):.4f}, max|v|={np.max(np.abs(v)):.4f}, max|div|={div_max:.3e}")
    total_time = time.time() - start_time
    print(f"Simulation completed in {total_time:.2f} s")

    # Postprocess: interpolate u and v to cell centers for plotting
    u_c, v_c = compute_cell_center_velocities(u, v)
    # Plot contours of u, v, p on cell-centered grid
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)
    levels = 50
    cf0 = axes[0].contourf(Xp, Yp, u_c, levels=levels, cmap='RdBu_r')
    axes[0].set_title('u (cell-centered) at t=1.0')
    axes[0].set_xlabel('x'); axes[0].set_ylabel('y')
    fig.colorbar(cf0, ax=axes[0])
    cf1 = axes[1].contourf(Xp, Yp, v_c, levels=levels, cmap='RdBu_r')
    axes[1].set_title('v (cell-centered) at t=1.0')
    axes[1].set_xlabel('x'); axes[1].set_ylabel('y')
    fig.colorbar(cf1, ax=axes[1])
    cf2 = axes[2].contourf(Xp, Yp, p, levels=levels, cmap='RdBu_r')
    axes[2].set_title('p (cell-centered) at t=1.0')
    axes[2].set_xlabel('x'); axes[2].set_ylabel('y')
    fig.colorbar(cf2, ax=axes[2])
    plt.suptitle('2D Incompressible Navier-Stokes (MAC, SSP-RK3 projection)')
    plt.savefig('ns_mac_results.png', dpi=150)

if __name__ == '__main__':
    main()



#### Script block2:

#!/usr/bin/env python
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time

# Global viscous coefficient
Re = 100.0
visc = 1.0 / Re

# Safety threshold for runaway values
SAFETY_MAX = 1e6

# ----------------------- Utility / Grid functions -----------------------
def create_grid(Nx, Ny, Lx=2.0, Ly=1.0):
    dx = Lx / float(Nx)
    dy = Ly / float(Ny)
    xu = np.linspace(0.0, Lx, Nx+1)
    yu = (np.arange(Ny) + 0.5) * dy
    xv = (np.arange(Nx) + 0.5) * dx
    yv = np.linspace(0.0, Ly, Ny+1)
    xp = (np.arange(Nx) + 0.5) * dx
    yp = (np.arange(Ny) + 0.5) * dy
    Xu, Yu = np.meshgrid(xu, yu, indexing='ij')  # (Nx+1, Ny)
    Xv, Yv = np.meshgrid(xv, yv, indexing='ij')  # (Nx, Ny+1)
    Xp, Yp = np.meshgrid(xp, yp, indexing='ij')  # (Nx, Ny)
    # ensure Fortran order for coordinate arrays used in vector ops (not strictly necessary but consistent)
    return dx, dy, np.asfortranarray(Xu), np.asfortranarray(Yu), np.asfortranarray(Xv), np.asfortranarray(Yv), np.asfortranarray(Xp), np.asfortranarray(Yp)

# ----------------------- BC and field helpers -----------------------
def inlet_u_profile(y, t):
    return np.sin(np.pi * y) * (np.sin(np.pi * t) + np.sin(3*np.pi*t) + np.sin(5*np.pi*t))


def apply_velocity_bcs(u, v, Xu, Yu, Xv, Yv, t):
    # u: (Nx+1, Ny), v: (Nx, Ny+1)
    # Enforce inlet u (Dirichlet) at i=0
    u[0, :] = inlet_u_profile(Yu[0, :], t)
    # Enforce no-slip at horizontal walls (approximate by setting tangential u at adjacent faces to zero)
    u[:, 0] = 0.0
    u[:, -1] = 0.0
    # Outlet (rightmost face) Neumann for u: approximate by copying interior neighbor
    u[-1, :] = u[-2, :]
    # v: enforce inlet v ~ 0 at leftmost face column i=0 (v is defined at x centers, leftmost is near inlet)
    v[0, :] = 0.0
    # bottom and top walls: v at j=0 and j=Ny are exactly on the walls -> Dirichlet zero
    v[:, 0] = 0.0
    v[:, -1] = 0.0
    # outlet Neumann for v (right boundary): copy neighbor
    v[-1, :] = v[-2, :]
    return

# ----------------------- Discrete operators -----------------------

def compute_cell_center_velocities(u, v):
    # u: (Nx+1, Ny) -> u_c: (Nx, Ny) = avg of u[i] and u[i+1]
    # v: (Nx, Ny+1) -> v_c: (Nx, Ny) = avg of v[:, j] and v[:, j+1]
    u_c = 0.5 * (u[0:-1, :] + u[1:, :])
    v_c = 0.5 * (v[:, 0:-1] + v[:, 1:])
    return np.asfortranarray(u_c), np.asfortranarray(v_c)


def compute_convective_terms(u, v, dx, dy):
    # Vectorized first-order upwind on cell-centered velocities then map back to faces by averaging
    Nx = u.shape[0] - 1
    Ny = u.shape[1]
    uc, vc = compute_cell_center_velocities(u, v)  # (Nx,Ny)
    # Prepare shifted arrays
    # For x-derivatives: left (i-1), center (i), right (i+1)
    uc_left = np.empty_like(uc)
    uc_right = np.empty_like(uc)
    # left shift: for i=0, left neighbor use same value (one-sided) -> will be replaced by forward diff
    uc_left[1:, :] = uc[0:-1, :]
    uc_left[0, :] = uc[0, :]
    uc_right[0:-1, :] = uc[1:, :]
    uc_right[-1, :] = uc[-1, :]
    # For y-derivatives: down (j-1), up (j+1)
    uc_down = np.empty_like(uc)
    uc_up = np.empty_like(uc)
    uc_down[:, 1:] = uc[:, 0:-1]
    uc_down[:, 0] = uc[:, 0]
    uc_up[:, 0:-1] = uc[:, 1:]
    uc_up[:, -1] = uc[:, -1]
    # dv arrays
    vc_left = np.empty_like(vc)
    vc_right = np.empty_like(vc)
    vc_left[1:, :] = vc[0:-1, :]
    vc_left[0, :] = vc[0, :]
    vc_right[0:-1, :] = vc[1:, :]
    vc_right[-1, :] = vc[-1, :]
    vc_down = np.empty_like(vc)
    vc_up = np.empty_like(vc)
    vc_down[:, 1:] = vc[:, 0:-1]
    vc_down[:, 0] = vc[:, 0]
    vc_up[:, 0:-1] = vc[:, 1:]
    vc_up[:, -1] = vc[:, -1]
    # compute derivatives using upwind logic
    du_dx = np.where(uc > 0.0, (uc - uc_left) / dx, (uc_right - uc) / dx)
    du_dy = np.where(vc > 0.0, (uc - uc_down) / dy, (uc_up - uc) / dy)
    dv_dx = np.where(uc > 0.0, (vc - vc_left) / dx, (vc_right - vc) / dx)
    dv_dy = np.where(vc > 0.0, (vc - vc_down) / dy, (vc_up - vc) / dy)
    conv_uc = uc * du_dx + vc * du_dy
    conv_vc = uc * dv_dx + vc * dv_dy
    # map cell-centered convective terms back to faces by averaging
    conv_u = np.zeros_like(u)
    # interior faces
    conv_u[1:-1, :] = 0.5 * (conv_uc[0:-1, :] + conv_uc[1:, :])
    conv_u[0, :] = conv_uc[0, :]
    conv_u[-1, :] = conv_uc[-1, :]
    conv_v = np.zeros_like(v)
    conv_v[:, 1:-1] = 0.5 * (conv_vc[:, 0:-1] + conv_vc[:, 1:])
    conv_v[:, 0] = conv_vc[:, 0]
    conv_v[:, -1] = conv_vc[:, -1]
    return np.asfortranarray(conv_u), np.asfortranarray(conv_v)


def laplacian_u(u, dx, dy):
    Nx_u, Ny = u.shape
    Lu = np.zeros_like(u)
    # x-second derivative vectorized
    # interior i=1..Nx_u-2
    Lu[1:-1, :] += (u[2:, :] - 2.0 * u[1:-1, :] + u[0:-2, :]) / dx**2
    # rightmost face i = Nx_u-1 -> Neumann mirrored ghost: d2 = 2*(u[i-1]-u[i])/dx^2
    Lu[-1, :] += 2.0 * (u[-2, :] - u[-1, :]) / dx**2
    # leftmost face i = 0 -> Dirichlet inlet, approximate second derivative using one-sided (u[1]-u[0]) / dx^2
    Lu[0, :] += (u[1, :] - u[0, :]) / dx**2
    # y-second derivative with antisymmetric ghost for no-slip walls: ghost bottom = -u[:,0], ghost top = -u[:,-1]
    # interior j=1..Ny-2
    Lu[:, 1:-1] += (u[:, 2:] - 2.0 * u[:, 1:-1] + u[:, 0:-2]) / dy**2
    # bottom j=0: ghost = -u[:,0]
    Lu[:, 0] += (u[:, 1] - 2.0 * u[:, 0] + (-u[:, 0])) / dy**2
    # top j=Ny-1
    Lu[:, -1] += ((-u[:, -1]) - 2.0 * u[:, -1] + u[:, -2]) / dy**2
    return np.asfortranarray(Lu)


def laplacian_v(v, dx, dy):
    Nx, Ny_v = v.shape
    Lv = np.zeros_like(v)
    # y-second derivative interior j=1..Ny_v-2
    Lv[:, 1:-1] += (v[:, 2:] - 2.0 * v[:, 1:-1] + v[:, 0:-2]) / dy**2
    # bottom j=0 Dirichlet v=0 -> one-sided reuse
    Lv[:, 0] += (v[:, 1] - v[:, 0]) / dy**2
    # top j=Ny_v-1 Dirichlet v=0
    Lv[:, -1] += ( - v[:, -1] - v[:, -1]) / dy**2  # (ghost= -v[:, -1])
    # x-second derivative interior i=1..Nx-2
    Lv[1:-1, :] += (v[2:, :] - 2.0 * v[1:-1, :] + v[0:-2, :]) / dx**2
    # leftmost i=0 Dirichlet v=0 -> one-sided
    Lv[0, :] += (v[1, :] - v[0, :]) / dx**2
    # rightmost i=Nx-1 Neumann mirrored
    Lv[-1, :] += 2.0 * (v[-2, :] - v[-1, :]) / dx**2
    return np.asfortranarray(Lv)

# ----------------------- Pressure Poisson assembly -----------------------

def assemble_poisson(Nx, Ny, dx, dy):
    # Assemble standard discrete Laplacian L (negative definite): diag = -sum, off-diagonals = +1/dx^2 or +1/dy^2
    N = Nx * Ny
    data = []
    rows = []
    cols = []
    def idx(i, j):
        return i + j * Nx  # Fortran-order mapping
    for j in range(Ny):
        for i in range(Nx):
            row = idx(i, j)
            diag = 0.0
            # left neighbor
            if i - 1 >= 0:
                rows.append(row); cols.append(idx(i-1, j)); data.append(1.0 / dx**2); diag -= 1.0 / dx**2
            else:
                # inlet Neumann (skip)
                pass
            # right neighbor
            if i + 1 <= Nx - 1:
                rows.append(row); cols.append(idx(i+1, j)); data.append(1.0 / dx**2); diag -= 1.0 / dx**2
            else:
                # outlet Dirichlet p=0 -> implement via diag term so that solution enforces p=0 at boundary cell center beyond domain
                diag -= 1.0 / dx**2
            # down neighbor
            if j - 1 >= 0:
                rows.append(row); cols.append(idx(i, j-1)); data.append(1.0 / dy**2); diag -= 1.0 / dy**2
            else:
                # bottom Neumann skip
                pass
            # up neighbor
            if j + 1 <= Ny - 1:
                rows.append(row); cols.append(idx(i, j+1)); data.append(1.0 / dy**2); diag -= 1.0 / dy**2
            else:
                # top Neumann skip
                pass
            # diagonal
            rows.append(row); cols.append(row); data.append(diag)
    A = sp.csr_matrix((data, (rows, cols)), shape=(N, N))
    return A

# ----------------------- Divergence and projection -----------------------

def divergence(u, v, dx, dy):
    div = (u[1:, :] - u[0:-1, :]) / dx + (v[:, 1:] - v[:, 0:-1]) / dy
    return div


def project(u, v, dx, dy, dt, poisson_solver, Nx, Ny):
    # Ensure array contiguity and shapes
    assert u.flags['F_CONTIGUOUS'] and v.flags['F_CONTIGUOUS']
    div = divergence(u, v, dx, dy)
    rhs = (div / dt).ravel(order='F')  # using assembled L, so rhs = div/dt
    # solve A * phi = rhs
    phi_vec = poisson_solver(rhs)
    phi = phi_vec.reshape((Nx, Ny), order='F')
    # Correct u
    gradx = np.zeros_like(u)
    gradx[1:-1, :] = (phi[1:, :] - phi[0:-1, :]) / dx
    # rightmost face uses outlet Dirichlet p=0
    gradx[-1, :] = (0.0 - phi[-1, :]) / dx
    # leftmost inlet face not corrected (Dirichlet)
    u_corr = u - dt * gradx
    # Correct v
    grady = np.zeros_like(v)
    grady[:, 1:-1] = (phi[:, 1:] - phi[:, 0:-1]) / dy
    v_corr = v - dt * grady
    p_inc = phi.copy()
    return np.asfortranarray(u_corr), np.asfortranarray(v_corr), np.asfortranarray(p_inc)

# ----------------------- Main solver -----------------------

def main():
    # Grid and time parameters
    Nx = 80
    Ny = 40
    Lx = 2.0
    Ly = 1.0
    dx, dy, Xu, Yu, Xv, Yv, Xp, Yp = create_grid(Nx, Ny, Lx, Ly)
    # CFL-based dt
    maxU_est = 3.0
    CFL = 0.4
    dt_conv = CFL * min(dx / maxU_est, dy / maxU_est)
    dt_diff = 0.25 * min(dx*dx, dy*dy) / visc
    dt = min(dt_conv, dt_diff)
    dt = min(dt, 0.008)
    t_final = 1.0
    Nt = int(np.ceil(t_final / dt))
    dt = t_final / Nt
    print(f"Grid: Nx={Nx}, Ny={Ny}, dx={dx:.4e}, dy={dy:.4e}, dt={dt:.4e}, steps={Nt}")

    # Initialize fields (Fortran-order)
    u = np.asfortranarray(np.zeros((Nx+1, Ny)))
    v = np.asfortranarray(np.zeros((Nx, Ny+1)))
    p = np.asfortranarray(np.zeros((Nx, Ny)))

    # Preassemble Poisson matrix and factorize
    A = assemble_poisson(Nx, Ny, dx, dy)
    A_csc = A.tocsc()
    lu = spla.splu(A_csc)
    def poisson_solver(rhs):
        return lu.solve(rhs)

    # Time stepping SSP-RK3 with projection each substage
    t = 0.0
    print_interval = max(1, Nt // 10)
    start_time = time.time()

    apply_velocity_bcs(u, v, Xu, Yu, Xv, Yv, t)

    for step in range(1, Nt+1):
        u0 = u.copy(order='F')
        v0 = v.copy(order='F')
        p_accum = np.zeros_like(p)
        u_stage = u.copy(order='F')
        v_stage = v.copy(order='F')
        # 3-stage SSP-RK3
        for stage in range(3):
            if stage == 0:
                a = 1.0; b = 0.0; t_sub = t
            elif stage == 1:
                a = 3.0/4.0; b = 1.0/4.0; t_sub = t + dt
            else:
                a = 1.0/3.0; b = 2.0/3.0; t_sub = t + dt
            # enforce BCs consistently before operators
            apply_velocity_bcs(u_stage, v_stage, Xu, Yu, Xv, Yv, t_sub)
            conv_u, conv_v = compute_convective_terms(u_stage, v_stage, dx, dy)
            Lu = laplacian_u(u_stage, dx, dy)
            Lv = laplacian_v(v_stage, dx, dy)
            fx = np.zeros_like(u_stage)
            fy = -np.sin(np.pi * Xv) * np.sin(np.pi * Yv) * np.sin(np.pi * t_sub)
            # explicit update
            u_tent = u_stage + dt * ( - conv_u + visc * Lu + fx )
            v_tent = v_stage + dt * ( - conv_v + visc * Lv + fy )
            # apply BCs to tentative
            apply_velocity_bcs(u_tent, v_tent, Xu, Yu, Xv, Yv, t_sub)
            # safety checks before projection
            if not np.isfinite(u_tent).all() or not np.isfinite(v_tent).all():
                raise RuntimeError(f"Non-finite values detected before projection at step {step}, stage {stage}")
            if np.max(np.abs(u_tent)) > SAFETY_MAX or np.max(np.abs(v_tent)) > SAFETY_MAX:
                raise RuntimeError(f"Values exceeded safety threshold before projection at step {step}, stage {stage}")
            # project
            u_proj, v_proj, phi = project(u_tent.copy(order='F'), v_tent.copy(order='F'), dx, dy, dt, poisson_solver, Nx, Ny)
            apply_velocity_bcs(u_proj, v_proj, Xu, Yu, Xv, Yv, t_sub)
            # safety checks after projection
            if not np.isfinite(u_proj).all() or not np.isfinite(v_proj).all():
                raise RuntimeError(f"Non-finite values detected after projection at step {step}, stage {stage}")
            if np.max(np.abs(u_proj)) > SAFETY_MAX or np.max(np.abs(v_proj)) > SAFETY_MAX:
                raise RuntimeError(f"Values exceeded safety threshold after projection at step {step}, stage {stage}")
            # SSP combination
            if stage == 0:
                u_stage = u_proj.copy(order='F')
                v_stage = v_proj.copy(order='F')
                p_accum = phi.copy(order='F')
            else:
                u_stage = a * u0 + b * u_proj
                v_stage = a * v0 + b * v_proj
                p_accum = a * p_accum + b * phi
        # accept stage
        u[:] = u_stage
        v[:] = v_stage
        p[:] = p_accum
        t += dt
        # diagnostics
        if step % print_interval == 0 or step == Nt:
            div_max = np.max(np.abs(divergence(u, v, dx, dy)))
            print(f"Step {step}/{Nt}, t={t:.4f}, max|u|={np.max(np.abs(u)):.6f}, max|v|={np.max(np.abs(v)):.6f}, max|div|={div_max:.3e}")
    total_time = time.time() - start_time
    print(f"Simulation completed in {total_time:.2f} s")

    # Postprocess
    u_c, v_c = compute_cell_center_velocities(u, v)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)
    levels = 50
    cf0 = axes[0].contourf(Xp, Yp, u_c, levels=levels, cmap='RdBu_r')
    axes[0].set_title('u (cell-centered) at t=1.0')
    axes[0].set_xlabel('x'); axes[0].set_ylabel('y')
    fig.colorbar(cf0, ax=axes[0])
    cf1 = axes[1].contourf(Xp, Yp, v_c, levels=levels, cmap='RdBu_r')
    axes[1].set_title('v (cell-centered) at t=1.0')
    axes[1].set_xlabel('x'); axes[1].set_ylabel('y')
    fig.colorbar(cf1, ax=axes[1])
    cf2 = axes[2].contourf(Xp, Yp, p, levels=levels, cmap='RdBu_r')
    axes[2].set_title('p (cell-centered) at t=1.0')
    axes[2].set_xlabel('x'); axes[2].set_ylabel('y')
    fig.colorbar(cf2, ax=axes[2])
    plt.suptitle('2D Incompressible Navier-Stokes (MAC, SSP-RK3 projection)')
    plt.savefig('ns_mac_results.png', dpi=150)

if __name__ == '__main__':
    main()



#### Script block3:

```python
#!/usr/bin/env python
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time

# Global viscous coefficient
Re = 100.0
visc = 1.0 / Re

# Safety threshold for runaway values
SAFETY_MAX = 1e6

# ----------------------- Utility / Grid functions -----------------------
def create_grid(Nx, Ny, Lx=2.0, Ly=1.0):
    dx = Lx / float(Nx)
    dy = Ly / float(Ny)
    xu = np.linspace(0.0, Lx, Nx + 1)
    yu = (np.arange(Ny) + 0.5) * dy
    xv = (np.arange(Nx) + 0.5) * dx
    yv = np.linspace(0.0, Ly, Ny + 1)
    xp = (np.arange(Nx) + 0.5) * dx
    yp = (np.arange(Ny) + 0.5) * dy
    Xu, Yu = np.meshgrid(xu, yu, indexing='ij')  # (Nx+1, Ny)
    Xv, Yv = np.meshgrid(xv, yv, indexing='ij')  # (Nx, Ny+1)
    Xp, Yp = np.meshgrid(xp, yp, indexing='ij')  # (Nx, Ny)
    return dx, dy, np.asfortranarray(Xu), np.asfortranarray(Yu), np.asfortranarray(Xv), np.asfortranarray(Yv), np.asfortranarray(Xp), np.asfortranarray(Yp)

# ----------------------- BC and field helpers -----------------------
def inlet_u_profile(y, t):
    return np.sin(np.pi * y) * (np.sin(np.pi * t) + np.sin(3 * np.pi * t) + np.sin(5 * np.pi * t))


def apply_velocity_bcs(u, v, Xu, Yu, Xv, Yv, t):
    # u: (Nx+1, Ny), v: (Nx, Ny+1)
    # Enforce inlet u (Dirichlet) at i=0
    u[0, :] = inlet_u_profile(Yu[0, :], t)
    # Enforce no-slip at horizontal walls (approximate by setting tangential u at adjacent faces to zero)
    u[:, 0] = 0.0
    u[:, -1] = 0.0
    # Outlet (rightmost face) Neumann for u: approximate by copying interior neighbor
    u[-1, :] = u[-2, :]
    # v: enforce inlet v ~ 0 at leftmost face column i=0
    v[0, :] = 0.0
    # bottom and top walls: v at j=0 and j=Ny are exactly on the walls -> Dirichlet zero
    v[:, 0] = 0.0
    v[:, -1] = 0.0
    # outlet Neumann for v (right boundary): copy neighbor
    v[-1, :] = v[-2, :]
    return

# Build Dirichlet masks (which faces are prescribed and thus NOT correctable by projection)
def build_dirichlet_masks(Nx, Ny):
    # u: shape (Nx+1, Ny); v: shape (Nx, Ny+1); p: shape (Nx, Ny)
    u_dir = np.zeros((Nx + 1, Ny), dtype=bool)
    v_dir = np.zeros((Nx, Ny + 1), dtype=bool)
    p_dir = np.zeros((Nx, Ny), dtype=bool)
    # u Dirichlet faces (inlet and walls): i=0 inlet, j=0 and j=Ny-1 walls
    u_dir[0, :] = True
    u_dir[:, 0] = True
    u_dir[:, -1] = True
    # v Dirichlet faces (inlet vertical faces and horizontal walls): v[0,:] inlet, v[:,0], v[:,-1]
    v_dir[:, 0] = True
    v_dir[:, -1] = True
    # pressure Dirichlet at outlet cell centers: i = Nx-1 (rightmost column of cell centers)
    p_dir[-1, :] = True
    # correctable masks are the logical NOT of dirichlet masks
    u_correctable = ~u_dir
    v_correctable = ~v_dir
    return u_dir, v_dir, u_correctable, v_correctable, p_dir

# ----------------------- Discrete operators -----------------------

def compute_cell_center_velocities(u, v):
    u_c = 0.5 * (u[0:-1, :] + u[1:, :])
    v_c = 0.5 * (v[:, 0:-1] + v[:, 1:])
    return np.asfortranarray(u_c), np.asfortranarray(v_c)


def compute_convective_terms(u, v, dx, dy):
    Nx = u.shape[0] - 1
    Ny = u.shape[1]
    uc, vc = compute_cell_center_velocities(u, v)  # (Nx,Ny)
    # Prepare shifted arrays with one-sided fill
    uc_left = np.empty_like(uc)
    uc_right = np.empty_like(uc)
    uc_left[1:, :] = uc[0:-1, :]
    uc_left[0, :] = uc[0, :]
    uc_right[0:-1, :] = uc[1:, :]
    uc_right[-1, :] = uc[-1, :]
    uc_down = np.empty_like(uc)
    uc_up = np.empty_like(uc)
    uc_down[:, 1:] = uc[:, 0:-1]
    uc_down[:, 0] = uc[:, 0]
    uc_up[:, 0:-1] = uc[:, 1:]
    uc_up[:, -1] = uc[:, -1]
    vc_left = np.empty_like(vc)
    vc_right = np.empty_like(vc)
    vc_left[1:, :] = vc[0:-1, :]
    vc_left[0, :] = vc[0, :]
    vc_right[0:-1, :] = vc[1:, :]
    vc_right[-1, :] = vc[-1, :]
    vc_down = np.empty_like(vc)
    vc_up = np.empty_like(vc)
    vc_down[:, 1:] = vc[:, 0:-1]
    vc_down[:, 0] = vc[:, 0]
    vc_up[:, 0:-1] = vc[:, 1:]
    vc_up[:, -1] = vc[:, -1]
    du_dx = np.where(uc > 0.0, (uc - uc_left) / dx, (uc_right - uc) / dx)
    du_dy = np.where(vc > 0.0, (uc - uc_down) / dy, (uc_up - uc) / dy)
    dv_dx = np.where(uc > 0.0, (vc - vc_left) / dx, (vc_right - vc) / dx)
    dv_dy = np.where(vc > 0.0, (vc - vc_down) / dy, (vc_up - vc) / dy)
    conv_uc = uc * du_dx + vc * du_dy
    conv_vc = uc * dv_dx + vc * dv_dy
    conv_u = np.zeros_like(u)
    conv_u[1:-1, :] = 0.5 * (conv_uc[0:-1, :] + conv_uc[1:, :])
    conv_u[0, :] = conv_uc[0, :]
    conv_u[-1, :] = conv_uc[-1, :]
    conv_v = np.zeros_like(v)
    conv_v[:, 1:-1] = 0.5 * (conv_vc[:, 0:-1] + conv_vc[:, 1:])
    conv_v[:, 0] = conv_vc[:, 0]
    conv_v[:, -1] = conv_vc[:, -1]
    return np.asfortranarray(conv_u), np.asfortranarray(conv_v)


def laplacian_u(u, dx, dy):
    Nx_u, Ny = u.shape
    Lu = np.zeros_like(u)
    Lu[1:-1, :] += (u[2:, :] - 2.0 * u[1:-1, :] + u[0:-2, :]) / dx ** 2
    Lu[-1, :] += 2.0 * (u[-2, :] - u[-1, :]) / dx ** 2
    Lu[0, :] += (u[1, :] - u[0, :]) / dx ** 2
    Lu[:, 1:-1] += (u[:, 2:] - 2.0 * u[:, 1:-1] + u[:, 0:-2]) / dy ** 2
    Lu[:, 0] += (u[:, 1] - 2.0 * u[:, 0] + (-u[:, 0])) / dy ** 2
    Lu[:, -1] += ((-u[:, -1]) - 2.0 * u[:, -1] + u[:, -2]) / dy ** 2
    return np.asfortranarray(Lu)


def laplacian_v(v, dx, dy):
    Nx, Ny_v = v.shape
    Lv = np.zeros_like(v)
    Lv[:, 1:-1] += (v[:, 2:] - 2.0 * v[:, 1:-1] + v[:, 0:-2]) / dy ** 2
    Lv[:, 0] += (v[:, 1] - v[:, 0]) / dy ** 2
    Lv[:, -1] += (-v[:, -1] - v[:, -1]) / dy ** 2
    Lv[1:-1, :] += (v[2:, :] - 2.0 * v[1:-1, :] + v[0:-2, :]) / dx ** 2
    Lv[0, :] += (v[1, :] - v[0, :]) / dx ** 2
    Lv[-1, :] += 2.0 * (v[-2, :] - v[-1, :]) / dx ** 2
    return np.asfortranarray(Lv)

# ----------------------- Pressure Poisson assembly (masked) -----------------------

def assemble_poisson_masked(Nx, Ny, dx, dy, u_correctable, v_correctable, p_dirichlet_mask):
    # Assemble discrete operator corresponding to div( mask * grad phi ) on cell centers
    N = Nx * Ny
    data = []
    rows = []
    cols = []
    def idx(i, j):
        return i + j * Nx  # Fortran-order mapping
    for j in range(Ny):
        for i in range(Nx):
            row = idx(i, j)
            if p_dirichlet_mask[i, j]:
                # enforce Dirichlet pressure: set row to identity
                rows.append(row); cols.append(row); data.append(1.0)
                continue
            diag = 0.0
            # left neighbor: face between (i-1,j) and (i,j) is u[i, j]
            if i - 1 >= 0:
                if u_correctable[i, j]:
                    rows.append(row); cols.append(idx(i - 1, j)); data.append(1.0 / dx ** 2); diag -= 1.0 / dx ** 2
                else:
                    # face not correctable -> no coupling across it
                    pass
            else:
                # domain boundary left: if face is correctable include; else ignore (Neumann-like handled by masks)
                if u_correctable[i, j]:
                    # left ghost handled by p_outside = p_inside -> contributes -1/dx^2 on diag
                    diag -= 1.0 / dx ** 2
            # right neighbor: face between (i,j) and (i+1,j) is u[i+1, j]
            if i + 1 <= Nx - 1:
                if u_correctable[i + 1, j]:
                    rows.append(row); cols.append(idx(i + 1, j)); data.append(1.0 / dx ** 2); diag -= 1.0 / dx ** 2
                else:
                    pass
            else:
                # right boundary (outside cell centers). If face correctable it would couple to outside p (treated via known Dirichlet p), else ignore
                if u_correctable[i + 1 - 0, j]:
                    # unreachable in typical layouts but keep for completeness
                    diag -= 1.0 / dx ** 2
            # down neighbor: face between (i,j-1) and (i,j) is v[i, j]
            if j - 1 >= 0:
                if v_correctable[i, j]:
                    rows.append(row); cols.append(idx(i, j - 1)); data.append(1.0 / dy ** 2); diag -= 1.0 / dy ** 2
                else:
                    pass
            else:
                if v_correctable[i, j]:
                    diag -= 1.0 / dy ** 2
            # up neighbor: face between (i,j) and (i,j+1) is v[i, j+1]
            if j + 1 <= Ny - 1:
                if v_correctable[i, j + 1]:
                    rows.append(row); cols.append(idx(i, j + 1)); data.append(1.0 / dy ** 2); diag -= 1.0 / dy ** 2
                else:
                    pass
            else:
                if v_correctable[i, j + 1 - 0]:
                    diag -= 1.0 / dy ** 2
            # diagonal
            rows.append(row); cols.append(row); data.append(diag)
    A = sp.csr_matrix((data, (rows, cols)), shape=(N, N))
    return A

# ----------------------- Divergence and projection -----------------------

def divergence(u, v, dx, dy):
    div = (u[1:, :] - u[0:-1, :]) / dx + (v[:, 1:] - v[:, 0:-1]) / dy
    return div


def project(u, v, dx, dy, dt, poisson_solver, Nx, Ny, u_correctable, v_correctable, p_dirichlet_mask):
    assert u.flags['F_CONTIGUOUS'] and v.flags['F_CONTIGUOUS']
    # full divergence from tentative velocity
    div_full = divergence(u, v, dx, dy)
    rhs = (div_full / dt).ravel(order='F')
    # enforce rhs entries where p is Dirichlet (set RHS to prescribed pressure value which is 0 here)
    if np.any(p_dirichlet_mask):
        rhs[p_dirichlet_mask.ravel(order='F')] = 0.0
    # solve masked Poisson
    phi_vec = poisson_solver(rhs)
    phi = phi_vec.reshape((Nx, Ny), order='F')
    # compute gradient on faces but only apply corrections on correctable faces
    gradx = np.zeros_like(u)
    # interior faces i=1..Nx-1
    gradx[1:-1, :] = (phi[1:, :] - phi[0:-1, :]) / dx
    # rightmost face uses p_outside = 0.0 (pressure Dirichlet at outlet cell centers ensures phi[-1,:] already 0 if enforced)
    gradx[-1, :] = (0.0 - phi[-1, :]) / dx
    # leftmost face grad is not used (Dirichlet inlet), gradx[0,:] stays zero
    # but we must not apply grad where face is not correctable
    u_corr = u.copy(order='F')
    mask_u = u_correctable
    u_corr[mask_u] = u_corr[mask_u] - dt * gradx[mask_u]
    # v corrections
    grady = np.zeros_like(v)
    grady[:, 1:-1] = (phi[:, 1:] - phi[:, 0:-1]) / dy
    # top/bottom: grady at top uses ghost p_outside=0? but v faces at boundaries are often Dirichlet (walls); we will not correct them if not correctable
    v_corr = v.copy(order='F')
    mask_v = v_correctable
    v_corr[mask_v] = v_corr[mask_v] - dt * grady[mask_v]
    return np.asfortranarray(u_corr), np.asfortranarray(v_corr), np.asfortranarray(phi)

# ----------------------- Main solver -----------------------

def main():
    # Grid and time parameters
    Nx = 80
    Ny = 40
    Lx = 2.0
    Ly = 1.0
    dx, dy, Xu, Yu, Xv, Yv, Xp, Yp = create_grid(Nx, Ny, Lx, Ly)
    # CFL-based dt
    maxU_est = 3.0
    CFL = 0.4
    dt_conv = CFL * min(dx / maxU_est, dy / maxU_est)
    dt_diff = 0.25 * min(dx * dx, dy * dy) / visc
    dt = min(dt_conv, dt_diff)
    dt = min(dt, 0.008)
    t_final = 1.0
    Nt = int(np.ceil(t_final / dt))
    dt = t_final / Nt
    print(f"Grid: Nx={Nx}, Ny={Ny}, dx={dx:.4e}, dy={dy:.4e}, dt={dt:.4e}, steps={Nt}")

    # Initialize fields (Fortran-order)
    u = np.asfortranarray(np.zeros((Nx + 1, Ny)))
    v = np.asfortranarray(np.zeros((Nx, Ny + 1)))
    p = np.asfortranarray(np.zeros((Nx, Ny)))

    # Build Dirichlet / correctable masks (static topology)
    u_dir, v_dir, u_correctable, v_correctable, p_dir = build_dirichlet_masks(Nx, Ny)

    # Preassemble masked Poisson matrix and factorize
    A_masked = assemble_poisson_masked(Nx, Ny, dx, dy, u_correctable, v_correctable, p_dir)
    A_csc = A_masked.tocsc()
    lu = spla.splu(A_csc)

    def poisson_solver(rhs):
        return lu.solve(rhs)

    # Time stepping SSP-RK3 with projection each substage
    t = 0.0
    print_interval = max(1, Nt // 10)
    start_time = time.time()

    apply_velocity_bcs(u, v, Xu, Yu, Xv, Yv, t)

    for step in range(1, Nt + 1):
        u0 = u.copy(order='F')
        v0 = v.copy(order='F')
        p_accum = np.zeros_like(p)
        u_stage = u.copy(order='F')
        v_stage = v.copy(order='F')
        # 3-stage SSP-RK3
        for stage in range(3):
            if stage == 0:
                a = 1.0; b = 0.0; t_sub = t
            elif stage == 1:
                a = 3.0 / 4.0; b = 1.0 / 4.0; t_sub = t + dt
            else:
                a = 1.0 / 3.0; b = 2.0 / 3.0; t_sub = t + dt
            # enforce BCs consistently before operators
            apply_velocity_bcs(u_stage, v_stage, Xu, Yu, Xv, Yv, t_sub)
            conv_u, conv_v = compute_convective_terms(u_stage, v_stage, dx, dy)
            Lu = laplacian_u(u_stage, dx, dy)
            Lv = laplacian_v(v_stage, dx, dy)
            fx = np.zeros_like(u_stage)
            fy = -np.sin(np.pi * Xv) * np.sin(np.pi * Yv) * np.sin(np.pi * t_sub)
            # explicit update
            u_tent = u_stage + dt * ( - conv_u + visc * Lu + fx )
            v_tent = v_stage + dt * ( - conv_v + visc * Lv + fy )
            # apply BCs to tentative (ensure prescribed faces carry prescribed values)
            apply_velocity_bcs(u_tent, v_tent, Xu, Yu, Xv, Yv, t_sub)
            # safety checks before projection
            if not np.isfinite(u_tent).all() or not np.isfinite(v_tent).all():
                raise RuntimeError(f"Non-finite values detected before projection at step {step}, stage {stage}")
            if np.max(np.abs(u_tent)) > SAFETY_MAX or np.max(np.abs(v_tent)) > SAFETY_MAX:
                raise RuntimeError(f"Values exceeded safety threshold before projection at step {step}, stage {stage}")
            # project using masked Poisson operator and only update correctable faces
            u_proj, v_proj, phi = project(u_tent.copy(order='F'), v_tent.copy(order='F'), dx, dy, dt, poisson_solver, Nx, Ny, u_correctable, v_correctable, p_dir)
            # re-apply BCs after projection to ensure exact Dirichlet values on prescribed faces
            apply_velocity_bcs(u_proj, v_proj, Xu, Yu, Xv, Yv, t_sub)
            # safety checks after projection
            if not np.isfinite(u_proj).all() or not np.isfinite(v_proj).all():
                raise RuntimeError(f"Non-finite values detected after projection at step {step}, stage {stage}")
            if np.max(np.abs(u_proj)) > SAFETY_MAX or np.max(np.abs(v_proj)) > SAFETY_MAX:
                raise RuntimeError(f"Values exceeded safety threshold after projection at step {step}, stage {stage}")
            # SSP combination
            if stage == 0:
                u_stage = u_proj.copy(order='F')
                v_stage = v_proj.copy(order='F')
                p_accum = phi.copy(order='F')
            else:
                u_stage = a * u0 + b * u_proj
                v_stage = a * v0 + b * v_proj
                p_accum = a * p_accum + b * phi
        # accept stage
        u[:] = u_stage
        v[:] = v_stage
        p[:] = p_accum
        t += dt
        # diagnostics
        if step % print_interval == 0 or step == Nt:
            div_max = np.max(np.abs(divergence(u, v, dx, dy)))
            div_l2 = np.linalg.norm(divergence(u, v, dx, dy).ravel(order='F'))
            print(f"Step {step}/{Nt}, t={t:.4f}, max|u|={np.max(np.abs(u)):.6f}, max|v|={np.max(np.abs(v)):.6f}, max|div|={div_max:.3e}, L2|div|={div_l2:.3e}")
    total_time = time.time() - start_time
    print(f"Simulation completed in {total_time:.2f} s")

    # Postprocess
    u_c, v_c = compute_cell_center_velocities(u, v)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)
    levels = 50
    cf0 = axes[0].contourf(Xp, Yp, u_c, levels=levels, cmap='RdBu_r')
    axes[0].set_title('u (cell-centered) at t=1.0')
    axes[0].set_xlabel('x'); axes[0].set_ylabel('y')
    fig.colorbar(cf0, ax=axes[0])
    cf1 = axes[1].contourf(Xp, Yp, v_c, levels=levels, cmap='RdBu_r')
    axes[1].set_title('v (cell-centered) at t=1.0')
    axes[1].set_xlabel('x'); axes[1].set_ylabel('y')
    fig.colorbar(cf1, ax=axes[1])
    cf2 = axes[2].contourf(Xp, Yp, p, levels=levels, cmap='RdBu_r')
    axes[2].set_title('p (cell-centered) at t=1.0')
    axes[2].set_xlabel('x'); axes[2].set_ylabel('y')
    fig.colorbar(cf2, ax=axes[2])
    plt.suptitle('2D Incompressible Navier-Stokes (MAC, SSP-RK3 projection - masked Poisson)')
    plt.savefig('ns_mac_results.png', dpi=150)

if __name__ == '__main__':
    main()
```



### runtime_outputs

#### Output block1

Grid: Nx=80, Ny=40, dx=2.5000e-02, dy=2.5000e-02, dt=3.3333e-03, steps=300
<string>:126: RuntimeWarning: overflow encountered in multiply
<string>:126: RuntimeWarning: invalid value encountered in add
<string>:127: RuntimeWarning: overflow encountered in multiply
<string>:127: RuntimeWarning: invalid value encountered in add
<string>:132: RuntimeWarning: invalid value encountered in add
<string>:139: RuntimeWarning: invalid value encountered in add
<string>:248: RuntimeWarning: invalid value encountered in subtract
Step 30/300, t=0.1000, max|u|=nan, max|v|=nan, max|div|=nan
Step 60/300, t=0.2000, max|u|=nan, max|v|=nan, max|div|=nan
Step 90/300, t=0.3000, max|u|=nan, max|v|=nan, max|div|=nan
Step 120/300, t=0.4000, max|u|=nan, max|v|=nan, max|div|=nan
Step 150/300, t=0.5000, max|u|=nan, max|v|=nan, max|div|=nan
Step 180/300, t=0.6000, max|u|=nan, max|v|=nan, max|div|=nan
Step 210/300, t=0.7000, max|u|=nan, max|v|=nan, max|div|=nan
Step 240/300, t=0.8000, max|u|=nan, max|v|=nan, max|div|=nan
Step 270/300, t=0.9000, max|u|=nan, max|v|=nan, max|div|=nan
Step 300/300, t=1.0000, max|u|=nan, max|v|=nan, max|div|=nan
Simulation completed in 6.03 s



#### Output block2

Grid: Nx=80, Ny=40, dx=2.5000e-02, dy=2.5000e-02, dt=3.3333e-03, steps=300
Step 30/300, t=0.1000, max|u|=2.100426, max|v|=0.414427, max|div|=2.036e-01
Step 60/300, t=0.2000, max|u|=1.620120, max|v|=0.247832, max|div|=5.007e-01
Step 90/300, t=0.3000, max|u|=0.477387, max|v|=0.394262, max|div|=1.972e-01
Step 120/300, t=0.4000, max|u|=0.589887, max|v|=0.279896, max|div|=2.455e-01
Step 150/300, t=0.5000, max|u|=1.040095, max|v|=0.187189, max|div|=1.288e-01
Step 180/300, t=0.6000, max|u|=0.581553, max|v|=0.189404, max|div|=2.393e-01
Step 210/300, t=0.7000, max|u|=0.357434, max|v|=0.211153, max|div|=2.093e-01
Step 240/300, t=0.8000, max|u|=1.510447, max|v|=0.261031, max|div|=3.085e-01
Step 270/300, t=0.9000, max|u|=2.130208, max|v|=0.251786, max|div|=4.695e-01
Step 300/300, t=1.0000, max|u|=0.622437, max|v|=0.462823, max|div|=6.295e-01
Simulation completed in 0.48 s



#### Output block3

Grid: Nx=80, Ny=40, dx=2.5000e-02, dy=2.5000e-02, dt=3.3333e-03, steps=300
Step 30/300, t=0.1000, max|u|=2.100426, max|v|=0.413012, max|div|=2.193e-01, L2|div|=5.731e-01
Step 60/300, t=0.2000, max|u|=1.619338, max|v|=0.248429, max|div|=5.145e-01, L2|div|=9.149e-01
Step 90/300, t=0.3000, max|u|=0.476440, max|v|=0.393666, max|div|=1.913e-01, L2|div|=6.049e-01
Step 120/300, t=0.4000, max|u|=0.590489, max|v|=0.280990, max|div|=2.491e-01, L2|div|=7.145e-01
Step 150/300, t=0.5000, max|u|=1.040102, max|v|=0.186843, max|div|=2.075e-01, L2|div|=9.061e-01
Step 180/300, t=0.6000, max|u|=0.580382, max|v|=0.189048, max|div|=3.167e-01, L2|div|=1.422e+00
Step 210/300, t=0.7000, max|u|=0.357880, max|v|=0.211947, max|div|=3.817e-01, L2|div|=1.607e+00
Step 240/300, t=0.8000, max|u|=1.510447, max|v|=0.259177, max|div|=4.884e-01, L2|div|=2.137e+00
Step 270/300, t=0.9000, max|u|=2.130149, max|v|=0.252389, max|div|=7.014e-01, L2|div|=2.914e+00
Step 300/300, t=1.0000, max|u|=0.621769, max|v|=0.459681, max|div|=7.938e-01, L2|div|=3.570e+00
Simulation completed in 0.52 s




### review_decision
revise

### review_comments
	Current Stage [A/3]
1. MAJOR

Category: accuracy
Issue: Wrong sign convention between assembled Poisson matrix and RHS in projection: A was assembled as -Laplace (A = -L), but the code solves A * phi = +div/dt instead of A * phi = -div/dt (or assembling L instead). This produces an incorrect pressure increment and can drive the velocity update unstable, causing blow-up / NaNs.

Recommendation: Fix the sign mismatch. Either (a) change assemble_poisson to assemble the standard discrete Laplacian L (neighbors +1/dx^2 on off-diagonals, diagonal negative) and keep rhs = div/dt, or (b) keep A = -L but set rhs = -div/dt before solving. After fixing, re-run and confirm that div ~ O(dt) and no overflow/NaN occurs.


2. MAJOR

Category: runtime
Issue: Numerical overflow and NaNs observed early in the run (RuntimeWarnings), indicating an instability likely originating from the projection sign error and/or inconsistent boundary treatments.

Recommendation: After fixing the sign issue above, add safety checks and short-circuit early: monitor max(|u|,|v|,|div|) every sub-step and abort (with informative error) if values exceed a large threshold (e.g. 1e6). Run a very coarse test (small Nx,Ny) and very small dt to verify stable behavior before full production runs.


3. MAJOR

Category: structure
Issue: Inconsistent / implicit handling of velocity boundary conditions (especially no-slip at top/bottom for u faces and outlet Neumann) across the convective, Laplacian, and projection operators. Some face values are 'skipped' in Laplacian (left as zeros) or treated via ghost approximations, which can lead to inconsistencies and mass conservation errors.

Recommendation: Make the BC enforcement explicit and consistent for all operators: (1) clearly define ghost values for all velocity faces (u and v) so that Laplacian and convective discretizations see the same boundary treatment; (2) set u at wall-adjacent faces to enforce no-slip (use mirrored/antisymmetric ghost values where appropriate) before computing convective terms and Laplacians; (3) document and unit-test BC handling on small grids (e.g. analytic checks for constant fields).


4. MINOR

Category: accuracy
Issue: The convective discretization is implemented with nested Python loops and a cell-centered upwind approximation that is hard to reason about and easy to get index offsets wrong. It is also inefficient and likely a source of subtle mistakes (e.g. sign/offset errors near boundaries).

Recommendation: Vectorize the upwind discretization (use array operations) or at least isolate and unit-test the stencil on small examples. Verify that the convective fluxes satisfy expected sign and limit tests (e.g. advect a passive scalar). Ensure mapping between face/cell-centered velocities and derivatives is consistent with MAC staggering. Consider clearer naming: u_c and v_c -> uc, vc and document shapes explicitly.


5. MINOR

Category: structure
Issue: Mix of array memory orders and some arrays created without explicit order='F' may still work but risks subtle alignment bugs with ravel/reshape operations (code tries to use order='F' sometimes). Also ghost treatment and indexing comments are incomplete and ambiguous.

Recommendation: Be consistent: create major arrays (u,v,p) and coordinate grids with order='F' if you rely on Fortran-order ravel/reshape. Add explicit asserts on shapes and memory order where you call ravel/reshape with order='F', or always use numpy.asfortranarray when needed. Improve docstrings / comments on array shapes (which index is i/j) and unit-test reshape/ravel round-trips.



	Current Stage [B/3]
1. MAJOR

Category: accuracy
Issue: Projection / Poisson assembly inconsistent with velocity Dirichlet BCs leading to poor divergence control.

Recommendation: The projection step is not enforcing mass conservation — the reported max|div| up to O(1e-0) at the final time indicates the Poisson solve / velocity correction are inconsistent with the BC treatment. Fix the Poisson assembly and RHS so that contributions from prescribed (Dirichlet) face velocities are moved to the RHS, and impose pressure Dirichlet rows exactly (or eliminate corresponding DOFs). Concretely:
- When a face velocity is Dirichlet (e.g. inlet u[0,:]), its flux contribution (−u_prescribed/dx or similar sign depending on your divergence sign) must be removed from the unknown part and added to the solve RHS prior to solving A phi = rhs. Currently the code forms rhs = div/dt with the full discrete divergence (including prescribed faces) and does not account for known face values, so phi is incorrect and the correction does not produce divergence-free velocity.
- Implement Dirichlet p on the outlet cell-centers by changing the corresponding matrix rows to identity and setting rhs to the prescribed pressure value (0). Do not merely adjust the diagonal in a heuristic way (the current assemble_poisson's 'else: diag -= 1/dx**2' for the right boundary is not a correct enforcement of Dirichlet).


2. MAJOR

Category: structure
Issue: assemble_poisson handles boundaries incorrectly (inconsistent/skipped neighbors and heuristic diag modification).

Recommendation: Correct the discrete Laplacian boundary treatments in assemble_poisson. The current assembly mixes 'skip neighbor' semantics for Neumann and a heuristic diagonal tweak for Dirichlet, which is not a correct finite-difference imposition. Two robust approaches:
1) Use ghost-point elimination: explicitly construct the finite-difference stencil including ghost values and fold known boundary conditions into the RHS and diagonal entries. For Neumann, use one-sided differences to represent the derivative and update diagonal/RHS accordingly; for Dirichlet, set the matrix row to enforce phi = value.
2) Alternatively, build the Laplacian for interior unknowns only and handle boundary cells by elimination (remove Dirichlet unknowns). For mixed BCs this is more straightforward and less error-prone.
After fixing assembly, recompute and test that projection reduces divergence to machine precision (or very small values) each step.


3. MAJOR

Category: accuracy
Issue: Velocity BC re-application after projection without corresponding RHS adjustment invalidates projection consistency.

Recommendation: Ensure consistency between the velocity BC enforcement and the projection correction: if you re-impose Dirichlet velocity after projection (apply_velocity_bcs called after project), you must still include that Dirichlet information in the Poisson RHS (see the first comment). Otherwise the projection solves for a phi incompatible with the final enforced velocity field, undoing divergence-free enforcement.
Practical fix: compute the divergence used for rhs using only the unknown face values (treat known faces as sources), or equivalently compute rhs = (div_interior − known_face_flux_contributions)/dt before solving A phi = rhs.


4. MINOR

Category: structure
Issue: Poisson rows for Dirichlet p are not set explicitly to identity; current heuristic may produce an incorrect or poorly conditioned matrix.

Recommendation: Improve the Poisson assembly to explicitly and clearly set rows for Dirichlet pressure nodes (e.g. set row to identity). This both ensures a nonsingular well-conditioned system and makes LU/iterative solves more robust. For performance, consider assembling the positive-definite -A and use Cholesky (scipy.sparse.linalg.splu is fine but LU is overkill for symmetric positive definite problems).


5. MINOR

Category: structure
Issue: Code quality / diagnostics and clarity around discrete ghost/BC implementations can be improved.

Recommendation: A few other improvements / checks that will help robustness and reproducibility:
- After fixing the projection, report a divergence norm (L2 and max) immediately after projection to verify enforcement each step. Use small prints or asserts in debug mode.
- The laplacian_{u,v} stencils use heuristic ghost/antisymmetric assignments; verify these produce the intended order of accuracy and consistency with boundary conditions (it may be fine, but cross-check using manufactured solution tests). Consider documenting the discrete ghost-point formulas used.
- Consider making the Poisson matrix assembly and boundary handling modular and well-documented so it is straightforward to unit-test (e.g. test that a constant pressure field with zero RHS returns zero phi).



	Current Stage [C/3]
1. MAJOR

Category: accuracy
Issue: Projection / masked Poisson formulation is inconsistent: the assembled 'masked' matrix removes coupling across non-correctable faces but the RHS uses the full divergence (div_full/dt). This does not move the known fluxes from Dirichlet faces into the RHS and therefore the linear system cannot in general enforce continuity. The result is a substantial divergence residual (max|div| ~ O(1)) at the end of the run.

Recommendation: Reformulate the masked Poisson correctly. Two robust approaches: (A) Assemble the full discrete divergence-of-gradient operator (standard cell-centered Laplacian) and treat Dirichlet face velocities by moving their contribution to the RHS: when a face is Dirichlet, subtract (face_value / dx or dy) term from RHS so the unknown only corrects the remaining (correctable) faces; or (B) keep a masked operator but explicitly include the known flux contribution of non-correctable faces into the RHS (i.e. the RHS should be div_correctable/dt = (div_full - div_from_prescribed_faces)/dt). Concretely: for each cell, compute fluxes across all four faces but for faces that are prescribed, evaluate their contribution and transfer it to RHS; assemble A as if all faces were unknown. This will drive post-projection divergence down to numerical tolerance. Add a diagnostic that computes the Poisson residual (A*phi - RHS) and the divergence before/after projection to confirm correctness.


2. MAJOR

Category: accuracy
Issue: Inconsistent and likely incorrect finite-difference boundary treatments in laplacian_u / laplacian_v. Several lines use odd ghost-value substitutions (e.g. Lu[:,0] uses (-u[:,0]) as ghost, Lv[:, -1] += (-v[:, -1] - v[:, -1]) / dy ** 2) that appear wrong and will produce incorrect viscous terms and break symmetry of the operator.

Recommendation: Replace the ad-hoc boundary formulas with clear, documented implementations that match the physical BCs. For Dirichlet faces (u or v prescribed) impose ghost values by reflection for homogeneous Neumann (ghost = interior) or odd reflection for homogeneous Dirichlet (ghost = -interior) depending on the variable and boundary type. Prefer writing explicit loops or vectorized formulas that handle each boundary side with comments stating the boundary condition used. Verify the Laplacian by checking matrix-free consistency: compare finite-difference Laplacian to sparse matrix application for a random field and ensure symmetry where appropriate.


3. MAJOR

Category: structure
Issue: Boundary conditions and mask definitions are inconsistent in places: e.g., apply_velocity_bcs enforces u[-1,:] = u[-2,:] (outlet Neumann for u) while build_dirichlet_masks marks some top/bottom faces as Dirichlet but not the outlet faces. The mask and the BC enforcement must be consistent (faces you treat as Dirichlet should be enforced and removed from correction; faces you want Neumann should be treated accordingly).

Recommendation: Make a single authoritative specification of which faces are Dirichlet/Neumann and ensure both apply_velocity_bcs and build_dirichlet_masks follow it. If the outlet is Neumann for velocities, do NOT set p Dirichlet at the same location without treating the coupling correctly; if you want pressure Dirichlet on the outlet (p=0) keep pressure Dirichlet but then ensure velocity faces there are considered correctable (or include their prescribed Neumann flux into RHS). Add unit tests for simple steady cases (e.g. constant inlet) to ensure BCs and masks produce expected flux balance.


4. MAJOR

Category: accuracy
Issue: Sign convention / operator definiteness is unclear in assemble_poisson_masked: diagonal entries are created as negative of neighbor sums, producing an operator that is negative-definite (diagonal negative, off-diagonals positive). This confuses interpretation and makes it hard to use SPD solvers or reason about conditioning.

Recommendation: Adopt the standard convention consistently: either assemble A = -Laplacian (so A is symmetric positive-definite and you can use efficient SPD solvers) and set RHS = -div/dt, or assemble A = Laplacian and use matching RHS sign. Document the convention at the top of the assembly routine. This also helps diagnosing residuals and using iterative solvers (CG/PCG) if desired. Ensure the assembled matrix is symmetric (check A - A.T norm) and that row-sums and boundary rows behave as intended.


5. MINOR

Category: accuracy
Issue: Diagnostics indicate divergence remains large (max|div| up to ~0.8), and no residual checks are performed for the Poisson solve. Without these checks it's difficult to know whether the linear solve and projection are working as intended.

Recommendation: Add diagnostics: compute the Poisson residual r = A*phi - rhs and report norm(r) occasionally; compute divergence before and after projection and assert it is significantly reduced (e.g., L2 reduction by many orders). If you intentionally allow nonzero divergence near prescribed faces, compute and report divergence only over the subdomain of correctable faces or cells so it's clear what's controlled. This will make debugging much easier.





### iteration_count
3

### rev_image_description
	Current Stage [A/3]


	Current Stage [B/3]


	Current Stage [C/3]




