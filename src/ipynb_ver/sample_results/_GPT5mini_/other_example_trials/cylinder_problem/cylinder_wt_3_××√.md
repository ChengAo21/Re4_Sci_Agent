### 🛌🏻 Re4gent with Multi-Modal Review

#### Configs 🏷️:
MaxCount: 2, MultiModalReview: False, AlternativeNum: 2, MaxComments: 5, ConsModule_ON: True, CUT_OUTEXT: 2000

### prob_todo

The governing PDEs for 2D incompressible flow around a cylinder are given by:
\begin{cases}
\frac{\partial u }{\partial t} + u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} + \frac{\partial p}{\partial x} - \frac{1}{\mathrm{Re}} \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right) = 0, (x,y) \in \Omega,\\
\frac{\partial v }{\partial t} + u \frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} + \frac{\partial p}{\partial y} - \frac{1}{\mathrm{Re}} \left( \frac{\partial^2 v}{\partial x^2} + \frac{\partial^2 v}{\partial y^2} \right) = 0, (x,y) \in \Omega,\\
\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0, (x,y) \in \Omega,\\
\end{cases}

The computational domain \Omega is a rectangle with:
\begin{cases}
Cylinder diameter d = 1.0 (Normalized length scale)\
Distance from cylinder center to left boundary: x_1 = 3.0\
Distance from cylinder center to right boundary: x_2 = 8.0\
Distance from cylinder center to top/bottom boundaries: y_1 = y_2 = 3.0\
Total domain size: (x_1 + x_2) \times (y_1 + y_2)\
\end{cases}
The boundary conditions are:
\begin{cases}
Inlet \Gamma_{in} (left): (u,v) = (1, 0) constant inflow\
Walls \Gamma_{top} \cup \Gamma_{bottom}: (u,v) with \frac{\partial u}{\partial y}=0, v=0 (Free-slip condition)\
Cylinder Surface: No-slip condition\
Outlet \Gamma_{out}: zero pressure gradient\
\end{cases}
The Reynolds number Re = 100 (Based on diameter and inlet velocity). 

Please implement a stable and efficient method to solve this problem.
Implement reasonable acceleration strategies to reduce computational cost.
Define visc = 1.0 / Re at the GLOBAL scope (outside main), and strictly use visc instead of nu throughout the code to avoid confusion with u.
Ensure all constants are passed explicitly to functions or defined globally to prevent NameError.
Simulate until t=20.0. Plot the contour of the velocity magnitude at the final step in one figure using 'RdBu_r' colormap, and mark the cylinder in the plot. 
Just save figs do not use plt.show() in the code.

[HINTS]:
Never use the keyword argument 'tol' in any SciPy solver. Always use 'atol'.
Do not use 'np.trapz' as it is removed in NumPy 2.0.
Fix the pressure at one single node to remove the nullspace singularity.
Use the Immersed Boundary Method with volume penalization on a Staggered (MAC) Grid. Treat the penalization term implicitly in the momentum matrix diagonal to avoid restrictive time-steps.
Rely solely on the penalization term to enforce the no-slip condition on the cylinder surface.
Focus on generating runnable code with reasonable physical results, do not over-optimize.
Print concise progress information ONLY every 10% of total steps.


### expanded_prob
We must compute an incompressible 2D flow (Navier–Stokes) around a circular cylinder embedded in a rectangular domain. The unknown fields are the velocity components (u,v) and pressure p. The domain is rectangular with the cylinder of diameter d=1 centered at a known position; distances to the boundaries are given (left 3.0, right 8.0, top/bottom 3.0), so the computational box is (x1+x2) by (y1+y2). The Reynolds number is Re=100, so visc = 1.0/Re is small but finite. Boundary conditions: prescribed uniform inflow (u=1,v=0) at the left boundary, free-slip (v=0, ∂u/∂y=0) at top and bottom, outlet condition of zero pressure gradient at the right boundary, and a no-slip cylinder surface. The problem requires a stable and efficient time-accurate simulation until t = 20.0 and a final contour plot of velocity magnitude (colormap 'RdBu_r') marking the cylinder. Implementation constraints/hints to respect: use an immersed boundary method via volume penalization on a staggered (MAC) grid; treat the penalization term implicitly by adding its contribution to the momentum linear system diagonal; rely on the penalization term alone to enforce the cylinder no-slip; fix pressure at a single node to remove the Poisson nullspace; define visc globally and pass constants explicitly; do not use SciPy solver keyword 'tol' (use 'atol'); do not use np.trapz; print progress only every 10% of total steps; save figure(s) but do not call plt.show().

Primary mathematical and numerical challenges:
- Enforcing incompressibility: discrete divergence-free condition requires a pressure-Poisson solve (or equivalent) which has a nullspace (constant pressure mode) unless fixed; boundary conditions for pressure are mixed (Neumann at outlet) which must be handled carefully.
- Irregular immersed boundary (cylinder) inside a rectangular grid: the cylinder boundary is not aligned with the Cartesian grid, so a method that embeds the geometry (here volume-penalization) is required; volume penalization approximates solid enforcement with a stiff Darcy/penalty term scaled by χ/η, where χ is the mask and η a small penalization parameter; this term is stiff and can impose a restrictive time-step if treated explicitly, hence the need to treat it implicitly (fold into momentum diagonal) to avoid undue time-step restriction.
- Temporal discretization: convective term is nonlinear and can be destabilizing if treated explicitly; a balance between stability, cost, and accuracy suggests semi-implicit schemes (explicit advection + implicit diffusion and penalization) or fully implicit (Newton/Krylov) strategies. Explicit treatment of advection imposes a CFL restriction on Δt that scales with U/Δx; viscous stability constraints scale with visc/Δx^2 if viscous terms are explicit.
- Spatial discretization: use a staggered MAC grid to avoid pressure-velocity decoupling (checkerboarding). Implement free-slip and inflow/outflow BCs correctly on the staggered grid, taking care of ghost/edge treatments for derivatives and the Poisson operator.
- Linear algebra: each time-step requires solving elliptic linear systems (momentum blocks for implicit viscous+penalization solves, and pressure Poisson). Efficiency demands reuse of factorizations or good preconditioners (AMG, multigrid, or sparse LU/ILU) because matrices are constant in time (for fixed grid and constant visc and χ mask) allowing precomputation.
- Accuracy vs. penalty parameter and grid resolution: volume penalization introduces an intrinsic modeling error O(η) and a boundary-layer thickness ~O(√η) or grid-scale smearing; controlling this while keeping η small enough to enforce no-slip but not so small as to produce ill-conditioning is a trade-off that interacts with Δt and solver robustness.
- Practical implementation cautions: fix a single pressure degree of freedom; define visc in global scope and pass constants explicitly; do not use removed/forbidden APIs (np.trapz, SciPy tol kwarg). Print progress sparsely (every 10% of steps). Save final figure but do not call plt.show().

### solution_plans
	Current Stage [A/2]
solu_name='Semi-implicit second-order projection on MAC grid with implicit penalization (recommended)' content="Governing idea:\nUse a fractional-step (projection) scheme on a staggered MAC grid. Treat convection explicitly (Adams–Bashforth 2), treat viscous diffusion implicitly (Crank–Nicolson or BDF2) and fold the volume-penalization term implicitly into the diagonal of the momentum operator so the penalization does not shrink Δt. Solve two decoupled velocity linear systems (u and v MAC components) with a constant sparse matrix (factorize once), then enforce incompressibility by solving a pressure Poisson equation and projecting velocities. Fix pressure at a single node to remove nullspace.\n\nAlgorithmic steps (numbered):\n1) Preliminaries and constants:\n   a) Define global visc = 1.0 / Re outside main and pass all constants explicitly to functions (grid sizes Nx, Ny, dt, final_time, η, mask χ, cylinder center/radius, BC flags). Choose η (e.g. 1e-3–1e-4) consistent with grid resolution.\n   b) Build staggered MAC grid: u at vertical faces, v at horizontal faces, p at cell centers. Compute grid spacings dx, dy.\n   c) Build mask chi on the cell-centered grid and interpolate to u/v locations as needed (chi_u, chi_v). Define penalty coefficient α = chi/η evaluated at velocity DOFs.\n2) Discrete operators and constant matrix assembly:\n   a) Assemble the discrete Laplacian for u and v at their staggered locations. Include appropriate BC enforcement for free-slip top/bottom, Dirichlet inlet u=1 on left, and outlet pressure BC handled separately.\n   b) Form constant momentum operator for each velocity component M = (I/dt) + (1/2)*visc*L + diag(α) if using Crank–Nicolson (or (I/dt)+visc*L+B for BDF2 variant). Because visc, dx, dy, eta, and chi are constant the sparse matrix M is constant in time; compute a sparse factorization (sparse LU) or precompute a good preconditioner (ILU or AMG) for fast solves.\n   c) Assemble discrete divergence and gradient operators and the pressure Poisson matrix A_p = div * grad. Impose BCs consistent with velocity BCs and outlet Neumann condition. Fix one pressure node to eliminate nullspace (set that equation to p=0 or add small diagonal penalty only for solver stabilization).\n3) Initial condition and time integration variables:\n   a) Initialize u^0, v^0 (fill inlet area with u=1), and compute an initial pressure (e.g., zero).\n   b) Choose dt based on CFL and viscous limits: dt_CFL = CFL * min(dx/|u_max|, dy/|v_max|), dt_visc = C_visc * min(dx^2,dy^2)/visc; choose dt = safety_factor * min(dt_CFL, dt_visc). Print chosen dt and number of steps Nsteps = ceil(T/dt). Ensure progress prints every 10% of Nsteps.\n4) Time-stepping loop (for n=0..Nsteps-1):\n   a) Compute advective terms explicitly using AB2 (for n=0 use forward-Euler or AB1 to bootstrap): Adv(u) = u*du/dx + v*du/dy, similarly for v. Use conservative interpolation consistent with MAC staggering to compute face-centered advection.\n   b) Form RHS for implicit momentum solve: RHS_u = (I/dt)*u^n - explicit_advection_terms + any explicit pressure-grad term (scheme choice: incremental projection uses old pressure on RHS). Include Dirichlet inflow/ghost enforcement in RHS.\n   c) Solve the linear systems M * u_star = RHS_u and M * v_star = RHS_v using the precomputed LU/solver. Because diag(α) is included, solid velocity relaxes toward zero implicitly.\n   d) Solve the pressure Poisson equation A_p * p^{n+1} = (1/dt) * div(u_star, v_star). Enforce the single-node pressure fix (set p at that cell to 0) and use a robust iterative solver (CG or BiCGStab) with an efficient preconditioner (AMG or ILU). When calling SciPy solvers, do not pass 'tol' but pass 'atol' for the desired absolute tolerance.\n   e) Correct velocities: u^{n+1} = u_star - dt * grad_x(p^{n+1}), v^{n+1} = v_star - dt * grad_y(p^{n+1}). Re-enforce inflow Dirichlet u=1 at inlet ghost/face nodes if necessary to maintain exact inflow.\n   f) Apply any filters or divergence-cleaning if small discrete divergence remains. Compute diagnostics (max divergence, lift/drag if desired using penalization force integral, kinetic energy) occasionally.\n   g) Print progress only every 10% of total steps (concise line with current time, step index, max|u|, max divergence).\n5) Postprocessing and plotting:\n   a) At final time compute velocity magnitude on cell centers (interpolate u/v from MAC faces to centers). Create a contour plot of |u| using 'RdBu_r' colormap; plot the cylinder outline (circle) on top. Save the figure to disk using savefig; do NOT call plt.show().\n\nStability / accuracy / complexity / efficiency limitations:\n- Stability: AB2 explicit treatment of advection imposes a CFL limit on dt that scales with U/max(dx,dy). Implicit viscous plus implicit penalization relieve diffusive and penalization stiffness, but convective CFL still binds. If very large dt are desired, move to fully implicit advection (see Plan B) at higher implementation cost.\n- Accuracy: Projection + AB2 + Crank–Nicolson yields overall second-order accuracy in time (given consistent initialization) and second-order in space with central differences. Volume-penalization introduces an error O(η) and boundary smearing; reduce η and refine grid to improve boundary accuracy at the cost of conditioning.\n- Complexity: Per time step requires two velocity solves (one for each component) plus one scalar Poisson solve. With precomputed factorization of the velocity matrices the per-step cost is dominated by the Poisson solve. Memory and assembly cost are O(N) with sparse storage; CPU cost per time-step is roughly O(N) for multigrid-preconditioned solvers or O(N log N)–O(N^1.5) for direct or poorly preconditioned iterative solvers depending on solver choice.\n- Efficiency: Reusing factorizations or using AMG/multigrid preconditioners for A_p reduces runtime substantially. The penalty parameter η choice trades enforcement accuracy vs. conditioning: too small η forces ill-conditioning; folding α into the diagonal helps but does not remove ill-conditioning entirely. Projection schemes avoid solving fully coupled saddle-point systems at each step, saving complexity.\n\n(End of plan — no code provided.)"

	Current Stage [B/2]
solu_name='Fully-coupled implicit Newton–Krylov (monolithic) with block preconditioning (higher robustness, fewer steps)' content="Governing idea:\nAssemble and solve the coupled saddle-point system for (u,v,p) at each time step using a semi-implicit backward-differentiation (BDF2) time discretization and Newton linearization of nonlinear convection, treating viscous and penalization terms implicitly. Use a Krylov solver (GMRES) for the linearized system with a physics-informed block preconditioner (approximate inverse of the velocity block via multigrid/AMG, and Poisson-like Schur complement approximation for pressure) so larger dt can be used and the number of time steps reduced. Fix pressure at a single node to remove nullspace.\n\nAlgorithmic steps (numbered):\n1) Problem setup and constants:\n   a) Define global visc = 1.0/Re and pass all other constants explicitly. Build MAC grid and masks chi_u, chi_v for penalization. Choose η moderately small (e.g. 1e-4–1e-5) compatible with solver tolerances and grid resolution.\n   b) Discretize spatial operators (convective, viscous Laplacian, gradient, divergence) on the staggered grid. Form the nonlinear residual R(U,P) corresponding to the time-discrete momentum and continuity equations.\n2) Time discretization and nonlinear solve:\n   a) Choose a time-step dt larger than in Plan A if desired, bounded by physical considerations and Newton solver convergence. Use BDF2 for second-order accuracy in time.\n   b) For each time-step perform a Newton (or inexact Newton) iteration: for k=0..K\n       i) Evaluate residual vector R at the current Newton iterate [u_k,v_k,p_k].\n      ii) Form or approximate Jacobian J (a block 3x3 structure [[Auu, Auv, Gx],[Avu, Avv, Gy],[Dx, Dy, 0]]). Include diag(χ/η) in the velocity block A to implicitly enforce penalization.\n     iii) Solve the linear system J * delta = -R using GMRES with a block preconditioner P^{-1} approximating J^{-1}. A standard choice: approximate inverse of velocity block via AMG or multigrid solves, and Schur complement S ≈ -D * A^{-1} * G approximated by a scaled Poisson operator; apply block preconditioner in left or right application inside GMRES.\n      iv) Update the solution [u_{k+1},v_{k+1},p_{k+1}] = [u_k,v_k,p_k] + delta and check Newton convergence. Use an inexact Newton stopping criteria tied to the nonlinear residual reduction.\n3) Preconditioning and linear solver specifics:\n   a) Build or provide AMG/multigrid preconditioner for the velocity-block operator A (including viscous + diag(χ/η)). If a direct solver is used for A as preconditioner, assemble factorizations once and reuse them when coefficients are constant.\n   b) Approximate Schur complement action using a scalable Poisson solve (AMG) for the pressure block inside the preconditioner. Fix one pressure node to remove nullspace; ensure pressure BCs consistent with outlet Neumann by modifying the Poisson operator appropriately.\n   c) When calling GMRES or SciPy iterative solvers, avoid 'tol' keyword and supply 'atol' and a reasonable restart/iteration cap. Optionally use flexible GMRES if preconditioner changes.\n4) Post-Newton step and corrections:\n   a) After Newton convergence, the solution fields satisfy momentum and continuity up to the chosen tolerances. Compute diagnostics (divergence norm, energy, drag via integral of penalization force χ*(u/η)).\n   b) Print progress every 10% of total expected steps. Optionally adjust dt adaptively if Newton struggles or converges rapidly.\n5) Postprocessing:\n   a) Interpolate velocities to cell centers and contour plot velocity magnitude using 'RdBu_r'; overplot cylinder and save the figure. Do not call plt.show().\n\nStability / accuracy / complexity / efficiency limitations:\n- Stability: Fully implicit Newton solves nonlinear stability issues and allows much larger dt than explicit/AB2 projection schemes; however Newton convergence can fail or be slow if dt is too large or initial guess is poor. Good globalization (line search) or time-step control may be necessary.\n- Accuracy: BDF2 + consistent Jacobian gives second-order time accuracy and the staggered MAC discretization retains second-order spatial accuracy. Volume penalization still introduces an O(η) modeling error near the immersed boundary.\n- Complexity: Each time-step requires solving a large coupled linear system via Krylov iterations; preconditioning cost and memory for AMG/block preconditioners are significant. Per-step cost is higher than projection schemes, but fewer time-steps may be needed. Building effective block preconditioners and AMG hierarchies is nontrivial and may require external libraries (PyAMG, PETSc) for performance.\n- Efficiency: For moderate problem sizes a well-designed block preconditioner yields near-optimal performance (iteration counts independent of mesh size); however, implementation and tuning are more complex than the projection approach. If AMG is unavailable, performance may degrade substantially. The implicit approach reduces CFL restrictions but increases solver complexity and memory.\n\n(End of plan — no code provided.)"



### technical_spec
	Current Stage [A/3]
This script implements a semi-implicit projection solver for 2D incompressible Navier–Stokes around a circular cylinder on a staggered MAC grid with volume-penalization (immersed boundary). The main data flow:
- Build MAC grid and cylinder mask (χ) on cell centers; interpolate χ to u/v grid locations to form implicit penalization diagonals.
- Assemble constant sparse matrices once: two momentum matrices (u and v) with implicit viscous Crank–Nicolson-like operator and penalization on the diagonal; one pressure Poisson matrix (with one pressure DOF fixed). Pre-factorize these matrices for reuse.
- Time loop: compute advective terms explicitly (AB2/@Euler for first step), form RHSs, solve two linear systems for provisional velocities (u*, v*), solve Poisson for pressure, project velocities to be divergence free, apply Dirichlet inflow and no-penetration walls, compute diagnostics and print sparse progress every 10% of steps.
- Postprocess: interpolate velocity magnitude to cell centers, plot contour with 'RdBu_r', overlay cylinder outline, and save the figure (no plt.show()).

Key algorithmic choices for balance of stability and simplicity:
- Semi-implicit scheme: explicit AB2 advection, implicit viscous Crank–Nicolson-like operator (matrix assembled once), penalization χ/η treated implicitly (added to momentum diagonal). This avoids penalization-induced timestep restriction.
- Projection method on MAC grid: pressure Poisson built from discrete divergence/gradient; single pressure DOF fixed to remove nullspace.
- Conservative, central-difference spatial discretizations and careful boundary treatment: Dirichlet on inlet (u=1,v=0), free-slip top/bottom via Neumann (∂u/∂y=0) implemented in stencils, v=0 on top/bottom, outflow handled with zero-gradient lightly via Neumann. Cylinder no-slip enforced solely via penalization.

	Current Stage [B/3]
This script implements a semi-implicit projection solver for 2D incompressible Navier–Stokes around a circular cylinder on a staggered (MAC) grid using volume penalization (immersed boundary). Key architecture and data flow:

- Grid and masks: build_grid returns MAC grid sizes and coordinates. Cylinder mask chi is cell-centered and interpolated to u/v faces once. Assertions check interpolation sizes match Nu/Nv.

- Operators assembled once: L_u and L_v (staggered Laplacians) and the pressure Poisson matrix A_p. The momentum matrices A_u and A_v are assembled with implicit viscous and penalization terms and the (1/dt) term.

- Dirichlet enforcement: After assembling A_u/A_v, Dirichlet rows are replaced by identity rows (A[idx,:]=0; A[idx,idx]=1). RHS entries for those rows are set to the physical boundary values (e.g., u=U_in), and the matrices are factorized once (spla.splu). A small unit-check verifies enforced values are returned by the linear solves.

- Time stepping: explicit AB2 for advection, implicit viscous/penalization solved via pre-factorized momentum matrices to obtain provisional velocities; pressure Poisson (pre-factorized) is solved once per step; velocities are projected to be divergence-free. No matrices are rebuilt inside the time loop.

- Boundary consistency: projection gradients at domain faces use zero normal pressure-gradient (consistent with Neumann treatment in Poisson) rather than arbitrary zeros outside the domain.

- Diagnostics & output: progress printed every 10% of steps, final mean/max speed printed, and a contour plot of the velocity magnitude saved (RdBu_r colormap) with the cylinder outline overlayed (no plt.show()).

Notes on robustness/performance: arrays and mask shapes are asserted, Dirichlet enforcement occurs prior to factorization (prevents baked-in inconsistencies), and penalization is included implicitly on momentum diagonals to avoid timestep restrictions.

	Current Stage [C/3]
Semi-implicit projection solver on a staggered MAC grid using volume penalization (immersed boundary).

Architecture summary:
- Grid build: build_grid returns cell/face coordinates and spacings. Cylinder mask built on p-grid and interpolated to u/v faces once.
- Operators: build_laplacian_staggered_u/v and build_poisson_p assemble sparse matrices with the correct sign convention (diagonal = -sum(neighbors)). These are assembled once before time-stepping.
- Momentum matrices: A_u/A_v = (1/dt)I - (visc/2)*L + diag(alpha) (penalization implicit on diagonal). Dirichlet rows are replaced by identity prior to factorization. Matrices are factorized once (spla.splu).
- Advection: explicit, stabilized first-order upwind discretization (for robustness). AB2 time-stepping for advection terms is used to keep second-order temporal accuracy for advection.
- Poisson: pressure Poisson matrix assembled once, nullspace removed by fixing a single pressure cell chosen in the fluid region; LU factorization performed once.
- Time loop: at each step, compute advective terms (upwind), form RHS for momentum, solve for provisional velocities (using pre-factorized LU), compute divergence, solve Poisson for pressure, project velocities. No sparse matrices are rebuilt inside loop.
- Diagnostics: numeric checks for NaN/Inf after major operations; progress prints every 10% of steps; final mean/max speeds and max velocity deviation inside the solid are printed.
- Output: final velocity magnitude contour saved as 'velocity_magnitude.png' with cylinder outline overlayed (RdBu_r colormap).



### python_codes

#### Script block1:

```python
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import ceil

# Global viscosity (use this name exactly as requested)
Re = 100.0
visc = 1.0 / Re

# --------------------------------------------------------------------------------
# Helper functions for grid/mask/operators
# --------------------------------------------------------------------------------

def build_grid(x1, x2, y1, y2, points_per_diameter):
    Lx = x1 + x2
    Ly = y1 + y2
    dx = 1.0 / points_per_diameter  # diameter normalized to 1
    Nx = int(round(Lx / dx))
    Ny = int(round(Ly / dx))
    dx = Lx / Nx
    dy = Ly / Ny
    x_u = np.linspace(0.0, Lx, Nx + 1)  # u at vertical faces
    y_u = (np.arange(Ny) + 0.5) * dy
    x_v = (np.arange(Nx) + 0.5) * dx
    y_v = np.linspace(0.0, Ly, Ny + 1)  # v at horizontal faces
    x_p = (np.arange(Nx) + 0.5) * dx
    y_p = (np.arange(Ny) + 0.5) * dy
    return {
        'Lx': Lx, 'Ly': Ly, 'dx': dx, 'dy': dy,
        'Nx': Nx, 'Ny': Ny,
        'x_u': x_u, 'y_u': y_u,
        'x_v': x_v, 'y_v': y_v,
        'x_p': x_p, 'y_p': y_p
    }


def build_cylinder_mask(Nx, Ny, dx, dy, x1, y1, x2, y2, center_x, center_y, radius):
    # cell-centered coordinates
    x = (np.arange(Nx) + 0.5) * dx
    y = (np.arange(Ny) + 0.5) * dy
    Xc, Yc = np.meshgrid(x, y)
    r2 = (Xc - center_x) ** 2 + (Yc - center_y) ** 2
    chi = (r2 <= radius ** 2).astype(float)
    return chi


def interp_chi_to_u(chi):
    # chi shape (Ny, Nx); u shape (Ny, Nx+1)
    Ny, Nx = chi.shape
    chi_u = np.zeros((Ny, Nx + 1))
    # interior faces: average left and right cell chi
    chi_u[:, 1:Nx] = 0.5 * (chi[:, 0:Nx - 1] + chi[:, 0:Nx])
    # leftmost face i=0: only right cell
    chi_u[:, 0] = chi[:, 0]
    # rightmost face i=Nx: only left cell
    chi_u[:, Nx] = chi[:, Nx - 1]
    return chi_u


def interp_chi_to_v(chi):
    # v shape (Ny+1, Nx)
    Ny, Nx = chi.shape
    chi_v = np.zeros((Ny + 1, Nx))
    # interior faces: average top and bottom cells
    chi_v[1:Ny, :] = 0.5 * (chi[0:Ny - 1, :] + chi[1:Ny, :])
    # bottom j=0: only top cell j=0
    chi_v[0, :] = chi[0, :]
    # top j=Ny: only bottom cell j=Ny-1
    chi_v[Ny, :] = chi[Ny - 1, :]
    return chi_v


def node_index_u(i, j, Nx, Ny):
    # u index for (j,i) where i in 0..Nx, j in 0..Ny-1
    return j * (Nx + 1) + i


def node_index_v(i, j, Nx, Ny):
    # v index for (j,i) where i in 0..Nx-1, j in 0..Ny
    return j * Nx + i


def node_index_p(i, j, Nx, Ny):
    # p index for (j,i) where i in 0..Nx-1, j in 0..Ny-1
    return j * Nx + i


def build_laplacian_staggered_u(Nx, Ny, dx, dy, bcs):
    # Build Laplacian operator for u-grid (Nu x Nu)
    Nu = Ny * (Nx + 1)
    rows = []
    cols = []
    data = []

    for j in range(Ny):
        for i in range(Nx + 1):
            idx = node_index_u(i, j, Nx, Ny)
            # Dirichlet at inlet (i==0) if requested
            if bcs['u_dirichlet_mask'][j, i]:
                rows.append(idx); cols.append(idx); data.append(1.0)
                continue
            # otherwise build finite-difference with Neumann at missing neighbor
            center_coef = 0.0
            # x-direction
            # left neighbor
            if i - 1 >= 0:
                nidx = node_index_u(i - 1, j, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dx ** 2)
                center_coef -= 1.0 / dx ** 2
            else:
                # missing neighbor: if Neumann, emulate ghost = center -> contribute center -1/dx^2
                if bcs['u_x_neumann_left']:
                    center_coef -= 1.0 / dx ** 2
                else:
                    # if Dirichlet and not captured, leave center zero
                    pass
            # right neighbor
            if i + 1 <= Nx:
                nidx = node_index_u(i + 1, j, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dx ** 2)
                center_coef -= 1.0 / dx ** 2
            else:
                if bcs['u_x_neumann_right']:
                    center_coef -= 1.0 / dx ** 2
            # y-direction (use u nodes j-1 and j+1 exist between u vertical faces)
            if j - 1 >= 0:
                nidx = node_index_u(i, j - 1, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dy ** 2)
                center_coef -= 1.0 / dy ** 2
            else:
                # bottom Neumann free-slip for u: ∂u/∂y = 0
                if bcs['u_y_neumann_bottom']:
                    center_coef -= 1.0 / dy ** 2
            if j + 1 <= Ny - 1:
                nidx = node_index_u(i, j + 1, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dy ** 2)
                center_coef -= 1.0 / dy ** 2
            else:
                if bcs['u_y_neumann_top']:
                    center_coef -= 1.0 / dy ** 2
            # center
            rows.append(idx); cols.append(idx); data.append(-center_coef)
    L = sp.csr_matrix((data, (rows, cols)), shape=(Nu, Nu))
    return L


def build_laplacian_staggered_v(Nx, Ny, dx, dy, bcs):
    Nv = (Ny + 1) * Nx
    rows = []
    cols = []
    data = []
    for j in range(Ny + 1):
        for i in range(Nx):
            idx = node_index_v(i, j, Nx, Ny)
            if bcs['v_dirichlet_mask'][j, i]:
                rows.append(idx); cols.append(idx); data.append(1.0)
                continue
            center_coef = 0.0
            # x-direction neighbors
            if i - 1 >= 0:
                nidx = node_index_v(i - 1, j, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dx ** 2)
                center_coef -= 1.0 / dx ** 2
            else:
                if bcs['v_x_neumann_left']:
                    center_coef -= 1.0 / dx ** 2
            if i + 1 <= Nx - 1:
                nidx = node_index_v(i + 1, j, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dx ** 2)
                center_coef -= 1.0 / dx ** 2
            else:
                if bcs['v_x_neumann_right']:
                    center_coef -= 1.0 / dx ** 2
            # y-direction neighbors
            if j - 1 >= 0:
                nidx = node_index_v(i, j - 1, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dy ** 2)
                center_coef -= 1.0 / dy ** 2
            else:
                if bcs['v_y_neumann_bottom']:
                    center_coef -= 1.0 / dy ** 2
            if j + 1 <= Ny:
                nidx = node_index_v(i, j + 1, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dy ** 2)
                center_coef -= 1.0 / dy ** 2
            else:
                if bcs['v_y_neumann_top']:
                    center_coef -= 1.0 / dy ** 2
            rows.append(idx); cols.append(idx); data.append(-center_coef)
    L = sp.csr_matrix((data, (rows, cols)), shape=(Nv, Nv))
    return L


def build_poisson_p(Nx, Ny, dx, dy):
    Np = Nx * Ny
    rows = []
    cols = []
    data = []
    for j in range(Ny):
        for i in range(Nx):
            idx = node_index_p(i, j, Nx, Ny)
            center_coef = 0.0
            # left
            if i - 1 >= 0:
                nidx = node_index_p(i - 1, j, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dx ** 2)
                center_coef -= 1.0 / dx ** 2
            else:
                # Neumann (zero normal derivative) -> one-sided contributes -1/dx^2
                center_coef -= 1.0 / dx ** 2
            # right
            if i + 1 <= Nx - 1:
                nidx = node_index_p(i + 1, j, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dx ** 2)
                center_coef -= 1.0 / dx ** 2
            else:
                center_coef -= 1.0 / dx ** 2
            # bottom
            if j - 1 >= 0:
                nidx = node_index_p(i, j - 1, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dy ** 2)
                center_coef -= 1.0 / dy ** 2
            else:
                center_coef -= 1.0 / dy ** 2
            # top
            if j + 1 <= Ny - 1:
                nidx = node_index_p(i, j + 1, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dy ** 2)
                center_coef -= 1.0 / dy ** 2
            else:
                center_coef -= 1.0 / dy ** 2
            rows.append(idx); cols.append(idx); data.append(-center_coef)
    A = sp.csr_matrix((data, (rows, cols)), shape=(Np, Np))
    return A

# --------------------------------------------------------------------------------
# Main solver
# --------------------------------------------------------------------------------

def main():
    # Problem geometry and parameters (explicitly passed)
    d = 1.0
    x1 = 3.0
    x2 = 8.0
    y1 = 3.0
    y2 = 3.0
    Lx = x1 + x2
    Ly = y1 + y2
    cx = x1  # cylinder center x
    cy = Ly / 2.0
    radius = d / 2.0

    # grid resolution (points per diameter) - moderate
    ppd = 12  # points per diameter -> reasonable runtime and visible wake
    grid = build_grid(x1, x2, y1, y2, ppd)
    Nx = grid['Nx']; Ny = grid['Ny']; dx = grid['dx']; dy = grid['dy']

    # penalization parameter eta (keep not too small to avoid ill-conditioning)
    eta = 1e-3

    # Build masks and interpolations
    chi = build_cylinder_mask(Nx, Ny, dx, dy, x1, y1, x2, y2, cx, cy, radius)
    chi_u = interp_chi_to_u(chi)
    chi_v = interp_chi_to_v(chi)
    alpha_u = (chi_u / eta).ravel()
    alpha_v = (chi_v / eta).ravel()

    # boundary condition masks for assembling Laplacians
    # For u-grid (Ny, Nx+1): Dirichlet at inlet (i==0) set to u=1
    u_dirichlet_mask = np.zeros((Ny, Nx + 1), dtype=bool)
    u_dirichlet_mask[:, 0] = True
    # for u free-slip top/bottom -> Neumann in y
    bcs_u = {
        'u_dirichlet_mask': u_dirichlet_mask,
        'u_x_neumann_left': False,
        'u_x_neumann_right': True,  # outflow ~ zero gradient
        'u_y_neumann_bottom': True,
        'u_y_neumann_top': True
    }
    # For v-grid (Ny+1, Nx): Dirichlet v=0 at top/bottom and v=0 at inlet (i==0)
    v_dirichlet_mask = np.zeros((Ny + 1, Nx), dtype=bool)
    v_dirichlet_mask[0, :] = True  # bottom wall v=0
    v_dirichlet_mask[-1, :] = True  # top wall v=0
    v_dirichlet_mask[:, 0] = True  # enforce v=0 at inlet
    bcs_v = {
        'v_dirichlet_mask': v_dirichlet_mask,
        'v_x_neumann_left': True,
        'v_x_neumann_right': True,
        'v_y_neumann_bottom': False,
        'v_y_neumann_top': False
    }

    # Build Laplacians once
    print('Assembling matrices...')
    L_u = build_laplacian_staggered_u(Nx, Ny, dx, dy, bcs_u)
    L_v = build_laplacian_staggered_v(Nx, Ny, dx, dy, bcs_v)

    Nu = Ny * (Nx + 1)
    Nv = (Ny + 1) * Nx

    # Momentum operator: A = (1/dt)*I - (visc/2)*L + diag(alpha)
    # Since dt unknown at assembly time, we'll assemble base matrix components and add (1/dt)I at solve-time by adding diagonal.
    # But we must not rebuild matrices inside time loop; instead we will assemble full A using chosen dt below and then factorize.

    # Poisson matrix
    A_p = build_poisson_p(Nx, Ny, dx, dy)

    # Initial condition: u=1 everywhere (inflow) and v=0
    u = np.ones((Ny, Nx + 1))
    v = np.zeros((Ny + 1, Nx))
    # Inside cylinder we let penalization enforce u,v->0; initial guess keeps inlet

    # Pressure initial guess
    p = np.zeros((Ny, Nx))

    # Time stepping parameters
    U_in = 1.0
    umax = U_in
    CFL = 0.5
    dt_CFL = CFL * min(dx / max(1e-8, umax), dy / max(1e-8, umax))
    dt_visc = 0.5 * min(dx * dx, dy * dy) / max(1e-12, visc)
    dt = min(dt_CFL, dt_visc)
    # Be conservative
    dt = min(dt, 0.02)
    T = 20.0
    Nsteps = int(ceil(T / dt))
    dt = T / Nsteps
    print(f'Grid: Nx={Nx}, Ny={Ny}, dx={dx:.4f}, dy={dy:.4f}, dt={dt:.5f}, steps={Nsteps}')

    # Build final momentum matrices with dt and penalization included
    # A_u = (1/dt)*I - (visc/2)*L_u + diag(alpha_u)
    Iu = sp.eye(Nu, format='csr')
    Iv = sp.eye(Nv, format='csr')
    A_u = (1.0 / dt) * Iu - (visc / 2.0) * L_u + sp.diags(alpha_u)
    A_v = (1.0 / dt) * Iv - (visc / 2.0) * L_v + sp.diags(alpha_v)

    # Enforce Dirichlet rows exactly: for any dirichlet mask node, replace row by identity (already handled in Laplacian) and ensure diagonal large
    # Factorize momentum matrices
    print('Factorizing momentum matrices (one-time)...')
    A_u_csc = A_u.tocsc()
    A_v_csc = A_v.tocsc()
    lu_u = spla.splu(A_u_csc)
    lu_v = spla.splu(A_v_csc)

    # Factorize Poisson after pressure fix; but first modify A_p to fix one pressure DOF to remove nullspace
    Np = Nx * Ny
    # Choose pressure fix node index (e.g., cell nearest to outlet center)
    fix_i = Nx - 1
    fix_j = Ny // 2
    fix_idx = node_index_p(fix_i, fix_j, Nx, Ny)
    A_p = A_p.tolil()
    A_p[fix_idx, :] = 0.0
    A_p[fix_idx, fix_idx] = 1.0
    A_p = A_p.tocsc()
    print('Factorizing Poisson matrix (one-time)...')
    lu_p = spla.splu(A_p)

    # Prepare AB2 storage
    adv_u_old = np.zeros_like(u).ravel()
    adv_v_old = np.zeros_like(v).ravel()

    # diagnostics
    def divergence_field(u, v, dx, dy, Nx, Ny):
        div = np.zeros((Ny, Nx))
        # div at p centers: (u[i+1]-u[i])/dx + (v[j+1]-v[j])/dy
        for j in range(Ny):
            for i in range(Nx):
                div[j, i] = (u[j, i + 1] - u[j, i]) / dx + (v[j + 1, i] - v[j, i]) / dy
        return div

    # Time loop
    print('Starting time-stepping...')
    report_every = max(1, Nsteps // 10)
    for step in range(Nsteps):
        t = (step + 1) * dt

        # Compute advective terms for u and v
        # adv_u = u * du/dx + v_at_u * du/dy
        adv_u = np.zeros_like(u)
        adv_v = np.zeros_like(v)

        # compute du/dx and du/dy on u-grid
        for j in range(Ny):
            for i in range(Nx + 1):
                ui = u[j, i]
                # du/dx
                if i == 0:
                    du_dx = (u[j, i + 1] - u[j, i]) / dx
                elif i == Nx:
                    du_dx = (u[j, i] - u[j, i - 1]) / dx
                else:
                    du_dx = (u[j, i + 1] - u[j, i - 1]) / (2.0 * dx)
                # du/dy
                if j == 0:
                    du_dy = (u[j + 1, i] - u[j, i]) / dy
                elif j == Ny - 1:
                    du_dy = (u[j, i] - u[j - 1, i]) / dy
                else:
                    du_dy = (u[j + 1, i] - u[j - 1, i]) / (2.0 * dy)
                # v at u location: average surrounding v nodes (bilinear-like)
                # v indices go j_v=0..Ny, i_v=0..Nx-1. For u at (j,i) average v at i-1 and i columns
                vvals = []
                for jj in [j, j + 1]:
                    for ii in [i - 1, i]:
                        if 0 <= jj <= Ny and 0 <= ii <= Nx - 1:
                            vvals.append(v[jj, ii])
                v_at_u = np.mean(vvals) if vvals else 0.0
                adv_u[j, i] = ui * du_dx + v_at_u * du_dy

        # adv_v
        for j in range(Ny + 1):
            for i in range(Nx):
                vi = v[j, i]
                # dv/dx
                if i == 0:
                    dv_dx = (v[j, i + 1] - v[j, i]) / dx
                elif i == Nx - 1:
                    dv_dx = (v[j, i] - v[j, i - 1]) / dx
                else:
                    dv_dx = (v[j, i + 1] - v[j, i - 1]) / (2.0 * dx)
                # dv/dy
                if j == 0:
                    dv_dy = (v[j + 1, i] - v[j, i]) / dy
                elif j == Ny:
                    dv_dy = (v[j, i] - v[j - 1, i]) / dy
                else:
                    dv_dy = (v[j + 1, i] - v[j - 1, i]) / (2.0 * dy)
                # u at v location: average surrounding u nodes
                uvals = []
                for jj in [j - 1, j]:
                    for ii in [i, i + 1]:
                        if 0 <= jj <= Ny - 1 and 0 <= ii <= Nx:
                            uvals.append(u[jj, ii])
                u_at_v = np.mean(uvals) if uvals else 0.0
                adv_v[j, i] = u_at_v * dv_dx + vi * dv_dy

        # Flatten adv arrays
        adv_u_flat = adv_u.ravel()
        adv_v_flat = adv_v.ravel()

        # AB2 for advection: adv_term = 1.5*adv^n - 0.5*adv^{n-1}
        if step == 0:
            adv_u_eff = adv_u_flat
            adv_v_eff = adv_v_flat
        else:
            adv_u_eff = 1.5 * adv_u_flat - 0.5 * adv_u_old
            adv_v_eff = 1.5 * adv_v_flat - 0.5 * adv_v_old

        adv_u_old = adv_u_flat.copy()
        adv_v_old = adv_v_flat.copy()

        # Build RHS for u: (1/dt)*u^n + (visc/2)*L_u*u^n - adv_eff (no pressure grad for intermediate step)
        u_flat = u.ravel()
        v_flat = v.ravel()
        RHS_u = (1.0 / dt) * u_flat + (visc / 2.0) * (L_u.dot(u_flat)) - adv_u_eff
        RHS_v = (1.0 / dt) * v_flat + (visc / 2.0) * (L_v.dot(v_flat)) - adv_v_eff

        # Enforce Dirichlet values in RHS for u and v (inlet and walls)
        # For u inlet i==0: set RHS to 1.0
        # obtain list of dirichlet indices for u and v
        u_dir_indices = np.where(u_dirichlet_mask.ravel())[0]
        for idx in u_dir_indices:
            RHS_u[idx] = U_in / dt  # since matrix row is identity, RHS must be value*(1/dt)
        v_dir_indices = np.where(v_dirichlet_mask.ravel())[0]
        for idx in v_dir_indices:
            RHS_v[idx] = 0.0

        # Solve for provisional velocities u_star, v_star
        u_star_flat = lu_u.solve(RHS_u)
        v_star_flat = lu_v.solve(RHS_v)

        u_star = u_star_flat.reshape(u.shape)
        v_star = v_star_flat.reshape(v.shape)

        # Solve Poisson: A_p * p_flat = (1/dt) * div(u_star, v_star) with pressure fix already built into A_p
        div_star = divergence_field(u_star, v_star, dx, dy, Nx, Ny)
        rhs_p = (1.0 / dt) * div_star.ravel()
        # enforce pressure fix RHS at fix_idx
        rhs_p[fix_idx] = 0.0
        p_flat = lu_p.solve(rhs_p)
        p = p_flat.reshape((Ny, Nx))

        # Project velocities: u^{n+1} = u_star - dt * grad_x(p)
        u_new = u_star.copy()
        v_new = v_star.copy()

        # grad_x at u faces: (p[i] - p[i-1]) / dx
        for j in range(Ny):
            for i in range(Nx + 1):
                if u_dirichlet_mask[j, i]:
                    # enforce inlet exactly
                    u_new[j, i] = U_in
                    continue
                if i == 0:
                    # left face (outside p domain): use one-sided
                    dpdx = (p[j, 0] - 0.0) / dx  # assume p outside ~ 0, but pressure at inlet not fixed; small impact
                elif i == Nx:
                    dpdx = (0.0 - p[j, Nx - 1]) / dx
                else:
                    dpdx = (p[j, i] - p[j, i - 1]) / dx
                u_new[j, i] = u_star[j, i] - dt * dpdx
        # grad_y at v faces: (p[j] - p[j-1]) / dy
        for j in range(Ny + 1):
            for i in range(Nx):
                if v_dirichlet_mask[j, i]:
                    v_new[j, i] = 0.0
                    continue
                if j == 0:
                    dpdy = (p[0, i] - 0.0) / dy
                elif j == Ny:
                    dpdy = (0.0 - p[Ny - 1, i]) / dy
                else:
                    dpdy = (p[j, i] - p[j - 1, i]) / dy
                v_new[j, i] = v_star[j, i] - dt * dpdy

        # Overwrite velocities
        u = u_new
        v = v_new

        # Diagnostics and progress prints only every 10% of steps
        if (step + 1) % report_every == 0 or (step == Nsteps - 1):
            div_norm = np.max(np.abs(divergence_field(u, v, dx, dy, Nx, Ny)))
            umax = np.max(np.abs(u))
            vmax = np.max(np.abs(v))
            print(f'Step {step+1}/{Nsteps}, time {t:.3f}, umax {umax:.4f}, vmax {vmax:.4f}, max|div| {div_norm:.3e}')

    # Postprocessing: compute velocity magnitude at cell centers by averaging surrounding u/v
    Uc = np.zeros((Ny, Nx))
    Vc = np.zeros((Ny, Nx))
    for j in range(Ny):
        for i in range(Nx):
            # u at left and right faces
            Uc[j, i] = 0.5 * (u[j, i] + u[j, i + 1])
            # v at bottom and top faces
            Vc[j, i] = 0.5 * (v[j, i] + v[j + 1, i])
    speed = np.sqrt(Uc ** 2 + Vc ** 2)

    # Plot contour of velocity magnitude with cylinder
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = (np.arange(Nx) + 0.5) * dx
    y = (np.arange(Ny) + 0.5) * dy
    X, Y = np.meshgrid(x, y)
    cs = ax.contourf(X, Y, speed, levels=50, cmap='RdBu_r')
    cbar = fig.colorbar(cs)
    cbar.set_label('Velocity magnitude')
    # plot cylinder outline
    theta = np.linspace(0, 2 * np.pi, 200)
    xc = cx + radius * np.cos(theta)
    yc = cy + radius * np.sin(theta)
    ax.plot(xc, yc, 'k', linewidth=1.2)
    ax.set_aspect('equal')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_title(f'Velocity magnitude at t={T:.2f} (Re={Re:.0f})')
    plt.tight_layout()
    fig.savefig('velocity_magnitude.png', dpi=150)
    print('Saved figure to velocity_magnitude.png')


if __name__ == '__main__':
    main()
```


#### Script block2:

```python
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import ceil

# Global viscosity (use this name exactly as requested)
Re = 100.0
visc = 1.0 / Re

# --------------------------------------------------------------------------------
# Helper functions for grid/mask/operators
# --------------------------------------------------------------------------------

def build_grid(x1, x2, y1, y2, points_per_diameter):
    Lx = x1 + x2
    Ly = y1 + y2
    dx = 1.0 / points_per_diameter  # diameter normalized to 1
    Nx = int(round(Lx / dx))
    Ny = int(round(Ly / dx))
    dx = Lx / Nx
    dy = Ly / Ny
    x_u = np.linspace(0.0, Lx, Nx + 1)  # u at vertical faces
    y_u = (np.arange(Ny) + 0.5) * dy
    x_v = (np.arange(Nx) + 0.5) * dx
    y_v = np.linspace(0.0, Ly, Ny + 1)  # v at horizontal faces
    x_p = (np.arange(Nx) + 0.5) * dx
    y_p = (np.arange(Ny) + 0.5) * dy
    return {
        'Lx': Lx, 'Ly': Ly, 'dx': dx, 'dy': dy,
        'Nx': Nx, 'Ny': Ny,
        'x_u': x_u, 'y_u': y_u,
        'x_v': x_v, 'y_v': y_v,
        'x_p': x_p, 'y_p': y_p
    }


def build_cylinder_mask(Nx, Ny, dx, dy, x1, y1, x2, y2, center_x, center_y, radius):
    # cell-centered coordinates
    x = (np.arange(Nx) + 0.5) * dx
    y = (np.arange(Ny) + 0.5) * dy
    Xc, Yc = np.meshgrid(x, y)
    r2 = (Xc - center_x) ** 2 + (Yc - center_y) ** 2
    chi = (r2 <= radius ** 2).astype(float)
    return chi


def interp_chi_to_u(chi):
    # chi shape (Ny, Nx); u shape (Ny, Nx+1)
    Ny, Nx = chi.shape
    chi_u = np.zeros((Ny, Nx + 1))
    # interior faces: average left and right cell chi
    # fix: use correct right-cell slice (1:Nx) to match 0:Nx-1
    chi_u[:, 1:Nx] = 0.5 * (chi[:, 0:Nx - 1] + chi[:, 1:Nx])
    # leftmost face i=0: only right cell
    chi_u[:, 0] = chi[:, 0]
    # rightmost face i=Nx: only left cell
    chi_u[:, Nx] = chi[:, Nx - 1]
    return chi_u


def interp_chi_to_v(chi):
    # v shape (Ny+1, Nx)
    Ny, Nx = chi.shape
    chi_v = np.zeros((Ny + 1, Nx))
    # interior faces: average top and bottom cells
    chi_v[1:Ny, :] = 0.5 * (chi[0:Ny - 1, :] + chi[1:Ny, :])
    # bottom j=0: only top cell j=0
    chi_v[0, :] = chi[0, :]
    # top j=Ny: only bottom cell j=Ny-1
    chi_v[Ny, :] = chi[Ny - 1, :]
    return chi_v


def node_index_u(i, j, Nx, Ny):
    # u index for (j,i) where i in 0..Nx, j in 0..Ny-1
    return j * (Nx + 1) + i


def node_index_v(i, j, Nx, Ny):
    # v index for (j,i) where i in 0..Nx-1, j in 0..Ny
    return j * Nx + i


def node_index_p(i, j, Nx, Ny):
    # p index for (j,i) where i in 0..Nx-1, j in 0..Ny-1
    return j * Nx + i


def build_laplacian_staggered_u(Nx, Ny, dx, dy, bcs):
    # Build Laplacian operator for u-grid (Nu x Nu)
    Nu = Ny * (Nx + 1)
    rows = []
    cols = []
    data = []

    for j in range(Ny):
        for i in range(Nx + 1):
            idx = node_index_u(i, j, Nx, Ny)
            # Dirichlet at inlet (i==0) if requested
            if bcs['u_dirichlet_mask'][j, i]:
                rows.append(idx); cols.append(idx); data.append(1.0)
                continue
            center_coef = 0.0
            # x-direction
            # left neighbor
            if i - 1 >= 0:
                nidx = node_index_u(i - 1, j, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dx ** 2)
                center_coef -= 1.0 / dx ** 2
            else:
                if bcs['u_x_neumann_left']:
                    center_coef -= 1.0 / dx ** 2
            # right neighbor
            if i + 1 <= Nx:
                nidx = node_index_u(i + 1, j, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dx ** 2)
                center_coef -= 1.0 / dx ** 2
            else:
                if bcs['u_x_neumann_right']:
                    center_coef -= 1.0 / dx ** 2
            # y-direction
            if j - 1 >= 0:
                nidx = node_index_u(i, j - 1, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dy ** 2)
                center_coef -= 1.0 / dy ** 2
            else:
                if bcs['u_y_neumann_bottom']:
                    center_coef -= 1.0 / dy ** 2
            if j + 1 <= Ny - 1:
                nidx = node_index_u(i, j + 1, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dy ** 2)
                center_coef -= 1.0 / dy ** 2
            else:
                if bcs['u_y_neumann_top']:
                    center_coef -= 1.0 / dy ** 2
            rows.append(idx); cols.append(idx); data.append(-center_coef)
    L = sp.csr_matrix((data, (rows, cols)), shape=(Nu, Nu))
    return L


def build_laplacian_staggered_v(Nx, Ny, dx, dy, bcs):
    Nv = (Ny + 1) * Nx
    rows = []
    cols = []
    data = []
    for j in range(Ny + 1):
        for i in range(Nx):
            idx = node_index_v(i, j, Nx, Ny)
            if bcs['v_dirichlet_mask'][j, i]:
                rows.append(idx); cols.append(idx); data.append(1.0)
                continue
            center_coef = 0.0
            # x-direction neighbors
            if i - 1 >= 0:
                nidx = node_index_v(i - 1, j, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dx ** 2)
                center_coef -= 1.0 / dx ** 2
            else:
                if bcs['v_x_neumann_left']:
                    center_coef -= 1.0 / dx ** 2
            if i + 1 <= Nx - 1:
                nidx = node_index_v(i + 1, j, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dx ** 2)
                center_coef -= 1.0 / dx ** 2
            else:
                if bcs['v_x_neumann_right']:
                    center_coef -= 1.0 / dx ** 2
            # y-direction neighbors
            if j - 1 >= 0:
                nidx = node_index_v(i, j - 1, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dy ** 2)
                center_coef -= 1.0 / dy ** 2
            else:
                if bcs['v_y_neumann_bottom']:
                    center_coef -= 1.0 / dy ** 2
            if j + 1 <= Ny:
                nidx = node_index_v(i, j + 1, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dy ** 2)
                center_coef -= 1.0 / dy ** 2
            else:
                if bcs['v_y_neumann_top']:
                    center_coef -= 1.0 / dy ** 2
            rows.append(idx); cols.append(idx); data.append(-center_coef)
    L = sp.csr_matrix((data, (rows, cols)), shape=(Nv, Nv))
    return L


def build_poisson_p(Nx, Ny, dx, dy):
    Np = Nx * Ny
    rows = []
    cols = []
    data = []
    for j in range(Ny):
        for i in range(Nx):
            idx = node_index_p(i, j, Nx, Ny)
            center_coef = 0.0
            # left
            if i - 1 >= 0:
                nidx = node_index_p(i - 1, j, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dx ** 2)
                center_coef -= 1.0 / dx ** 2
            else:
                # Neumann (zero normal derivative) -> one-sided contributes -1/dx^2
                center_coef -= 1.0 / dx ** 2
            # right
            if i + 1 <= Nx - 1:
                nidx = node_index_p(i + 1, j, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dx ** 2)
                center_coef -= 1.0 / dx ** 2
            else:
                center_coef -= 1.0 / dx ** 2
            # bottom
            if j - 1 >= 0:
                nidx = node_index_p(i, j - 1, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dy ** 2)
                center_coef -= 1.0 / dy ** 2
            else:
                center_coef -= 1.0 / dy ** 2
            # top
            if j + 1 <= Ny - 1:
                nidx = node_index_p(i, j + 1, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dy ** 2)
                center_coef -= 1.0 / dy ** 2
            else:
                center_coef -= 1.0 / dy ** 2
            rows.append(idx); cols.append(idx); data.append(-center_coef)
    A = sp.csr_matrix((data, (rows, cols)), shape=(Np, Np))
    return A


# --------------------------------------------------------------------------------
# Main solver
# --------------------------------------------------------------------------------

def main():
    # Problem geometry and parameters (explicitly passed)
    d = 1.0
    x1 = 3.0
    x2 = 8.0
    y1 = 3.0
    y2 = 3.0
    Lx = x1 + x2
    Ly = y1 + y2
    cx = x1  # cylinder center x
    cy = Ly / 2.0
    radius = d / 2.0

    # grid resolution (points per diameter) - moderate
    ppd = 12  # points per diameter -> reasonable runtime and visible wake
    grid = build_grid(x1, x2, y1, y2, ppd)
    Nx = grid['Nx']; Ny = grid['Ny']; dx = grid['dx']; dy = grid['dy']

    # penalization parameter eta (keep not too small to avoid ill-conditioning)
    eta = 1e-3

    # Build masks and interpolations
    chi = build_cylinder_mask(Nx, Ny, dx, dy, x1, y1, x2, y2, cx, cy, radius)
    chi_u = interp_chi_to_u(chi)
    chi_v = interp_chi_to_v(chi)

    Nu = Ny * (Nx + 1)
    Nv = (Ny + 1) * Nx
    assert chi_u.size == Nu, f"chi_u.size {chi_u.size} != Nu {Nu}"
    assert chi_v.size == Nv, f"chi_v.size {chi_v.size} != Nv {Nv}"

    alpha_u = (chi_u / eta).ravel()
    alpha_v = (chi_v / eta).ravel()

    # boundary condition masks for assembling Laplacians
    # For u-grid (Ny, Nx+1): Dirichlet at inlet (i==0) set to u=1
    u_dirichlet_mask = np.zeros((Ny, Nx + 1), dtype=bool)
    u_dirichlet_mask[:, 0] = True
    # for u free-slip top/bottom -> Neumann in y
    bcs_u = {
        'u_dirichlet_mask': u_dirichlet_mask,
        'u_x_neumann_left': False,
        'u_x_neumann_right': True,  # outflow ~ zero gradient
        'u_y_neumann_bottom': True,
        'u_y_neumann_top': True
    }
    # For v-grid (Ny+1, Nx): Dirichlet v=0 at top/bottom and v=0 at inlet (i==0)
    v_dirichlet_mask = np.zeros((Ny + 1, Nx), dtype=bool)
    v_dirichlet_mask[0, :] = True  # bottom wall v=0
    v_dirichlet_mask[-1, :] = True  # top wall v=0
    v_dirichlet_mask[:, 0] = True  # enforce v=0 at inlet
    # Because we set Dirichlet at inlet for v, the corresponding Neumann flag should be False
    bcs_v = {
        'v_dirichlet_mask': v_dirichlet_mask,
        'v_x_neumann_left': False,
        'v_x_neumann_right': True,
        'v_y_neumann_bottom': False,
        'v_y_neumann_top': False
    }

    # Build Laplacians once
    print('Assembling matrices...')
    L_u = build_laplacian_staggered_u(Nx, Ny, dx, dy, bcs_u)
    L_v = build_laplacian_staggered_v(Nx, Ny, dx, dy, bcs_v)

    # Poisson matrix
    A_p = build_poisson_p(Nx, Ny, dx, dy)

    # Initial condition: u=1 everywhere (inflow) and v=0
    u = np.ones((Ny, Nx + 1))
    v = np.zeros((Ny + 1, Nx))

    # Pressure initial guess
    p = np.zeros((Ny, Nx))

    # Time stepping parameters
    U_in = 1.0
    umax = U_in
    CFL = 0.5
    dt_CFL = CFL * min(dx / max(1e-8, umax), dy / max(1e-8, umax))
    dt_visc = 0.5 * min(dx * dx, dy * dy) / max(1e-12, visc)
    dt = min(dt_CFL, dt_visc)
    # Be conservative
    dt = min(dt, 0.02)
    T = 20.0
    Nsteps = int(ceil(T / dt))
    dt = T / Nsteps
    print(f'Grid: Nx={Nx}, Ny={Ny}, dx={dx:.4f}, dy={dy:.4f}, dt={dt:.5f}, steps={Nsteps}')

    # Build final momentum matrices with dt and penalization included
    Iu = sp.eye(Nu, format='csr')
    Iv = sp.eye(Nv, format='csr')
    A_u = (1.0 / dt) * Iu - (visc / 2.0) * L_u + sp.diags(alpha_u)
    A_v = (1.0 / dt) * Iv - (visc / 2.0) * L_v + sp.diags(alpha_v)

    # Enforce Dirichlet rows on full momentum matrices before factorization
    print('Applying Dirichlet enforcement to momentum matrices...')
    A_u = A_u.tolil()
    A_v = A_v.tolil()

    u_dir_indices = np.where(u_dirichlet_mask.ravel())[0]
    v_dir_indices = np.where(v_dirichlet_mask.ravel())[0]
    assert u_dir_indices.size > 0, 'Expected some u Dirichlet indices'
    assert v_dir_indices.size > 0, 'Expected some v Dirichlet indices'

    # Replace rows with identity to enforce u = U_in and v = 0 exactly
    for idx in u_dir_indices:
        A_u.rows[idx] = [idx]
        A_u.data[idx] = [1.0]
    for idx in v_dir_indices:
        A_v.rows[idx] = [idx]
        A_v.data[idx] = [1.0]

    A_u = A_u.tocsc()
    A_v = A_v.tocsc()

    # Factorize momentum matrices
    print('Factorizing momentum matrices (one-time)...')
    lu_u = spla.splu(A_u)
    lu_v = spla.splu(A_v)

    # Quick unit-check: solving with RHS equal to Dirichlet targets yields those values at masked indices
    # Build test RHS
    test_rhs_u = np.zeros(Nu)
    test_rhs_v = np.zeros(Nv)
    test_rhs_u[u_dir_indices] = U_in
    test_rhs_v[v_dir_indices] = 0.0
    sol_u_test = lu_u.solve(test_rhs_u)
    sol_v_test = lu_v.solve(test_rhs_v)
    # verify enforcement
    assert np.allclose(sol_u_test[u_dir_indices], U_in, atol=1e-12), 'Dirichlet u enforcement failed'
    assert np.allclose(sol_v_test[v_dir_indices], 0.0, atol=1e-12), 'Dirichlet v enforcement failed'

    # Factorize Poisson after pressure fix; modify A_p to fix one pressure DOF to remove nullspace
    Np = Nx * Ny
    # Choose pressure fix node index (e.g., cell nearest to outlet center)
    fix_i = Nx - 1
    fix_j = Ny // 2
    fix_idx = node_index_p(fix_i, fix_j, Nx, Ny)
    A_p = A_p.tolil()
    A_p[fix_idx, :] = 0.0
    A_p[fix_idx, fix_idx] = 1.0
    A_p = A_p.tocsc()
    print('Factorizing Poisson matrix (one-time)...')
    lu_p = spla.splu(A_p)

    # Prepare AB2 storage
    adv_u_old = np.zeros_like(u).ravel()
    adv_v_old = np.zeros_like(v).ravel()

    # diagnostics
    def divergence_field(u_field, v_field, dx_loc, dy_loc, Nx_loc, Ny_loc):
        div = np.zeros((Ny_loc, Nx_loc))
        for j in range(Ny_loc):
            for i in range(Nx_loc):
                div[j, i] = (u_field[j, i + 1] - u_field[j, i]) / dx_loc + (v_field[j + 1, i] - v_field[j, i]) / dy_loc
        return div

    # Time loop
    print('Starting time-stepping...')
    report_every = max(1, Nsteps // 10)
    for step in range(Nsteps):
        t = (step + 1) * dt

        # Compute advective terms for u and v (explicit). Keep loops for clarity; small grid.
        adv_u = np.zeros_like(u)
        adv_v = np.zeros_like(v)

        # adv_u
        for j in range(Ny):
            for i in range(Nx + 1):
                ui = u[j, i]
                # du/dx
                if i == 0:
                    du_dx = (u[j, i + 1] - u[j, i]) / dx
                elif i == Nx:
                    du_dx = (u[j, i] - u[j, i - 1]) / dx
                else:
                    du_dx = (u[j, i + 1] - u[j, i - 1]) / (2.0 * dx)
                # du/dy
                if j == 0:
                    du_dy = (u[j + 1, i] - u[j, i]) / dy
                elif j == Ny - 1:
                    du_dy = (u[j, i] - u[j - 1, i]) / dy
                else:
                    du_dy = (u[j + 1, i] - u[j - 1, i]) / (2.0 * dy)
                # v at u location: average surrounding v nodes
                vvals = []
                for jj in [j, j + 1]:
                    for ii in [i - 1, i]:
                        if 0 <= jj <= Ny and 0 <= ii <= Nx - 1:
                            vvals.append(v[jj, ii])
                v_at_u = np.mean(vvals) if vvals else 0.0
                adv_u[j, i] = ui * du_dx + v_at_u * du_dy

        # adv_v
        for j in range(Ny + 1):
            for i in range(Nx):
                vi = v[j, i]
                # dv/dx
                if i == 0:
                    dv_dx = (v[j, i + 1] - v[j, i]) / dx
                elif i == Nx - 1:
                    dv_dx = (v[j, i] - v[j, i - 1]) / dx
                else:
                    dv_dx = (v[j, i + 1] - v[j, i - 1]) / (2.0 * dx)
                # dv/dy
                if j == 0:
                    dv_dy = (v[j + 1, i] - v[j, i]) / dy
                elif j == Ny:
                    dv_dy = (v[j, i] - v[j - 1, i]) / dy
                else:
                    dv_dy = (v[j + 1, i] - v[j - 1, i]) / (2.0 * dy)
                # u at v location: average surrounding u nodes
                uvals = []
                for jj in [j - 1, j]:
                    for ii in [i, i + 1]:
                        if 0 <= jj <= Ny - 1 and 0 <= ii <= Nx:
                            uvals.append(u[jj, ii])
                u_at_v = np.mean(uvals) if uvals else 0.0
                adv_v[j, i] = u_at_v * dv_dx + vi * dv_dy

        # Flatten adv arrays
        adv_u_flat = adv_u.ravel()
        adv_v_flat = adv_v.ravel()

        # AB2 for advection: adv_term = 1.5*adv^n - 0.5*adv^{n-1}
        if step == 0:
            adv_u_eff = adv_u_flat
            adv_v_eff = adv_v_flat
        else:
            adv_u_eff = 1.5 * adv_u_flat - 0.5 * adv_u_old
            adv_v_eff = 1.5 * adv_v_flat - 0.5 * adv_v_old

        adv_u_old = adv_u_flat.copy()
        adv_v_old = adv_v_flat.copy()

        # Build RHS for u and v
        u_flat = u.ravel()
        v_flat = v.ravel()
        RHS_u = (1.0 / dt) * u_flat + (visc / 2.0) * (L_u.dot(u_flat)) - adv_u_eff
        RHS_v = (1.0 / dt) * v_flat + (visc / 2.0) * (L_v.dot(v_flat)) - adv_v_eff

        # Enforce Dirichlet values in RHS for u and v consistent with identity rows in A
        RHS_u[u_dir_indices] = U_in
        RHS_v[v_dir_indices] = 0.0

        # Solve for provisional velocities u_star, v_star
        u_star_flat = lu_u.solve(RHS_u)
        v_star_flat = lu_v.solve(RHS_v)

        u_star = u_star_flat.reshape(u.shape)
        v_star = v_star_flat.reshape(v.shape)

        # Solve Poisson: A_p * p_flat = (1/dt) * div(u_star, v_star) with pressure fix already built into A_p
        div_star = divergence_field(u_star, v_star, dx, dy, Nx, Ny)
        rhs_p = (1.0 / dt) * div_star.ravel()
        rhs_p[fix_idx] = 0.0
        p_flat = lu_p.solve(rhs_p)
        p = p_flat.reshape((Ny, Nx))

        # Project velocities: u^{n+1} = u_star - dt * grad_x(p)
        u_new = u_star.copy()
        v_new = v_star.copy()

        # grad_x at u faces: use one-sided zero-gradient consistent with Poisson Neumann treatment
        for j in range(Ny):
            for i in range(Nx + 1):
                if u_dirichlet_mask[j, i]:
                    u_new[j, i] = U_in
                    continue
                if i == 0 or i == Nx:
                    # zero normal pressure gradient -> dpdx = 0
                    dpdx = 0.0
                else:
                    dpdx = (p[j, i] - p[j, i - 1]) / dx
                u_new[j, i] = u_star[j, i] - dt * dpdx

        # grad_y at v faces: similar treatment
        for j in range(Ny + 1):
            for i in range(Nx):
                if v_dirichlet_mask[j, i]:
                    v_new[j, i] = 0.0
                    continue
                if j == 0 or j == Ny:
                    dpdy = 0.0
                else:
                    dpdy = (p[j, i] - p[j - 1, i]) / dy
                v_new[j, i] = v_star[j, i] - dt * dpdy

        # Overwrite velocities
        u = u_new
        v = v_new

        # Diagnostics and progress prints only every 10% of steps
        if (step + 1) % report_every == 0 or (step == Nsteps - 1):
            div_norm = np.max(np.abs(divergence_field(u, v, dx, dy, Nx, Ny)))
            umax_now = np.max(np.abs(u))
            vmax_now = np.max(np.abs(v))
            print(f'Step {step+1}/{Nsteps}, time {t:.3f}, umax {umax_now:.4f}, vmax {vmax_now:.4f}, max|div| {div_norm:.3e}')

    # Postprocessing: compute velocity magnitude at cell centers by averaging surrounding u/v
    Uc = np.zeros((Ny, Nx))
    Vc = np.zeros((Ny, Nx))
    for j in range(Ny):
        for i in range(Nx):
            Uc[j, i] = 0.5 * (u[j, i] + u[j, i + 1])
            Vc[j, i] = 0.5 * (v[j, i] + v[j + 1, i])
    speed = np.sqrt(Uc ** 2 + Vc ** 2)

    # Print a simple quantitative metric
    mean_speed = np.mean(speed)
    max_speed = np.max(speed)
    print(f'Postprocess: mean speed={mean_speed:.5f}, max speed={max_speed:.5f}')

    # Plot contour of velocity magnitude with cylinder
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = (np.arange(Nx) + 0.5) * dx
    y = (np.arange(Ny) + 0.5) * dy
    X, Y = np.meshgrid(x, y)
    cs = ax.contourf(X, Y, speed, levels=50, cmap='RdBu_r')
    cbar = fig.colorbar(cs)
    cbar.set_label('Velocity magnitude')
    # plot cylinder outline
    theta = np.linspace(0, 2 * np.pi, 200)
    xc = cx + radius * np.cos(theta)
    yc = cy + radius * np.sin(theta)
    ax.plot(xc, yc, 'k', linewidth=1.2)
    ax.set_aspect('equal')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_title(f'Velocity magnitude at t={T:.2f} (Re={Re:.0f})')
    plt.tight_layout()
    fig.savefig('velocity_magnitude.png', dpi=150)
    print('Saved figure to velocity_magnitude.png')


if __name__ == '__main__':
    main()
```


#### Script block3:

```python
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import ceil

# Global viscosity (use this name exactly as requested)
Re = 100.0
visc = 1.0 / Re

# --------------------------------------------------------------------------------
# Helper functions for grid/mask/operators
# --------------------------------------------------------------------------------

def build_grid(x1, x2, y1, y2, points_per_diameter):
    """Build MAC grid geometry and coordinates.
    Returns a dict with grid spacings and coordinates.
    """
    Lx = x1 + x2
    Ly = y1 + y2
    dx = 1.0 / points_per_diameter  # diameter normalized to 1
    Nx = int(round(Lx / dx))
    Ny = int(round(Ly / dx))
    dx = Lx / Nx
    dy = Ly / Ny
    x_u = np.linspace(0.0, Lx, Nx + 1)  # u at vertical faces
    y_u = (np.arange(Ny) + 0.5) * dy
    x_v = (np.arange(Nx) + 0.5) * dx
    y_v = np.linspace(0.0, Ly, Ny + 1)  # v at horizontal faces
    x_p = (np.arange(Nx) + 0.5) * dx
    y_p = (np.arange(Ny) + 0.5) * dy
    return {
        'Lx': Lx, 'Ly': Ly, 'dx': dx, 'dy': dy,
        'Nx': Nx, 'Ny': Ny,
        'x_u': x_u, 'y_u': y_u,
        'x_v': x_v, 'y_v': y_v,
        'x_p': x_p, 'y_p': y_p
    }


def build_cylinder_mask(Nx, Ny, dx, dy, center_x, center_y, radius):
    """Cell-centered solid mask chi on p-grid (Ny, Nx)."""
    x = (np.arange(Nx) + 0.5) * dx
    y = (np.arange(Ny) + 0.5) * dy
    Xc, Yc = np.meshgrid(x, y)
    r2 = (Xc - center_x) ** 2 + (Yc - center_y) ** 2
    chi = (r2 <= radius ** 2).astype(float)
    return chi


def interp_chi_to_u(chi):
    # chi shape (Ny, Nx); u shape (Ny, Nx+1)
    Ny, Nx = chi.shape
    chi_u = np.zeros((Ny, Nx + 1))
    # interior faces: average left and right cell chi
    chi_u[:, 1:Nx] = 0.5 * (chi[:, 0:Nx - 1] + chi[:, 1:Nx])
    # leftmost face i=0: only right cell
    chi_u[:, 0] = chi[:, 0]
    # rightmost face i=Nx: only left cell
    chi_u[:, Nx] = chi[:, Nx - 1]
    return chi_u


def interp_chi_to_v(chi):
    # v shape (Ny+1, Nx)
    Ny, Nx = chi.shape
    chi_v = np.zeros((Ny + 1, Nx))
    # interior faces: average top and bottom cells
    chi_v[1:Ny, :] = 0.5 * (chi[0:Ny - 1, :] + chi[1:Ny, :])
    # bottom j=0: only top cell j=0
    chi_v[0, :] = chi[0, :]
    # top j=Ny: only bottom cell j=Ny-1
    chi_v[Ny, :] = chi[Ny - 1, :]
    return chi_v


def node_index_u(i, j, Nx, Ny):
    return j * (Nx + 1) + i


def node_index_v(i, j, Nx, Ny):
    return j * Nx + i


def node_index_p(i, j, Nx, Ny):
    return j * Nx + i


def build_laplacian_staggered_u(Nx, Ny, dx, dy, bcs):
    """Assemble 2D Laplacian on u-grid. Correct sign convention: diag = -sum(neighbors).
    Dirichlet rows are set as identity (1) here to simplify enforcement; caller will
    replace rows with identity semantics as needed.
    """
    Nu = Ny * (Nx + 1)
    rows = []
    cols = []
    data = []

    for j in range(Ny):
        for i in range(Nx + 1):
            idx = node_index_u(i, j, Nx, Ny)
            if bcs['u_dirichlet_mask'][j, i]:
                rows.append(idx); cols.append(idx); data.append(1.0)
                continue
            center_coef = 0.0
            # x-direction neighbors
            if i - 1 >= 0:
                nidx = node_index_u(i - 1, j, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dx ** 2)
                center_coef += 1.0 / dx ** 2
            else:
                # Neumann treated as one-sided -> contributes same center accumulation
                if bcs['u_x_neumann_left']:
                    center_coef += 1.0 / dx ** 2
            if i + 1 <= Nx:
                nidx = node_index_u(i + 1, j, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dx ** 2)
                center_coef += 1.0 / dx ** 2
            else:
                if bcs['u_x_neumann_right']:
                    center_coef += 1.0 / dx ** 2
            # y-direction neighbors
            if j - 1 >= 0:
                nidx = node_index_u(i, j - 1, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dy ** 2)
                center_coef += 1.0 / dy ** 2
            else:
                if bcs['u_y_neumann_bottom']:
                    center_coef += 1.0 / dy ** 2
            if j + 1 <= Ny - 1:
                nidx = node_index_u(i, j + 1, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dy ** 2)
                center_coef += 1.0 / dy ** 2
            else:
                if bcs['u_y_neumann_top']:
                    center_coef += 1.0 / dy ** 2
            # diagonal is negative sum of neighbor coefficients
            rows.append(idx); cols.append(idx); data.append(-center_coef)
    L = sp.csr_matrix((data, (rows, cols)), shape=(Nu, Nu))
    return L


def build_laplacian_staggered_v(Nx, Ny, dx, dy, bcs):
    Nv = (Ny + 1) * Nx
    rows = []
    cols = []
    data = []
    for j in range(Ny + 1):
        for i in range(Nx):
            idx = node_index_v(i, j, Nx, Ny)
            if bcs['v_dirichlet_mask'][j, i]:
                rows.append(idx); cols.append(idx); data.append(1.0)
                continue
            center_coef = 0.0
            # x-direction neighbors
            if i - 1 >= 0:
                nidx = node_index_v(i - 1, j, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dx ** 2)
                center_coef += 1.0 / dx ** 2
            else:
                if bcs['v_x_neumann_left']:
                    center_coef += 1.0 / dx ** 2
            if i + 1 <= Nx - 1:
                nidx = node_index_v(i + 1, j, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dx ** 2)
                center_coef += 1.0 / dx ** 2
            else:
                if bcs['v_x_neumann_right']:
                    center_coef += 1.0 / dx ** 2
            # y-direction neighbors
            if j - 1 >= 0:
                nidx = node_index_v(i, j - 1, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dy ** 2)
                center_coef += 1.0 / dy ** 2
            else:
                if bcs['v_y_neumann_bottom']:
                    center_coef += 1.0 / dy ** 2
            if j + 1 <= Ny:
                nidx = node_index_v(i, j + 1, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dy ** 2)
                center_coef += 1.0 / dy ** 2
            else:
                if bcs['v_y_neumann_top']:
                    center_coef += 1.0 / dy ** 2
            rows.append(idx); cols.append(idx); data.append(-center_coef)
    L = sp.csr_matrix((data, (rows, cols)), shape=(Nv, Nv))
    return L


def build_poisson_p(Nx, Ny, dx, dy):
    """Assemble pressure Poisson matrix with correct sign convention."""
    Np = Nx * Ny
    rows = []
    cols = []
    data = []
    for j in range(Ny):
        for i in range(Nx):
            idx = node_index_p(i, j, Nx, Ny)
            center_coef = 0.0
            # left
            if i - 1 >= 0:
                nidx = node_index_p(i - 1, j, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dx ** 2)
                center_coef += 1.0 / dx ** 2
            else:
                # Neumann -> one-sided contributes same accumulation
                center_coef += 1.0 / dx ** 2
            # right
            if i + 1 <= Nx - 1:
                nidx = node_index_p(i + 1, j, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dx ** 2)
                center_coef += 1.0 / dx ** 2
            else:
                center_coef += 1.0 / dx ** 2
            # bottom
            if j - 1 >= 0:
                nidx = node_index_p(i, j - 1, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dy ** 2)
                center_coef += 1.0 / dy ** 2
            else:
                center_coef += 1.0 / dy ** 2
            # top
            if j + 1 <= Ny - 1:
                nidx = node_index_p(i, j + 1, Nx, Ny)
                rows.append(idx); cols.append(nidx); data.append(1.0 / dy ** 2)
                center_coef += 1.0 / dy ** 2
            else:
                center_coef += 1.0 / dy ** 2
            rows.append(idx); cols.append(idx); data.append(-center_coef)
    A = sp.csr_matrix((data, (rows, cols)), shape=(Np, Np))
    return A


# --------------------------------------------------------------------------------
# Main solver
# --------------------------------------------------------------------------------

def main():
    # Problem geometry and parameters (explicitly passed)
    d = 1.0
    x1 = 3.0
    x2 = 8.0
    y1 = 3.0
    y2 = 3.0
    Lx = x1 + x2
    Ly = y1 + y2
    cx = x1  # cylinder center x
    cy = Ly / 2.0
    radius = d / 2.0

    # grid resolution (points per diameter) - moderate
    ppd = 12  # points per diameter -> reasonable runtime and visible wake
    grid = build_grid(x1, x2, y1, y2, ppd)
    Nx = grid['Nx']; Ny = grid['Ny']; dx = grid['dx']; dy = grid['dy']

    # penalization parameter eta (start with a modestly larger value for stability)
    eta = 1e-2

    # Build masks and interpolations
    chi = build_cylinder_mask(Nx, Ny, dx, dy, cx, cy, radius)
    chi_u = interp_chi_to_u(chi)
    chi_v = interp_chi_to_v(chi)

    Nu = Ny * (Nx + 1)
    Nv = (Ny + 1) * Nx
    assert chi_u.size == Nu, f"chi_u.size {chi_u.size} != Nu {Nu}"
    assert chi_v.size == Nv, f"chi_v.size {chi_v.size} != Nv {Nv}"

    alpha_u = (chi_u / eta).ravel()
    alpha_v = (chi_v / eta).ravel()

    # boundary condition masks for assembling Laplacians
    u_dirichlet_mask = np.zeros((Ny, Nx + 1), dtype=bool)
    u_dirichlet_mask[:, 0] = True  # inlet u=U_in
    bcs_u = {
        'u_dirichlet_mask': u_dirichlet_mask,
        'u_x_neumann_left': False,
        'u_x_neumann_right': True,
        'u_y_neumann_bottom': True,
        'u_y_neumann_top': True
    }

    v_dirichlet_mask = np.zeros((Ny + 1, Nx), dtype=bool)
    v_dirichlet_mask[0, :] = True  # bottom wall v=0
    v_dirichlet_mask[-1, :] = True  # top wall v=0
    v_dirichlet_mask[:, 0] = True  # inlet v=0
    bcs_v = {
        'v_dirichlet_mask': v_dirichlet_mask,
        'v_x_neumann_left': False,
        'v_x_neumann_right': True,
        'v_y_neumann_bottom': False,
        'v_y_neumann_top': False
    }

    # Build Laplacians once
    print('Assembling matrices...')
    L_u = build_laplacian_staggered_u(Nx, Ny, dx, dy, bcs_u)
    L_v = build_laplacian_staggered_v(Nx, Ny, dx, dy, bcs_v)

    # Poisson matrix
    A_p = build_poisson_p(Nx, Ny, dx, dy)

    # Initial condition: u=1 everywhere (inlet) and v=0
    u = np.ones((Ny, Nx + 1))
    v = np.zeros((Ny + 1, Nx))

    # Pressure initial guess
    p = np.zeros((Ny, Nx))

    # Time stepping parameters
    U_in = 1.0
    umax = U_in
    CFL = 0.5
    dt_CFL = CFL * min(dx / max(1e-8, umax), dy / max(1e-8, umax))
    dt_visc = 0.5 * min(dx * dx, dy * dy) / max(1e-12, visc)
    dt = min(dt_CFL, dt_visc)
    dt = min(dt, 0.02)
    T = 20.0
    Nsteps = int(ceil(T / dt))
    dt = T / Nsteps
    print(f'Grid: Nx={Nx}, Ny={Ny}, dx={dx:.4f}, dy={dy:.4f}, dt={dt:.5f}, steps={Nsteps}')

    # Build final momentum matrices with dt and penalization included
    Iu = sp.eye(Nu, format='csr')
    Iv = sp.eye(Nv, format='csr')
    # Use Crank-Nicolson for viscous (semi-implicit) and implicit penalization on diagonal
    A_u = (1.0 / dt) * Iu - (visc / 2.0) * L_u + sp.diags(alpha_u)
    A_v = (1.0 / dt) * Iv - (visc / 2.0) * L_v + sp.diags(alpha_v)

    # Enforce Dirichlet rows on full momentum matrices before factorization
    print('Applying Dirichlet enforcement to momentum matrices...')
    A_u = A_u.tolil()
    A_v = A_v.tolil()

    u_dir_indices = np.where(u_dirichlet_mask.ravel())[0]
    v_dir_indices = np.where(v_dirichlet_mask.ravel())[0]
    assert u_dir_indices.size > 0, 'Expected some u Dirichlet indices'
    assert v_dir_indices.size > 0, 'Expected some v Dirichlet indices'

    for idx in u_dir_indices:
        A_u.rows[idx] = [idx]
        A_u.data[idx] = [1.0]
    for idx in v_dir_indices:
        A_v.rows[idx] = [idx]
        A_v.data[idx] = [1.0]

    A_u = A_u.tocsc()
    A_v = A_v.tocsc()

    # Factorize momentum matrices
    print('Factorizing momentum matrices (one-time)...')
    try:
        lu_u = spla.splu(A_u)
        lu_v = spla.splu(A_v)
    except Exception as e:
        print('Momentum factorization failed:', e)
        return

    # Quick unit-check
    test_rhs_u = np.zeros(Nu)
    test_rhs_v = np.zeros(Nv)
    test_rhs_u[u_dir_indices] = U_in
    test_rhs_v[v_dir_indices] = 0.0
    sol_u_test = lu_u.solve(test_rhs_u)
    sol_v_test = lu_v.solve(test_rhs_v)
    if not (np.allclose(sol_u_test[u_dir_indices], U_in, atol=1e-12) and np.allclose(sol_v_test[v_dir_indices], 0.0, atol=1e-12)):
        print('Dirichlet enforcement unit-check failed; inspecting norms: ',
              np.max(np.abs(sol_u_test[u_dir_indices] - U_in)), np.max(np.abs(sol_v_test[v_dir_indices] - 0.0)))
        return

    # Factorize Poisson after pressure fix; pick a pressure fix in the fluid region near outlet
    Np = Nx * Ny
    # search for a fluid cell near outlet center (right side)
    fix_idx = None
    for ii in range(Nx - 1, max(Nx - 6, 0), -1):
        jj = Ny // 2
        if chi[jj, ii] < 0.5:
            fix_idx = node_index_p(ii, jj, Nx, Ny)
            break
    if fix_idx is None:
        # fallback: first fluid cell found
        found = False
        for jj in range(Ny):
            for ii in range(Nx):
                if chi[jj, ii] < 0.5:
                    fix_idx = node_index_p(ii, jj, Nx, Ny)
                    found = True
                    break
            if found:
                break
    assert fix_idx is not None, 'Could not find a fluid cell for pressure fixing.'

    A_p = A_p.tolil()
    A_p[fix_idx, :] = 0.0
    A_p[fix_idx, fix_idx] = 1.0
    A_p = A_p.tocsc()
    print('Factorizing Poisson matrix (one-time)...')
    try:
        lu_p = spla.splu(A_p)
    except Exception as e:
        print('Poisson factorization failed:', e)
        return

    # Prepare AB2 storage
    adv_u_old = np.zeros_like(u).ravel()
    adv_v_old = np.zeros_like(v).ravel()

    # diagnostics
    def divergence_field(u_field, v_field, dx_loc, dy_loc, Nx_loc, Ny_loc):
        div = np.zeros((Ny_loc, Nx_loc))
        for j in range(Ny_loc):
            for i in range(Nx_loc):
                div[j, i] = (u_field[j, i + 1] - u_field[j, i]) / dx_loc + (v_field[j + 1, i] - v_field[j, i]) / dy_loc
        return div

    def check_numerics(name, *arrays):
        for arr in arrays:
            if np.any(np.isnan(arr)) or np.any(np.isinf(arr)):
                print(f'Numeric problem detected in {name}: NaN/Inf found. Stopping.')
                return False
        return True

    # Time loop
    print('Starting time-stepping...')
    report_every = max(1, Nsteps // 10)
    for step in range(Nsteps):
        t = (step + 1) * dt

        # Compute advective terms with first-order upwind for stability
        adv_u = np.zeros_like(u)
        adv_v = np.zeros_like(v)

        # adv_u (u at faces)
        for j in range(Ny):
            for i in range(Nx + 1):
                ui = u[j, i]
                # du/dx: upwind based on ui
                if i == 0:
                    du_dx = (u[j, i + 1] - u[j, i]) / dx
                elif i == Nx:
                    du_dx = (u[j, i] - u[j, i - 1]) / dx
                else:
                    if ui >= 0.0:
                        du_dx = (u[j, i] - u[j, i - 1]) / dx
                    else:
                        du_dx = (u[j, i + 1] - u[j, i]) / dx
                # du/dy: upwind based on v_at_u
                if j == 0:
                    du_dy = (u[j + 1, i] - u[j, i]) / dy
                elif j == Ny - 1:
                    du_dy = (u[j, i] - u[j - 1, i]) / dy
                else:
                    # compute v at u location
                    vvals = []
                    for jj in [j, j + 1]:
                        for ii in [i - 1, i]:
                            if 0 <= jj <= Ny and 0 <= ii <= Nx - 1:
                                vvals.append(v[jj, ii])
                    v_at_u = np.mean(vvals) if vvals else 0.0
                    if v_at_u >= 0.0:
                        du_dy = (u[j, i] - u[j - 1, i]) / dy
                    else:
                        du_dy = (u[j + 1, i] - u[j, i]) / dy
                # v at u location
                vvals = []
                for jj in [j, j + 1]:
                    for ii in [i - 1, i]:
                        if 0 <= jj <= Ny and 0 <= ii <= Nx - 1:
                            vvals.append(v[jj, ii])
                v_at_u = np.mean(vvals) if vvals else 0.0
                adv_u[j, i] = ui * du_dx + v_at_u * du_dy

        # adv_v (v at faces)
        for j in range(Ny + 1):
            for i in range(Nx):
                vi = v[j, i]
                # dv/dx: upwind based on u_at_v
                # u at v location: average surrounding u nodes
                uvals = []
                for jj in [j - 1, j]:
                    for ii in [i, i + 1]:
                        if 0 <= jj <= Ny - 1 and 0 <= ii <= Nx:
                            uvals.append(u[jj, ii])
                u_at_v = np.mean(uvals) if uvals else 0.0
                if i == 0:
                    dv_dx = (v[j, i + 1] - v[j, i]) / dx
                elif i == Nx - 1:
                    dv_dx = (v[j, i] - v[j, i - 1]) / dx
                else:
                    if u_at_v >= 0.0:
                        dv_dx = (v[j, i] - v[j, i - 1]) / dx
                    else:
                        dv_dx = (v[j, i + 1] - v[j, i]) / dx
                # dv/dy: upwind based on vi
                if j == 0:
                    dv_dy = (v[j + 1, i] - v[j, i]) / dy
                elif j == Ny:
                    dv_dy = (v[j, i] - v[j - 1, i]) / dy
                else:
                    if vi >= 0.0:
                        dv_dy = (v[j, i] - v[j - 1, i]) / dy
                    else:
                        dv_dy = (v[j + 1, i] - v[j, i]) / dy
                adv_v[j, i] = u_at_v * dv_dx + vi * dv_dy

        # Flatten adv arrays
        adv_u_flat = adv_u.ravel()
        adv_v_flat = adv_v.ravel()

        # AB2 for advection
        if step == 0:
            adv_u_eff = adv_u_flat
            adv_v_eff = adv_v_flat
        else:
            adv_u_eff = 1.5 * adv_u_flat - 0.5 * adv_u_old
            adv_v_eff = 1.5 * adv_v_flat - 0.5 * adv_v_old

        adv_u_old = adv_u_flat.copy()
        adv_v_old = adv_v_flat.copy()

        # Basic numeric checks
        if not check_numerics('adv', adv_u_flat, adv_v_flat):
            break

        # Build RHS for u and v: note penalization is implicit in A (diagonal), so RHS uses (1/dt)*u + visc/2 * L u - adv
        u_flat = u.ravel()
        v_flat = v.ravel()
        RHS_u = (1.0 / dt) * u_flat + (visc / 2.0) * (L_u.dot(u_flat)) - adv_u_eff
        RHS_v = (1.0 / dt) * v_flat + (visc / 2.0) * (L_v.dot(v_flat)) - adv_v_eff

        # Enforce Dirichlet values in RHS consistent with identity rows in A
        RHS_u[u_dir_indices] = U_in
        RHS_v[v_dir_indices] = 0.0

        if not check_numerics('RHS_momentum', RHS_u, RHS_v):
            break

        # Solve for provisional velocities u_star, v_star
        try:
            u_star_flat = lu_u.solve(RHS_u)
            v_star_flat = lu_v.solve(RHS_v)
        except Exception as e:
            print('Momentum solve failed at step', step, 'error:', e)
            break

        if not check_numerics('u_star/v_star', u_star_flat, v_star_flat):
            break

        u_star = u_star_flat.reshape(u.shape)
        v_star = v_star_flat.reshape(v.shape)

        # Solve Poisson: A_p * p_flat = (1/dt) * div(u_star, v_star)
        div_star = divergence_field(u_star, v_star, dx, dy, Nx, Ny)
        rhs_p = (1.0 / dt) * div_star.ravel()
        # Enforce pressure-fix rhs value
        rhs_p[fix_idx] = 0.0

        if not check_numerics('rhs_p', rhs_p):
            break
        try:
            p_flat = lu_p.solve(rhs_p)
        except Exception as e:
            print('Poisson solve failed at step', step, 'error:', e)
            break
        if not check_numerics('p', p_flat):
            break

        p = p_flat.reshape((Ny, Nx))

        # Project velocities
        u_new = u_star.copy()
        v_new = v_star.copy()

        for j in range(Ny):
            for i in range(Nx + 1):
                if u_dirichlet_mask[j, i]:
                    u_new[j, i] = U_in
                    continue
                if i == 0 or i == Nx:
                    dpdx = 0.0
                else:
                    dpdx = (p[j, i] - p[j, i - 1]) / dx
                u_new[j, i] = u_star[j, i] - dt * dpdx

        for j in range(Ny + 1):
            for i in range(Nx):
                if v_dirichlet_mask[j, i]:
                    v_new[j, i] = 0.0
                    continue
                if j == 0 or j == Ny:
                    dpdy = 0.0
                else:
                    dpdy = (p[j, i] - p[j - 1, i]) / dy
                v_new[j, i] = v_star[j, i] - dt * dpdy

        # Overwrite velocities and check numerics
        u = u_new
        v = v_new
        if not check_numerics('u/v after proj', u, v):
            break

        # Diagnostics and progress prints only every 10% of steps
        if (step + 1) % report_every == 0 or (step == Nsteps - 1):
            div_norm = np.max(np.abs(divergence_field(u, v, dx, dy, Nx, Ny)))
            umax_now = np.max(np.abs(u))
            vmax_now = np.max(np.abs(v))
            print(f'Step {step+1}/{Nsteps}, time {t:.3f}, umax {umax_now:.4f}, vmax {vmax_now:.4f}, max|div| {div_norm:.3e}')

    # Postprocessing: compute velocity magnitude at cell centers by averaging surrounding u/v
    Uc = np.zeros((Ny, Nx))
    Vc = np.zeros((Ny, Nx))
    for j in range(Ny):
        for i in range(Nx):
            Uc[j, i] = 0.5 * (u[j, i] + u[j, i + 1])
            Vc[j, i] = 0.5 * (v[j, i] + v[j + 1, i])
    speed = np.sqrt(Uc ** 2 + Vc ** 2)

    # Print a simple quantitative metric, and check speeds within solid
    mean_speed = np.nanmean(speed)
    max_speed = np.nanmax(speed)
    # verify penalization effect: max deviation inside solid from solid velocity (u=0, v=0)
    solid_cells = chi > 0.5
    if np.any(solid_cells):
        dev_in_solid = np.sqrt((Uc[solid_cells] - 0.0) ** 2 + (Vc[solid_cells] - 0.0) ** 2)
        max_dev_solid = np.max(dev_in_solid)
    else:
        max_dev_solid = 0.0
    print(f'Postprocess: mean speed={mean_speed:.5f}, max speed={max_speed:.5f}, max_dev_in_solid={max_dev_solid:.5e}')

    # Plot contour of velocity magnitude with cylinder
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = (np.arange(Nx) + 0.5) * dx
    y = (np.arange(Ny) + 0.5) * dy
    X, Y = np.meshgrid(x, y)
    cs = ax.contourf(X, Y, speed, levels=50, cmap='RdBu_r')
    cbar = fig.colorbar(cs)
    cbar.set_label('Velocity magnitude')
    # plot cylinder outline
    theta = np.linspace(0, 2 * np.pi, 200)
    xc = cx + radius * np.cos(theta)
    yc = cy + radius * np.sin(theta)
    ax.plot(xc, yc, 'k', linewidth=1.2)
    ax.set_aspect('equal')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_title(f'Velocity magnitude at t={T:.2f} (Re={Re:.0f})')
    plt.tight_layout()
    fig.savefig('velocity_magnitude.png', dpi=150)
    print('Saved figure to velocity_magnitude.png')


if __name__ == '__main__':
    main()
```




### runtime_outputs

#### Output block1


Runtime Error:
Traceback (most recent call last):
  File "C:\Users\admin\AppData\Local\Temp\ipykernel_49500\2663542217.py", line 138, in execute_code_tool
    exec(code, global_scope, global_scope)
  File "<string>", line 563, in <module>
  File "<string>", line 266, in main
  File "<string>", line 55, in interp_chi_to_u
ValueError: operands could not be broadcast together with shapes (72,131) (72,132) 



#### Output block2

Assembling matrices...
Grid: Nx=132, Ny=72, dx=0.0833, dy=0.0833, dt=0.02000, steps=1000
Applying Dirichlet enforcement to momentum matrices...
Factorizing momentum matrices (one-time)...
Factorizing Poisson matrix (one-time)...
Starting time-stepping...
<string>:431: RuntimeWarning: overflow encountered in scalar multiply
<string>:431: RuntimeWarning: invalid value encountered in scalar add
<string>:458: RuntimeWarning: overflow encountered in scalar multiply
<string>:458: RuntimeWarning: invalid value encountered in scalar add
Step 100/1000, time 2.000, umax nan, vmax nan, max|div| nan
Step 200/1000, time 4.000, umax nan, vmax nan, max|div| nan
Step 300/1000, time 6.000, umax nan, vmax nan, max|div| nan
Step 400/1000, time 8.000, umax nan, vmax nan, max|div| nan
Step 500/1000, time 10.000, umax nan, vmax nan, max|div| nan
Step 600/1000, time 12.000, umax nan, vmax nan, max|div| nan
Step 700/1000, time 14.000, umax nan, vmax nan, max|div| nan
Step 800/1000, time 16.000, umax nan, vmax nan, max|div| nan
Step 900/1000, time 18.000, umax nan, vmax nan, max|div| nan
Step 1000/1000, time 20.000, umax nan, vmax nan, max|div| nan
Postprocess: mean speed=nan, max speed=nan
Saved figure to velocity_magnitude.png



#### Output block3

Assembling matrices...
Grid: Nx=132, Ny=72, dx=0.0833, dy=0.0833, dt=0.02000, steps=1000
Applying Dirichlet enforcement to momentum matrices...
Factorizing momentum matrices (one-time)...
Factorizing Poisson matrix (one-time)...
Starting time-stepping...
Step 100/1000, time 2.000, umax 1.3260, vmax 0.5878, max|div| 6.398e+00
Step 200/1000, time 4.000, umax 1.2965, vmax 0.5748, max|div| 6.605e+00
Step 300/1000, time 6.000, umax 1.2757, vmax 0.5655, max|div| 6.704e+00
Step 400/1000, time 8.000, umax 1.2616, vmax 0.5600, max|div| 6.820e+00
Step 500/1000, time 10.000, umax 1.2515, vmax 0.5567, max|div| 7.014e+00
Step 600/1000, time 12.000, umax 1.2436, vmax 0.5545, max|div| 7.198e+00
Step 700/1000, time 14.000, umax 1.2369, vmax 0.5530, max|div| 7.359e+00
Step 800/1000, time 16.000, umax 1.2312, vmax 0.5519, max|div| 7.502e+00
Step 900/1000, time 18.000, umax 1.2264, vmax 0.5510, max|div| 7.622e+00
Step 1000/1000, time 20.000, umax 1.2226, vmax 0.5503, max|div| 7.711e+00
Postprocess: mean speed=0.81858, max speed=1.22308, max_dev_in_solid=9.08114e-02
Saved figure to velocity_magnitude.png




### review_decision
revise

### review_comments
	Current Stage [A/3]
1. MAJOR

Category: runtime
Issue: Broadcast error in interp_chi_to_u: the interior-face averaging uses mismatched slices 0:Nx-1 and 0:Nx leading to shapes (Ny,Nx-1) and (Ny,Nx) and a ValueError at runtime.

Recommendation: Fix the slice for the right-hand cell when averaging. Change the line
  chi_u[:, 1:Nx] = 0.5 * (chi[:, 0:Nx - 1] + chi[:, 0:Nx])
to
  chi_u[:, 1:Nx] = 0.5 * (chi[:, 0:Nx - 1] + chi[:, 1:Nx])
Then verify that chi_u.ravel() length equals Nu = Ny*(Nx+1) before using alpha_u. Add an assertion after interpolation, e.g. assert chi_u.size == Nu, to catch similar mismatches early.


2. MAJOR

Category: accuracy
Issue: Dirichlet boundary rows are only set in the Laplacian L_u/L_v but not enforced after assembling the full momentum matrices A_u/A_v. Because you add (1/dt)I and diag(alpha) to L_u/L_v, those Dirichlet rows are no longer exact identity rows in A_u/A_v. Later you treat their RHS as if the row were identity (RHS = value/dt), which is inconsistent and will produce incorrect BC enforcement and wrong solutions.

Recommendation: Enforce Dirichlet rows on the full momentum matrices before factorizing (and before solving) by zeroing the off-diagonal entries and setting the diagonal to 1.0 (or to a large coefficient consistent with RHS scaling), then set the corresponding RHS to the target value (e.g. U_in/dt if you keep (1/dt)*I structure). Concretely:
 - After building A_u and A_v, get the list of dirichlet indices from the masks.
 - For each dirichlet index idx: set A[idx, :] = 0 and A[idx, idx] = 1.0.
 - When building RHS, set RHS[idx] = desired_value (not scaled by 1/dt if you set diagonal=1).
 - Factorize the modified A matrices afterwards.
Also ensure you apply the same enforcement consistently for v-dirichlet nodes. Do this prior to calling splu(). This is a functional blocker for correctness.


3. MAJOR

Category: structure
Issue: Dirichlet enforcement is applied inconsistently between the Laplacian assembly and the momentum-assembly/factorization timing. In the current code you factorize A_u/A_v before you've guaranteed Dirichlet rows have been made exact in A, which will bake incorrect rows into LU factors.

Recommendation: Move the Dirichlet-row replacement step to immediately after A_u/A_v assembly and before calling spla.splu. Replace full rows with identity rows and then call splu. Also add unit tests (small grid) to check that solving with zero RHS returns the enforced Dirichlet values at the masked indices.


4. MINOR

Category: accuracy
Issue: Projection / pressure gradient at domain faces uses arbitrary 0.0 outside-the-domain values (e.g. dpdx at i==0 uses p[j,0] - 0.0). That is inconsistent and can inject spurious pressure forcing near boundaries. The Poisson operator was assembled assuming Neumann at boundaries, but the gradient step should be consistent with those boundary conditions.

Recommendation: Replace the '0.0 outside' approximations with ghost values consistent with Neumann BCs used in the Poisson assembly. For example, for left/right faces use one-sided gradients using the adjacent cell pressure (i==0: dpdx = (p[j,0] - p[j,0]) / dx = 0) if you intend zero normal pressure gradient, or use extrapolated ghost values if you want to maintain compatibility with the Poisson implementation. Ideally compute gradient at faces via the discrete pressure-to-velocity coupling implied by the divergence stencil (use the same indexing conventions). This will reduce boundary artifacts in the projection.


5. MINOR

Category: structure
Issue: Performance and maintainability: advective term assembly is implemented with Python-level nested loops, and some BC flags/masks are inconsistent (e.g. v_x_neumann_left=True while v_dirichlet_mask[:,0]=True). The code lacks assertions on array shapes and sizes after interpolations and before factorization.

Recommendation: Vectorize the finite-difference and interpolation operations with numpy where possible (this will massively speed up the time loop). Clean up BC flag logic so flags are not contradictory (or document precedence of Dirichlet masks over Neumann flags). Add assertions after building arrays: assert chi_u.size==Nu, assert chi_v.size==Nv, and check that the lists of Dirichlet indices are non-empty where expected. These changes are optional for correctness but will improve runtime and debuggability.



	Current Stage [B/3]
1. MAJOR

Category: accuracy
Issue: Discrete Laplacian / Poisson operators are assembled with the wrong sign (center coefficient ends up positive instead of the negative sum of neighbor coefficients). This makes the momentum and Poisson matrices incorrect (likely indefinite) and is the primary cause of the solver blow-up (overflow / NaN) during time stepping.

Recommendation: Fix the Laplacian/Poisson assembly sign. In every build_*laplacian* and build_poisson_p: accumulate neighbor coefficients with a positive accumulator (e.g. center_coef += 1/dx**2 or +=1/dy**2) and set the diagonal to -center_coef (so that a typical 1D stencil gives diag = -2/dx**2 and off-diagonals = +1/dx**2). After this change, re-check matrix spectra/conditioning and re-run factorization; most of the runtime overflow/NaN issues will disappear once the elliptic operators have the correct sign. Example change pattern:
- replace "center_coef = 0.0; ... center_coef -= 1.0/dx**2; ... rows.append(idx); cols.append(idx); data.append(-center_coef)" with
  "center_coef = 0.0; ... center_coef += 1.0/dx**2; ... rows.append(idx); cols.append(idx); data.append(-center_coef)".
Also apply the same fix to build_poisson_p and to both staggered Laplacian builders.


2. MAJOR

Category: runtime
Issue: The time stepping produces overflow/NaN during the run (as seen in runtime warnings). Even after fixing the Laplacian sign, the explicit central-difference treatment of convection (with AB2) and the immersed-boundary/penalization coupling can produce local instability or ill-conditioned algebraic systems.

Recommendation: (a) Add runtime checks to stop early and report location of NaNs/Infs (e.g. check u_flat/v_flat/adv arrays and the RHS after each major operation). (b) Use a more robust discretization for advection (donor-cell / first-order upwind or a limited higher-order scheme) instead of pure central differences — this greatly improves stability for explicit convective terms. Example: for du_dx use upwind based on face velocity sign. (c) Temporarily increase penalization parameter eta (e.g. 1e-2 to 1e-1) to check solver behavior; once Laplacian sign and advection discretization are fixed, reduce eta to the desired value and monitor conditioning. (d) Consider reducing dt to confirm stability, and add a CFL check that includes effects of penalization/solid region if needed.


3. MAJOR

Category: accuracy
Issue: Inconsistent / incomplete treatment of boundary conditions for pressure (Poisson) and projection step. The Poisson assembly assumes homogeneous Neumann at all boundaries (by subtracting neighbor contributions at boundaries), but the velocity Dirichlet at the inflow and the free-slip BCs must be consistently reflected in the Poisson RHS and in the discrete divergence operator; also the pressure-fix node should be outside the solid mask.

Recommendation: Ensure the Poisson RHS and boundary handling are consistent: (a) incorporate contribution of velocity Dirichlet BCs into the discrete divergence and Poisson RHS (normal velocity at Dirichlet boundaries must be inserted correctly). (b) Verify that the chosen pressure-fix cell (fix_idx) is in the fluid (not inside the penalized/solid region). If it lies inside the cylinder mask the pressure fix is nonphysical and can cause solver issues; pick a cell in the downstream fluid near the outlet. (c) After fixing sign convention of the Poisson matrix, verify the sign of RHS used in lu_p.solve (i.e. correct sign following your definition of Laplacian); add a small validation test (e.g. steady Couette or manufactured solution) to validate Poisson/Projection correctness.


4. MINOR

Category: structure
Issue: Implementation details that hamper debugging and performance: heavy nested Python loops for advection, lack of NaN/infinite checks, and some unused/unclear function arguments (e.g. build_cylinder_mask signature includes x1,y1,x2,y2 but they are unused inside).

Recommendation: (a) Add quick diagnostic checks after each major linear algebra operation (e.g. after factorization, after solve) to detect Inf/NaN early and print the offending norm. (b) Vectorize the advection assembly or use numba to accelerate the tight loops — on larger grids the current Python loops will be slow. (c) Clean up function signatures (remove unused arguments) and add docstrings to clarify indexing conventions and grid origins. (d) Keep the one-time factorization approach (good), but add an option to use permutation reordering (use splu with specified permc_spec) if factorization is slow or fails.


5. MINOR

Category: accuracy
Issue: Penalization strategy and values: using a very small eta (1e-3) and relying on central differencing for convection can make the flow near the cylinder stiff; also the penalization term is added only to the diagonal but the projected velocity update may not fully respect the no-slip if projection and penalization signs are inconsistent.

Recommendation: (a) After fixing Laplacian sign, keep penalization implicit (diagonal) but test with a larger eta first to ensure conditioning; reduce eta once correctness is verified. (b) Add verification that velocities inside the penalized mask tend to the solid velocity (e.g. compute max|u-Usolid| within chi region); (c) consider treating penalization in RHS too (i.e., include alpha * u^n contribution appropriately) so that the temporal discretization is consistent; (d) document the chosen eta and its expected effect on boundary-layer thickness and modeling error.



	Current Stage [C/3]
1. MAJOR

Category: accuracy
Issue: The solver does not enforce incompressibility: reported max|div| is O(1e+0..1e+1) (6–8) throughout the run instead of being driven to near zero after projection.

Recommendation: This is a fundamental accuracy issue. The discrete projection step is inconsistent with the discrete divergence/Poisson operator. Replace the ad-hoc divergence/gradient handling with a consistent pair of discrete operators D (divergence on the staggered grid) and G (face gradient built from cell-centered pressure), and assemble the Poisson matrix as A_p = D * G. Build D and G as sparse matrices and use them consistently for (1) RHS_p = (1/dt) * D * u_star and (2) velocity correction u_new = u_star - dt * G * p. Verify that after one projection the L2 (or max) norm of divergence drops to near numerical zero. Use a simple manufactured / analytical test (e.g. divergence-free prescribed field or lid-driven cavity) to validate D,G,A_p before running the cylinder case.


2. MAJOR

Category: accuracy
Issue: Boundary treatments for the pressure gradient in the projection are inconsistent / ad-hoc (dpdx/dpdy forced to 0 at some face indices). This breaks the discrete compatibility between D and G and prevents mass conservation.

Recommendation: Do not set dpdx/dpdy manually to zero at inlet/outlet faces. Instead, impose physical boundary conditions consistently in the matrix operators (D and G) and in the Poisson problem (Neumann at outlet, Dirichlet/Neumann elsewhere as required). Implement boundary rows/columns in A_p and in G so the projection respects the same boundary treatment used when assembling D. If you want to enforce a pressure reference, only modify one row/column of A_p (and corresponding RHS) as you already do for the nullspace — do not zero gradients manually in the projection loop.


3. MAJOR

Category: accuracy
Issue: Poisson assembly uses local one-sided accumulation on boundaries rather than explicit matching of the discrete divergence operator, risking an inconsistent operator and incorrect Neumann treatment.

Recommendation: Assemble A_p from the discrete D and G operators (A_p = D * G). This automatically ensures the correct boundary treatment and adjoint consistency. If you keep the current direct assembler, validate it by checking that A_p equals D*G (within round-off) and that applying A_p to a pressure field produces the same result as D(G(p)).


4. MINOR

Category: accuracy
Issue: Penalization parameter eta is rather large (1e-2), producing non-negligible slip inside the solid (reported max_dev_in_solid ~ 9e-2).

Recommendation: Reduce eta (e.g. 1e-4 or smaller) to enforce no-slip more tightly. Because penalization is treated implicitly (folded into momentum diagonal), smaller eta is feasible. Monitor conditioning and, if needed, use a more robust preconditioner or increase solver tolerance. Report the norm of velocity in the solid after change; expect it to drop substantially. Consider local grid refinement near the cylinder or increase points_per_diameter to reduce smearing of the immersed boundary.


5. MINOR

Category: structure
Issue: Laplacian matrices L_u / L_v have Dirichlet rows replaced by identity before they are used both in the momentum matrix A and in the RHS viscous term (visc/2 * L_u.dot(u)). This mixes Dirichlet enforcement into the Laplacian used on the RHS and can create inconsistencies.

Recommendation: Keep the pure Laplacian operators (with consistent ghost/Neumann handling) for use in RHS viscous terms, and apply Dirichlet enforcement only on the final A_u/A_v matrices and RHS vectors used in the factorization/solve. Concretely: (1) preserve an untouched L_u/L_v for L.dot(u) in RHS, (2) assemble A_u/A_v from I, L and penalization, and then overwrite rows for Dirichlet nodes as you do. This avoids pollution of the RHS by identity rows and makes the discretization clearer and easier to debug.





### iteration_count
3

### rev_image_description
	Current Stage [A/3]


	Current Stage [B/3]


	Current Stage [C/3]




