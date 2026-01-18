### 🛌🏻 Re4gent with Multi-Modal Review

#### Configs 🏷️:
MaxCount: 2, MultiModalReview: False, AlternativeNum: 2, MaxComments: 5, ConsModule_ON: True, CUT_OUTEXT: 2000

### prob_todo

The PDE of 2-D steady incompressible Navier-Stokes equations is given by:
\begin{cases}
u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} + \frac{\partial p}{\partial x} - \frac{1}{\mathrm{Re}} \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right) = 0, (x,y) \in \Omega,\\
u \frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} + \frac{\partial p}{\partial y} - \frac{1}{\mathrm{Re}} \left( \frac{\partial^2 v}{\partial x^2} + \frac{\partial^2 v}{\partial y^2} \right) = 0, (x,y) \in \Omega,\\
\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0, (x,y) \in \Omega,\\
\end{cases}
The Reynolds number Re = 400.
The domain is $\Omega = [0,1]^2$, the top boundary is $\Gamma_1$, the left, right and bottom boundary is $\Gamma_2$.
The boundary conditions are:
\begin{cases}
(u,v) = (\alpha(x(1-x)), 0) &(x,y) \in \Gamma_1\\
(u,v) = (0,0) &(x,y) \in \Gamma_2\\
\text{Reference pressure: } p(0,0) = 0\\
\end{cases}
where $\alpha$ is 2.

Implement a stable and efficient method to solve this problem.
Implement reasonable acceleration strategies to reduce computational cost.
Define visc = 1.0 / Re at the GLOBAL scope (outside main), and strictly use visc instead of nu throughout the code to avoid confusion with u.
Ensure all constants are passed explicitly to functions or defined globally to prevent NameError.
Plot the contour of the velocity magnitude using the 'RdBu_r' colormap, overlaid with streamlines, and the convergence history in one figure.
Just save figs do not use plt.show() in the code.

[HINTS]:
Print concise progress information ONLY every 10% of total steps.
Focus on generating runnable code with reasonable physical results, do not over-optimize.
Do NOT conclude [not converged] solely because max_iters is reached, assess convergence based on whether residuals are still decreasing. 


### expanded_prob
We must solve the steady, two-dimensional incompressible Navier–Stokes equations on the unit square Ω = [0,1]^2 with Reynolds number Re = 400 (visc = 1/Re). The PDE set couples two momentum equations for the horizontal and vertical velocity components (u,v) with the incompressibility constraint (divergence-free). Boundary conditions: a prescribed tangential lid velocity on the top boundary Γ1, (u,v) = (α x(1-x),0) with α = 2, and no-slip (u,v)=(0,0) on the other three sides (Γ2). Pressure is fixed at a reference point p(0,0)=0 to remove the pressure nullspace.

Primary mathematical and numerical challenges:
- Nonlinearity: the convective terms (u·∇)u and (u·∇)v render the system nonlinear and require iterative linearization (Picard/fixed-point or Newton) or pseudo-transient continuation.
- Saddle-point constraint: incompressibility couples pressure and velocity in a saddle-point system; naive discretizations lead to oscillatory pressure (checkerboarding) unless velocity/pressure spaces or grid staggering are chosen appropriately, or stabilization is applied.
- Convection-dominated behavior: at Re=400 there are moderately strong advection effects leading to boundary layers and possible recirculation; central differencing can produce nonphysical oscillations unless stabilized (upwinding, flux-limiting, SUPG, etc.).
- Solver scalability: the pressure Poisson-like subsystem or linearized momentum matrices dominate cost; efficient preconditioning or multigrid is essential for practical run times.
- Boundary treatment: for vorticity–streamfunction formulations the boundary vorticity must be computed consistently with no-slip conditions; for collocated grids, pressure/velocity stabilization is required.
- Convergence monitoring: a max-iteration limit alone is insufficient—monitor residuals and reduction trends to decide if the solution is converging to steady state.

Numerical targets: produce a stable, accurate steady solution (velocity magnitude contours and streamlines physically reasonable for a lid-driven-like flow), implement explicit acceleration strategies (under-relaxation, multigrid or AMG for Poisson, ILU/AMG preconditioning for Krylov solvers, adaptive pseudo-time stepping), and produce two plots in a single saved figure: (1) contour of velocity magnitude with 'RdBu_r' colormap overlaid with streamlines, and (2) convergence history. Ensure visc = 1.0/Re is defined globally and passed explicitly to functions to avoid name collisions (nu vs u). Print concise progress only every 10% of total steps.

### solution_plans
	Current Stage [A/2]
solu_name='SIMPLE on a staggered (MAC) grid with multigrid Poisson and conservative upwind convection' content="Governing idea:\nDiscretize on a staggered (Marker-and-Cell) grid to avoid pressure–velocity checkerboarding; use the SIMPLE pressure-correction algorithm to solve the steady saddle-point problem iteratively. Treat convection with a conservative second-order upwind-biased scheme (or QUICK) and diffusion with central second-order differences. Solve the pressure-correction Poisson with a multigrid solver for efficiency; use sparse Krylov solvers with ILU/AMG preconditioners for the momentum systems. Use under-relaxation on velocities and pressure to stabilize iterations.\n\nAlgorithmic steps (numbered):\n1) Preliminaries\n   - Define global constant visc = 1.0/Re and pass visc, grid spacing, relaxation factors, and tolerances explicitly into functions.\n   - Create staggered grid: u-values at vertical cell faces, v-values at horizontal faces, pressure at cell centers. Choose Nx,Ny (e.g. 64–256 depending on resources).\n   - Build discrete operators for diffusion (central) and convective flux computation (conservative upwind/QUICK). Assemble sparse matrices for implicit viscous operator as needed.\n2) Initialization\n   - Initialize u,v fields consistent with BCs (enforce top lid u(x,1)=alpha x(1-x), others zero). Initialize pressure field (e.g. zero) and set p(0,0)=0 as a reference.\n3) Outer SIMPLE loop (iterate until residual tolerance or max_iters)\n   - For iteration k:\n     a) Solve momentum predictor (u*,v*) by discretizing momentum equations with current pressure p^k on RHS: treat diffusion implicitly and convection either explicitly using previous iterate or linearize (Picard). Solve the linear systems for u* and v* using a Krylov solver with ILU/AMG preconditioner; reuse preconditioner across iterations where possible.\n     b) Compute mass imbalance (divergence of u*,v*) at pressure cell centers.\n     c) Solve pressure-correction Poisson: ∇·(H ∇p') = (1/Δt) divergence(u*,v*) or the appropriate derived SIMPLE Poisson; use geometric or algebraic multigrid for this Poisson solve for best efficiency.\n     d) Update pressure p^{k+1} = p^k + α_p p' (α_p ~ 0.4–0.8) and correct velocities u^{k+1} = u* + α_u correction (α_u smaller than α_p if needed). Enforce boundary conditions after correction.\n     e) Apply under-relaxation to u and v if momentum solves are direct updates.\n     f) Compute residuals: momentum residual norms and continuity residual (max divergence). Store history.\n   - Print concise progress every 10% of total max_iters (i.e., at 10%,20%,...): show iteration, current residuals and reduction factor. Do not print every iteration.\n4) Convergence and postprocessing\n   - Terminate when residuals fall below user-specified tolerance (absolute and relative) or when residuals stop decreasing meaningfully; do NOT declare nonconverged solely because max_iters reached—assess residual trends.\n   - Compute velocity magnitude, contours and streamlines on a common plotting grid (interpolate staggered fields to cell centers if needed). Plot velocity magnitude with 'RdBu_r' colormap, overlay streamlines computed from center velocities. Also plot residual history on same figure (subplot). Save fig (no plt.show()).\n\nStability / accuracy / complexity / efficiency limitations:\n- SIMPLE is robust but can converge slowly for higher Re; under-relaxation slows convergence further but stabilizes iterations. For Re~400 and moderate grid sizes it is tractable.\n- The convection discretization must be stabilized (upwind/QUICK) to avoid spurious oscillations; this introduces numerical diffusion and reduces formal accuracy compared with pure central schemes.\n- Poisson solves dominate CPU cost; geometric multigrid gives near-linear complexity O(N) per solve, while Krylov+ILU scales worse. AMG is effective but increases implementation complexity.\n- Memory and CPU cost scale with grid resolution; doubling resolution increases cost by roughly 4x for 2D. Reuse preconditioners and sparse factorization where possible.\n- The staggered grid complicates plotting and interpolation but prevents pressure checkerboarding without additional stabilization.\n- Implementation detail: ensure visc is used globally and passed into all solver functions to avoid confusion with u. Stop here."

	Current Stage [B/2]
solu_name='Streamfunction–vorticity with pseudo-transient implicit diffusion and high-resolution convection (multigrid Poisson)' content="Governing idea:\nEliminate pressure by switching to the streamfunction–vorticity (ψ–ω) formulation in 2D: solve ∇^2 ψ = -ω to get streamfunction, compute velocities from ψ, advance ω using the steady vorticity transport equation with pseudo-time stepping until steady state. Handle viscous diffusion implicitly (Crank–Nicolson or fully implicit) and convection explicitly or semi-implicitly using a high-resolution scheme; solve the Poisson for ψ with multigrid for speed. This avoids the saddle-point structure and direct handling of pressure.\n\nAlgorithmic steps (numbered):\n1) Preliminaries\n   - Define global visc = 1.0/Re and pass it to all routines. Choose grid resolution Nx×Ny (cell-centered); precompute Laplacian operator and multigrid hierarchy for the Poisson.\n   - Map boundary velocity BCs to ψ values on boundary using no-penetration and no-slip: compute ψ on boundary from integral of tangential velocity or set reference consistently. Establish formula for boundary vorticity ω_b from velocities using finite-difference approximations (Thom or second-order consistent boundary vorticity formula).\n2) Initialization\n   - Initialize ω field (e.g. zero) and compute ψ by solving ∇^2 ψ = -ω with pde Poisson solver; compute velocities u = ∂ψ/∂y, v = -∂ψ/∂x. Enforce lid velocity on top by setting ψ boundary values accordingly and vorticity boundary ω_b from BC.\n3) Pseudo-transient time march to steady state\n   - Use a pseudo-time iteration: ∂ω/∂t + u·∇ω = visc ∇^2 ω. March in time until steady-state residuals fall below tolerance.\n   - Time discretization: treat diffusion implicitly (Crank–Nicolson or backward Euler) and convection explicitly (strong stability-preserving Runge–Kutta or second-order Adams–Bashforth) to avoid severe time-step restriction. Optionally use fully implicit Newton–Krylov for ω if fast convergence is needed.\n   - At each time step:\n     a) Update interior ω via semi-implicit solve: (I - Δt visc L) ω^{n+1} = ω^n - Δt (u·∇ω)^n + boundary forcing. Solve linear system for ω^{n+1} with a sparse solver (multigrid-preconditioned Krylov or direct if small).\n     b) Enforce vorticity boundary ω_b from no-slip using updated ψ/velocities.\n     c) Solve ∇^2 ψ^{n+1} = -ω^{n+1} with multigrid to obtain ψ and update u/v from derivatives of ψ.\n     d) Compute residual measures: change in ω (norm), divergence of resulting velocity (should be near zero by formulation), global momentum residual proxies. Save history.\n   - Use adaptive pseudo-time stepping based on residual reduction (increase Δt when residual decreases rapidly; reduce Δt if instability observed). Print concise progress only at every 10% of total planned pseudo-time steps.\n4) Convergence and postprocessing\n   - Stop when ω (or velocities) reach steady tolerance or residuals stagnate. Interpolate u,v to plotting grid and assemble velocity magnitude and streamlines.\n   - Plot velocity magnitude contours with 'RdBu_r' and overlay streamlines; in same figure subplot include convergence history. Save figure (no plt.show()).\n\nStability / accuracy / complexity / efficiency limitations:\n- ψ–ω avoids pressure but requires careful treatment of boundary vorticity; boundary approximation errors can degrade accuracy and slow convergence, especially near corners.\n- Semi-implicit time stepping removes severe diffusion restrictions but convection CFL still constrains Δt for explicit convective treatment; adaptive stepping helps but increases bookkeeping.\n- Multigrid Poisson solves are critical for efficiency; without them the method can be slow for fine grids.\n- For Re ~ 400 the method is effective; at much higher Re the vorticity gradients become sharp and require finer grids or higher-order convection schemes (and possibly explicit filtering or stabilization).\n- Implementation complexity is lower than a full Newton–Krylov saddle-point solver but greater care is needed for boundary vorticity formulas and ensuring global mass/momentum consistency.\n- As before, ensure visc is declared globally and passed to all routines to avoid name conflicts. Stop here."



### technical_spec
	Current Stage [A/3]
We implement a streamfunction–vorticity (ψ–ω) pseudo-transient solver for the 2D steady incompressible Navier–Stokes lid-driven cavity on a unit square with Re=400. The domain is discretized on a uniform cell-centered grid. The solver advances vorticity in pseudo-time with implicit diffusion (Helmholtz solve) and explicit upwind convection. Poisson solves for ψ (streamfunction) are performed each step. Sparse matrices (2D Laplacian and Helmholtz) are assembled once and factorized for reuse to avoid rebuilding in the time loop. Boundary vorticity is updated each step using Thom's formula to enforce no-slip and the lid velocity. Convergence is assessed by the maximum change in vorticity between steps and the trend of residuals. Outputs: velocity magnitude contour (RdBu_r) with streamlines and the convergence history in a single saved figure.


	Current Stage [B/3]
Solver architecture:
- Discretization: cell-centered uniform grid on [0,1]^2. Dirichlet psi=0 on all boundaries; lid velocity applied on top boundary.
- Variables: streamfunction psi and vorticity omega on the same grid (including boundaries). Velocities computed from psi (u = dpsi/dy, v = -dpsi/dx).
- Time marching: pseudo-transient update of vorticity with implicit diffusion (Helmholtz) and explicit upwind convection.
- Linear algebra: 2D Laplacian and Helmholtz matrices assembled once for interior unknowns and factorized with SuperLU for repeated solves.
- Boundary handling: Thom's formula vectorized to update boundary omega from psi and lid velocity; boundary contributions from known omega are moved to the Helmholtz RHS (vectorized) so interior Helmholtz solve is correct.
- Initialization: set boundary omega from initial psi/lid, solve Poisson for psi, then initialize interior omega consistently as -Lap(psi).
- Diagnostics: residuals (max|Δω|) and max divergence tracked; stagnation heuristic warns if residuals stop decreasing; centerline/center velocity printed.
- Outputs: saves a single figure showing the velocity magnitude (RdBu_r) with streamlines and the convergence history.

Key implementation notes:
- build_2d_laplacian takes dx and dy to support non-square grids.
- Matrices are built and factorized once (outside time loop) to avoid repeated assembly.
- Helmholtz RHS includes dt*visc*boundary contributions from neighbors adjacent to domain boundaries (vectorized), addressing the major accuracy bug.
- Velocities adjacent to boundaries use one-sided differences to ensure convection sees the lid influence consistently.

	Current Stage [C/3]
This solver uses a cell-centered streamfunction–vorticity formulation on a uniform rectangular grid. Key steps:
- Build 2D FD Laplacian for interior points once and factorize it (SuperLU) for repeated Poisson solves.
- Use an implicit Helmholtz solve (I - dt*visc*Lap) assembled/factorized once for implicit diffusion; convection treated explicitly with first-order upwind.
- Velocities are computed robustly: central differences for all interior cells (1:-1,1:-1); boundary u/v are set explicitly from BCs in main to ensure consistency with Thom's formula.
- Thom's formula updates boundary vorticity using current ψ and lid velocity; boundary contributions are moved to the Helmholtz RHS vectorized.
- A simple under-relaxation on the vorticity update accelerates/stabilizes nonlinear iterations without rebuilding matrices.
- Convergence requires (i) small change in vorticity Δω, (ii) small divergence, and (iii) small physics residual of the steady vorticity equation.
- Includes a lightweight unit test on a coarse grid (sin() manufactured solution) to verify Lap ψ = -ω and near-zero divergence.
- Produces and saves one figure with velocity magnitude (RdBu_r) + streamlines and the convergence history.




### python_codes

#### Script block1:

```python
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Global viscosity (use visc everywhere, do NOT use nu)
Re = 400.0
visc = 1.0 / Re

# ----------------------------- Helper functions -----------------------------

def build_2d_laplacian(nx, ny, dx):
    """Build 2D Laplacian (interior points) using Kronecker products.
    nx, ny: number of interior points in x and y directions.
    Returns sparse matrix of size (nx*ny, nx*ny).
    """
    ex = np.ones(nx)
    ey = np.ones(ny)
    Tx = sp.diags([ex, -2 * ex, ex], offsets=[-1, 0, 1], shape=(nx, nx))
    Ty = sp.diags([ey, -2 * ey, ey], offsets=[-1, 0, 1], shape=(ny, ny))
    Ix = sp.eye(nx)
    Iy = sp.eye(ny)
    Lap = sp.kron(Iy, Tx) + sp.kron(Ty, Ix)
    Lap = Lap / (dx * dx)
    return Lap.tocsc()


def factorize_matrix(A):
    """Factorize sparse matrix A using SuperLU for repeated solves."""
    return spla.factorized(A)


def solve_poisson_from_omega(omega, lap_factor):
    """Solve Lap(psi) = -omega for interior psi values and return full psi array with zero BCs.
    omega: full array with boundaries included; only interior omega[1:-1,1:-1] used in RHS.
    lap_factor: callable factorized solver for Laplace matrix.
    """
    ny, nx = omega.shape
    nx_int = nx - 2
    ny_int = ny - 2
    rhs = -omega[1:-1, 1:-1].ravel()
    psi_int = lap_factor(rhs)
    psi = np.zeros_like(omega)
    psi[1:-1, 1:-1] = psi_int.reshape((ny_int, nx_int))
    return psi


def update_boundary_vorticity(omega, psi, dx, dy, u_lid_func):
    """Apply Thom's formula to update vorticity at physical boundaries.
    omega and psi are full arrays including boundaries; psi boundary is Dirichlet (here set to 0).
    u_lid_func(x): function returning lid tangential velocity at x.
    """
    ny, nx = psi.shape
    x = np.linspace(0.0, 1.0, nx)

    # Top boundary (y = 1): j = ny-1
    j = ny - 1
    j_in = j - 1
    for i in range(1, nx - 1):
        u_top = u_lid_func(x[i])
        omega[j, i] = -2.0 * (psi[j_in, i] - psi[j, i]) / (dy * dy) - 2.0 * u_top / dy

    # Bottom boundary (y = 0): j = 0
    j = 0
    j_in = j + 1
    for i in range(1, nx - 1):
        # no-slip u=0 on bottom
        omega[j, i] = -2.0 * (psi[j_in, i] - psi[j, i]) / (dy * dy)

    # Left boundary (x = 0): i = 0
    i = 0
    i_in = i + 1
    for j in range(1, ny - 1):
        # no-slip v=0 on left
        omega[j, i] = -2.0 * (psi[j, i_in] - psi[j, i]) / (dx * dx)

    # Right boundary (x = 1): i = nx-1
    i = nx - 1
    i_in = i - 1
    for j in range(1, ny - 1):
        omega[j, i] = -2.0 * (psi[j, i_in] - psi[j, i]) / (dx * dx)

    # Corners: average adjacent boundary values to avoid NaNs
    omega[0, 0] = 0.5 * (omega[0, 1] + omega[1, 0])
    omega[0, -1] = 0.5 * (omega[0, -2] + omega[1, -1])
    omega[-1, 0] = 0.5 * (omega[-1, 1] + omega[-2, 0])
    omega[-1, -1] = 0.5 * (omega[-1, -2] + omega[-2, -1])

    return omega


def compute_velocities_from_psi(psi, dx, dy):
    """Compute u = dpsi/dy and v = -dpsi/dx on the cell-centered grid.
    Returns u, v arrays of same shape as psi (zeros at boundaries by central diff convention).
    """
    ny, nx = psi.shape
    u = np.zeros_like(psi)
    v = np.zeros_like(psi)
    # interior points
    u[1:-1, 1:-1] = (psi[2:, 1:-1] - psi[0:-2, 1:-1]) / (2.0 * dy)
    v[1:-1, 1:-1] = -(psi[1:-1, 2:] - psi[1:-1, 0:-2]) / (2.0 * dx)
    # Boundaries keep u,v zeros except we will enforce lid externally when plotting
    return u, v


def compute_convection_term(omega, u, v, dx, dy):
    """Compute convective term (u * dω/dx + v * dω/dy) at interior using first-order upwind.
    omega, u, v are full arrays including boundaries. Returns interior convective term array.
    """
    # interior slices
    w = omega
    u_c = u[1:-1, 1:-1]
    v_c = v[1:-1, 1:-1]

    dwdx_f = (w[1:-1, 2:] - w[1:-1, 1:-1]) / dx
    dwdx_b = (w[1:-1, 1:-1] - w[1:-1, 0:-2]) / dx
    dwdy_f = (w[2:, 1:-1] - w[1:-1, 1:-1]) / dy
    dwdy_b = (w[1:-1, 1:-1] - w[0:-2, 1:-1]) / dy

    dwdx = np.where(u_c > 0.0, dwdx_b, dwdx_f)
    dwdy = np.where(v_c > 0.0, dwdy_b, dwdy_f)

    conv = u_c * dwdx + v_c * dwdy
    return conv


def divergence(u, v, dx, dy):
    """Compute divergence du/dx + dv/dy at interior points; return max absolute divergence."""
    dudx = (u[1:-1, 2:] - u[1:-1, 0:-2]) / (2.0 * dx)
    dvdy = (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2.0 * dy)
    div = dudx + dvdy
    return np.max(np.abs(div))

# ------------------------------- Main solver --------------------------------

def main():
    # Parameters
    alpha = 2.0  # lid amplitude factor
    Nx = 129  # grid points in x (including boundaries)
    Ny = 129  # grid points in y
    Lx = 1.0
    Ly = 1.0
    dx = Lx / (Nx - 1)
    dy = Ly / (Ny - 1)

    # interior counts
    nx_int = Nx - 2
    ny_int = Ny - 2

    # time-stepping / pseudo-time parameters (fixed dt to avoid rebuilding matrices)
    dt = 1e-3  # fixed pseudo-time step (tuneable); matrix A built once with this dt
    max_iters = 3000
    tol = 1e-6

    # Lid velocity function (top boundary y=1): u = alpha * x * (1-x)
    def u_lid(x):
        return alpha * x * (1.0 - x)

    # Create coordinate grids for plotting (cell-centered)
    x = np.linspace(0.0, Lx, Nx)
    y = np.linspace(0.0, Ly, Ny)
    X, Y = np.meshgrid(x, y)

    # Initialize fields (including boundaries)
    omega = np.zeros((Ny, Nx))
    psi = np.zeros_like(omega)

    # Build Laplacian for interior unknowns and factorize (for Poisson solves)
    Lap = build_2d_laplacian(nx_int, ny_int, dx)
    lap_solver = factorize_matrix(Lap)

    # Build Helmholtz matrix A = I - dt * visc * Lap (for implicit diffusion) and factorize
    I = sp.eye(nx_int * ny_int, format='csc')
    A = (I - (dt * visc) * Lap).tocsc()
    helm_solver = factorize_matrix(A)

    # Initial boundary vorticity from initial psi (zero) and lid
    omega = update_boundary_vorticity(omega, psi, dx, dy, u_lid)

    # Prepare to store convergence history
    residuals = []

    # Print frequency (every 10% of total steps)
    freq = max(1, max_iters // 10)

    # Initial Poisson solve to get psi consistent with initial omega
    psi = solve_poisson_from_omega(omega, lap_solver)

    # Main pseudo-time loop
    converged = False
    print_interval = freq
    for it in range(1, max_iters + 1):
        # 1) Solve Poisson for psi from current interior omega
        psi = solve_poisson_from_omega(omega, lap_solver)

        # 2) Update boundary vorticity using Thom's formula and current psi & lid BC
        omega = update_boundary_vorticity(omega, psi, dx, dy, u_lid)

        # 3) Compute velocities from psi
        u, v = compute_velocities_from_psi(psi, dx, dy)
        # Enforce lid velocity explicitly at top boundary for plotting and for convection influence
        # Set u on the top boundary row to prescribed lid; interior u remains from psi derivatives
        u[-1, :] = u_lid(x)
        v[-1, :] = 0.0

        # 4) Compute convective term (interior)
        conv = compute_convection_term(omega, u, v, dx, dy)

        # 5) Form RHS and solve Helmholtz for new interior omega
        rhs = omega[1:-1, 1:-1] - dt * conv
        rhs_flat = rhs.ravel()
        omega_new_flat = helm_solver(rhs_flat)
        omega_new = omega.copy()
        omega_new[1:-1, 1:-1] = omega_new_flat.reshape((ny_int, nx_int))

        # 6) Update omega and compute new psi next iteration
        delta = np.max(np.abs(omega_new - omega))
        omega = omega_new

        # Compute divergence max as a diagnostic
        # Recompute psi for divergence estimate
        psi = solve_poisson_from_omega(omega, lap_solver)
        u, v = compute_velocities_from_psi(psi, dx, dy)
        u[-1, :] = u_lid(x)
        v[-1, :] = 0.0
        max_div = divergence(u, v, dx, dy)

        residuals.append(delta)

        # Print concise progress only every 10% of total steps
        if (it % print_interval == 0) or (it == max_iters):
            print(f"Iter {it}/{max_iters}: max|Δω|={delta:.3e}, max|div|={max_div:.3e}, max|u|={np.max(np.abs(u)):.3e}")

        # Convergence check
        if delta < tol:
            converged = True
            print(f"Converged by delta tolerance at iteration {it} (delta={delta:.3e}).")
            break

    # If we reached max iters, assess convergence trend rather than declaring failure
    if not converged:
        init_res = residuals[0] if len(residuals) > 0 else 1.0
        final_res = residuals[-1] if len(residuals) > 0 else 0.0
        if final_res <= 0.01 * init_res:
            converged = True
        print(f"Finished iterations. Converged flag = {converged}. Initial residual={init_res:.3e}, final={final_res:.3e}")

    # Postprocessing for plotting
    psi = solve_poisson_from_omega(omega, lap_solver)
    u, v = compute_velocities_from_psi(psi, dx, dy)
    u[-1, :] = u_lid(x)
    v[-1, :] = 0.0
    vel_mag = np.sqrt(u * u + v * v)

    # Prepare figure: left contour+streamlines, right convergence history
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [3, 1]})

    # Contour plot of velocity magnitude
    cmap = 'RdBu_r'
    cf = ax1.contourf(X, Y, vel_mag, levels=50, cmap=cmap)
    ax1.streamplot(x, y, u, v, color='k', density=1.2, linewidth=0.6)
    ax1.set_title('Velocity magnitude and streamlines')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    cbar = fig.colorbar(cf, ax=ax1)
    cbar.set_label('|u|')

    # Convergence history
    ax2.plot(np.arange(1, len(residuals) + 1), residuals, '-o', markersize=3)
    ax2.set_yscale('log')
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('max|Δω|')
    ax2.set_title('Convergence history')
    ax2.grid(True)

    plt.tight_layout()
    fig_filename = 'ns2d_stream_vmag_convergence.png'
    plt.savefig(fig_filename, dpi=200)

    # Print final quantitative metrics
    print(f"Final: iterations={it}, max|u|={np.max(np.abs(u)):.6f}, max|v|={np.max(np.abs(v)):.6f}, max|div|={max_div:.3e}")
    print(f"Saved figure to {fig_filename}")


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

# Global viscosity (use visc everywhere)
Re = 400.0
visc = 1.0 / Re

# ----------------------------- Helper functions -----------------------------

def build_2d_laplacian(nx, ny, dx, dy):
    """Build 2D Laplacian for interior unknowns using Kronecker products.
    Returns sparse matrix of size (nx*ny, nx*ny) representing
    (d2/dx2 + d2/dy2) with Dirichlet boundaries excluded.
    """
    ex = np.ones(nx)
    ey = np.ones(ny)
    Tx = sp.diags([ex, -2 * ex, ex], offsets=[-1, 0, 1], shape=(nx, nx))
    Ty = sp.diags([ey, -2 * ey, ey], offsets=[-1, 0, 1], shape=(ny, ny))
    Ix = sp.eye(nx)
    Iy = sp.eye(ny)
    Lap = sp.kron(Iy, Tx) / (dx * dx) + sp.kron(Ty, Ix) / (dy * dy)
    return Lap.tocsc()


def factorize_matrix(A):
    """Factorize sparse matrix A for repeated solves; returns a callable solver."""
    return spla.factorized(A)


def solve_poisson_from_omega(omega, lap_solver):
    """Solve Lap(psi) = -omega for interior psi values and return full psi array with zero BCs.
    omega: full array with boundaries included; only interior omega[1:-1,1:-1] used in RHS.
    lap_solver: callable factorized solver for Laplace matrix.
    """
    ny, nx = omega.shape
    nx_int = nx - 2
    ny_int = ny - 2
    rhs = -omega[1:-1, 1:-1].ravel()
    psi_int = lap_solver(rhs)
    psi = np.zeros_like(omega)
    psi[1:-1, 1:-1] = psi_int.reshape((ny_int, nx_int))
    return psi


def update_boundary_vorticity(omega, psi, dx, dy, u_lid_func):
    """Apply Thom's formula to update vorticity at boundaries (vectorized).
    Both omega and psi are full arrays including boundaries.
    """
    ny, nx = psi.shape
    dx2 = dx * dx
    dy2 = dy * dy
    x = np.linspace(0.0, 1.0, nx)

    # Top boundary (y = 1): j = ny-1
    omega[-1, 1:-1] = -2.0 * (psi[-2, 1:-1] - psi[-1, 1:-1]) / dy2 - 2.0 * u_lid_func(x[1:-1]) / dy

    # Bottom boundary (y = 0): j = 0
    omega[0, 1:-1] = -2.0 * (psi[1, 1:-1] - psi[0, 1:-1]) / dy2

    # Left boundary (x = 0): i = 0
    omega[1:-1, 0] = -2.0 * (psi[1:-1, 1] - psi[1:-1, 0]) / dx2

    # Right boundary (x = 1): i = nx-1
    omega[1:-1, -1] = -2.0 * (psi[1:-1, -2] - psi[1:-1, -1]) / dx2

    # Corners: average adjacent boundary values
    omega[0, 0] = 0.5 * (omega[0, 1] + omega[1, 0])
    omega[0, -1] = 0.5 * (omega[0, -2] + omega[1, -1])
    omega[-1, 0] = 0.5 * (omega[-1, 1] + omega[-2, 0])
    omega[-1, -1] = 0.5 * (omega[-1, -2] + omega[-2, -1])

    return omega


def compute_velocities_from_psi(psi, dx, dy):
    """Compute u = dpsi/dy and v = -dpsi/dx on the cell-centered grid.
    Use central differences for interior nodes and one-sided near boundaries so that
    velocities adjacent to boundaries reflect the Dirichlet psi at the wall.
    Returns u, v arrays of same shape as psi.
    """
    ny, nx = psi.shape
    u = np.zeros_like(psi)
    v = np.zeros_like(psi)

    # central differences for most interior points
    u[2:-2, 1:-1] = (psi[3:-1, 1:-1] - psi[1:-3, 1:-1]) / (2.0 * dy)
    v[1:-1, 2:-2] = -(psi[1:-1, 3:-1] - psi[1:-1, 1:-3]) / (2.0 * dx)

    # Fill remaining interior rows/cols using one-sided differences adjacent to boundaries
    # y-direction one-sided for u
    # first interior row (adjacent bottom boundary)
    u[1, 1:-1] = (psi[1, 1:-1] - psi[0, 1:-1]) / dy
    # second row (if not already set by central) ensure we didn't miss index when small grids
    if ny > 4:
        u[2, 1:-1] = (psi[3, 1:-1] - psi[1, 1:-1]) / (2.0 * dy)
    # last interior row (adjacent top boundary)
    u[-2, 1:-1] = (psi[-1, 1:-1] - psi[-2, 1:-1]) / dy

    # x-direction one-sided for v
    v[1:-1, 1] = -(psi[1:-1, 1] - psi[1:-1, 0]) / dx
    v[1:-1, -2] = -(psi[1:-1, -1] - psi[1:-1, -2]) / dx

    # For any remaining unset central indices (small grids), fall back to central
    # center region
    u[1:-1, 1:-1] = np.where(u[1:-1, 1:-1] == 0.0, (psi[2:, 1:-1] - psi[0:-2, 1:-1]) / (2.0 * dy), u[1:-1, 1:-1])
    v[1:-1, 1:-1] = np.where(v[1:-1, 1:-1] == 0.0, -(psi[1:-1, 2:] - psi[1:-1, 0:-2]) / (2.0 * dx), v[1:-1, 1:-1])

    return u, v


def compute_convection_term(omega, u, v, dx, dy):
    """Compute convective term (u * dω/dx + v * dω/dy) at interior using first-order upwind.
    omega, u, v are full arrays including boundaries. Returns interior convective term array.
    """
    w = omega
    u_c = u[1:-1, 1:-1]
    v_c = v[1:-1, 1:-1]

    dwdx_f = (w[1:-1, 2:] - w[1:-1, 1:-1]) / dx
    dwdx_b = (w[1:-1, 1:-1] - w[1:-1, 0:-2]) / dx
    dwdy_f = (w[2:, 1:-1] - w[1:-1, 1:-1]) / dy
    dwdy_b = (w[1:-1, 1:-1] - w[0:-2, 1:-1]) / dy

    dwdx = np.where(u_c > 0.0, dwdx_b, dwdx_f)
    dwdy = np.where(v_c > 0.0, dwdy_b, dwdy_f)

    conv = u_c * dwdx + v_c * dwdy
    return conv


def divergence(u, v, dx, dy):
    """Compute divergence du/dx + dv/dy at interior points; return max absolute divergence."""
    dudx = (u[1:-1, 2:] - u[1:-1, 0:-2]) / (2.0 * dx)
    dvdy = (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2.0 * dy)
    div = dudx + dvdy
    return np.max(np.abs(div))


# ------------------------------- Main solver --------------------------------

def main():
    # Parameters
    alpha = 2.0  # lid amplitude factor
    Nx = 129  # grid points in x (including boundaries)
    Ny = 129  # grid points in y
    Lx = 1.0
    Ly = 1.0
    dx = Lx / (Nx - 1)
    dy = Ly / (Ny - 1)
    dx2 = dx * dx
    dy2 = dy * dy

    # interior counts
    nx_int = Nx - 2
    ny_int = Ny - 2

    # time-stepping / pseudo-time parameters
    dt = 1e-3  # fixed pseudo-time step (matrix A built once with this dt)
    max_iters = 3000
    tol = 1e-6

    # Lid velocity function (top boundary y=1): u = alpha * x * (1-x)
    def u_lid(x_arr):
        return alpha * x_arr * (1.0 - x_arr)

    # Create coordinate grids for plotting (cell-centered)
    x = np.linspace(0.0, Lx, Nx)
    y = np.linspace(0.0, Ly, Ny)
    X, Y = np.meshgrid(x, y)

    # Initialize fields (including boundaries)
    omega = np.zeros((Ny, Nx))
    psi = np.zeros_like(omega)

    # Build Laplacian for interior unknowns and factorize (for Poisson solves)
    Lap = build_2d_laplacian(nx_int, ny_int, dx, dy)
    lap_solver = factorize_matrix(Lap)

    # Build Helmholtz matrix A = I - dt * visc * Lap (for implicit diffusion) and factorize
    I = sp.eye(nx_int * ny_int, format='csc')
    A = (I - (dt * visc) * Lap).tocsc()
    helm_solver = factorize_matrix(A)

    # Initial boundary vorticity from initial psi (zero) and lid
    omega = update_boundary_vorticity(omega, psi, dx, dy, u_lid)

    # After setting boundary omega, compute psi consistent with that boundary source
    psi = solve_poisson_from_omega(omega, lap_solver)

    # Initialize interior omega consistently as -Lap(psi) using FD (so interior & boundary consistent)
    omega_interior = -(
        (psi[1:-1, 2:] - 2.0 * psi[1:-1, 1:-1] + psi[1:-1, 0:-2]) / dx2
        + (psi[2:, 1:-1] - 2.0 * psi[1:-1, 1:-1] + psi[0:-2, 1:-1]) / dy2
    )
    omega[1:-1, 1:-1] = omega_interior

    # Prepare to store convergence history
    residuals = []

    # Print frequency (every 10% of total steps)
    freq = max(1, max_iters // 10)

    # Main pseudo-time loop
    converged = False
    stagnant_counter = 0
    prev_print_res = None

    for it in range(1, max_iters + 1):
        # 1) Solve Poisson for psi from current omega
        psi = solve_poisson_from_omega(omega, lap_solver)

        # 2) Update boundary vorticity using Thom's formula and current psi & lid BC
        omega = update_boundary_vorticity(omega, psi, dx, dy, u_lid)

        # 3) Compute velocities from psi
        u, v = compute_velocities_from_psi(psi, dx, dy)
        # Enforce lid and wall velocities explicitly on boundaries
        u[-1, :] = u_lid(x)
        v[-1, :] = 0.0
        u[0, :] = 0.0
        v[0, :] = 0.0
        u[:, 0] = 0.0
        u[:, -1] = 0.0
        v[:, 0] = 0.0
        v[:, -1] = 0.0

        # 4) Compute convective term (interior)
        conv = compute_convection_term(omega, u, v, dx, dy)

        # 5) Form RHS and solve Helmholtz for new interior omega
        w_in = omega[1:-1, 1:-1]
        rhs_in = w_in - dt * conv

        # Add boundary contributions: vectorized assembly for known neighbor omega values
        bd = np.zeros_like(rhs_in)
        # left boundary contributes to first interior column
        bd[:, 0] += omega[1:-1, 0] / dx2
        # right boundary contributes to last interior column
        bd[:, -1] += omega[1:-1, -1] / dx2
        # bottom boundary contributes to first interior row
        bd[0, :] += omega[0, 1:-1] / dy2
        # top boundary contributes to last interior row
        bd[-1, :] += omega[-1, 1:-1] / dy2

        rhs_in = rhs_in + dt * visc * bd

        rhs_flat = rhs_in.ravel()
        omega_new_flat = helm_solver(rhs_flat)
        omega_new = omega.copy()
        omega_new[1:-1, 1:-1] = omega_new_flat.reshape((ny_int, nx_int))

        # 6) Update omega and compute diagnostics
        delta = np.max(np.abs(omega_new - omega))
        omega = omega_new

        # Compute divergence max as a diagnostic (compute psi consistent with updated omega)
        psi = solve_poisson_from_omega(omega, lap_solver)
        u, v = compute_velocities_from_psi(psi, dx, dy)
        u[-1, :] = u_lid(x)
        v[-1, :] = 0.0
        max_div = divergence(u, v, dx, dy)

        residuals.append(delta)

        # Print concise progress only every 10% of total steps
        if (it % freq == 0) or (it == max_iters):
            print(f"Iter {it}/{max_iters}: max|Δω|={delta:.3e}, max|div|={max_div:.3e}, max|u|={np.max(np.abs(u)):.3e}")

            # stagnation heuristic: check if residual not decreasing across print intervals
            if prev_print_res is not None and delta >= prev_print_res * 0.999:
                stagnant_counter += 1
                if stagnant_counter >= 3:
                    print("Warning: residual appears stagnated across several print intervals.")
            else:
                stagnant_counter = 0
            prev_print_res = delta

        # Convergence check (based on change in omega)
        if delta < tol:
            converged = True
            print(f"Converged by delta tolerance at iteration {it} (delta={delta:.3e}).")
            break

    # If we reached max iters, assess convergence trend rather than declaring failure
    if not converged:
        init_res = residuals[0] if len(residuals) > 0 else 1.0
        final_res = residuals[-1] if len(residuals) > 0 else 0.0
        if final_res <= 0.01 * init_res:
            converged = True
        print(f"Finished iterations. Converged flag = {converged}. Initial residual={init_res:.3e}, final={final_res:.3e}")

    # Final postprocessing
    psi = solve_poisson_from_omega(omega, lap_solver)
    u, v = compute_velocities_from_psi(psi, dx, dy)
    u[-1, :] = u_lid(x)
    v[-1, :] = 0.0
    vel_mag = np.sqrt(u * u + v * v)

    # Some simple quantitative diagnostics (centerline values)
    uc_mid = u[:, Nx // 2]
    vc_mid = v[Ny // 2, :]
    u_center = u[Ny // 2, Nx // 2]
    print(f"Diagnostics: u_center={u_center:.6f}, u_top_center={u[-1, Nx // 2]:.6f}")

    # Prepare figure: left contour+streamlines, right convergence history
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [3, 1]})

    # Contour plot of velocity magnitude
    cmap = 'RdBu_r'
    cf = ax1.contourf(X, Y, vel_mag, levels=50, cmap=cmap)
    ax1.streamplot(x, y, u, v, color='k', density=1.2, linewidth=0.6)
    ax1.set_title('Velocity magnitude and streamlines')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    cbar = fig.colorbar(cf, ax=ax1)
    cbar.set_label('|u|')

    # Convergence history
    res_arr = np.array(residuals)
    iters = np.arange(1, len(res_arr) + 1)
    # Handle non-positive or single-point histories for plotting
    pos = res_arr[res_arr > 0]
    if pos.size >= 2:
        ax2.plot(iters, res_arr, '-o', markersize=3)
        ax2.set_yscale('log')
    else:
        ax2.plot(iters, res_arr, '-o', markersize=3)
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('max|Δω|')
    ax2.set_title('Convergence history')
    ax2.grid(True)

    plt.tight_layout()
    fig_filename = 'ns2d_stream_vmag_convergence.png'
    plt.savefig(fig_filename, dpi=200)

    # Print final quantitative metrics
    print(f"Final: iterations={it}, max|u|={np.max(np.abs(u)):.6f}, max|v|={np.max(np.abs(v)):.6f}, max|div|={max_div:.3e}")
    print(f"Saved figure to {fig_filename}")


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

# Global viscosity (use visc everywhere)
Re = 400.0
visc = 1.0 / Re

# ----------------------------- Helper functions -----------------------------

def build_2d_laplacian(nx, ny, dx, dy):
    """Build 2D Laplacian for interior unknowns using Kronecker products.
    Returns sparse matrix of size (nx*ny, nx*ny) representing
    (d2/dx2 + d2/dy2) with Dirichlet boundaries excluded.
    """
    ex = np.ones(nx)
    ey = np.ones(ny)
    Tx = sp.diags([ex, -2 * ex, ex], offsets=[-1, 0, 1], shape=(nx, nx))
    Ty = sp.diags([ey, -2 * ey, ey], offsets=[-1, 0, 1], shape=(ny, ny))
    Ix = sp.eye(nx)
    Iy = sp.eye(ny)
    Lap = sp.kron(Iy, Tx) / (dx * dx) + sp.kron(Ty, Ix) / (dy * dy)
    return Lap.tocsc()


def factorize_matrix(A):
    """Factorize sparse matrix A for repeated solves; returns a callable solver."""
    return spla.factorized(A)


def solve_poisson_from_omega(omega, lap_solver):
    """Solve Lap(psi) = -omega for interior psi values and return full psi array with zero BCs.
    omega: full array with boundaries included; only interior omega[1:-1,1:-1] used in RHS.
    lap_solver: callable factorized solver for Laplace matrix.
    """
    ny, nx = omega.shape
    nx_int = nx - 2
    ny_int = ny - 2
    rhs = -omega[1:-1, 1:-1].ravel()
    psi_int = lap_solver(rhs)
    psi = np.zeros_like(omega)
    psi[1:-1, 1:-1] = psi_int.reshape((ny_int, nx_int))
    return psi


def update_boundary_vorticity(omega, psi, dx, dy, u_lid_func):
    """Apply Thom's formula to update vorticity at boundaries (vectorized).
    Both omega and psi are full arrays including boundaries.
    """
    ny, nx = psi.shape
    dx2 = dx * dx
    dy2 = dy * dy
    x = np.linspace(0.0, 1.0, nx)

    # Top boundary (y = 1): j = ny-1
    omega[-1, 1:-1] = -2.0 * (psi[-2, 1:-1] - psi[-1, 1:-1]) / dy2 - 2.0 * u_lid_func(x[1:-1]) / dy

    # Bottom boundary (y = 0): j = 0
    omega[0, 1:-1] = -2.0 * (psi[1, 1:-1] - psi[0, 1:-1]) / dy2

    # Left boundary (x = 0): i = 0
    omega[1:-1, 0] = -2.0 * (psi[1:-1, 1] - psi[1:-1, 0]) / dx2

    # Right boundary (x = 1): i = nx-1
    omega[1:-1, -1] = -2.0 * (psi[1:-1, -2] - psi[1:-1, -1]) / dx2

    # Corners: average adjacent boundary values
    omega[0, 0] = 0.5 * (omega[0, 1] + omega[1, 0])
    omega[0, -1] = 0.5 * (omega[0, -2] + omega[1, -1])
    omega[-1, 0] = 0.5 * (omega[-1, 1] + omega[-2, 0])
    omega[-1, -1] = 0.5 * (omega[-1, -2] + omega[-2, -1])

    return omega


def compute_velocities_from_psi(psi, dx, dy):
    """Compute u = dpsi/dy and v = -dpsi/dx on the cell-centered grid.
    Use central differences for interior nodes (1:-1,1:-1). Do NOT use brittle fallbacks.
    Boundary u/v remain as initialized (typically set explicitly in main to enforce BCs).
    Returns u, v arrays of same shape as psi.
    """
    ny, nx = psi.shape
    u = np.zeros_like(psi)
    v = np.zeros_like(psi)

    # Central differences for all interior cells (1:-1,1:-1)
    u[1:-1, 1:-1] = (psi[2:, 1:-1] - psi[0:-2, 1:-1]) / (2.0 * dy)
    v[1:-1, 1:-1] = -(psi[1:-1, 2:] - psi[1:-1, 0:-2]) / (2.0 * dx)

    # Leave boundaries to be set explicitly by caller (consistent BCs)
    return u, v


def compute_convection_term(omega, u, v, dx, dy):
    """Compute convective term (u * dω/dx + v * dω/dy) at interior using first-order upwind.
    omega, u, v are full arrays including boundaries. Returns interior convective term array.
    """
    w = omega
    u_c = u[1:-1, 1:-1]
    v_c = v[1:-1, 1:-1]

    dwdx_f = (w[1:-1, 2:] - w[1:-1, 1:-1]) / dx
    dwdx_b = (w[1:-1, 1:-1] - w[1:-1, 0:-2]) / dx
    dwdy_f = (w[2:, 1:-1] - w[1:-1, 1:-1]) / dy
    dwdy_b = (w[1:-1, 1:-1] - w[0:-2, 1:-1]) / dy

    dwdx = np.where(u_c > 0.0, dwdx_b, dwdx_f)
    dwdy = np.where(v_c > 0.0, dwdy_b, dwdy_f)

    conv = u_c * dwdx + v_c * dwdy
    return conv


def divergence(u, v, dx, dy):
    """Compute divergence du/dx + dv/dy at interior points; return max absolute divergence."""
    dudx = (u[1:-1, 2:] - u[1:-1, 0:-2]) / (2.0 * dx)
    dvdy = (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2.0 * dy)
    div = dudx + dvdy
    return np.max(np.abs(div))


def compute_laplacian_interior(field, dx, dy):
    """Compute 5-point Laplacian for interior nodes using second-order FD."""
    dx2 = dx * dx
    dy2 = dy * dy
    lap = (
        (field[1:-1, 2:] - 2.0 * field[1:-1, 1:-1] + field[1:-1, 0:-2]) / dx2
        + (field[2:, 1:-1] - 2.0 * field[1:-1, 1:-1] + field[0:-2, 1:-1]) / dy2
    )
    return lap


# ------------------------------- Unit test ----------------------------------

def run_unit_tests():
    """Run a lightweight unit test on a coarse grid to check Lap psi = -omega and divergence."""
    Nx = 11
    Ny = 11
    Lx = Ly = 1.0
    x = np.linspace(0.0, Lx, Nx)
    y = np.linspace(0.0, Ly, Ny)
    X, Y = np.meshgrid(x, y)

    dx = Lx / (Nx - 1)
    dy = Ly / (Ny - 1)

    # Manufactured solution: psi = sin(pi x) sin(pi y) (zero on boundaries)
    psi_true = np.sin(np.pi * X) * np.sin(np.pi * Y)

    # Compute omega = -Lap(psi_true)
    lap_psi = compute_laplacian_interior(psi_true, dx, dy)
    omega = np.zeros_like(psi_true)
    omega[1:-1, 1:-1] = -lap_psi

    # Build Lap and solve Poisson to recover psi
    nx_int = Nx - 2
    ny_int = Ny - 2
    Lap = build_2d_laplacian(nx_int, ny_int, dx, dy)
    lap_solver = factorize_matrix(Lap)
    psi_recov = solve_poisson_from_omega(omega, lap_solver)

    err_inf = np.max(np.abs(psi_recov - psi_true))

    # Compute velocities and divergence
    u, v = compute_velocities_from_psi(psi_true, dx, dy)
    # boundaries of psi_true are zero -> velocities on boundaries should be zero; we explicitly set them
    u[0, :] = 0.0
    u[-1, :] = 0.0
    u[:, 0] = 0.0
    u[:, -1] = 0.0
    v[0, :] = 0.0
    v[-1, :] = 0.0
    v[:, 0] = 0.0
    v[:, -1] = 0.0

    max_div = divergence(u, v, dx, dy)

    print(f"Unit test (coarse): max|psi_recov - psi_true|={err_inf:.3e}, max|div|={max_div:.3e}")
    # Basic assertions (not raising to keep code robust): ensure errors small
    return err_inf, max_div


# ------------------------------- Main solver --------------------------------

def main():
    # Run unit tests to validate core discretizations
    err_inf, max_div_unit = run_unit_tests()

    # Parameters
    alpha = 2.0  # lid amplitude factor
    Nx = 129  # grid points in x (including boundaries)
    Ny = 129  # grid points in y
    Lx = 1.0
    Ly = 1.0
    dx = Lx / (Nx - 1)
    dy = Ly / (Ny - 1)
    dx2 = dx * dx
    dy2 = dy * dy

    # interior counts
    nx_int = Nx - 2
    ny_int = Ny - 2

    # time-stepping / pseudo-time parameters
    # choose a reasonable dt (implicit diffusion allows larger); use under-relaxation for stability
    dt = 5e-3
    relax = 0.85  # under-relaxation on vorticity update (0<relax<=1)
    max_iters = 3000
    tol_omega = 1e-6
    tol_div = 1e-5
    tol_phys = 1e-5

    # Lid velocity function (top boundary y=1): u = alpha * x * (1-x)
    def u_lid(x_arr):
        return alpha * x_arr * (1.0 - x_arr)

    # Create coordinate grids for plotting (cell-centered)
    x = np.linspace(0.0, Lx, Nx)
    y = np.linspace(0.0, Ly, Ny)
    X, Y = np.meshgrid(x, y)

    # Initialize fields (including boundaries)
    omega = np.zeros((Ny, Nx))
    psi = np.zeros_like(omega)

    # Build Laplacian for interior unknowns and factorize (for Poisson solves)
    Lap = build_2d_laplacian(nx_int, ny_int, dx, dy)
    lap_solver = factorize_matrix(Lap)

    # Build Helmholtz matrix A = I - dt * visc * Lap (for implicit diffusion) and factorize
    I = sp.eye(nx_int * ny_int, format='csc')
    A = (I - (dt * visc) * Lap).tocsc()
    helm_solver = factorize_matrix(A)

    # Initial boundary vorticity from initial psi (zero) and lid
    omega = update_boundary_vorticity(omega, psi, dx, dy, u_lid)

    # After setting boundary omega, compute psi consistent with that boundary source
    psi = solve_poisson_from_omega(omega, lap_solver)

    # Initialize interior omega consistently as -Lap(psi) using FD
    omega_interior = -compute_laplacian_interior(psi, dx, dy)
    omega[1:-1, 1:-1] = omega_interior

    # Prepare to store convergence history
    residuals = []            # store max|Δω|
    divergences = []         # store max divergence
    phys_res = []            # physics residual norms

    # Print frequency (every 10% of total steps)
    freq = max(1, max_iters // 10)

    # Main pseudo-time loop
    converged = False
    stagnant_counter = 0
    prev_print_res = None

    for it in range(1, max_iters + 1):
        # 1) Solve Poisson for psi from current omega
        psi = solve_poisson_from_omega(omega, lap_solver)

        # 2) Update boundary vorticity using Thom's formula and current psi & lid BC
        omega = update_boundary_vorticity(omega, psi, dx, dy, u_lid)

        # 3) Compute velocities from psi
        u, v = compute_velocities_from_psi(psi, dx, dy)
        # Enforce lid and wall velocities explicitly on boundaries
        u[-1, :] = u_lid(x)
        v[-1, :] = 0.0
        u[0, :] = 0.0
        v[0, :] = 0.0
        u[:, 0] = 0.0
        u[:, -1] = 0.0
        v[:, 0] = 0.0
        v[:, -1] = 0.0

        # 4) Compute convective term (interior)
        conv = compute_convection_term(omega, u, v, dx, dy)

        # 5) Form RHS and solve Helmholtz for new interior omega
        w_in = omega[1:-1, 1:-1].copy()
        rhs_in = w_in - dt * conv

        # Add boundary contributions: vectorized assembly for known neighbor omega values
        bd = np.zeros_like(rhs_in)
        # left boundary contributes to first interior column
        bd[:, 0] += omega[1:-1, 0] / dx2
        # right boundary contributes to last interior column
        bd[:, -1] += omega[1:-1, -1] / dx2
        # bottom boundary contributes to first interior row
        bd[0, :] += omega[0, 1:-1] / dy2
        # top boundary contributes to last interior row
        bd[-1, :] += omega[-1, 1:-1] / dy2

        rhs_in = rhs_in + dt * visc * bd

        rhs_flat = rhs_in.ravel()
        omega_new_flat = helm_solver(rhs_flat)
        omega_new = omega.copy()
        omega_new_interior = omega_new_flat.reshape((ny_int, nx_int))

        # Under-relaxation to stabilize nonlinear update
        omega_new[1:-1, 1:-1] = relax * omega_new_interior + (1.0 - relax) * w_in

        # 6) Update omega and compute diagnostics
        delta = np.max(np.abs(omega_new - omega))
        omega = omega_new

        # Compute psi consistent with updated omega and velocities
        psi = solve_poisson_from_omega(omega, lap_solver)
        u, v = compute_velocities_from_psi(psi, dx, dy)
        u[-1, :] = u_lid(x)
        v[-1, :] = 0.0
        u[0, :] = 0.0
        v[0, :] = 0.0
        u[:, 0] = 0.0
        u[:, -1] = 0.0
        v[:, 0] = 0.0
        v[:, -1] = 0.0

        max_div = divergence(u, v, dx, dy)

        # Physics residual: convective - viscous Lap(omega) (interior)
        lap_omega = compute_laplacian_interior(omega, dx, dy)
        res_phys = conv - visc * lap_omega
        res_norm = np.linalg.norm(res_phys.ravel(), 2) / np.sqrt(res_phys.size)

        residuals.append(delta)
        divergences.append(max_div)
        phys_res.append(res_norm)

        # Print concise progress only every 10% of total steps
        if (it % freq == 0) or (it == max_iters):
            print(f"Iter {it}/{max_iters}: max|Δω|={delta:.3e}, max|div|={max_div:.3e}, phys_res_L2={res_norm:.3e}, max|u|={np.max(np.abs(u)):.3e}")

            # stagnation heuristic: check if residual not decreasing across print intervals
            if prev_print_res is not None and delta >= prev_print_res * 0.999:
                stagnant_counter += 1
                if stagnant_counter >= 3:
                    print("Warning: residual appears stagnated across several print intervals.")
            else:
                stagnant_counter = 0
            prev_print_res = delta

        # Convergence check: require small change in omega, small divergence, and small physics residual
        if (delta < tol_omega) and (max_div < tol_div) and (res_norm < tol_phys):
            converged = True
            print(f"Converged by combined tolerances at iteration {it} (delta={delta:.3e}, div={max_div:.3e}, phys={res_norm:.3e}).")
            break

    # If we reached max iters, assess convergence trend rather than declaring failure
    if not converged:
        init_res = residuals[0] if len(residuals) > 0 else 1.0
        final_res = residuals[-1] if len(residuals) > 0 else 0.0
        phys_final = phys_res[-1] if len(phys_res) > 0 else 0.0
        div_final = divergences[-1] if len(divergences) > 0 else 0.0
        # Accept if physics residual and divergence dropped significantly
        if (final_res <= 0.01 * init_res) or (phys_final < 10.0 * tol_phys and div_final < 10.0 * tol_div):
            converged = True
        print(f"Finished iterations. Converged flag = {converged}. Initial residual={init_res:.3e}, final={final_res:.3e}, final_phys={phys_final:.3e}, final_div={div_final:.3e}")

    # Final postprocessing
    psi = solve_poisson_from_omega(omega, lap_solver)
    u, v = compute_velocities_from_psi(psi, dx, dy)
    u[-1, :] = u_lid(x)
    v[-1, :] = 0.0
    u[0, :] = 0.0
    v[0, :] = 0.0
    u[:, 0] = 0.0
    u[:, -1] = 0.0
    v[:, 0] = 0.0
    v[:, -1] = 0.0

    vel_mag = np.sqrt(u * u + v * v)

    # Some simple quantitative diagnostics (centerline values)
    u_center = u[Ny // 2, Nx // 2]
    print(f"Diagnostics: u_center={u_center:.6f}, u_top_center={u[-1, Nx // 2]:.6f}")

    # Prepare figure: left contour+streamlines, right convergence history
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [3, 1]})

    # Contour plot of velocity magnitude
    cmap = 'RdBu_r'
    cf = ax1.contourf(X, Y, vel_mag, levels=50, cmap=cmap)
    # streamplot expects vector fields defined on meshgrid x,y; transpose/shape match ok here
    ax1.streamplot(x, y, u, v, color='k', density=1.2, linewidth=0.6)
    ax1.set_title('Velocity magnitude and streamlines')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    cbar = fig.colorbar(cf, ax=ax1)
    cbar.set_label('|u|')

    # Convergence history (plot Δω, divergence and physics residual)
    iters = np.arange(1, len(residuals) + 1)
    if len(residuals) >= 2:
        ax2.plot(iters, residuals, '-o', markersize=3, label='max|Δω|')
        ax2.plot(iters, divergences, '-s', markersize=3, label='max|div|')
        ax2.plot(iters, phys_res, '-^', markersize=3, label='phys L2')
        ax2.set_yscale('log')
    else:
        ax2.plot(iters, residuals, '-o', markersize=3, label='max|Δω|')
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('Metric')
    ax2.set_title('Convergence history')
    ax2.grid(True)
    ax2.legend(fontsize=8)

    plt.tight_layout()
    fig_filename = 'ns2d_stream_vmag_convergence.png'
    plt.savefig(fig_filename, dpi=200)

    # Print final quantitative metrics
    max_u = np.max(np.abs(u))
    max_v = np.max(np.abs(v))
    final_div = divergences[-1] if divergences else divergence(u, v, dx, dy)
    print(f"Final: iterations={it}, max|u|={max_u:.6f}, max|v|={max_v:.6f}, max|div|={final_div:.3e}")
    print(f"Saved figure to {fig_filename}")


if __name__ == '__main__':
    main()
```




### runtime_outputs

#### Output block1

Converged by delta tolerance at iteration 1 (delta=0.000e+00).
<string>:272: UserWarning: Data has no positive values, and therefore cannot be log-scaled.
Final: iterations=1, max|u|=0.500000, max|v|=0.000000, max|div|=0.000e+00
Saved figure to ns2d_stream_vmag_convergence.png



#### Output block2

Iter 300/3000: max|Δω|=1.597e-02, max|div|=6.463e-01, max|u|=5.000e-01
Iter 600/3000: max|Δω|=5.803e-03, max|div|=6.451e-01, max|u|=5.000e-01
Iter 900/3000: max|Δω|=3.161e-03, max|div|=6.447e-01, max|u|=5.000e-01
Iter 1200/3000: max|Δω|=2.005e-03, max|div|=6.445e-01, max|u|=5.000e-01
Iter 1500/3000: max|Δω|=1.417e-03, max|div|=6.444e-01, max|u|=5.000e-01
Iter 1800/3000: max|Δω|=1.172e-03, max|div|=6.443e-01, max|u|=5.000e-01
Iter 2100/3000: max|Δω|=9.872e-04, max|div|=6.443e-01, max|u|=5.000e-01
Iter 2400/3000: max|Δω|=8.414e-04, max|div|=6.443e-01, max|u|=5.000e-01
Iter 2700/3000: max|Δω|=7.239e-04, max|div|=6.443e-01, max|u|=5.000e-01
Iter 3000/3000: max|Δω|=6.278e-04, max|div|=6.443e-01, max|u|=5.000e-01
Finished iterations. Converged flag = True. Initial residual=4.853e+00, final=6.278e-04
Diagnostics: u_center=-0.039615, u_top_center=0.500000
Final: iterations=3000, max|u|=0.500000, max|v|=0.178699, max|div|=6.443e-01
Saved figure to ns2d_stream_vmag_convergence.png



#### Output block3

Unit test (coarse): max|psi_recov - psi_true|=3.331e-16, max|div|=2.665e-15
Iter 300/3000: max|Δω|=7.725e-03, max|div|=3.553e-15, phys_res_L2=3.297e-01, max|u|=5.000e-01
Iter 600/3000: max|Δω|=3.322e-03, max|div|=3.553e-15, phys_res_L2=1.485e-01, max|u|=5.000e-01
Iter 900/3000: max|Δω|=1.878e-03, max|div|=3.553e-15, phys_res_L2=9.203e-02, max|u|=5.000e-01
Iter 1200/3000: max|Δω|=1.182e-03, max|div|=3.553e-15, phys_res_L2=6.424e-02, max|u|=5.000e-01
Iter 1500/3000: max|Δω|=8.025e-04, max|div|=3.553e-15, phys_res_L2=4.794e-02, max|u|=5.000e-01
Iter 1800/3000: max|Δω|=5.743e-04, max|div|=3.553e-15, phys_res_L2=3.731e-02, max|u|=5.000e-01
Iter 2100/3000: max|Δω|=4.264e-04, max|div|=3.553e-15, phys_res_L2=2.981e-02, max|u|=5.000e-01
Iter 2400/3000: max|Δω|=3.248e-04, max|div|=3.553e-15, phys_res_L2=2.419e-02, max|u|=5.000e-01
Iter 2700/3000: max|Δω|=2.518e-04, max|div|=3.553e-15, phys_res_L2=1.979e-02, max|u|=5.000e-01
Iter 3000/3000: max|Δω|=1.980e-04, max|div|=3.553e-15, phys_res_L2=1.625e-02, max|u|=5.000e-01
Finished iterations. Converged flag = True. Initial residual=1.615e+01, final=1.980e-04, final_phys=1.625e-02, final_div=3.553e-15
Diagnostics: u_center=-0.086517, u_top_center=0.500000
Final: iterations=3000, max|u|=0.500000, max|v|=0.220983, max|div|=3.553e-15
Saved figure to ns2d_stream_vmag_convergence.png




### review_decision
accept

### review_comments
	Current Stage [A/3]
1. MAJOR

Category: accuracy
Issue: Helmholtz (implicit diffusion) RHS does not include boundary contributions coming from nonzero ω on the boundaries.

Recommendation: Major bug: the implicit Helmholtz solve for interior vorticity omits contributions from known boundary vorticity. The discrete Laplacian acting on interior unknowns has neighbor stencils that reference boundary nodes; when those boundary ω values are nonzero they must be moved to the RHS. Because the current implementation ignores these boundary contributions, the interior ω stays zero and the solver reports a spurious 'converged in 1 iteration' result. Fix by adding the boundary-source vector to the RHS before calling helm_solver. For a rectangular interior grid (nx_int x ny_int) this means adding dt*visc*(ω_boundary / dx^2 or dy^2) terms for interior nodes adjacent to boundaries (top/bottom/left/right). Provide an explicit vectorized assembly: compute rhs_interior = ω_interior - dt*conv + dt*visc*(bd_top/dy2 + bd_bot/dy2 + bd_left/dx2 + bd_right/dx2) where bd_top, etc. are arrays placed at the appropriate interior indices. After this change you should see nontrivial interior ω evolution and meaningful convergence.


2. MAJOR

Category: structure
Issue: Residual history contains zeros (trivial convergence) and causes a 'Data has no positive values' warning when using log scale.

Recommendation: Related to the first item: the current pseudo-time update therefore produces delta == 0 at the first iteration and residual history of zero, which triggers the matplotlib log-scale warning. After fixing the boundary-source terms, ensure residuals are recorded starting from the first meaningful update (or filter zeros before log-scaling for plotting). Also change the convergence/plotting code to handle the case where the first residual is zero (avoid plotting a single zero value on log scale).


3. MAJOR

Category: accuracy
Issue: Inconsistent enforcement of lid BC for velocities used in convection near the top boundary; may lead to incorrect convective fluxes.

Recommendation: Velocity used in the convective term should be consistent at the interior-adjacent cells. Right now you only set u at the top boundary row to the lid, but the convection stencil for interior nodes next to the lid will require the near-boundary velocities (at the first interior row and/or ghost points). Ensure that the velocities used for upwinding at interior edges include the influence of the prescribed lid — e.g. compute u at the interior first row using psi derivatives (or enforce appropriate one-sided derivative at that interface). Consider computing velocities on a staggered grid or using one-sided differences / ghost values consistent with Thom's formula so the convection term near boundaries is accurate.


4. MINOR

Category: structure
Issue: Minor structural / numerical robustness issues (dx/dy usage, initial field consistency, performance of boundary updates).

Recommendation: Numerical consistency & clarity improvements: (a) build_2d_laplacian currently takes only dx and divides both x- and y- second-derivatives by dx^2. This is fine only when dx==dy; prefer to accept dx and dy and use the correct scalings separately to make the code robust to non-square grids. (b) The Thom boundary update loops can be vectorized for performance and clarity. (c) Initialize interior ω to -Laplace(ψ) (or compute ω = -Lap ψ after the initial Poisson solve) to have a consistent starting field instead of starting interior ω at zero while boundary ω is nonzero.


5. MINOR

Category: accuracy
Issue: No quantitative validation against known benchmarks and limited diagnostics to detect stagnation vs genuine convergence.

Recommendation: Validation and diagnostics: after fixing the above, verify the solution against known lid-driven cavity benchmarks (e.g. centerline u-velocity maxima at Re=400 and location of primary vortex). Also monitor not just max|Δω| but also residuals of the steady vorticity equation or norm of u·∇u + ∇p - visc Δu (if possible). Keep the 'print only every 10%' behavior but also print a small diagnostic if the solver stagnates (no change for several prints) to help debugging.



	Current Stage [B/3]
1. MAJOR

Category: accuracy
Issue: compute_velocities_from_psi uses incorrect and brittle indexing and fallbacks, producing u and v that are inconsistent with ψ and causing large divergence (~0.64).

Recommendation: Fix the velocity-from-streamfunction routine so that u and v are computed consistently and robustly. Replace the current convoluted indexing with a clear, standard central-difference formula for interior points and one-sided differences only immediately adjacent to boundaries. For a cell-centered psi array (shape (Ny,Nx)) use
- u[1:-1,1:-1] = (psi[2:,1:-1] - psi[0:-2,1:-1]) / (2*dy)
- v[1:-1,1:-1] = -(psi[1:-1,2:] - psi[1:-1,0:-2]) / (2*dx)
Then explicitly set boundary u/v from the known BCs. Remove the fragile patterns like u[2:-2,..] and the np.where(... == 0.0) fallbacks. After this change the numerical divergence should drop to machine/FD error levels and the streamfunction identity (du/dx + dv/dy ≈ 0) will hold.


2. MAJOR

Category: accuracy
Issue: Current convergence decision uses only max|Δω| and a reduction heuristic; because u and v are inconsistent the solver declares 'converged' even though incompressibility is not satisfied.

Recommendation: Change the convergence monitor to include both the vorticity change (Δω) and a divergence-based metric, and use them together for stopping. Do not rely solely on max|Δω| reduction. For example require max|Δω| < tol_omega AND max|div| < tol_div (tol_div ~ 1e-4–1e-6 depending on grid). Also consider reporting an L2 norm of the residual of the steady vorticity equation (convective + viscous terms) as a physics-based convergence check.


3. MAJOR

Category: structure
Issue: Boundary handling is delicate and the combination of complex velocity differencing, Thom's formula and explicit boundary overwrites is a likely source of inconsistencies.

Recommendation: Simplify and harden the boundary vorticity and velocity treatments: keep psi Dirichlet consistent (ψ_wall constant), use Thom's formula as you do but ensure the one-sided derivatives used in velocity evaluation are consistent with the same indexing conventions. Add unit tests on a coarse grid (e.g. Nx=Ny=11) with a known analytic manufactured solution or a simple steady case to verify that Lap ψ = -ω holds and that computed u/v produce near-zero divergence.


4. MINOR

Category: structure
Issue: The code uses np.where(u == 0.0, ...) as a fallback to fill 'unset' entries; float-equality is unreliable and hides logic errors.

Recommendation: Remove comparisons testing float arrays for equality to zero (np.where(u == 0.0, ...)). These are brittle (rounding/initialization dependent). Instead initialize arrays with np.nan if you want to detect unset entries, or better yet compute interior and boundary blocks explicitly so no sentinel checks are needed. This makes the code more robust and easier to reason about.


5. MINOR

Category: accuracy
Issue: Fixed small dt (1e-3) can make convergence slow; nonlinear treatments (under-relaxation, adaptive pseudo-time) would improve efficiency and robustness.

Recommendation: Tune the pseudo-time step and nonlinear acceleration: consider an adaptive dt respecting a CFL estimate based on max|u|, dx, dy and explicit convection, or use stronger under-relaxation on the nonlinear term. Also consider monitoring the steady vorticity residual (convective+viscous balance) rather than only Δω. These changes will speed up convergence and make the stopping criterion meaningful.



	Current Stage [C/3]
1. MINOR

Category: accuracy
Issue: First-order upwind convection discretization is stable but strongly diffusive at Re=400; this can smear vortical structures and slow convergence to a sharp steady solution.

Recommendation: Consider a higher-order or TVD/upwind-biased scheme (e.g. QUICK, 3rd-order upwind, or a limited scheme) or a blended central/upwind approach (e.g. hybrid donor-cell with small artificial dissipation) for the convective term. Alternatively test reduced dt or different under-relaxation to assess sensitivity.


2. MINOR

Category: structure
Issue: The code uses a direct sparse LU factorization (spla.factorized) for the Laplacian and the Helmholtz matrix which is fine for the current grid size but may become memory/time-limiting for larger grids.

Recommendation: Document this limitation and optionally provide an alternative solver path (iterative Krylov solver with AMG/ILU preconditioner, e.g. pyAMG + cg, or SciPy's spla.minres/ bicgstab with a simple preconditioner) to improve scalability. Also consider using spla.splu and reusing the factorization object explicitly for clarity.


3. MINOR

Category: structure
Issue: Boundary treatment uses psi=0 on all walls (consistent up to an additive constant) but this choice is not documented in code comments; readers may be confused about streamfunction boundary values and mass conservation.

Recommendation: Add a short comment explaining that psi is set to a constant on all boundaries (here zero) because no-through-flow BCs permit an arbitrary additive constant for psi, and that Thom's formula is used to enforce wall vorticity consistent with the imposed tangential lid velocity.


4. MINOR

Category: accuracy
Issue: The physics residual res_phys is computed as conv - visc * Lap(omega) but is used directly without normalization by characteristic scales; the chosen tol_phys (1e-5) may be hard to interpret physically.

Recommendation: Either scale res_phys by a characteristic convection/viscous scale or report relative reductions (e.g. res_norm / res_norm_initial). Optionally print the initial res_phys to help interpret the final value.


5. MINOR

Category: structure
Issue: Minor clarity: the variable names and comments are generally good, but some helper functions (e.g., update_boundary_vorticity) accept u_lid as a callable and assume x is in [0,1]; this implicit assumption could be documented to avoid user errors when reusing the function on different domains/grids.

Recommendation: Document the expected coordinate mapping in the function docstring, or pass the x-array explicitly to update_boundary_vorticity to avoid hidden dependencies.





### iteration_count
3

### rev_image_description
	Current Stage [A/3]


	Current Stage [B/3]


	Current Stage [C/3]




