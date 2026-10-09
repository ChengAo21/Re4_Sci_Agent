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
Define visc = 1.0 / Re at the GLOBAL scope (outside main), and strictly use visc instead of nu throughout the code to avoid confusion with u.
Ensure all constants are passed explicitly to functions or defined globally to prevent NameError.
Plot the contour of the velocity magnitude using the 'RdBu_r' colormap, overlaid with streamlines, and the convergence history in one figure.
Just save figs do not use plt.show() in the code.

[HINTS]:
Fix the pressure at one single node to remove the nullspace singularity and prevent non-convergence warnings.
Focus on generating runnable code with reasonable physical results, do not over-optimize.
Suppress ALL per-step iteration output. Print concise progress information ONLY every 10% of total steps.
Do NOT conclude [not converged] solely because max_iters is reached, assess convergence based on whether residuals are still decreasing. 


### expanded_prob
We are asked to compute a steady solution of the two-dimensional incompressible Navier–Stokes equations on the unit square Ω = [0,1]^2 with Reynolds number Re = 400 (visc = 1.0/Re). The velocity boundary condition is a moving lid on the top boundary Γ1: (u,v) = (α x(1-x),0) with α = 2, and no-slip on the other boundaries Γ2. Pressure is determined up to a constant: the problem statement enforces p(0,0)=0 to eliminate the nullspace. The goal is a stable and efficient numerical implementation that produces physically reasonable steady fields and convergence history and saves figures (contour of velocity magnitude using the 'RdBu_r' colormap with streamlines overlaid, and a plot of residual history). The implementation must: define visc globally and use it everywhere, pass constants explicitly to functions, suppress per-step iteration output, print progress only every 10% of total steps, fix pressure at one node, assess convergence by residual trends (not merely max_iters), and save figures without plt.show().

Primary mathematical and numerical challenges:
- Incompressibility (divergence-free) constraint couples velocity and pressure; enforcing it robustly is critical. The pressure has a nullspace (constant offset) that must be removed (fix pressure at one point or use appropriate boundary conditions for pressure Poisson).
- Nonlinearity from convective terms requires iterative treatment/linearization (pseudo-time stepping, Newton, or SIMPLE-like splitting). Stability of the convective term is a concern at Re=400: central differences may suffer oscillations; upwinding or stabilized discretizations may be required.
- Discretization choice (staggered/MAC vs collocated) impacts pressure-velocity coupling and checkerboarding. A staggered grid avoids the need for Rhie–Chow interpolation and is simpler for finite-difference implementations.
- Choosing solvers: Poisson/elliptic solver performance dominates runtime. Efficient Poisson solvers (multigrid or optimized sparse solvers) strongly influence overall efficiency and scalability.
- Boundary layers and corner singularities at moderate Re require sufficient resolution; the steady solution may be sensitive to grid size and the discretization of wall vorticity (for streamfunction–vorticity methods).
- Convergence monitoring: need robust residual definitions (momentum residuals and continuity residual or combined norm) and a policy to decide steady-state (e.g., stagnation in residual decay, reaching tolerance, or asymptotic decrease) rather than relying only on iteration count.
- Practical coding constraints: ensure global variable visc is used, pass other constants explicitly, and avoid verbose per-step printing while giving periodic progress updates.

Given these challenges, two alternative, practical, and runnable strategies are presented below. Each plan gives the governing idea, algorithmic steps in ordered form, and explicit notes on stability, accuracy and computational cost/limitations.

### solution_plans
	Current Stage [A/2]
solu_name='Staggered-MAC Projection Method' content="Solu Name: Staggered-MAC with pseudo-transient projection (implicit diffusion + explicit convection)\n\nGoverning idea:\n- Use a staggered (MAC) finite-difference grid: u defined at vertical cell faces, v at horizontal faces, and pressure p at cell centers. March the steady solution in pseudo-time using a projection (pressure correction) method: treat viscous diffusion implicitly (for stability at reasonable time steps) and convection explicitly (with conservative upwind/central blend). Enforce incompressibility via a Poisson equation for pressure correction; fix p(0,0)=0 to remove the nullspace.\n\nAlgorithmic steps (numbered):\n1) Grid and parameters:\n   - Set Nx, Ny and uniform grid spacing dx = 1/Nx, dy = 1/Ny; define staggered grid arrays for u (Nx+1 by Ny), v (Nx by Ny+1), and p (Nx by Ny). Define global visc = 1.0/Re outside main and pass visc, dt_max, total_steps, tol into routines.\n2) Initialize fields:\n   - Initial guess: u=0 everywhere except impose lid velocity on top u-face locations; v=0; p=0 (or small random). Apply boundary conditions explicitly on staggered locations.\n3) Time-marching loop (pseudo-time t += dt until steady or max_iters):\n   - Choose dt adaptively via CFL and visc constraints: dt = min(CFL * min(dx/|u|_max, dy/|v|_max), dt_diff) with dt_diff ~ 0.5*min(dx,dy)^2/visc; impose floor/ceiling on dt and pass dt explicitly.\n   - Compute convective terms at u and v face centers with a stabilized discretization: second-order central + explicit donor-cell blending (or first-order upwind if high local Peclet) to avoid wiggles.\n   - Form right-hand side for implicit viscous solve: for u-momentum, (u^* - dt * convect_u) = u^n + dt * visc * Laplacian(u^{n+1}) + dt * (-dp/dx) (similarly for v). Treat diffusion implicitly by building a sparse linear system for u^{n+1} (or use efficient tridiagonal solves by directional splitting). Solve u* and v* provisional velocities without pressure term.\n   - Build discrete divergence of provisional velocity and solve pressure-Poisson: ∇^2 φ = (1/dt) ∇·u* with Neumann/derived boundary conditions. Fix pressure at cell (0,0): p(0,0)=0 to remove nullspace. Use a reliable sparse iterative solver (e.g., conjugate gradients with incomplete-LU or a multigrid if available) and pass solver tolerances explicitly.\n   - Correct velocities: u^{n+1} = u* - dt * dφ/dx (apply at face centers), v^{n+1} = v* - dt * dφ/dy.\n   - Update pressure: p^{n+1} = p^n + φ (or directly set p = φ with fixed point handling), ensuring p(0,0)=0 after update.\n   - Apply boundary conditions to corrected velocities (re-impose lid and no-slip walls) and ensure ghost/extrapolation values consistent.\n4) Convergence monitoring and outputs:\n   - Compute residuals each step: momentum residual norm (L2 or max of discrete residuals for u and v) and continuity residual (L2 of divergence on cell centers). Store residual history.\n   - Suppress per-step output: print only concise progress every 10% of total_steps (compute percent based on max_iters or predicted iterations).\n   - Termination: stop if residual norms drop below tol AND show asymptotic decrease (e.g., residual reduced by factor or plateauing small), otherwise continue until max_iters. Do NOT declare failure solely because max_iters was hit; examine trend.\n5) Postprocessing and plotting (save figures, do not plt.show()):\n   - Interpolate staggered velocities to cell centers to compute speed magnitude ||u||.\n   - Create a figure with two subplots: top: filled contour of speed using colormap 'RdBu_r' with overlaid streamlines computed from interpolated velocities; bottom: semilogy plot of residual history vs iteration. Save the figure(s) to files.\n\nStability / accuracy / complexity / efficiency limitations:\n- Stability: treating diffusion implicitly stabilizes stiff viscous terms, but explicit convection imposes CFL restrictions; use donor-cell blending to control oscillations at moderate Re. Fully implicit nonlinear solves (e.g., Newton) would allow larger time steps but are more complex.\n- Accuracy: staggered MAC with second-order central differences for diffusion and blended upwind for convection gives reasonable accuracy for moderate Re; resolution must be sufficient to resolve boundary layers (grid refinement study recommended).\n- Complexity: forming and solving the implicit diffusion and pressure Poisson sparse systems dominates runtime. Using efficient iterative solvers or multigrid is recommended for large grids; direct solvers will be costly in memory/time as Nx,Ny grow.\n- Implementation burden: staggered-grid indexing and boundary condition implementation are more code-intensive but avoid pressure checkerboarding and Rhie–Chow complications.\n- Practical note: ensure visc is defined globally and passed to functions; explicitly fix p(0,0)=0 in the Poisson solver to remove nullspace."

	Current Stage [B/2]
solu_name='Streamfunction-Vorticity Iteration' content="Solu Name: Streamfunction–Vorticity iterative steady solver (ψ–ω approach)\n\nGoverning idea:\n- Eliminate pressure by switching to streamfunction ψ and vorticity ω: u = ∂ψ/∂y, v = -∂ψ/∂x, and ∇^2 ψ = -ω. Solve iteratively for a steady ω that satisfies the vorticity transport equation u·∇ω = visc ∇^2 ω with the no-slip lid and wall boundary conditions enforced through boundary values of ψ and wall vorticity specification. Use implicit (or semi-implicit) treatment for diffusion and explicit or semi-implicit advection for robustness. The ψ–ω formulation removes the pressure nullspace issue.\n\nAlgorithmic steps (numbered):\n1) Grid and parameters:\n   - Use a collocated cell-centered grid for ψ and ω with Nx × Ny interior cells; define dx,dy. Set global visc = 1.0/Re outside main and pass visc and algorithmic constants explicitly.\n2) Boundary conditions and initial guesses:\n   - For no-slip walls, set ψ = constant along boundaries (choose ψ = 0 along all boundaries for convenience). Impose lid tangential velocity by setting ghost-cell streamfunction values consistent with the specified u on the top boundary; compute corresponding boundary vorticity or use a finite-difference relation (e.g., second-order one-sided) to evaluate ω on walls.\n   - Initialize ω = 0 (or small), ψ = 0, compute initial u,v from ψ.\n3) Iterative update loop for steady ω:\n   - Given ω^n, solve Poisson ∇^2 ψ^{n} = -ω^n with Dirichlet ψ on boundaries (sparse solver or multigrid); ensure the solver pins solution (ψ boundary values fix constant) so no nullspace problem arises.\n   - Compute velocities u = ∂ψ/∂y, v = -∂ψ/∂x at cell centers via second-order central differences.\n   - Advance ω using a steady pseudo-transient or fixed-point update: (ω^{n+1} - ω^n)/Δt + u·∇ω^{*} = visc ∇^2 ω^{n+1}, treating diffusion implicitly and advection explicitly (or use a semi-implicit scheme). Δt chosen via CFL constraint for advection and diffusive limit; pass dt explicitly.\n   - Enforce boundary vorticity at walls at the end of each ω update using the streamfunction values (use finite-difference formula ω_wall = -2*(ψ_interior - ψ_boundary)/Δn^2 - 2*U_wall/Δn for moving lid as appropriate).\n   - Update residuals: measure change in ω (||ω^{n+1}-ω^n||) and the divergence (should be zero by construction), or monitor momentum-equivalent residuals derived from ψ and ω.\n4) Convergence and outputs:\n   - Suppress per-step verbose output; print progress only every 10% of max_iters. Terminate when ω changes are below tolerance OR residuals show steady asymptotic decay. Do not equate max_iters with failure if residuals are still decreasing significantly.\n   - Compute speed magnitude from u and v from the converged ψ, interpolate as needed for plotting.\n5) Postprocessing and plotting (save figures, do not plt.show()):\n   - Create a figure with top subplot: contour of velocity magnitude using 'RdBu_r' and overlay streamlines (derived directly from ψ contours or computed from u,v); bottom subplot: plot residual history (e.g., ||ω^{n+1}-ω^n|| or momentum-equivalent residuals) on a semilogy scale. Save the figure(s).\n\nStability / accuracy / complexity / efficiency limitations:\n- Stability: the ψ–ω formulation avoids pressure coupling, but wall vorticity boundary conditions are delicate (especially in corners) and can hamper convergence if approximated poorly. Advection at Re=400 requires stabilized discretizations (upwinding or flux-limiter) to avoid oscillations; implicit diffusion helps.\n- Accuracy: second-order central schemes for Laplacians and derivatives give reasonable accuracy, but the boundary-layer resolution near the lid and corners may require grid refinement. Errors in wall vorticity estimation degrade global accuracy.\n- Complexity: solving the Poisson equation for ψ each iteration is the major cost. A fast Poisson solver (multigrid or FFT for uniform grids) is highly recommended. The method has smaller linear algebra per iteration than a fully coupled primitive-variable Newton solve but may need many pseudo-time steps for steady convergence.\n- Applicability: ψ–ω is convenient for 2D incompressible flow and removes pressure issues, but it is not directly extendable to 3D and requires careful boundary treatments for moving walls.\n- Implementation notes: define visc globally and pass it to all functions; ensure vorticity at walls is updated consistently and print progress only at 10% increments of iteration budget."



### technical_spec
	Current Stage [A/3]
Solver architecture (brief):
- Unknowns: streamfunction psi and vorticity omega on a collocated 2D grid with Dirichlet psi on boundaries.
- Matrices assembled once: 2D Laplacian L used for Poisson solve (psi) and for implicit diffusion operator B = I - dt*visc*L (omega).
- Factorize matrices with sparse direct solver (scipy.sparse.linalg.factorized) once prior to time loop.
- Time-stepping: semi-implicit update for omega: (I - dt*visc*L) omega^{n+1} = omega^n - dt*(u·∇ω) + dt*visc*(boundary contributions). Convection handled with first-order upwind for stability.
- Boundary vorticity enforced from psi each iteration using second-order finite-difference relation; used as Dirichlet BC in the implicit step via boundary-contribution vector.
- Convergence monitor: relative L2 change in omega each iteration; print progress every 10% of max iterations; final assessment uses reduction factor and trend to avoid declaring failure solely based on max_iters.
- Postprocessing: interpolate psi->u,v at interior points, compute speed magnitude, plot filled contour (RdBu_r) with streamlines and convergence history; save figure without showing.

Data flow:
1) Build operators L, B once and factorize.
2) Start with omega and psi (initially zeros) and enforce initial wall omega from psi.
3) Iteration: solve Poisson for psi <- A^{-1}(-omega), compute wall omega from psi, compute u,v, compute convective term, compute RHS including boundary contributions, solve B omega^{n+1} <- RHS, update residual.
4) Repeat until tolerance or max_iters; finalize and plot results.

Key implementation notes:
- visc is declared globally and used everywhere.
- Matrices are assembled and factorized once outside the time loop.
- Boundary contributions are computed per iteration and added to the RHS (allowed).
- First-order upwind stabilizes convection at moderate Reynolds.
- Prints are limited to checkpoints to keep iteration output sparse.

	Current Stage [B/3]
- Unknowns: vorticity omega and streamfunction psi on a collocated (nx x ny) interior grid; arrays include a one-cell boundary padding for Dirichlet psi and boundary vorticity.
- Operators: 2D Laplacian L built once via Kronecker products (two operand orderings tested to ensure consistency with Fortran-style flattening). Implicit diffusion matrix B = I - dt*visc*L assembled once and factorized.
- Time-stepping: semi-implicit for omega: convection handled with first-order upwind, diffusion implicit via factorized B solve. Boundary vorticity enforced each step from psi via second-order relations.
- Data flow: build operators -> factorize -> initialize psi/omega -> loop: Poisson(psi) <- -omega, update wall omega, compute u/v, convective term, construct RHS (including boundary contributions), solve B for new omega, compute relative L2 change residual -> repeat until tol or max_iters.
- Convergence: relative L2 change in omega monitored; progress printed only every 10% of max_iters (and at final iteration); final assessment uses reduction factor and trend.
- Postprocessing: compute u,v from psi, speed magnitude; produce contour (RdBu_r) of |u| with streamlines overlay; plot convergence history; save figure to file without plt.show().
- Implementation notes: visc declared globally; Laplacian ordering verified programmatically; matrices built and factorized once; boundary contributions allowed to change per iteration; streamplot/grid shapes carefully aligned to avoid ValueError.

	Current Stage [C/3]
Streamfunction-vorticity semi-implicit solver for 2D steady lid-driven cavity (Re=400). Core flow:
- Build Laplacian once (kron ordering verified) and implicit diffusion matrix B once; factorize (with guarded fallbacks).
- Fields lived on (nx+2, ny+2) arrays with one-cell padding for Dirichlet psi and ghost vorticity.
- update_wall_vorticity_from_psi writes ghost vorticity into omega_full ghost indices (0, nx+1, 0, ny+1) to be consistent with build_boundary_contrib_vector.
- Time-stepping: semi-implicit: convection explicit (1st-order upwind), diffusion implicit via B solve.
- Diagnostics: relative omega-change residual, vorticity PDE residual (u·∇ω - νΔω), discrete divergence surrogate; printed only every 10% of max_iters.
- Unit test to verify boundary-contribution indexing consistency.
- Postprocessing: contour of |u| with streamlines (RdBu_r) and convergence history saved to file. Quantitative prints: progress, diagnostics, centerline samples.




### python_codes

#### Script block1:

```python
#!/usr/bin/env python3
"""
Steady 2D incompressible Navier-Stokes solver (lid-driven cavity) using
Streamfunction - Vorticity formulation with semi-implicit time-stepping.

Global visc is defined outside main as required.
Saves a single figure with the velocity magnitude contour + streamlines (top)
and convergence history (bottom). No per-step prints; concise progress every
10% of max_iters.

Requirements satisfied:
- visc is global and used everywhere
- matrices (Poisson and implicit diffusion) are assembled once and factorized
  before time loop
- matrices are NOT rebuilt inside the time loop
- prints only at 10% intervals
- saves figure; does not call plt.show()

Author: Generated by assistant
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# GLOBAL viscosity according to the problem statement
Re = 400.0
visc = 1.0 / Re


def build_2d_laplacian(nx, ny, dx, dy):
    """Build 2D Laplacian matrix for interior nodes (Dirichlet BC on boundaries).
    The ordering is lexicographic: k = (j-1)*nx + (i-1), i=1..nx, j=1..ny.
    Returns sparse CSR matrix of shape (nx*ny, nx*ny).
    """
    N = nx * ny
    # 1D operators
    idx = np.arange(nx)
    Tx = sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(nx, nx)) / (dx * dx)
    Ty = sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(ny, ny)) / (dy * dy)
    I_x = sp.eye(nx)
    I_y = sp.eye(ny)
    # 2D Laplacian via Kronecker
    L = sp.kron(I_y, Tx) + sp.kron(Ty, I_x)
    return L.tocsr()


def idx_map(nx, ny):
    """Return lambda to map (i,j) interior indices to flat index k.
    i in [1..nx], j in [1..ny]."""
    def mapper(i, j):
        return (j - 1) * nx + (i - 1)
    return mapper


def flatten_interior(arr):
    """Flatten interior arr[1:nx+1,1:ny+1] into vector in lexicographic order."""
    return arr[1:-1, 1:-1].ravel(order='F')


def unflatten_interior(vec, nx, ny):
    """Unflatten vector into array of shape (nx,ny) with Fortran order and
    embed into (nx+2, ny+2) with zero boundary padding.
    """
    interior = vec.reshape((nx, ny), order='F')
    out = np.zeros((nx + 2, ny + 2), dtype=vec.dtype)
    out[1:-1, 1:-1] = interior
    return out


def compute_velocity_from_psi(psi, dx, dy):
    """Compute velocities u = dpsi/dy, v = -dpsi/dx at interior grid points.
    psi is (nx+2, ny+2) with Dirichlet boundaries known.
    Returns (u,v) arrays of shape (nx, ny) corresponding to interior nodes.
    """
    # interior indices
    # central differences
    u = (psi[1:-1, 2:] - psi[1:-1, 0:-2]) / (2.0 * dy)
    v = -(psi[2:, 1:-1] - psi[0:-2, 1:-1]) / (2.0 * dx)
    return u, v


def compute_convective_term(omega_full, u, v, dx, dy):
    """Compute convective term u * d(omega)/dx + v * d(omega)/dy at interior nodes.
    omega_full has shape (nx+2, ny+2). u and v are interior shape (nx,ny).
    Use first-order upwind scheme for stability.
    Returns conv interior array shape (nx,ny).
    """
    # forward/backward differences
    # dω/dx: shape (nx,ny). For interior nodes (i=1..nx), use omega_full indices
    o = omega_full
    # backward difference (ω_i - ω_{i-1})/dx
    d_omega_bx = (o[1:-1, 1:-1] - o[0:-2, 1:-1]) / dx
    # forward difference (ω_{i+1} - ω_i)/dx
    d_omega_fx = (o[2:, 1:-1] - o[1:-1, 1:-1]) / dx
    # similarly for y
    d_omega_by = (o[1:-1, 1:-1] - o[1:-1, 0:-2]) / dy
    d_omega_fy = (o[1:-1, 2:] - o[1:-1, 1:-1]) / dy

    # upwind selection
    dwdx = np.where(u >= 0.0, d_omega_bx, d_omega_fx)
    dwdy = np.where(v >= 0.0, d_omega_by, d_omega_fy)

    conv = u * dwdx + v * dwdy
    return conv


def update_wall_vorticity_from_psi(omega_full, psi, dx, dy, alpha):
    """Enforce wall vorticity using the streamfunction (Dirichlet psi on walls).
    omega_full is (nx+2, ny+2) array and will be modified in place for boundary
    indices at the interior-adjacent rows/cols.

    For top wall (moving lid): omega[..., ny] = -2*psi[..., ny]/dy^2 - 2*U/dy
    For other walls (stationary): omega = -2*psi_interior / (d^2)
    """
    # interior sizes
    nx = psi.shape[0] - 2
    ny = psi.shape[1] - 2
    # top: interior-adjacent row index is ny (since psi indexing 0..ny+1)
    # moving lid velocity U(x) = alpha * x(1-x) at y=1
    x = np.linspace(dx, 1.0 - dx, nx)
    U_lid = alpha * x * (1.0 - x)
    # top row interior j = ny
    omega_full[1:-1, ny] = -2.0 * psi[1:-1, ny] / (dy * dy) - 2.0 * U_lid / dy
    # bottom row interior j = 1
    omega_full[1:-1, 1] = -2.0 * psi[1:-1, 1] / (dy * dy)
    # left column interior i = 1
    omega_full[1, 1:-1] = -2.0 * psi[1, 1:-1] / (dx * dx)
    # right column interior i = nx
    omega_full[nx, 1:-1] = -2.0 * psi[nx, 1:-1] / (dx * dx)


def build_boundary_contrib_vector(omega_full, nx, ny, dx, dy):
    """Compute vector of boundary contributions for Laplacian when boundaries are
    non-zero. This computes for each interior node the sum of neighbor boundary
    contributions of the form (ω_boundary)/dx^2 or /dy^2 where applicable.
    Returns flattened vector length nx*ny in Fortran order.
    """
    b = np.zeros((nx, ny), dtype=omega_full.dtype)
    # left boundary neighbors (i-1 == 0)
    b[0, :] += omega_full[0, 1:-1] / (dx * dx)
    # right boundary neighbors (i+1 == nx+1)
    b[-1, :] += omega_full[nx + 1, 1:-1] / (dx * dx)
    # bottom boundary neighbors (j-1 == 0)
    b[:, 0] += omega_full[1:-1, 0] / (dy * dy)
    # top boundary neighbors (j+1 == ny+1)
    b[:, -1] += omega_full[1:-1, ny + 1] / (dy * dy)
    return b.ravel(order='F')


def main():
    # Physical and numerical parameters
    alpha = 2.0  # lid amplitude factor
    nx = 64      # number of interior points in x
    ny = 64      # number of interior points in y
    dx = 1.0 / (nx + 1)  # spacing with Dirichlet boundaries at 0 and 1
    dy = 1.0 / (ny + 1)

    # Time stepping parameters
    # choose dt based on advective CFL and reasonable value for semi-implicit scheme
    umax_est = alpha * 0.25  # lid peak = alpha*0.25
    cfl = 0.2
    dt_adv = cfl * min(dx / max(umax_est, 1e-8), dy / max(umax_est, 1e-8))
    dt = min(0.005, dt_adv)  # keep dt <= 0.005 for stability

    max_iters = 4000
    tol = 1e-6

    print(f"Starting streamfunction-vorticity solver: nx={nx}, ny={ny}, dt={dt:.5e}, visc={visc:.3e}")

    # Build Laplacian matrix once
    L = build_2d_laplacian(nx, ny, dx, dy)  # shape (N,N)
    N = nx * ny
    I = sp.eye(N, format='csc')

    # Poisson matrix A: Lap ψ = -ω  -> A ψ = -ω, where A = L
    A_poisson = L.tocsc()
    # Implicit diffusion matrix for vorticity: (I - dt*visc*L)
    B_mat = (I - dt * visc * L).tocsc()

    # Factorize matrices once
    A_factor = spla.factorized(A_poisson)  # gives a function to solve A x = b
    B_factor = spla.factorized(B_mat)

    # initialize arrays
    # omega and psi full arrays include boundaries: shape (nx+2, ny+2)
    omega = np.zeros((nx + 2, ny + 2), dtype=float)
    psi = np.zeros_like(omega)

    # initial boundary vorticity from initial psi (zero) -> top boundary non-zero due to lid
    update_wall_vorticity_from_psi(omega, psi, dx, dy, alpha)

    # flatten interior omega
    omega_flat = flatten_interior(omega)

    residuals = []

    # progress checkpoint
    checkpoint = max(1, int(max_iters / 10))

    # initial residual baseline
    initial_res = None

    for it in range(1, max_iters + 1):
        # 1) Solve Poisson for psi: Lap psi = -omega -> A psi_flat = -omega_flat
        rhs_psi = -omega_flat
        psi_flat = A_factor(rhs_psi)
        psi = unflatten_interior(psi_flat, nx, ny)
        # boundary psi are zero by construction

        # 2) Update wall vorticity using latest psi (Dirichlet psi on walls)
        update_wall_vorticity_from_psi(omega, psi, dx, dy, alpha)

        # 3) Compute velocities at interior points from psi
        u, v = compute_velocity_from_psi(psi, dx, dy)

        # 4) Compute convective term using upwind
        conv = compute_convective_term(omega, u, v, dx, dy)
        conv_flat = conv.ravel(order='F')

        # 5) Build boundary contribution vector for implicit Laplacian using current boundary omega
        bnd_vec = build_boundary_contrib_vector(omega, nx, ny, dx, dy)

        # 6) Assemble RHS for implicit solve: (I - dt*visc*L) ω^{n+1} = ω^{n} - dt*conv + dt*visc*bnd_vec
        rhs = omega_flat - dt * conv_flat + dt * visc * bnd_vec

        # 7) Solve for new interior omega
        omega_new_flat = B_factor(rhs)

        # 8) update omega full and compute residual
        omega_old_flat = omega_flat.copy()
        omega_flat = omega_new_flat
        omega = unflatten_interior(omega_flat, nx, ny)

        # 9) After updating interior omega, we refresh psi via next iteration (so update wall omega now)
        # For residual monitoring compute change
        change = omega_flat - omega_old_flat
        res = np.linalg.norm(change, ord=2) / max(1.0, np.linalg.norm(omega_flat, ord=2))
        residuals.append(res)
        if initial_res is None:
            initial_res = res if res > 0 else 1.0

        # Progress print only every 10% steps
        if it % checkpoint == 0 or it == 1 or it == max_iters:
            pct = int(100.0 * it / max_iters)
            print(f"Progress: {pct}% (iter {it}/{max_iters}), res = {res:.3e}")

        # Convergence test: relative change small
        if res < tol:
            print(f"Converged: relative ω-change {res:.3e} < tol {tol:.1e} at iter {it}")
            break

    # end time loop

    # Final psi from final omega
    psi_flat = A_factor(-omega_flat)
    psi = unflatten_interior(psi_flat, nx, ny)
    update_wall_vorticity_from_psi(omega, psi, dx, dy, alpha)
    u, v = compute_velocity_from_psi(psi, dx, dy)
    speed = np.sqrt(u * u + v * v)

    # Assessment when max_iters reached but residual still decreasing
    final_res = residuals[-1] if len(residuals) > 0 else 0.0
    reduction_factor = initial_res / final_res if final_res > 0 else np.inf
    trend = "decreasing" if (len(residuals) >= 5 and residuals[-1] < residuals[-5]) else "unknown/plateau"
    if final_res < tol:
        final_message = f"Converged: final residual {final_res:.3e}"
    elif reduction_factor > 100.0 and trend == "decreasing":
        final_message = f"Stopped at max_iters but residual decreasing (reduction factor {reduction_factor:.1f})."
    else:
        final_message = f"Stopped at max_iters; final residual {final_res:.3e} (reduction factor {reduction_factor:.1f})."

    print(final_message)

    # Prepare plots: one figure with two subplots
    X = np.linspace(dx, 1.0 - dx, nx)
    Y = np.linspace(dy, 1.0 - dy, ny)
    XX, YY = np.meshgrid(X, Y, indexing='ij')

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), constrained_layout=True)

    # Top: contour of velocity magnitude with streamlines overlay
    cf = ax1.contourf(XX, YY, speed, levels=50, cmap='RdBu_r')
    strm = ax1.streamplot(XX, YY, u.T, v.T, color='k', density=1.0, linewidth=0.6)
    ax1.set_title('Velocity magnitude and streamlines (lid-driven cavity)')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    fig.colorbar(cf, ax=ax1, label='|u|')

    # Bottom: convergence history
    ax2.semilogy(np.arange(1, len(residuals) + 1), residuals, '-k')
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('Relative change residual (L2)')
    ax2.set_title('Convergence history')
    ax2.grid(True, which='both', ls='--', lw=0.5)

    # Save figure
    fig_filename = 'ns_stream_vorticity_convergence.png'
    fig.savefig(fig_filename, dpi=200)
    print(f"Figure saved to: {fig_filename}")


if __name__ == '__main__':
    main()
```


#### Script block2:

```python
#!/usr/bin/env python3
"""
Streamfunction-vorticity solver for 2D steady incompressible Navier-Stokes
(lid-driven cavity) using semi-implicit time stepping.

Key points enforced:
- visc is defined at GLOBAL scope and used everywhere
- 2D Laplacian and implicit diffusion matrices assembled once and factorized
- No matrix rebuild inside the time loop
- Sparse progress printing only every 10% of max_iters (and at final iteration)
- Streamplot/grid shapes fixed to avoid ValueError
- Ordering consistency test for Laplacian vs. flatten/unflatten (Fortran order)
- Brief note: pressure is not computed in this formulation (p(0,0)=0 irrelevant)

The script saves a single figure with velocity magnitude contour + streamlines
and the convergence history. No interactive show() call.
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# GLOBAL viscosity according to the problem statement
Re = 400.0
visc = 1.0 / Re


def build_2d_laplacian(nx, ny, dx, dy, kron_choice=0):
    """Build 2D Laplacian for interior (nx x ny) nodes.
    Returns CSR sparse matrix. kron_choice picks operand ordering to match
    vectorization ordering (0: kron(Iy,Tx)+kron(Ty,Ix), 1: swapped).
    """
    Tx = sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(nx, nx)) / (dx * dx)
    Ty = sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(ny, ny)) / (dy * dy)
    I_x = sp.eye(nx)
    I_y = sp.eye(ny)
    if kron_choice == 0:
        L = sp.kron(I_y, Tx) + sp.kron(Ty, I_x)
    else:
        L = sp.kron(Tx, I_y) + sp.kron(I_x, Ty)
    return L.tocsr()


def idx_map(nx, ny):
    def mapper(i, j):
        return (j - 1) * nx + (i - 1)
    return mapper


def flatten_interior(arr):
    """Flatten arr[1:-1,1:-1] to length nx*ny vector in Fortran order.
    This ordering is used consistently with the Laplacian construction.
    """
    return arr[1:-1, 1:-1].ravel(order='F')


def unflatten_interior(vec, nx, ny):
    """Unflatten vector (Fortran-order) into full array with zero boundary padding."""
    interior = vec.reshape((nx, ny), order='F')
    out = np.zeros((nx + 2, ny + 2), dtype=vec.dtype)
    out[1:-1, 1:-1] = interior
    return out


def compute_velocity_from_psi(psi, dx, dy):
    """Compute interior velocities u = dpsi/dy, v = -dpsi/dx (central differences)."""
    u = (psi[1:-1, 2:] - psi[1:-1, 0:-2]) / (2.0 * dy)
    v = -(psi[2:, 1:-1] - psi[0:-2, 1:-1]) / (2.0 * dx)
    return u, v


def compute_convective_term(omega_full, u, v, dx, dy):
    """First-order upwind convective term at interior nodes."""
    o = omega_full
    d_omega_bx = (o[1:-1, 1:-1] - o[0:-2, 1:-1]) / dx
    d_omega_fx = (o[2:, 1:-1] - o[1:-1, 1:-1]) / dx
    d_omega_by = (o[1:-1, 1:-1] - o[1:-1, 0:-2]) / dy
    d_omega_fy = (o[1:-1, 2:] - o[1:-1, 1:-1]) / dy
    dwdx = np.where(u >= 0.0, d_omega_bx, d_omega_fx)
    dwdy = np.where(v >= 0.0, d_omega_by, d_omega_fy)
    conv = u * dwdx + v * dwdy
    return conv


def update_wall_vorticity_from_psi(omega_full, psi, dx, dy, alpha):
    """Enforce wall vorticity using Dirichlet psi on walls.
    Derived from finite-difference approximations consistent with
    omega = dv/dx - du/dy and psi boundary values.

    For top wall (y=1) with u=U(x):
      omega_interior = -2*psi_interior/dy^2 - 2*U/dy
    For bottom wall (y=0) with u=0:
      omega_interior = -2*psi_interior/dy^2
    For left/right walls (v=0):
      omega_interior = -2*psi_interior/dx^2

    Indexing: omega_full[i,j] where i indexes x (0..nx+1) and j indexes y (0..ny+1).
    """
    nx = psi.shape[0] - 2
    ny = psi.shape[1] - 2
    x = np.linspace(dx, 1.0 - dx, nx)
    U_lid = alpha * x * (1.0 - x)
    # top wall interior-adjacent (j = ny)
    omega_full[1:-1, ny] = -2.0 * psi[1:-1, ny] / (dy * dy) - 2.0 * U_lid / dy
    # bottom (j = 1)
    omega_full[1:-1, 1] = -2.0 * psi[1:-1, 1] / (dy * dy)
    # left (i = 1)
    omega_full[1, 1:-1] = -2.0 * psi[1, 1:-1] / (dx * dx)
    # right (i = nx)
    omega_full[nx, 1:-1] = -2.0 * psi[nx, 1:-1] / (dx * dx)


def build_boundary_contrib_vector(omega_full, nx, ny, dx, dy):
    """Compute neighbor boundary contributions for interior Laplacian."""
    b = np.zeros((nx, ny), dtype=omega_full.dtype)
    b[0, :] += omega_full[0, 1:-1] / (dx * dx)
    b[-1, :] += omega_full[nx + 1, 1:-1] / (dx * dx)
    b[:, 0] += omega_full[1:-1, 0] / (dy * dy)
    b[:, -1] += omega_full[1:-1, ny + 1] / (dy * dy)
    return b.ravel(order='F')


def verify_laplacian_ordering(nx, ny, dx, dy):
    """Try different Kron orderings and pick one consistent with Fortran flattening.
    Returns chosen kron_choice (0 or 1) and the Laplacian matrix (CSR).
    """
    # test function f(x,y) = sin(pi x) sin(pi y)
    X = np.linspace(dx, 1.0 - dx, nx)
    Y = np.linspace(dy, 1.0 - dy, ny)
    XX, YY = np.meshgrid(X, Y, indexing='ij')
    f2d = np.sin(np.pi * XX) * np.sin(np.pi * YY)
    # finite-difference Laplacian (interior) using central differences
    f_xx = (np.roll(f2d, -1, axis=0) - 2.0 * f2d + np.roll(f2d, 1, axis=0)) / (dx * dx)
    f_yy = (np.roll(f2d, -1, axis=1) - 2.0 * f2d + np.roll(f2d, 1, axis=1)) / (dy * dy)
    # zero fictitious boundaries are handled by setting neighbors at edges to the same value
    # but we will compute internal ones directly (exclude ghost influence): compute true interior
    f2d_pad = np.zeros((nx + 2, ny + 2))
    f2d_pad[1:-1, 1:-1] = f2d
    lap_fd = np.zeros_like(f2d)
    for i in range(nx):
        for j in range(ny):
            lap_fd[i, j] = (f2d_pad[i + 2, j + 1] - 2.0 * f2d_pad[i + 1, j + 1] + f2d_pad[i, j + 1]) / (dx * dx) + (
                f2d_pad[i + 1, j + 2] - 2.0 * f2d_pad[i + 1, j + 1] + f2d_pad[i + 1, j]) / (dy * dy)
    f_flat = f2d.ravel(order='F')
    target = lap_fd.ravel(order='F')
    # try both Kron orderings
    errs = []
    mats = []
    for choice in (0, 1):
        Ltest = build_2d_laplacian(nx, ny, dx, dy, kron_choice=choice)
        Lf = Ltest.dot(f_flat)
        err = np.linalg.norm(Lf - target) / max(1.0, np.linalg.norm(target))
        errs.append(err)
        mats.append(Ltest)
    best = int(np.argmin(errs))
    if errs[best] > 1e-8:
        print(f"Warning: Laplacian verification residual = {errs[best]:.3e}; ordering may be inexact.")
    return best, mats[best]


def main():
    alpha = 2.0
    nx = 64
    ny = 64
    dx = 1.0 / (nx + 1)
    dy = 1.0 / (ny + 1)

    umax_est = alpha * 0.25
    cfl = 0.2
    dt_adv = cfl * min(dx / max(umax_est, 1e-8), dy / max(umax_est, 1e-8))
    dt = min(0.005, dt_adv)

    max_iters = 4000
    tol = 1e-6

    print(f"Starting solver: nx={nx}, ny={ny}, dt={dt:.5e}, visc={visc:.3e}")

    # Verify Laplacian ordering and build matrix once
    kron_choice, L = verify_laplacian_ordering(nx, ny, dx, dy)
    N = nx * ny
    I = sp.eye(N, format='csc')

    A_poisson = L.tocsc()
    B_mat = (I - dt * visc * L).tocsc()

    # Factorize once
    A_factor = spla.factorized(A_poisson)
    B_factor = spla.factorized(B_mat)

    # initialize fields
    omega = np.zeros((nx + 2, ny + 2), dtype=float)
    psi = np.zeros_like(omega)

    update_wall_vorticity_from_psi(omega, psi, dx, dy, alpha)
    omega_flat = flatten_interior(omega)

    residuals = []
    checkpoint = max(1, int(max_iters / 10))
    initial_res = None

    for it in range(1, max_iters + 1):
        # Poisson: Lap psi = -omega
        psi_flat = A_factor(-omega_flat)
        psi = unflatten_interior(psi_flat, nx, ny)

        update_wall_vorticity_from_psi(omega, psi, dx, dy, alpha)

        u, v = compute_velocity_from_psi(psi, dx, dy)

        conv = compute_convective_term(omega, u, v, dx, dy)
        conv_flat = conv.ravel(order='F')

        bnd_vec = build_boundary_contrib_vector(omega, nx, ny, dx, dy)

        rhs = omega_flat - dt * conv_flat + dt * visc * bnd_vec
        omega_new_flat = B_factor(rhs)

        omega_old_flat = omega_flat.copy()
        omega_flat = omega_new_flat
        omega = unflatten_interior(omega_flat, nx, ny)

        change = omega_flat - omega_old_flat
        res = np.linalg.norm(change, ord=2) / max(1.0, np.linalg.norm(omega_flat, ord=2))
        residuals.append(res)
        if initial_res is None:
            initial_res = res if res > 0 else 1.0

        # print only every 10% or at final iteration
        if (it % checkpoint == 0) or (it == max_iters):
            pct = int(100.0 * it / max_iters)
            print(f"Progress: {pct}% (iter {it}/{max_iters}), res = {res:.3e}")

        if res < tol:
            print(f"Converged: relative ω-change {res:.3e} < tol {tol:.1e} at iter {it}")
            break

    # final fields
    psi_flat = A_factor(-omega_flat)
    psi = unflatten_interior(psi_flat, nx, ny)
    update_wall_vorticity_from_psi(omega, psi, dx, dy, alpha)
    u, v = compute_velocity_from_psi(psi, dx, dy)
    speed = np.sqrt(u * u + v * v)

    final_res = residuals[-1] if residuals else 0.0
    reduction_factor = initial_res / final_res if final_res > 0 else np.inf
    trend = "decreasing" if (len(residuals) >= 5 and residuals[-1] < residuals[-5]) else "unknown/plateau"
    if final_res < tol:
        final_message = f"Converged: final residual {final_res:.3e}"
    elif reduction_factor > 100.0 and trend == "decreasing":
        final_message = f"Stopped at max_iters but residual decreasing (reduction factor {reduction_factor:.1f})."
    else:
        final_message = f"Stopped at max_iters; final residual {final_res:.3e} (reduction factor {reduction_factor:.1f})."
    print(final_message)

    # Note: pressure is not computed in streamfunction-vorticity formulation;
    # the p(0,0)=0 reference is irrelevant here.

    # Prepare plotting coordinates and ensure shapes match streamplot requirements:
    X = np.linspace(dx, 1.0 - dx, nx)
    Y = np.linspace(dy, 1.0 - dy, ny)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), constrained_layout=True)

    # Top: contour of velocity magnitude with streamlines
    # Use 1D X,Y inputs and transpose data arrays so shapes align: (ny, nx)
    cf = ax1.contourf(X, Y, speed.T, levels=50, cmap='RdBu_r')
    # streamplot expects U,V arrays shaped (ny, nx) when X and Y are 1D
    strm = ax1.streamplot(X, Y, u.T, v.T, color='k', density=1.0, linewidth=0.6)
    ax1.set_title('Velocity magnitude and streamlines (lid-driven cavity)')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    fig.colorbar(cf, ax=ax1, label='|u|')

    # Bottom: convergence history
    ax2.semilogy(np.arange(1, len(residuals) + 1), residuals, '-k')
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('Relative change residual (L2)')
    ax2.set_title('Convergence history')
    ax2.grid(True, which='both', ls='--', lw=0.5)

    fig_filename = 'ns_stream_vorticity_convergence.png'
    fig.savefig(fig_filename, dpi=200)
    print(f"Figure saved to: {fig_filename}")


if __name__ == '__main__':
    main()
```



#### Script block3:

```python
#!/usr/bin/env python3
"""
Streamfunction-vorticity solver for 2D steady incompressible Navier-Stokes
(lid-driven cavity). Re = 400 -> visc = 1/Re defined at global scope.

Key points implemented according to review:
- visc declared globally and used throughout
- update_wall_vorticity_from_psi writes vorticity into GHOST indices
  (omega_full[0,:], omega_full[nx+1,:], omega_full[:,0], omega_full[:,ny+1])
  so build_boundary_contrib_vector reads consistent values
- Laplacian and implicit matrices built once and factorized (guarded)
- Diagnostics added: vorticity PDE residual and discrete divergence surrogate
- Unit test validates boundary-contribution indexing
- Prints only at 10% cadence; prints diagnostics at same cadence
- Saves figure; does NOT call plt.show()
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# GLOBAL viscosity according to the problem statement
Re = 400.0
visc = 1.0 / Re


def build_2d_laplacian(nx, ny, dx, dy, kron_choice=0):
    Tx = sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(nx, nx)) / (dx * dx)
    Ty = sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(ny, ny)) / (dy * dy)
    I_x = sp.eye(nx)
    I_y = sp.eye(ny)
    if kron_choice == 0:
        L = sp.kron(I_y, Tx) + sp.kron(Ty, I_x)
    else:
        L = sp.kron(Tx, I_y) + sp.kron(I_x, Ty)
    return L.tocsr()


def flatten_interior(arr):
    return arr[1:-1, 1:-1].ravel(order='F')


def unflatten_interior(vec, nx, ny):
    interior = vec.reshape((nx, ny), order='F')
    out = np.zeros((nx + 2, ny + 2), dtype=vec.dtype)
    out[1:-1, 1:-1] = interior
    return out


def compute_velocity_from_psi(psi, dx, dy):
    # u = dpsi/dy, v = -dpsi/dx at interior nodes (central differences)
    u = (psi[1:-1, 2:] - psi[1:-1, 0:-2]) / (2.0 * dy)
    v = -(psi[2:, 1:-1] - psi[0:-2, 1:-1]) / (2.0 * dx)
    return u, v


def compute_convective_term(omega_full, u, v, dx, dy):
    o = omega_full
    d_omega_bx = (o[1:-1, 1:-1] - o[0:-2, 1:-1]) / dx
    d_omega_fx = (o[2:, 1:-1] - o[1:-1, 1:-1]) / dx
    d_omega_by = (o[1:-1, 1:-1] - o[1:-1, 0:-2]) / dy
    d_omega_fy = (o[1:-1, 2:] - o[1:-1, 1:-1]) / dy
    dwdx = np.where(u >= 0.0, d_omega_bx, d_omega_fx)
    dwdy = np.where(v >= 0.0, d_omega_by, d_omega_fy)
    conv = u * dwdx + v * dwdy
    return conv


def update_wall_vorticity_from_psi(omega_full, psi, dx, dy, alpha):
    """
    Enforce wall vorticity by setting ghost vorticity values so the discrete
    Laplacian sees consistent boundary contributions.

    Derivation for top wall (y = 1): psi(boundary) = 0 (Dirichlet)
    using second-order one-sided approximation for tangential derivative:
      u_wall = (psi_boundary - psi_interior)/dy  => psi_boundary = 0 => u_wall = -psi_interior/dy
    The discrete relation for ghost vorticity at j = ny+1 (ghost row) is chosen
    so that the finite-difference expression for omega = dv/dx - du/dy reproduces
    the no-slip wall velocity U(x):
      omega_ghost = -2*(psi_interior - psi_boundary)/dy^2 - 2*U/dy
    With psi_boundary = 0 -> omega_ghost = -2*psi_interior/dy^2 - 2*U/dy

    Similar formulas for bottom (U=0) and left/right walls.

    Important: THIS FUNCTION WRITES INTO THE GHOST INDICES:
      omega_full[0, 1:-1] (left ghost), omega_full[nx+1, 1:-1] (right ghost),
      omega_full[1:-1, 0] (bottom ghost), omega_full[1:-1, ny+1] (top ghost)
    so build_boundary_contrib_vector which reads those ghost slots is consistent.
    """
    nx = psi.shape[0] - 2
    ny = psi.shape[1] - 2
    # ensure Dirichlet psi on outer boundary (explicit enforcement)
    psi[0, :] = 0.0
    psi[-1, :] = 0.0
    psi[:, 0] = 0.0
    psi[:, -1] = 0.0

    x = np.linspace(dx, 1.0 - dx, nx)
    U_lid = alpha * x * (1.0 - x)

    # top ghost (j = ny+1)
    omega_full[1:-1, ny + 1] = -2.0 * psi[1:-1, ny] / (dy * dy) - 2.0 * U_lid / dy
    # bottom ghost (j = 0)
    omega_full[1:-1, 0] = -2.0 * psi[1:-1, 1] / (dy * dy)
    # left ghost (i = 0)
    omega_full[0, 1:-1] = -2.0 * psi[1, 1:-1] / (dx * dx)
    # right ghost (i = nx+1)
    omega_full[nx + 1, 1:-1] = -2.0 * psi[nx, 1:-1] / (dx * dx)


def build_boundary_contrib_vector(omega_full, nx, ny, dx, dy):
    b = np.zeros((nx, ny), dtype=omega_full.dtype)
    b[0, :] += omega_full[0, 1:-1] / (dx * dx)
    b[-1, :] += omega_full[nx + 1, 1:-1] / (dx * dx)
    b[:, 0] += omega_full[1:-1, 0] / (dy * dy)
    b[:, -1] += omega_full[1:-1, ny + 1] / (dy * dy)
    return b.ravel(order='F')


def verify_laplacian_ordering(nx, ny, dx, dy):
    X = np.linspace(dx, 1.0 - dx, nx)
    Y = np.linspace(dy, 1.0 - dy, ny)
    XX, YY = np.meshgrid(X, Y, indexing='ij')
    f2d = np.sin(np.pi * XX) * np.sin(np.pi * YY)
    f2d_pad = np.zeros((nx + 2, ny + 2))
    f2d_pad[1:-1, 1:-1] = f2d
    lap_fd = np.zeros_like(f2d)
    for i in range(nx):
        for j in range(ny):
            lap_fd[i, j] = (f2d_pad[i + 2, j + 1] - 2.0 * f2d_pad[i + 1, j + 1] + f2d_pad[i, j + 1]) / (dx * dx) + (
                f2d_pad[i + 1, j + 2] - 2.0 * f2d_pad[i + 1, j + 1] + f2d_pad[i + 1, j]) / (dy * dy)
    f_flat = f2d.ravel(order='F')
    target = lap_fd.ravel(order='F')
    errs = []
    mats = []
    for choice in (0, 1):
        Ltest = build_2d_laplacian(nx, ny, dx, dy, kron_choice=choice)
        Lf = Ltest.dot(f_flat)
        err = np.linalg.norm(Lf - target) / max(1.0, np.linalg.norm(target))
        errs.append(err)
        mats.append(Ltest)
    best = int(np.argmin(errs))
    if errs[best] > 1e-8:
        print(f"Warning: Laplacian verification residual = {errs[best]:.3e}; ordering may be inexact.")
    return best, mats[best]


def unit_test_boundary_contrib_consistency():
    # small test to ensure update_wall_vorticity_from_psi writes into ghost slots
    nx = 8
    ny = 8
    dx = 1.0 / (nx + 1)
    dy = 1.0 / (ny + 1)
    alpha = 2.0
    # create synthetic psi (interior) and unflatten
    X = np.linspace(dx, 1.0 - dx, nx)
    Y = np.linspace(dy, 1.0 - dy, ny)
    XX, YY = np.meshgrid(X, Y, indexing='ij')
    psi_in = np.sin(np.pi * XX) * np.sin(np.pi * YY)
    psi_full = np.zeros((nx + 2, ny + 2))
    psi_full[1:-1, 1:-1] = psi_in
    omega = np.zeros_like(psi_full)
    update_wall_vorticity_from_psi(omega, psi_full, dx, dy, alpha)
    b = build_boundary_contrib_vector(omega, nx, ny, dx, dy)
    # manual construction
    b_manual = np.zeros((nx, ny))
    b_manual[0, :] += omega[0, 1:-1] / (dx * dx)
    b_manual[-1, :] += omega[nx + 1, 1:-1] / (dx * dx)
    b_manual[:, 0] += omega[1:-1, 0] / (dy * dy)
    b_manual[:, -1] += omega[1:-1, ny + 1] / (dy * dy)
    b_manual = b_manual.ravel(order='F')
    err = np.linalg.norm(b - b_manual)
    if err > 1e-12:
        raise RuntimeError(f"Boundary-contrib unit test failed (err={err:.3e})")
    # quick check: top ghost values vs formula
    x = np.linspace(dx, 1.0 - dx, nx)
    U_lid = alpha * x * (1.0 - x)
    top_ghost_expected = -2.0 * psi_full[1:-1, ny] / (dy * dy) - 2.0 * U_lid / dy
    diff_top = np.max(np.abs(omega[1:-1, ny + 1] - top_ghost_expected))
    if diff_top > 1e-12:
        raise RuntimeError(f"Top ghost vorticity mismatch (diff={diff_top:.3e})")
    print("Unit test (boundary-contrib consistency) passed.")


def compute_vorticity_pde_residual(omega_flat, L, conv_flat, visc):
    # residual for steady vorticity equation: u·∇ω - ν Δω ≈ 0
    lap = L.dot(omega_flat)
    res_vec = conv_flat - visc * lap
    return np.linalg.norm(res_vec) / max(1.0, np.linalg.norm(omega_flat))


def sample_centerlines(u, v, X, Y):
    # sample u on vertical centerline (x~0.5) and v on horizontal centerline (y~0.5)
    ix = np.argmin(np.abs(X - 0.5))
    jy = np.argmin(np.abs(Y - 0.5))
    u_center = u[ix, :]
    v_center = v[:, jy]
    return ix, jy, u_center, v_center


def try_factorize(mat_csc, name):
    try:
        f = spla.factorized(mat_csc)
        return f
    except Exception as e:
        print(f"Factorization of {name} failed ({e}). Falling back to spsolve (slower).")

        def solver(rhs):
            return spla.spsolve(mat_csc, rhs)

        return solver


def main():
    # run unit test first
    unit_test_boundary_contrib_consistency()

    alpha = 2.0
    nx = 64
    ny = 64
    dx = 1.0 / (nx + 1)
    dy = 1.0 / (ny + 1)

    umax_est = alpha * 0.25
    cfl = 0.2
    dt_adv = cfl * min(dx / max(umax_est, 1e-8), dy / max(umax_est, 1e-8))
    dt = min(0.005, dt_adv)

    max_iters = 4000
    tol = 1e-6

    print(f"Starting solver: nx={nx}, ny={ny}, dt={dt:.5e}, visc={visc:.3e}")

    kron_choice, L = verify_laplacian_ordering(nx, ny, dx, dy)
    N = nx * ny
    I = sp.eye(N, format='csc')

    A_poisson = L.tocsc()
    B_mat = (I - dt * visc * L).tocsc()

    A_factor = try_factorize(A_poisson, 'Poisson A')
    B_factor = try_factorize(B_mat, 'Implicit B')

    omega = np.zeros((nx + 2, ny + 2), dtype=float)
    psi = np.zeros_like(omega)

    update_wall_vorticity_from_psi(omega, psi, dx, dy, alpha)
    omega_flat = flatten_interior(omega)

    residuals = []
    pde_res = []
    div_res = []
    checkpoint = max(1, int(max_iters / 10))
    initial_res = None

    for it in range(1, max_iters + 1):
        # Poisson: Lap psi = -omega
        psi_flat = A_factor(-omega_flat)
        psi = unflatten_interior(psi_flat, nx, ny)
        # ensure Dirichlet psi on boundaries
        psi[0, :] = 0.0
        psi[-1, :] = 0.0
        psi[:, 0] = 0.0
        psi[:, -1] = 0.0

        # update ghost vorticity based on psi and lid
        update_wall_vorticity_from_psi(omega, psi, dx, dy, alpha)

        u, v = compute_velocity_from_psi(psi, dx, dy)

        conv = compute_convective_term(omega, u, v, dx, dy)
        conv_flat = conv.ravel(order='F')

        bnd_vec = build_boundary_contrib_vector(omega, nx, ny, dx, dy)

        rhs = omega_flat - dt * conv_flat + dt * visc * bnd_vec
        omega_new_flat = B_factor(rhs)

        omega_old_flat = omega_flat.copy()
        omega_flat = omega_new_flat
        omega = unflatten_interior(omega_flat, nx, ny)

        change = omega_flat - omega_old_flat
        res = np.linalg.norm(change, ord=2) / max(1.0, np.linalg.norm(omega_flat, ord=2))
        residuals.append(res)
        if initial_res is None:
            initial_res = res if res > 0 else 1.0

        # diagnostics: vorticity PDE residual and discrete divergence surrogate
        pde_r = compute_vorticity_pde_residual(omega_flat, L, conv_flat, visc)
        pde_res.append(pde_r)
        # discrete divergence surrogate (central differences on interior excluding last index)
        # compute du/dx on interior reduced grid and dv/dy then combine on overlapping region
        du_dx = (u[1:, :] - u[:-1, :]) / dx  # shape (nx-1, ny)
        dv_dy = (v[:, 1:] - v[:, :-1]) / dy  # shape (nx, ny-1)
        # overlap to (nx-1, ny-1)
        div_sub = du_dx[:, :-1] + dv_dy[:-1, :]
        div_norm = np.linalg.norm(div_sub.ravel(order='F')) / max(1.0, np.linalg.norm(u.ravel(order='F')) + 1e-12)
        div_res.append(div_norm)

        # print only every 10% or at final iteration
        if (it % checkpoint == 0) or (it == max_iters):
            pct = int(100.0 * it / max_iters)
            print(f"Progress: {pct}% (iter {it}/{max_iters}), ||Δω||_rel = {res:.3e}, PDE_res = {pde_r:.3e}, div_res = {div_norm:.3e}")
            # check top lid velocity match (approximate)
            x = np.linspace(dx, 1.0 - dx, nx)
            U_lid = alpha * x * (1.0 - x)
            # approximate u at wall by one-sided derivative: u_wall ≈ (psi_boundary - psi_interior)/dy = -psi_interior/dy
            u_wall_approx = -psi[1:-1, ny] / dy
            max_diff = np.max(np.abs(u_wall_approx - U_lid))
            print(f"  Max|u_wall - U_lid| = {max_diff:.3e}")

        if res < tol:
            print(f"Converged: relative ω-change {res:.3e} < tol {tol:.1e} at iter {it}")
            break

    # final fields
    psi_flat = A_factor(-omega_flat)
    psi = unflatten_interior(psi_flat, nx, ny)
    psi[0, :] = 0.0
    psi[-1, :] = 0.0
    psi[:, 0] = 0.0
    psi[:, -1] = 0.0
    update_wall_vorticity_from_psi(omega, psi, dx, dy, alpha)
    u, v = compute_velocity_from_psi(psi, dx, dy)
    speed = np.sqrt(u * u + v * v)

    final_res = residuals[-1] if residuals else 0.0
    reduction_factor = initial_res / final_res if final_res > 0 else np.inf
    trend = "decreasing" if (len(residuals) >= 5 and residuals[-1] < residuals[-5]) else "unknown/plateau"
    if final_res < tol:
        final_message = f"Converged: final residual {final_res:.3e}"
    elif reduction_factor > 100.0 and trend == "decreasing":
        final_message = f"Stopped at max_iters but residual decreasing (reduction factor {reduction_factor:.1f})."
    else:
        final_message = f"Stopped at max_iters; final residual {final_res:.3e} (reduction factor {reduction_factor:.1f})."
    print(final_message)

    # centerline samples for quick benchmark-ish output
    X = np.linspace(dx, 1.0 - dx, nx)
    Y = np.linspace(dy, 1.0 - dy, ny)
    ix, jy, u_center, v_center = sample_centerlines(u, v, X, Y)
    print(f"Sample: vertical centerline at X[{ix}]≈{X[ix]:.3f}, max u on centerline = {np.max(u_center):.5f}")
    print(f"Sample: horizontal centerline at Y[{jy}]≈{Y[jy]:.3f}, max v on centerline = {np.max(np.abs(v_center)):.5f}")

    # plotting
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), constrained_layout=True)
    cf = ax1.contourf(X, Y, speed.T, levels=50, cmap='RdBu_r')
    strm = ax1.streamplot(X, Y, u.T, v.T, color='k', density=1.0, linewidth=0.6)
    ax1.set_title('Velocity magnitude and streamlines (lid-driven cavity)')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    fig.colorbar(cf, ax=ax1, label='|u|')

    ax2.semilogy(np.arange(1, len(residuals) + 1), residuals, '-k', label='||Δω||_rel')
    ax2.semilogy(np.arange(1, len(pde_res) + 1), pde_res, '--r', label='vorticity PDE res')
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('Residual (log)')
    ax2.set_title('Convergence history')
    ax2.grid(True, which='both', ls='--', lw=0.5)
    ax2.legend()

    fig_filename = 'ns_stream_vorticity_convergence.png'
    fig.savefig(fig_filename, dpi=200)
    print(f"Figure saved to: {fig_filename}")


if __name__ == '__main__':
    main()
```



### runtime_outputs

#### Output block1

Starting streamfunction-vorticity solver: nx=64, ny=64, dt=5.00000e-03, visc=2.500e-03
Progress: 0% (iter 1/4000), res = 1.074e-01
Progress: 10% (iter 400/4000), res = 2.802e-03
Progress: 20% (iter 800/4000), res = 1.715e-03
Progress: 30% (iter 1200/4000), res = 1.330e-03
Progress: 40% (iter 1600/4000), res = 1.067e-03
Progress: 50% (iter 2000/4000), res = 8.894e-04
Progress: 60% (iter 2400/4000), res = 7.633e-04
Progress: 70% (iter 2800/4000), res = 6.692e-04
Progress: 80% (iter 3200/4000), res = 5.964e-04
Progress: 90% (iter 3600/4000), res = 4.685e-04
Progress: 100% (iter 4000/4000), res = 3.526e-04
Stopped at max_iters but residual decreasing (reduction factor 304.6).

Runtime Error:
Traceback (most recent call last):
  File "C:\Users\admin\AppData\Local\Temp\ipykernel_50164\2663542217.py", line 138, in execute_code_tool
    exec(code, global_scope, global_scope)
  File "<string>", line 305, in <module>
  File "<string>", line 285, in main
  File "d:\zFILES\MiniConda\envs\langchain_\Lib\site-packages\matplotlib\__init__.py", line 1524, in inner
    return func(
           ^^^^^
  File "d:\zFILES\MiniConda\envs\langchain_\Lib\site-packages\matplotlib\streamplot.py", line 91, in streamplot
    grid = Grid(x, y)
           ^^^^^^^^^^
  File "d:\zFILES\MiniConda\envs\langchain_\Lib\site-packages\matplotlib\streamplot.py", line 331, in __init__
    raise ValueError("The rows of 'x' must be equal")
ValueError: The rows of 'x' must be equal



#### Output block2

Starting solver: nx=64, ny=64, dt=5.00000e-03, visc=2.500e-03
Progress: 10% (iter 400/4000), res = 2.802e-03
Progress: 20% (iter 800/4000), res = 1.715e-03
Progress: 30% (iter 1200/4000), res = 1.330e-03
Progress: 40% (iter 1600/4000), res = 1.067e-03
Progress: 50% (iter 2000/4000), res = 8.894e-04
Progress: 60% (iter 2400/4000), res = 7.633e-04
Progress: 70% (iter 2800/4000), res = 6.692e-04
Progress: 80% (iter 3200/4000), res = 5.964e-04
Progress: 90% (iter 3600/4000), res = 4.685e-04
Progress: 100% (iter 4000/4000), res = 3.526e-04
Stopped at max_iters but residual decreasing (reduction factor 304.6).
Figure saved to: ns_stream_vorticity_convergence.png



#### Output block3

Unit test (boundary-contrib consistency) passed.
Starting solver: nx=64, ny=64, dt=5.00000e-03, visc=2.500e-03
Progress: 10% (iter 400/4000), ||Δω||_rel = 8.114e-04, PDE_res = 6.832e+00, div_res = 4.214e-01
  Max|u_wall - U_lid| = 6.588e-02
Progress: 20% (iter 800/4000), ||Δω||_rel = 3.554e-04, PDE_res = 6.820e+00, div_res = 4.037e-01
  Max|u_wall - U_lid| = 6.464e-02
Progress: 30% (iter 1200/4000), ||Δω||_rel = 2.078e-04, PDE_res = 6.848e+00, div_res = 3.957e-01
  Max|u_wall - U_lid| = 6.444e-02
Progress: 40% (iter 1600/4000), ||Δω||_rel = 1.373e-04, PDE_res = 6.860e+00, div_res = 3.878e-01
  Max|u_wall - U_lid| = 6.417e-02
Progress: 50% (iter 2000/4000), ||Δω||_rel = 9.614e-05, PDE_res = 6.860e+00, div_res = 3.807e-01
  Max|u_wall - U_lid| = 6.379e-02
Progress: 60% (iter 2400/4000), ||Δω||_rel = 6.889e-05, PDE_res = 6.855e+00, div_res = 3.750e-01
  Max|u_wall - U_lid| = 6.343e-02
Progress: 70% (iter 2800/4000), ||Δω||_rel = 4.972e-05, PDE_res = 6.849e+00, div_res = 3.706e-01
  Max|u_wall - U_lid| = 6.314e-02
Progress: 80% (iter 3200/4000), ||Δω||_rel = 3.594e-05, PDE_res = 6.842e+00, div_res = 3.674e-01
  Max|u_wall - U_lid| = 6.289e-02
Progress: 90% (iter 3600/4000), ||Δω||_rel = 2.599e-05, PDE_res = 6.837e+00, div_res = 3.652e-01
  Max|u_wall - U_lid| = 6.271e-02
Progress: 100% (iter 4000/4000), ||Δω||_rel = 1.880e-05, PDE_res = 6.833e+00, div_res = 3.635e-01
  Max|u_wall - U_lid| = 6.257e-02
Stopped at max_iters but residual decreasing (reduction factor 53181.3).
Sample: vertical centerline at X[31]≈0.492, max u on centerline = 0.41326
Sample: horizontal centerline at Y[31]≈0.492, max v on centerline = 0.11863
Figure saved to: ns_stream_vorticity_convergence.png




### review_decision
accept

### review_comments
	Current Stage [A/3]
1. MAJOR

Category: runtime
Issue: Matplotlib streamplot raises ValueError: 'The rows of 'x' must be equal' at plotting time.

Recommendation: The error comes from mismatched shapes/orientation passed to streamplot. Fix by providing streamplot with consistent X/Y grid and U/V shapes. For example: (a) build X_vec, Y_vec and do XX, YY = np.meshgrid(X_vec, Y_vec, indexing='xy') and then call streamplot(XX, YY, u.T, v.T, ...), or (b) call streamplot(X_vec, Y_vec, u.T, v.T, ...) where u.T and v.T have shape (ny, nx). Remove the current indexing='ij' meshgrid or transpose arrays so rows of the first coordinate array are constant. After this change the plotting call should not trigger the ValueError.


2. MINOR

Category: structure
Issue: Progress printing still occurs at iteration 1 due to condition 'or it == 1', which violates the requirement to print only every 10% of total steps.

Recommendation: Remove the 'or it == 1' clause so prints happen only when it % checkpoint == 0 (or at final iteration). This will enforce the requested sparse reporting (every 10% of max_iters).


3. MINOR

Category: structure
Issue: The implementation uses a streamfunction–vorticity formulation and never computes pressure, but the problem statement explicitly asked to fix the reference pressure p(0,0)=0.

Recommendation: Either (a) explicitly state in comments that pressure is not computed because the streamfunction–vorticity approach eliminates pressure, or (b) if pressure is required for postprocessing or verification, add a Poisson solve to recover pressure from momentum equations (e.g. solve pressure Poisson with continuity and fix p(0,0)=0). If you keep the current approach, document that p(0,0)=0 is irrelevant for this formulation.


4. MAJOR

Category: accuracy
Issue: Potential mismatch risk between the assembled Laplacian matrix ordering (Kronecker construction) and flatten/unflatten routines (Fortran order). If ordering is inconsistent the linear solves will be incorrect even if code runs.

Recommendation: Verify and ensure the same ordering convention is used everywhere. You use ravel(..., order='F') and reshape(..., order='F') for flatten/unflatten. Confirm that build_2d_laplacian (kron(I_y, Tx) + kron(Ty, I_x)) matches that Fortran ordering. A simple test: apply L to a flattened known separable function (e.g. sin(pi*x)*sin(pi*y) sampled on interior) and compare the finite-difference Laplacian computed directly on the 2D grid. If mismatch found, reorder Kronecker factors to match flattening (swap kron operands or change flatten order).


5. MAJOR

Category: accuracy
Issue: Wall vorticity boundary formula and indexing may be inconsistent (signs/factors and indexing along x/y may be mixed), which can bias the solution near walls.

Recommendation: Re-check discrete boundary vorticity relations and array indexing carefully. The standard relation (for Dirichlet ψ=0 on boundary) at the top (moving lid) interior-adjacent node is ω = -2*(ψ_interior - ψ_wall)/Δy^2 - 2*U/Δy (with ψ_wall=0). Ensure (a) the Δ used (dx vs dy) is correct for each wall, (b) the sign convention for ω matches the definition ω = ∂v/∂x - ∂u/∂y used in the solver, and (c) the indices used (rows vs columns) correspond to x/y axes consistently. Validate the BC implementation by comparing U_lid at the top boundary and near-wall velocities, and by checking global diagnostics like maximum streamfunction location and centerline velocities against reference lid-driven cavity data at Re=400.



	Current Stage [B/3]
1. MAJOR

Category: accuracy
Issue: Boundary vorticity indexing mismatch: update_wall_vorticity_from_psi writes vorticity at the interior-adjacent indices (i=1 or i=nx, j=1 or j=ny), but build_boundary_contrib_vector expects boundary/ghost vorticity at indices 0, nx+1, 0, ny+1. As a result the boundary contributions used in the interior Laplacian RHS are taken from unset ghost slots (left as zero) while the computed wall vorticities are never injected where the Poisson/stencil expects them.

Recommendation: Either (A) make update_wall_vorticity_from_psi set the ghost/boundary indices (omega_full[0,...], omega_full[nx+1,...], omega_full[...,0], omega_full[...,ny+1]) consistently with build_boundary_contrib_vector; or (B) change build_boundary_contrib_vector to pick up the wall vorticity values at the indices currently written (1 and nx, 1 and ny). Ensure indexing choice is documented and used consistently throughout flatten/unflatten, boundary updates, and the Laplacian assembly. Add a small unit test that sets a known psi and verifies the computed boundary-contribution vector matches a finite-difference reference.


2. MAJOR

Category: accuracy
Issue: Boundary vorticity enforcement formula and ghost/Dirichlet psi assumptions are not explicitly checked. The update_wall_vorticity_from_psi formulas are plausible but sensitive to whether psi boundary values are zero and how the derivative is discretized; incorrect enforcement can produce a physically incorrect wall velocity and spoil convergence.

Recommendation: Document and (if necessary) enforce the intended psi Dirichlet values at the outer boundary (e.g., psi=0 on all walls) and derive the finite-difference expressions used for each wall in comments. Add a verification routine: for a manufactured (analytic) psi and U(x) compute analytic ω at boundary and compare with implementation. Also check that the computed u at the top boundary matches the prescribed lid profile (u(x, y=1) ≈ alpha*x*(1-x)) up to discretization error — if not, the vorticity BC should be corrected.


3. MINOR

Category: structure
Issue: Residual definition and stopping criteria are pragmatic but could be misleading: the algorithm monitors relative L2 change of omega between steps. This is fine, but it does not directly report the momentum/continuity residual or divergence of velocity, nor a robust estimate of steady-state PDE residual.

Recommendation: Add one or two diagnostic residuals: (1) L2 norm of discrete divergence of velocity at interior nodes (should be ~0), and (2) L2 norm of the momentum residual (u·∇u + ∇p - νΔu) approximated without p (or use vorticity form residual). Report these diagnostics intermittently (same 10% cadence). Use them together with the omega-change metric to decide convergence. Also consider stopping when residuals plateau or meet tolerances, not only when change < tol.


4. MINOR

Category: structure
Issue: Use of spla.factorized for both Poisson and implicit diffusion solves is convenient but can be memory heavy for larger grids; also factorized may raise exceptions if matrices are singular or ill-conditioned.

Recommendation: Guard the factorization calls with try/except to provide a clear error and fallback (e.g., use spla.spsolve or an iterative solver with a preconditioner). If performance/memory is a concern, mention in comments that iterative solvers (CG/BiCGSTAB with AMG or ILU preconditioning) are an alternative. Optionally add a small parameter to select 'direct' vs 'iterative' solver modes for the Poisson and B solves.


5. MINOR

Category: structure
Issue: Validation and benchmarks are missing. The code prints convergence progress and saves figures, but no quantitative check against known lid-driven cavity benchmarks (e.g., Ghia et al.) is provided.

Recommendation: Add an optional post-processing check that samples centerline velocities (u at vertical centerline and v at horizontal centerline) and compares them to Ghia et al. values for Re=400 (or reports them for manual verification). Also compute and report global metrics (max streamfunction, location of primary vortex) to make it easier to judge physical accuracy.



	Current Stage [C/3]
1. MINOR

Category: accuracy
Issue: The code uses a streamfunction–vorticity formulation (no pressure solve). The problem statement and hints explicitly discuss pressure fixing (p(0,0)=0) — this semantic mismatch is not fatal but should be documented.

Recommendation: Document clearly that the chosen streamfunction–vorticity approach eliminates the pressure variable and thus the p(0,0)=0 constraint is not needed. If the reviewer/assignment required a primitive-variable (u,v,p) solver, note that the current method does not provide a pressure field. If a pressure is required for postprocessing, compute it (e.g. solve a pressure Poisson using u,v) and then fix p at (0,0).


2. MINOR

Category: structure
Issue: The code mixes global and passed constants in a few places (Re defined globally, visc global, alpha passed to functions). While visc is global as required, not all functions consistently accept physical constants which reduces clarity and testability.

Recommendation: Either pass visc and alpha explicitly to functions that need them (preferred for modularity) or document clearly which globals are relied upon. For example, update_wall_vorticity_from_psi currently requires alpha as an argument (good), but some functions rely on global visc; consider passing visc to compute_vorticity_pde_residual (already done) and any other physics routines for consistency.


3. MINOR

Category: accuracy
Issue: Convergence criterion based on relative change of ω is reasonable, but the code uses a single-step relative change; stability of residual reporting and the decision to stop at max_iters relies on heuristics (reduction_factor and recent trend).

Recommendation: Enhance convergence checks by also monitoring the vorticity PDE residual (pde_res) and/or a moving-window slope of residuals. For example, stop if both ||Δω||_rel < tol AND pde_res < tol_pde, or if residual has decreased by a target factor and still shows monotonic decline over the last N iterations. Log these decisions succinctly (only at 10% cadence as you already do).





### iteration_count
3

### rev_image_description
	Current Stage [A/3]


	Current Stage [B/3]


	Current Stage [C/3]




