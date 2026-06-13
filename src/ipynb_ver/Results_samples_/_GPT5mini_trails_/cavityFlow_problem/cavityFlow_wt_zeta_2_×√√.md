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
Never use the keyword argument 'tol' in any SciPy solver. Always use 'atol'. 
Print concise progress information ONLY every 10% of total steps.
Focus on generating runnable code with reasonable physical results, do not over-optimize.
Do NOT conclude [not converged] solely because max_iters is reached, assess convergence based on whether residuals are still decreasing. 


### expanded_prob
We must compute a steady 2-D incompressible Navier–Stokes flow on the unit square Omega = [0,1]^2 with Reynolds number Re = 400 (visc = 1.0/Re). The PDE system is nonlinear (convective terms u·grad u) and constrained (divergence-free velocity). Boundary conditions: a moving lid on the top (u = alpha x(1-x), v = 0, alpha = 2) and no-slip elsewhere; pressure can be fixed by p(0,0)=0. The numerical task is to produce a stable and efficient solver that converges to the steady solution and produces useful diagnostics and plots: a contour of velocity magnitude with the RdBu_r colormap overlaid by streamlines, and a convergence-history plot; both saved to files (no plt.show()).

Primary mathematical and numerical challenges:
- Nonlinearity: convective terms create multiple difficulties — they require linearization (Picard or Newton) or explicit/semi-implicit time stepping to converge to a steady state. At Re=400 convective transport is significant and may produce boundary-layer-like gradients and corner singularities.
- Incompressibility (saddle-point structure): coupled velocity-pressure equations form a saddle-point system that requires either careful discretization (inf-sup stable elements or a staggered grid) or stabilization (Rhie–Chow, SUPG, pressure-stabilization) on collocated meshes.
- Stability vs. accuracy tradeoffs: central differences for convection yield high accuracy but can produce spurious oscillations at higher Re; upwinding or flux-limited schemes increase stability at some accuracy cost.
- Linear algebraic solvers and preconditioning: discretization produces large sparse linear systems; obtaining acceptable runtimes requires iterative solvers with good preconditioners (multigrid, ILU, block-Schur approximations). Memory and CPU costs depend strongly on solver choices.
- Convergence diagnostics: for a steady solver we must monitor residual norms and divergence; do not prematurely declare failure when iteration counts reach max_iters — check trend of residuals instead.
- Implementation constraints: define visc = 1.0/Re at global scope and use that identifier everywhere; pass all constants explicitly into functions; when using SciPy iterative solvers use atol (not tol); print progress only every 10% of total steps. Plotting must save figures and not call plt.show().

### solution_plans
	Current Stage [A/2]
solu_name='Fractional-step (projection) method on a collocated grid with semi-implicit time-marching and multigrid-accelerated Poisson solves' content="Governing idea:\n- March pseudo-time to steady state with a fractional-step (projection) method: solve momentum for an intermediate velocity (with implicit viscous term and linearized convection), solve a Poisson equation for pressure correction to enforce divergence-free constraint, correct velocity, iterate until steady.\n\nAlgorithmic steps (numbered):\n1) Discretize Omega on a uniform collocated grid (or slightly staggered if preferred). Choose second-order central differences for diffusion and a stabilized convective discretization: either QUICK/second-order upwind or central + small artificial dissipation. Use Rhie–Chow interpolation if collocated.\n2) Set global visc = 1.0/Re. Initialize u, v, p (e.g., zero velocity and p=0), apply BCs (moving lid top, no-slip elsewhere). Enforce p(0,0)=0 at each pressure solve to fix constant.\n3) Use pseudo-time stepping (pseudo-t) with a semi-implicit scheme: treat viscous terms implicitly to avoid stiff timestep constraints and treat convective terms explicitly (or linearize convective term with previous iterate for Picard-type treatment). Use backward-Euler (or BDF2 for second-order) in pseudo-time.\n4) At each pseudo-time step k:\n   a) Assemble/linearize the momentum operator A_k = M/dt + visc * Laplacian + linearized convection terms. Right-hand side = M/dt * u^n + explicit convective contributions + body terms (if any).\n   b) Solve the two decoupled momentum equations for intermediate velocity u* and v* (use sparse GMRES with ILU or an AMG-accelerated solver). Use SciPy iterative solvers with atol set and appropriate restart.\n   c) Form and solve the pressure-Poisson equation: Laplacian(phi) = div(u*)/dt with Neumann/Dirichlet consistent BCs; enforce phi(0,0) = 0 to fix reference pressure. Solve the Poisson using a multigrid solver (PyAMG) or SciPy CG preconditioned by AMG/ILU. Reuse preconditioner between steps.\n   d) Correct velocity: u^{n+1} = u* - dt * grad(phi). Update p^{n+1} = p^n + phi.\n   e) Apply under-relaxation to velocity/pressure updates if oscillatory (u := omega_u * u^{n+1} + (1-omega_u)*u^n). Use omega in (0.6,1.0).\n5) Time stepping acceleration:\n   - Use adaptive pseudo-time step dt based on CFL (max |u| * dt / dx <= CFL_target) and increase dt as residuals reduce; cap dt to avoid instability.\n   - Stop when L2 norm of momentum residual and max divergence both drop below user tolerances or show stagnation. Print progress every 10% of total planned iterations (or every 10% of elapsed pseudo-time if adaptively stepping).\n6) Post-process: compute velocity magnitude, generate contour with 'RdBu_r' colormap and overlay streamlines sampled from the steady velocity field. Also create a subplot with convergence history (residual vs iteration). Save figure(s) to disk — do not call plt.show().\n\nStability / accuracy / complexity / efficiency limitations:\n- Stability depends on convective discretization and pseudo-time step control. Explicit treatment of convection limits dt by CFL; linearization (Picard) relaxes that but slows convergence. Using fully implicit convection (Newton linearization) increases per-step cost.\n- Collocated grids require Rhie–Chow to avoid checkerboarding; staggered grids alleviate this at some implementation cost.\n- Multigrid/AMG for the Poisson is critical for efficiency; if unavailable, direct sparse solvers may be necessary but will be slower and memory intensive for fine grids.\n- Method is robust and simple to implement, but convergence to steady can be slow for higher Re unless under-relaxation, adaptive dt, and good preconditioners are used. Memory/CPU cost dominated by repeated momentum and Poisson solves; reuse of preconditioners and matrix assembly amortizes cost.\n- Ensure all constants (visc, dt, omega, grid params, max_iters, atol) are passed explicitly to functions; when calling SciPy solvers use atol (never tol). Print progress only at the 10% cadence required."

	Current Stage [B/2]
solu_name='Monolithic steady Newton–Krylov with block-preconditioned GMRES using a finite-element discretization (Taylor–Hood) and continuation in Re' content='Governing idea:\n- Formulate the steady nonlinear residual R(U) for unknown U = [u,v,p]. Solve R(U)=0 via inexact Newton iterations. Each Newton step requires solving a linear saddle-point system J * delta = -R, where J is the Jacobian. Use GMRES with a block preconditioner approximating the Schur complement (block triangular preconditioner using AMG for the velocity block and a mass/Laplacian approximation for the Schur complement). Use continuation in Re (start at low Re and ramp to 400) to enhance robustness.\n\nAlgorithmic steps (numbered):\n1) Use an inf-sup stable finite-element pair (Taylor–Hood: P2 velocity, P1 pressure) and assemble system matrices: stiffness (viscous Laplacian), convection linearization contributions, divergence B and gradient operators, and mass matrices. Define visc = 1.0/Re globally and pass to assembly routines.\n2) Initialize U^0 (e.g., Stokes solution or low-Re steady solution). Optionally perform continuation: start at Re_small (e.g., 10), solve to steady, then increase Re in stages until 400, using previous solution as initial guess.\n3) Newton loop until residual norm is below tolerance or stagnates:\n   a) Form residual R(U^k) and assemble Jacobian J(U^k) — analytic Jacobian of convection and diffusion terms (or approximate Jacobian using Picard linearization for lower cost).\n   b) Solve the linearized system J * delta = -R using GMRES with a block preconditioner P ≈ [[A, 0],[B, -S]] where A is the velocity block and S ≈ B A^{-1} B^T the Schur complement. Implement P as a block-triangular preconditioner using:\n       - AMG (or ILU) for approximate inverse of A (velocity block). Reuse AMG hierarchy across Newton steps where possible.\n       - Approximate S by visc-scaled pressure mass matrix or by B * (diag(A)^{-1}) * B^T; invert S approximately with AMG or a few CG iterations.\n   c) Use an inexact Newton stopping criterion for inner GMRES (tighter solves as Newton progresses). Use atol in SciPy solvers (never tol).\n   d) Apply a globalization strategy (Armijo/backtracking) on the Newton update if norm does not reduce.\n4) After convergence, compute velocity magnitude and produce the required plot (RdBu_r contour + streamlines) and convergence history. Save figures; do not call plt.show().\n5) Logging: print concise progress only every 10% of total Newton steps (or during continuation stages). Record residual norms for plotting.\n\nStability / accuracy / complexity / efficiency limitations:\n- Newton–Krylov is fast in iteration count (quadratic near solution) but each iteration is expensive because it solves a large linear system; cost dominated by GMRES and preconditioning.\n- Building and applying an effective block preconditioner requires access to robust AMG/ILU libraries and careful tuning; implementation complexity is higher than fractional-step schemes.\n- Memory demands are significant (matrix assembly, AMG hierarchies). For very fine meshes direct application without good preconditioning becomes intractable.\n- Continuation in Re mitigates nonlinearity but adds extra solves as Re is ramped up — overall extra cost but improves robustness.\n- Requires a finite-element library (assembly of consistent matrices) and linear-algebra backends (PETSc, hypre, PyAMG) for practical efficiency; pure SciPy may be usable for moderate grids but will scale worse.\n- As always, pass visc and other constants explicitly to assembly/solver functions and use atol in iterative solver calls. Print progress at the specified 10% cadence.'



### technical_spec
	Current Stage [A/3]
We implement a semi-implicit fractional-step (projection) solver on a uniform collocated finite-difference grid for the 2D steady incompressible Navier–Stokes lid-driven cavity. The algorithm marches pseudo-time to steady state: convection is treated explicitly, viscous terms implicitly (Helmholtz solve reused via LU factorization built once), and incompressibility enforced by solving a pressure-Poisson equation with Neumann discretization (also assembled and factorized once). Boundary velocities are Dirichlet (moving lid on top), and pressure reference is enforced at one corner. All sparse matrices are assembled once before the time-loop. Convergence is monitored via velocity-change norms and maximum divergence, printed only every 10% of the maximum iterations. Final outputs: saved figure with velocity magnitude contour (RdBu_r) overlaid by streamlines and convergence history, plus printed quantitative metrics.

	Current Stage [B/3]
- Build discrete Laplacians once (interior Dirichlet Laplacian and full-grid Neumann Laplacian) with correct sign convention: center entries = -sum(off-diagonals), off-diagonals = +1/dx^2. Assemble Helmholtz A = I - dt*visc*L_int once and LU-factorize.
- Precompute boundary-contribution vectors (bc_vec_u, bc_vec_v) once; these store the NEGATIVE of neighbor coefficients times known boundary values so RHS can add dt*visc*bc_vec directly.
- Time-march pseudo-time to steady state: explicit convective term (central differences), implicit viscous via Helmholtz solve, and pressure projection via Poisson with Neumann BC (reference pressure fixed at a corner). Pressure Laplacian is factorized once (with reference row replacement) and used to solve for divergence correction each step.
- Stability: choose dt conservatively using CFL based on lid velocity; under-relax velocities each step with omega to improve robustness. Check for NaNs and abort with diagnostics if they occur.
- Outputs: save a single figure with (left) contour of velocity magnitude (RdBu_r) plus streamlines and (right) convergence histories; print concise iteration progress only every 10% of max_iters and final quantitative metrics (velocity-change norm, max divergence, mean |u|).

	Current Stage [C/3]
This script solves the steady 2D incompressible Navier–Stokes lid-driven cavity using a pseudo-time projection method. Key steps: central-difference convective terms (explicit), implicit viscous solve via a pre-factorized Helmholtz matrix built once for interior unknowns, pressure projection using a pre-factorized Neumann Poisson operator (with one reference pressure row replaced), and velocity correction. Boundary contributions for the implicit viscous RHS are precomputed once. Convergence diagnostics focus on interior divergence and interior velocity-change norms (L2 and Linf). The script reports concise progress every 10% of iterations and verifies that the Poisson solve reduces the interior divergence by checking the Poisson residual. The final figure (saved) shows the velocity magnitude with streamlines and the convergence histories (interior norms).



### python_codes

#### Script block1:

```python
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# Global viscosity from Reynolds number
Re = 400.0
visc = 1.0 / Re  # MUST use 'visc' everywhere


def build_interior_laplacian_dirichlet(Nx, Ny, dx):
    """
    Build 2D Laplacian matrix for interior nodes only (Dirichlet BC on boundaries).
    Returns A (sparse) and bc_coeffs mapping to boundary values for each interior index.
    Interior nodes are those with 1 <= i <= Nx-2, 1 <= j <= Ny-2 (0-based indexing)
    Flattening order: row-major over j (y) then i (x).
    bc_coeffs_u/v should be applied as sum(c * boundary_value_at_neighbor) for each interior node.
    We'll return a single bc_coefficient mapping (same stencil) and user can multiply by boundary arrays when needed.
    """
    nx_in = Nx - 2
    ny_in = Ny - 2
    N = nx_in * ny_in
    dx2 = dx * dx

    rows = []
    cols = []
    data = []
    bc_coeffs = np.zeros(N)  # coefficients to multiply boundary values for RHS

    def idx(i, j):
        return j * nx_in + i

    for j in range(ny_in):
        for i in range(nx_in):
            row = idx(i, j)
            center = 0.0
            # neighbors: left (i-1,j), right (i+1,j), down (i,j-1), up (i,j+1)
            # map to global indices: global_i = i+1, global_j = j+1
            gi = i + 1
            gj = j + 1

            # left neighbor
            if gi - 1 >= 1:
                # interior neighbor
                rows.append(row); cols.append(idx(i - 1, j)); data.append(1.0 / dx2)
                center -= 1.0 / dx2
            else:
                # neighbor is boundary -> Dirichlet: contribution moves to RHS via bc_coeffs
                bc_coeffs[row] += 1.0 / dx2
                center -= 1.0 / dx2

            # right neighbor
            if gi + 1 <= Nx - 2 + 0:
                if gi + 1 <= Nx - 2:
                    # interior neighbor exists
                    rows.append(row); cols.append(idx(i + 1, j)); data.append(1.0 / dx2)
                    center -= 1.0 / dx2
                else:
                    # this condition won't happen because gi+1 <= Nx-2 check suffices
                    pass
            else:
                # right boundary
                bc_coeffs[row] += 1.0 / dx2
                center -= 1.0 / dx2

            # down neighbor (j-1)
            if gj - 1 >= 1:
                rows.append(row); cols.append(idx(i, j - 1)); data.append(1.0 / dx2)
                center -= 1.0 / dx2
            else:
                bc_coeffs[row] += 1.0 / dx2
                center -= 1.0 / dx2

            # up neighbor (j+1)
            if gj + 1 <= Ny - 2:
                rows.append(row); cols.append(idx(i, j + 1)); data.append(1.0 / dx2)
                center -= 1.0 / dx2
            else:
                bc_coeffs[row] += 1.0 / dx2
                center -= 1.0 / dx2

            # center
            rows.append(row); cols.append(row); data.append(-center)

    A = sp.csr_matrix((data, (rows, cols)), shape=(N, N), dtype=np.float64)
    return A.tocsc(), bc_coeffs


def build_pressure_laplacian_neumann(Nx, Ny, dx):
    """
    Build 2D Laplacian matrix representing Neumann BCs via mirrored ghost nodes (zero normal derivative).
    All nodes included; flattening order row-major over j (y) then i (x).
    We will later fix p at reference node by replacing the corresponding row.
    """
    N = Nx * Ny
    dx2 = dx * dx

    rows = []
    cols = []
    data = []

    def idx(i, j):
        return j * Nx + i

    for j in range(Ny):
        for i in range(Nx):
            row = idx(i, j)
            center = 0.0
            # x-direction
            left = i - 1
            right = i + 1
            if left >= 0 and right <= Nx - 1:
                # interior in x-direction
                rows.append(row); cols.append(idx(left, j)); data.append(1.0 / dx2)
                rows.append(row); cols.append(idx(right, j)); data.append(1.0 / dx2)
                center -= 2.0 / dx2
            elif left < 0 and right <= Nx - 1:
                # left missing (boundary), apply Neumann: ghost = right neighbor -> coefficient 2/right
                rows.append(row); cols.append(idx(right, j)); data.append(2.0 / dx2)
                center -= 2.0 / dx2
            elif right > Nx - 1 and left >= 0:
                rows.append(row); cols.append(idx(left, j)); data.append(2.0 / dx2)
                center -= 2.0 / dx2
            # else both missing won't happen unless Nx==1

            # y-direction
            down = j - 1
            up = j + 1
            if down >= 0 and up <= Ny - 1:
                rows.append(row); cols.append(idx(i, down)); data.append(1.0 / dx2)
                rows.append(row); cols.append(idx(i, up)); data.append(1.0 / dx2)
                center -= 2.0 / dx2
            elif down < 0 and up <= Ny - 1:
                rows.append(row); cols.append(idx(i, up)); data.append(2.0 / dx2)
                center -= 2.0 / dx2
            elif up > Ny - 1 and down >= 0:
                rows.append(row); cols.append(idx(i, down)); data.append(2.0 / dx2)
                center -= 2.0 / dx2

            rows.append(row); cols.append(row); data.append(-center)

    Lp = sp.csr_matrix((data, (rows, cols)), shape=(N, N), dtype=np.float64)
    return Lp.tocsc()


def build_helmholtz_matrix_for_vel(Nx, Ny, dx, dt, visc):
    """
    Using the interior Laplacian for Dirichlet BCs, build Helmholtz matrix A = I - dt*visc*L
    Returns A (csc) and bc_coeffs (to apply Dirichlet boundary values into RHS).
    """
    L_int, bc_coeffs = build_interior_laplacian_dirichlet(Nx, Ny, dx)
    A = sp.eye(L_int.shape[0], format='csc') - (dt * visc) * L_int
    return A, bc_coeffs


def idx_interior(i, j, Nx):
    # interior indexing function: interior grid indices flattened
    nx_in = Nx - 2
    return (j - 1) * nx_in + (i - 1)


def flatten_interior(field):
    # field is full Nx x Ny array; extract interior nodes into vector
    Ny, Nx = field.shape
    nx_in = Nx - 2
    ny_in = Ny - 2
    vec = np.empty(nx_in * ny_in, dtype=np.float64)
    k = 0
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            vec[k] = field[j, i]
            k += 1
    return vec


def unflatten_interior(vec, Nx, Ny):
    nx_in = Nx - 2
    ny_in = Ny - 2
    field = np.zeros((Ny, Nx), dtype=np.float64)
    k = 0
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            field[j, i] = vec[k]
            k += 1
    return field


def compute_convective(u, v, dx):
    """
    Compute convective terms (u·grad u) and (u·grad v) at interior nodes using central differences.
    u,v are full arrays (Ny,Nx) including boundary values.
    Returns conv_u and conv_v at interior nodes flattened.
    """
    Ny, Nx = u.shape
    nx_in = Nx - 2
    ny_in = Ny - 2
    conv_u = np.empty(nx_in * ny_in, dtype=np.float64)
    conv_v = np.empty(nx_in * ny_in, dtype=np.float64)
    dx2 = 2.0 * dx
    k = 0
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            du_dx = (u[j, i + 1] - u[j, i - 1]) / dx2
            du_dy = (u[j + 1, i] - u[j - 1, i]) / dx2
            dv_dx = (v[j, i + 1] - v[j, i - 1]) / dx2
            dv_dy = (v[j + 1, i] - v[j - 1, i]) / dx2
            ui = u[j, i]
            vi = v[j, i]
            conv_u[k] = ui * du_dx + vi * du_dy
            conv_v[k] = ui * dv_dx + vi * dv_dy
            k += 1
    return conv_u, conv_v


def gradient_flattened_scalar(phi, dx):
    """
    Compute gradient of scalar phi (full array size Ny x Nx) and return dp/dx and dp/dy at interior nodes flattened.
    """
    Ny, Nx = phi.shape
    nx_in = Nx - 2
    ny_in = Ny - 2
    dphidx = np.empty(nx_in * ny_in, dtype=np.float64)
    dphidy = np.empty(nx_in * ny_in, dtype=np.float64)
    dx2 = 2.0 * dx
    k = 0
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            dphidx[k] = (phi[j, i + 1] - phi[j, i - 1]) / dx2
            dphidy[k] = (phi[j + 1, i] - phi[j - 1, i]) / dx2
            k += 1
    return dphidx, dphidy


def divergence_full(u, v, dx):
    """
    Compute divergence at all nodes using central differences for interior and one-sided for boundaries.
    Returns div array shape (Ny, Nx).
    """
    Ny, Nx = u.shape
    div = np.zeros_like(u)
    # interior
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            dudx = (u[j, i + 1] - u[j, i - 1]) / (2.0 * dx)
            dvdy = (v[j + 1, i] - v[j - 1, i]) / (2.0 * dx)
            div[j, i] = dudx + dvdy
    # boundaries: use one-sided differences
    # left and right
    for j in range(Ny):
        # left i=0
        dudx = (u[j, 1] - u[j, 0]) / dx
        # dvdy using neighbors if possible
        if j == 0:
            dvdy = (v[1, 0] - v[0, 0]) / dx
        elif j == Ny - 1:
            dvdy = (v[-1, 0] - v[-2, 0]) / dx
        else:
            dvdy = (v[j + 1, 0] - v[j - 1, 0]) / (2.0 * dx)
        div[j, 0] = dudx + dvdy
        # right i=Nx-1
        dudx = (u[j, -1] - u[j, -2]) / dx
        if j == 0:
            dvdy = (v[1, -1] - v[0, -1]) / dx
        elif j == Ny - 1:
            dvdy = (v[-1, -1] - v[-2, -1]) / dx
        else:
            dvdy = (v[j + 1, -1] - v[j - 1, -1]) / (2.0 * dx)
        div[j, -1] = dudx + dvdy
    # top and bottom
    for i in range(Nx):
        # bottom j=0
        if i == 0:
            dudx = (u[0, 1] - u[0, 0]) / dx
        elif i == Nx - 1:
            dudx = (u[0, -1] - u[0, -2]) / dx
        else:
            dudx = (u[0, i + 1] - u[0, i - 1]) / (2.0 * dx)
        dvdy = (v[1, i] - v[0, i]) / dx
        div[0, i] = dudx + dvdy
        # top j=Ny-1
        if i == 0:
            dudx = (u[-1, 1] - u[-1, 0]) / dx
        elif i == Nx - 1:
            dudx = (u[-1, -1] - u[-1, -2]) / dx
        else:
            dudx = (u[-1, i + 1] - u[-1, i - 1]) / (2.0 * dx)
        dvdy = (v[-1, i] - v[-2, i]) / dx
        div[-1, i] = dudx + dvdy
    return div


def apply_velocity_bcs(u, v, alpha, Nx, Ny, x_coords):
    # Apply Dirichlet BCs: top boundary (moving lid) and no-slip elsewhere
    # u, v are arrays shape (Ny, Nx)
    # top boundary j = Ny-1
    jtop = Ny - 1
    u[jtop, :] = alpha * x_coords * (1.0 - x_coords)
    v[jtop, :] = 0.0
    # other boundaries set to zero
    u[0, :] = 0.0
    v[0, :] = 0.0
    u[:, 0] = 0.0
    v[:, 0] = 0.0
    u[:, -1] = 0.0
    v[:, -1] = 0.0


def main():
    # Problem parameters and grid
    alpha = 2.0
    Nx = 64  # points in x
    Ny = 64  # points in y
    Lx = 1.0
    Ly = 1.0
    dx = Lx / (Nx - 1)
    dy = Ly / (Ny - 1)
    assert abs(dx - dy) < 1e-12, "Currently assumes uniform square grid"

    # Time-stepping / pseudo-time parameters
    max_iters = 1200
    dt = 1e-3  # fixed pseudo-time step; chosen conservatively for stability
    omega = 0.75  # under-relaxation for velocity

    # Convergence tolerances
    tol_u = 1e-6
    tol_div = 1e-6

    # Build grid coordinates
    x = np.linspace(0.0, Lx, Nx)
    y = np.linspace(0.0, Ly, Ny)
    X, Y = np.meshgrid(x, y)

    # Initialize fields
    u = np.zeros((Ny, Nx), dtype=np.float64)
    v = np.zeros((Ny, Nx), dtype=np.float64)
    p = np.zeros((Ny, Nx), dtype=np.float64)

    # Apply BCs initially
    apply_velocity_bcs(u, v, alpha, Nx, Ny, x)

    # Pre-assemble matrices and LU factorizations
    # Helmholtz for velocities (interior unknowns)
    A_helm, bc_coeffs = build_helmholtz_matrix_for_vel(Nx, Ny, dx, dt, visc)
    # pre-factorize
    A_lu = spla.splu(A_helm)

    # Pressure Laplacian with Neumann
    Lp = build_pressure_laplacian_neumann(Nx, Ny, dx)
    # Enforce reference p(0,0)=0 by replacing first row with identity; create modified copy
    Lp_mod = Lp.tolil()
    ref_idx = 0  # (i=0,j=0) flattened at idx 0
    Lp_mod[ref_idx, :] = 0.0
    Lp_mod[ref_idx, ref_idx] = 1.0
    Lp_mod = Lp_mod.tocsc()
    Lp_lu = spla.splu(Lp_mod)

    # Precompute mapping of boundary contributions for velocities: bc_coeffs holds sum(coeff * boundary_value)
    # But boundary values differ for u and v (top lid nonzero in u). We precompute static boundary-value contributions vectors.
    # For each interior node, bc_coeffs is the sum of neighbor-coeffs (1/dx2) for each missing neighbor; However in our Dirichlet assembly we added neighbors that are on the boundary as coefficients in bc_coeffs. To compute bc contribution we must weight by actual boundary values of the particular boundary neighbor positions.
    # For simplicity, recompute bc contribution explicitly using the stencil used in assembly: iterate neighbors and if neighbor is boundary multiply coefficient by boundary value.

    nx_in = Nx - 2
    ny_in = Ny - 2
    N_in = nx_in * ny_in

    # We will build a matrix B such that bc_vec = B @ boundary_values_flat; but since BCs are simple and fixed, easier to compute bc_vec_u and bc_vec_v now.
    bc_vec_u = np.zeros(N_in, dtype=np.float64)
    bc_vec_v = np.zeros(N_in, dtype=np.float64)

    # For each interior node, check its four neighbors; if neighbor is boundary, add coeff*(boundary value)
    dx2 = dx * dx
    def bidx(i, j):
        return (j - 1) * nx_in + (i - 1)

    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            k = bidx(i, j)
            # left neighbor (i-1)
            if i - 1 == 0:
                # boundary at (i-1,j)
                val_u = 0.0
                val_v = 0.0
                bc_vec_u[k] += (1.0 / dx2) * val_u
                bc_vec_v[k] += (1.0 / dx2) * val_v
            # right neighbor (i+1)
            if i + 1 == Nx - 1:
                val_u = 0.0
                val_v = 0.0
                bc_vec_u[k] += (1.0 / dx2) * val_u
                bc_vec_v[k] += (1.0 / dx2) * val_v
            # bottom neighbor (j-1)
            if j - 1 == 0:
                val_u = 0.0
                val_v = 0.0
                bc_vec_u[k] += (1.0 / dx2) * val_u
                bc_vec_v[k] += (1.0 / dx2) * val_v
            # top neighbor (j+1)
            if j + 1 == Ny - 1:
                # boundary top: moving lid for u
                val_u = alpha * x[i] * (1.0 - x[i])
                val_v = 0.0
                bc_vec_u[k] += (1.0 / dx2) * val_u
                bc_vec_v[k] += (1.0 / dx2) * val_v

    # Note: Because our earlier assembly for interior Laplacian placed coefficient -center etc, we also already included center coefficient via A matrix. The bc_vec built here accounts for the neighbor contributions from Dirichlet boundaries.

    # Precompute flattened index mapping helpers
    def flatten(field):
        return flatten_interior(field)

    # Residual history
    res_u_hist = []
    res_div_hist = []

    # Progress reporting schedule
    report_every = max(1, max_iters // 10)

    converged = False
    for it in range(1, max_iters + 1):
        u_old = u.copy()
        v_old = v.copy()
        p_old = p.copy()

        # Compute convective terms evaluated at u_old, v_old
        conv_u, conv_v = compute_convective(u_old, v_old, dx)

        # Pressure gradient at interior (from current p)
        dpdx, dpdy = gradient_flattened_scalar(p_old, dx)

        # RHS for Helmholtz solve for u and v (interior unknowns)
        u_in = flatten(u_old)
        v_in = flatten(v_old)

        RHS_u = u_in - dt * (conv_u + dpdx) + dt * visc * bc_vec_u
        RHS_v = v_in - dt * (conv_v + dpdy) + dt * visc * bc_vec_v

        # Solve Helmholtz: A * u_star = RHS_u
        u_star_in = A_lu.solve(RHS_u)
        v_star_in = A_lu.solve(RHS_v)

        # Reconstruct full arrays with boundary values
        u_star = np.zeros_like(u)
        v_star = np.zeros_like(v)
        # place interior values
        u_star += unflatten_interior(u_star_in, Nx, Ny)
        v_star += unflatten_interior(v_star_in, Nx, Ny)
        # apply Dirichlet BCs (fixed)
        apply_velocity_bcs(u_star, v_star, alpha, Nx, Ny, x)

        # Compute divergence of u_star over full domain
        div_star = divergence_full(u_star, v_star, dx)
        rhs_p = (div_star / dt).ravel()

        # Enforce reference pressure right side zero at ref_idx
        rhs_p_mod = rhs_p.copy()
        rhs_p_mod[ref_idx] = 0.0

        # Solve Poisson for phi
        phi_flat = Lp_lu.solve(rhs_p_mod)
        phi = phi_flat.reshape((Ny, Nx))

        # Update pressure
        p_new = p_old + phi

        # Correct velocities: u^{n+1} = u_star - dt * grad(phi)
        dphidx_flat, dphidy_flat = gradient_flattened_scalar(phi, dx)
        u_corr_in = u_star_in - dt * dphidx_flat
        v_corr_in = v_star_in - dt * dphidy_flat

        # Unflatten and rebuild full arrays
        u_new = np.zeros_like(u)
        v_new = np.zeros_like(v)
        u_new += unflatten_interior(u_corr_in, Nx, Ny)
        v_new += unflatten_interior(v_corr_in, Nx, Ny)
        # apply Dirichlet velocity BCs
        apply_velocity_bcs(u_new, v_new, alpha, Nx, Ny, x)

        # Under-relaxation
        u = omega * u_new + (1.0 - omega) * u_old
        v = omega * v_new + (1.0 - omega) * v_old
        p = p_new  # pressure not relaxed here

        # Compute convergence metrics
        du = np.linalg.norm((u - u_old).ravel(), 2)
        dv = np.linalg.norm((v - v_old).ravel(), 2)
        res_u = max(du, dv)
        div_field = divergence_full(u, v, dx)
        max_div = np.max(np.abs(div_field))

        res_u_hist.append(res_u)
        res_div_hist.append(max_div)

        # Progress printing every 10% of max_iters
        if it % report_every == 0 or it == 1 or it == max_iters:
            print(f"Iter {it}/{max_iters}: vel_change_norm={res_u:.3e}, max_div={max_div:.3e}")

        # Convergence criteria
        if res_u < tol_u and max_div < tol_div:
            converged = True
            print(f"Converged at iteration {it}: vel_change_norm={res_u:.3e}, max_div={max_div:.3e}")
            break

    if not converged:
        # assess whether residuals have been decreasing
        dec = np.less(res_u_hist[-1], res_u_hist[0]) if len(res_u_hist) > 1 else False
        if dec:
            print(f"Reached max_iters={max_iters} but residuals decreased from {res_u_hist[0]:.2e} to {res_u_hist[-1]:.2e} - accept last iterate.")
        else:
            print(f"Reached max_iters={max_iters} and residuals did not decrease sufficiently. Last vel_change={res_u_hist[-1]:.3e}")

    # Post-processing: compute velocity magnitude
    vel_mag = np.sqrt(u ** 2 + v ** 2)

    # Plotting: left contour of velocity magnitude with streamlines, right convergence history
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    cf = ax.contourf(X, Y, vel_mag, levels=50, cmap='RdBu_r')
    ax.streamplot(X, Y, u, v, color='k', density=1.2, linewidth=0.5)
    ax.set_title('Velocity magnitude and streamlines')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    fig.colorbar(cf, ax=ax, label='|u|')

    ax2 = axes[1]
    iters = np.arange(1, len(res_u_hist) + 1)
    ax2.semilogy(iters, res_u_hist, '-o', markersize=3, label='vel-change norm')
    ax2.semilogy(iters, res_div_hist, '-s', markersize=3, label='max divergence')
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('Residual (log scale)')
    ax2.set_title('Convergence history')
    ax2.legend()
    ax2.grid(True, which='both', ls='--')

    fig.tight_layout()
    fig.savefig('navier_stokes_lid_cavity.png', dpi=200)

    # Print final quantitative results
    print(f"Final iteration: {it}")
    print(f"Final vel_change_norm = {res_u_hist[-1]:.3e}")
    print(f"Final max_divergence = {res_div_hist[-1]:.3e}")
    print(f"Mean velocity magnitude = {vel_mag.mean():.6f}")


if __name__ == '__main__':
    main()
```}


#### Script block2:

```python
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# Global viscosity from Reynolds number
Re = 400.0
visc = 1.0 / Re  # MUST use 'visc' everywhere


def build_interior_laplacian_dirichlet(Nx, Ny, dx):
    """
    Build 2D Laplacian matrix for interior nodes only (Dirichlet BC on boundaries).
    Discrete Laplacian L such that L(u)_i = (sum neighbors - 4*u_i)/dx^2.
    Interior nodes: i=1..Nx-2, j=1..Ny-2 (0-based global grid indexing).
    Returns L_int in CSC format (shape = N_in x N_in).
    """
    nx_in = Nx - 2
    ny_in = Ny - 2
    N = nx_in * ny_in
    dx2 = dx * dx

    rows = []
    cols = []
    data = []

    def idx(i, j):
        # i,j are interior-local indices in [0,nx_in-1]x[0,ny_in-1]
        return j * nx_in + i

    # neighbor offsets (left, right, down, up)
    offsets = [(-1, 0), (1, 0), (0, -1), (0, 1)]

    for j in range(ny_in):
        for i in range(nx_in):
            row = idx(i, j)
            center = 0.0
            # For each neighbor, if neighbor is interior add +1/dx2 off-diagonal
            for di, dj in offsets:
                ni = i + di
                nj = j + dj
                if 0 <= ni < nx_in and 0 <= nj < ny_in:
                    # interior neighbor
                    rows.append(row)
                    cols.append(idx(ni, nj))
                    data.append(1.0 / dx2)
                    center -= 1.0 / dx2
                else:
                    # neighbor is boundary: its contribution will be moved to RHS (handled outside)
                    center -= 1.0 / dx2
            # append center (which is negative: -4/dx2 for interior node)
            rows.append(row)
            cols.append(row)
            data.append(center)

    L = sp.csr_matrix((data, (rows, cols)), shape=(N, N), dtype=np.float64)
    return L.tocsc()


def build_pressure_laplacian_neumann(Nx, Ny, dx):
    """
    Build 2D Laplacian matrix for full grid with Neumann BC (zero normal derivative) handled by mirroring ghost nodes.
    Assembles Lp such that Lp(p) approximates Laplacian p.
    """
    N = Nx * Ny
    dx2 = dx * dx

    rows = []
    cols = []
    data = []

    def idx(i, j):
        return j * Nx + i

    # offsets: left,right,down,up in global indices
    offsets = [(-1, 0), (1, 0), (0, -1), (0, 1)]

    for j in range(Ny):
        for i in range(Nx):
            row = idx(i, j)
            center = 0.0
            for di, dj in offsets:
                ni = i + di
                nj = j + dj
                if 0 <= ni < Nx and 0 <= nj < Ny:
                    # regular neighbor
                    rows.append(row)
                    cols.append(idx(ni, nj))
                    data.append(1.0 / dx2)
                    center -= 1.0 / dx2
                else:
                    # Neumann mirror: ghost equals mirror interior neighbor
                    # Find the mirror interior neighbor (reflect across boundary): mirror by reversing offset
                    mi = i - di
                    mj = j - dj
                    # mi,mj should be a valid interior neighbor (since grid size > 1)
                    if 0 <= mi < Nx and 0 <= mj < Ny:
                        rows.append(row)
                        cols.append(idx(mi, mj))
                        data.append(1.0 * (1.0 / dx2))
                        center -= 1.0 * (1.0 / dx2)
                    else:
                        # This should not happen for standard Nx,Ny > 1, but handle gracefully
                        pass
            rows.append(row)
            cols.append(row)
            data.append(center)

    Lp = sp.csr_matrix((data, (rows, cols)), shape=(N, N), dtype=np.float64)
    return Lp.tocsc()


def build_helmholtz_matrix_for_vel(Nx, Ny, dx, dt, visc):
    L_int = build_interior_laplacian_dirichlet(Nx, Ny, dx)
    A = sp.eye(L_int.shape[0], format='csc') - (dt * visc) * L_int
    return A, L_int


def flatten_interior(field):
    Ny, Nx = field.shape
    nx_in = Nx - 2
    ny_in = Ny - 2
    vec = np.empty(nx_in * ny_in, dtype=np.float64)
    k = 0
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            vec[k] = field[j, i]
            k += 1
    return vec


def unflatten_interior(vec, Nx, Ny):
    nx_in = Nx - 2
    ny_in = Ny - 2
    field = np.zeros((Ny, Nx), dtype=np.float64)
    k = 0
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            field[j, i] = vec[k]
            k += 1
    return field


def compute_convective(u, v, dx):
    """
    Central differences for convective term at interior nodes.
    Returns flattened convective terms for u and v (interior only).
    """
    Ny, Nx = u.shape
    nx_in = Nx - 2
    ny_in = Ny - 2
    conv_u = np.empty(nx_in * ny_in, dtype=np.float64)
    conv_v = np.empty(nx_in * ny_in, dtype=np.float64)
    two_dx = 2.0 * dx
    k = 0
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            du_dx = (u[j, i + 1] - u[j, i - 1]) / two_dx
            du_dy = (u[j + 1, i] - u[j - 1, i]) / two_dx
            dv_dx = (v[j, i + 1] - v[j, i - 1]) / two_dx
            dv_dy = (v[j + 1, i] - v[j - 1, i]) / two_dx
            ui = u[j, i]
            vi = v[j, i]
            conv_u[k] = ui * du_dx + vi * du_dy
            conv_v[k] = ui * dv_dx + vi * dv_dy
            k += 1
    return conv_u, conv_v


def gradient_flattened_scalar(phi, dx):
    """
    Compute gradient of scalar phi (full array size Ny x Nx) and return dp/dx and dp/dy at interior nodes flattened.
    """
    Ny, Nx = phi.shape
    nx_in = Nx - 2
    ny_in = Ny - 2
    dphidx = np.empty(nx_in * ny_in, dtype=np.float64)
    dphidy = np.empty(nx_in * ny_in, dtype=np.float64)
    two_dx = 2.0 * dx
    k = 0
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            dphidx[k] = (phi[j, i + 1] - phi[j, i - 1]) / two_dx
            dphidy[k] = (phi[j + 1, i] - phi[j - 1, i]) / two_dx
            k += 1
    return dphidx, dphidy


def divergence_full(u, v, dx):
    """
    Compute divergence at all nodes using central differences for interior and one-sided for boundaries.
    Returns div array shape (Ny, Nx).
    """
    Ny, Nx = u.shape
    div = np.zeros_like(u)
    # interior
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            dudx = (u[j, i + 1] - u[j, i - 1]) / (2.0 * dx)
            dvdy = (v[j + 1, i] - v[j - 1, i]) / (2.0 * dx)
            div[j, i] = dudx + dvdy
    # boundaries: use one-sided differences
    for j in range(Ny):
        # left i=0
        dudx = (u[j, 1] - u[j, 0]) / dx
        if j == 0:
            dvdy = (v[1, 0] - v[0, 0]) / dx
        elif j == Ny - 1:
            dvdy = (v[-1, 0] - v[-2, 0]) / dx
        else:
            dvdy = (v[j + 1, 0] - v[j - 1, 0]) / (2.0 * dx)
        div[j, 0] = dudx + dvdy
        # right i=Nx-1
        dudx = (u[j, -1] - u[j, -2]) / dx
        if j == 0:
            dvdy = (v[1, -1] - v[0, -1]) / dx
        elif j == Ny - 1:
            dvdy = (v[-1, -1] - v[-2, -1]) / dx
        else:
            dvdy = (v[j + 1, -1] - v[j - 1, -1]) / (2.0 * dx)
        div[j, -1] = dudx + dvdy
    for i in range(Nx):
        # bottom j=0
        if i == 0:
            dudx = (u[0, 1] - u[0, 0]) / dx
        elif i == Nx - 1:
            dudx = (u[0, -1] - u[0, -2]) / dx
        else:
            dudx = (u[0, i + 1] - u[0, i - 1]) / (2.0 * dx)
        dvdy = (v[1, i] - v[0, i]) / dx
        div[0, i] = dudx + dvdy
        # top j=Ny-1
        if i == 0:
            dudx = (u[-1, 1] - u[-1, 0]) / dx
        elif i == Nx - 1:
            dudx = (u[-1, -1] - u[-1, -2]) / dx
        else:
            dudx = (u[-1, i + 1] - u[-1, i - 1]) / (2.0 * dx)
        dvdy = (v[-1, i] - v[-2, i]) / dx
        div[-1, i] = dudx + dvdy
    return div


def apply_velocity_bcs(u, v, alpha, Nx, Ny, x_coords):
    # Apply Dirichlet BCs: top boundary (moving lid) and no-slip elsewhere
    jtop = Ny - 1
    u[jtop, :] = alpha * x_coords * (1.0 - x_coords)
    v[jtop, :] = 0.0
    # other boundaries set to zero
    u[0, :] = 0.0
    v[0, :] = 0.0
    u[:, 0] = 0.0
    v[:, 0] = 0.0
    u[:, -1] = 0.0
    v[:, -1] = 0.0


def compute_bc_vecs(Nx, Ny, dx, alpha, x_coords):
    """
    Compute boundary contribution vectors for u and v interior RHSs.
    Returns bc_vec_u and bc_vec_v of length N_in where each entry is the NEGATIVE of (1/dx^2 * boundary_value)
    accumulated for missing neighbors (so RHS can add dt*visc*bc_vec directly).
    """
    nx_in = Nx - 2
    ny_in = Ny - 2
    N_in = nx_in * ny_in
    bc_vec_u = np.zeros(N_in, dtype=np.float64)
    bc_vec_v = np.zeros(N_in, dtype=np.float64)
    dx2 = dx * dx

    def bidx(i, j):
        return (j - 1) * nx_in + (i - 1)

    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            k = bidx(i, j)
            # left neighbor
            if i - 1 == 0:
                val_u = 0.0
                val_v = 0.0
                bc_vec_u[k] -= (1.0 / dx2) * val_u
                bc_vec_v[k] -= (1.0 / dx2) * val_v
            # right neighbor
            if i + 1 == Nx - 1:
                val_u = 0.0
                val_v = 0.0
                bc_vec_u[k] -= (1.0 / dx2) * val_u
                bc_vec_v[k] -= (1.0 / dx2) * val_v
            # bottom neighbor
            if j - 1 == 0:
                val_u = 0.0
                val_v = 0.0
                bc_vec_u[k] -= (1.0 / dx2) * val_u
                bc_vec_v[k] -= (1.0 / dx2) * val_v
            # top neighbor
            if j + 1 == Ny - 1:
                val_u = alpha * x_coords[i] * (1.0 - x_coords[i])
                val_v = 0.0
                bc_vec_u[k] -= (1.0 / dx2) * val_u
                bc_vec_v[k] -= (1.0 / dx2) * val_v
    return bc_vec_u, bc_vec_v


def main():
    # Problem parameters and grid
    alpha = 2.0
    Nx = 64  # points in x
    Ny = 64  # points in y
    Lx = 1.0
    Ly = 1.0
    dx = Lx / (Nx - 1)
    dy = Ly / (Ny - 1)
    assert abs(dx - dy) < 1e-12, "Currently assumes uniform square grid"

    # Time-stepping / pseudo-time parameters
    max_iters = 1200
    dt_user = 1e-3  # user's preferred dt; will be limited by CFL below
    cfl = 0.4
    omega = 0.75  # under-relaxation for velocity

    # Convergence tolerances
    tol_u = 1e-6
    tol_div = 1e-6

    # Build grid coordinates
    x = np.linspace(0.0, Lx, Nx)
    y = np.linspace(0.0, Ly, Ny)
    X, Y = np.meshgrid(x, y)

    # Initialize fields
    u = np.zeros((Ny, Nx), dtype=np.float64)
    v = np.zeros((Ny, Nx), dtype=np.float64)
    p = np.zeros((Ny, Nx), dtype=np.float64)

    # Apply BCs initially
    apply_velocity_bcs(u, v, alpha, Nx, Ny, x)

    # Choose dt conservatively based on lid velocity (max expected ~ alpha*0.25)
    u_lid_max = alpha * 0.25
    dt_cfl = cfl * dx / (u_lid_max + 1e-12)
    dt = min(dt_user, dt_cfl)

    # Pre-assemble matrices and LU factorizations
    A_helm, L_int = build_helmholtz_matrix_for_vel(Nx, Ny, dx, dt, visc)
    A_lu = spla.splu(A_helm)

    # Pressure Laplacian with Neumann
    Lp = build_pressure_laplacian_neumann(Nx, Ny, dx)

    # Quick consistency check: Lp applied to constant should be ~0 (Neumann nullspace)
    ones = np.ones(Nx * Ny)
    max_abs_Lp_ones = np.max(np.abs(Lp.dot(ones)))
    if max_abs_Lp_ones > 1e-10:
        print(f"Warning: pressure Laplacian applied to constant has max abs = {max_abs_Lp_ones:.3e}")

    # Enforce reference p(0,0)=0 by replacing first row with identity; create modified copy
    Lp_mod = Lp.tolil()
    ref_idx = 0  # (i=0,j=0) flattened at idx 0
    Lp_mod[ref_idx, :] = 0.0
    Lp_mod[ref_idx, ref_idx] = 1.0
    Lp_mod = Lp_mod.tocsc()
    Lp_lu = spla.splu(Lp_mod)

    # Precompute boundary contributions for velocity Helmholtz RHS (constant BCs)
    bc_vec_u, bc_vec_v = compute_bc_vecs(Nx, Ny, dx, alpha, x)

    # Residual history
    res_u_hist = []
    res_div_hist = []

    # Progress reporting schedule
    report_every = max(1, max_iters // 10)

    converged = False
    for it in range(1, max_iters + 1):
        u_old = u.copy()
        v_old = v.copy()
        p_old = p.copy()

        # Compute convective terms evaluated at u_old, v_old
        conv_u, conv_v = compute_convective(u_old, v_old, dx)

        # Pressure gradient at interior (from current p)
        dpdx, dpdy = gradient_flattened_scalar(p_old, dx)

        # Flatten interior velocities
        u_in = flatten_interior(u_old)
        v_in = flatten_interior(v_old)

        # RHS for Helmholtz solve for u and v (interior unknowns)
        # Note: bc_vec_u and bc_vec_v were computed as negative sums of boundary contributions -> add dt*visc*bc_vec
        RHS_u = u_in - dt * (conv_u + dpdx) + dt * visc * bc_vec_u
        RHS_v = v_in - dt * (conv_v + dpdy) + dt * visc * bc_vec_v

        # Solve Helmholtz: A * u_star = RHS_u
        u_star_in = A_lu.solve(RHS_u)
        v_star_in = A_lu.solve(RHS_v)

        # Reconstruct full arrays with boundary values
        u_star = np.zeros_like(u)
        v_star = np.zeros_like(v)
        u_star += unflatten_interior(u_star_in, Nx, Ny)
        v_star += unflatten_interior(v_star_in, Nx, Ny)
        apply_velocity_bcs(u_star, v_star, alpha, Nx, Ny, x)

        # Compute divergence of u_star over full domain
        div_star = divergence_full(u_star, v_star, dx)
        rhs_p = (div_star / dt).ravel()

        # Enforce reference pressure right side zero at ref_idx
        rhs_p_mod = rhs_p.copy()
        rhs_p_mod[ref_idx] = 0.0

        # Solve Poisson for phi
        phi_flat = Lp_lu.solve(rhs_p_mod)
        phi = phi_flat.reshape((Ny, Nx))

        # Update pressure
        p_new = p_old + phi

        # Correct velocities: u^{n+1} = u_star - dt * grad(phi)
        dphidx_flat, dphidy_flat = gradient_flattened_scalar(phi, dx)
        u_corr_in = u_star_in - dt * dphidx_flat
        v_corr_in = v_star_in - dt * dphidy_flat

        # Unflatten and rebuild full arrays
        u_new = np.zeros_like(u)
        v_new = np.zeros_like(v)
        u_new += unflatten_interior(u_corr_in, Nx, Ny)
        v_new += unflatten_interior(v_corr_in, Nx, Ny)
        apply_velocity_bcs(u_new, v_new, alpha, Nx, Ny, x)

        # Under-relaxation
        u = omega * u_new + (1.0 - omega) * u_old
        v = omega * v_new + (1.0 - omega) * v_old
        p = p_new  # pressure not relaxed here

        # Convergence metrics
        du = np.linalg.norm((u - u_old).ravel(), 2)
        dv = np.linalg.norm((v - v_old).ravel(), 2)
        res_u = max(du, dv)
        div_field = divergence_full(u, v, dx)
        max_div = np.max(np.abs(div_field))

        res_u_hist.append(res_u)
        res_div_hist.append(max_div)

        # Progress printing every 10% of max_iters (and first/last)
        if it % report_every == 0 or it == 1 or it == max_iters:
            print(f"Iter {it}/{max_iters}: vel_change_norm={res_u:.3e}, max_div={max_div:.3e}")

        # Check for NaNs or Infs and abort with diagnostics
        if not np.isfinite(res_u) or not np.isfinite(max_div):
            print("Numerical blow-up detected (NaN/Inf). Diagnostics:")
            print(f"  u_star_in: min={np.min(u_star_in):.3e}, max={np.max(u_star_in):.3e}")
            print(f"  phi_flat: min={np.min(phi_flat):.3e}, max={np.max(phi_flat):.3e}")
            break

        # Convergence criteria
        if res_u < tol_u and max_div < tol_div:
            converged = True
            print(f"Converged at iteration {it}: vel_change_norm={res_u:.3e}, max_div={max_div:.3e}")
            break

    if not converged:
        if len(res_u_hist) > 1 and res_u_hist[-1] < res_u_hist[0]:
            print(f"Reached max_iters={max_iters} but residuals decreased from {res_u_hist[0]:.2e} to {res_u_hist[-1]:.2e} - accept last iterate.")
        else:
            print(f"Reached max_iters={max_iters} and residuals did not decrease sufficiently. Last vel_change={res_u_hist[-1]:.3e}")

    # Post-processing: compute velocity magnitude
    vel_mag = np.sqrt(u ** 2 + v ** 2)

    # Plotting: left contour of velocity magnitude with streamlines, right convergence history
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    cf = ax.contourf(X, Y, vel_mag, levels=50, cmap='RdBu_r')
    ax.streamplot(X, Y, u, v, color='k', density=1.2, linewidth=0.5)
    ax.set_title('Velocity magnitude and streamlines')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    fig.colorbar(cf, ax=ax, label='|u|')

    ax2 = axes[1]
    iters = np.arange(1, len(res_u_hist) + 1)
    if len(iters) > 0:
        ax2.semilogy(iters, res_u_hist, '-o', markersize=3, label='vel-change norm')
        ax2.semilogy(iters, res_div_hist, '-s', markersize=3, label='max divergence')
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('Residual (log scale)')
    ax2.set_title('Convergence history')
    ax2.legend()
    ax2.grid(True, which='both', ls='--')

    fig.tight_layout()
    fig.savefig('navier_stokes_lid_cavity.png', dpi=200)

    # Print final quantitative results
    print(f"Final iteration: {it}")
    if len(res_u_hist) > 0:
        print(f"Final vel_change_norm = {res_u_hist[-1]:.3e}")
        print(f"Final max_divergence = {res_div_hist[-1]:.3e}")
    print(f"Mean velocity magnitude = {np.mean(vel_mag):.6f}")


if __name__ == '__main__':
    main()
```



#### Script block3:

```python
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# Global viscosity from Reynolds number
Re = 400.0
visc = 1.0 / Re  # MUST use 'visc' everywhere


def build_interior_laplacian_dirichlet(Nx, Ny, dx):
    """
    Build 2D Laplacian for interior nodes only (Dirichlet BC on boundaries).
    Discrete Laplacian L such that center entries = -sum(off-diagonals), off-diagonals = +1/dx^2.
    Ordering: row-major over interior grid (j then i), matching flatten_interior/unflatten_interior.
    Returns L_int in CSC format (shape = N_in x N_in).
    """
    nx_in = Nx - 2
    ny_in = Ny - 2
    N = nx_in * ny_in
    dx2 = dx * dx

    rows = []
    cols = []
    data = []

    def idx(i, j):
        return (j - 1) * nx_in + (i - 1)

    # iterate over interior nodes in global indices i=1..Nx-2, j=1..Ny-2
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            row = idx(i, j)
            center = 0.0
            # left
            if i - 1 >= 1:
                rows.append(row); cols.append(idx(i - 1, j)); data.append(1.0 / dx2)
                center -= 1.0 / dx2
            else:
                center -= 1.0 / dx2
            # right
            if i + 1 <= Nx - 2:
                rows.append(row); cols.append(idx(i + 1, j)); data.append(1.0 / dx2)
                center -= 1.0 / dx2
            else:
                center -= 1.0 / dx2
            # down
            if j - 1 >= 1:
                rows.append(row); cols.append(idx(i, j - 1)); data.append(1.0 / dx2)
                center -= 1.0 / dx2
            else:
                center -= 1.0 / dx2
            # up
            if j + 1 <= Ny - 2:
                rows.append(row); cols.append(idx(i, j + 1)); data.append(1.0 / dx2)
                center -= 1.0 / dx2
            else:
                center -= 1.0 / dx2
            rows.append(row); cols.append(row); data.append(center)

    L = sp.csr_matrix((data, (rows, cols)), shape=(N, N), dtype=np.float64)
    return L.tocsc()


def build_pressure_laplacian_neumann(Nx, Ny, dx):
    """
    Build 2D Laplacian on full grid with Neumann boundary treatment by mirroring.
    Uses row-major flattening idx = j*Nx + i, consistent with full-grid raveling.
    """
    N = Nx * Ny
    dx2 = dx * dx

    rows = []
    cols = []
    data = []

    def idx(i, j):
        return j * Nx + i

    for j in range(Ny):
        for i in range(Nx):
            row = idx(i, j)
            center = 0.0
            # left
            if i - 1 >= 0:
                rows.append(row); cols.append(idx(i - 1, j)); data.append(1.0 / dx2)
                center -= 1.0 / dx2
            else:
                # mirror to i+1
                rows.append(row); cols.append(idx(i + 1, j)); data.append(1.0 / dx2)
                center -= 1.0 / dx2
            # right
            if i + 1 < Nx:
                rows.append(row); cols.append(idx(i + 1, j)); data.append(1.0 / dx2)
                center -= 1.0 / dx2
            else:
                # mirror to i-1
                rows.append(row); cols.append(idx(i - 1, j)); data.append(1.0 / dx2)
                center -= 1.0 / dx2
            # down
            if j - 1 >= 0:
                rows.append(row); cols.append(idx(i, j - 1)); data.append(1.0 / dx2)
                center -= 1.0 / dx2
            else:
                # mirror to j+1
                rows.append(row); cols.append(idx(i, j + 1)); data.append(1.0 / dx2)
                center -= 1.0 / dx2
            # up
            if j + 1 < Ny:
                rows.append(row); cols.append(idx(i, j + 1)); data.append(1.0 / dx2)
                center -= 1.0 / dx2
            else:
                # mirror to j-1
                rows.append(row); cols.append(idx(i, j - 1)); data.append(1.0 / dx2)
                center -= 1.0 / dx2

            rows.append(row); cols.append(row); data.append(center)

    Lp = sp.csr_matrix((data, (rows, cols)), shape=(N, N), dtype=np.float64)
    return Lp.tocsc()


def build_helmholtz_matrix_for_vel(Nx, Ny, dx, dt, visc):
    L_int = build_interior_laplacian_dirichlet(Nx, Ny, dx)
    A = sp.eye(L_int.shape[0], format='csc') - (dt * visc) * L_int
    return A, L_int


def flatten_interior(field):
    Ny, Nx = field.shape
    nx_in = Nx - 2
    ny_in = Ny - 2
    vec = np.empty(nx_in * ny_in, dtype=np.float64)
    k = 0
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            vec[k] = field[j, i]
            k += 1
    return vec


def unflatten_interior(vec, Nx, Ny):
    nx_in = Nx - 2
    ny_in = Ny - 2
    field = np.zeros((Ny, Nx), dtype=np.float64)
    k = 0
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            field[j, i] = vec[k]
            k += 1
    return field


def compute_convective(u, v, dx):
    """
    Central differences for convective term at interior nodes; returns flattened convective terms.
    """
    Ny, Nx = u.shape
    nx_in = Nx - 2
    ny_in = Ny - 2
    conv_u = np.empty(nx_in * ny_in, dtype=np.float64)
    conv_v = np.empty(nx_in * ny_in, dtype=np.float64)
    two_dx = 2.0 * dx
    k = 0
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            du_dx = (u[j, i + 1] - u[j, i - 1]) / two_dx
            du_dy = (u[j + 1, i] - u[j - 1, i]) / two_dx
            dv_dx = (v[j, i + 1] - v[j, i - 1]) / two_dx
            dv_dy = (v[j + 1, i] - v[j - 1, i]) / two_dx
            ui = u[j, i]
            vi = v[j, i]
            conv_u[k] = ui * du_dx + vi * du_dy
            conv_v[k] = ui * dv_dx + vi * dv_dy
            k += 1
    return conv_u, conv_v


def gradient_flattened_scalar(phi, dx):
    """
    Compute gradient of scalar phi (full array Ny x Nx) and return dp/dx and dp/dy at interior nodes flattened.
    Ordering matches flatten_interior/unflatten_interior.
    """
    Ny, Nx = phi.shape
    nx_in = Nx - 2
    ny_in = Ny - 2
    dphidx = np.empty(nx_in * ny_in, dtype=np.float64)
    dphidy = np.empty(nx_in * ny_in, dtype=np.float64)
    two_dx = 2.0 * dx
    k = 0
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            dphidx[k] = (phi[j, i + 1] - phi[j, i - 1]) / two_dx
            dphidy[k] = (phi[j + 1, i] - phi[j - 1, i]) / two_dx
            k += 1
    return dphidx, dphidy


def divergence_full(u, v, dx):
    """
    Compute divergence at all nodes using central differences for interior and one-sided for boundaries.
    Returns div array shape (Ny, Nx).
    """
    Ny, Nx = u.shape
    div = np.zeros_like(u)
    # interior
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            dudx = (u[j, i + 1] - u[j, i - 1]) / (2.0 * dx)
            dvdy = (v[j + 1, i] - v[j - 1, i]) / (2.0 * dx)
            div[j, i] = dudx + dvdy
    # boundaries: use one-sided differences
    for j in range(Ny):
        # left i=0
        dudx = (u[j, 1] - u[j, 0]) / dx
        if j == 0:
            dvdy = (v[1, 0] - v[0, 0]) / dx
        elif j == Ny - 1:
            dvdy = (v[-1, 0] - v[-2, 0]) / dx
        else:
            dvdy = (v[j + 1, 0] - v[j - 1, 0]) / (2.0 * dx)
        div[j, 0] = dudx + dvdy
        # right i=Nx-1
        dudx = (u[j, -1] - u[j, -2]) / dx
        if j == 0:
            dvdy = (v[1, -1] - v[0, -1]) / dx
        elif j == Ny - 1:
            dvdy = (v[-1, -1] - v[-2, -1]) / dx
        else:
            dvdy = (v[j + 1, -1] - v[j - 1, -1]) / (2.0 * dx)
        div[j, -1] = dudx + dvdy
    for i in range(Nx):
        # bottom j=0
        if i == 0:
            dudx = (u[0, 1] - u[0, 0]) / dx
        elif i == Nx - 1:
            dudx = (u[0, -1] - u[0, -2]) / dx
        else:
            dudx = (u[0, i + 1] - u[0, i - 1]) / (2.0 * dx)
        dvdy = (v[1, i] - v[0, i]) / dx
        div[0, i] = dudx + dvdy
        # top j=Ny-1
        if i == 0:
            dudx = (u[-1, 1] - u[-1, 0]) / dx
        elif i == Nx - 1:
            dudx = (u[-1, -1] - u[-1, -2]) / dx
        else:
            dudx = (u[-1, i + 1] - u[-1, i - 1]) / (2.0 * dx)
        dvdy = (v[-1, i] - v[-2, i]) / dx
        div[-1, i] = dudx + dvdy
    return div


def apply_velocity_bcs(u, v, alpha, Nx, Ny, x_coords):
    # Apply Dirichlet BCs: top boundary (moving lid) and no-slip elsewhere
    jtop = Ny - 1
    u[jtop, :] = alpha * x_coords * (1.0 - x_coords)
    v[jtop, :] = 0.0
    # other boundaries set to zero (no-slip)
    u[0, :] = 0.0
    v[0, :] = 0.0
    u[:, 0] = 0.0
    v[:, 0] = 0.0
    u[:, -1] = 0.0
    v[:, -1] = 0.0


def compute_bc_vecs(Nx, Ny, dx, alpha, x_coords):
    """
    Compute boundary contribution vectors for u and v interior RHSs.
    Each entry equals the NEGATIVE of (1/dx^2 * boundary_value) accumulated for missing neighbors.
    Note: currently only top boundary has nonzero Dirichlet values; other boundaries are zero explicitly.
    Returns bc_vec_u and bc_vec_v of length N_in.
    """
    nx_in = Nx - 2
    ny_in = Ny - 2
    N_in = nx_in * ny_in
    bc_vec_u = np.zeros(N_in, dtype=np.float64)
    bc_vec_v = np.zeros(N_in, dtype=np.float64)
    dx2 = dx * dx

    def bidx(i, j):
        return (j - 1) * nx_in + (i - 1)

    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            k = bidx(i, j)
            # left neighbor (i-1 == 0) -> boundary
            if i - 1 == 0:
                val_u = 0.0
                val_v = 0.0
                bc_vec_u[k] -= (1.0 / dx2) * val_u
                bc_vec_v[k] -= (1.0 / dx2) * val_v
            # right neighbor (i+1 == Nx-1)
            if i + 1 == Nx - 1:
                val_u = 0.0
                val_v = 0.0
                bc_vec_u[k] -= (1.0 / dx2) * val_u
                bc_vec_v[k] -= (1.0 / dx2) * val_v
            # bottom neighbor (j-1 == 0)
            if j - 1 == 0:
                val_u = 0.0
                val_v = 0.0
                bc_vec_u[k] -= (1.0 / dx2) * val_u
                bc_vec_v[k] -= (1.0 / dx2) * val_v
            # top neighbor (j+1 == Ny-1)
            if j + 1 == Ny - 1:
                val_u = alpha * x_coords[i] * (1.0 - x_coords[i])
                val_v = 0.0
                bc_vec_u[k] -= (1.0 / dx2) * val_u
                bc_vec_v[k] -= (1.0 / dx2) * val_v
    return bc_vec_u, bc_vec_v


def main():
    # Problem parameters and grid
    alpha = 2.0
    Nx = 64  # points in x
    Ny = 64  # points in y
    Lx = 1.0
    Ly = 1.0
    dx = Lx / (Nx - 1)
    dy = Ly / (Ny - 1)
    assert abs(dx - dy) < 1e-12, "Currently assumes uniform square grid"

    # Time-stepping / pseudo-time parameters
    max_iters = 1200
    dt_user = 1e-3  # user's preferred dt; will be limited by CFL below
    cfl = 0.4
    omega = 0.75  # under-relaxation for velocity

    # Convergence tolerances (interior-based)
    tol_u = 1e-6
    tol_div = 1e-6

    # Build grid coordinates
    x = np.linspace(0.0, Lx, Nx)
    y = np.linspace(0.0, Ly, Ny)
    X, Y = np.meshgrid(x, y)

    # Initialize fields
    u = np.zeros((Ny, Nx), dtype=np.float64)
    v = np.zeros((Ny, Nx), dtype=np.float64)
    p = np.zeros((Ny, Nx), dtype=np.float64)

    # Apply BCs initially
    apply_velocity_bcs(u, v, alpha, Nx, Ny, x)

    # Choose dt conservatively based on lid velocity (max expected ~ alpha*0.25)
    u_lid_max = alpha * 0.25
    dt_cfl = cfl * dx / (u_lid_max + 1e-12)
    dt = min(dt_user, dt_cfl)

    # Pre-assemble matrices and LU factorizations (done once)
    A_helm, L_int = build_helmholtz_matrix_for_vel(Nx, Ny, dx, dt, visc)
    A_lu = spla.splu(A_helm)

    # Pressure Laplacian with Neumann
    Lp = build_pressure_laplacian_neumann(Nx, Ny, dx)

    # Enforce reference p(0,0)=0 by replacing first row with identity; create modified copy
    Lp_mod = Lp.tolil()
    ref_idx = 0  # (i=0,j=0) flattened at idx 0 in row-major
    Lp_mod[ref_idx, :] = 0.0
    Lp_mod[ref_idx, ref_idx] = 1.0
    Lp_mod = Lp_mod.tocsc()
    Lp_lu = spla.splu(Lp_mod)

    # Precompute boundary contributions for velocity Helmholtz RHS (constant BCs)
    bc_vec_u, bc_vec_v = compute_bc_vecs(Nx, Ny, dx, alpha, x)

    # Residual history (interior-based)
    vel_change_L2_hist = []
    vel_change_Linf_hist = []
    div_int_L2_hist = []
    div_int_Linf_hist = []
    poisson_res_hist = []

    # Progress reporting schedule
    report_every = max(1, max_iters // 10)

    converged = False
    # Compute initial divergence diagnostics
    div0 = divergence_full(u, v, dx)
    div0_int = div0[1:-1, 1:-1]
    div0_int_L2 = np.linalg.norm(div0_int.ravel(), 2)
    div0_int_Linf = np.max(np.abs(div0_int))

    for it in range(1, max_iters + 1):
        u_old = u.copy()
        v_old = v.copy()
        p_old = p.copy()

        # Compute convective terms evaluated at u_old, v_old
        conv_u, conv_v = compute_convective(u_old, v_old, dx)

        # Pressure gradient at interior (from current p)
        dpdx, dpdy = gradient_flattened_scalar(p_old, dx)

        # Flatten interior velocities
        u_in = flatten_interior(u_old)
        v_in = flatten_interior(v_old)

        # RHS for Helmholtz solve for u and v (interior unknowns)
        RHS_u = u_in - dt * (conv_u + dpdx) + dt * visc * bc_vec_u
        RHS_v = v_in - dt * (conv_v + dpdy) + dt * visc * bc_vec_v

        # Solve Helmholtz: A * u_star = RHS_u
        u_star_in = A_lu.solve(RHS_u)
        v_star_in = A_lu.solve(RHS_v)

        # Reconstruct full arrays with boundary values
        u_star = unflatten_interior(u_star_in, Nx, Ny)
        v_star = unflatten_interior(v_star_in, Nx, Ny)
        apply_velocity_bcs(u_star, v_star, alpha, Nx, Ny, x)

        # Compute divergence of u_star over full domain and interior diagnostics before projection
        div_star = divergence_full(u_star, v_star, dx)
        div_star_int = div_star[1:-1, 1:-1]
        div_star_int_L2 = np.linalg.norm(div_star_int.ravel(), 2)
        div_star_int_Linf = np.max(np.abs(div_star_int))

        # RHS for Poisson (full-grid ordering): div_star / dt flattened row-major
        rhs_p = (div_star / dt).ravel()
        rhs_p_mod = rhs_p.copy()
        rhs_p_mod[ref_idx] = 0.0

        # Solve Poisson for phi (with modified Lp that pins reference)
        phi_flat = Lp_lu.solve(rhs_p_mod)

        # Check Poisson residual using Lp_mod to be consistent with the solve
        poisson_res = Lp_mod.dot(phi_flat) - rhs_p_mod
        poisson_res_norm = np.linalg.norm(poisson_res, 2)

        # Reshape phi to full grid
        phi = phi_flat.reshape((Ny, Nx))

        # Update pressure and correct velocities by gradient of phi (interior flattened grads)
        p_new = p_old + phi
        dphidx_flat, dphidy_flat = gradient_flattened_scalar(phi, dx)
        u_corr_in = u_star_in - dt * dphidx_flat
        v_corr_in = v_star_in - dt * dphidy_flat

        # Unflatten and apply BCs
        u_new = unflatten_interior(u_corr_in, Nx, Ny)
        v_new = unflatten_interior(v_corr_in, Nx, Ny)
        apply_velocity_bcs(u_new, v_new, alpha, Nx, Ny, x)

        # Under-relaxation
        u = omega * u_new + (1.0 - omega) * u_old
        v = omega * v_new + (1.0 - omega) * v_old
        p = p_new

        # Convergence metrics (interior-focused)
        du_field = u - u_old
        dv_field = v - v_old
        # interior velocity change
        du_int = du_field[1:-1, 1:-1].ravel()
        dv_int = dv_field[1:-1, 1:-1].ravel()
        vel_change_L2 = np.linalg.norm(np.hstack((du_int, dv_int)), 2)
        vel_change_Linf = max(np.max(np.abs(du_int)), np.max(np.abs(dv_int)))

        div_field = divergence_full(u, v, dx)
        div_int = div_field[1:-1, 1:-1]
        div_int_L2 = np.linalg.norm(div_int.ravel(), 2)
        div_int_Linf = np.max(np.abs(div_int))

        # Record histories
        vel_change_L2_hist.append(vel_change_L2)
        vel_change_Linf_hist.append(vel_change_Linf)
        div_int_L2_hist.append(div_int_L2)
        div_int_Linf_hist.append(div_int_Linf)
        poisson_res_hist.append(poisson_res_norm)

        # Progress printing every 10% of max_iters (plus first and last)
        if (it % report_every == 0) or (it == 1) or (it == max_iters):
            print(f"Iter {it}/{max_iters}: vel_change_L2={vel_change_L2:.3e}, div_int_Linf={div_int_Linf:.3e}, poisson_res={poisson_res_norm:.3e}")

        # Check Poisson residual for consistency
        if poisson_res_norm > 1e-6:
            # If residual unexpectedly large, print a concise warning once
            if it == 1 or (it % (5 * report_every) == 0):
                print(f"  Warning: Poisson residual norm = {poisson_res_norm:.3e} (it={it})")

        # Check for NaNs or Infs and abort with diagnostics
        if not np.isfinite(vel_change_L2) or not np.isfinite(div_int_Linf):
            print("Numerical blow-up detected (NaN/Inf). Diagnostics:")
            print(f"  u_star_in: min={np.min(u_star_in):.3e}, max={np.max(u_star_in):.3e}")
            print(f"  phi_flat: min={np.min(phi_flat):.3e}, max={np.max(phi_flat):.3e}")
            break

        # Convergence criteria (interior-based)
        if (vel_change_L2 < tol_u) and (div_int_Linf < tol_div):
            converged = True
            print(f"Converged at iteration {it}: vel_change_L2={vel_change_L2:.3e}, div_int_Linf={div_int_Linf:.3e}")
            break

    # Post-loop acceptance logic
    if not converged:
        if len(vel_change_L2_hist) > 1 and vel_change_L2_hist[-1] < vel_change_L2_hist[0]:
            print(f"Reached max_iters={max_iters} but interior residuals decreased from {vel_change_L2_hist[0]:.2e} to {vel_change_L2_hist[-1]:.2e} - accept last iterate.")
        else:
            print(f"Reached max_iters={max_iters} and residuals did not decrease sufficiently. Last vel_change_L2={vel_change_L2_hist[-1]:.3e}")

    # Final diagnostic summary (print interior initial vs final)
    if len(div_int_L2_hist) > 0:
        print(f"Initial interior div L2 = {div0_int_L2:.3e}, final interior div L2 = {div_int_L2_hist[-1]:.3e}")
        print(f"Initial interior div Linf = {div0_int_Linf:.3e}, final interior div Linf = {div_int_Linf_hist[-1]:.3e}")
    if len(vel_change_L2_hist) > 0:
        print(f"Final vel_change_L2 = {vel_change_L2_hist[-1]:.3e}, Final vel_change_Linf = {vel_change_Linf_hist[-1]:.3e}")
    # Show Poisson residual history final
    if len(poisson_res_hist) > 0:
        print(f"Final Poisson residual norm = {poisson_res_hist[-1]:.3e}")

    # Post-processing: compute velocity magnitude
    vel_mag = np.sqrt(u ** 2 + v ** 2)

    # Plotting: left contour of velocity magnitude with streamlines, right convergence history
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    cf = ax.contourf(X, Y, vel_mag, levels=50, cmap='RdBu_r')
    ax.streamplot(X, Y, u, v, color='k', density=1.2, linewidth=0.5)
    ax.set_title('Velocity magnitude and streamlines')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    fig.colorbar(cf, ax=ax, label='|u|')

    ax2 = axes[1]
    iters = np.arange(1, len(vel_change_L2_hist) + 1)
    if len(iters) > 0:
        ax2.semilogy(iters, vel_change_L2_hist, '-o', markersize=3, label='vel-change L2')
        ax2.semilogy(iters, vel_change_Linf_hist, '-s', markersize=3, label='vel-change Linf')
        ax2.semilogy(iters, div_int_L2_hist, '-^', markersize=3, label='div interior L2')
        ax2.semilogy(iters, div_int_Linf_hist, '-v', markersize=3, label='div interior Linf')
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('Residual (log scale)')
    ax2.set_title('Convergence history (interior metrics)')
    ax2.legend()
    ax2.grid(True, which='both', ls='--')

    fig.tight_layout()
    fig.savefig('navier_stokes_lid_cavity.png', dpi=200)

    # Print final quantitative results
    print(f"Final iteration: {it}")
    if len(vel_change_L2_hist) > 0:
        print(f"Final vel_change_L2 = {vel_change_L2_hist[-1]:.3e}")
        print(f"Final div interior Linf = {div_int_Linf_hist[-1]:.3e}")
    print(f"Mean velocity magnitude = {np.mean(vel_mag):.6f}")


if __name__ == '__main__':
    main()
```



### runtime_outputs

#### Output block1

Iter 1/1200: vel_change_norm=2.278e-02, max_div=2.061e+00
<string>:209: RuntimeWarning: overflow encountered in scalar multiply
<string>:209: RuntimeWarning: invalid value encountered in scalar add
<string>:210: RuntimeWarning: overflow encountered in scalar multiply
<string>:210: RuntimeWarning: invalid value encountered in scalar add
Iter 120/1200: vel_change_norm=nan, max_div=nan
Iter 240/1200: vel_change_norm=nan, max_div=nan
Iter 360/1200: vel_change_norm=nan, max_div=nan
Iter 480/1200: vel_change_norm=nan, max_div=nan
Iter 600/1200: vel_change_norm=nan, max_div=nan
Iter 720/1200: vel_change_norm=nan, max_div=nan
Iter 840/1200: vel_change_norm=nan, max_div=nan
Iter 960/1200: vel_change_norm=nan, max_div=nan
Iter 1080/1200: vel_change_norm=nan, max_div=nan
Iter 1200/1200: vel_change_norm=nan, max_div=nan
Reached max_iters=1200 and residuals did not decrease sufficiently. Last vel_change=nan
Final iteration: 1200
Final vel_change_norm = nan
Final max_divergence = nan
Mean velocity magnitude = nan



#### Output block2

Iter 1/1200: vel_change_norm=1.741e-01, max_div=1.968e+00
Iter 120/1200: vel_change_norm=6.199e-03, max_div=1.968e+00
Iter 240/1200: vel_change_norm=3.619e-03, max_div=1.968e+00
Iter 360/1200: vel_change_norm=2.648e-03, max_div=1.968e+00
Iter 480/1200: vel_change_norm=2.108e-03, max_div=1.968e+00
Iter 600/1200: vel_change_norm=1.759e-03, max_div=1.968e+00
Iter 720/1200: vel_change_norm=1.513e-03, max_div=1.968e+00
Iter 840/1200: vel_change_norm=1.329e-03, max_div=1.968e+00
Iter 960/1200: vel_change_norm=1.186e-03, max_div=1.968e+00
Iter 1080/1200: vel_change_norm=1.070e-03, max_div=1.968e+00
Iter 1200/1200: vel_change_norm=9.753e-04, max_div=1.968e+00
Reached max_iters=1200 but residuals decreased from 1.74e-01 to 9.75e-04 - accept last iterate.
Final iteration: 1200
Final vel_change_norm = 9.753e-04
Final max_divergence = 1.968e+00
Mean velocity magnitude = 0.028449



#### Output block3

Iter 1/1200: vel_change_L2=2.362e-01, div_int_Linf=3.686e-01, poisson_res=1.816e-10
Iter 120/1200: vel_change_L2=6.520e-03, div_int_Linf=1.365e-02, poisson_res=9.281e-14
Iter 240/1200: vel_change_L2=3.813e-03, div_int_Linf=1.708e-03, poisson_res=1.199e-14
Iter 360/1200: vel_change_L2=2.825e-03, div_int_Linf=1.384e-03, poisson_res=7.678e-15
Iter 480/1200: vel_change_L2=2.273e-03, div_int_Linf=1.365e-03, poisson_res=2.479e-13
Iter 600/1200: vel_change_L2=1.915e-03, div_int_Linf=1.174e-03, poisson_res=5.742e-14
Iter 720/1200: vel_change_L2=1.663e-03, div_int_Linf=1.074e-03, poisson_res=8.081e-14
Iter 840/1200: vel_change_L2=1.474e-03, div_int_Linf=1.010e-03, poisson_res=8.077e-14
Iter 960/1200: vel_change_L2=1.327e-03, div_int_Linf=9.649e-04, poisson_res=8.077e-14
Iter 1080/1200: vel_change_L2=1.209e-03, div_int_Linf=9.325e-04, poisson_res=5.740e-14
Iter 1200/1200: vel_change_L2=1.111e-03, div_int_Linf=9.081e-04, poisson_res=5.744e-14
Reached max_iters=1200 but interior residuals decreased from 2.36e-01 to 1.11e-03 - accept last iterate.
Initial interior div L2 = 0.000e+00, final interior div L2 = 3.444e-02
Initial interior div Linf = 0.000e+00, final interior div Linf = 9.081e-04
Final vel_change_L2 = 1.111e-03, Final vel_change_Linf = 8.280e-05
Final Poisson residual norm = 5.744e-14
Final iteration: 1200
Final vel_change_L2 = 1.111e-03
Final div interior Linf = 9.081e-04
Mean velocity magnitude = 0.028449




### review_decision
accept

### review_comments
	Current Stage [A/3]
1. MAJOR

Category: accuracy
Issue: Discrete Laplacian matrices have incorrect sign convention (diagonal sign flipped). build_interior_laplacian_dirichlet and build_pressure_laplacian_neumann append -center instead of center, producing a matrix with positive diagonal and positive off-diagonals instead of the standard discrete Laplacian (diag = -4/dx2, off-diag = +1/dx2).

Recommendation: Fix the assembly: append the accumulated 'center' value directly (rows.append(row); cols.append(row); data.append(center)) instead of '-center'. Ensure the matrix implements L(u) = (u_{i-1}+u_{i+1}+u_{j-1}+u_{j+1}-4u_i)/dx^2. After correcting the sign, re-run and verify eigenvalues / positive-definiteness where expected (pressure Poisson should be SPD up to the reference row).


2. MAJOR

Category: accuracy
Issue: Boundary contribution sign is inconsistent with Laplacian assembly. bc_vec_u and bc_vec_v are computed as positive sums of boundary values ( +1/dx^2 * val ) but when the Laplacian is assembled in the standard form the RHS should receive the negative of these boundary contributions.

Recommendation: Change bc_vec computation to accumulate the negative boundary contribution, e.g. bc_vec[k] -= (1.0/dx2) * boundary_value for each missing neighbor. After fixing the Laplacian sign (previous item), ensure RHS for Helmholtz is RHS = u_in - dt*(conv + dpdx) - dt*visc*(bc_vec) if bc_vec holds positive sum, or equivalently compute bc_vec as negative and keep +dt*visc*bc_vec. Be consistent and test with a manufactured solution or simple Dirichlet test to ensure boundary terms are applied correctly.


3. MAJOR

Category: accuracy
Issue: Pressure Poisson (Lp) assembly is likely sign-inverted and the Neumann 'mirror' treatment is fragile. With the current sign errors and boundary treatment the Poisson solve produces incorrect pressure corrections leading to blow-up (overflow/NaN).

Recommendation: Use the same sign convention fix as for the interior Laplacian (append center, not -center). Re-check the Neumann implementation: simpler approach is to assemble the standard 5-point Laplacian for interior nodes and for boundary nodes apply one-sided second-derivative consistent discretization (or explicitly set ghost = neighbor to get doubled coefficient). After assembly, verify that the nullspace is handled correctly: you already fix the reference row which is OK. Validate Lp by applying it to a smooth test phi and compare to finite-difference Laplacian (unit test).


4. MAJOR

Category: accuracy
Issue: The algorithm became numerically unstable (overflow/NaN) after a few iterations — symptomatic of the sign/linear-system bugs and also of the explicit treatment of convection at moderate Re. The projection / time stepping formulation and RHS assembly must be checked for correct signs and consistency with the discrete operators.

Recommendation: After fixing matrix signs and BC handling, reduce dt or add an adaptive CFL-based dt selection initially (e.g. dt <= CFL * dx / max(|u|+|v|) with CFL ~ 0.2-0.5) to avoid blow-up from explicit advection. Consider using a skew-symmetric or upwind discretization for the convective term to improve stability at Re=400. Keep under-relaxation (omega) but consider applying mild relaxation to pressure as well if needed. Re-run and monitor residual curves (they should monotonically decrease if stable).


5. MINOR

Category: structure
Issue: Several implementation/structure issues increase fragility and complicate debugging: (a) build_interior_laplacian_dirichlet has awkward/incorrect neighbor condition branches (right-neighbor check is confusing), (b) bc_vec assembly logic mixes index arithmetic and equality checks (i-1==0 etc.) which is easy to get wrong, (c) comments and variable names are occasionally ambiguous (dx2 used for dx*dx in places and for 2*dx in others).

Recommendation: Refactor the Laplacian assembly to a clearer, symmetric stencil loop: iterate neighbors with explicit offsets [(−1,0),(+1,0),(0,−1),(0,+1)] and handle interior vs boundary with a small helper. Compute dx2 = dx*dx for second-derivative denominators and use a separate variable for 2*dx in central-difference denominators. Add a small unit-test routine (2D analytical solution) or a consistency check that A @ ones (or numerical Laplace) behaves as expected before the time loop. Keep prints at the requested frequency (you already did), and add short checks to stop early if NaNs appear and print diagnostics (max/min values of rhs, matrix diagonal range).



	Current Stage [B/3]
1. MAJOR

Category: accuracy
Issue: Max divergence is reported over the full grid and remains ~1.968 constant — this is caused by boundary nodes (Dirichlet BCs) dominating the max-norm, so the convergence diagnostic incorrectly indicates no improvement.

Recommendation: When monitoring divergence and reporting convergence you must exclude Dirichlet boundary nodes (or else the max-norm will be dominated by the enforced boundary values and not reflect interior mass conservation). Change the diagnostic to compute divergence only on interior (1:-1,1:-1) nodes, e.g. max_div = np.max(np.abs(div_field[1:-1,1:-1])) or use an L2-norm over the interior. After that, re-run and confirm the projection step is actually reducing the interior divergence.


2. MAJOR

Category: accuracy
Issue: Projection step appears ineffective (divergence not reduced). This can be caused by inconsistencies between discretizations/orderings of Lp, the flattened gradients, or the flattened RHS.

Recommendation: Verify that the projection solve and velocity correction use exactly the same discretization and ordering. Concretely: (a) ensure the pressure-Poisson operator Lp uses the same finite-difference stencil and index ordering as gradient_flattened_scalar and flatten_interior/unflatten_interior; (b) check that the RHS ordering (rhs_p.ravel()) matches Lp row ordering (you use row-major, so keep consistent); (c) after fixing the divergence metric per above, test that the Poisson solve reduces interior divergence by measuring ||div_star - dt*L(phi)|| or directly computing Lp.dot(phi_flat) - rhs_p_mod.


3. MINOR

Category: structure
Issue: Convergence diagnostics are currently limited to global max-div (which is misleading) and the global vel-change norm; richer interior metrics will help justify acceptance at max_iters.

Recommendation: When reporting convergence accept/reject logic is reasonable but improve diagnostics: record and print interior divergence norms (L2 and Linf) and also show initial and final numbers in the summary. Also consider using interior L2 of velocity-change for a smoother convergence indicator. This will make the 'accept at max_iters' decision more robust.


4. MINOR

Category: structure
Issue: Code contains some minor clarity and micro-performance issues (index arithmetic in hot loops, awkward additions) that do not change correctness but reduce readability/maintainability.

Recommendation: Small robustness/clarity fixes: (a) compute_bc_vecs currently treats all non-top boundaries as zero — add comments or explicit handling if boundary conditions change; (b) consider precomputing the interior index mapping (i.e. arrays of i/j indices) once to avoid repeated index arithmetic in flatten/unflatten/convective loops for readability and slight speed-up; (c) check and remove the redundant addition pattern u_star = np.zeros_like(u); u_star += unflatten_interior(...) (direct assignment is clearer).


5. MINOR

Category: accuracy
Issue: Central-difference convection is marginally stable at higher Reynolds; current small dt and under-relaxation mask potential instabilities but adding mild stabilization or more advanced solvers would improve efficiency and robustness.

Recommendation: Optional numerical improvements if you want faster/stabler convergence at Re=400: use a small amount of upwinding (or flux-limiter) for convective term, increase under-relaxation if instability observed, or use a larger dt until residuals stop decreasing. Also consider replacing splu with a direct multigrid or AMG Poisson solver if planning larger grids.



	Current Stage [C/3]
1. MINOR

Category: accuracy
Issue: Convergence criteria and acceptance logic are weak: the run ended by accepting the last iterate simply because the final velocity-change L2 was smaller than the initial one, even though neither vel-change nor divergence reached the stated tolerances.

Recommendation: Tighten stopping criteria or make acceptance logic more diagnostic: require a relative reduction (e.g. factor 1e-3) or monotonic decrease over a window of iterations, monitor both velocity-change and divergence trends, and optionally allow adaptive under-relaxation or increase max_iters if trends are still improving. Report a clear flag when the solution is only 'partially converged'.


2. MINOR

Category: structure
Issue: Progress logging prints at iteration 1 in addition to every 10% of max_iters, which violates the stated 'ONLY every 10%' requirement.

Recommendation: Restrict printing strictly to the schedule (report_every) unless printing the first/last iteration is explicitly desired. If keeping first/last prints is useful, document that exception in comments.


3. MINOR

Category: accuracy
Issue: Pressure Poisson boundary treatment uses simple mirroring (zero normal derivative) on all walls. This is a reasonable default but not strictly consistent with the projection method's Neumann BCs derived from the momentum equations and the imposed Dirichlet velocity BCs.

Recommendation: Improve physical consistency by assembling the Poisson RHS including boundary flux terms (phi Neumann boundary values computed from normal components of u_star and prescribed wall velocities) or apply the correct discrete normal-gradient conditions. This typically reduces residual divergence and yields a cleaner pressure field.


4. MINOR

Category: structure
Issue: Several array/matrix assembly routines and finite-difference operations are implemented with Python loops which is simple and clear but could be a performance bottleneck for larger grids.

Recommendation: Consider vectorizing derivative computations (numpy slicing) and using sparse matrix assembly utilities (e.g. kronecker products or scipy.sparse.diags) to speed setup. Keep the current, readable implementation for clarity but profile and replace hotspots if you scale to finer meshes.


5. MINOR

Category: runtime
Issue: Post-processing does not close the matplotlib figure after saving; in long batch runs this can accumulate memory. Also the single saved figure name is hard-coded.

Recommendation: After saving call plt.close(fig) to release memory. Consider parameterizing output filenames and adding a short header comment in the saved figure (or separate text file) describing grid and solver parameters (Nx, Ny, Re, dt, omega, iterations) for reproducibility.





### iteration_count
3

### rev_image_description
	Current Stage [A/3]


	Current Stage [B/3]


	Current Stage [C/3]




