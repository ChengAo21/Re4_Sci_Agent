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
We are to solve the 2D unsteady incompressible Navier–Stokes equations on the rectangular domain Omega = [0,2] x [0,1] with a body force f = (0, -sin(pi x) sin(pi y) sin(pi t)) and Reynolds number Re = 100 (visc = 1/Re). The unknowns are the velocity vector (u,v) and pressure p, subject to no-slip on top and bottom walls, a time-dependent Dirichlet inlet profile at x=0 for u (v=0), homogeneous Neumann velocity at the outlet x=2 (du/dx = dv/dx = 0), and Dirichlet pressure p=0 at the outlet. Initial velocities are zero and simulation runs until t = 1.0.

Primary mathematical and numerical challenges:
- Incompressibility constraint: the divergence-free condition couples u, v and p elliptically. A suitable algorithm must enforce divergence-free condition every time-step (projection, pressure correction, or saddle-point solve). Projection methods on staggered (MAC) grids are widely used for robustness and conservation.  
- Pressure-velocity coupling and boundary conditions: consistent imposition of the outlet pressure Dirichlet and velocity Neumann BCs on a staggered grid is delicate; the Poisson problem must be posed to respect those BCs and maintain solvability.  
- Advection nonlinearity: convective terms are nonlinear and can create instabilities or spurious oscillations if discretized improperly; using an upwind flux (donor-cell) or flux-limited schemes stabilizes advection at moderate Reynolds numbers.  
- Time-step selection: explicit treatment of convection imposes a CFL-type condition (dt <= O(dx/|u|)), while implicit diffusion relaxes the viscous constraint (dt <= O(dx^2/visc)). A semi-implicit IMEX strategy or fully implicit treatment of diffusion is recommended to allow reasonable dt.  
- Solver bottleneck: the pressure Poisson (elliptic) linear system dominates the cost per time-step. Efficient solvers (geometric multigrid or AMG) and reuse of factorizations/preconditioners are essential for reasonable runtimes.  
- Accuracy vs stability: tradeoffs between upwind (stable, lower order) and central/higher-order schemes (accurate but less stable) must be balanced. Using second-order temporal schemes (AB2+CN) reduces time discretization error, but implementation is more complex.  
- Implementation robustness: the staggered grid indexing, ghost cell treatment for boundary conditions, flattening order for linear solvers, and explicit declaration of global parameters (visc) all have to be handled carefully to avoid subtle bugs and shape mismatches.  

We should therefore implement a Chorin/projection (MAC) solver with upwind advection and semi-implicit diffusion or a second-order incremental projection (Van Kan) with AB2+Crank–Nicolson, accelerate the Poisson solves by multigrid or preconditioned Krylov, reuse matrix factorizations, and strictly follow the provided coding hints (use visc global, pass constants explicitly, use indexing='ij' and order='F', use atol in SciPy solvers, and print progress only at 10% increments). The final output must save a figure with contours of u, v and p at t=1.0 using the 'RdBu_r' colormap without calling plt.show().

### solution_plans
	Current Stage [A/3]
solu_name='Chorin (classical) projection on a MAC grid with semi‑implicit diffusion and upwind advection' content="Governing idea:\n1) Use a staggered (MAC) grid: u at vertical faces, v at horizontal faces, p at cell centers. 2) March in time with a fractional step: (a) compute a tentative velocity u* by advancing momentum without the pressure term (explicit convection + implicit diffusion + explicit forcing); (b) solve a Poisson equation for the pressure increment (phi) to enforce divergence-free velocity; (c) project u* to obtain u^{n+1} and update pressure. 3) Use upwind (donor‑cell) discretization for convective fluxes for stability and second‑order central for diffusion/pressure gradient. 4) Use sparse, reusable linear solvers for the Helmholtz (viscous) solves and an efficient Poisson solver (multigrid or preconditioned Krylov).\n\nAlgorithmic steps (pseudocode / numbered steps):\n1. Define global constants (visc = 1.0/Re) and pass all other constants explicitly into functions. 2. Build MAC grid with indexing='ij'. Place arrays with Fortran ordering when flattening (order='F'). 3. Discretize operators on MAC grid: divergence (cell centers), gradient (staggered), Laplacian for u/v (face-centered) and for p (cell-centered). Use second‑order central differences for diffusion and pressure gradient; use first‑order upwind for convective fluxes to ensure monotonicity at Re=100. 4. Choose dt by CFL: dt = CFL*min(dx/umax, dy/vmax, 0.5*min(dx^2,dy^2)/visc); set initial small dt and optionally adapt (bounded by accuracy). 5. Preassemble sparse matrices that do not change in time:\n   - Helmholtz matrix A = I - dt*visc*L (for u and v solves) if using implicit diffusion (or I + dt*visc*L depending on discretization sign convention).\n   - Poisson matrix for cell-centered Laplacian (for pressure correction). 6. Pre-factorize (or create preconditioners) for the time-invariant matrices where possible (e.g., sparse LU via splu or ILU + GMRES). Always call SciPy solvers with atol (not tol). 7. Time loop for n = 0..Nsteps up to t_final:\n   a. Compute convective terms C_x(u^n,v^n), C_y(u^n,v^n) using upwind fluxes on MAC grid.\n   b. Form RHS for tentative velocity: rhs_u = u^n + dt*( - C_x + f_x ); similarly for v. If using semi-implicit diffusion, move viscous term to LHS.\n   c. Solve Helmholtz systems A u* = rhs_u and A v* = rhs_v using the pre-factorized solver or an iterative method with atol. Ensure correct boundary conditions on faces (no-slip walls, inlet Dirichlet for u, homogeneous Neumann at outlet implemented by ghost/one-sided stencils). Reuse factorization for all time-steps to accelerate.\n   d. Compute divergence of tentative velocity div_u* at cell centers.\n   e. Solve Poisson: Lap phi = (1/dt) * div_u* with boundary conditions consistent with velocity BCs. Use a fast multigrid solver (pyamg / geometric multigrid) or preconditioned CG/BiCGStab with appropriate preconditioner. Use 'atol' for solver stopping. Enforce p( x=2 ) = 0 as Dirichlet on the appropriate cell-centered nodes (so the Poisson problem is well-posed).\n   f. Correct velocities: u^{n+1} = u* - dt * dphi/dx (apply at face locations), v^{n+1} = v* - dt * dphi/dy.\n   g. Update pressure: p^{n+1} = phi (incremental vs non-incremental variant: choose p^{n+1} = phi or p^{n+1} = p^n + phi depending on chosen projection variant and consistency with BCs).\n   h. Impose physical boundary conditions exactly on staggered grid after projection (inlet u, walls no-slip, outlet pressure Dirichlet p=0, outlet velocity Neumann consistency via ghost values). 8. Print concise progress every 10% of nsteps only. 9. At t_final save contours of u, v, p (interpolate from staggered to cell centers for plotting) in one figure using 'RdBu_r' colormap; save figure(s) without plt.show().\n\nStability / accuracy / complexity / efficiency limitations:\n- Temporal accuracy: classical Chorin projection with first-order splitting is only first-order accurate in time; can be improved by using incremental or higher-order projection variants (at extra implementation cost).  - Spatial accuracy: upwind advection is stable but only first-order accurate for advection; central advection would be higher-order but can produce oscillations at moderate Re.  - Stability: explicit treatment of convection imposes a CFL constraint on dt (dt <= O(min(dx/|u|,dy/|v|))). Semi-implicit diffusion relaxes viscous constraint.  - Linear solver cost: Poisson solve per time-step dominates computational cost; using a good multigrid or reusing factorizations is crucial.  - Implementation complexity: correct staggered BC implementation, Poisson BC consistency and constructing face/center stencils is delicate and error-prone.  - Memory: storing sparse factorizations or AMG hierarchies increases memory usage.  - Robustness: if AMG isn't available, iterative solvers may stagnate without good preconditioning; always use atol in SciPy calls and check convergence."

	Current Stage [B/3]
solu_name='Second‑order incremental projection (Van Kan) on MAC grid with AB2 convective treatment and Crank‑Nicolson diffusion + geometric multigrid Poisson' content="Governing idea:\n1) Improve temporal accuracy by using an incremental pressure‑correction (Van Kan) / second‑order projection: use Adams‑Bashforth 2 (AB2) for nonlinear convective terms (explicit) and Crank‑Nicolson (CN) for viscous term (implicit second order). 2) Solve two Helmholtz-like systems per velocity component for implicit diffusion at each step; then solve a cell-centered Poisson for the pressure increment and correct velocities. 3) Use a geometric multigrid Poisson solver tuned for the MAC Laplacian (fast O(N)) to accelerate the dominant solve. Use flux-limited higher-order upwind (MUSCL or limited linear reconstruction) if improved advection accuracy is desired.\n\nAlgorithmic steps (pseudocode / numbered steps):\n1. Define global constant visc = 1.0/Re and pass domain and solver parameters explicitly to all functions. Use indexing='ij' for grid generation and order='F' for flattening. 2. Build staggered MAC grid and assemble discrete operators as in Plan 1; precompute constant sparse matrices for CN Helmholtz solves and for the Poisson operator (cell-centered Laplacian). 3. Initialize velocity fields u^0, v^0 = 0 and compute a first sub-step using e.g., forward Euler (or Crank–Nicolson with AB1) to seed AB2. 4. Time stepping (for each n >= 1):\n   a. Compute convective terms C^n and C^{n-1} using upwind or a limited reconstruction; form AB2 predictor for convection: C_AB2 = 1.5*C^n - 0.5*C^{n-1}.\n   b. Form left-hand side for CN: (I - 0.5*dt*visc*L) u_tent = (I + 0.5*dt*visc*L) u^n + dt*( - C_AB2 + f^{n+1/2} )  (and similarly for v). This produces Helmholtz-like systems with the same LHS each step; pre-factorize LHS or build a multigrid/AMG preconditioner.\n   c. Solve for tentative velocities u*, v* with the chosen linear solver (reuse factorization or preconditioner). Always pass atol and appropriate solver options.  \n   d. Compute divergence div(u*) and solve the Poisson problem for pressure increment phi: Lap phi = (1/dt) div(u*). Use a geometric multigrid solver specialized to the MAC/cell layout (fast and scalable). Enforce p( x=2 ) = 0 as the Dirichlet reference. \n   e. Correct velocities: u^{n+1} = u* - dt * dphi/dx; v^{n+1} = v* - dt * dphi/dy. Update pressure incrementally: p^{n+1} = p^n + phi (this is the incremental variant improving pressure accuracy). \n   f. Enforce boundary conditions explicitly on staggered velocities and pressure after correction. 5. Acceleration strategies: reuse factorization/ILU of Helmholtz LHS (constant coefficients); use geometric multigrid for Poisson (optimal O(N) work); vectorize flux computations; optionally use coarse tolerance for intermediate iterative solves and tighten tolerances near final time or when reporting. Print concise progress every 10% of total steps.\n\nStability / accuracy / complexity / efficiency limitations:\n- Temporal accuracy: scheme is formally second order in time (AB2+CN) for velocity, and incremental projection improves pressure accuracy; projection splitting still introduces splitting error but smaller than first-order projection.  - Spatial accuracy: central diffusion/pressure gradients give second-order; advection with AB2 + upwind is at least second-order in time but spatial upwind reconstruction may limit spatial order; using higher-order limited reconstructions raises complexity.  - Stability: AB2 for advection is conditionally stable — dt must satisfy convective CFL (less restrictive than explicit Euler but still limited). CN implicit diffusion removes strong viscous restriction.  - Implementation complexity: building a robust geometric multigrid Poisson solver and ensuring correct staggered BCs and coarse-grid transfers is more complex than using black-box solvers.  - Computational cost: each time step involves two Helmholtz solves (u and v) plus one Poisson; multigrid and reuse of factorizations mitigate cost.  - Robustness: higher-order temporal discretization can be more sensitive to boundary-condition inconsistency; careful initialization and consistent treatment of inlet/outlet BCs are required.  - Practical caveat: ensure all SciPy solver calls use atol, use indexing='ij' and order='F' consistently, and define visc globally; print only every 10% of progress to satisfy output requirements."

	Current Stage [C/3]
solu_name='Implementation and practical notes (common to both strategies)' content="Governing idea and practical enforcement:\n1) Always declare visc = 1.0 / Re at global scope and use visc everywhere in numerical formulae. 2) Use a staggered MAC grid and be meticulous about where each discrete variable lives (u: Nx+1 by Ny, v: Nx by Ny+1, p: Nx by Ny). 3) Use indexing='ij' when creating coordinate arrays with numpy.meshgrid, and use order='F' when flattening arrays before passing them to sparse solvers. 4) Use upwind differencing for convection on faces to guarantee stability at Re=100 and to avoid spurious oscillations. 5) For Poisson solves prefer a geometric multigrid implementation or pyamg if available; otherwise use preconditioned Krylov solvers (CG/BiCGStab) with ILU or AMG preconditioners. Always pass atol rather than tol to SciPy iterative solvers. 6) Reuse sparse factorizations or preconditioners for time‑invariant matrices to accelerate the time loop. 7) Compute time step dt via a CFL formula and optionally adapt it but cap the change per step to avoid instability. 8) Interpolate staggered u,v to cell centers for plotting pcolor/contour and create a single figure with subplots of u, v, p using cmap='RdBu_r'; call plt.savefig(...) and do not call plt.show(). 9) Print concise progress only at 10%, 20%, ... of total steps. 10) Test solver on coarse grids first and verify mass conservation (L2 norm of divergence should be near machine tolerance after projection) and monitor residuals of linear solves. \n\nLimitations / warnings:\n- The main implementation difficulties are correct enforcement of staggered BCs and Poisson boundary conditions (consistency with velocity BC). - Using higher-order advection schemes or unsplit Navier–Stokes solvers can give better accuracy but are more complex to implement robustly. - Ensure solver tolerances are chosen to balance accuracy and compute time; using too loose tolerances yields mass conservation errors, too strict tolerances waste CPU. - Follow the provided hints: never use tol keyword in SciPy solvers, do not use np.trapz, use indexing='ij' and order='F', and explicitly pass constants to avoid NameError."



### technical_spec
	Current Stage [A/3]
This script implements a Chorin projection solver on a staggered (MAC) grid for 2D unsteady incompressible Navier–Stokes with explicit body force. It uses first-order upwind-like advection approximations (sign-based upwind differences), central diffusion and pressure gradient discretization, and a cell-centered Poisson solver (sparse direct factorization) assembled once and reused. Velocity boundary conditions (no-slip walls, time-dependent inlet Dirichlet, outlet Neumann) are enforced explicitly each time-step; pressure uses Neumann everywhere except a Dirichlet reference at the outlet. The solver is explicit in time for convective and viscous terms (stable with chosen dt), and uses projection to enforce divergence-free velocities. Final fields (u,v,p) are interpolated to cell centers and saved as a figure with three contours using the RdBu_r colormap. Progress prints are emitted only at 10% intervals.


	Current Stage [B/3]
This script implements a Chorin projection solver on a staggered (MAC) grid for 2D unsteady incompressible Navier–Stokes with explicit body force. We explicitly assemble the discrete divergence (D) and gradient (G) operators on the MAC grid so the pressure Poisson matrix A = D*G is algebraically consistent with the velocity update and divergence operator. Dirichlet pressure at the outlet (x=2) is handled by moving the known ghost-face contributions into a constant RHS vector (g_const). The Poisson matrix A is assembled once and factorized (direct solver). Convective terms and face interpolations are vectorized using numpy array slicing. Boundary conditions are applied explicitly each step. The velocity correction uses the same G operator (and g_const) that was used to create A which yields a projection that enforces discrete divergence to solver tolerance. Quantitative outputs printed: progress with divergence L2 every 10% and final u/v/p min/max and final divergence L2. A three-panel contour figure (u,v,p at cell centers) is saved to disk as ns_mac_u_v_p.png using RdBu_r colormap.

	Current Stage [C/3]
This script implements a Chorin projection solver for 2D unsteady incompressible Navier–Stokes on a MAC (staggered) grid. Core data flow:
- Grid and MAC indexing are created once.
- Discrete divergence D and gradient G operators are assembled once and Poisson matrix A = D*G is factorized.
- Each time step: apply consistent BCs, compute convective terms (upwind, elementwise decisions), diffusion (Laplacian), explicit forcing, form tentative velocities, project via Poisson solve to enforce discrete divergence, correct velocities, and update pressure.
Key design choices to ensure consistency and accuracy:
- All flatten/reshape use Fortran ordering (order='F') via helpers flatten_F/reshape_F.
- Boundary conditions are centralized in apply_bcs and all differential stencils and interpolations use the same BC assumptions (walls Dirichlet u/v=0, inlet Dirichlet u-profile & v=0, outlet Neumann copying interior faces).
- Convective upwind uses elementwise sign tests (no .any()), and boundary one-sided stencils use the same ghost values implied by apply_bcs (explicitly using known zeros or copied values).
- Operator consistency tests (D*(G*p + g_const) == A*p + D*g_const) and divergence-reduction unit tests are performed once after operator assembly to detect algebraic inconsistencies.
- Poisson matrix A is assembled once and factorized (direct LU) for efficiency.
Diagnostics:
- Progress printed every 10% with both L2 and Linf divergence norms.
- After operator assembly: printed algebraic consistency residual and divergence-reduction test norm.
- Final printed minima/maxima for u,v,p and final divergence norms.
Plots:
- u, v, p interpolated to cell centers and saved as ns_mac_u_v_p.png using RdBu_r colormap (no plt.show()).



### python_codes

#### Script block1:

```python
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# Global viscosity defined as visc = 1.0 / Re
Re = 100.0
visc = 1.0 / Re

# -------------------- Utility / Discretization functions --------------------

def build_poisson_matrix(Nx, Ny, dx, dy, bc_right='dirichlet'):
    """
    Build cell-centered Poisson matrix with Neumann on left/top/bottom and
    Dirichlet or Neumann on right depending on bc_right.
    Uses finite difference: for interior neighbors adds +1/dx2 off-diagonal and
    diagonal = -sum_off_diagonals. For Dirichlet neighbor (known value g) the
    neighbor contributes to diag and to RHS (here g=0 for outlet), for Neumann
    missing neighbor contributes nothing.
    Returns sparse matrix A (Np x Np) and a mask of Dirichlet nodes (False here)
    """
    Np = Nx * Ny
    dx2 = dx * dx
    dy2 = dy * dy
    rows = []
    cols = []
    data = []

    def idx(i, j):
        return i + j * Nx  # Fortran order flattening (i fastest)

    for j in range(Ny):
        for i in range(Nx):
            k = idx(i, j)
            sum_coef = 0.0
            # left neighbor
            if i - 1 >= 0:
                kr = idx(i - 1, j)
                rows.append(k); cols.append(kr); data.append(1.0 / dx2)
                sum_coef += 1.0 / dx2
            else:
                # left boundary: Neumann -> do nothing
                pass
            # right neighbor
            if i + 1 < Nx:
                kr = idx(i + 1, j)
                rows.append(k); cols.append(kr); data.append(1.0 / dx2)
                sum_coef += 1.0 / dx2
            else:
                # at right external boundary: either Dirichlet (p=0) or Neumann
                if bc_right == 'dirichlet':
                    # known boundary value g=0 -> contributes to diag but RHS remains unchanged
                    sum_coef += 1.0 / dx2
                else:
                    # Neumann: do nothing
                    pass
            # down neighbor
            if j - 1 >= 0:
                kr = idx(i, j - 1)
                rows.append(k); cols.append(kr); data.append(1.0 / dy2)
                sum_coef += 1.0 / dy2
            else:
                # bottom boundary: Neumann -> do nothing
                pass
            # up neighbor
            if j + 1 < Ny:
                kr = idx(i, j + 1)
                rows.append(k); cols.append(kr); data.append(1.0 / dy2)
                sum_coef += 1.0 / dy2
            else:
                # top boundary: Neumann -> do nothing
                pass
            # diagonal
            rows.append(k); cols.append(k); data.append(-sum_coef)
    A = sp.csr_matrix((data, (rows, cols)), shape=(Np, Np))
    return A


def flatten_F(arr):
    return arr.ravel(order='F')


def inv_idx(k, Nx):
    # returns (i,j) from flattened Fortran-order index
    i = k % Nx
    j = k // Nx
    return i, j

# -------------------- Grid and problem setup --------------------

def create_mac_grid(Lx, Ly, Nx, Ny):
    dx = Lx / Nx
    dy = Ly / Ny
    # u at vertical faces: shape (Nx+1, Ny), x=0..Nx, y centers
    u_shape = (Nx + 1, Ny)
    # v at horizontal faces: shape (Nx, Ny+1), x centers, y=0..Ny
    v_shape = (Nx, Ny + 1)
    # p at cell centers: shape (Nx, Ny)
    p_shape = (Nx, Ny)

    # cell center coordinates
    xc = (np.arange(Nx) + 0.5) * dx
    yc = (np.arange(Ny) + 0.5) * dy
    Xc, Yc = np.meshgrid(xc, yc, indexing='ij')  # shape (Nx,Ny)
    return dx, dy, Xc, Yc, u_shape, v_shape, p_shape

# -------------------- Boundary conditions and forcing --------------------

def inlet_u_profile(y, t):
    return np.sin(np.pi * y) * (np.sin(np.pi * t) + np.sin(3 * np.pi * t) + np.sin(5 * np.pi * t))


def apply_bcs(u, v, dx, dy, Xc, Yc, t):
    # u shape (Nx+1, Ny), v shape (Nx, Ny+1)
    Nx_plus1, Ny = u.shape
    Nx, Ny_plus1 = v.shape
    # top/bottom no-slip: set u adjacent rows and v boundary rows to zero
    # For u, indices j=0..Ny-1 correspond to y centers; enforce u[:,0]=0 and u[:,-1]=0 to approximate no-slip
    u[:, 0] = 0.0
    u[:, -1] = 0.0
    # For v, j=0 and j=Ny are exactly on bottom and top walls
    v[:, 0] = 0.0
    v[:, -1] = 0.0

    # inlet left x=0: u(0,y,t) is prescribed at vertical face i=0
    # Evaluate at u face y positions which are at yc = (j+0.5)*dy
    y_u = (np.arange(Ny) + 0.5) * dy
    u[0, :] = inlet_u_profile(y_u, t)
    # inlet v at left boundary: enforce zero normal component at boundary approximately
    # set v[0,:] = 0
    v[0, :] = 0.0

    # outlet right x=Lx: homogeneous Neumann for velocities -> du/dx = 0, dv/dx = 0
    # approximate by copying interior values to outer face: u[-1,:] = u[-2,:]
    u[-1, :] = u[-2, :].copy()
    # For v, rightmost v index is i=Nx-1 (last interior center), enforce v[-1,:] = v[-2,:]
    v[-1, :] = v[-2, :].copy()

    return u, v


def compute_divergence(u, v, dx, dy):
    # u shape (Nx+1, Ny), v shape (Nx, Ny+1)
    Nx_plus1, Ny = u.shape
    Nx, Ny_plus1 = v.shape
    Nx_cells = Nx
    Ny_cells = Ny
    div = np.zeros((Nx_cells, Ny_cells))
    # div[i,j] = (u[i+1,j]-u[i,j])/dx + (v[i,j+1]-v[i,j])/dy
    div[:, :] = (u[1:, :] - u[:-1, :]) / dx + (v[:, 1:] - v[:, :-1]) / dy
    return div

# -------------------- Differential operators --------------------

def lap_u(u, dx, dy):
    # u shape (Nx+1, Ny)
    lap = np.zeros_like(u)
    # interior in i: 1..Nx-1, j: 1..Ny-2 (since u has Ny points in y)
    Nx_plus1, Ny = u.shape
    # x second derivative
    lap[1:-1, :] += (u[2:, :] - 2.0 * u[1:-1, :] + u[:-2, :]) / (dx * dx)
    # y second derivative (handle j indices 0..Ny-1 -> interior j=1..Ny-2)
    lap[:, 1:-1] += (u[:, 2:] - 2.0 * u[:, 1:-1] + u[:, :-2]) / (dy * dy)
    # For j=0 and j=Ny-1 (adjacent to walls), assume Dirichlet u=0 at wall -> u[:, -] already set to 0 by BC
    lap[:, 0] += (u[:, 1] - 2.0 * u[:, 0] + 0.0) / (dy * dy)
    lap[:, -1] += (0.0 - 2.0 * u[:, -1] + u[:, -2]) / (dy * dy)
    return lap


def lap_v(v, dx, dy):
    lap = np.zeros_like(v)
    Nx, Ny_plus1 = v.shape
    lap[1:-1, :] += (v[2:, :] - 2.0 * v[1:-1, :] + v[:-2, :]) / (dx * dx)
    lap[:, 1:-1] += (v[:, 2:] - 2.0 * v[:, 1:-1] + v[:, :-2]) / (dy * dy)
    # For i edges (left/right) where Neumann is used for v at rightmost and v[0,:] set to 0 at inlet
    lap[0, :] += (v[1, :] - 2.0 * v[0, :] + 0.0) / (dx * dx)
    lap[-1, :] += (0.0 - 2.0 * v[-1, :] + v[-2, :]) / (dx * dx)
    return lap

# -------------------- Convective terms (upwind-ish) --------------------

def compute_convective_u(u, v, dx, dy):
    # returns conv term u * du/dx + v_on_u * du/dy at u-locations (shape u)
    Nx_plus1, Ny = u.shape
    conv = np.zeros_like(u)
    # compute du/dx using upwind based on u sign
    # interior in i: 1..Nx-1
    # du/dx
    du_dx = np.zeros_like(u)
    # for i=1..Nx-1
    pos = u >= 0.0
    # backward difference where positive
    du_dx[1:-1, :] = np.where(pos[1:-1, :], (u[1:-1, :] - u[:-2, :]) / dx, (u[2:, :] - u[1:-1, :]) / dx)
    # boundaries: left face (i=0) and right face (i=Nx)
    du_dx[0, :] = (u[0, :] - u[1, :]) / dx  # approximate
    du_dx[-1, :] = (u[-1, :] - u[-2, :]) / dx

    # interpolate v to u locations (bilinear average): v has shape (Nx, Ny+1)
    v_on_u = np.zeros_like(u)
    # for interior u at i=1..Nx-1, j=0..Ny-1
    for i in range(1, Nx_plus1 - 1):
        for j in range(0, Ny):
            # v indices to average: v[i-1,j] v[i-1,j+1] v[i,j] v[i,j+1]
            v00 = v[i - 1, j]
            v01 = v[i - 1, j + 1]
            v10 = v[i, j]
            v11 = v[i, j + 1]
            v_on_u[i, j] = 0.25 * (v00 + v01 + v10 + v11)
    # edges: i=0
    i = 0
    for j in range(0, Ny):
        v_on_u[i, j] = 0.5 * (v[0, j] + v[0, j + 1])
    # i = Nx
    i = Nx_plus1 - 1
    for j in range(0, Ny):
        v_on_u[i, j] = 0.5 * (v[-1, j] + v[-1, j + 1])

    # du/dy using upwind based on v_on_u
    du_dy = np.zeros_like(u)
    pos_v = v_on_u >= 0.0
    # interior j indices 1..Ny-2
    du_dy[:, 1:-1] = np.where(pos_v[:, 1:-1], (u[:, 1:-1] - u[:, :-2]) / dy, (u[:, 2:] - u[:, 1:-1]) / dy)
    # boundaries in y
    du_dy[:, 0] = (u[:, 0] - u[:, 1]) / dy
    du_dy[:, -1] = (u[:, -1] - u[:, -2]) / dy

    conv = u * du_dx + v_on_u * du_dy
    return conv


def compute_convective_v(u, v, dx, dy):
    # returns conv term u_on_v * dv/dx + v * dv/dy at v-locations (shape v)
    Nx, Ny_plus1 = v.shape
    conv = np.zeros_like(v)
    # dv/dy using upwind based on v sign
    dv_dy = np.zeros_like(v)
    pos = v >= 0.0
    dv_dy[:, 1:-1] = np.where(pos[:, 1:-1], (v[:, 1:-1] - v[:, :-2]) / dy, (v[:, 2:] - v[:, 1:-1]) / dy)
    dv_dy[:, 0] = (v[:, 0] - v[:, 1]) / dy
    dv_dy[:, -1] = (v[:, -1] - v[:, -2]) / dy

    # interpolate u to v locations
    u_on_v = np.zeros_like(v)
    for i in range(0, Nx):
        for j in range(1, Ny_plus1 - 1):
            # u indices around v: u[i, j-1], u[i+1, j-1], u[i, j], u[i+1, j]
            u00 = u[i, j - 1]
            u10 = u[i + 1, j - 1]
            u01 = u[i, j]
            u11 = u[i + 1, j]
            u_on_v[i, j] = 0.25 * (u00 + u10 + u01 + u11)
    # j=0 (bottom boundary)
    j = 0
    for i in range(0, Nx):
        u_on_v[i, j] = 0.5 * (u[i, 0] + u[i + 1, 0])
    # j=Ny (top boundary)
    j = Ny_plus1 - 1
    for i in range(0, Nx):
        u_on_v[i, j] = 0.5 * (u[i, -1] + u[i + 1, -1])

    # dv/dx using upwind based on u_on_v
    dv_dx = np.zeros_like(v)
    pos_u = u_on_v >= 0.0
    dv_dx[1:-1, :] = np.where(pos_u[1:-1, :], (v[1:-1, :] - v[:-2, :]) / dx, (v[2:, :] - v[1:-1, :]) / dx)
    dv_dx[0, :] = (v[0, :] - v[1, :]) / dx
    dv_dx[-1, :] = (v[-1, :] - v[-2, :]) / dx

    conv = u_on_v * dv_dx + v * dv_dy
    return conv

# -------------------- Main solver --------------------

def run_simulation(Lx=2.0, Ly=1.0, Nx=64, Ny=32, t_final=1.0):
    dx, dy, Xc, Yc, u_shape, v_shape, p_shape = create_mac_grid(Lx, Ly, Nx, Ny)

    # Initialize fields
    u = np.zeros(u_shape)  # (Nx+1, Ny)
    v = np.zeros(v_shape)  # (Nx, Ny+1)
    p = np.zeros(p_shape)  # (Nx, Ny)

    # Poisson matrix assembly (cell-centered)
    A = build_poisson_matrix(Nx, Ny, dx, dy, bc_right='dirichlet')
    A_csc = A.tocsc()
    # Factorize once (direct solver) for repeated solves
    lu = spla.splu(A_csc)

    # Time stepping parameters (explicit conv + explicit visc used here)
    # Choose dt by CFL: convection and diffusion
    umax_est = 3.0  # rough max inlet amplitude
    cfl = 0.2
    dt_conv = cfl * min(dx, dy) / max(umax_est, 1e-6)
    dt_diff = 0.25 * min(dx * dx, dy * dy) / visc
    dt = min(dt_conv, dt_diff, 0.002)
    Nsteps = int(np.ceil(t_final / dt))
    dt = t_final / Nsteps  # adjust dt to fit exactly

    print(f"Grid: Nx={Nx}, Ny={Ny}, dx={dx:.4e}, dy={dy:.4e}")
    print(f"Time stepping: dt={dt:.6e}, Nsteps={Nsteps}")

    # Progress printing every 10%
    progress_step = max(1, Nsteps // 10)

    t = 0.0
    for n in range(1, Nsteps + 1):
        t_new = t + dt
        # Enforce BCs before computing conv/diff
        u, v = apply_bcs(u, v, dx, dy, Xc, Yc, t)

        # Compute convective terms
        conv_u = compute_convective_u(u, v, dx, dy)
        conv_v = compute_convective_v(u, v, dx, dy)

        # Diffusion terms via laplacian
        lapu = lap_u(u, dx, dy)
        lapv = lap_v(v, dx, dy)

        # External forcing (only fy nonzero)
        # f_x = 0, f_y = -sin(pi x) sin(pi y) sin(pi t)
        # Evaluate forcing at face centers for u/v properly: for v we need f_y at v locations
        # Approximate f_y at v-locations by sampling at v-location coords
        # create coords for v location
        xv = (np.arange(Nx) + 0.5) * dx
        yv = np.arange(Ny + 1) * dy
        XV, YV = np.meshgrid(xv, yv, indexing='ij')  # careful indexing
        f_x = np.zeros_like(u)
        f_y = -np.sin(np.pi * XV) * np.sin(np.pi * YV) * np.sin(np.pi * t_new)
        # f_y shape (Nx, Ny+1) matches v

        # Tentative velocities (explicit time stepping)
        u_star = u + dt * (-conv_u + visc * lapu + f_x)
        v_star = v + dt * (-conv_v + visc * lapv + f_y)

        # Enforce BCs on tentative as well
        u_star, v_star = apply_bcs(u_star, v_star, dx, dy, Xc, Yc, t_new)

        # Compute divergence of tentative velocity at cell centers
        div_star = compute_divergence(u_star, v_star, dx, dy)

        # Solve Poisson: Lap(phi) = (1/dt) * div_star
        rhs = (flatten_F(div_star) / dt)
        # Adjust RHS for Dirichlet neighbors contribution: at right boundary Dirichlet p=0, neighbor contributes nothing
        # Solve A phi = rhs
        phi_flat = lu.solve(rhs)
        phi = phi_flat.reshape((Nx, Ny), order='F')

        # Correct velocities: u_{i} = u_star - dt * dphi/dx at u faces (i=1..Nx-1)
        # dp/dx at u faces interior i=1..Nx-1
        dpdx = np.zeros_like(u)
        dpdx[1:-1, :] = (phi[1:, :] - phi[:-1, :]) / dx
        # For faces at i=0 (inlet) and i=Nx (outlet) don't correct since Dirichlet/Neumann applied

        dpdy = np.zeros_like(v)
        dpdy[:, 1:-1] = (phi[:, 1:] - phi[:, :-1]) / dy

        u = u_star - dt * dpdx
        v = v_star - dt * dpdy

        # Enforce BCs after projection
        u, v = apply_bcs(u, v, dx, dy, Xc, Yc, t_new)

        # Update pressure (use phi as pressure approximation)
        p = phi.copy()

        t = t_new

        # Print progress only at 10% intervals
        if (n % progress_step) == 0 or n == Nsteps:
            div_after = compute_divergence(u, v, dx, dy)
            div_L2 = np.linalg.norm(div_after.ravel(order='F')) * np.sqrt(dx * dy)
            print(f"Step {n}/{Nsteps}, t={t:.4f}, ||div||_2={div_L2:.3e}")

    # Interpolate u,v to cell centers for plotting
    u_center = 0.5 * (u[:-1, :] + u[1:, :])
    v_center = 0.5 * (v[:, :-1] + v[:, 1:])

    # Save figure with three contours using RdBu_r colormap
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    levels = 50
    cf0 = axes[0].contourf(Xc, Yc, u_center, levels=levels, cmap='RdBu_r')
    axes[0].set_title('u at t=%.3f' % t)
    axes[0].set_xlabel('x'); axes[0].set_ylabel('y')
    fig.colorbar(cf0, ax=axes[0])

    cf1 = axes[1].contourf(Xc, Yc, v_center, levels=levels, cmap='RdBu_r')
    axes[1].set_title('v at t=%.3f' % t)
    axes[1].set_xlabel('x'); axes[1].set_ylabel('y')
    fig.colorbar(cf1, ax=axes[1])

    cf2 = axes[2].contourf(Xc, Yc, p, levels=levels, cmap='RdBu_r')
    axes[2].set_title('p at t=%.3f' % t)
    axes[2].set_xlabel('x'); axes[2].set_ylabel('y')
    fig.colorbar(cf2, ax=axes[2])

    plt.tight_layout()
    plt.savefig('ns_mac_u_v_p.png', dpi=200)

    # Print concise summary metrics
    print('Final statistics:')
    print(f"u: min={u_center.min():.4e}, max={u_center.max():.4e}")
    print(f"v: min={v_center.min():.4e}, max={v_center.max():.4e}")
    print(f"p: min={p.min():.4e}, max={p.max():.4e}")
    div_final = compute_divergence(u, v, dx, dy)
    div_final_norm = np.linalg.norm(div_final.ravel(order='F')) * np.sqrt(dx * dy)
    print(f"Final divergence L2 norm: {div_final_norm:.3e}")

    return u, v, p, Xc, Yc

# -------------------- Run if main --------------------

if __name__ == '__main__':
    # Keep grid modest to ensure reasonable run-time
    u, v, p, Xc, Yc = run_simulation(Nx=64, Ny=32, t_final=1.0)
```


#### Script block2:

```python
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# Global viscosity defined as visc = 1.0 / Re
Re = 100.0
visc = 1.0 / Re

# Utility flatten/reshape helpers use Fortran order everywhere
def flatten_F(arr):
    return arr.ravel(order='F')

# -------------------- Grid and problem setup --------------------

def create_mac_grid(Lx, Ly, Nx, Ny):
    dx = Lx / Nx
    dy = Ly / Ny
    # u at vertical faces: shape (Nx+1, Ny)
    u_shape = (Nx + 1, Ny)
    # v at horizontal faces: shape (Nx, Ny+1)
    v_shape = (Nx, Ny + 1)
    # p at cell centers: shape (Nx, Ny)
    p_shape = (Nx, Ny)

    # cell center coordinates (Fortran ordering for arrays)
    xc = (np.arange(Nx) + 0.5) * dx
    yc = (np.arange(Ny) + 0.5) * dy
    Xc, Yc = np.meshgrid(xc, yc, indexing='ij')  # shape (Nx,Ny)
    return dx, dy, Xc, Yc, u_shape, v_shape, p_shape

# -------------------- Operators assembly (MAC-consistent) --------------------

def build_mac_operators(Nx, Ny, dx, dy, bc_right_dirichlet_value=0.0):
    """
    Build discrete divergence D (Np x Nfaces) and gradient G (Nfaces x Np) operators
    consistent with MAC staggering. Faces ordering: u-faces first (Fortran flattening
    of u array), then v-faces. Returns A = D*G (cell-centered Poisson matrix),
    sparse D, G, and g_const (length Nfaces) containing contributions from Dirichlet
    ghost-face pressures (so that face gradients = G*p + g_const).

    Boundary treatment:
    - Pressure Dirichlet applied at outlet (x = Lx) is represented by ghost face
      pressure equal to bc_right_dirichlet_value and its contribution moved into g_const.
    - Other external faces assume Neumann zero-gradient for pressure (implemented by
      treating missing neighbor as equal to adjacent interior cell -> no G contribution).
    """
    Np = Nx * Ny
    Nuf = (Nx + 1) * Ny
    Nvf = Nx * (Ny + 1)
    Nfaces = Nuf + Nvf

    # index helpers (Fortran order flattening)
    def idx_p(i, j):
        return i + j * Nx

    def idx_u(i, j):
        return i + j * (Nx + 1)

    def idx_v(i, j):
        return Nuf + (i + j * Nx)

    # Build D: maps face values (u and v) -> cell divergence
    rows_D = []
    cols_D = []
    data_D = []
    for j in range(Ny):
        for i in range(Nx):
            row = idx_p(i, j)
            # u contributions: (u_{i+1,j} - u_{i,j})/dx
            up = idx_u(i + 1, j)
            um = idx_u(i, j)
            rows_D.extend([row, row])
            cols_D.extend([up, um])
            data_D.extend([1.0 / dx, -1.0 / dx])
            # v contributions: (v_{i,j+1} - v_{i,j})/dy
            vp = idx_v(i, j + 1)
            vm = idx_v(i, j)
            rows_D.extend([row, row])
            cols_D.extend([vp, vm])
            data_D.extend([1.0 / dy, -1.0 / dy])

    D = sp.csr_matrix((data_D, (rows_D, cols_D)), shape=(Np, Nfaces))

    # Build G: maps cell pressures -> face gradients (dpdx at u faces, dpdy at v faces)
    rows_G = []
    cols_G = []
    data_G = []
    # g_const holds constant contributions from Dirichlet ghost pressures
    g_const = np.zeros(Nfaces)

    # u-faces
    for j in range(Ny):
        for i in range(Nx + 1):
            row = idx_u(i, j)
            # interior u-face between cell (i-1,j) and (i,j)
            if 1 <= i <= Nx - 1:
                # dpdx = (p_{i,j} - p_{i-1,j})/dx
                pR = idx_p(i, j)
                pL = idx_p(i - 1, j)
                rows_G.extend([row, row])
                cols_G.extend([pR, pL])
                data_G.extend([1.0 / dx, -1.0 / dx])
            elif i == 0:
                # left boundary face: assume Neumann (zero gradient): dpdx=0 -> no contribution
                pass
            else:  # i == Nx, rightmost face touching outlet boundary
                # dpdx = (p_ghost - p_{Nx-1,j})/dx
                # p_ghost = bc_right_dirichlet_value
                pL = idx_p(Nx - 1, j)
                rows_G.append(row)
                cols_G.append(pL)
                data_G.append(-1.0 / dx)
                g_const[row] = bc_right_dirichlet_value / dx

    # v-faces
    for j in range(Ny + 1):
        for i in range(Nx):
            row = idx_v(i, j)
            # interior v-face between cell (i,j-1) and (i,j)
            if 1 <= j <= Ny - 1:
                pT = idx_p(i, j)
                pB = idx_p(i, j - 1)
                rows_G.extend([row, row])
                cols_G.extend([pT, pB])
                data_G.extend([1.0 / dy, -1.0 / dy])
            else:
                # bottom (j==0) or top (j==Ny) boundaries: assume Neumann (zero gradient)
                pass

    G = sp.csr_matrix((data_G, (rows_G, cols_G)), shape=(Nfaces, Np))

    # Build A = D * G (cell-centered Poisson matrix consistent with D and G)
    A = (D @ G).tocsr()

    return D, G, A, g_const

# -------------------- Boundary conditions and forcing --------------------

def inlet_u_profile(y, t):
    return np.sin(np.pi * y) * (np.sin(np.pi * t) + np.sin(3 * np.pi * t) + np.sin(5 * np.pi * t))


def apply_bcs(u, v, dx, dy, t, Lx=2.0):
    """
    Apply velocity boundary conditions in-place on u (Nx+1,Ny) and v (Nx,Ny+1).
    - top/bottom walls (y=0,1): no-slip (u and v zero at walls)
    - inlet x=0: u face i=0 prescribed by inlet profile, v normal component zero at inlet
    - outlet x=Lx: homogeneous Neumann for velocities approximated by copying interior
      values into boundary face (i.e., du/dx=0 -> u[-1,:]=u[-2,:])

    Note: Xc/Yc removed since not needed. This function documents BC approximation.
    """
    Nx_plus1, Ny_local = u.shape
    Nx_local, Ny_plus1 = v.shape

    # Top/bottom walls no-slip
    # For u faces, j=0..Ny-1 correspond to centers in y; adjacent to walls: set u[:,0]=0 and u[:,-1]=0
    u[:, 0] = 0.0
    u[:, -1] = 0.0
    # For v faces, j=0 and j=Ny correspond to walls, set zero
    v[:, 0] = 0.0
    v[:, -1] = 0.0

    # Inlet at x=0: prescribe u at face i=0
    y_u = (np.arange(Ny_local) + 0.5) * dy
    u[0, :] = inlet_u_profile(y_u, t)
    # inlet normal v component approximate zero
    v[0, :] = 0.0

    # Outlet at x=Lx: homogeneous Neumann du/dx=dv/dx=0 approx by copying interior
    u[-1, :] = u[-2, :].copy()
    v[-1, :] = v[-2, :].copy()

    return u, v

# -------------------- Differential operators (vectorized) --------------------

def lap_u(u, dx, dy):
    lap = np.zeros_like(u)
    # x second derivative interior
    lap[1:-1, :] += (u[2:, :] - 2.0 * u[1:-1, :] + u[:-2, :]) / (dx * dx)
    # y second derivative interior
    lap[:, 1:-1] += (u[:, 2:] - 2.0 * u[:, 1:-1] + u[:, :-2]) / (dy * dy)
    # one-sided in y next to walls (Dirichlet u=0 at walls)
    lap[:, 0] += (u[:, 1] - 2.0 * u[:, 0] + 0.0) / (dy * dy)
    lap[:, -1] += (0.0 - 2.0 * u[:, -1] + u[:, -2]) / (dy * dy)
    # x boundaries handled by available faces (Neumann copy) so above interior stencil is sufficient
    return lap


def lap_v(v, dx, dy):
    lap = np.zeros_like(v)
    lap[1:-1, :] += (v[2:, :] - 2.0 * v[1:-1, :] + v[:-2, :]) / (dx * dx)
    lap[:, 1:-1] += (v[:, 2:] - 2.0 * v[:, 1:-1] + v[:, :-2]) / (dy * dy)
    lap[0, :] += (v[1, :] - 2.0 * v[0, :] + 0.0) / (dx * dx)
    lap[-1, :] += (0.0 - 2.0 * v[-1, :] + v[-2, :]) / (dx * dx)
    return lap

# -------------------- Convective terms (vectorized) --------------------

def compute_convective_u(u, v, dx, dy):
    Nx_plus1, Ny = u.shape
    conv = np.zeros_like(u)

    # du/dx (upwind based on u sign)
    du_dx = np.zeros_like(u)
    pos = u >= 0.0
    # interior faces
    du_dx[1:-1, :] = np.where(pos[1:-1, :], (u[1:-1, :] - u[:-2, :]) / dx, (u[2:, :] - u[1:-1, :]) / dx)
    # boundaries: use one-sided consistent with BCs (forward at left, backward at right)
    du_dx[0, :] = (u[1, :] - u[0, :]) / dx
    du_dx[-1, :] = (u[-1, :] - u[-2, :]) / dx

    # interpolate v to u-locations using vectorized averaging
    v_on_u = np.zeros_like(u)
    # interior u faces i=1..Nx-1
    v_on_u[1:-1, :] = 0.25 * (v[:-1, :-1] + v[:-1, 1:] + v[1:, :-1] + v[1:, 1:])
    # left boundary u-face i=0: average v[0,:] and v[0,1:]
    v_on_u[0, :] = 0.5 * (v[0, :-1] + v[0, 1:])
    # right boundary u-face i=Nx
    v_on_u[-1, :] = 0.5 * (v[-1, :-1] + v[-1, 1:])

    # du/dy using upwind based on v_on_u
    du_dy = np.zeros_like(u)
    pos_v = v_on_u >= 0.0
    du_dy[:, 1:-1] = np.where(pos_v[:, 1:-1], (u[:, 1:-1] - u[:, :-2]) / dy, (u[:, 2:] - u[:, 1:-1]) / dy)
    # one-sided at walls
    du_dy[:, 0] = (u[:, 1] - u[:, 0]) / dy * (-1.0 if pos_v[:, 0].any() else 1.0)
    du_dy[:, -1] = (u[:, -1] - u[:, -2]) / dy

    conv = u * du_dx + v_on_u * du_dy
    return conv


def compute_convective_v(u, v, dx, dy):
    Nx, Ny_plus1 = v.shape
    conv = np.zeros_like(v)

    # dv/dy (upwind based on v sign)
    dv_dy = np.zeros_like(v)
    pos = v >= 0.0
    dv_dy[:, 1:-1] = np.where(pos[:, 1:-1], (v[:, 1:-1] - v[:, :-2]) / dy, (v[:, 2:] - v[:, 1:-1]) / dy)
    dv_dy[:, 0] = (v[:, 1] - v[:, 0]) / dy * (-1.0 if pos[:, 0].any() else 1.0)
    dv_dy[:, -1] = (v[:, -1] - v[:, -2]) / dy

    # interpolate u to v locations vectorized
    u_on_v = np.zeros_like(v)
    # interior v j=1..Ny-1
    u_on_v[:, 1:-1] = 0.25 * (u[:-1, :-1] + u[1:, :-1] + u[:-1, 1:] + u[1:, 1:])
    # bottom v-face j=0
    u_on_v[:, 0] = 0.5 * (u[:-1, 0] + u[1:, 0])
    # top v-face j=Ny
    u_on_v[:, -1] = 0.5 * (u[:-1, -1] + u[1:, -1])

    # dv/dx using upwind based on u_on_v
    dv_dx = np.zeros_like(v)
    pos_u = u_on_v >= 0.0
    dv_dx[1:-1, :] = np.where(pos_u[1:-1, :], (v[1:-1, :] - v[:-2, :]) / dx, (v[2:, :] - v[1:-1, :]) / dx)
    dv_dx[0, :] = (v[1, :] - v[0, :]) / dx * (-1.0 if pos_u[0, :].any() else 1.0)
    dv_dx[-1, :] = (v[-1, :] - v[-2, :]) / dx

    conv = u_on_v * dv_dx + v * dv_dy
    return conv

# -------------------- Divergence computation --------------------

def compute_divergence(u, v, dx, dy):
    # u shape (Nx+1, Ny), v shape (Nx, Ny+1)
    div = (u[1:, :] - u[:-1, :]) / dx + (v[:, 1:] - v[:, :-1]) / dy
    return div

# -------------------- Main solver --------------------

def run_simulation(Lx=2.0, Ly=1.0, Nx=64, Ny=32, t_final=1.0):
    dx, dy, Xc, Yc, u_shape, v_shape, p_shape = create_mac_grid(Lx, Ly, Nx, Ny)

    # Initialize fields
    u = np.zeros(u_shape)  # (Nx+1, Ny)
    v = np.zeros(v_shape)  # (Nx, Ny+1)
    p = np.zeros(p_shape)  # (Nx, Ny)

    # Build D,G,A operators consistent with MAC discretization; include Dirichlet p=0 at outlet
    D, G, A, g_const = build_mac_operators(Nx, Ny, dx, dy, bc_right_dirichlet_value=0.0)
    # Factorize A once
    A_csc = A.tocsc()
    lu = spla.splu(A_csc)

    # Time stepping parameters (explicit conv + explicit visc used here)
    umax_est = 3.0  # rough max inlet amplitude
    cfl = 0.2
    dt_conv = cfl * min(dx, dy) / max(umax_est, 1e-6)
    dt_diff = 0.25 * min(dx * dx, dy * dy) / visc
    dt = min(dt_conv, dt_diff, 0.002)
    Nsteps = int(np.ceil(t_final / dt))
    dt = t_final / Nsteps  # adjust dt to fit exactly

    print(f"Grid: Nx={Nx}, Ny={Ny}, dx={dx:.4e}, dy={dy:.4e}")
    print(f"Time stepping: dt={dt:.6e}, Nsteps={Nsteps}")

    # Progress printing every 10%
    progress_step = max(1, Nsteps // 10)

    t = 0.0
    Nuf = (Nx + 1) * Ny
    Nvf = Nx * (Ny + 1)

    # Precompute face coordinates used for forcing at v faces
    xv = (np.arange(Nx) + 0.5) * dx
    yv = np.arange(Ny + 1) * dy
    XV_v, YV_v = np.meshgrid(xv, yv, indexing='ij')  # (Nx, Ny+1)

    for n in range(1, Nsteps + 1):
        t_new = t + dt
        # Enforce BCs before computing conv/diff
        u, v = apply_bcs(u, v, dx, dy, t)

        # Compute convective terms
        conv_u = compute_convective_u(u, v, dx, dy)
        conv_v = compute_convective_v(u, v, dx, dy)

        # Diffusion terms via laplacian
        lapu = lap_u(u, dx, dy)
        lapv = lap_v(v, dx, dy)

        # External forcing (only fy nonzero). Evaluate at v-locations and at new time.
        f_x = np.zeros_like(u)
        f_y = -np.sin(np.pi * XV_v) * np.sin(np.pi * YV_v) * np.sin(np.pi * t_new)

        # Tentative velocities (explicit time stepping)
        u_star = u + dt * (-conv_u + visc * lapu + f_x)
        v_star = v + dt * (-conv_v + visc * lapv + f_y)

        # Enforce BCs on tentative as well (important for consistent divergence)
        u_star, v_star = apply_bcs(u_star, v_star, dx, dy, t_new)

        # Compute divergence of tentative velocity at cell centers
        div_star = compute_divergence(u_star, v_star, dx, dy)

        # Solve Poisson: A * phi = (1/dt) * div_star - D * g_const
        rhs = (flatten_F(div_star) / dt) - (D @ g_const)
        phi_flat = lu.solve(rhs)

        # Compute face gradients consistently: face_grad = G * phi + g_const
        face_grad = G @ phi_flat + g_const
        # split face_grad into dpdx (u faces) and dpdy (v faces)
        dpdx_flat = face_grad[:Nuf]
        dpdy_flat = face_grad[Nuf:]
        dpdx = dpdx_flat.reshape(u.shape, order='F')
        dpdy = dpdy_flat.reshape(v.shape, order='F')

        # Correct velocities
        u = u_star - dt * dpdx
        v = v_star - dt * dpdy

        # Enforce BCs after projection
        u, v = apply_bcs(u, v, dx, dy, t_new)

        # Update pressure (cell-centered)
        p = phi_flat.reshape(p_shape, order='F')

        t = t_new

        # Print progress only at 10% intervals
        if (n % progress_step) == 0 or n == Nsteps:
            div_after = compute_divergence(u, v, dx, dy)
            div_L2 = np.linalg.norm(flatten_F(div_after)) * np.sqrt(dx * dy)
            print(f"Step {n}/{Nsteps}, t={t:.4f}, ||div||_2={div_L2:.3e}")

    # Interpolate u,v to cell centers for plotting
    u_center = 0.5 * (u[:-1, :] + u[1:, :])
    v_center = 0.5 * (v[:, :-1] + v[:, 1:])

    # Save figure with three contours using RdBu_r colormap
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    levels = 50
    cf0 = axes[0].contourf(Xc, Yc, u_center, levels=levels, cmap='RdBu_r')
    axes[0].set_title('u at t=%.3f' % t)
    axes[0].set_xlabel('x'); axes[0].set_ylabel('y')
    fig.colorbar(cf0, ax=axes[0])

    cf1 = axes[1].contourf(Xc, Yc, v_center, levels=levels, cmap='RdBu_r')
    axes[1].set_title('v at t=%.3f' % t)
    axes[1].set_xlabel('x'); axes[1].set_ylabel('y')
    fig.colorbar(cf1, ax=axes[1])

    cf2 = axes[2].contourf(Xc, Yc, p, levels=levels, cmap='RdBu_r')
    axes[2].set_title('p at t=%.3f' % t)
    axes[2].set_xlabel('x'); axes[2].set_ylabel('y')
    fig.colorbar(cf2, ax=axes[2])

    plt.tight_layout()
    plt.savefig('ns_mac_u_v_p.png', dpi=200)

    # Print concise summary metrics
    print('Final statistics:')
    print(f"u: min={u_center.min():.4e}, max={u_center.max():.4e}")
    print(f"v: min={v_center.min():.4e}, max={v_center.max():.4e}")
    print(f"p: min={p.min():.4e}, max={p.max():.4e}")
    div_final = compute_divergence(u, v, dx, dy)
    div_final_norm = np.linalg.norm(flatten_F(div_final)) * np.sqrt(dx * dy)
    print(f"Final divergence L2 norm: {div_final_norm:.3e}")

    return u, v, p, Xc, Yc

# -------------------- Run if main --------------------

if __name__ == '__main__':
    # Keep grid modest to ensure reasonable run-time
    u, v, p, Xc, Yc = run_simulation(Nx=64, Ny=32, t_final=1.0)
```


#### Script block3:

```python
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# Global viscosity defined as visc = 1.0 / Re
Re = 100.0
visc = 1.0 / Re

# Flatten/reshape helpers using Fortran order consistently
def flatten_F(arr):
    return arr.ravel(order='F')

def reshape_F(flat, shape):
    return np.reshape(flat, shape, order='F')

# -------------------- Grid and problem setup --------------------

def create_mac_grid(Lx, Ly, Nx, Ny):
    dx = Lx / Nx
    dy = Ly / Ny
    # u at vertical faces: shape (Nx+1, Ny)
    u_shape = (Nx + 1, Ny)
    # v at horizontal faces: shape (Nx, Ny+1)
    v_shape = (Nx, Ny + 1)
    # p at cell centers: shape (Nx, Ny)
    p_shape = (Nx, Ny)

    # cell center coordinates (Fortran ordering for arrays)
    xc = (np.arange(Nx) + 0.5) * dx
    yc = (np.arange(Ny) + 0.5) * dy
    Xc, Yc = np.meshgrid(xc, yc, indexing='ij')  # shape (Nx,Ny)
    return dx, dy, Xc, Yc, u_shape, v_shape, p_shape

# -------------------- Operators assembly (MAC-consistent) --------------------

def build_mac_operators(Nx, Ny, dx, dy, bc_right_dirichlet_value=0.0):
    Np = Nx * Ny
    Nuf = (Nx + 1) * Ny
    Nvf = Nx * (Ny + 1)
    Nfaces = Nuf + Nvf

    # index helpers (Fortran order flattening)
    def idx_p(i, j):
        return i + j * Nx

    def idx_u(i, j):
        return i + j * (Nx + 1)

    def idx_v(i, j):
        return Nuf + (i + j * Nx)

    # Build D: maps face values (u and v) -> cell divergence
    rows_D = []
    cols_D = []
    data_D = []
    for j in range(Ny):
        for i in range(Nx):
            row = idx_p(i, j)
            # u contributions: (u_{i+1,j} - u_{i,j})/dx
            up = idx_u(i + 1, j)
            um = idx_u(i, j)
            rows_D.extend([row, row])
            cols_D.extend([up, um])
            data_D.extend([1.0 / dx, -1.0 / dx])
            # v contributions: (v_{i,j+1} - v_{i,j})/dy
            vp = idx_v(i, j + 1)
            vm = idx_v(i, j)
            rows_D.extend([row, row])
            cols_D.extend([vp, vm])
            data_D.extend([1.0 / dy, -1.0 / dy])

    D = sp.csr_matrix((data_D, (rows_D, cols_D)), shape=(Np, Nfaces))

    # Build G: maps cell pressures -> face gradients (dpdx at u faces, dpdy at v faces)
    rows_G = []
    cols_G = []
    data_G = []
    # g_const holds constant contributions from Dirichlet ghost pressures
    g_const = np.zeros(Nfaces)

    # u-faces
    for j in range(Ny):
        for i in range(Nx + 1):
            row = idx_u(i, j)
            # interior u-face between cell (i-1,j) and (i,j)
            if 1 <= i <= Nx - 1:
                # dpdx = (p_{i,j} - p_{i-1,j})/dx
                pR = idx_p(i, j)
                pL = idx_p(i - 1, j)
                rows_G.extend([row, row])
                cols_G.extend([pR, pL])
                data_G.extend([1.0 / dx, -1.0 / dx])
            elif i == 0:
                # left boundary face: Neumann zero-pressure-gradient -> no contribution
                pass
            else:  # i == Nx, rightmost face touching outlet boundary (Dirichlet p at ghost)
                # dpdx = (p_ghost - p_{Nx-1,j})/dx
                # p_ghost = bc_right_dirichlet_value
                pL = idx_p(Nx - 1, j)
                rows_G.append(row)
                cols_G.append(pL)
                data_G.append(-1.0 / dx)
                g_const[row] = bc_right_dirichlet_value / dx

    # v-faces
    for j in range(Ny + 1):
        for i in range(Nx):
            row = idx_v(i, j)
            # interior v-face between cell (i,j-1) and (i,j)
            if 1 <= j <= Ny - 1:
                pT = idx_p(i, j)
                pB = idx_p(i, j - 1)
                rows_G.extend([row, row])
                cols_G.extend([pT, pB])
                data_G.extend([1.0 / dy, -1.0 / dy])
            else:
                # bottom (j==0) or top (j==Ny) boundaries: Neumann (zero gradient) -> no contribution
                pass

    G = sp.csr_matrix((data_G, (rows_G, cols_G)), shape=(Nfaces, Np))

    # Build A = D * G (cell-centered Poisson matrix consistent with D and G)
    A = (D @ G).tocsr()

    return D, G, A, g_const

# -------------------- Boundary conditions and forcing --------------------

def inlet_u_profile(y, t):
    return np.sin(np.pi * y) * (np.sin(np.pi * t) + np.sin(3 * np.pi * t) + np.sin(5 * np.pi * t))


def apply_bcs(u, v, dx, dy, t, Lx=2.0):
    """
    Apply velocity boundary conditions in-place on u (Nx+1,Ny) and v (Nx,Ny+1).
    Centralized BC application so that all stencils/interpolations use the same implied
    ghost-values: walls -> Dirichlet zero (u/v=0), inlet -> prescribed u face & v=0,
    outlet -> Neumann approximated by copying interior face values (du/dx=0).
    """
    Nx_plus1, Ny_local = u.shape
    Nx_local, Ny_plus1 = v.shape

    # Top/bottom walls no-slip (Dirichlet)
    u[:, 0] = 0.0
    u[:, -1] = 0.0
    v[:, 0] = 0.0
    v[:, -1] = 0.0

    # Inlet at x=0: prescribe u at face i=0
    y_u = (np.arange(Ny_local) + 0.5) * dy
    u[0, :] = inlet_u_profile(y_u, t)
    # inlet normal v component zero
    v[0, :] = 0.0

    # Outlet at x=Lx: homogeneous Neumann du/dx=dv/dx=0 approx by copying interior
    u[-1, :] = u[-2, :].copy()
    v[-1, :] = v[-2, :].copy()

    return u, v

# -------------------- Interpolations helpers --------------------

def interp_v_to_u(v):
    # v shape (Nx, Ny+1) -> v_on_u shape (Nx+1, Ny)
    Nx, Nyp1 = v.shape
    Ny = Nyp1 - 1
    v_on_u = np.zeros((Nx + 1, Ny))
    # interior u faces
    v_on_u[1:-1, :] = 0.25 * (v[:-1, :-1] + v[:-1, 1:] + v[1:, :-1] + v[1:, 1:])
    # left u-face i=0: average v[0, :-1] and v[0, 1:]
    v_on_u[0, :] = 0.5 * (v[0, :-1] + v[0, 1:])
    # right u-face i=Nx: average v[-1, :-1] and v[-1, 1:]
    v_on_u[-1, :] = 0.5 * (v[-1, :-1] + v[-1, 1:])
    assert v_on_u.shape == (Nx + 1, Ny)
    return v_on_u


def interp_u_to_v(u):
    # u shape (Nx+1, Ny) -> u_on_v shape (Nx, Ny+1)
    Nxp1, Ny = u.shape
    Nx = Nxp1 - 1
    u_on_v = np.zeros((Nx, Ny + 1))
    # interior v faces j=1..Ny-1
    u_on_v[:, 1:-1] = 0.25 * (u[:-1, :-1] + u[1:, :-1] + u[:-1, 1:] + u[1:, 1:])
    # bottom v-face j=0
    u_on_v[:, 0] = 0.5 * (u[:-1, 0] + u[1:, 0])
    # top v-face j=Ny
    u_on_v[:, -1] = 0.5 * (u[:-1, -1] + u[1:, -1])
    assert u_on_v.shape == (Nx, Ny + 1)
    return u_on_v

# -------------------- Differential operators (vectorized) --------------------

def lap_u(u, dx, dy):
    # u shape (Nx+1, Ny)
    lap = np.zeros_like(u)
    dx2 = dx * dx
    dy2 = dy * dy
    # x second derivative interior
    lap[1:-1, :] += (u[2:, :] - 2.0 * u[1:-1, :] + u[:-2, :]) / dx2
    # y second derivative interior
    lap[:, 1:-1] += (u[:, 2:] - 2.0 * u[:, 1:-1] + u[:, :-2]) / dy2
    # one-sided in y next to walls (use Dirichlet wall value 0)
    # bottom: ghost below = 0
    lap[:, 0] += (u[:, 1] - 2.0 * u[:, 0] + 0.0) / dy2
    # top: ghost above = 0
    lap[:, -1] += (0.0 - 2.0 * u[:, -1] + u[:, -2]) / dy2
    # x boundaries: inlet (i=0) has neighbor i=-1 ghost (for laplacian we use available faces)
    # left: use one-sided with ghost value consistent with apply_bcs (inlet u[0,:] prescribed, ghost left assumed 0 for du/dx?)
    # here we approximate x-boundary lap with available neighbor values; the outlet uses copy so interior stencil is OK
    lap[0, :] += (u[1, :] - 2.0 * u[0, :] + 0.0) / dx2
    lap[-1, :] += (0.0 - 2.0 * u[-1, :] + u[-2, :]) / dx2
    return lap


def lap_v(v, dx, dy):
    # v shape (Nx, Ny+1)
    lap = np.zeros_like(v)
    dx2 = dx * dx
    dy2 = dy * dy
    # x second derivative interior
    lap[1:-1, :] += (v[2:, :] - 2.0 * v[1:-1, :] + v[:-2, :]) / dx2
    # y second derivative interior
    lap[:, 1:-1] += (v[:, 2:] - 2.0 * v[:, 1:-1] + v[:, :-2]) / dy2
    # For v, left boundary (inlet) has v[0,:]=0 (Dirichlet), right boundary copied for Neumann
    lap[0, :] += (v[1, :] - 2.0 * v[0, :] + 0.0) / dx2
    lap[-1, :] += (0.0 - 2.0 * v[-1, :] + v[-2, :]) / dx2
    # y boundaries (walls) have Dirichlet v=0
    lap[:, 0] += (v[:, 1] - 2.0 * v[:, 0] + 0.0) / dy2
    lap[:, -1] += (0.0 - 2.0 * v[:, -1] + v[:, -2]) / dy2
    return lap

# -------------------- Convective terms (vectorized, elementwise upwind) --------------------

def compute_convective_u(u, v, dx, dy):
    # u shape (Nx+1, Ny), v shape (Nx, Ny+1)
    conv = np.zeros_like(u)

    # du/dx using upwind based on u sign (elementwise)
    du_dx = np.zeros_like(u)
    pos = u >= 0.0
    # interior faces
    du_dx[1:-1, :] = np.where(pos[1:-1, :], (u[1:-1, :] - u[:-2, :]) / dx, (u[2:, :] - u[1:-1, :]) / dx)
    # left boundary i=0: backward uses ghost at left (inlet) assume ghost 0 for upwind? use one-sided forward/backward consistent with apply_bcs
    du_dx[0, :] = np.where(pos[0, :], (u[0, :] - 0.0) / dx, (u[1, :] - u[0, :]) / dx)
    # right boundary i=-1: rightmost face uses interior copy for Neumann -> ghost equal to u[-1]
    du_dx[-1, :] = np.where(pos[-1, :], (u[-1, :] - u[-2, :]) / dx, (0.0 - u[-1, :]) / dx)

    # interpolate v to u-locations
    v_on_u = interp_v_to_u(v)

    # du/dy using upwind based on v_on_u (elementwise) with boundary one-sided using wall Dirichlet=0
    du_dy = np.zeros_like(u)
    pos_v = v_on_u >= 0.0
    # interior in y
    du_dy[:, 1:-1] = np.where(pos_v[:, 1:-1], (u[:, 1:-1] - u[:, :-2]) / dy, (u[:, 2:] - u[:, 1:-1]) / dy)
    # bottom j=0: ghost below is wall (u_ghost=0)
    du_dy[:, 0] = np.where(pos_v[:, 0], (u[:, 0] - 0.0) / dy, (u[:, 1] - u[:, 0]) / dy)
    # top j=-1: ghost above is wall (0)
    du_dy[:, -1] = np.where(pos_v[:, -1], (u[:, -1] - u[:, -2]) / dy, (0.0 - u[:, -1]) / dy)

    conv = u * du_dx + v_on_u * du_dy
    return conv


def compute_convective_v(u, v, dx, dy):
    # v shape (Nx, Ny+1), u shape (Nx+1, Ny)
    conv = np.zeros_like(v)

    # dv/dy using upwind based on v sign (elementwise)
    dv_dy = np.zeros_like(v)
    pos = v >= 0.0
    dv_dy[:, 1:-1] = np.where(pos[:, 1:-1], (v[:, 1:-1] - v[:, :-2]) / dy, (v[:, 2:] - v[:, 1:-1]) / dy)
    # bottom j=0: ghost below is wall (v[:,0]=0)
    dv_dy[:, 0] = np.where(pos[:, 0], (v[:, 0] - 0.0) / dy, (v[:, 1] - v[:, 0]) / dy)
    # top j=-1: ghost above is wall (0)
    dv_dy[:, -1] = np.where(pos[:, -1], (v[:, -1] - v[:, -2]) / dy, (0.0 - v[:, -1]) / dy)

    # interpolate u to v locations
    u_on_v = interp_u_to_v(u)

    # dv/dx using upwind based on u_on_v (elementwise)
    dv_dx = np.zeros_like(v)
    pos_u = u_on_v >= 0.0
    # interior in x
    dv_dx[1:-1, :] = np.where(pos_u[1:-1, :], (v[1:-1, :] - v[:-2, :]) / dx, (v[2:, :] - v[1:-1, :]) / dx)
    # left i=0: ghost left for v is inlet v=0
    dv_dx[0, :] = np.where(pos_u[0, :], (v[0, :] - 0.0) / dx, (v[1, :] - v[0, :]) / dx)
    # right i=-1: ghost right equals copied interior (Neumann) so forward difference gives zero
    dv_dx[-1, :] = np.where(pos_u[-1, :], (v[-1, :] - v[-2, :]) / dx, 0.0)

    conv = u_on_v * dv_dx + v * dv_dy
    return conv

# -------------------- Divergence computation --------------------

def compute_divergence(u, v, dx, dy):
    # u shape (Nx+1, Ny), v shape (Nx, Ny+1)
    div = (u[1:, :] - u[:-1, :]) / dx + (v[:, 1:] - v[:, :-1]) / dy
    return div

# -------------------- Main solver --------------------

def run_simulation(Lx=2.0, Ly=1.0, Nx=64, Ny=32, t_final=1.0):
    dx, dy, Xc, Yc, u_shape, v_shape, p_shape = create_mac_grid(Lx, Ly, Nx, Ny)

    # Initialize fields
    u = np.zeros(u_shape)  # (Nx+1, Ny)
    v = np.zeros(v_shape)  # (Nx, Ny+1)
    p = np.zeros(p_shape)  # (Nx, Ny)

    # Build D,G,A operators consistent with MAC discretization; include Dirichlet p=0 at outlet
    D, G, A, g_const = build_mac_operators(Nx, Ny, dx, dy, bc_right_dirichlet_value=0.0)

    # Quick algebraic consistency tests (unit tests) to detect assembly/ordering mismatches
    # Random test for identity: D*(G*p + g_const) == A*p + D*g_const
    p_rand = np.random.randn(Nx * Ny)
    left = D @ (G @ p_rand + g_const)
    right = (A @ p_rand) + (D @ g_const)
    resid_AG = np.linalg.norm(left - right)
    print(f"Operator consistency residual ||D*(G*p+g_const) - (A*p + D*g_const)|| = {resid_AG:.3e}")

    # Factorize A once
    A_csc = A.tocsc()
    lu = spla.splu(A_csc)

    # Divergence-reduction test with random tentative velocities
    # Build random u_star/v_star consistent shapes
    u_test = np.random.randn(*u_shape)
    v_test = np.random.randn(*v_shape)
    # enforce BCs on these test velocities to be consistent
    u_test, v_test = apply_bcs(u_test, v_test, dx, dy, 0.0)
    div_star_test = compute_divergence(u_test, v_test, dx, dy)
    rhs_test = flatten_F(div_star_test) - (D @ g_const)  # dt=1 here just for unit test
    phi_flat_test = lu.solve(rhs_test)
    face_grad_test = G @ phi_flat_test + g_const
    Nuf = (Nx + 1) * Ny
    dpdx_flat_test = face_grad_test[:Nuf]
    dpdy_flat_test = face_grad_test[Nuf:]
    dpdx_test = reshape_F(dpdx_flat_test, u_shape)
    dpdy_test = reshape_F(dpdy_flat_test, v_shape)
    u_corr = u_test - 1.0 * dpdx_test
    v_corr = v_test - 1.0 * dpdy_test
    div_corr = compute_divergence(u_corr, v_corr, dx, dy)
    div_corr_norm = np.linalg.norm(flatten_F(div_corr))
    print(f"Divergence-reduction test ||div(after correction)|| = {div_corr_norm:.3e} (should be ~1e-12)")

    # Time stepping parameters (explicit conv + explicit visc used here)
    umax_est = 3.0  # rough max inlet amplitude
    cfl = 0.2
    dt_conv = cfl * min(dx, dy) / max(umax_est, 1e-6)
    dt_diff = 0.25 * min(dx * dx, dy * dy) / visc
    dt = min(dt_conv, dt_diff, 0.002)
    Nsteps = int(np.ceil(t_final / dt))
    dt = t_final / Nsteps  # adjust dt to fit exactly

    print(f"Grid: Nx={Nx}, Ny={Ny}, dx={dx:.4e}, dy={dy:.4e}")
    print(f"Time stepping: dt={dt:.6e}, Nsteps={Nsteps}")

    # Progress printing every 10%
    progress_step = max(1, Nsteps // 10)

    t = 0.0
    Nuf = (Nx + 1) * Ny
    Nvf = Nx * (Ny + 1)

    # Precompute face coordinates used for forcing at v faces
    xv = (np.arange(Nx) + 0.5) * dx
    yv = np.arange(Ny + 1) * dy
    XV_v, YV_v = np.meshgrid(xv, yv, indexing='ij')  # (Nx, Ny+1)

    for n in range(1, Nsteps + 1):
        t_new = t + dt
        # Enforce BCs before computing conv/diff
        u, v = apply_bcs(u, v, dx, dy, t)

        # Compute convective terms
        conv_u = compute_convective_u(u, v, dx, dy)
        conv_v = compute_convective_v(u, v, dx, dy)

        # Diffusion terms via laplacian
        lapu = lap_u(u, dx, dy)
        lapv = lap_v(v, dx, dy)

        # External forcing (only fy nonzero). Evaluate at v-locations and at new time.
        f_x = np.zeros_like(u)
        f_y = -np.sin(np.pi * XV_v) * np.sin(np.pi * YV_v) * np.sin(np.pi * t_new)

        # Tentative velocities (explicit time stepping)
        u_star = u + dt * (-conv_u + visc * lapu + f_x)
        v_star = v + dt * (-conv_v + visc * lapv + f_y)

        # Enforce BCs on tentative as well (important for consistent divergence)
        u_star, v_star = apply_bcs(u_star, v_star, dx, dy, t_new)

        # Compute divergence of tentative velocity at cell centers
        div_star = compute_divergence(u_star, v_star, dx, dy)

        # Solve Poisson: A * phi = (1/dt) * div_star - D * g_const
        rhs = (flatten_F(div_star) / dt) - (D @ g_const)
        phi_flat = lu.solve(rhs)

        # Optional: small residual check for Poisson solve (diagnostic)
        # Compute residual r = A*phi - rhs
        res_poisson = A @ phi_flat - rhs
        res_norm = np.linalg.norm(res_poisson)
        if n == 1:
            print(f"Poisson residual norm after solve (first step): {res_norm:.3e}")

        # Compute face gradients consistently: face_grad = G * phi + g_const
        face_grad = G @ phi_flat + g_const
        dpdx_flat = face_grad[:Nuf]
        dpdy_flat = face_grad[Nuf:]
        dpdx = reshape_F(dpdx_flat, u.shape)
        dpdy = reshape_F(dpdy_flat, v.shape)

        # Correct velocities
        u = u_star - dt * dpdx
        v = v_star - dt * dpdy

        # Enforce BCs after projection
        u, v = apply_bcs(u, v, dx, dy, t_new)

        # Update pressure (cell-centered)
        p = reshape_F(phi_flat, p_shape)

        t = t_new

        # Print progress only at 10% intervals
        if (n % progress_step) == 0 or n == Nsteps:
            div_after = compute_divergence(u, v, dx, dy)
            div_flat = flatten_F(div_after)
            div_L2 = np.linalg.norm(div_flat)
            div_Linf = np.max(np.abs(div_flat))
            print(f"Step {n}/{Nsteps}, t={t:.4f}, ||div||_2={div_L2:.3e}, ||div||_inf={div_Linf:.3e}")

    # Interpolate u,v to cell centers for plotting
    u_center = 0.5 * (u[:-1, :] + u[1:, :])
    v_center = 0.5 * (v[:, :-1] + v[:, 1:])

    # Save figure with three contours using RdBu_r colormap
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    levels = 50
    cf0 = axes[0].contourf(Xc, Yc, u_center, levels=levels, cmap='RdBu_r')
    axes[0].set_title('u at t=%.3f' % t)
    axes[0].set_xlabel('x'); axes[0].set_ylabel('y')
    fig.colorbar(cf0, ax=axes[0])

    cf1 = axes[1].contourf(Xc, Yc, v_center, levels=levels, cmap='RdBu_r')
    axes[1].set_title('v at t=%.3f' % t)
    axes[1].set_xlabel('x'); axes[1].set_ylabel('y')
    fig.colorbar(cf1, ax=axes[1])

    cf2 = axes[2].contourf(Xc, Yc, p, levels=levels, cmap='RdBu_r')
    axes[2].set_title('p at t=%.3f' % t)
    axes[2].set_xlabel('x'); axes[2].set_ylabel('y')
    fig.colorbar(cf2, ax=axes[2])

    plt.tight_layout()
    plt.savefig('ns_mac_u_v_p.png', dpi=200)

    # Print concise summary metrics
    print('Final statistics:')
    print(f"u: min={u_center.min():.4e}, max={u_center.max():.4e}")
    print(f"v: min={v_center.min():.4e}, max={v_center.max():.4e}")
    print(f"p: min={p.min():.4e}, max={p.max():.4e}")
    div_final = compute_divergence(u, v, dx, dy)
    div_flat_final = flatten_F(div_final)
    div_final_L2 = np.linalg.norm(div_flat_final)
    div_final_Linf = np.max(np.abs(div_flat_final))
    print(f"Final divergence L2 norm: {div_final_L2:.3e}, Linf: {div_final_Linf:.3e}")

    return u, v, p, Xc, Yc

# -------------------- Run if main --------------------

if __name__ == '__main__':
    # Keep grid modest to ensure reasonable run-time
    u, v, p, Xc, Yc = run_simulation(Nx=64, Ny=32, t_final=1.0)
```




### runtime_outputs

#### Output block1

Grid: Nx=64, Ny=32, dx=3.1250e-02, dy=3.1250e-02
Time stepping: dt=2.000000e-03, Nsteps=500
Step 50/500, t=0.1000, ||div||_2=1.709e-01
Step 100/500, t=0.2000, ||div||_2=1.269e-01
Step 150/500, t=0.3000, ||div||_2=1.506e-02
Step 200/500, t=0.4000, ||div||_2=3.086e-02
Step 250/500, t=0.5000, ||div||_2=8.221e-02
Step 300/500, t=0.6000, ||div||_2=3.729e-02
Step 350/500, t=0.7000, ||div||_2=2.118e-02
Step 400/500, t=0.8000, ||div||_2=1.266e-01
Step 450/500, t=0.9000, ||div||_2=1.780e-01
Step 500/500, t=1.0000, ||div||_2=5.422e-02
Final statistics:
u: min=-6.4535e-01, max=4.2238e-01
v: min=-4.5177e-01, max=4.5455e-01
p: min=-4.0691e+01, max=-4.5430e-01
Final divergence L2 norm: 5.422e-02



#### Output block2

Grid: Nx=64, Ny=32, dx=3.1250e-02, dy=3.1250e-02
Time stepping: dt=2.000000e-03, Nsteps=500
Step 50/500, t=0.1000, ||div||_2=1.709e-01
Step 100/500, t=0.2000, ||div||_2=1.269e-01
Step 150/500, t=0.3000, ||div||_2=1.506e-02
Step 200/500, t=0.4000, ||div||_2=3.086e-02
Step 250/500, t=0.5000, ||div||_2=8.221e-02
Step 300/500, t=0.6000, ||div||_2=3.729e-02
Step 350/500, t=0.7000, ||div||_2=2.118e-02
Step 400/500, t=0.8000, ||div||_2=1.266e-01
Step 450/500, t=0.9000, ||div||_2=1.780e-01
Step 500/500, t=1.0000, ||div||_2=5.422e-02
Final statistics:
u: min=-6.4535e-01, max=4.2238e-01
v: min=-4.5177e-01, max=4.5455e-01
p: min=-4.0691e+01, max=-4.5430e-01
Final divergence L2 norm: 5.422e-02



#### Output block3

Operator consistency residual ||D*(G*p+g_const) - (A*p + D*g_const)|| = 2.125e-11
Divergence-reduction test ||div(after correction)|| = 8.107e-13 (should be ~1e-12)
Grid: Nx=64, Ny=32, dx=3.1250e-02, dy=3.1250e-02
Time stepping: dt=2.000000e-03, Nsteps=500
Poisson residual norm after solve (first step): 5.598e-10
Step 50/500, t=0.1000, ||div||_2=5.470e+00, ||div||_inf=3.326e+00
Step 100/500, t=0.2000, ||div||_2=4.061e+00, ||div||_inf=2.416e+00
Step 150/500, t=0.3000, ||div||_2=4.820e-01, ||div||_inf=1.972e-01
Step 200/500, t=0.4000, ||div||_2=9.875e-01, ||div||_inf=5.704e-01
Step 250/500, t=0.5000, ||div||_2=2.631e+00, ||div||_inf=1.570e+00
Step 300/500, t=0.6000, ||div||_2=1.193e+00, ||div||_inf=5.704e-01
Step 350/500, t=0.7000, ||div||_2=6.776e-01, ||div||_inf=1.853e-01
Step 400/500, t=0.8000, ||div||_2=4.050e+00, ||div||_inf=2.416e+00
Step 450/500, t=0.9000, ||div||_2=5.695e+00, ||div||_inf=3.326e+00
Step 500/500, t=1.0000, ||div||_2=1.735e+00, ||div||_inf=4.176e-01
Final statistics:
u: min=-6.4535e-01, max=4.2238e-01
v: min=-4.5177e-01, max=4.5455e-01
p: min=-4.0691e+01, max=-4.5430e-01
Final divergence L2 norm: 1.735e+00, Linf: 4.176e-01




### review_decision
revise

### review_comments
	Current Stage [A/3]
1. MAJOR

Category: accuracy
Issue: Discrete Poisson solve is inconsistent with the staggered (MAC) divergence/gradient operators, so the projection does not enforce a near-zero discrete divergence (final ||div||_2 ~ 5e-2).

Recommendation: Construct the pressure Poisson matrix using the same discrete divergence and gradient operators as used for the velocity update (i.e. assemble D and G on the staggered grid and form A = D*G or explicitly build the cell-centered Laplacian from face-centered flux differences). Include boundary-condition contributions in the RHS consistently (ghost/Dirichlet values or explicit RHS terms). This ensures the projection reduces divergence to machine precision and yields physically consistent pressure/velocity fields.


2. MAJOR

Category: accuracy
Issue: Pressure field magnitude and sign look suspicious (p_min ~ -40) and pressure may be shifted or scaled incorrectly due to sign/BC/RHS inconsistencies in the Poisson assembly.

Recommendation: Verify the sign convention used for the Poisson equation and gradient correction. Ensure you solve the equation A phi = (1/dt) * div_star using A assembled consistently (see previous comment). Check boundary condition treatment: if p=0 at the outlet, ensure contributions from that Dirichlet boundary are applied to the linear system RHS where appropriate (or eliminate Dirichlet dofs). Add a small diagnostic test (manufactured solution or simple steady-state case) where the analytical pressure is known to validate sign and scaling.


3. MAJOR

Category: structure
Issue: Poisson matrix assembly in build_poisson_matrix uses ad-hoc handling for Neumann/Dirichlet neighbors (adding to diagonal for Dirichlet g=0) — this is fragile and not tied to the discrete divergence operator used in compute_divergence/gradient.

Recommendation: Either (a) explicitly assemble A by computing D and G operators on the MAC grid (recommended), or (b) when using a cell-centered Laplacian, derive it from the face flux discretization so boundary conditions seen by A exactly match those used to correct u and v. Avoid relying on special-casing g=0; make the BC treatment explicit and documented in the code.


4. MINOR

Category: structure
Issue: Several loops for interpolation (u_on_v, v_on_u) are implemented in Python with explicit for-loops which will be slow on larger grids; some boundary stencils and boundary derivative approximations are also ad-hoc (e.g. du_dx at i=0 and i=-1).

Recommendation: Vectorize interpolation using array slicing and averaging (e.g. u_on_v = 0.25*(u[:-1,:-1] + u[1:,:-1] + u[:-1,1:] + u[1:,1:]) with careful indexing and order='F' or 'C' consistency). Replace ad-hoc boundary one-sided differences with well-documented one-sided stencils that match the physical BCs. This will improve performance and reduce chance of index errors.


5. MINOR

Category: structure
Issue: apply_bcs signature takes Xc,Yc but does not use them; some BC treatments (copying interior values to set Neumann) are valid but should be documented; compute_convective_* and lap_* functions make assumptions on ghost/Dirichlet values that should be explicit.

Recommendation: Remove unused arguments or use Xc,Yc where needed for clear BC evaluation. Add comments documenting the BC approximation (how Neumann is implemented via copying, which face indices correspond to physical walls/inlet/outlet). Add unit tests for boundary treatments (single step with known input) to catch indexing/BC mistakes early.



	Current Stage [B/3]
1. MAJOR

Category: accuracy
Issue: Incorrect/ambiguous one-sided upwind handling at boundaries in compute_convective_u and compute_convective_v (uses aggregated boolean .any() to choose sign), producing wrong elementwise upwind stencils.

Recommendation: Replace the scalar .any() logic with elementwise upwind decisions using the local normal velocity sign (pos_v, pos_u) for each face. At physical walls use the known ghost values from Dirichlet BCs (u=v=0) to form one-sided differences. Example fixes (conceptual):
- For du_dy at bottom face j=0: du_dy[:,0] = np.where(v_on_u[:,0] >= 0.0, (u[:,0] - 0.0)/dy, (u[:,1] - u[:,0])/dy)
- For du_dy at top face j=-1: du_dy[:,-1] = np.where(v_on_u[:,-1] >= 0.0, (u[:,-1] - u[:,-2])/dy, (0.0 - u[:,-1])/dy)
Perform analogous elementwise fixes for dv_dy/dv_dx boundaries. After changes, verify upwind behaviour with simple manufactured velocity fields.


2. MAJOR

Category: accuracy
Issue: Discrete projection not reducing divergence to solver tolerance; reported ||div||_2 is O(1e-1) rather than near numerical zero, indicating an inconsistency between D, G, A assembly, BC treatment, or ordering/flattening conventions.

Recommendation: Perform targeted consistency/unit tests to locate the mismatch:
- Test algebraic identity: for random cell-centered p, compute face_grad = G @ p + g_const and check that (D @ face_grad) equals (A @ p) + (D @ g_const) (should match to numerical precision).
- Do divergence-reduction test: pick random u_star,v_star, compute div_star, solve for phi as in code and then apply correction; the divergence of corrected velocity must be ~machine epsilon. If not, inspect ordering/flattening (order='F') usage and face ordering indexers (idx_u/idx_v) for consistency.
- Check sign conventions: ensure the Poisson RHS and correction use consistent signs. The derivation in the code expects A = D*G and u_new = u_star - dt * dpdx. Any sign flip in G or in RHS will prevent divergence cancellation.
Fix the inconsistency (most likely in either boundary ghost-value treatment in G/g_const or in the way g_const is subtracted in the RHS). After fixes, divergence L2 should fall to solver tolerance (~1e-10 for direct solve).


3. MAJOR

Category: structure
Issue: Ad-hoc/fragile handling of boundary stencils in lap_u/lap_v and in convective interpolation at boundaries (use of zeros as ghost values and mixed one-sided formulae). This leads to inconsistent treatment of BCs across operators (apply_bcs vs. laplacian/convective stencils).

Recommendation: Unify BC treatment via explicit ghost layers or consistent formulas: either (A) store ghost cells for u and v and fill them from BCs in apply_bcs, then write stencils uniformly using ghost values, or (B) always derive one-sided stencils using the same boundary values (0 at walls, copy at outlet). Avoid mixing 'copy interior' for Neumann and directly inserting 0.0 in laplacian code unless that matches apply_bcs ghost values. Document and centralize BC application (preferably with small helper functions to set ghost-values) so D, G, laplacian and convection all use the same boundary state.


4. MINOR

Category: accuracy
Issue: Convective interpolation and indexing are complex and error-prone (many hand-written averaged stencils for v_on_u and u_on_v). Potential shape/offset mistakes may create subtle errors and aliasing.

Recommendation: Simplify and clearly document interpolation logic. For example compute face-centered quantities using explicit averaged expressions with clear index ranges (and assert shapes). Consider building simple helper functions: interp_v_to_u(v) and interp_u_to_v(u) that (a) explicitly handle interior and boundary indices and (b) include assert checks for shapes/Fortran ordering. Add small unit tests comparing interpolation against a known analytic function to ensure correct placement and parity (e.g., linear field should interpolate exactly).


5. MINOR

Category: structure
Issue: Runtime/diagnostic niceties and numerical output interpretation: large pressure magnitudes and the way divergence norm is reported may obscure debugging. Also some flattening/reshape operations are used in some places but not others.

Recommendation: (a) Report both L2 and Linf divergence norms without multiplication by sqrt(dx*dy) to make magnitude comparisons easier. (b) Add small checks after Poisson solve: print residual ||A*phi - rhs|| or assert it is small when debugging. (c) Ensure flattening/reshaping always uses order='F' (already used in several places) — consider centralizing flatten_F and reshape_F helpers and use them everywhere to avoid subtle order mismatches. (d) Consider plotting or printing a few diagnostic snapshots (max/min of phi, RHS, face_grad) during development to isolate the source of large pressure values.



	Current Stage [C/3]
1. MAJOR

Category: accuracy
Issue: Projection step is inconsistent with velocity Dirichlet boundary enforcement; re-applying BCs after projection reintroduces divergence and leads to poor mass conservation.

Recommendation: Do not re-apply Dirichlet velocity boundary values after the projection correction. Either (A) treat Dirichlet faces as fixed and only update non-Dirichlet faces during the correction (i.e. mask-out faces where velocity is prescribed), or (B) incorporate the velocity Dirichlet conditions into the Poisson/RHS assembly so that the pressure correction respects the prescribed face values. Concretely, create boolean masks u_fixed_mask and v_fixed_mask that mark inlet and wall faces; after computing dpdx/dpdy do:

    u = u_star.copy()
    u[~u_fixed_mask] -= dt * dpdx[~u_fixed_mask]

and analogous for v; do NOT overwrite u/v by apply_bcs() after this step. Alternatively modify A/RHS to enforce the projection consistently for prescribed faces. This fixes the large, time-growing divergence seen in the logs (final ||div||_2 ~ 1.735).


2. MAJOR

Category: accuracy
Issue: Inconsistent treatment of ghost values/one-sided stencils across Laplacian, convective terms, and BC application leads to inconsistency and can degrade accuracy and conservation.

Recommendation: Make the boundary-handling consistent across all operators (convective, diffusive, and projection). Right now different functions assume different ghost values (e.g., lap_u/lap_v use zeros for some ghosts while apply_bcs copies interior values for outlet). This leads to inconsistent stencils and nonphysical behavior. I recommend one of the following consistent strategies:

 - Explicit ghost layers: store and update ghost values that reflect the BCs, and use those ghosts in laplacian and upwind stencils.
 - Or implement stencils that explicitly use the known boundary values (e.g., for Dirichlet walls use the prescribed value, for outlet Neumann use one-sided derivative matched to copying interior values) and document the ghost assumptions.

Specifically, make left/inlet, right/outlet and top/bottom ghost choices consistent in lap_u, lap_v, compute_convective_u, compute_convective_v, and apply_bcs. Add unit tests that check consistency of divergence after one projection when using these BC assumptions.


3. MAJOR

Category: structure
Issue: Poisson solve and face gradient correction do not respect which faces are fixed by velocity Dirichlet BCs; reapplying BCs afterward reintroduces divergence.

Recommendation: The projection routine must not change faces where velocity is prescribed. Currently the code solves the Poisson globally and subtracts dt * (G phi) on all faces and then re-applies BCs — this forces the algorithm to overwrite the projection on Dirichlet faces and thereby breaks mass conservation. Implement face masks for Dirichlet faces and only update free faces (see first suggestion). Also ensure that g_const and D are assembled consistently for the pressure Dirichlet at outlet and that the RHS of Poisson includes contributions from fixed face gradients (if you choose to keep projection changing all faces).


4. MINOR

Category: structure
Issue: BC / ghost-value handling is insufficiently documented and inconsistent between functions which makes maintenance and debugging difficult.

Recommendation: Improve and simplify boundary-related comments and documenting assumptions. For example, clarify the assumed ghost values used in lap_u/lap_v and convective stencils. Replace approximations such as using 0.0 as ghost for left/right in laplacian with explicit code that reads ghost cells or uses the same apply_bcs convention. This reduces subtle bugs and makes it easier to verify correctness.


5. MINOR

Category: runtime
Issue: Poisson factorization with splu is fine for modest grids but will not scale; explicit diffusion makes dt restrictive.

Recommendation: For performance and larger grids consider using an iterative solver (e.g. CG/GMRES) with an appropriate preconditioner (e.g. algebraic multigrid via pyamg or SciPy ILU/AMG) instead of a direct LU for A. Also consider implicit/semi-implicit diffusion (CN or implicit solve for viscous term) to relax dt restriction and reduce Nsteps. Finally, keep the progress printing requirement (every 10%) as implemented, and prefer atol over tol for SciPy solvers as noted in hints.





### iteration_count
3

### rev_image_description
	Current Stage [A/3]


	Current Stage [B/3]


	Current Stage [C/3]




