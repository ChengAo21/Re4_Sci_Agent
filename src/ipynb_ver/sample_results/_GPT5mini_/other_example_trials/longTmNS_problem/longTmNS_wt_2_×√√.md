### 🛌🏻 Re4gent with Multi-Modal Review

#### Configs 🏷️:
MaxCount: 2, MultiModalReview: False, AlternativeNum: 2, MaxComments: 5, ConsModule_ON: False, CUT_OUTEXT: 2000, FROM_Seed: True

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
Use Chorin's projection method on a Staggered (MAC) Grid using the upwind scheme to ensure stability.
Focus on generating runnable code with reasonable physical results, do not over-optimize.
Print concise progress information ONLY every 10% of total steps.


### expanded_prob


### solution_plans


### technical_spec
	Current Stage [A/3]

We implement a semi-implicit Chorin projection solver on a staggered MAC grid. The code builds sparse Laplacian matrices once for u-faces, v-faces, and pressure cells including Dirichlet/Neumann boundary handling. Each time step uses an explicit upwind-ish convective discretization, implicit viscous solves (direct LU factorizations reused), projection via pressure Poisson solve, and velocity correction. Velocity BCs (inlet Dirichlet, top/bottom no-slip, outlet Neumann) and pressure outlet Dirichlet are enforced consistently. Diagnostics (max divergence and kinetic energy) are printed at 10% progress intervals. Final u, v (interpolated to cell centers) and p are contoured in a single figure saved to disk using the 'RdBu_r' colormap. All constants (including visc defined globally) are passed or referenced explicitly to avoid NameError; matrices are assembled only once and reused.


	Current Stage [B/3]
- Build pure Laplacian operators once for u-faces, v-faces, and pressure cells (no Dirichlet rows altered there).
- Form implicit diffusion matrices M = I - dt*visc*L for u and v; then replace rows corresponding to Dirichlet DOFs with identity rows (matrix-level enforcement). Factorize these matrices once.
- For pressure Poisson, replace Dirichlet rows by identity (p=0 at outlet) and factorize once.
- At each time step: enforce instantaneous Dirichlet values into RHS vectors BEFORE solving (u inlet at t+dt, zero for walls), check RHS finite, then solve via re-used LU factors.
- Use upwind-like convective computations with robust boundary handling; check convective terms for finiteness.
- Perform Chorin projection: compute intermediate velocities, solve Poisson for pressure correction (with Dirichlet RHS zeros), correct velocities, enforce BCs, compute diagnostics.
- Save a contour figure of u, v, p at t=1.0 using 'RdBu_r'.
- All matrices assembled ONCE at program start and reused; Dirichlet enforcement on matrices done once as identity rows; RHS is adjusted each step to impose boundary values.

	Current Stage [C/3]
We implement a stable 2D MAC-grid Chorin-projection solver that:
- Assembles discrete divergence operator D once and uses G = -D^T for the discrete gradient, ensuring adjointness (div = -G^T).
- Builds implicit diffusion matrices for u and v once (M = I - dt*visc*L) and enforces Dirichlet rows by setting identity rows and zeroing column coupling; original column entries are stored so that, at each time step, RHS vectors receive the correct contributions from time-varying Dirichlet values.
- Builds the pressure Poisson matrix consistently as A_p = D * G and imposes Dirichlet pressure rows similarly.
- Uses a moderately vectorized upwind convective operator with separate, tested interpolation helpers.
- Performs unit-tests at startup: verifies A_p * phi residuals for manufactured phi and verifies projection reduces divergence to near machine precision.
- Reuses factorizations (spla.splu) for u, v, and p solves; enforces BC-values into RHS using stored column maps before each solve.
- Prints concise progress at 10% intervals, and writes a single contour figure (u, v, p) saved to file (no plt.show()).

Data flow summary:
1) Create grid and MAC indexing functions.
2) Assemble D (velocity-to-divergence) and derive G = -D.T, then A_p = D*G.
3) Assemble diffusion Laplacians Lu, Lv (standard 5-pt with Neumann as ghost=interior) and form M = I - dt*visc*L.
4) Prepare matrices for Dirichlet enforcement (zero columns, set identity rows); store original column entries for RHS adjustments.
5) Factorize matrices once.
6) Time loop: compute forcing, convective terms, build predictor RHS, add Dirichlet contributions to RHS using stored columns and prescribed boundary values at t+dt; solve implicit velocity; compute divergence, build Poisson RHS, add Dirichlet contributions for pressure; solve Poisson; correct velocities via G; re-enforce BCs; diagnostics and progress printing.

Quantitative outputs printed: operator test residual norms and projection divergence check; progress printing every 10% with max_div and KE; final max_div and KE. Plot saved to 'ns_mac_solution.png'.



### python_codes

#### Script block1:

```python
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Global viscosity (use this name exactly as requested)
Re = 100.0
visc = 1.0 / Re

# ---------------------- Utility & Assembler Functions ----------------------
def build_laplace_matrix(nx, ny, dx, dy, dirichlet_mask, neumann_x_mask=None, neumann_y_mask=None):
    """
    Build 2D Laplacian matrix for a grid of size nx * ny with spacing dx, dy.
    dirichlet_mask: boolean array (nx, ny) True where variable is Dirichlet (row set to identity)
    neumann_x_mask: boolean array (nx, ny) True where +x neighbor replaced by mirror (zero normal derivative)
    neumann_y_mask: boolean array (nx, ny) True where +y neighbor replaced by mirror
    Returns sparse CSC matrix.
    """
    N = nx * ny
    A = sp.lil_matrix((N, N))
    idx = lambda i, j: i + j * nx
    invdx2 = 1.0 / dx ** 2
    invdy2 = 1.0 / dy ** 2
    for j in range(ny):
        for i in range(nx):
            k = idx(i, j)
            if dirichlet_mask[i, j]:
                A[k, k] = 1.0
                continue
            diag = 0.0
            # x- neighbors
            # left neighbor exists?
            if i - 1 >= 0:
                A[k, idx(i - 1, j)] = invdx2
                diag -= invdx2
            else:
                # left boundary: assume Neumann (mirror) by default -> contributes like + neighbor
                # treat as ghost mirrored: adds invdx2 to left neighbor which is the right neighbor: handled below if right exists
                # but simpler: no left neighbor => subtract invdx2 from diag (equivalent to zero second derivative)
                diag -= invdx2
            # right neighbor
            if i + 1 < nx:
                if neumann_x_mask is not None and neumann_x_mask[i, j]:
                    # Neumann on +x side: ghost mirrored -> coefficient to left neighbor weighted 2*invdx2
                    # Implemented by setting coefficient for right neighbor as effectively invdx2 but distributing to left
                    # Here we use a standard mirrored approximation: second derivative_x -> 2*(u_{i-1} - u_i)/dx^2
                    # Translate to coefficients: left neighbor gets 2*invdx2, diag gets -2*invdx2
                    # Since left neighbor already added above if exists, adjust accordingly
                    if i - 1 >= 0:
                        A[k, idx(i - 1, j)] = A[k, idx(i - 1, j)] + invdx2  # extra invdx2 to make it 2*invdx2
                        diag -= invdx2  # adjust extra negative diag
                    else:
                        # degenerate: single column, treat as zero operator
                        diag -= invdx2
                else:
                    A[k, idx(i + 1, j)] = invdx2
                    diag -= invdx2
            else:
                # right boundary missing: assume Neumann by default
                diag -= invdx2
            # y- neighbors
            if j - 1 >= 0:
                A[k, idx(i, j - 1)] = invdy2
                diag -= invdy2
            else:
                diag -= invdy2
            if j + 1 < ny:
                if neumann_y_mask is not None and neumann_y_mask[i, j]:
                    # mirror at +y -> add extra invdy2 to lower neighbor if exists
                    if j - 1 >= 0:
                        A[k, idx(i, j - 1)] = A[k, idx(i, j - 1)] + invdy2
                        diag -= invdy2
                    else:
                        diag -= invdy2
                else:
                    A[k, idx(i, j + 1)] = invdy2
                    diag -= invdy2
            else:
                diag -= invdy2
            A[k, k] = diag
    return A.tocsc()


def map_uv_to_1d_indices(shape):
    nx, ny = shape
    n = nx * ny
    return lambda i, j: i + j * nx


# ---------------------- Boundary Conditions & Forcing ----------------------

def inlet_u_profile(y, t, Ly):
    # u(0, y, t) = sin(pi y) * (sin(pi t) + sin(3 pi t) + sin(5 pi t))
    return np.sin(np.pi * y / Ly * Ly) * (np.sin(np.pi * t) + np.sin(3 * np.pi * t) + np.sin(5 * np.pi * t))


# ---------------------- Convective Terms (Upwind-like with simple interpolation) ----------------------

def compute_convective_u(u, v, dx, dy):
    # u shape: (nx_u, ny_u) where nx_u = Nx+1, ny_u = Ny
    nx_u, ny_u = u.shape
    conv = np.zeros_like(u)
    # helper for v interpolation to u-face
    # v shape: (nx_v, ny_v) where nx_v = Nx, ny_v = Ny+1
    nx_v, ny_v = v.shape
    for j in range(ny_u):
        for i in range(nx_u):
            ui = u[i, j]
            # du/dx upwind
            if i == 0:
                du_dx = (u[i, j] - u[i, j]) / dx  # forced 0 (inlet Dirichlet handled separately)
            elif i == nx_u - 1:
                # rightmost: use backward
                du_dx = (u[i, j] - u[i - 1, j]) / dx
            else:
                if ui >= 0:
                    du_dx = (u[i, j] - u[i - 1, j]) / dx
                else:
                    du_dx = (u[i + 1, j] - u[i, j]) / dx
            # interpolate v to u-face (average of nearby v values)
            # approximate v_at_u by averaging up to four surrounding v-values
            # corresponding v indices: v has i_v in [0,nx_v-1], j_v in [0,ny_v-1]
            # u at i corresponds between v indices i-1 and i
            vi_vals = []
            for ii in [i - 1, i]:
                for jj in [j, j + 1]:
                    if 0 <= ii < nx_v and 0 <= jj < ny_v:
                        vi_vals.append(v[ii, jj])
            v_at = np.mean(vi_vals) if vi_vals else 0.0
            # du/dy upwind for v component
            if v_at >= 0:
                # use downward difference
                if j - 1 >= 0:
                    du_dy = (u[i, j] - u[i, j - 1]) / dy
                else:
                    du_dy = (u[i, j + 1] - u[i, j]) / dy if j + 1 < ny_u else 0.0
            else:
                if j + 1 < ny_u:
                    du_dy = (u[i, j + 1] - u[i, j]) / dy
                else:
                    du_dy = (u[i, j] - u[i, j - 1]) / dy if j - 1 >= 0 else 0.0
            conv[i, j] = ui * du_dx + v_at * du_dy
    return conv


def compute_convective_v(u, v, dx, dy):
    nx_v, ny_v = v.shape
    conv = np.zeros_like(v)
    nx_u, ny_u = u.shape
    for j in range(ny_v):
        for i in range(nx_v):
            vi = v[i, j]
            # dv/dy upwind
            if j == 0:
                dv_dy = (v[i, j] - v[i, j]) / dy
            elif j == ny_v - 1:
                dv_dy = (v[i, j] - v[i, j - 1]) / dy
            else:
                if vi >= 0:
                    dv_dy = (v[i, j] - v[i, j - 1]) / dy
                else:
                    dv_dy = (v[i, j + 1] - v[i, j]) / dy
            # interpolate u to v-face
            ui_vals = []
            for ii in [i, i + 1]:
                for jj in [j - 1, j]:
                    if 0 <= ii < nx_u and 0 <= jj < ny_u:
                        ui_vals.append(u[ii, jj])
            u_at = np.mean(ui_vals) if ui_vals else 0.0
            # dv/dx upwind
            if u_at >= 0:
                if i - 1 >= 0:
                    dv_dx = (v[i, j] - v[i - 1, j]) / dx
                else:
                    dv_dx = (v[i + 1, j] - v[i, j]) / dx if i + 1 < nx_v else 0.0
            else:
                if i + 1 < nx_v:
                    dv_dx = (v[i + 1, j] - v[i, j]) / dx
                else:
                    dv_dx = (v[i, j] - v[i - 1, j]) / dx if i - 1 >= 0 else 0.0
            conv[i, j] = u_at * dv_dx + vi * dv_dy
    return conv


# ---------------------- Divergence & Pressure Gradient ----------------------

def divergence(u, v, dx, dy):
    # returns cell-centered divergence of shape (Nx, Ny)
    nx_u, ny_u = u.shape
    nx_v, ny_v = v.shape
    Nx = nx_u - 1
    Ny = ny_v - 1
    div = np.zeros((Nx, Ny))
    for j in range(Ny):
        for i in range(Nx):
            div_x = (u[i + 1, j] - u[i, j]) / dx
            div_y = (v[i, j + 1] - v[i, j]) / dy
            div[i, j] = div_x + div_y
    return div


def grad_p_to_u(p, dx, dy):
    # p shape (Nx, Ny), returns dp/dx at u-faces shape (Nx+1, Ny)
    Nx, Ny = p.shape
    dpdx = np.zeros((Nx + 1, Ny))
    # interior faces
    for j in range(Ny):
        for i in range(1, Nx):
            dpdx[i, j] = (p[i, j] - p[i - 1, j]) / dx
    # left boundary face i=0: use one-sided
    for j in range(Ny):
        dpdx[0, j] = (p[0, j] - 0.0) / dx  # p at left ghost assumed 0? actually p unknown but u at left Dirichlet so will be overwritten
    # right boundary face i=Nx: one-sided
    for j in range(Ny):
        dpdx[Nx, j] = (0.0 - p[Nx - 1, j]) / dx  # p at right ghost is Dirichlet p=0 at outlet
    return dpdx


def grad_p_to_v(p, dx, dy):
    # returns dp/dy at v-faces shape (Nx, Ny+1)
    Nx, Ny = p.shape
    dpdy = np.zeros((Nx, Ny + 1))
    for i in range(Nx):
        for j in range(1, Ny):
            dpdy[i, j] = (p[i, j] - p[i, j - 1]) / dy
    # bottom j=0 one-sided
    for i in range(Nx):
        dpdy[i, 0] = (p[i, 0] - 0.0) / dy
    # top j=Ny one-sided
    for i in range(Nx):
        dpdy[i, Ny] = (0.0 - p[i, Ny - 1]) / dy
    return dpdy


# ---------------------- Main Solver ----------------------

def main():
    # Domain and grid
    Lx = 2.0
    Ly = 1.0
    Nx = 64  # number of pressure cells in x
    Ny = 32  # number of pressure cells in y
    dx = Lx / Nx
    dy = Ly / Ny

    # MAC grid shapes
    nx_u = Nx + 1
    ny_u = Ny
    nx_v = Nx
    ny_v = Ny + 1

    # Time stepping
    t_final = 1.0
    dt = 0.005  # fixed dt for simplicity and consistent progress reporting
    Nt = int(np.ceil(t_final / dt))
    dt = t_final / Nt

    # Initialize fields
    u = np.zeros((nx_u, ny_u))
    v = np.zeros((nx_v, ny_v))
    p = np.zeros((Nx, Ny))

    # Build BC masks
    # u: Dirichlet at left inlet (i==0), top/bottom (j==0 and j==Ny-1)
    dirichlet_u = np.zeros((nx_u, ny_u), dtype=bool)
    dirichlet_u[0, :] = True
    dirichlet_u[:, 0] = True
    dirichlet_u[:, ny_u - 1] = True
    # Neumann for u at right boundary (i==nx_u-1)
    neumann_u_x = np.zeros((nx_u, ny_u), dtype=bool)
    neumann_u_x[nx_u - 1, :] = True

    # v: Dirichlet at left inlet for v (i==0), top/bottom (j==0,j==ny_v-1)
    dirichlet_v = np.zeros((nx_v, ny_v), dtype=bool)
    dirichlet_v[0, :] = True
    dirichlet_v[:, 0] = True
    dirichlet_v[:, ny_v - 1] = True
    # Neumann for v at right side i==nx_v-1
    neumann_v_x = np.zeros((nx_v, ny_v), dtype=bool)
    neumann_v_x[nx_v - 1, :] = True

    # Pressure: Dirichlet at outlet cell centers (i==Nx-1)
    dirichlet_p = np.zeros((Nx, Ny), dtype=bool)
    dirichlet_p[Nx - 1, :] = True
    # Neumann masks for pressure (left/top/bottom) -> handle via mirrored ghost on +x/+y if needed
    neumann_p_x = np.zeros((Nx, Ny), dtype=bool)
    neumann_p_y = np.zeros((Nx, Ny), dtype=bool)
    # set neumann on +x for pressure at i==Nx-1? Actually we have Dirichlet there. For other boundaries default mirrored behavior suffices.

    # Assemble Laplacian matrices once and build (I - dt*visc*L) matrices for u and v
    A_u = build_laplace_matrix(nx_u, ny_u, dx, dy, dirichlet_u, neumann_x_mask=neumann_u_x)
    A_v = build_laplace_matrix(nx_v, ny_v, dx, dy, dirichlet_v, neumann_x_mask=neumann_v_x)
    A_p = build_laplace_matrix(Nx, Ny, dx, dy, dirichlet_p)

    # Build diffusion operators for implicit solve: (I - dt*visc*L)
    I_u = sp.identity(A_u.shape[0], format='csc')
    I_v = sp.identity(A_v.shape[0], format='csc')
    # We'll factorize the matrices (constant coefficients)
    Lu_mat = (I_u - dt * visc * A_u).tocsc()
    Lv_mat = (I_v - dt * visc * A_v).tocsc()
    Lp_mat = A_p.tocsc()  # Poisson: A_p * phi = rhs

    # Factorize using sparse LU for speed/reuse
    print('Factorizing matrices (one-time) ...')
    lu_u = spla.splu(Lu_mat)
    lu_v = spla.splu(Lv_mat)
    lu_p = spla.splu(Lp_mat)

    # Precompute index mapping
    idx_u = map_uv_to_1d_indices((nx_u, ny_u))
    idx_v = map_uv_to_1d_indices((nx_v, ny_v))
    idx_p = map_uv_to_1d_indices((Nx, Ny))

    # Coordinates for inlet profile
    y_u = (np.arange(ny_u) + 0.5) * dy  # u faces y-locations

    # Progress print schedule
    progress_steps = set([int(Nt * frac / 10) for frac in range(1, 11)])
    if 0 in progress_steps:
        progress_steps.remove(0)

    t = 0.0
    print('Starting time loop: Nt =', Nt, ', dt =', dt)
    for n in range(Nt):
        # compute forcing at time t (explicit) at face centers
        # f_x = 0 everywhere (u-equation)
        f_u = np.zeros_like(u)
        # f_y = -sin(pi x) sin(pi y) sin(pi t) applied at v-face centers
        xv = (np.arange(nx_v) + 0.5) * dx
        yv = np.arange(ny_v) * dy  # v faces y positions go from 0..Ny
        # build mesh
        Xv, Yv = np.meshgrid(xv, yv, indexing='xy')
        # note: sin(pi x) sin(pi y) sin(pi t)
        f_v = -np.sin(np.pi * Xv) * np.sin(np.pi * Yv) * np.sin(np.pi * t)
        f_v = f_v.T  # match shape (nx_v, ny_v)

        # Apply boundary conditions to current u,v explicitly before convective computation
        # Inlet u at i=0
        u[0, :] = np.sin(np.pi * y_u) * (np.sin(np.pi * t) + np.sin(3 * np.pi * t) + np.sin(5 * np.pi * t))
        # No-slip top/bottom for u
        u[:, 0] = 0.0
        u[:, ny_u - 1] = 0.0
        # For outlet Neumann u[:, -1] is already free; we can mirror for ghost via setting last face equal to previous
        u[-1, :] = u[-2, :]
        # v BCs
        v[0, :] = 0.0
        v[:, 0] = 0.0
        v[:, ny_v - 1] = 0.0
        v[-1, :] = v[-2, :]

        # Convective terms
        conv_u = compute_convective_u(u, v, dx, dy)
        conv_v = compute_convective_v(u, v, dx, dy)

        # Predictor step: RHS for implicit diffusion solve (I - dt*visc*L) u* = u^n - dt * conv + dt * f
        # Flatten fields to 1D matching matrix ordering
        rhs_u = (u - dt * conv_u + dt * f_u).ravel(order='F')
        rhs_v = (v - dt * conv_v + dt * f_v).ravel(order='F')

        # Solve for intermediate velocities
        u_star_flat = lu_u.solve(rhs_u)
        v_star_flat = lu_v.solve(rhs_v)
        u_star = u_star_flat.reshape(u.shape, order='F')
        v_star = v_star_flat.reshape(v.shape, order='F')

        # Re-apply Dirichlet BCs exactly to u_star and v_star
        u_star[0, :] = np.sin(np.pi * y_u) * (np.sin(np.pi * t) + np.sin(3 * np.pi * t) + np.sin(5 * np.pi * t))
        u_star[:, 0] = 0.0
        u_star[:, ny_u - 1] = 0.0
        u_star[-1, :] = u_star[-2, :]

        v_star[0, :] = 0.0
        v_star[:, 0] = 0.0
        v_star[:, ny_v - 1] = 0.0
        v_star[-1, :] = v_star[-2, :]

        # Compute divergence of u_star
        div_star = divergence(u_star, v_star, dx, dy)
        rhs_p_vec = (div_star / dt).ravel(order='F')

        # Solve Poisson for pressure correction (phi). Dirichlet p at outlet already enforced in A_p assembly.
        phi_flat = lu_p.solve(rhs_p_vec)
        phi = phi_flat.reshape((Nx, Ny), order='F')

        # Update pressure (standard Chorin: p^{n+1} = phi)
        p = phi.copy()

        # Correct velocities: u_new = u_star - dt * dp/dx, v_new = v_star - dt * dp/dy
        dpdx = grad_p_to_u(p, dx, dy)
        dpdy = grad_p_to_v(p, dx, dy)
        u = u_star - dt * dpdx
        v = v_star - dt * dpdy

        # Enforce velocity BCs again
        u[0, :] = np.sin(np.pi * y_u) * (np.sin(np.pi * (t + dt)) + np.sin(3 * np.pi * (t + dt)) + np.sin(5 * np.pi * (t + dt)))
        u[:, 0] = 0.0
        u[:, ny_u - 1] = 0.0
        u[-1, :] = u[-2, :]

        v[0, :] = 0.0
        v[:, 0] = 0.0
        v[:, ny_v - 1] = 0.0
        v[-1, :] = v[-2, :]

        # Diagnostics
        div_after = divergence(u, v, dx, dy)
        max_div = np.max(np.abs(div_after))
        # Kinetic energy approximate integral
        u_center = 0.5 * (u[0:Nx, :] + u[1:Nx + 1, :])
        v_center = 0.5 * (v[:, 0:Ny] + v[:, 1:Ny + 1])
        KE = 0.5 * np.sum(u_center ** 2 + v_center ** 2) * dx * dy

        t += dt
        # Print progress every 10% of total steps
        if (n + 1) in progress_steps:
            pct = int(100 * (n + 1) / Nt)
            print(f'Progress: {pct}% | step {n+1}/{Nt} | t={t:.3f} | max_div={max_div:.3e} | KE={KE:.6e}')

    # Final diagnostics
    print(f'Finished. t={t:.6f} | max_div_final={max_div:.3e} | KE_final={KE:.6e}')

    # Interpolate u and v to cell centers for plotting
    u_center = 0.5 * (u[0:Nx, :] + u[1:Nx + 1, :])
    v_center = 0.5 * (v[:, 0:Ny] + v[:, 1:Ny + 1])

    # Coordinates for cell centers
    xc = (np.arange(Nx) + 0.5) * dx
    yc = (np.arange(Ny) + 0.5) * dy
    Xc, Yc = np.meshgrid(xc, yc, indexing='xy')

    # Create one figure with 3 contours: u, v, p
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    cs0 = axes[0].contourf(Xc, Yc, u_center.T, levels=50, cmap='RdBu_r')
    axes[0].set_title('u at t=1.0')
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    fig.colorbar(cs0, ax=axes[0])

    cs1 = axes[1].contourf(Xc, Yc, v_center.T, levels=50, cmap='RdBu_r')
    axes[1].set_title('v at t=1.0')
    axes[1].set_xlabel('x')
    fig.colorbar(cs1, ax=axes[1])

    cs2 = axes[2].contourf(Xc, Yc, p.T, levels=50, cmap='RdBu_r')
    axes[2].set_title('p at t=1.0')
    axes[2].set_xlabel('x')
    fig.colorbar(cs2, ax=axes[2])

    plt.tight_layout()
    plt.savefig('ns_mac_solution.png', dpi=200)
    print('Saved figure to ns_mac_solution.png')


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

# Global viscosity (use this name exactly as requested)
Re = 100.0
visc = 1.0 / Re

# ---------------------- Utility & Assembler Functions ----------------------
def build_laplace_matrix(nx, ny, dx, dy, neumann_x_mask=None, neumann_y_mask=None):
    """
    Build 2D Laplacian matrix (standard 5-point) on an nx-by-ny grid with spacing dx,dy.
    Returns sparse CSC matrix with ordering k = i + j*nx (Fortran-style / column-major).
    Handles mirror-type Neumann on +x / +y sides via masks (optional).
    Does NOT apply Dirichlet identity rows (this must be done on the final operator).
    """
    N = nx * ny
    A = sp.lil_matrix((N, N))
    idx = lambda i, j: i + j * nx
    invdx2 = 1.0 / dx ** 2
    invdy2 = 1.0 / dy ** 2

    for j in range(ny):
        for i in range(nx):
            k = idx(i, j)
            diag = 0.0
            # left neighbor
            if i - 1 >= 0:
                A[k, idx(i - 1, j)] = invdx2
                diag -= invdx2
            else:
                # left physical boundary: assume mirror if neumann on -x would be used,
                # but default behavior: one-sided (second derivative uses only right neighbor)
                diag -= invdx2
            # right neighbor
            if i + 1 < nx:
                # handle +x Neumann (mirror) at this cell if requested
                if neumann_x_mask is not None and neumann_x_mask[i, j]:
                    # mirror: approximate ghost = interior neighbor -> gives extra weight to left neighbor
                    # effectively add invdx2 to left neighbor if exists, else leave diag
                    if i - 1 >= 0:
                        A[k, idx(i - 1, j)] = A[k, idx(i - 1, j)] + invdx2
                        diag -= invdx2
                    else:
                        diag -= invdx2
                else:
                    A[k, idx(i + 1, j)] = invdx2
                    diag -= invdx2
            else:
                # missing right neighbor: one-sided
                diag -= invdx2
            # bottom neighbor
            if j - 1 >= 0:
                A[k, idx(i, j - 1)] = invdy2
                diag -= invdy2
            else:
                diag -= invdy2
            # top neighbor
            if j + 1 < ny:
                if neumann_y_mask is not None and neumann_y_mask[i, j]:
                    if j - 1 >= 0:
                        A[k, idx(i, j - 1)] = A[k, idx(i, j - 1)] + invdy2
                        diag -= invdy2
                    else:
                        diag -= invdy2
                else:
                    A[k, idx(i, j + 1)] = invdy2
                    diag -= invdy2
            else:
                diag -= invdy2
            A[k, k] = diag
    return A.tocsc()


def map_uv_to_1d_indices(shape):
    nx, ny = shape
    return lambda i, j: i + j * nx

# ---------------------- Boundary Conditions & Forcing ----------------------

def inlet_u_profile(y, t):
    # u(0, y, t) = sin(pi y) * (sin(pi t) + sin(3 pi t) + sin(5 pi t))
    return np.sin(np.pi * y) * (np.sin(np.pi * t) + np.sin(3 * np.pi * t) + np.sin(5 * np.pi * t))

# ---------------------- Convective Terms (safer, with checks) ----------------------

def compute_convective_u(u, v, dx, dy):
    # Conservative-ish upwind for u on u-faces (slow but robust).
    nx_u, ny_u = u.shape
    conv = np.zeros_like(u)
    nx_v, ny_v = v.shape
    for j in range(ny_u):
        for i in range(nx_u):
            ui = u[i, j]
            # du/dx upwind: prefer one-sided stencils that exist
            if i == 0:
                # use forward difference
                du_dx = (u[i + 1, j] - u[i, j]) / dx if nx_u > 1 else 0.0
            elif i == nx_u - 1:
                du_dx = (u[i, j] - u[i - 1, j]) / dx
            else:
                if ui >= 0:
                    du_dx = (u[i, j] - u[i - 1, j]) / dx
                else:
                    du_dx = (u[i + 1, j] - u[i, j]) / dx
            # interpolate v to u-face robustly
            vi_vals = []
            for ii in (i - 1, i):
                for jj in (j, j + 1):
                    if 0 <= ii < nx_v and 0 <= jj < ny_v:
                        vi_vals.append(v[ii, jj])
            v_at = float(np.mean(vi_vals)) if vi_vals else 0.0
            # du/dy upwind
            if v_at >= 0:
                if j - 1 >= 0:
                    du_dy = (u[i, j] - u[i, j - 1]) / dy
                elif j + 1 < ny_u:
                    du_dy = (u[i, j + 1] - u[i, j]) / dy
                else:
                    du_dy = 0.0
            else:
                if j + 1 < ny_u:
                    du_dy = (u[i, j + 1] - u[i, j]) / dy
                elif j - 1 >= 0:
                    du_dy = (u[i, j] - u[i, j - 1]) / dy
                else:
                    du_dy = 0.0
            conv[i, j] = ui * du_dx + v_at * du_dy
    # sanity
    if not np.all(np.isfinite(conv)):
        raise RuntimeError('Non-finite convective term for u encountered')
    return conv


def compute_convective_v(u, v, dx, dy):
    nx_v, ny_v = v.shape
    conv = np.zeros_like(v)
    nx_u, ny_u = u.shape
    for j in range(ny_v):
        for i in range(nx_v):
            vi = v[i, j]
            # dv/dy upwind
            if j == 0:
                dv_dy = (v[i, j + 1] - v[i, j]) / dy if ny_v > 1 else 0.0
            elif j == ny_v - 1:
                dv_dy = (v[i, j] - v[i, j - 1]) / dy
            else:
                if vi >= 0:
                    dv_dy = (v[i, j] - v[i, j - 1]) / dy
                else:
                    dv_dy = (v[i, j + 1] - v[i, j]) / dy
            # interpolate u to v-face robustly
            ui_vals = []
            for ii in (i, i + 1):
                for jj in (j - 1, j):
                    if 0 <= ii < nx_u and 0 <= jj < ny_u:
                        ui_vals.append(u[ii, jj])
            u_at = float(np.mean(ui_vals)) if ui_vals else 0.0
            # dv/dx upwind
            if u_at >= 0:
                if i - 1 >= 0:
                    dv_dx = (v[i, j] - v[i - 1, j]) / dx
                elif i + 1 < nx_v:
                    dv_dx = (v[i + 1, j] - v[i, j]) / dx
                else:
                    dv_dx = 0.0
            else:
                if i + 1 < nx_v:
                    dv_dx = (v[i + 1, j] - v[i, j]) / dx
                elif i - 1 >= 0:
                    dv_dx = (v[i, j] - v[i - 1, j]) / dx
                else:
                    dv_dx = 0.0
            conv[i, j] = u_at * dv_dx + vi * dv_dy
    if not np.all(np.isfinite(conv)):
        raise RuntimeError('Non-finite convective term for v encountered')
    return conv

# ---------------------- Divergence & Pressure Gradient ----------------------

def divergence(u, v, dx, dy):
    nx_u, ny_u = u.shape
    nx_v, ny_v = v.shape
    Nx = nx_u - 1
    Ny = ny_v - 1
    div = np.zeros((Nx, Ny))
    for j in range(Ny):
        for i in range(Nx):
            div_x = (u[i + 1, j] - u[i, j]) / dx
            div_y = (v[i, j + 1] - v[i, j]) / dy
            div[i, j] = div_x + div_y
    return div


def grad_p_to_u(p, dx, dy):
    Nx, Ny = p.shape
    dpdx = np.zeros((Nx + 1, Ny))
    for j in range(Ny):
        for i in range(1, Nx):
            dpdx[i, j] = (p[i, j] - p[i - 1, j]) / dx
    # left face: use one-sided (assume ghost p= p[0]) -> zero gradient approximation
    for j in range(Ny):
        dpdx[0, j] = (p[0, j] - p[0, j]) / dx
    # right face: p at outlet Dirichlet = 0, use (0 - p[Nx-1])/dx
    for j in range(Ny):
        dpdx[Nx, j] = (0.0 - p[Nx - 1, j]) / dx
    return dpdx


def grad_p_to_v(p, dx, dy):
    Nx, Ny = p.shape
    dpdy = np.zeros((Nx, Ny + 1))
    for i in range(Nx):
        for j in range(1, Ny):
            dpdy[i, j] = (p[i, j] - p[i, j - 1]) / dy
    for i in range(Nx):
        dpdy[i, 0] = (p[i, 0] - p[i, 0]) / dy
        dpdy[i, Ny] = (0.0 - p[i, Ny - 1]) / dy
    return dpdy

# ---------------------- Main Solver ----------------------

def main():
    # Domain and grid
    Lx = 2.0
    Ly = 1.0
    Nx = 64  # number of pressure cells in x
    Ny = 32  # number of pressure cells in y
    dx = Lx / Nx
    dy = Ly / Ny

    # MAC grid shapes
    nx_u = Nx + 1
    ny_u = Ny
    nx_v = Nx
    ny_v = Ny + 1

    # Time stepping
    t_final = 1.0
    dt = 0.005
    Nt = int(np.ceil(t_final / dt))
    dt = t_final / Nt

    # Initialize fields
    u = np.zeros((nx_u, ny_u))
    v = np.zeros((nx_v, ny_v))
    p = np.zeros((Nx, Ny))

    # Build BC masks (used for enforcing Dirichlet rows)
    dirichlet_u = np.zeros((nx_u, ny_u), dtype=bool)
    dirichlet_u[0, :] = True  # inlet
    dirichlet_u[:, 0] = True  # bottom no-slip
    dirichlet_u[:, ny_u - 1] = True  # top no-slip

    neumann_u_x = np.zeros((nx_u, ny_u), dtype=bool)
    neumann_u_x[nx_u - 1, :] = True  # outlet Neumann on rightmost u-faces

    dirichlet_v = np.zeros((nx_v, ny_v), dtype=bool)
    dirichlet_v[0, :] = True
    dirichlet_v[:, 0] = True
    dirichlet_v[:, ny_v - 1] = True

    neumann_v_x = np.zeros((nx_v, ny_v), dtype=bool)
    neumann_v_x[nx_v - 1, :] = True

    dirichlet_p = np.zeros((Nx, Ny), dtype=bool)
    dirichlet_p[Nx - 1, :] = True  # pressure Dirichlet at outlet cells

    # Assemble Laplacian matrices (pure Laplace operator; no Dirichlet rows)
    A_u = build_laplace_matrix(nx_u, ny_u, dx, dy, neumann_x_mask=neumann_u_x)
    A_v = build_laplace_matrix(nx_v, ny_v, dx, dy, neumann_x_mask=neumann_v_x)
    A_p = build_laplace_matrix(Nx, Ny, dx, dy)

    # Build implicit diffusion matrices: M = I - dt*visc*L
    I_u = sp.identity(A_u.shape[0], format='csc')
    I_v = sp.identity(A_v.shape[0], format='csc')
    Lu_mat = (I_u - dt * visc * A_u).tolil()  # modify rows easily
    Lv_mat = (I_v - dt * visc * A_v).tolil()

    # Impose identity rows on Lu_mat/Lv_mat for Dirichlet DOFs (rows -> identity). Do this once.
    # Create index lists
    idx_u = map_uv_to_1d_indices((nx_u, ny_u))
    idx_v = map_uv_to_1d_indices((nx_v, ny_v))
    dirichlet_u_idx = [idx_u(i, j) for i in range(nx_u) for j in range(ny_u) if dirichlet_u[i, j]]
    dirichlet_v_idx = [idx_v(i, j) for i in range(nx_v) for j in range(ny_v) if dirichlet_v[i, j]]

    for k in dirichlet_u_idx:
        Lu_mat.rows[k] = [k]
        Lu_mat.data[k] = [1.0]
    for k in dirichlet_v_idx:
        Lv_mat.rows[k] = [k]
        Lv_mat.data[k] = [1.0]

    Lu_mat = Lu_mat.tocsc()
    Lv_mat = Lv_mat.tocsc()

    # Poisson matrix: we must enforce Dirichlet rows for pressure and factorize once
    Lp_mat = A_p.tolil()
    idx_p = map_uv_to_1d_indices((Nx, Ny))
    dirichlet_p_idx = [idx_p(i, j) for i in range(Nx) for j in range(Ny) if dirichlet_p[i, j]]
    if len(dirichlet_p_idx) == 0:
        raise RuntimeError('At least one Dirichlet pressure cell is required to fix pressure zero.')
    for k in dirichlet_p_idx:
        Lp_mat.rows[k] = [k]
        Lp_mat.data[k] = [1.0]
    Lp_mat = Lp_mat.tocsc()

    # Factorize matrices once (Dirichlet rows already set in matrices)
    print('Factorizing matrices (one-time) ...')
    lu_u = spla.splu(Lu_mat)
    lu_v = spla.splu(Lv_mat)
    lu_p = spla.splu(Lp_mat)

    # Coordinates for inlet profile (u faces y-locations)
    y_u = (np.arange(ny_u) + 0.5) * dy

    # Progress print schedule: print at exact 10% increments
    progress_steps = set([int(Nt * frac / 10) for frac in range(1, 11)])
    if 0 in progress_steps:
        progress_steps.remove(0)

    t = 0.0
    print('Starting time loop: Nt =', Nt, ', dt =', dt)

    for n in range(Nt):
        # Forcing at current time t (explicit in predictor)
        f_u = np.zeros_like(u)
        xv = (np.arange(nx_v) + 0.5) * dx
        yv = np.arange(ny_v) * dy
        Xv, Yv = np.meshgrid(xv, yv, indexing='xy')
        f_v = (-np.sin(np.pi * Xv) * np.sin(np.pi * Yv) * np.sin(np.pi * t)).T

        # Enforce Dirichlet BCs on current fields for convective evaluation (use current t)
        u[0, :] = inlet_u_profile(y_u, t)
        u[:, 0] = 0.0
        u[:, ny_u - 1] = 0.0
        u[-1, :] = u[-2, :]

        v[0, :] = 0.0
        v[:, 0] = 0.0
        v[:, ny_v - 1] = 0.0
        v[-1, :] = v[-2, :]

        # Convective terms (explicit)
        conv_u = compute_convective_u(u, v, dx, dy)
        conv_v = compute_convective_v(u, v, dx, dy)

        # Predictor RHS (flatten Fortran-order to match matrix ordering)
        rhs_u = (u - dt * conv_u + dt * f_u).ravel(order='F')
        rhs_v = (v - dt * conv_v + dt * f_v).ravel(order='F')

        # Before solving, impose Dirichlet values into RHS for Dirichlet DOFs (use boundary at t+dt)
        u_inlet_tnp1 = inlet_u_profile(y_u, t + dt)
        # set u dirichlet RHS entries to value at t+dt
        for j in range(ny_u):
            k = idx_u(0, j)
            rhs_u[k] = u_inlet_tnp1[j]
        # top/bottom u are zero
        for i in range(nx_u):
            rhs_u[idx_u(i, 0)] = 0.0
            rhs_u[idx_u(i, ny_u - 1)] = 0.0
        # v Dirichlet RHS zeros
        for j in range(ny_v):
            rhs_v[idx_v(0, j)] = 0.0
        for i in range(nx_v):
            rhs_v[idx_v(i, 0)] = 0.0
            rhs_v[idx_v(i, ny_v - 1)] = 0.0

        # Sanity check before solves
        if not (np.all(np.isfinite(rhs_u)) and np.all(np.isfinite(rhs_v))):
            raise RuntimeError('Non-finite RHS for velocity solve encountered')

        # Solve implicit diffusion systems (Dirichlet enforced via modified matrices)
        u_star_flat = lu_u.solve(rhs_u)
        v_star_flat = lu_v.solve(rhs_v)
        u_star = u_star_flat.reshape(u.shape, order='F')
        v_star = v_star_flat.reshape(v.shape, order='F')

        # Ensure Dirichlet BCs satisfied (numerical round-off)
        u_star[0, :] = u_inlet_tnp1
        u_star[:, 0] = 0.0
        u_star[:, ny_u - 1] = 0.0
        v_star[0, :] = 0.0
        v_star[:, 0] = 0.0
        v_star[:, ny_v - 1] = 0.0

        # Compute divergence of intermediate velocity
        div_star = divergence(u_star, v_star, dx, dy)
        rhs_p_vec = (div_star / dt).ravel(order='F')
        # Enforce pressure Dirichlet RHS entries to zero
        for k in dirichlet_p_idx:
            rhs_p_vec[k] = 0.0
        if not np.all(np.isfinite(rhs_p_vec)):
            raise RuntimeError('Non-finite RHS for pressure solve encountered')

        # Solve Poisson (Dirichlet rows in Lp_mat enforce p=0 at outlet)
        phi_flat = lu_p.solve(rhs_p_vec)
        phi = phi_flat.reshape((Nx, Ny), order='F')

        # Pressure is phi (Chorin projection)
        p = phi.copy()

        # Correct velocities
        dpdx = grad_p_to_u(p, dx, dy)
        dpdy = grad_p_to_v(p, dx, dy)
        u = u_star - dt * dpdx
        v = v_star - dt * dpdy

        # Re-enforce velocity Dirichlet boundary values at t+dt exactly
        u[0, :] = u_inlet_tnp1
        u[:, 0] = 0.0
        u[:, ny_u - 1] = 0.0
        u[-1, :] = u[-2, :]
        v[0, :] = 0.0
        v[:, 0] = 0.0
        v[:, ny_v - 1] = 0.0
        v[-1, :] = v[-2, :]

        # Diagnostics
        div_after = divergence(u, v, dx, dy)
        max_div = np.max(np.abs(div_after))
        u_center = 0.5 * (u[0:Nx, :] + u[1:Nx + 1, :])
        v_center = 0.5 * (v[:, 0:Ny] + v[:, 1:Ny + 1])
        KE = 0.5 * np.sum(u_center ** 2 + v_center ** 2) * dx * dy

        t += dt
        if (n + 1) in progress_steps:
            pct = int(100 * (n + 1) / Nt)
            print(f'Progress: {pct}% | step {n+1}/{Nt} | t={t:.3f} | max_div={max_div:.3e} | KE={KE:.6e}')

    # Final diagnostics
    print(f'Finished. t={t:.6f} | max_div_final={max_div:.3e} | KE_final={KE:.6e}')

    # Interpolate u and v to cell centers for plotting
    u_center = 0.5 * (u[0:Nx, :] + u[1:Nx + 1, :])
    v_center = 0.5 * (v[:, 0:Ny] + v[:, 1:Ny + 1])

    # Coordinates for cell centers
    xc = (np.arange(Nx) + 0.5) * dx
    yc = (np.arange(Ny) + 0.5) * dy
    Xc, Yc = np.meshgrid(xc, yc, indexing='xy')

    # Save contour figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    cs0 = axes[0].contourf(Xc, Yc, u_center.T, levels=50, cmap='RdBu_r')
    axes[0].set_title('u at t=1.0')
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    fig.colorbar(cs0, ax=axes[0])

    cs1 = axes[1].contourf(Xc, Yc, v_center.T, levels=50, cmap='RdBu_r')
    axes[1].set_title('v at t=1.0')
    axes[1].set_xlabel('x')
    fig.colorbar(cs1, ax=axes[1])

    cs2 = axes[2].contourf(Xc, Yc, p.T, levels=50, cmap='RdBu_r')
    axes[2].set_title('p at t=1.0')
    axes[2].set_xlabel('x')
    fig.colorbar(cs2, ax=axes[2])

    plt.tight_layout()
    plt.savefig('ns_mac_solution.png', dpi=200)
    print('Saved figure to ns_mac_solution.png')


if __name__ == '__main__':
    main()
```



#### Script block3:

#!/usr/bin/env python3
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Global viscosity (use this name exactly as requested)
Re = 100.0
visc = 1.0 / Re

# ---------------------- Index helpers ----------------------
def idx_map(nx, ny):
    return lambda i, j: i + j * nx

# ---------------------- Discrete operators assembly ----------------------

def build_divergence_matrix(Nx, Ny, dx, dy):
    """
    Build sparse divergence operator D mapping [u_flat; v_flat] -> div_flat
    Ordering: u flattened Fortran on (Nx+1, Ny), then v flattened Fortran on (Nx, Ny+1)
    Returns D (Np x (Nu+Nv)) and index functions idx_u, idx_v, idx_p
    """
    nx_u = Nx + 1
    ny_u = Ny
    nx_v = Nx
    ny_v = Ny + 1
    Nu = nx_u * ny_u
    Nv = nx_v * ny_v
    Np = Nx * Ny

    idx_u = idx_map(nx_u, ny_u)
    idx_v = idx_map(nx_v, ny_v)
    idx_p = idx_map(Nx, Ny)

    rows = []
    cols = []
    data = []

    for j in range(Ny):
        for i in range(Nx):
            p_k = idx_p(i, j)
            # u contribution: +1/dx at face i+1, -1/dx at face i
            k_u_right = idx_u(i + 1, j)
            rows.append(p_k); cols.append(k_u_right); data.append(1.0 / dx)
            k_u_left = idx_u(i, j)
            rows.append(p_k); cols.append(k_u_left); data.append(-1.0 / dx)
            # v contribution: +1/dy at face j+1, -1/dy at face j
            k_v_top = idx_v(i, j + 1)
            rows.append(p_k); cols.append(Nu + k_v_top); data.append(1.0 / dy)
            k_v_bottom = idx_v(i, j)
            rows.append(p_k); cols.append(Nu + k_v_bottom); data.append(-1.0 / dy)

    D = sp.csr_matrix((data, (rows, cols)), shape=(Np, Nu + Nv))
    return D, idx_u, idx_v, idx_p


def build_laplace_matrix_standard(nx, ny, dx, dy):
    """
    Standard 5-point Laplacian on a regular grid of nx-by-ny unknowns.
    Missing neighbors (domain boundary) are handled with a ghost=interior (Neumann zero) assumption
    so that second derivatives reduce naturally at boundaries. This assembly does NOT apply Dirichlet rows.
    Returns CSC sparse matrix.
    """
    N = nx * ny
    idx = idx_map(nx, ny)
    invdx2 = 1.0 / dx ** 2
    invdy2 = 1.0 / dy ** 2
    A = sp.lil_matrix((N, N))
    for j in range(ny):
        for i in range(nx):
            k = idx(i, j)
            diag = 0.0
            # left
            if i - 1 >= 0:
                A[k, idx(i - 1, j)] = invdx2
                diag -= invdx2
            else:
                # ghost=interior -> contributes -invdx2 to diag
                diag -= invdx2
            # right
            if i + 1 < nx:
                A[k, idx(i + 1, j)] = invdx2
                diag -= invdx2
            else:
                diag -= invdx2
            # bottom
            if j - 1 >= 0:
                A[k, idx(i, j - 1)] = invdy2
                diag -= invdy2
            else:
                diag -= invdy2
            # top
            if j + 1 < ny:
                A[k, idx(i, j + 1)] = invdy2
                diag -= invdy2
            else:
                diag -= invdy2
            A[k, k] = diag
    return A.tocsc()

# ---------------------- Dirichlet enforcement utilities ----------------------

def prepare_matrix_with_dirichlet(A_csc, dirichlet_idx):
    """
    Modify A so that rows corresponding to dirichlet_idx become identity rows
    and the columns corresponding to these dofs are zeroed (except diagonal).
    Store original column off-diagonal entries so RHS can be adjusted at solve-time when prescribed values vary.
    Returns modified matrix (CSC) and a dict col_contribs: k -> (rows_arr, vals_arr)
    where original A[row, k] values are given (for row != k).
    Note: It is assumed dirichlet_idx are indices in 0..N-1.
    """
    A = A_csc.tolil()
    N = A.shape[0]
    col_contribs = {}
    for k in dirichlet_idx:
        # extract original column entries
        col = A[:, k].toarray().ravel()
        mask = np.nonzero(col)[0]
        mask = mask[mask != k]
        if mask.size > 0:
            rows = mask.copy()
            vals = col[mask].copy()
            col_contribs[k] = (rows, vals)
        else:
            col_contribs[k] = (np.array([], dtype=int), np.array([], dtype=float))
        # zero column except diagonal
        for i in range(A.shape[0]):
            if i != k and A[i, k] != 0:
                A[i, k] = 0.0
        # set row k to identity
        A.rows[k] = [k]
        A.data[k] = [1.0]
    return A.tocsc(), col_contribs

# ---------------------- Convective operator (vectorized-ish) ----------------------

def interp_v_to_u(v):
    # v has shape (Nx, Ny+1), want v at u-faces of shape (Nx+1, Ny)
    nx_v, ny_v = v.shape
    nx_u = nx_v + 1
    ny_u = ny_v - 1
    v_at_u = np.zeros((nx_u, ny_u))
    # interior u faces i=1..nx_u-2 average of four nearby v values
    # v indices used: v[i-1,j], v[i,j], v[i-1,j+1], v[i,j+1]
    v_at_u[1:-1, :] = 0.25 * (
        v[0:-1, 0:ny_u] + v[1:, 0:ny_u] + v[0:-1, 1:ny_u+1] + v[1:, 1:ny_u+1]
    )
    # leftmost u-face (inlet): average of v[0, j] and v[0, j+1]
    v_at_u[0, :] = 0.5 * (v[0, 0:ny_u] + v[0, 1:ny_u+1])
    # rightmost u-face (outlet): average of v[-1, j] and v[-1, j+1]
    v_at_u[-1, :] = 0.5 * (v[-1, 0:ny_u] + v[-1, 1:ny_u+1])
    return v_at_u


def interp_u_to_v(u):
    # u has shape (Nx+1, Ny), want u at v-faces of shape (Nx, Ny+1)
    nx_u, ny_u = u.shape
    nx_v = nx_u - 1
    ny_v = ny_u + 1
    u_at_v = np.zeros((nx_v, ny_v))
    # interior v faces j=1..ny_v-2 average of four nearby u values
    # u indices: u[i, j-1], u[i+1, j-1], u[i, j], u[i+1, j]
    u_at_v[:, 1:-1] = 0.25 * (
        u[0:-1, 0:ny_u-1] + u[1:, 0:ny_u-1] + u[0:-1, 1:ny_u] + u[1:, 1:ny_u]
    )
    # bottom v-face j=0
    u_at_v[:, 0] = 0.5 * (u[0:-1, 0] + u[1:, 0])
    # top v-face j=ny_v-1
    u_at_v[:, -1] = 0.5 * (u[0:-1, -1] + u[1:, -1])
    return u_at_v


def compute_convective(u, v, dx, dy):
    """
    Compute convective terms for u and v on their staggered faces using an upwind-like scheme.
    Returns conv_u (shape u) and conv_v (shape v).
    """
    # u convective
    v_at_u = interp_v_to_u(v)
    conv_u = np.zeros_like(u)
    # du/dx upwind
    du_forward = (np.roll(u, -1, axis=0) - u) / dx
    du_backward = (u - np.roll(u, 1, axis=0)) / dx
    # at left boundary backward uses same-index, at right boundary forward wraps; fix boundaries below
    # choose upwind based on u sign
    mask_pos = u >= 0
    mask_neg = ~mask_pos
    du_dx = np.where(mask_pos, du_backward, du_forward)
    # fix boundaries explicitly to one-sided where roll wrapped
    du_dx[0, :] = (u[1, :] - u[0, :]) / dx  # forward
    du_dx[-1, :] = (u[-1, :] - u[-2, :]) / dx  # backward
    # du/dy upwind using v_at_u sign
    du_up = (u - np.roll(u, 1, axis=1)) / dy
    du_down = (np.roll(u, -1, axis=1) - u) / dy
    du_dy = np.where(v_at_u >= 0, du_up, du_down)
    # fix top/bottom boundaries
    du_dy[:, 0] = (u[:, 1] - u[:, 0]) / dy
    du_dy[:, -1] = (u[:, -1] - u[:, -2]) / dy
    conv_u = u * du_dx + v_at_u * du_dy

    # v convective
    u_at_v = interp_u_to_v(u)
    conv_v = np.zeros_like(v)
    dv_forward = (np.roll(v, -1, axis=1) - v) / dy
    dv_backward = (v - np.roll(v, 1, axis=1)) / dy
    mask_pos_v = v >= 0
    dv_dy = np.where(mask_pos_v, dv_backward, dv_forward)
    dv_dy[:, 0] = (v[:, 1] - v[:, 0]) / dy
    dv_dy[:, -1] = (v[:, -1] - v[:, -2]) / dy

    dv_forward_x = (np.roll(v, -1, axis=0) - v) / dx
    dv_backward_x = (v - np.roll(v, 1, axis=0)) / dx
    dv_dx = np.where(u_at_v >= 0, dv_backward_x, dv_forward_x)
    dv_dx[0, :] = (v[0, :] - v[0, :]) / dx  # leftmost uses zero-gradient approx due to inlet Dirichlet
    dv_dx[-1, :] = (v[-1, :] - v[-2, :]) / dx

    conv_v = u_at_v * dv_dx + v * dv_dy

    if not (np.all(np.isfinite(conv_u)) and np.all(np.isfinite(conv_v))):
        raise RuntimeError('Non-finite convective term encountered')
    return conv_u, conv_v

# ---------------------- Inlet profile ----------------------

def inlet_u_profile(y, t):
    return np.sin(np.pi * y) * (np.sin(np.pi * t) + np.sin(3 * np.pi * t) + np.sin(5 * np.pi * t))

# ---------------------- Main Solver ----------------------

def main():
    # Domain and grid (passed explicitly)
    Lx = 2.0
    Ly = 1.0
    Nx = 64  # pressure cells
    Ny = 32
    dx = Lx / Nx
    dy = Ly / Ny

    # MAC shapes
    nx_u = Nx + 1
    ny_u = Ny
    nx_v = Nx
    ny_v = Ny + 1
    Nu = nx_u * ny_u
    Nv = nx_v * ny_v
    Np = Nx * Ny

    # Time stepping
    t_final = 1.0
    dt = 0.005
    Nt = int(np.ceil(t_final / dt))
    dt = t_final / Nt

    # Initialize fields
    u = np.zeros((nx_u, ny_u))
    v = np.zeros((nx_v, ny_v))
    p = np.zeros((Nx, Ny))

    # Build divergence operator and grad (adjoint) operator
    D, idx_u, idx_v, idx_p = build_divergence_matrix(Nx, Ny, dx, dy)
    G = -D.T  # discrete gradient operator mapping p_flat -> [dpdx; dpdy]

    # Build pressure Poisson matrix as A_p = D * G
    A_p = (D @ G).tocsc()

    # Build diffusion Laplacians for u and v
    A_u = build_laplace_matrix_standard(nx_u, ny_u, dx, dy)
    A_v = build_laplace_matrix_standard(nx_v, ny_v, dx, dy)

    # Build implicit diffusion matrices M = I - dt * visc * L
    Iu = sp.identity(A_u.shape[0], format='csc')
    Iv = sp.identity(A_v.shape[0], format='csc')
    Lu_mat = (Iu - dt * visc * A_u).tocsc()
    Lv_mat = (Iv - dt * visc * A_v).tocsc()

    # Build Dirichlet masks (True where Dirichlet is enforced)
    dirichlet_u = np.zeros((nx_u, ny_u), dtype=bool)
    dirichlet_u[0, :] = True  # inlet Dirichlet
    dirichlet_u[:, 0] = True  # bottom no-slip
    dirichlet_u[:, ny_u - 1] = True  # top no-slip

    dirichlet_v = np.zeros((nx_v, ny_v), dtype=bool)
    dirichlet_v[0, :] = True
    dirichlet_v[:, 0] = True
    dirichlet_v[:, ny_v - 1] = True

    dirichlet_p = np.zeros((Nx, Ny), dtype=bool)
    dirichlet_p[Nx - 1, :] = True  # pressure Dirichlet at outlet cells

    dirichlet_u_idx = [idx_u(i, j) for i in range(nx_u) for j in range(ny_u) if dirichlet_u[i, j]]
    dirichlet_v_idx = [idx_v(i, j) for i in range(nx_v) for j in range(ny_v) if dirichlet_v[i, j]]
    dirichlet_p_idx = [idx_p(i, j) for i in range(Nx) for j in range(Ny) if dirichlet_p[i, j]]

    # Prepare matrices by zeroing columns for Dirichlet DOFs and making identity rows; store column maps
    Lu_mat_mod, Lu_colmaps = prepare_matrix_with_dirichlet(Lu_mat, dirichlet_u_idx)
    Lv_mat_mod, Lv_colmaps = prepare_matrix_with_dirichlet(Lv_mat, dirichlet_v_idx)
    Lp_mat_mod, Lp_colmaps = prepare_matrix_with_dirichlet(A_p, dirichlet_p_idx)

    # Factorize once
    print('Factorizing matrices (one-time) ...')
    lu_u = spla.splu(Lu_mat_mod)
    lu_v = spla.splu(Lv_mat_mod)
    lu_p = spla.splu(Lp_mat_mod)

    # Unit-tests to verify Poisson operator and projection adjointness
    # Test 1: manufactured phi = x -> Laplacian zero ideally
    xc = (np.arange(Nx) + 0.5) * dx
    yc = (np.arange(Ny) + 0.5) * dy
    Xc, Yc = np.meshgrid(xc, yc, indexing='xy')
    phi_test1 = Xc  # linear in x
    phi_flat = phi_test1.ravel(order='F')
    res1 = A_p.dot(phi_flat)
    norm1 = np.linalg.norm(res1)

    # Test 2: phi = x^2 + y^2 -> Laplacian = 2 + 2 = 4 constant * cell area? Our discrete operator approximates Laplacian
    phi_test2 = Xc ** 2 + Yc ** 2
    phi2_flat = phi_test2.ravel(order='F')
    res2 = A_p.dot(phi2_flat)
    # Compute mean of res2 to inspect magnitude
    mean_res2 = np.mean(res2)
    norm2 = np.linalg.norm(res2 - mean_res2)

    print(f'Operator tests: ||A_p*phi(x)||_2 = {norm1:.3e} (should be small), mean(A_p*phi(x^2+y^2))={mean_res2:.3e}, ||res-mean||={norm2:.3e}')

    # Small projection test: create u_star and ensure projection reduces divergence
    # Random small u_star, compute rhs = D*u_star / dt, solve for phi_proj and compute divergence after correction
    np.random.seed(42)
    u_star_rand = 1e-3 * np.random.randn(nx_u, ny_u)
    v_star_rand = 1e-3 * np.random.randn(nx_v, ny_v)
    ustar_flat = u_star_rand.ravel(order='F')
    vstar_flat = v_star_rand.ravel(order='F')
    uv_flat = np.concatenate([ustar_flat, vstar_flat])
    div_star = D.dot(uv_flat).reshape((Nx, Ny), order='F')
    rhs_p_test = (div_star / dt).ravel(order='F')
    # no dirichlet p contributions since p_dirichlet values are zero -> just solve
    phi_proj_flat = lu_p.solve(rhs_p_test)
    g_flat = G.dot(phi_proj_flat)
    # corrected uv
    uv_corr = uv_flat - dt * g_flat
    div_after_proj = D.dot(uv_corr)
    max_div_proj = np.max(np.abs(div_after_proj))
    print(f'Projection test: max|div_after_proj| = {max_div_proj:.3e} (should be near machine precision)')

    # Coordinates for inlet profile (u faces y-locations)
    y_u = (np.arange(ny_u) + 0.5) * dy

    # Progress schedule
    progress_steps = set([int(Nt * frac / 10) for frac in range(1, 11)])
    if 0 in progress_steps:
        progress_steps.remove(0)

    t = 0.0
    print('Starting time loop: Nt =', Nt, ', dt =', dt)

    # Precompute indices arrays for adding Dirichlet contributions into RHS quickly
    # Each colmap: k -> (rows_arr, vals_arr) where original A[row,k] = vals

    for n in range(Nt):
        # Forcing at current time t on v-cells
        xv = (np.arange(nx_v) + 0.5) * dx
        yv = np.arange(ny_v) * dy
        Xv, Yv = np.meshgrid(xv, yv, indexing='xy')
        f_v = (-np.sin(np.pi * Xv) * np.sin(np.pi * Yv) * np.sin(np.pi * t)).T
        f_u = np.zeros_like(u)

        # Enforce Dirichlet BCs on current fields for convective evaluation
        u[0, :] = inlet_u_profile(y_u, t)
        u[:, 0] = 0.0
        u[:, ny_u - 1] = 0.0
        # outlet Neumann for u/v: keep ghost value same as interior for convective eval
        u[-1, :] = u[-2, :]

        v[0, :] = 0.0
        v[:, 0] = 0.0
        v[:, ny_v - 1] = 0.0
        v[-1, :] = v[-2, :]

        # Convective terms
        conv_u, conv_v = compute_convective(u, v, dx, dy)

        # Predictor RHS for velocities (flatten Fortran-order)
        rhs_u = (u - dt * conv_u + dt * f_u).ravel(order='F')
        rhs_v = (v - dt * conv_v + dt * f_v).ravel(order='F')

        # Before solving, impose Dirichlet values into RHS for Dirichlet DOFs (use boundary at t+dt)
        u_inlet_tnp1 = inlet_u_profile(y_u, t + dt)
        # Prepare a combined rhs vector for u and v solves separately
        # For Lu: we must add contributions from columns that were zeroed earlier: for each Dirichlet k with value val,
        # add original col entries A[i,k] * val to rhs[i]
        # First set the prescribed values array for u and v
        prescribed_u_vals = {idx_u(0, j): u_inlet_tnp1[j] for j in range(ny_u)}
        # top/bottom u zeros
        for i in range(nx_u):
            prescribed_u_vals[idx_u(i, 0)] = 0.0
            prescribed_u_vals[idx_u(i, ny_u - 1)] = 0.0
        prescribed_v_vals = {}
        for j in range(ny_v):
            prescribed_v_vals[idx_v(0, j)] = 0.0
        for i in range(nx_v):
            prescribed_v_vals[idx_v(i, 0)] = 0.0
            prescribed_v_vals[idx_v(i, ny_v - 1)] = 0.0

        # Add column contributions to rhs_u
        for k, val in prescribed_u_vals.items():
            rows_arr, vals_arr = Lu_colmaps[k]
            if rows_arr.size > 0 and val != 0.0:
                rhs_u[rows_arr] -= vals_arr * val  # move to RHS (note sign consistent with row-zeroing)
            rhs_u[k] = val  # enforced identity row RHS

        # Add column contributions to rhs_v
        for k, val in prescribed_v_vals.items():
            rows_arr, vals_arr = Lv_colmaps[k]
            if rows_arr.size > 0 and val != 0.0:
                rhs_v[rows_arr] -= vals_arr * val
            rhs_v[k] = val

        # Sanity check
        if not (np.all(np.isfinite(rhs_u)) and np.all(np.isfinite(rhs_v))):
            raise RuntimeError('Non-finite RHS for velocity solve encountered')

        # Solve implicit diffusion systems
        u_star_flat = lu_u.solve(rhs_u)
        v_star_flat = lu_v.solve(rhs_v)
        u_star = u_star_flat.reshape(u.shape, order='F')
        v_star = v_star_flat.reshape(v.shape, order='F')

        # Ensure Dirichlet BCs satisfied numerically
        u_star[0, :] = u_inlet_tnp1
        u_star[:, 0] = 0.0
        u_star[:, ny_u - 1] = 0.0
        v_star[0, :] = 0.0
        v_star[:, 0] = 0.0
        v_star[:, ny_v - 1] = 0.0

        # Compute divergence of intermediate velocity
        ustar_flat = u_star.ravel(order='F')
        vstar_flat = v_star.ravel(order='F')
        uvstar = np.concatenate([ustar_flat, vstar_flat])
        div_star = D.dot(uvstar).reshape((Nx, Ny), order='F')

        # Pressure Poisson RHS = div_star / dt (flattened)
        rhs_p = (div_star / dt).ravel(order='F')

        # Enforce Dirichlet pressure RHS entries to prescribed values contributions (pressure outlet is zero always)
        # For generality, get prescribed p vals (here zeros)
        prescribed_p_vals = {k: 0.0 for k in dirichlet_p_idx}
        # Add column contributions from Lp_colmaps
        for k, val in prescribed_p_vals.items():
            rows_arr, vals_arr = Lp_colmaps[k]
            if rows_arr.size > 0 and val != 0.0:
                rhs_p[rows_arr] -= vals_arr * val
            rhs_p[k] = val

        if not np.all(np.isfinite(rhs_p)):
            raise RuntimeError('Non-finite RHS for pressure solve encountered')

        # Solve Poisson
        phi_flat = lu_p.solve(rhs_p)
        phi = phi_flat.reshape((Nx, Ny), order='F')

        # Pressure is phi
        p = phi.copy()

        # Correct velocities using discrete gradient G (consistent with divergence operator)
        g_flat = G.dot(phi_flat)
        dpdx_flat = g_flat[:Nu]
        dpdy_flat = g_flat[Nu:]
        dpdx = dpdx_flat.reshape((nx_u, ny_u), order='F')
        dpdy = dpdy_flat.reshape((nx_v, ny_v), order='F')

        u = u_star - dt * dpdx
        v = v_star - dt * dpdy

        # Re-enforce velocity Dirichlet boundary values at t+dt exactly
        u[0, :] = u_inlet_tnp1
        u[:, 0] = 0.0
        u[:, ny_u - 1] = 0.0
        u[-1, :] = u[-2, :]
        v[0, :] = 0.0
        v[:, 0] = 0.0
        v[:, ny_v - 1] = 0.0
        v[-1, :] = v[-2, :]

        # Diagnostics
        uv_flat = np.concatenate([u.ravel(order='F'), v.ravel(order='F')])
        div_after = D.dot(uv_flat).reshape((Nx, Ny), order='F')
        max_div = np.max(np.abs(div_after))
        # interpolate to cell centers
        u_center = 0.5 * (u[0:Nx, :] + u[1:Nx + 1, :])
        v_center = 0.5 * (v[:, 0:Ny] + v[:, 1:Ny + 1])
        KE = 0.5 * np.sum(u_center ** 2 + v_center ** 2) * dx * dy

        t += dt
        if (n + 1) in progress_steps:
            pct = int(100 * (n + 1) / Nt)
            print(f'Progress: {pct}% | step {n+1}/{Nt} | t={t:.3f} | max_div={max_div:.3e} | KE={KE:.6e}')

    # Final diagnostics
    print(f'Finished. t={t:.6f} | max_div_final={max_div:.3e} | KE_final={KE:.6e}')

    # Interpolate u and v to cell centers for plotting
    u_center = 0.5 * (u[0:Nx, :] + u[1:Nx + 1, :])
    v_center = 0.5 * (v[:, 0:Ny] + v[:, 1:Ny + 1])

    # Coordinates for cell centers
    xc = (np.arange(Nx) + 0.5) * dx
    yc = (np.arange(Ny) + 0.5) * dy
    Xc, Yc = np.meshgrid(xc, yc, indexing='xy')

    # Save contour figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    cs0 = axes[0].contourf(Xc, Yc, u_center.T, levels=50, cmap='RdBu_r')
    axes[0].set_title('u at t=1.0')
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    fig.colorbar(cs0, ax=axes[0])

    cs1 = axes[1].contourf(Xc, Yc, v_center.T, levels=50, cmap='RdBu_r')
    axes[1].set_title('v at t=1.0')
    axes[1].set_xlabel('x')
    fig.colorbar(cs1, ax=axes[1])

    cs2 = axes[2].contourf(Xc, Yc, p.T, levels=50, cmap='RdBu_r')
    axes[2].set_title('p at t=1.0')
    axes[2].set_xlabel('x')
    fig.colorbar(cs2, ax=axes[2])

    plt.tight_layout()
    plt.savefig('ns_mac_solution.png', dpi=200)
    print('Saved figure to ns_mac_solution.png')


if __name__ == '__main__':
    main()




### runtime_outputs

#### Output block1

Factorizing matrices (one-time) ...
Starting time loop: Nt = 200 , dt = 0.005
<string>:413: RuntimeWarning: overflow encountered in square
<string>:144: RuntimeWarning: overflow encountered in scalar multiply
<string>:144: RuntimeWarning: invalid value encountered in scalar add
<string>:183: RuntimeWarning: overflow encountered in scalar multiply
<string>:183: RuntimeWarning: invalid value encountered in scalar add
Progress: 10% | step 20/200 | t=0.100 | max_div=nan | KE=nan
Progress: 20% | step 40/200 | t=0.200 | max_div=nan | KE=nan
Progress: 30% | step 60/200 | t=0.300 | max_div=nan | KE=nan
Progress: 40% | step 80/200 | t=0.400 | max_div=nan | KE=nan
Progress: 50% | step 100/200 | t=0.500 | max_div=nan | KE=nan
Progress: 60% | step 120/200 | t=0.600 | max_div=nan | KE=nan
Progress: 70% | step 140/200 | t=0.700 | max_div=nan | KE=nan
Progress: 80% | step 160/200 | t=0.800 | max_div=nan | KE=nan
Progress: 90% | step 180/200 | t=0.900 | max_div=nan | KE=nan
Progress: 100% | step 200/200 | t=1.000 | max_div=nan | KE=nan
Finished. t=1.000000 | max_div_final=nan | KE_final=nan
Saved figure to ns_mac_solution.png



#### Output block2

Factorizing matrices (one-time) ...
Starting time loop: Nt = 200 , dt = 0.005
Progress: 10% | step 20/200 | t=0.100 | max_div=2.847e+01 | KE=1.335428e-01
Progress: 20% | step 40/200 | t=0.200 | max_div=1.208e+01 | KE=2.819059e-01
Progress: 30% | step 60/200 | t=0.300 | max_div=1.377e+01 | KE=1.339708e-01
Progress: 40% | step 80/200 | t=0.400 | max_div=3.075e+00 | KE=8.851328e-02
Progress: 50% | step 100/200 | t=0.500 | max_div=7.703e+00 | KE=1.365910e-01
Progress: 60% | step 120/200 | t=0.600 | max_div=7.423e+00 | KE=1.384910e-01
Progress: 70% | step 140/200 | t=0.700 | max_div=4.407e+00 | KE=1.004279e-01
Progress: 80% | step 160/200 | t=0.800 | max_div=2.093e+01 | KE=1.682454e-01
Progress: 90% | step 180/200 | t=0.900 | max_div=1.352e+01 | KE=4.616995e-01
Progress: 100% | step 200/200 | t=1.000 | max_div=2.611e+01 | KE=3.344876e-01
Finished. t=1.000000 | max_div_final=2.611e+01 | KE_final=3.344876e-01
Saved figure to ns_mac_solution.png



#### Output block3

Factorizing matrices (one-time) ...
Operator tests: ||A_p*phi(x)||_2 = 1.958e+04 (should be small), mean(A_p*phi(x^2+y^2))=-2.036e+02, ||res-mean||=4.192e+04
Projection test: max|div_after_proj| = 6.509e+02 (should be near machine precision)
Starting time loop: Nt = 200 , dt = 0.005
Progress: 10% | step 20/200 | t=0.100 | max_div=2.847e+01 | KE=1.335428e-01
Progress: 20% | step 40/200 | t=0.200 | max_div=1.208e+01 | KE=2.819059e-01
Progress: 30% | step 60/200 | t=0.300 | max_div=1.377e+01 | KE=1.339708e-01
Progress: 40% | step 80/200 | t=0.400 | max_div=3.075e+00 | KE=8.851328e-02
Progress: 50% | step 100/200 | t=0.500 | max_div=7.703e+00 | KE=1.365910e-01
Progress: 60% | step 120/200 | t=0.600 | max_div=7.423e+00 | KE=1.384910e-01
Progress: 70% | step 140/200 | t=0.700 | max_div=4.407e+00 | KE=1.004279e-01
Progress: 80% | step 160/200 | t=0.800 | max_div=2.093e+01 | KE=1.682454e-01
Progress: 90% | step 180/200 | t=0.900 | max_div=1.352e+01 | KE=4.616995e-01
Progress: 100% | step 200/200 | t=1.000 | max_div=2.611e+01 | KE=3.344876e-01
Finished. t=1.000000 | max_div_final=2.611e+01 | KE_final=3.344876e-01
Saved figure to ns_mac_solution.png




### review_decision
revise

### review_comments
	Current Stage [A/3]
1. MAJOR

Category: runtime
Issue: Dirichlet rows are set in A (Laplacian) before forming implicit matrices; RHS not adjusted; leads to inconsistent linear systems and blow-up.

Recommendation: The root cause of the NaNs/overflows is incorrect handling of Dirichlet boundary rows when building and using the discrete Laplace operators. Do not permanently replace rows of the discrete Laplacian (A) with identity before forming the time-stepping matrix. Instead: (1) build the Laplacian L with the correct finite-difference stencil everywhere (no special identity rows), (2) form the implicit matrix M = I - dt * visc * L, and only then impose Dirichlet BCs on M by replacing the corresponding rows with identity and setting the RHS entries to the Dirichlet values. Similarly, for the pressure Poisson solve, enforce p=0 by replacing rows of the Poisson matrix with identity and setting the RHS to 0 at Dirichlet cell indices before solving. This will prevent spurious coupling and the blow-up you see. Also ensure the Laplacian sign is consistent with your discretization so that M is well-conditioned. (category: runtime, severity: MAJOR)


2. MAJOR

Category: accuracy
Issue: Poisson matrix rows for Dirichlet pressure exist but RHS is not overwritten to enforce p=0; leads to incorrect pressure and velocity corrections.

Recommendation: The current Poisson solve uses the assembled A_p (which had Dirichlet rows set earlier) but you never set the RHS to enforce p=0 on outlet cells — the solver therefore returns phi values at Dirichlet locations equal to whatever RHS was there (div/dt). Fix: after forming the Poisson matrix (A_p) and before factorizing/solving, set rows corresponding to Dirichlet pressure cells to identity and set the corresponding RHS entries to the prescribed pressure (zero). If you keep factorization, rebuild / re-factorize the Poisson matrix after enforcing Dirichlet rows (this is necessary once, before the time loop). Also check the solvability: ensure at least one Dirichlet cell exists to remove singularity. (category: accuracy, severity: MAJOR)


3. MAJOR

Category: structure
Issue: Inconsistent enforcement of Dirichlet BCs: rows altered incorrectly and RHS not adjusted; re-applying BCs after solve masks the real problem.

Recommendation: The way you enforce Dirichlet BCs by solving with modified matrices then re-applying BCs afterwards is inconsistent. For each implicit solve (u and v), enforce Dirichlet by: (a) forming RHS including any contributions from known Dirichlet values (or simply set RHS entries at Dirichlet DOFs to the Dirichlet values), (b) set corresponding rows in the matrix to identity, (c) solve, which yields values that satisfy BCs exactly without needing to reapply them. Your current approach (modify A early and not adjust RHS) results in incorrect linear systems and is a primary source of the instability. Also do not rely on setting the last face equal to the previous face to realize Neumann BCs — handle Neumann consistently in the operator or in the RHS. (category: structure, severity: MAJOR)


4. MINOR

Category: accuracy
Issue: Convective terms have fragile boundary-handling and are slow (nested Python loops). This can produce incorrect boundary derivatives and contribute to instability.

Recommendation: The convective term routines are correct in spirit but fragile at boundaries and very slow (pure Python loops). Two actions: (1) audit and fix boundary interpolations carefully — several fall-back branches (e.g. du_dx at i==0) use expressions that may produce zeros or spurious large slopes; ensure upwind stencils use well-defined neighbors and do not read out-of-bounds or take means of empty arrays; add assertions/checks that convective fields are finite. (2) For performance and stability, vectorize these computations (NumPy) or at least precompute interpolated face velocities in arrays. Finally, the explicit convective step imposes a CFL constraint — if you keep it explicit, reduce dt or use an IMEX scheme for convection if you need larger dt. (category: accuracy, severity: MINOR)


5. MINOR

Category: structure
Issue: Lack of sanity checks and improper matrix factorization timing; missing enforcement of BCs in RHS before factorization/solve.

Recommendation: Add runtime sanity checks and propagate BCs into RHS prior to solves: after creating RHS vectors but before calling LU solve, check np.isfinite(rhs). If any non-finite entries are found, raise an informative error. Also, when replacing rows for Dirichlet BCs, set the diagonal to 1.0 and RHS to the desired boundary value. Re-factorize the implicit matrices after the rows are modified (do this once outside the time loop) so solves use the correct operator. For the pressure Poisson, either factorize the final Dirichlet-enforced matrix once or switch to an iterative solver (e.g. conjugate gradient with proper preconditioning) to avoid expensive sparse LU for large grids. Finally, test the solver with a manufactured solution or small grid where analytic/expected behavior (e.g. divergence reduction after projection) can be verified. (category: structure, severity: MINOR)



	Current Stage [B/3]
1. MAJOR

Category: accuracy
Issue: Projection fails to enforce incompressibility: reported max |div| at final time is O(10)--O(10^1), not near machine precision. The Chorin projection is therefore not removing divergence and the velocity/pressure fields are not physically consistent.

Recommendation: Investigate and fix the Poisson/Laplacian assembly and pressure-gradient discretization. As immediate checks: (1) verify that A_p * phi = rhs works for simple manufactured phi (e.g. phi = x, phi = x^2+y^2) by computing residuals; (2) after solving the Poisson, compute the discrete divergence of u_star - dt*grad(phi) analytically for a simple manufactured u_star to ensure divergence goes to zero. The goal is to get max|div| ~ 1e-10. Do not proceed until the projection properly enforces incompressibility.


2. MAJOR

Category: accuracy
Issue: build_laplace_matrix implements boundary and Neumann handling in an ad-hoc way (mirror-approximation via adding coefficients to other neighbors and subtracting stencil weights into diag). This nonstandard treatment likely yields an incorrect discrete Laplacian, breaking both diffusion solves and the pressure Poisson.

Recommendation: Replace build_laplace_matrix with a clean, standard assembly: for interior nodes use the 5-point stencil; for Dirichlet nodes either (a) remove row/column and modify RHS (preferred for symmetric solvers), or (b) set row to identity and zero corresponding column entries (to preserve factorization behavior). For Neumann BCs, implement ghost-point elimination or explicitly modify the stencil to enforce ∂φ/∂n = g by shifting the finite-difference formula (i.e. use one-sided discretization or substitute ghost value from Neumann condition). Document and unit-test these cases (apply operator to known fields).


3. MAJOR

Category: accuracy
Issue: Imposition of Dirichlet boundary conditions is inconsistent: the code only replaces rows with identity but does not zero corresponding columns or update RHS of neighboring rows to account for fixed values. This can produce incorrect coupling in the linear systems and contaminate interior equations.

Recommendation: When imposing Dirichlet BCs in a matrix solve, either: (1) eliminate DOFs analytically (reduce system and modify RHS), or (2) set row to identity and also zero the column entries (except the diagonal) and set RHS to the prescribed value. If using (2) confirm that the RHS of all affected neighbor rows has been adjusted so the fixed value's contribution is not double-counted or lost. Add unit tests to confirm operator+BC correctness.


4. MAJOR

Category: accuracy
Issue: Pressure gradient to face interpolation (grad_p_to_u, grad_p_to_v) uses one-sided/zero-gradient approximations at inlet faces (dpdx[0]=0) and constructs rightmost face gradient using a hard-coded 0 ghost value. This is inconsistent with how pressure Dirichlet cells were defined and can introduce large errors in the correction step.

Recommendation: Re-derive the mapping from cell-centered pressure to face-centered gradient for the exact MAC arrangement used. Ensure consistency: if p is defined at cell centers and outlet pressure is given at a cell adjacent to the boundary, use the correct ghost/one-sided formula that matches the Poisson discretization. Verify the discrete identity div(u^{n+1}) = 0 by plugging in the discrete grad and div operators; make them adjoint where required (discrete grad = -discrete div^T up to scaling).


5. MINOR

Category: structure
Issue: Convective term functions are implemented as explicit nested Python loops with many per-cell Python checks and branch logic. This is very slow and complicates correctness audits (and may hide indexing errors).

Recommendation: Vectorize convective term calculations with NumPy where possible, or at least isolate interpolation formulas into small, tested helper functions. Add small unit tests (constant velocity / linear velocity fields) to verify convective operator behavior and sign convention. While performance is secondary, cleaner vectorized code will make it easier to verify correctness.



	Current Stage [C/3]
1. MAJOR

Category: accuracy
Issue: Projection test and operator tests fail: A_p*phi and the projection residual are large (Operator test ||A_p*phi(x)||_2 ~ 2e4, projection max|div_after_proj| ~ 6.5e2). The solver does not enforce incompressibility and produces large divergence (max_div ~ O(10)-O(10^1) at runtime).

Recommendation: Stop further refinement until the discrete Poisson/gradient/divergence operators are verified. Suggested checks and fixes: (1) Verify assembly/sign convention of G and A_p: compute A_p explicitly and check symmetry: A_p + A_p.T should be near zero if using the expected sign convention; otherwise adjust G = +D.T or remove the extra negative sign so that A_p approximates the standard positive-definite (or sign-consistent) discrete Laplacian. (2) Compare A_p with a separately assembled 5-point Laplacian on pressure cells (use same spacing and the same flattening order) — their action on polynomials (phi=x, phi=y, phi=x^2+y^2) should match theoretical expectations (linear -> zero). (3) After fixes, re-run the projection unit test: build a random small u_star/v_star, solve for phi and confirm max|div_after_proj| is near machine precision (e.g. <= 1e-10). Do not proceed until these pass.


2. MAJOR

Category: structure
Issue: build_laplace_matrix_standard uses an unconventional 'ghost = interior' handling and constructs diagonals by subtracting neighbor contributions in a way that yields negative diagonal entries and possibly wrong sign conventions for the Laplacian. This makes reasoning about signs in implicit solves and in Poisson assembly hard.

Recommendation: Replace build_laplace_matrix_standard with an explicit, standard finite-difference 5-point stencil assembly: set A[k,k] = -2*(1/dx^2 + 1/dy^2) (or the appropriate sign consistent with your diffusion convention) and A[k, neighbor] = 1/dx^2 or 1/dy^2 for valid neighbors. Handle boundary conditions explicitly (Dirichlet rows or Neumann stencil modifications) instead of hiding them via 'ghost=interior' logic. After re-writing, check eigenvalues of the diffusion operator to verify sign and scale before using it in forming I - dt*visc*L.


3. MAJOR

Category: structure
Issue: Convective term implementation uses np.roll which implicitly wraps data and creates periodic ghosting; then the code attempts to patch boundary rows. This pattern is fragile and can silently introduce incorrect differences at boundaries (especially on staggered grid), and the upwind logic mixes arrays of different staggering with several transposes/averages that are not validated.

Recommendation: Avoid np.roll for boundary differencing. Compute one-sided differences explicitly at domain boundaries and use properly-interpolated face velocities without roll-induced wrapping. Add small unit tests for convective operators: feed a known polynomial velocity field and verify convective operator respects expected behavior (e.g. for linear fields convective derivative should be zero). Also check that interp_u_to_v and interp_v_to_u produce expected values on simple analytic fields (constant and linear tests).


4. MAJOR

Category: accuracy
Issue: Dirichlet enforcement machinery (prepare_matrix_with_dirichlet plus stored column contributions) is complex and error-prone. It is not confirmed that the stored column maps use the same indexing and sign conventions as the RHS adjustments; misuse will corrupt RHS leading to incorrect solves (this may explain large residuals).

Recommendation: Simplify and validate Dirichlet enforcement: (a) build matrices with known boundary-handling rules (set Dirichlet rows to identity and explicitly subtract the full original column contribution from RHS using a vectorized approach), (b) add unit tests: for a matrix A and a vector x with known Dirichlet values on mask, check that solving the modified system returns the original x. Also print norms of the original and modified matrix rows/columns for a few representative Dirichlet DOFs to verify you saved the correct column entries. Consider assembling velocity diffusion matrices while applying BCs directly (i.e., never assemble ghost couplings) to avoid needing to store and reapply column contributions.


5. MINOR

Category: structure
Issue: Several indexing/shape operations are fragile (use of .T on forcing f_v, mixed meshgrid indexing, and many ravel(order='F') usages). These can easily introduce shape mismatches or transposed fields that silently break operator adjointness.

Recommendation: Make indexing and array shapes explicit and consistent: prefer indexing='ij' for meshgrid when using Fortran ordering; remove unnecessary .T usages and validate shapes with assertions (e.g. assert f_v.shape == v.shape). Add small consistency assertions after building D and G: for a random phi compute finite-difference gradient and compare to G*phi for a few entries; likewise ensure D and G are adjoints up to sign: np.linalg.norm(D.toarray() + G.T.toarray()) should be small (up to floating rounding) given your sign convention.





### iteration_count
3

### rev_image_description
	Current Stage [A/3]


	Current Stage [B/3]


	Current Stage [C/3]




