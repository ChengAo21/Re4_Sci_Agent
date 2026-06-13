### 🛌🏻 Re4gent with Multi-Modal Review

#### Configs 🏷️:
MaxCount: 2, MultiModalReview: False, AlternativeNum: 2, MaxComments: 5, ConsModule_ON: False, CUT_OUTEXT: 2000, FROM_Seed: True

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
Use indexing='ij' for meshgrid and order='F' for all ravel/reshape operations to ensure consistent shapes and correct coordinate alignment.
Fix the pressure at one single node to remove the nullspace singularity.
Use the Immersed Boundary Method with volume penalization on a Staggered (MAC) Grid. Treat the penalization term implicitly in the momentum matrix diagonal to avoid restrictive time-steps.
Rely solely on the penalization term to enforce the no-slip condition on the cylinder surface.
Focus on generating runnable code with reasonable physical results, do not over-optimize.
Print concise progress information ONLY every 10% of total steps.


### expanded_prob


### solution_plans


### technical_spec
	Current Stage [A/2]

Grid & discretization:
- Uniform rectangular MAC grid: pressure at cell centers (Nx x Ny), u on vertical faces ((Nx+1) x Ny), v on horizontal faces (Nx x (Ny+1)).
- Cylinder centered at (0,0) with radius 0.5 represented by binary masks chi_u and chi_v at face locations.

Operators & linear solves:
- Discrete Laplacians constructed for u, v and p (Neumann approximations at outer boundaries). Dirichlet inlet (u=1) and v=0 at top/bottom enforced by replacing rows with identity.
- Velocity linear systems: Au = I - dt*visc*L_u + dt*(chi/eta) diag, Av likewise; penalization diag added implicitly and matrices factorized once (scipy.splu) and reused.
- Pressure Poisson: L_p built and pinned at one node to remove nullspace; factorized once.

Time-stepping & nonlinearity:
- Semi-implicit scheme: explicit AB2 for convection (N_extrap = 1.5*N^n - 0.5*N^{n-1}).
- Implicit backward-Euler for viscous + penalization via Au/Av solves to get tentative velocities u*, v*.
- Projection: solve Poisson for p correction and project velocities to divergence-free field.
- Boundary enforcement: Dirichlet inlet u and v top/bottom, zero-gradient (copying) used for outlet and for free-slip u boundary.

Acceleration & implementation choices:
- Matrices Au, Av, Lp built and factorized once outside time loop; solves reuse LU factors to reduce cost.
- Penalization treated implicitly by including dt*(chi/eta) on diagonal.
- Convection computed with central differences and simple interpolations between staggered grids (keeps code compact and robust).

Diagnostics & output:
- Prints progress every 10% of total steps with concise KE and max-velocity metrics.
- Saves final velocity magnitude contour (RdBu_r) to 'velocity_magnitude.png' and overlays cylinder boundary.

Limitations & notes:
- Boundary treatments use simple mirrored/one-sided approximations for Neumann conditions; this favors robustness over strict conservation.
- Penalization parameter eta and grid resolution influence enforcement accuracy near the cylinder; values chosen provide stable, qualitative results for Re=100.
- The solver aims at clarity and reproducibility: matrices are assembled once and reused.


	Current Stage [B/2]
Semi-implicit MAC solver with volume-penalized immersed boundary
- Grid: staggered MAC (u on vertical faces (Nx+1,Ny), v on horizontal faces (Nx,Ny+1), p at centers (Nx,Ny)).
- Ordering: column-major / Fortran (order='F') used for all flatten/reshape and meshgrid indexing='ij'. This ensures DOF indexing (idx_u/idx_v/idx_p) matches sparse matrix assembly and penalization diagonals.
- Operators assembled once: Laplacians for u, v, p. Pressure pinned at interior cell to remove nullspace.
- Time stepping: AB2 extrapolation for convection (explicit), backward-Euler implicit viscous + implicit diagonal penalization (Au/Av solves), projection via Poisson (Lp).
- Penalization: included implicitly on diagonal (dt * chi / eta) and the matrices factorized once (splu) and reused.
- Convection interpolation vectorized using array slices (no nested Python loops) for better performance.
- Outputs: concise progress messages every 10% of steps, final velocity magnitude contour saved to 'velocity_magnitude.png'.




### python_codes

#### Script block1:

```python
#!/usr/bin/env python3
"""
2D incompressible Navier-Stokes around a cylinder using a MAC grid
with immersed boundary volume penalization. Semi-implicit time stepping:
- explicit AB2 for convection
- implicit (backward-Euler) for viscous + implicit diagonal penalization
- pressure projection (Poisson with one pinned node)

Requirements from problem statement are followed:
- global visc defined (1/Re)
- penalization implicit via diagonal term
- staggered MAC grid
- precompute matrices and reuse factorizations
- print progress every 10% of steps
- save final contour plot (RdBu_r) and mark cylinder

This code is written for clarity and reasonable computational cost; it favors
robustness and simplicity over maximum accuracy/efficiency.

"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import floor

# Global physical parameter
Re = 100.0
visc = 1.0 / Re   # explicit global visc as required

# Penalization parameter (smaller = stronger enforcement)
eta = 1e-3


def build_grid(x1=3.0, x2=8.0, y1=3.0, y2=3.0, Nx=220, Ny=120):
    """Build domain and MAC grid coordinates.
    Domain: x in [-x1, x2], y in [-y1, y2] with cylinder centered at (0,0).
    Returns dx, dy and coordinate arrays for u (vertical faces), v (horizontal faces), p (cells).
    """
    Lx = x1 + x2
    Ly = y1 + y2
    dx = Lx / Nx
    dy = Ly / Ny

    # pressure cell centers
    x_p = np.linspace(-x1 + dx/2, x2 - dx/2, Nx)
    y_p = np.linspace(-y1 + dy/2, y2 - dy/2, Ny)

    # u at vertical faces: x faces i=0..Nx, y centers j=0..Ny-1
    x_u = np.linspace(-x1, x2, Nx+1)
    y_u = np.linspace(-y1 + dy/2, y2 - dy/2, Ny)

    # v at horizontal faces: x centers i=0..Nx-1, y faces j=0..Ny
    x_v = np.linspace(-x1 + dx/2, x2 - dx/2, Nx)
    y_v = np.linspace(-y1, y2, Ny+1)

    return dx, dy, x_u, y_u, x_v, y_v, x_p, y_p


def idx_u(i, j, Nx, Ny):
    return j * (Nx + 1) + i


def idx_v(i, j, Nx, Ny):
    return j * Nx + i


def idx_p(i, j, Nx, Ny):
    return j * Nx + i


def build_cylinder_masks(x_u, y_u, x_v, y_v, radius=0.5):
    """Build boolean masks chi_u and chi_v for penalization (1 inside cylinder).
    Cylinder is centered at (0,0).
    """
    Xu, Yu = np.meshgrid(x_u, y_u, indexing='xy')
    Xv, Yv = np.meshgrid(x_v, y_v, indexing='xy')

    chi_u = ((Xu**2 + Yu**2) <= radius**2).astype(float)
    chi_v = ((Xv**2 + Yv**2) <= radius**2).astype(float)

    return chi_u, chi_v


def build_laplacian_u(Nx, Ny, dx, dy, inlet_u_indices, chi_u_shape):
    """Build sparse Laplacian (second-order) for u-grid (Nu = (Nx+1)*Ny).
    Neumann boundaries (zero normal derivative) are approximated using mirrored ghost points.
    Rows corresponding to inlet Dirichlet indices will be overwritten later by identity rows.
    """
    Nu = (Nx + 1) * Ny
    data = []
    rows = []
    cols = []

    for j in range(Ny):
        for i in range(Nx + 1):
            row = idx_u(i, j, Nx, Ny)
            # Dirichlet inlet (we will set rows to identity later outside)
            if row in inlet_u_indices:
                rows.append(row); cols.append(row); data.append(1.0)
                continue

            center = 0.0
            # x-direction neighbors
            # left neighbor
            if i - 1 >= 0:
                left = idx_u(i - 1, j, Nx, Ny)
                rows.append(row); cols.append(left); data.append(1.0 / dx**2)
                center += -1.0 / dx**2
            else:
                # left boundary: Neumann mirrored ghost -> contributes as if coefficient 2 for right neighbor
                # This case corresponds to the leftmost u-face, but leftmost is inlet Dirichlet usually.
                pass

            # right neighbor
            if i + 1 <= Nx:
                right = idx_u(i + 1, j, Nx, Ny)
                rows.append(row); cols.append(right); data.append(1.0 / dx**2)
                center += -1.0 / dx**2
            else:
                pass

            # y-direction neighbors
            if j - 1 >= 0:
                down = idx_u(i, j - 1, Nx, Ny)
                rows.append(row); cols.append(down); data.append(1.0 / dy**2)
                center += -1.0 / dy**2
            else:
                # Neumann at bottom (mirror)
                # treat by adding an extra -1/dy^2 and doubling top neighbor when exists
                center += 0.0
            if j + 1 <= Ny - 1:
                up = idx_u(i, j + 1, Nx, Ny)
                rows.append(row); cols.append(up); data.append(1.0 / dy**2)
                center += -1.0 / dy**2
            else:
                center += 0.0

            # center
            rows.append(row); cols.append(row); data.append(center)

    L = sp.csr_matrix((data, (rows, cols)), shape=(Nu, Nu))
    return L


def build_laplacian_v(Nx, Ny, dx, dy, v_dirichlet_indices):
    """Build sparse Laplacian for v-grid (Nv = Nx*(Ny+1)).
    Dirichlet rows (top/bottom/left) will be identity overwritten as needed.
    Neumann approximated by mirrored ghost.
    """
    Nv = Nx * (Ny + 1)
    data = []
    rows = []
    cols = []

    for j in range(Ny + 1):
        for i in range(Nx):
            row = idx_v(i, j, Nx, Ny)
            if row in v_dirichlet_indices:
                rows.append(row); cols.append(row); data.append(1.0)
                continue

            center = 0.0
            # x-direction
            if i - 1 >= 0:
                left = idx_v(i - 1, j, Nx, Ny)
                rows.append(row); cols.append(left); data.append(1.0 / dx**2)
                center += -1.0 / dx**2
            if i + 1 <= Nx - 1:
                right = idx_v(i + 1, j, Nx, Ny)
                rows.append(row); cols.append(right); data.append(1.0 / dx**2)
                center += -1.0 / dx**2

            # y-direction
            if j - 1 >= 0:
                down = idx_v(i, j - 1, Nx, Ny)
                rows.append(row); cols.append(down); data.append(1.0 / dy**2)
                center += -1.0 / dy**2
            if j + 1 <= Ny:
                up = idx_v(i, j + 1, Nx, Ny)
                rows.append(row); cols.append(up); data.append(1.0 / dy**2)
                center += -1.0 / dy**2

            rows.append(row); cols.append(row); data.append(center)

    L = sp.csr_matrix((data, (rows, cols)), shape=(Nv, Nv))
    return L


def build_laplacian_p(Nx, Ny, dx, dy, pinned_idx=0):
    """Pressure Poisson Laplacian on pressure cell centers (Np = Nx*Ny).
    Homogeneous Neumann approximated by mirrored ghost, then one pressure node is pinned to remove nullspace.
    """
    Np = Nx * Ny
    data = []
    rows = []
    cols = []

    for j in range(Ny):
        for i in range(Nx):
            row = idx_p(i, j, Nx, Ny)
            if row == pinned_idx:
                rows.append(row); cols.append(row); data.append(1.0)
                continue

            center = 0.0
            # left
            if i - 1 >= 0:
                left = idx_p(i - 1, j, Nx, Ny)
                rows.append(row); cols.append(left); data.append(1.0 / dx**2)
                center += -1.0 / dx**2
            else:
                # Neumann mirror -> add contribution to right neighbor effectively
                center += 0.0
            # right
            if i + 1 <= Nx - 1:
                right = idx_p(i + 1, j, Nx, Ny)
                rows.append(row); cols.append(right); data.append(1.0 / dx**2)
                center += -1.0 / dx**2
            else:
                center += 0.0
            # down
            if j - 1 >= 0:
                down = idx_p(i, j - 1, Nx, Ny)
                rows.append(row); cols.append(down); data.append(1.0 / dy**2)
                center += -1.0 / dy**2
            else:
                center += 0.0
            # up
            if j + 1 <= Ny - 1:
                up = idx_p(i, j + 1, Nx, Ny)
                rows.append(row); cols.append(up); data.append(1.0 / dy**2)
                center += -1.0 / dy**2
            else:
                center += 0.0

            rows.append(row); cols.append(row); data.append(center)

    L = sp.csr_matrix((data, (rows, cols)), shape=(Np, Np))
    return L


def divergence(u, v, Nx, Ny, dx, dy):
    """Compute discrete divergence at pressure cell centers from face-centered u and v.
    u shape: (Nx+1, Ny), v shape: (Nx, Ny+1)
    returns div shape (Nx, Ny)
    """
    # divergence = (u_{i+1,j} - u_{i,j})/dx + (v_{i,j+1} - v_{i,j})/dy
    div = (u[1:Nx+1, :] - u[0:Nx, :]) / dx + (v[:, 1:Ny+1] - v[:, 0:Ny]) / dy
    return div


def gradient_p_to_u(p, Nx, Ny, dx, dy):
    """Compute dp/dx at u face locations from cell-centered pressure p (shape Nx x Ny).
    Returns array of shape (Nx+1, Ny)
    For faces at i=0..Nx: dp/dx at face between cells i-1 and i: (p_i - p_{i-1})/dx
    For boundaries where one side missing use one-sided difference (Neumann->0 gradient assumed if pinned)
    """
    dpdx = np.zeros((Nx + 1, Ny))
    # interior faces i=1..Nx-1
    dpdx[1:Nx, :] = (p[1:Nx, :] - p[0:Nx-1, :]) / dx
    # left boundary face i=0: use one-sided (p0 - p0) -> 0
    dpdx[0, :] = (p[0, :] - p[0, :]) / dx
    # right face i=Nx: use one-sided
    dpdx[Nx, :] = (p[Nx - 1, :] - p[Nx - 1, :]) / dx
    return dpdx


def gradient_p_to_v(p, Nx, Ny, dx, dy):
    """Compute dp/dy at v face locations from p.
    Returns array shape (Nx, Ny+1)
    """
    dpdy = np.zeros((Nx, Ny + 1))
    dpdy[:, 1:Ny] = (p[:, 1:Ny] - p[:, 0:Ny-1]) / dy
    dpdy[:, 0] = (p[:, 0] - p[:, 0]) / dy
    dpdy[:, Ny] = (p[:, Ny - 1] - p[:, Ny - 1]) / dy
    return dpdy


def interp_to_u_for_v(u, i, j):
    # placeholder if more advanced interpolation needed
    return None


def compute_convective(u, v, Nx, Ny, dx, dy):
    """Compute convective terms N_u and N_v on u and v grids respectively.
    Uses central differences with simple interpolation between grids.
    Returns N_u shape (Nx+1, Ny), N_v shape (Nx, Ny+1).
    """
    # N_u = u * du/dx + v_at_u * du/dy
    # du/dx at u faces
    du_dx = np.zeros_like(u)
    # central differences in x
    du_dx[1:Nx, :] = (u[2:Nx+1, :] - u[0:Nx-1, :]) / (2 * dx)
    # boundaries: use one-sided
    du_dx[0, :] = (u[1, :] - u[0, :]) / dx
    du_dx[Nx, :] = (u[Nx, :] - u[Nx - 1, :]) / dx

    # du/dy at u faces
    du_dy = np.zeros_like(u)
    du_dy[:, 1:Ny - 1] = (u[:, 2:Ny] - u[:, 0:Ny - 2]) / (2 * dy)
    # one-sided near boundaries
    if Ny >= 2:
        du_dy[:, 0] = (u[:, 1] - u[:, 0]) / dy
        du_dy[:, Ny - 1] = (u[:, Ny - 1] - u[:, Ny - 2]) / dy

    # interpolate v to u locations: v at (i-0.5,j+0.0) approximated by average of surrounding v
    v_at_u = np.zeros_like(u)
    # v shape (Nx, Ny+1), u shape (Nx+1, Ny)
    # for u at i=0..Nx, j=0..Ny-1, take average of v[i-1,j+1] & v[i,j+1] ... approximate carefully
    # We'll average the two nearest v-values when available
    for j in range(Ny):
        for i in range(Nx + 1):
            vals = []
            if i - 1 >= 0:
                vals.append(v[i - 1, j])
            if i <= Nx - 1:
                vals.append(v[i, j])
            if vals:
                v_at_u[i, j] = 0.5 * sum(vals)
            else:
                v_at_u[i, j] = 0.0

    N_u = u * du_dx + v_at_u * du_dy

    # N_v = u_at_v * dv/dx + v * dv/dy
    dv_dy = np.zeros_like(v)
    dv_dy[:, 1:Ny] = (v[:, 2:Ny + 1] - v[:, 0:Ny - 1]) / (2 * dy)
    dv_dy[:, 0] = (v[:, 1] - v[:, 0]) / dy
    dv_dy[:, Ny] = (v[:, Ny] - v[:, Ny - 1]) / dy

    dv_dx = np.zeros_like(v)
    dv_dx[1:Nx - 1, :] = (v[2:Nx, :] - v[0:Nx - 2, :]) / (2 * dx)
    # one-sided
    if Nx >= 2:
        dv_dx[0, :] = (v[1, :] - v[0, :]) / dx
        dv_dx[Nx - 1, :] = (v[Nx - 1, :] - v[Nx - 2, :]) / dx

    # interpolate u to v locations
    u_at_v = np.zeros_like(v)
    for j in range(Ny + 1):
        for i in range(Nx):
            vals = []
            if j - 1 >= 0:
                vals.append(u[i + 1, j - 1])
                vals.append(u[i, j - 1])
            if j <= Ny - 1:
                vals.append(u[i + 1, j])
                vals.append(u[i, j])
            if vals:
                u_at_v[i, j] = sum(vals) / len(vals)
            else:
                u_at_v[i, j] = 0.0

    N_v = u_at_v * dv_dx + v * dv_dy

    return N_u, N_v


def flatten(arr):
    return arr.ravel()


def unflatten_u(vec, Nx, Ny):
    return vec.reshape((Ny, Nx + 1)).T


def unflatten_v(vec, Nx, Ny):
    return vec.reshape((Ny + 1, Nx)).T


def unflatten_p(vec, Nx, Ny):
    return vec.reshape((Ny, Nx)).T


def apply_dirichlet_u(vec_u, u_inlet_value, Nx, Ny, inlet_indices):
    """Set inlet u-values to prescribed u_inlet_value in flattened vector.
    inlet_indices: set/list of flattened DOF indices corresponding to inlet faces.
    """
    for idx in inlet_indices:
        vec_u[idx] = u_inlet_value
    return vec_u


def apply_dirichlet_v(vec_v, Nx, Ny, v_zero_indices):
    for idx in v_zero_indices:
        vec_v[idx] = 0.0
    return vec_v


def main():
    # Domain and grid
    x1 = 3.0
    x2 = 8.0
    y1 = 3.0
    y2 = 3.0
    # grid resolution (moderate)
    Nx = 220
    Ny = 120

    dx, dy, x_u, y_u, x_v, y_v, x_p, y_p = build_grid(x1, x2, y1, y2, Nx=Nx, Ny=Ny)

    # construct masks for cylinder penalization
    chi_u, chi_v = build_cylinder_masks(x_u, y_u, x_v, y_v, radius=0.5)

    # Flattened sizes
    Nu = (Nx + 1) * Ny
    Nv = Nx * (Ny + 1)
    Np = Nx * Ny

    # inlet u-face indices (left boundary) flattened
    inlet_u_indices = [idx_u(0, j, Nx, Ny) for j in range(Ny)]
    inlet_u_set = set(inlet_u_indices)

    # v Dirichlet indices (top/bottom and inlet v=0)
    v_zero_indices = set()
    # top j = Ny -> v index j = Ny
    for i in range(Nx):
        v_zero_indices.add(idx_v(i, 0, Nx, Ny))   # bottom v=0
        v_zero_indices.add(idx_v(i, Ny, Nx, Ny))  # top v=0
    # inlet vertical velocity v at left boundary is naturally zero (inflow horizontal)
    # For safety we won't force all left/right v to zero; we used only top/bottom

    # Build Laplacians
    print('Building operators...')
    L_u = build_laplacian_u(Nx, Ny, dx, dy, inlet_u_set, chi_u.shape)
    L_v = build_laplacian_v(Nx, Ny, dx, dy, v_zero_indices)
    # build pressure Laplacian and pin one node at (0,0)
    pinned_idx = 0
    L_p = build_laplacian_p(Nx, Ny, dx, dy, pinned_idx=pinned_idx)

    # Time-stepping parameters
    T = 20.0
    # choose dt via convective CFL: dt <= CFL * min(dx,dy) / U_max
    Umax = 1.0
    CFL = 0.45
    dt = CFL * min(dx, dy) / Umax
    nsteps = int(floor(T / dt))
    dt = T / nsteps  # adjust dt to fit exactly T

    print(f'Nx={Nx} Ny={Ny} dx={dx:.4f} dy={dy:.4f} dt={dt:.5f} steps={nsteps}')

    # Build constant matrices for velocity solves: Au = I - dt*visc*L_u + dt*(chi/eta) diag
    # Au and Av are sparse (Nu x Nu) and (Nv x Nv)
    # Build identity
    Iu = sp.eye(Nu, format='csr')
    Iv = sp.eye(Nv, format='csr')

    # diag penalization arrays
    chi_u_flat = chi_u.ravel()
    chi_v_flat = chi_v.ravel()

    # Build penalization diagonal matrices
    Pen_u = sp.diags(dt * (chi_u_flat / eta), format='csr')
    Pen_v = sp.diags(dt * (chi_v_flat / eta), format='csr')

    Au = Iu - dt * visc * L_u + Pen_u
    Av = Iv - dt * visc * L_v + Pen_v

    # Enforce Dirichlet (inlet) rows in Au to identity to fix u=1 at inlet faces
    Au = Au.tolil()
    for idx in inlet_u_indices:
        Au.rows[idx] = [idx]
        Au.data[idx] = [1.0]
    Au = Au.tocsc()

    # Enforce Dirichlet rows for v at top/bottom
    Av = Av.tolil()
    for idx in v_zero_indices:
        Av.rows[idx] = [idx]
        Av.data[idx] = [1.0]
    Av = Av.tocsc()

    # Factorize matrices for fast repeated solves
    print('Factorizing velocity matrices...')
    Au_fac = spla.splu(Au)
    Av_fac = spla.splu(Av)

    # Factorize pressure Laplacian (pinning included) once
    print('Factorizing pressure Poisson...')
    Lp_csc = L_p.tocsc()
    Lp_fac = spla.splu(Lp_csc)

    # Initialize fields
    u = np.ones((Nx + 1, Ny)) * 1.0  # horizontal velocity faces
    v = np.zeros((Nx, Ny + 1))       # vertical velocity faces
    p = np.zeros((Nx, Ny))

    # apply penalization to enforce zero inside cylinder initially
    u = u * (1 - chi_u)  # inside cylinder initial 0
    v = v * (1 - chi_v)

    # History for AB2 convective term
    N_u, N_v = compute_convective(u, v, Nx, Ny, dx, dy)
    N_u_prev = N_u.copy()
    N_v_prev = N_v.copy()

    # Arrays for flattened solves
    u_flat = flatten(u)
    v_flat = flatten(v)

    # Precompute location sets for Dirichlet enforcement on flattened vectors
    inlet_u_flat_indices = inlet_u_indices
    v_zero_flat_indices = list(v_zero_indices)

    # Time loop
    print('Starting time-stepping...')
    progress_interval = max(1, nsteps // 10)
    for step in range(1, nsteps + 1):
        # Extrapolate convective term with AB2
        N_u, N_v = compute_convective(u, v, Nx, Ny, dx, dy)
        N_u_ex = 1.5 * N_u - 0.5 * N_u_prev
        N_v_ex = 1.5 * N_v - 0.5 * N_v_prev

        # RHS for tentative velocity solves: u^n - dt * N_ex + dt*(chi/eta)*u_b
        # u_b = 0 on cylinder
        rhs_u = flatten(u) - dt * flatten(N_u_ex)
        rhs_v = flatten(v) - dt * flatten(N_v_ex)

        # enforce inlet Dirichlet in RHS
        for idx in inlet_u_flat_indices:
            rhs_u[idx] = 1.0

        # enforce v Dirichlet rows in RHS
        for idx in v_zero_flat_indices:
            rhs_v[idx] = 0.0

        # Solve implicit systems (Au * u_star = rhs_u) and (Av * v_star = rhs_v)
        u_star_flat = Au_fac.solve(rhs_u)
        v_star_flat = Av_fac.solve(rhs_v)

        # reshape
        u_star = unflatten_u(u_star_flat, Nx, Ny)
        v_star = unflatten_v(v_star_flat, Nx, Ny)

        # compute divergence of tentative velocity
        div_u_star = divergence(u_star, v_star, Nx, Ny, dx, dy)
        rhs_p = (1.0 / dt) * flatten(div_u_star)

        # enforce pinned pressure node (index pinned_idx)
        rhs_p[pinned_idx] = 0.0
        # solve Poisson L_p * p_corr = rhs_p
        p_corr_flat = Lp_fac.solve(rhs_p)
        p_corr = unflatten_p(p_corr_flat, Nx, Ny)

        # compute pressure gradients at faces
        dpdx = gradient_p_to_u(p_corr, Nx, Ny, dx, dy)
        dpdy = gradient_p_to_v(p_corr, Nx, Ny, dx, dy)

        # update velocities: u^{n+1} = u_star - dt * dpdx ; v^{n+1} = v_star - dt * dpdy
        u_new = u_star - dt * dpdx
        v_new = v_star - dt * dpdy

        # enforce boundary conditions explicitly
        # inlet u
        u_new[0, :] = 1.0
        # free-slip top/bottom: v=0 already via Dirichlet; enforce u derivative zero by copying neighbor
        u_new[:, 0] = u_new[:, 1]
        u_new[:, Ny - 1] = u_new[:, Ny - 2]
        # outlet zero-gradient (copy interior)
        u_new[Nx, :] = u_new[Nx - 1, :]
        v_new[:, 0] = 0.0
        v_new[:, Ny] = 0.0
        # optional: enforce small velocity inside cylinder via penalization effect; implicit penalization already included

        # update fields and history
        u = u_new.copy()
        v = v_new.copy()
        p = p_corr.copy()

        N_u_prev = N_u.copy()
        N_v_prev = N_v.copy()

        # progress printing every 10% of total steps
        if (step % progress_interval) == 0 or step == nsteps:
            perc = int(100 * step / nsteps)
            # compute kinetic energy and max velocity for diagnostics
            uc = 0.5 * (u[0:Nx, :] + u[1:Nx + 1, :]) if Nx >= 1 else u
            vc = 0.5 * (v[:, 0:Ny] + v[:, 1:Ny + 1]) if Ny >= 1 else v
            KE = 0.5 * np.sum(uc**2 + vc**2) * dx * dy
            umax = max(np.abs(u).max(), np.abs(v).max())
            print(f'Step {step}/{nsteps} ({perc}%) t={step*dt:.3f} umax={umax:.3f} KE={KE:.5f}')

    # After time loop, compute velocity magnitude at cell centers for plotting
    u_center = 0.5 * (u[0:Nx, :] + u[1:Nx + 1, :])  # shape (Nx, Ny)
    v_center = 0.5 * (v[:, 0:Ny] + v[:, 1:Ny + 1])  # shape (Nx, Ny)
    vel_mag = np.sqrt(u_center**2 + v_center**2)

    # build meshgrid for plotting (pressure cell centers)
    Xp, Yp = np.meshgrid(x_p, y_p, indexing='xy')

    # plot contour of velocity magnitude and mark cylinder
    fig, ax = plt.subplots(figsize=(10, 4.5))
    cp = ax.contourf(Xp, Yp, vel_mag.T, levels=50, cmap='RdBu_r')
    cbar = fig.colorbar(cp, ax=ax)
    cbar.set_label('Velocity magnitude')
    ax.set_aspect('equal')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_title(f'Velocity magnitude at t={T:.2f}, Re={Re:.0f}')

    # mark cylinder (circle centered at 0,0 radius=0.5)
    theta = np.linspace(0, 2 * np.pi, 200)
    xc = 0.5 * np.cos(theta)
    yc = 0.5 * np.sin(theta)
    ax.plot(xc, yc, 'k', linewidth=1.0)

    fig.tight_layout()
    fig.savefig('velocity_magnitude.png', dpi=200)
    print('Saved figure to velocity_magnitude.png')


if __name__ == '__main__':
    main()
```


#### Script block2:

```python
#!/usr/bin/env python3
"""
2D incompressible Navier-Stokes around a cylinder using MAC grid
with immersed boundary volume penalization. Semi-implicit time stepping.

Key fixes from prior run:
- Use consistent Fortran (column-major) ordering for all flatten/reshape operations
  (ravel(order='F'), reshape(..., order='F')).
- Use np.meshgrid(..., indexing='ij') everywhere so (i,j) ordering matches DOF layout.
- Pin pressure at an interior cell near cylinder center.
- Build penalization diagonals from chi.ravel(order='F').
- Vectorized convective interpolation to avoid slow Python loops.

Saves final contour to 'velocity_magnitude.png'.
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import floor

# Global physical parameter
Re = 100.0
visc = 1.0 / Re   # explicit global visc as required
# Penalization parameter (smaller = stronger enforcement)
eta = 1e-3


def build_grid(x1=3.0, x2=8.0, y1=3.0, y2=3.0, Nx=220, Ny=120):
    Lx = x1 + x2
    Ly = y1 + y2
    dx = Lx / Nx
    dy = Ly / Ny

    # pressure cell centers (Nx, Ny)
    x_p = np.linspace(-x1 + dx / 2, x2 - dx / 2, Nx)
    y_p = np.linspace(-y1 + dy / 2, y2 - dy / 2, Ny)

    # u at vertical faces: (Nx+1, Ny)
    x_u = np.linspace(-x1, x2, Nx + 1)
    y_u = np.linspace(-y1 + dy / 2, y2 - dy / 2, Ny)

    # v at horizontal faces: (Nx, Ny+1)
    x_v = np.linspace(-x1 + dx / 2, x2 - dx / 2, Nx)
    y_v = np.linspace(-y1, y2, Ny + 1)

    return dx, dy, x_u, y_u, x_v, y_v, x_p, y_p


# Indexing helpers consistent with Fortran flatten (order='F')
def idx_u(i, j, Nx, Ny):
    return j * (Nx + 1) + i


def idx_v(i, j, Nx, Ny):
    return j * Nx + i


def idx_p(i, j, Nx, Ny):
    return j * Nx + i


def build_cylinder_masks(x_u, y_u, x_v, y_v, radius=0.5):
    # Use indexing='ij' so Xu.shape == (Nx+1, Ny), Xv.shape == (Nx, Ny+1)
    Xu, Yu = np.meshgrid(x_u, y_u, indexing='ij')
    Xv, Yv = np.meshgrid(x_v, y_v, indexing='ij')

    chi_u = ((Xu ** 2 + Yu ** 2) <= radius ** 2).astype(float)
    chi_v = ((Xv ** 2 + Yv ** 2) <= radius ** 2).astype(float)

    return chi_u, chi_v


def build_laplacian_u(Nx, Ny, dx, dy, inlet_u_indices):
    Nu = (Nx + 1) * Ny
    data = []
    rows = []
    cols = []

    for j in range(Ny):
        for i in range(Nx + 1):
            row = idx_u(i, j, Nx, Ny)
            if row in inlet_u_indices:
                rows.append(row); cols.append(row); data.append(1.0)
                continue
            center = 0.0
            # left
            if i - 1 >= 0:
                left = idx_u(i - 1, j, Nx, Ny)
                rows.append(row); cols.append(left); data.append(1.0 / dx ** 2)
                center += -1.0 / dx ** 2
            # right
            if i + 1 <= Nx:
                right = idx_u(i + 1, j, Nx, Ny)
                rows.append(row); cols.append(right); data.append(1.0 / dx ** 2)
                center += -1.0 / dx ** 2
            # down
            if j - 1 >= 0:
                down = idx_u(i, j - 1, Nx, Ny)
                rows.append(row); cols.append(down); data.append(1.0 / dy ** 2)
                center += -1.0 / dy ** 2
            # up
            if j + 1 <= Ny - 1:
                up = idx_u(i, j + 1, Nx, Ny)
                rows.append(row); cols.append(up); data.append(1.0 / dy ** 2)
                center += -1.0 / dy ** 2
            rows.append(row); cols.append(row); data.append(center)

    L = sp.csr_matrix((data, (rows, cols)), shape=(Nu, Nu))
    return L


def build_laplacian_v(Nx, Ny, dx, dy, v_dirichlet_indices):
    Nv = Nx * (Ny + 1)
    data = []
    rows = []
    cols = []

    for j in range(Ny + 1):
        for i in range(Nx):
            row = idx_v(i, j, Nx, Ny)
            if row in v_dirichlet_indices:
                rows.append(row); cols.append(row); data.append(1.0)
                continue
            center = 0.0
            if i - 1 >= 0:
                left = idx_v(i - 1, j, Nx, Ny)
                rows.append(row); cols.append(left); data.append(1.0 / dx ** 2)
                center += -1.0 / dx ** 2
            if i + 1 <= Nx - 1:
                right = idx_v(i + 1, j, Nx, Ny)
                rows.append(row); cols.append(right); data.append(1.0 / dx ** 2)
                center += -1.0 / dx ** 2
            if j - 1 >= 0:
                down = idx_v(i, j - 1, Nx, Ny)
                rows.append(row); cols.append(down); data.append(1.0 / dy ** 2)
                center += -1.0 / dy ** 2
            if j + 1 <= Ny:
                up = idx_v(i, j + 1, Nx, Ny)
                rows.append(row); cols.append(up); data.append(1.0 / dy ** 2)
                center += -1.0 / dy ** 2
            rows.append(row); cols.append(row); data.append(center)

    L = sp.csr_matrix((data, (rows, cols)), shape=(Nv, Nv))
    return L


def build_laplacian_p(Nx, Ny, dx, dy, pinned_idx=0):
    Np = Nx * Ny
    data = []
    rows = []
    cols = []

    for j in range(Ny):
        for i in range(Nx):
            row = idx_p(i, j, Nx, Ny)
            if row == pinned_idx:
                rows.append(row); cols.append(row); data.append(1.0)
                continue
            center = 0.0
            if i - 1 >= 0:
                left = idx_p(i - 1, j, Nx, Ny)
                rows.append(row); cols.append(left); data.append(1.0 / dx ** 2)
                center += -1.0 / dx ** 2
            if i + 1 <= Nx - 1:
                right = idx_p(i + 1, j, Nx, Ny)
                rows.append(row); cols.append(right); data.append(1.0 / dx ** 2)
                center += -1.0 / dx ** 2
            if j - 1 >= 0:
                down = idx_p(i, j - 1, Nx, Ny)
                rows.append(row); cols.append(down); data.append(1.0 / dy ** 2)
                center += -1.0 / dy ** 2
            if j + 1 <= Ny - 1:
                up = idx_p(i, j + 1, Nx, Ny)
                rows.append(row); cols.append(up); data.append(1.0 / dy ** 2)
                center += -1.0 / dy ** 2
            rows.append(row); cols.append(row); data.append(center)

    L = sp.csr_matrix((data, (rows, cols)), shape=(Np, Np))
    return L


def divergence(u, v, Nx, Ny, dx, dy):
    # u: (Nx+1, Ny), v: (Nx, Ny+1)
    return (u[1:Nx + 1, :] - u[0:Nx, :]) / dx + (v[:, 1:Ny + 1] - v[:, 0:Ny]) / dy


def gradient_p_to_u(p, Nx, Ny, dx, dy):
    dpdx = np.zeros((Nx + 1, Ny))
    # interior faces
    if Nx > 1:
        dpdx[1:Nx, :] = (p[1:Nx, :] - p[0:Nx - 1, :]) / dx
    # boundaries: one-sided approx (Neumann ~ zero gradient assumed)
    dpdx[0, :] = 0.0
    dpdx[Nx, :] = 0.0
    return dpdx


def gradient_p_to_v(p, Nx, Ny, dx, dy):
    dpdy = np.zeros((Nx, Ny + 1))
    if Ny > 1:
        dpdy[:, 1:Ny] = (p[:, 1:Ny] - p[:, 0:Ny - 1]) / dy
    dpdy[:, 0] = 0.0
    dpdy[:, Ny] = 0.0
    return dpdy


def compute_convective(u, v, Nx, Ny, dx, dy):
    # compute du/dx at u faces (shape (Nx+1,Ny))
    du_dx = np.zeros_like(u)
    du_dx[1:Nx, :] = (u[2:Nx + 1, :] - u[0:Nx - 1, :]) / (2 * dx)
    du_dx[0, :] = (u[1, :] - u[0, :]) / dx
    du_dx[Nx, :] = (u[Nx, :] - u[Nx - 1, :]) / dx

    # du/dy at u faces
    du_dy = np.zeros_like(u)
    if Ny >= 3:
        du_dy[:, 1:Ny - 1] = (u[:, 2:Ny] - u[:, 0:Ny - 2]) / (2 * dy)
    if Ny >= 2:
        du_dy[:, 0] = (u[:, 1] - u[:, 0]) / dy
        du_dy[:, Ny - 1] = (u[:, Ny - 1] - u[:, Ny - 2]) / dy

    # interpolate v to u locations: v_at_u shape (Nx+1, Ny)
    v_at_u = np.zeros_like(u)
    # v[:, 0:Ny] corresponds to same j indices as u[:, :]
    if Nx >= 2:
        v_at_u[1:Nx, :] = 0.5 * (v[0:Nx - 1, 0:Ny] + v[1:Nx, 0:Ny])
    if Nx >= 1:
        v_at_u[0, :] = v[0, 0:Ny]
        v_at_u[Nx, :] = v[Nx - 1, 0:Ny]

    N_u = u * du_dx + v_at_u * du_dy

    # dv/dy at v faces
    dv_dy = np.zeros_like(v)
    if Ny >= 2:
        dv_dy[:, 1:Ny] = (v[:, 2:Ny + 1] - v[:, 0:Ny - 1]) / (2 * dy)
    dv_dy[:, 0] = (v[:, 1] - v[:, 0]) / dy
    dv_dy[:, Ny] = (v[:, Ny] - v[:, Ny - 1]) / dy

    # dv/dx at v faces
    dv_dx = np.zeros_like(v)
    if Nx >= 3:
        dv_dx[1:Nx - 1, :] = (v[2:Nx, :] - v[0:Nx - 2, :]) / (2 * dx)
    if Nx >= 2:
        dv_dx[0, :] = (v[1, :] - v[0, :]) / dx
        dv_dx[Nx - 1, :] = (v[Nx - 1, :] - v[Nx - 2, :]) / dx

    # interpolate u to v locations: u_at_v shape (Nx, Ny+1)
    u_at_v = np.zeros_like(v)
    # interior j
    if Ny >= 2:
        # j=1..Ny-1 average of four surrounding u
        u_at_v[:, 1:Ny] = 0.25 * (
            u[0:Nx, 0:Ny - 1] + u[1:Nx + 1, 0:Ny - 1] + u[0:Nx, 1:Ny] + u[1:Nx + 1, 1:Ny]
        )
    # boundaries j=0 and j=Ny
    if Ny >= 1:
        u_at_v[:, 0] = 0.5 * (u[0:Nx, 0] + u[1:Nx + 1, 0])
        u_at_v[:, Ny] = 0.5 * (u[0:Nx, Ny - 1] + u[1:Nx + 1, Ny - 1])

    N_v = u_at_v * dv_dx + v * dv_dy

    return N_u, N_v


# consistent Fortran-order flatten/reshape helpers
def flatten(arr):
    return arr.ravel(order='F')


def unflatten_u(vec, Nx, Ny):
    return vec.reshape((Nx + 1, Ny), order='F')


def unflatten_v(vec, Nx, Ny):
    return vec.reshape((Nx, Ny + 1), order='F')


def unflatten_p(vec, Nx, Ny):
    return vec.reshape((Nx, Ny), order='F')


def apply_dirichlet_u(vec_u, u_inlet_value, inlet_indices):
    for idx in inlet_indices:
        vec_u[idx] = u_inlet_value
    return vec_u


def apply_dirichlet_v(vec_v, v_zero_indices):
    for idx in v_zero_indices:
        vec_v[idx] = 0.0
    return vec_v


def main():
    # Domain and grid
    x1 = 3.0
    x2 = 8.0
    y1 = 3.0
    y2 = 3.0
    # moderate grid (reduce if debugging)
    Nx = 220
    Ny = 120

    dx, dy, x_u, y_u, x_v, y_v, x_p, y_p = build_grid(x1, x2, y1, y2, Nx=Nx, Ny=Ny)

    # masks for cylinder penalization (using consistent indexing)
    chi_u, chi_v = build_cylinder_masks(x_u, y_u, x_v, y_v, radius=0.5)

    Nu = (Nx + 1) * Ny
    Nv = Nx * (Ny + 1)
    Np = Nx * Ny

    # inlet u-face indices (left boundary) flattened using idx_u
    inlet_u_indices = [idx_u(0, j, Nx, Ny) for j in range(Ny)]
    inlet_u_set = set(inlet_u_indices)

    # v Dirichlet indices (top/bottom) flattened
    v_zero_indices = set()
    for i in range(Nx):
        v_zero_indices.add(idx_v(i, 0, Nx, Ny))
        v_zero_indices.add(idx_v(i, Ny, Nx, Ny))

    # Build Laplacians
    print('Building operators...')
    L_u = build_laplacian_u(Nx, Ny, dx, dy, inlet_u_set)
    L_v = build_laplacian_v(Nx, Ny, dx, dy, v_zero_indices)

    # pin pressure at interior cell near center
    pinned_idx = idx_p(Nx // 2, Ny // 2, Nx, Ny)
    L_p = build_laplacian_p(Nx, Ny, dx, dy, pinned_idx=pinned_idx)

    # Time-stepping parameters
    T = 20.0
    Umax = 1.0
    CFL = 0.45
    dt = CFL * min(dx, dy) / Umax
    nsteps = int(floor(T / dt))
    dt = T / nsteps

    print(f'Nx={Nx} Ny={Ny} dx={dx:.4f} dy={dy:.4f} dt={dt:.5f} steps={nsteps}')

    # Build constant matrices for velocity solves
    Iu = sp.eye(Nu, format='csr')
    Iv = sp.eye(Nv, format='csr')

    chi_u_flat = chi_u.ravel(order='F')
    chi_v_flat = chi_v.ravel(order='F')

    Pen_u = sp.diags(dt * (chi_u_flat / eta), format='csr')
    Pen_v = sp.diags(dt * (chi_v_flat / eta), format='csr')

    Au = Iu - dt * visc * L_u + Pen_u
    Av = Iv - dt * visc * L_v + Pen_v

    # Enforce Dirichlet rows
    Au = Au.tolil()
    for idx in inlet_u_indices:
        Au.rows[idx] = [idx]
        Au.data[idx] = [1.0]
    Au = Au.tocsc()

    Av = Av.tolil()
    for idx in v_zero_indices:
        Av.rows[idx] = [idx]
        Av.data[idx] = [1.0]
    Av = Av.tocsc()

    # Factorize matrices once
    print('Factorizing velocity matrices...')
    Au_fac = spla.splu(Au)
    Av_fac = spla.splu(Av)

    print('Factorizing pressure Poisson...')
    Lp_csc = L_p.tocsc()
    Lp_fac = spla.splu(Lp_csc)

    # Initialize fields with consistent shapes
    u = np.ones((Nx + 1, Ny)) * 1.0
    v = np.zeros((Nx, Ny + 1))
    p = np.zeros((Nx, Ny))

    # zero inside cylinder initially
    u = u * (1 - chi_u)
    v = v * (1 - chi_v)

    # initial convective terms
    N_u, N_v = compute_convective(u, v, Nx, Ny, dx, dy)
    N_u_prev = N_u.copy()
    N_v_prev = N_v.copy()

    # flattened arrays
    u_flat = flatten(u)
    v_flat = flatten(v)

    inlet_u_flat_indices = inlet_u_indices
    v_zero_flat_indices = list(v_zero_indices)

    print('Starting time-stepping...')
    progress_interval = max(1, nsteps // 10)
    for step in range(1, nsteps + 1):
        N_u, N_v = compute_convective(u, v, Nx, Ny, dx, dy)
        N_u_ex = 1.5 * N_u - 0.5 * N_u_prev
        N_v_ex = 1.5 * N_v - 0.5 * N_v_prev

        rhs_u = flatten(u) - dt * flatten(N_u_ex)
        rhs_v = flatten(v) - dt * flatten(N_v_ex)

        # enforce Dirichlet in RHS
        for idx in inlet_u_flat_indices:
            rhs_u[idx] = 1.0
        for idx in v_zero_flat_indices:
            rhs_v[idx] = 0.0

        # Solve implicit systems
        u_star_flat = Au_fac.solve(rhs_u)
        v_star_flat = Av_fac.solve(rhs_v)

        u_star = unflatten_u(u_star_flat, Nx, Ny)
        v_star = unflatten_v(v_star_flat, Nx, Ny)

        # pressure Poisson RHS
        div_u_star = divergence(u_star, v_star, Nx, Ny, dx, dy)
        rhs_p = (1.0 / dt) * flatten(div_u_star)
        rhs_p[pinned_idx] = 0.0

        p_corr_flat = Lp_fac.solve(rhs_p)
        p_corr = unflatten_p(p_corr_flat, Nx, Ny)

        dpdx = gradient_p_to_u(p_corr, Nx, Ny, dx, dy)
        dpdy = gradient_p_to_v(p_corr, Nx, Ny, dx, dy)

        u_new = u_star - dt * dpdx
        v_new = v_star - dt * dpdy

        # enforce boundary conditions explicitly
        u_new[0, :] = 1.0
        # free-slip: du/dy = 0 at top/bottom -> copy neighbor
        if Ny >= 2:
            u_new[:, 0] = u_new[:, 1]
            u_new[:, Ny - 1] = u_new[:, Ny - 2]
        # outlet zero-gradient
        u_new[Nx, :] = u_new[Nx - 1, :]
        v_new[:, 0] = 0.0
        v_new[:, Ny] = 0.0

        # update
        u = u_new.copy()
        v = v_new.copy()
        p = p_corr.copy()

        N_u_prev = N_u.copy()
        N_v_prev = N_v.copy()

        if (step % progress_interval) == 0 or step == nsteps:
            perc = int(100 * step / nsteps)
            uc = 0.5 * (u[0:Nx, :] + u[1:Nx + 1, :]) if Nx >= 1 else u
            vc = 0.5 * (v[:, 0:Ny] + v[:, 1:Ny + 1]) if Ny >= 1 else v
            KE = 0.5 * np.sum(uc ** 2 + vc ** 2) * dx * dy
            umax = max(np.abs(u).max(), np.abs(v).max())
            print(f'Step {step}/{nsteps} ({perc}%) t={step*dt:.3f} umax={umax:.3f} KE={KE:.5f}')

    # compute velocity magnitude at pressure cell centers
    u_center = 0.5 * (u[0:Nx, :] + u[1:Nx + 1, :])
    v_center = 0.5 * (v[:, 0:Ny] + v[:, 1:Ny + 1])
    vel_mag = np.sqrt(u_center ** 2 + v_center ** 2)

    # plotting meshgrid consistent with indexing='ij'
    Xp, Yp = np.meshgrid(x_p, y_p, indexing='ij')

    fig, ax = plt.subplots(figsize=(10, 4.5))
    cp = ax.contourf(Xp, Yp, vel_mag, levels=50, cmap='RdBu_r')
    cbar = fig.colorbar(cp, ax=ax)
    cbar.set_label('Velocity magnitude')
    ax.set_aspect('equal')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    ax.set_title(f'Velocity magnitude at t={T:.2f}, Re={Re:.0f}')

    theta = np.linspace(0, 2 * np.pi, 200)
    xc = 0.5 * np.cos(theta)
    yc = 0.5 * np.sin(theta)
    ax.plot(xc, yc, 'k', linewidth=1.0)

    fig.tight_layout()
    fig.savefig('velocity_magnitude.png', dpi=200)
    print('Saved figure to velocity_magnitude.png')


if __name__ == '__main__':
    main()
```




### runtime_outputs

#### Output block1

Building operators...
Nx=220 Ny=120 dx=0.0500 dy=0.0500 dt=0.02252 steps=888
Factorizing velocity matrices...
Factorizing pressure Poisson...

Runtime Error:
Traceback (most recent call last):
  File "C:\Users\admin\AppData\Local\Temp\ipykernel_49784\599011946.py", line 140, in execute_code_tool
    exec(code, global_scope, global_scope)
  File "<string>", line 616, in <module>
  File "<string>", line 493, in main
ValueError: operands could not be broadcast together with shapes (221,120) (120,221) 



#### Output block2

Building operators...
Nx=220 Ny=120 dx=0.0500 dy=0.0500 dt=0.02252 steps=888
Factorizing velocity matrices...
Factorizing pressure Poisson...
Starting time-stepping...
Step 88/888 (9%) t=1.982 umax=1.496 KE=34.53532
Step 176/888 (19%) t=3.964 umax=1.459 KE=35.30346
Step 264/888 (29%) t=5.946 umax=1.435 KE=35.95102
Step 352/888 (39%) t=7.928 umax=1.479 KE=36.52678
Step 440/888 (49%) t=9.910 umax=1.485 KE=37.11164
Step 528/888 (59%) t=11.892 umax=1.491 KE=37.77028
Step 616/888 (69%) t=13.874 umax=1.538 KE=38.24978
Step 704/888 (79%) t=15.856 umax=1.697 KE=39.21877
Step 792/888 (89%) t=17.838 umax=1.653 KE=39.16903
Step 880/888 (99%) t=19.820 umax=1.488 KE=37.91005
Step 888/888 (100%) t=20.000 umax=1.489 KE=38.12493
Saved figure to velocity_magnitude.png




### review_decision
accept

### review_comments
	Current Stage [A/2]
1. MAJOR

Category: runtime
Issue: Runtime broadcasting error: operands could not be broadcast together with shapes (221,120) (120,221).

Recommendation: This is caused by inconsistent memory/layout/orientation between flattened and reshaped arrays (flatten()/unflatten_*), combined with inconsistent meshgrid indexing. Fix by making the flatten/reshape operations use a single, consistent storage order and removing the ad-hoc transposes. Concretely:
- Use column-major (Fortran) order everywhere for flatten/reshape: flatten(arr) -> arr.ravel(order='F'); unflatten_u(vec) -> vec.reshape((Nx+1, Ny), order='F'), unflatten_v -> vec.reshape((Nx, Ny+1), order='F'), unflatten_p -> vec.reshape((Nx, Ny), order='F'). Remove the .T in those functions.
- Make all meshgrid calls use consistent indexing (prefer indexing='ij' for MAC grids) and generate masks with the same ordering.
- When building diagonal penalization and any flattened index lists, always use the same order='F' to match the flattened DOF ordering.
Applying these changes will resolve the (221,120) vs (120,221) transpose mismatch and stop the ValueError.


2. MINOR

Category: structure
Issue: Inconsistent use of meshgrid indexing and array layout (indexing='xy' used in some places, inconsistent ravel/reshape ordering).

Recommendation: Use a consistent convention across the code: for a MAC (i,j) ordering use np.meshgrid(x, y, indexing='ij') everywhere. Use order='F' for ravel/reshape so that the flattened index ordering matches the (i,j) ordering used to assemble matrices and masks. Update build_cylinder_masks, plotting meshgrid creation, and any place using meshgrid to indexing='ij'. This will harmonize coordinate arrays, matrix assembly and plotting (then you won't need to transpose vel_mag for plotting).


3. MINOR

Category: runtime
Issue: Pressure pinning index pinned_idx = 0 is at a corner cell which may produce non-physical reference and interacts badly with Neumann approximations.

Recommendation: Pin pressure at an interior cell (preferably nearest the cylinder center). E.g. pinned_idx = idx_p(Nx//2, Ny//2, Nx, Ny). This removes the nullspace while avoiding biasing boundary values. Ensure the same flattened ordering (order='F') is used when referencing pinned_idx.


4. MAJOR

Category: structure
Issue: Flatten/unflatten helpers and the way chi arrays are flattened for penalization are inconsistent and fragile.

Recommendation: Replace all occurrences of .ravel()/.reshape(...) and ad-hoc .T with explicit order='F' usage. Build Pen_u and Pen_v from chi_u.ravel(order='F') and chi_v.ravel(order='F'). Make sure inlet/v Dirichlet index lists are generated in the same order as flat ordering. Example replacements:
- def flatten(a): return a.ravel(order='F')
- def unflatten_u(vec,...): return vec.reshape((Nx+1, Ny), order='F')
- same pattern for v and p
This will make the penalization diagonal alignment and LU solves correct and reproducible.


5. MINOR

Category: accuracy
Issue: Numerical and performance concerns: compute_convective uses nested Python loops (slow) and boundary/Neumann treatments are ad-hoc and incomplete in a few Laplacian builders.

Recommendation: For correctness first fix the shape/layout bugs. For performance/accuracy (optional but recommended):
- Vectorize interpolation (avoid Python loops) to accelerate convective term (use array slicing and averaging). 
- Consider assembling Laplacians via Kronecker products (sp.kron) or using stencil-based sp.diags patterns for clarity and speed. 
- Improve Neumann/outlet handling (use ghost-value approximations consistent with projection) and consider pinning pressure in center as recommended. 
- Consider reducing Nx, Ny while debugging to speed up iteration (e.g. Nx=80, Ny=40) and only increase once correctness is confirmed.
These changes are not strictly blocking but will significantly improve runtime and robustness.



	Current Stage [B/2]
1. MINOR

Category: structure
Issue: Overall implementation and runtime behaviour

Recommendation: Acceptable as-is for a working research / prototyping solver. The implementation follows the requested MAC layout, uses Fortran ordering consistently, treats the penalization implicitly, pins pressure to remove the nullspace, prints progress only every 10% of steps, and saves the final figure without plt.show(). Numerical output in the runtime log shows a stable run to t=20.0. No runtime failures were observed in the provided run.


2. MINOR

Category: structure
Issue: Inconsistent enforcement of free-slip boundary conditions (operator vs post-processing)

Recommendation: Make the boundary-condition implementation consistent between linear-operator assembly and the time-stepping updates. Right now the free-slip top/bottom conditions for u are enforced explicitly after the projection step (u_new[:,0] = u_new[:,1], u_new[:,Ny-1] = u_new[:,Ny-2]) but the u-Laplacian matrix is built without embedding these Neumann-style BCs. Either incorporate the free-slip conditions into L_u / Au (modify rows/columns appropriately) or explicitly enforce them in the linear solves (set those rows as Dirichlet-like when appropriate). Doing so will avoid a small inconsistency between the implicit solve and the post-projection correction and improve accuracy / convergence.


3. MINOR

Category: accuracy
Issue: Pressure-pinning location is inside the cylinder (penalized region)

Recommendation: Pin the pressure at a point known to be in the fluid region rather than at a cell that may lie inside the penalized (solid) region. Currently pinned_idx = idx_p(Nx//2, Ny//2) is at the geometric center, which for this problem is inside the cylinder. Pinning inside the solid is not wrong per se, but it can hide pressure-conditioning issues or make the physical interpretation of the pressure field awkward. Consider selecting a pinned pressure index on the downstream boundary or at a guaranteed fluid cell.


4. MINOR

Category: accuracy
Issue: Penalization RHS for non-zero body velocity is not included

Recommendation: If the immersed body had a nonzero prescribed velocity u_b (moving body), the penalization RHS term dt*(chi/eta)*u_b must be included in rhs_u/rhs_v. Currently the code relies on u_b=0 (stationary no-slip), so omission is correct here; add a short comment and include the term in the RHS assembly if you ever switch to moving boundaries.


5. MINOR

Category: structure
Issue: Scalability: direct LU factorization of large matrices may become memory-bound

Recommendation: For larger grids or production runs consider replacing direct factorization of the three large sparse matrices with iterative solvers (GMRES/CG) and good preconditioners/AMG for the pressure Poisson. The current direct factorization works and is robust for the presented grid, but memory and CPU cost will grow quickly for finer meshes. Consider pyamg / AMG preconditioners for L_p and (block) preconditioning for the velocity solves.





### iteration_count
2

### rev_image_description
	Current Stage [A/2]


	Current Stage [B/2]




