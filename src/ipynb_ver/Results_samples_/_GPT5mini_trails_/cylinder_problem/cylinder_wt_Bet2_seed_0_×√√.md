### 🛌🏻 Re4gent with Multi-Modal Review

#### Configs 🏷️:
MaxCount: 2, MultiModalReview: False, AlternativeNum: 2, MaxComments: 5, ConsModule_ON: False, CUT_OUTEXT: 2000

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
	Current Stage [A/3]

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


	Current Stage [B/3]
This solver uses a staggered MAC grid (u on vertical faces, v on horizontal faces, p at cell centers) and an immersed-boundary volume-penalization to enforce the cylinder. Key choices made to satisfy the review and problem requirements:

- Consistent memory/flattening convention: meshgrid(..., indexing='ij') and Fortran-order flatten/reshape (order='F') everywhere so sparse matrix row/column indexing (which uses idx_* helpers) aligns with vector layouts used in linear solves (splu).
- Operators (velocity systems Au, Av and pressure Laplacian Lp) are assembled once, with the penalization term included implicitly as a diagonal in Au/Av. Dirichlet rows are overwritten once before factorization. Matrices are factorized once (scipy.splu) and reused each time step.
- Semi-implicit time stepping: AB2 extrapolation for convection, implicit viscous + penalization solves for tentative velocities, pressure Poisson projection, then velocity projection. Penalization is treated implicitly via diag additions to Au/Av.
- Vectorized convective term calculation (no Python i/j loops) for performance.
- Pressure pinning chosen at a physically neutral outflow-center cell and an additional zero-mean projection is applied to avoid nullspace bias. Divergence residual (L2) is computed and printed with progress diagnostics.
- Progress messages printed only each 10% of the total steps include umax, KE and divergence L2.
- Final velocity magnitude at cell centers is plotted using 'RdBu_r', cylinder boundary overplotted, and figure saved as 'velocity_magnitude.png' (no plt.show()).

Notes on consistency and correctness:
- All flattening/reshaping uses order='F' and idx_* functions are consistent with that layout.
- build_cylinder_masks and all meshgrid calls use indexing='ij' so arrays have shapes matching the MAC arrays (axis 0 -> x index, axis 1 -> y index).
- No use of removed numpy APIs or forbidden SciPy kwargs (e.g., no 'tol', use 'atol' if needed).

	Current Stage [C/3]
This solver implements a MAC-grid incompressible Navier–Stokes solver with immersed-boundary volume-penalization. Key algorithmic choices:
- Discrete divergence operator D is assembled once (sparse) mapping face-centered velocities to cell-centered divergence. The gradient operator G is its transpose (G = D.T) so L_p = D * G (positive semi-definite). This enforces algebraic consistency: projection uses G and L_p so D * (u_star - dt*G*p) = 0 when L_p p = (1/dt) D u_star.
- Momentum matrices Au and Av (implicit viscous + implicit penalization on diagonal) are assembled once and factorized for repeated solves (splu).
- Semi-implicit AB2 time stepping for convection, implicit viscous+penalization for velocity solve, projection using D/G consistent operators, and re-application of penalization (no-slip in cylinder) after projection to prevent velocities inside the solid.
- Pressure singularity removed by pinning a single pressure cell; no mean-subtraction is applied (pinned entry is explicitly enforced after solve as a sanity). Diagnostics include divergence L2 and L_inf along with kinetic energy and umax printed every 10% of steps.
- All arrays flattened/reshaped using Fortran-order and meshgrid(..., indexing='ij') to keep consistent layouts.




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

#!/usr/bin/env python3
"""
2D incompressible Navier-Stokes around a cylinder on a staggered MAC grid
with immersed-boundary volume penalization. The code follows the review
recommendations: consistent indexing (indexing='ij'), Fortran-order flattening
(order='F'), vectorized convective calculations and single factorization
of linear systems.
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
visc = 1.0 / Re   # use visc everywhere per instructions

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

    # pressure cell centers (Nx x Ny)
    x_p = np.linspace(-x1 + dx/2, x2 - dx/2, Nx)
    y_p = np.linspace(-y1 + dy/2, y2 - dy/2, Ny)

    # u at vertical faces: x faces i=0..Nx (Nx+1), y centers j=0..Ny-1 (Ny)
    x_u = np.linspace(-x1, x2, Nx+1)
    y_u = np.linspace(-y1 + dy/2, y2 - dy/2, Ny)

    # v at horizontal faces: x centers i=0..Nx-1 (Nx), y faces j=0..Ny (Ny+1)
    x_v = np.linspace(-x1 + dx/2, x2 - dx/2, Nx)
    y_v = np.linspace(-y1, y2, Ny+1)

    return dx, dy, x_u, y_u, x_v, y_v, x_p, y_p


# index helpers consistent with Fortran-order flattening (order='F')
def idx_u(i, j, Nx, Ny):
    return i + j * (Nx + 1)


def idx_v(i, j, Nx, Ny):
    return i + j * Nx


def idx_p(i, j, Nx, Ny):
    return i + j * Nx


def build_cylinder_masks(x_u, y_u, x_v, y_v, radius=0.5):
    """Build boolean masks chi_u and chi_v for penalization (1 inside cylinder).
    Use indexing='ij' so output shapes match unstaggered array shapes:
      chi_u: (Nx+1, Ny), chi_v: (Nx, Ny+1)
    """
    Xu, Yu = np.meshgrid(x_u, y_u, indexing='ij')
    Xv, Yv = np.meshgrid(x_v, y_v, indexing='ij')

    chi_u = ((Xu**2 + Yu**2) <= radius**2).astype(float)
    chi_v = ((Xv**2 + Yv**2) <= radius**2).astype(float)

    return chi_u, chi_v


def build_laplacian_u(Nx, Ny, dx, dy, inlet_u_indices):
    """Build sparse Laplacian (second-order) for u-grid (Nu = (Nx+1)*Ny).
    Dirichlet inlet rows are left as identity elsewhere we assemble a standard 5-point
    stencil with mirror/one-sided treatment for outer Neumann-like boundaries.
    """
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
            # left neighbor
            if i - 1 >= 0:
                left = idx_u(i - 1, j, Nx, Ny)
                rows.append(row); cols.append(left); data.append(1.0 / dx**2)
                center += -1.0 / dx**2
            # right neighbor
            if i + 1 <= Nx:
                right = idx_u(i + 1, j, Nx, Ny)
                rows.append(row); cols.append(right); data.append(1.0 / dx**2)
                center += -1.0 / dx**2

            # down (j-1)
            if j - 1 >= 0:
                down = idx_u(i, j - 1, Nx, Ny)
                rows.append(row); cols.append(down); data.append(1.0 / dy**2)
                center += -1.0 / dy**2
            # up (j+1)
            if j + 1 <= Ny - 1:
                up = idx_u(i, j + 1, Nx, Ny)
                rows.append(row); cols.append(up); data.append(1.0 / dy**2)
                center += -1.0 / dy**2

            rows.append(row); cols.append(row); data.append(center)

    L = sp.csr_matrix((data, (rows, cols)), shape=(Nu, Nu))
    return L


def build_laplacian_v(Nx, Ny, dx, dy, v_dirichlet_indices):
    """Build sparse Laplacian for v-grid (Nv = Nx*(Ny+1)).
    Dirichlet rows are set to identity; else assemble 5-point stencil.
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
    Pin one node to remove nullspace; choose pinned index at outflow center in main.
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
            # right
            if i + 1 <= Nx - 1:
                right = idx_p(i + 1, j, Nx, Ny)
                rows.append(row); cols.append(right); data.append(1.0 / dx**2)
                center += -1.0 / dx**2
            # down
            if j - 1 >= 0:
                down = idx_p(i, j - 1, Nx, Ny)
                rows.append(row); cols.append(down); data.append(1.0 / dy**2)
                center += -1.0 / dy**2
            # up
            if j + 1 <= Ny - 1:
                up = idx_p(i, j + 1, Nx, Ny)
                rows.append(row); cols.append(up); data.append(1.0 / dy**2)
                center += -1.0 / dy**2

            rows.append(row); cols.append(row); data.append(center)

    L = sp.csr_matrix((data, (rows, cols)), shape=(Np, Np))
    return L


def divergence(u, v, Nx, Ny, dx, dy):
    """Compute discrete divergence at pressure cell centers from face-centered u and v.
    u shape: (Nx+1, Ny), v shape: (Nx, Ny+1)
    returns div shape (Nx, Ny)
    """
    div = (u[1:Nx+1, :] - u[0:Nx, :]) / dx + (v[:, 1:Ny+1] - v[:, 0:Ny]) / dy
    return div


def gradient_p_to_u(p, Nx, Ny, dx, dy):
    """Compute dp/dx at u face locations from cell-centered pressure p (Nx x Ny).
    Returns array of shape (Nx+1, Ny)
    """
    dpdx = np.zeros((Nx + 1, Ny))
    # interior faces
    if Nx > 1:
        dpdx[1:Nx, :] = (p[1:Nx, :] - p[0:Nx-1, :]) / dx
    # boundaries: one-sided -> assume zero-gradient approximation
    dpdx[0, :] = 0.0
    dpdx[Nx, :] = 0.0
    return dpdx


def gradient_p_to_v(p, Nx, Ny, dx, dy):
    """Compute dp/dy at v face locations from p.
    Returns array shape (Nx, Ny+1)
    """
    dpdy = np.zeros((Nx, Ny + 1))
    if Ny > 1:
        dpdy[:, 1:Ny] = (p[:, 1:Ny] - p[:, 0:Ny-1]) / dy
    dpdy[:, 0] = 0.0
    dpdy[:, Ny] = 0.0
    return dpdy


def compute_convective(u, v, Nx, Ny, dx, dy):
    """Vectorized convective term computation on staggered grid.
    Returns N_u (Nx+1, Ny) and N_v (Nx, Ny+1).
    """
    # du/dx
    du_dx = np.zeros_like(u)
    if Nx >= 2:
        du_dx[1:Nx, :] = (u[2:Nx+1, :] - u[0:Nx-1, :]) / (2 * dx)
    du_dx[0, :] = (u[1, :] - u[0, :]) / dx
    du_dx[Nx, :] = (u[Nx, :] - u[Nx - 1, :]) / dx

    # du/dy
    du_dy = np.zeros_like(u)
    if Ny >= 3:
        du_dy[:, 1:Ny-1] = (u[:, 2:Ny] - u[:, 0:Ny-2]) / (2 * dy)
    if Ny >= 2:
        du_dy[:, 0] = (u[:, 1] - u[:, 0]) / dy
        du_dy[:, Ny-1] = (u[:, Ny-1] - u[:, Ny-2]) / dy

    # interpolate v to u locations (vectorized)
    v_at_u = np.zeros_like(u)
    if Nx >= 2:
        v_at_u[1:Nx, :] = 0.5 * (v[0:Nx-1, 0:Ny] + v[1:Nx, 0:Ny])
    # boundaries
    v_at_u[0, :] = v[0, 0:Ny]
    v_at_u[Nx, :] = v[Nx-1, 0:Ny]

    N_u = u * du_dx + v_at_u * du_dy

    # dv/dy
    dv_dy = np.zeros_like(v)
    if Ny >= 2:
        dv_dy[:, 1:Ny] = (v[:, 2:Ny+1] - v[:, 0:Ny-1]) / (2 * dy)
        dv_dy[:, 0] = (v[:, 1] - v[:, 0]) / dy
        dv_dy[:, Ny] = (v[:, Ny] - v[:, Ny-1]) / dy

    # dv/dx
    dv_dx = np.zeros_like(v)
    if Nx >= 3:
        dv_dx[1:Nx-1, :] = (v[2:Nx, :] - v[0:Nx-2, :]) / (2 * dx)
    if Nx >= 2:
        dv_dx[0, :] = (v[1, :] - v[0, :]) / dx
        dv_dx[Nx-1, :] = (v[Nx-1, :] - v[Nx-2, :]) / dx

    # interpolate u to v locations (vectorized)
    u_at_v = np.zeros_like(v)
    # interior v-j faces (1..Ny-1) use average of four surrounding u-values
    if Ny >= 2:
        # u left-bottom: u[0:Nx, 0:Ny-1]
        # u right-bottom: u[1:Nx+1, 0:Ny-1]
        # u left-top: u[0:Nx, 1:Ny]
        # u right-top: u[1:Nx+1, 1:Ny]
        u_at_v[:, 1:Ny] = 0.25 * (
            u[0:Nx, 0:Ny-1] + u[1:Nx+1, 0:Ny-1] + u[0:Nx, 1:Ny] + u[1:Nx+1, 1:Ny]
        )
    # bottom j=0: average bottom cell-centered u's
    u_at_v[:, 0] = 0.5 * (u[0:Nx, 0] + u[1:Nx+1, 0])
    # top j=Ny: average top cell-centered u's
    u_at_v[:, Ny] = 0.5 * (u[0:Nx, Ny-1] + u[1:Nx+1, Ny-1])

    N_v = u_at_v * dv_dx + v * dv_dy

    return N_u, N_v


# consistent flatten/unflatten using Fortran order
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

    # v Dirichlet indices (top/bottom)
    v_zero_indices = set()
    for i in range(Nx):
        v_zero_indices.add(idx_v(i, 0, Nx, Ny))   # bottom v=0
        v_zero_indices.add(idx_v(i, Ny, Nx, Ny))  # top v=0

    # Build Laplacians
    print('Building operators...')
    L_u = build_laplacian_u(Nx, Ny, dx, dy, inlet_u_set)
    L_v = build_laplacian_v(Nx, Ny, dx, dy, v_zero_indices)
    # choose a pinned pressure index at outflow center to minimize bias
    pinned_idx = idx_p(Nx - 1, Ny // 2, Nx, Ny)
    L_p = build_laplacian_p(Nx, Ny, dx, dy, pinned_idx=pinned_idx)

    # Time-stepping parameters
    T = 20.0
    Umax = 1.0
    CFL = 0.45
    dt = CFL * min(dx, dy) / Umax
    nsteps = int(floor(T / dt))
    dt = T / nsteps

    print(f'Nx={Nx} Ny={Ny} dx={dx:.4f} dy={dy:.4f} dt={dt:.5f} steps={nsteps}')

    # Build constant matrices for velocity solves: Au = I - dt*visc*L_u + dt*(chi/eta) diag
    Iu = sp.eye(Nu, format='csr')
    Iv = sp.eye(Nv, format='csr')

    # diag penalization arrays (Fortran order)
    chi_u_flat = flatten(chi_u)
    chi_v_flat = flatten(chi_v)

    Pen_u = sp.diags(dt * (chi_u_flat / eta), format='csr')
    Pen_v = sp.diags(dt * (chi_v_flat / eta), format='csr')

    Au = Iu - dt * visc * L_u + Pen_u
    Av = Iv - dt * visc * L_v + Pen_v

    # Enforce Dirichlet rows in Au and Av to identity
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

    # Factorize matrices for fast repeated solves
    print('Factorizing velocity matrices...')
    Au_fac = spla.splu(Au)
    Av_fac = spla.splu(Av)

    # Factorize pressure Laplacian (pinning included) once
    print('Factorizing pressure Poisson...')
    Lp_csc = L_p.tocsc()
    Lp_fac = spla.splu(Lp_csc)

    # Initialize fields (Fortran-order shapes)
    u = np.ones((Nx + 1, Ny), order='F') * 1.0
    v = np.zeros((Nx, Ny + 1), order='F')
    p = np.zeros((Nx, Ny), order='F')

    # apply penalization to enforce zero inside cylinder initially
    u = u * (1.0 - chi_u)
    v = v * (1.0 - chi_v)

    # History for AB2 convective term
    N_u, N_v = compute_convective(u, v, Nx, Ny, dx, dy)
    N_u_prev = N_u.copy()
    N_v_prev = N_v.copy()

    # Precompute index lists for enforcing Dirichlet in flattened vectors
    inlet_u_flat_indices = inlet_u_indices
    v_zero_flat_indices = list(v_zero_indices)

    # Time loop
    print('Starting time-stepping...')
    progress_interval = max(1, nsteps // 10)
    for step in range(1, nsteps + 1):
        # convective term and AB2 extrapolation
        N_u, N_v = compute_convective(u, v, Nx, Ny, dx, dy)
        N_u_ex = 1.5 * N_u - 0.5 * N_u_prev
        N_v_ex = 1.5 * N_v - 0.5 * N_v_prev

        # RHS for tentative velocity solves
        rhs_u = flatten(u) - dt * flatten(N_u_ex)
        rhs_v = flatten(v) - dt * flatten(N_v_ex)

        # enforce inlet Dirichlet in RHS
        for idx in inlet_u_flat_indices:
            rhs_u[idx] = 1.0
        # enforce v Dirichlet rows in RHS
        for idx in v_zero_flat_indices:
            rhs_v[idx] = 0.0

        # Solve implicit systems
        u_star_flat = Au_fac.solve(rhs_u)
        v_star_flat = Av_fac.solve(rhs_v)

        u_star = unflatten_u(u_star_flat, Nx, Ny)
        v_star = unflatten_v(v_star_flat, Nx, Ny)

        # compute divergence of tentative velocity
        div_u_star = divergence(u_star, v_star, Nx, Ny, dx, dy)
        rhs_p = (1.0 / dt) * flatten(div_u_star)

        # enforce pinned pressure RHS (matching pinned row in L_p)
        rhs_p[pinned_idx] = 0.0

        # solve Poisson L_p * p_corr = rhs_p
        p_corr_flat = Lp_fac.solve(rhs_p)
        # remove any residual mean to avoid location bias
        p_corr_flat = p_corr_flat - np.mean(p_corr_flat)
        p_corr = unflatten_p(p_corr_flat, Nx, Ny)

        # compute pressure gradients at faces
        dpdx = gradient_p_to_u(p_corr, Nx, Ny, dx, dy)
        dpdy = gradient_p_to_v(p_corr, Nx, Ny, dx, dy)

        # project velocities
        u_new = u_star - dt * dpdx
        v_new = v_star - dt * dpdy

        # enforce boundary conditions explicitly
        u_new[0, :] = 1.0
        # free-slip (du/dy=0) at top/bottom -> copy neighbor
        u_new[:, 0] = u_new[:, 1]
        u_new[:, Ny - 1] = u_new[:, Ny - 2]
        # outlet zero-gradient
        u_new[Nx, :] = u_new[Nx - 1, :]
        # v top/bottom Dirichlet
        v_new[:, 0] = 0.0
        v_new[:, Ny] = 0.0

        # update fields and history
        u = u_new.copy()
        v = v_new.copy()
        p = p_corr.copy()

        N_u_prev = N_u.copy()
        N_v_prev = N_v.copy()

        # compute divergence residual after projection for diagnostics
        div_res = divergence(u, v, Nx, Ny, dx, dy)
        div_L2 = np.sqrt(np.sum(div_res**2) * dx * dy)

        # progress printing every 10% of total steps
        if (step % progress_interval) == 0 or step == nsteps:
            perc = int(100 * step / nsteps)
            uc = 0.5 * (u[0:Nx, :] + u[1:Nx+1, :])
            vc = 0.5 * (v[:, 0:Ny] + v[:, 1:Ny+1])
            KE = 0.5 * np.sum(uc**2 + vc**2) * dx * dy
            umax = max(np.abs(u).max(), np.abs(v).max())
            print(f'Step {step}/{nsteps} ({perc}%) t={step*dt:.3f} umax={umax:.3f} KE={KE:.5f} div_L2={div_L2:.5e}')

    # After time loop, compute velocity magnitude at cell centers for plotting
    u_center = 0.5 * (u[0:Nx, :] + u[1:Nx + 1, :])  # shape (Nx, Ny)
    v_center = 0.5 * (v[:, 0:Ny] + v[:, 1:Ny + 1])  # shape (Nx, Ny)
    vel_mag = np.sqrt(u_center**2 + v_center**2)

    # build meshgrid for plotting (pressure cell centers) with indexing='ij'
    Xp, Yp = np.meshgrid(x_p, y_p, indexing='ij')

    fig, ax = plt.subplots(figsize=(10, 4.5))
    cp = ax.contourf(Xp, Yp, vel_mag, levels=50, cmap='RdBu_r')
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



#### Script block3:

```python
#!/usr/bin/env python3
"""
2D incompressible Navier-Stokes around a cylinder on a staggered MAC grid
with immersed-boundary volume penalization. Corrected discrete operators
so G = D.T, L_p = D*G, removed mean-subtraction that broke the pin, and
re-applied penalization after projection to ensure no-slip in the solid.
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
visc = 1.0 / Re   # use visc everywhere per instructions

# Penalization parameter (smaller = stronger enforcement)
eta = 1e-3


def build_grid(x1=3.0, x2=8.0, y1=3.0, y2=3.0, Nx=220, Ny=120):
    Lx = x1 + x2
    Ly = y1 + y2
    dx = Lx / Nx
    dy = Ly / Ny

    # pressure cell centers (Nx x Ny)
    x_p = np.linspace(-x1 + dx/2, x2 - dx/2, Nx)
    y_p = np.linspace(-y1 + dy/2, y2 - dy/2, Ny)

    # u at vertical faces: x faces i=0..Nx (Nx+1), y centers j=0..Ny-1 (Ny)
    x_u = np.linspace(-x1, x2, Nx+1)
    y_u = np.linspace(-y1 + dy/2, y2 - dy/2, Ny)

    # v at horizontal faces: x centers i=0..Nx-1 (Nx), y faces j=0..Ny (Ny+1)
    x_v = np.linspace(-x1 + dx/2, x2 - dx/2, Nx)
    y_v = np.linspace(-y1, y2, Ny+1)

    return dx, dy, x_u, y_u, x_v, y_v, x_p, y_p


# index helpers consistent with Fortran-order flattening (order='F')
def idx_u(i, j, Nx, Ny):
    return i + j * (Nx + 1)


def idx_v(i, j, Nx, Ny):
    return i + j * Nx


def idx_p(i, j, Nx, Ny):
    return i + j * Nx


def build_cylinder_masks(x_u, y_u, x_v, y_v, radius=0.5):
    Xu, Yu = np.meshgrid(x_u, y_u, indexing='ij')
    Xv, Yv = np.meshgrid(x_v, y_v, indexing='ij')

    chi_u = ((Xu**2 + Yu**2) <= radius**2).astype(float)
    chi_v = ((Xv**2 + Yv**2) <= radius**2).astype(float)

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
            # left neighbor
            if i - 1 >= 0:
                left = idx_u(i - 1, j, Nx, Ny)
                rows.append(row); cols.append(left); data.append(1.0 / dx**2)
                center += -1.0 / dx**2
            # right neighbor
            if i + 1 <= Nx:
                right = idx_u(i + 1, j, Nx, Ny)
                rows.append(row); cols.append(right); data.append(1.0 / dx**2)
                center += -1.0 / dx**2

            # down (j-1)
            if j - 1 >= 0:
                down = idx_u(i, j - 1, Nx, Ny)
                rows.append(row); cols.append(down); data.append(1.0 / dy**2)
                center += -1.0 / dy**2
            # up (j+1)
            if j + 1 <= Ny - 1:
                up = idx_u(i, j + 1, Nx, Ny)
                rows.append(row); cols.append(up); data.append(1.0 / dy**2)
                center += -1.0 / dy**2

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


def build_D(Nx, Ny, dx, dy):
    """Build divergence operator D mapping [u_flat; v_flat] -> p_flat (Np x (Nu+Nv)).
    Discretization matches divergence(u,v) = (u_{i+1,j}-u_{i,j})/dx + (v_{i,j+1}-v_{i,j})/dy
    """
    Nu = (Nx + 1) * Ny
    Nv = Nx * (Ny + 1)
    Np = Nx * Ny

    rows = []
    cols = []
    data = []

    for j in range(Ny):
        for i in range(Nx):
            row = idx_p(i, j, Nx, Ny)
            # u faces: left (i,j), right (i+1,j)
            u_left = idx_u(i, j, Nx, Ny)
            u_right = idx_u(i + 1, j, Nx, Ny)
            rows.extend([row, row]); cols.extend([u_right, u_left]); data.extend([1.0 / dx, -1.0 / dx])

            # v faces: bottom (i,j), top (i,j+1)
            v_bottom = idx_v(i, j, Nx, Ny)
            v_top = idx_v(i, j + 1, Nx, Ny)
            # v columns are offset by Nu in the combined vel vector
            rows.extend([row, row]); cols.extend([Nu + v_top, Nu + v_bottom]); data.extend([1.0 / dy, -1.0 / dy])

    D = sp.csr_matrix((data, (rows, cols)), shape=(Np, Nu + Nv))
    return D


def compute_convective(u, v, Nx, Ny, dx, dy):
    # identical convective computation as before (vectorized)
    du_dx = np.zeros_like(u)
    if Nx >= 2:
        du_dx[1:Nx, :] = (u[2:Nx+1, :] - u[0:Nx-1, :]) / (2 * dx)
    du_dx[0, :] = (u[1, :] - u[0, :]) / dx
    du_dx[Nx, :] = (u[Nx, :] - u[Nx - 1, :]) / dx

    du_dy = np.zeros_like(u)
    if Ny >= 3:
        du_dy[:, 1:Ny-1] = (u[:, 2:Ny] - u[:, 0:Ny-2]) / (2 * dy)
    if Ny >= 2:
        du_dy[:, 0] = (u[:, 1] - u[:, 0]) / dy
        du_dy[:, Ny-1] = (u[:, Ny-1] - u[:, Ny-2]) / dy

    v_at_u = np.zeros_like(u)
    if Nx >= 2:
        v_at_u[1:Nx, :] = 0.5 * (v[0:Nx-1, 0:Ny] + v[1:Nx, 0:Ny])
    v_at_u[0, :] = v[0, 0:Ny]
    v_at_u[Nx, :] = v[Nx-1, 0:Ny]

    N_u = u * du_dx + v_at_u * du_dy

    dv_dy = np.zeros_like(v)
    if Ny >= 2:
        dv_dy[:, 1:Ny] = (v[:, 2:Ny+1] - v[:, 0:Ny-1]) / (2 * dy)
        dv_dy[:, 0] = (v[:, 1] - v[:, 0]) / dy
        dv_dy[:, Ny] = (v[:, Ny] - v[:, Ny-1]) / dy

    dv_dx = np.zeros_like(v)
    if Nx >= 3:
        dv_dx[1:Nx-1, :] = (v[2:Nx, :] - v[0:Nx-2, :]) / (2 * dx)
    if Nx >= 2:
        dv_dx[0, :] = (v[1, :] - v[0, :]) / dx
        dv_dx[Nx-1, :] = (v[Nx-1, :] - v[Nx-2, :]) / dx

    u_at_v = np.zeros_like(v)
    if Ny >= 2:
        u_at_v[:, 1:Ny] = 0.25 * (
            u[0:Nx, 0:Ny-1] + u[1:Nx+1, 0:Ny-1] + u[0:Nx, 1:Ny] + u[1:Nx+1, 1:Ny]
        )
    u_at_v[:, 0] = 0.5 * (u[0:Nx, 0] + u[1:Nx+1, 0])
    u_at_v[:, Ny] = 0.5 * (u[0:Nx, Ny-1] + u[1:Nx+1, Ny-1])

    N_v = u_at_v * dv_dx + v * dv_dy

    return N_u, N_v


# consistent flatten/unflatten using Fortran order
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
    Nx = 220
    Ny = 120

    dx, dy, x_u, y_u, x_v, y_v, x_p, y_p = build_grid(x1, x2, y1, y2, Nx=Nx, Ny=Ny)

    # construct masks for cylinder penalization
    chi_u, chi_v = build_cylinder_masks(x_u, y_u, x_v, y_v, radius=0.5)

    # Flatten sizes
    Nu = (Nx + 1) * Ny
    Nv = Nx * (Ny + 1)
    Np = Nx * Ny

    # inlet u-face indices (left boundary)
    inlet_u_indices = [idx_u(0, j, Nx, Ny) for j in range(Ny)]
    inlet_u_set = set(inlet_u_indices)

    # v Dirichlet indices (top/bottom)
    v_zero_indices = set()
    for i in range(Nx):
        v_zero_indices.add(idx_v(i, 0, Nx, Ny))   # bottom v=0
        v_zero_indices.add(idx_v(i, Ny, Nx, Ny))  # top v=0

    # Build Laplacians for viscous terms
    print('Building operators...')
    L_u = build_laplacian_u(Nx, Ny, dx, dy, inlet_u_set)
    L_v = build_laplacian_v(Nx, Ny, dx, dy, v_zero_indices)

    # Build divergence matrix D and gradient G = D.T
    D = build_D(Nx, Ny, dx, dy)
    G = D.T  # gradient operator consistent with divergence

    # Build pressure Laplacian Lp = D * G
    L_p = (D.dot(G)).tocsr()

    # choose a pinned pressure index at outflow-center to minimize bias
    pinned_idx = idx_p(Nx - 1, Ny // 2, Nx, Ny)
    # enforce pin by replacing the pinned row with identity row
    Lp_lil = L_p.tolil()
    Lp_lil.rows[pinned_idx] = [pinned_idx]
    Lp_lil.data[pinned_idx] = [1.0]
    L_p = Lp_lil.tocsc()

    # Time-stepping parameters
    T = 20.0
    Umax = 1.0
    CFL = 0.45
    dt = CFL * min(dx, dy) / Umax
    nsteps = int(floor(T / dt))
    dt = T / nsteps

    print(f'Nx={Nx} Ny={Ny} dx={dx:.4f} dy={dy:.4f} dt={dt:.5f} steps={nsteps}')

    # Build constant matrices for velocity solves: Au = I - dt*visc*L_u + dt*(chi/eta) diag
    Iu = sp.eye(Nu, format='csr')
    Iv = sp.eye(Nv, format='csr')

    chi_u_flat = flatten(chi_u)
    chi_v_flat = flatten(chi_v)

    Pen_u = sp.diags(dt * (chi_u_flat / eta), format='csr')
    Pen_v = sp.diags(dt * (chi_v_flat / eta), format='csr')

    Au = Iu - dt * visc * L_u + Pen_u
    Av = Iv - dt * visc * L_v + Pen_v

    # Enforce Dirichlet rows in Au and Av to identity
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

    # Factorize matrices for fast repeated solves
    print('Factorizing velocity matrices...')
    Au_fac = spla.splu(Au)
    Av_fac = spla.splu(Av)

    # Factorize pressure Poisson (pinning included) once
    print('Factorizing pressure Poisson...')
    Lp_fac = spla.splu(L_p)

    # Initialize fields (Fortran-order shapes)
    u = np.ones((Nx + 1, Ny), order='F') * 1.0
    v = np.zeros((Nx, Ny + 1), order='F')
    p = np.zeros((Nx, Ny), order='F')

    # apply penalization to enforce zero inside cylinder initially
    u = u * (1.0 - chi_u)
    v = v * (1.0 - chi_v)

    # History for AB2 convective term
    N_u, N_v = compute_convective(u, v, Nx, Ny, dx, dy)
    N_u_prev = N_u.copy()
    N_v_prev = N_v.copy()

    inlet_u_flat_indices = inlet_u_indices
    v_zero_flat_indices = list(v_zero_indices)

    # Round-trip check for flatten/unflatten mapping (sanity)
    test_u = np.random.randn(Nx + 1, Ny)
    assert np.allclose(test_u, unflatten_u(flatten(test_u), Nx, Ny))

    print('Starting time-stepping...')
    progress_interval = max(1, nsteps // 10)
    for step in range(1, nsteps + 1):
        N_u, N_v = compute_convective(u, v, Nx, Ny, dx, dy)
        N_u_ex = 1.5 * N_u - 0.5 * N_u_prev
        N_v_ex = 1.5 * N_v - 0.5 * N_v_prev

        rhs_u = flatten(u) - dt * flatten(N_u_ex)
        rhs_v = flatten(v) - dt * flatten(N_v_ex)

        # enforce inlet Dirichlet in RHS
        for idx in inlet_u_flat_indices:
            rhs_u[idx] = 1.0
        for idx in v_zero_flat_indices:
            rhs_v[idx] = 0.0

        # Solve implicit systems
        u_star_flat = Au_fac.solve(rhs_u)
        v_star_flat = Av_fac.solve(rhs_v)

        u_star = unflatten_u(u_star_flat, Nx, Ny)
        v_star = unflatten_v(v_star_flat, Nx, Ny)

        # build combined velocity vector for divergence: [u; v]
        vel_star_combined = np.concatenate([u_star_flat, v_star_flat])

        # compute RHS for Poisson: (1/dt) * D * vel_star
        rhs_p = (1.0 / dt) * (D.dot(vel_star_combined))

        # enforce pinned pressure RHS
        rhs_p[pinned_idx] = 0.0

        # solve Poisson L_p * p_corr = rhs_p
        p_corr_flat = Lp_fac.solve(rhs_p)

        # ensure pinned node equals zero (sanity, no global mean subtraction)
        p_corr_flat[pinned_idx] = 0.0
        p_corr = unflatten_p(p_corr_flat, Nx, Ny)

        # compute pressure gradient at faces via G * p (G = D.T)
        grad_p_faces = G.dot(p_corr_flat)
        dpdx_flat = grad_p_faces[0:Nu]
        dpdy_flat = grad_p_faces[Nu:Nu + Nv]

        dpdx = unflatten_u(dpdx_flat, Nx, Ny)
        dpdy = unflatten_v(dpdy_flat, Nx, Ny)

        # project velocities
        u_new = u_star - dt * dpdx
        v_new = v_star - dt * dpdy

        # re-apply boundary conditions explicitly
        u_new[0, :] = 1.0
        # free-slip (du/dy=0) at top/bottom -> copy neighbor
        u_new[:, 0] = u_new[:, 1]
        u_new[:, Ny - 1] = u_new[:, Ny - 2]
        # outlet zero-gradient
        u_new[Nx, :] = u_new[Nx - 1, :]
        # v top/bottom Dirichlet
        v_new[:, 0] = 0.0
        v_new[:, Ny] = 0.0

        # re-apply penalization (no-slip inside cylinder) to avoid reintroduction of velocity
        u_new = u_new * (1.0 - chi_u)
        v_new = v_new * (1.0 - chi_v)

        # update fields and history
        u = u_new.copy()
        v = v_new.copy()
        p = p_corr.copy()

        N_u_prev = N_u.copy()
        N_v_prev = N_v.copy()

        # compute divergence residual after projection for diagnostics
        div_res = (u[1:Nx+1, :] - u[0:Nx, :]) / dx + (v[:, 1:Ny+1] - v[:, 0:Ny]) / dy
        div_L2 = np.sqrt(np.sum(div_res**2) * dx * dy)
        div_Linf = np.max(np.abs(div_res))

        if (step % progress_interval) == 0 or step == nsteps:
            perc = int(100 * step / nsteps)
            uc = 0.5 * (u[0:Nx, :] + u[1:Nx+1, :])
            vc = 0.5 * (v[:, 0:Ny] + v[:, 1:Ny+1])
            KE = 0.5 * np.sum(uc**2 + vc**2) * dx * dy
            umax = max(np.abs(u).max(), np.abs(v).max())
            print(f'Step {step}/{nsteps} ({perc}%) t={step*dt:.3f} umax={umax:.3f} KE={KE:.5f} div_L2={div_L2:.5e} div_Linf={div_Linf:.5e}')

    # After time loop, compute velocity magnitude at cell centers for plotting
    u_center = 0.5 * (u[0:Nx, :] + u[1:Nx + 1, :])  # shape (Nx, Ny)
    v_center = 0.5 * (v[:, 0:Ny] + v[:, 1:Ny + 1])  # shape (Nx, Ny)
    vel_mag = np.sqrt(u_center**2 + v_center**2)

    # build meshgrid for plotting (pressure cell centers) with indexing='ij'
    Xp, Yp = np.meshgrid(x_p, y_p, indexing='ij')

    fig, ax = plt.subplots(figsize=(10, 4.5))
    cp = ax.contourf(Xp, Yp, vel_mag, levels=50, cmap='RdBu_r')
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




### runtime_outputs

#### Output block1

Building operators...
Nx=220 Ny=120 dx=0.0500 dy=0.0500 dt=0.02252 steps=888
Factorizing velocity matrices...
Factorizing pressure Poisson...

Runtime Error:
Traceback (most recent call last):
  File "C:\Users\admin\AppData\Local\Temp\ipykernel_27508\2786029551.py", line 138, in execute_code_tool
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
Step 88/888 (9%) t=1.982 umax=1.496 KE=34.53532 div_L2=8.52052e-05
Step 176/888 (19%) t=3.964 umax=1.459 KE=35.30345 div_L2=1.95781e-04
Step 264/888 (29%) t=5.946 umax=1.435 KE=35.95091 div_L2=6.23056e-04
Step 352/888 (39%) t=7.928 umax=1.479 KE=36.52476 div_L2=3.17800e-03
Step 440/888 (49%) t=9.910 umax=1.485 KE=37.08953 div_L2=8.52232e-03
Step 528/888 (59%) t=11.892 umax=1.486 KE=37.61715 div_L2=2.21420e-02
Step 616/888 (69%) t=13.874 umax=1.528 KE=37.90682 div_L2=4.13568e-02
Step 704/888 (79%) t=15.856 umax=1.689 KE=38.08575 div_L2=8.53713e-02
Step 792/888 (89%) t=17.838 umax=1.621 KE=37.91326 div_L2=9.02908e-02
Step 880/888 (99%) t=19.820 umax=1.463 KE=37.34590 div_L2=7.61843e-02
Step 888/888 (100%) t=20.000 umax=1.461 KE=37.41149 div_L2=8.31657e-02
Saved figure to velocity_magnitude.png



#### Output block3

Building operators...
Nx=220 Ny=120 dx=0.0500 dy=0.0500 dt=0.02252 steps=888
Factorizing velocity matrices...
Factorizing pressure Poisson...
Starting time-stepping...
Step 88/888 (9%) t=1.982 umax=1.434 KE=32.03667 div_L2=7.62211e-01 div_Linf=4.44200e+00
Step 176/888 (19%) t=3.964 umax=1.401 KE=32.29324 div_L2=8.16461e-01 div_Linf=4.45055e+00
Step 264/888 (29%) t=5.946 umax=1.383 KE=32.67805 div_L2=9.10972e-01 div_Linf=4.63536e+00
Step 352/888 (39%) t=7.928 umax=1.424 KE=33.07357 div_L2=9.92096e-01 div_Linf=4.51905e+00
Step 440/888 (49%) t=9.910 umax=1.431 KE=33.42838 div_L2=1.06152e+00 div_Linf=4.48036e+00
Step 528/888 (59%) t=11.892 umax=1.442 KE=33.71274 div_L2=1.11784e+00 div_Linf=4.98673e+00
Step 616/888 (69%) t=13.874 umax=1.465 KE=33.81210 div_L2=1.14418e+00 div_Linf=5.06025e+00
Step 704/888 (79%) t=15.856 umax=1.484 KE=33.65972 div_L2=1.20230e+00 div_Linf=5.23463e+00
Step 792/888 (89%) t=17.838 umax=1.529 KE=33.27147 div_L2=1.29817e+00 div_Linf=5.83491e+00
Step 880/888 (99%) t=19.820 umax=1.511 KE=32.93176 div_L2=1.35217e+00 div_Linf=5.55204e+00
Step 888/888 (100%) t=20.000 umax=1.519 KE=32.92984 div_L2=1.35421e+00 div_Linf=5.54280e+00
Saved figure to velocity_magnitude.png




### review_decision
revise

### review_comments
	Current Stage [A/3]
1. MAJOR

Category: runtime
Issue: Broadcasting error (ValueError: operands could not be broadcast together with shapes (221,120) (120,221)) when applying the cylinder mask to the u field: u = u * (1 - chi_u).

Recommendation: The mask arrays are created with inconsistent axis ordering. In build_cylinder_masks you call np.meshgrid(..., indexing='xy'), which yields arrays shaped (Ny, Nx+1) whereas the code expects masks with shape (Nx+1, Ny). Change meshgrid calls to indexing='ij' so Xu, Yu (and Xv, Yv) have shape (len(x), len(y)). Also audit all meshgrid calls (including plotting) so all arrays use the same indexing convention.


2. MAJOR

Category: structure
Issue: Inconsistent flatten/reshape array ordering across the code (mix of C-order ravel/reshape and transposes). This makes the code fragile and will produce subtle shape bugs if any meshgrid/array creation uses a different ordering.

Recommendation: Adopt a single, explicit array ordering convention everywhere. Two safe options: (a) use meshgrid(..., indexing='ij') and keep C-order ravel/reshape (default) consistently; or (b) follow the problem hint and use order='F' for all ravel/reshape and build meshgrids with indexing='ij'. Concretely: update flatten/unflatten to use a single order argument (e.g. flatten(arr) -> arr.ravel(order='F') and reshape(..., order='F')), and update all .ravel()/.reshape() calls to specify order. Remove ad-hoc .T usage by making reshape shapes match the intended layout.


3. MINOR

Category: structure
Issue: compute_convective contains Python-level nested loops to interpolate staggered velocities (loops over i,j). For the chosen grid sizes (Nx=220, Ny=120) these loops are slow and will dominate runtime.

Recommendation: Vectorize the interpolation and derivative stencils using NumPy slicing and averaging. For example, compute v_at_u by averaging the two nearest v slices with appropriate bounds handling instead of looping. Similarly replace the u_at_v loop with vectorized operations. This will dramatically reduce runtime and is straightforward with the MAC grid indexing.


4. MINOR

Category: accuracy
Issue: Pressure nullspace handling: code pins a single pressure node (pinned_idx=0). While this removes the singularity, pinning a corner cell may introduce a small bias and interacts with the simple mirrored Neumann approximations used in L_p.

Recommendation: Either pin a pressure node in a physically neutral location (e.g. a cell on the outflow boundary center) or project the pressure to zero mean after solving (i.e., subtract mean(p) inside the solver) to avoid introducing a location-dependent bias. Also verify divergence residuals after projection to confirm incompressibility.


5. MINOR

Category: accuracy
Issue: Discrete operators and boundary treatment are simplistic: the Laplacian construction uses mirrored ghost approximations and Dirichlet enforcement by row replacement; the Poisson operator and gradient/divergence may not be strict discrete adjoints on the MAC mesh. This can degrade mass conservation and pressure-velocity coupling.

Recommendation: Revisit the discrete operator assembly to ensure consistency with MAC staggering: build divergence and gradient as discrete adjoints (transpose) and construct the Poisson operator as D * (M^{-1}) * G or directly the MAC-compatible discrete Laplacian. If time is limited, at minimum test and report L2-norm of divergence after projection and tune boundary approximations. If accuracy/stability issues remain, consider using a standard MAC Poisson assembly or a ready-made solver (e.g. PyAMG) for the pressure solve.



	Current Stage [B/3]
1. MAJOR

Category: accuracy
Issue: Pressure mean subtraction after solving pinned Poisson breaks the pinned constraint and can introduce inconsistency in the projection.

Recommendation: Do not subtract the pressure mean after solving the Poisson if you pin a pressure node (or if you must remove the mean, re‑enforce the pinned node afterwards). Currently you set L_p[pinned_idx,:]=I and rhs[pinned_idx]=0 so the direct solve returns p_corr[pinned_idx]=0. Immediately subtracting the global mean breaks that constraint and reintroduces an inconsistent pressure offset that corrupts the projection. Remove the mean-subtraction step or restore the pinned value (p_corr_flat[pinned_idx]=0) after subtracting the mean.


2. MAJOR

Category: accuracy
Issue: Inconsistent discrete divergence / gradient / Poisson discretizations causing a growing divergence residual (div_L2 ~ 1e-1).

Recommendation: Make the discrete divergence/gradient operators algebraically consistent and reflect the chosen boundary conditions in the Poisson assembly. The current approach builds the Poisson from a simple 5‑point stencil and computes gradients with ad‑hoc one‑sided zeros at faces; this mismatch is the likely reason the divergence L2 grows to O(1e-1). Fix by assembling discrete D (divergence) and G (gradient) sparse matrices on the MAC grid so that G = -D^T (up to the pinned row) and build the Poisson as D * G (or use the adjoint relation when imposing BCs). Also implement the correct boundary conditions for the pressure Poisson (Neumann at outlet): modify boundary rows appropriately or include ghost/one‑sided stencils so the RHS exactly matches net face fluxes.


3. MAJOR

Category: accuracy
Issue: Penalization is treated implicitly for tentative velocities but not for the projection step, which may reintroduce nonzero velocities inside the solid.

Recommendation: The projection step does not include the penalization term: you included penalization implicitly in the implicit momentum solve (Au, Av) but the subsequent pressure projection modifies u_star without penalization and can reintroduce velocity inside the cylinder. Enforce the no‑slip after projection by either (a) including penalization in the projection operator (modify the projection equation to account for the diagonal penalization), or (b) simply re-apply the penalization (set u_new[chi_u==1]=0 and v_new[chi_v==1]=0) immediately after projection as a pragmatic fix. The latter is cheap and consistent with the method requirement to rely on penalization for no‑slip.


4. MINOR

Category: structure
Issue: Mixed handling of pressure nullspace (pin + mean subtraction) and limited diagnostics make it hard to verify correctness of the projection.

Recommendation: Correct the pressure pinning / boundary handling decisions and improve numerical diagnostics: either (i) keep the current pinned row strategy but do not apply mean subtraction (or re-pin afterwards), or (ii) remove the explicit pin and instead enforce the pressure nullspace by projecting out the mean with the pinned row preserved. Add a check that the pinned entry after any post-processing is equal to the intended value (e.g. zero). For diagnostics, print the maximum absolute divergence (L_inf) in addition to the L2 norm to better monitor local divergence spikes.


5. MINOR

Category: structure
Issue: Opportunities to improve robustness, performance and maintainability (factorization choice, reapplication of no‑slip, small sanity checks).

Recommendation: Small robustness and performance notes: (a) After fixing correctness issues consider using sparse symmetric factorization (Cholesky) or iterative solvers with AMG for the Poisson to reduce memory/time for larger grids; (b) re-apply body velocities (no-slip) after each full step as a safety net until you add penalization to the projection; (c) explicitly document and assert the Fortran-order/idx mapping (you appear consistent, but an assert test on shapes / forward/backward mapping would help future maintenance).



	Current Stage [C/3]
1. MAJOR

Category: accuracy
Issue: Re-applying penalization (zeroing velocities inside the cylinder) after projection reintroduces divergence and breaks mass conservation.

Recommendation: Remove the explicit zeroing of velocities inside the cylinder after the projection step. That line (u_new = u_new * (1.0 - chi_u) and similarly for v) reintroduces a local jump in velocity inside the solid and therefore a large divergence immediately after the projection; this is the most likely cause of the large div_L2/div_Linf values reported (~1.3 / ~5.5). You already include the penalization implicitly in the momentum matrices (Pen_u, Pen_v) which is the correct way to enforce no-slip without breaking the projection. Either (a) do not re-apply the mask after projection and rely solely on the implicit penalization, or (b) if you must force exact zero inside the solid, incorporate that constraint consistently into the projection (e.g. modify the Poisson RHS to exclude solid cells or enforce mass-conservation only in the fluid region), or (c) include the penalization in the projection operator so that zeroing inside the solid does not produce divergence. After this fix re-check the divergence (report it in fluid-only cells) and the physical fields.


2. MAJOR

Category: structure
Issue: Pressure pinning replaces only a row; the column is left unchanged which yields an inconsistent asymmetric modification and can bias/complicate the Poisson solve.

Recommendation: Modify the pressure pinning step so the linear system remains consistent and numerically well-conditioned. Currently you replace only the pinned row of L_p with an identity row prior to factorization. This makes L_p non-symmetric (and not strictly the usual row+column pin), which may be acceptable for LU but is not the standard approach. Replace the pinned row and the pinned column with zeros and set the diagonal at pinned_idx to 1.0 (i.e. zero out row and column, set L_p[pinned_idx,pinned_idx]=1). Also set rhs[pinned_idx]=0. This keeps the system consistent and avoids subtle bias in the factorization/solve. After changing this, re-factorize L_p and verify the Poisson solve produces the pinned value exactly.


3. MINOR

Category: structure
Issue: Dirichlet enforcement is applied by setting rows to identity but columns are left intact; this produces non-symmetric matrices and potential inconsistencies.

Recommendation: When you enforce Dirichlet boundary conditions in Au and Av by replacing rows with identity, also clear the corresponding columns (set them to zero) so the matrix representation exactly corresponds to the algebraic elimination of those DOFs. Leaving columns intact makes the matrices non-symmetric and can produce non-physical couplings; clearing columns (or assemble reduced systems eliminating BC DOFs) yields cleaner and more robust factorization. Alternatively, assemble the linear systems with BCs eliminated (smaller systems) to save memory and improve stability.


4. MINOR

Category: accuracy
Issue: Divergence diagnostics are computed over the whole domain including solid cells, which is misleading for an immersed-boundary solver.

Recommendation: Compute and report divergence diagnostics only over the fluid cells (i.e. exclude/ignore cells fully inside the solid mask when reporting div_L2 and div_Linf). Because you are using an immersed method, diagnostics over the full domain are misleading: the sharp zeroing of velocities at the immersed boundary (or subtle penetration) will produce local large divergence that is not representative of mass conservation in the fluid. Create a fluid-cell mask (1 - chi_p where chi_p is the cell-centered indicator) and compute L2/Linf on that region for a physically meaningful diagnostic.


5. MINOR

Category: runtime
Issue: Opportunities for solver and performance improvements; small maintainability notes.

Recommendation: A few additional suggestions to improve robustness and maintainability: (1) For the Poisson solve, consider using a symmetric positive-definite Poisson solver (conjugate gradient with an AMG/preconditioner) for larger problems to reduce memory/time vs. full sparse LU. (2) Replace any remaining list membership checks in inner assembly loops with sets (you already use a set in one place) to guarantee O(1) checks; you already do this in places but verify consistency. (3) After implementing the changes above, run a short test and report divergence (fluid-only), energy behavior and check that velocity inside solid is sufficiently small without manual masking. These are optional but will improve code quality and performance.





### iteration_count
3

### rev_image_description
	Current Stage [A/3]


	Current Stage [B/3]


	Current Stage [C/3]




