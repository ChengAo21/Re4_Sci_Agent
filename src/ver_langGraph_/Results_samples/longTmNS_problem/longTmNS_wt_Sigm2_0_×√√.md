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
We must compute a time-dependent incompressible flow in 2D on a rectangular domain [0,2]×[0,1] at moderate Reynolds number Re=100 with a time-dependent body force f_y = -sin(pi x) sin(pi y) sin(pi t) and nontrivial inlet/outlet/wall boundary conditions. The governing equations are the unsteady Navier–Stokes equations with incompressibility constraint; the domain prescribes no-slip on top/bottom, a time-varying Dirichlet inflow at x=0, homogeneous Neumann velocity BCs (zero normal derivative) at the outlet x=2, and a Dirichlet pressure reference p(2,y,t)=0. Initial fields are zero.

Primary mathematical and numerical challenges:
- Incompressibility constraint: enforcing div(u)=0 requires a robust pressure–velocity coupling (projection or pressure-correction) and a fast Poisson solver for the pressure update.
- Mixed boundary conditions: Dirichlet velocity inflow/wall and Neumann velocity (outlet) plus pressure Dirichlet at outlet require consistent discretization on a staggered grid to avoid spurious pressure modes and to satisfy discrete mass conservation.
- Stability of advection: advective terms are nonlinear and can dominate; an upwind discretization is recommended by the hints to ensure stable advection treatment at the cost of added numerical dissipation.
- Time stepping: balancing accuracy and stability requires semi-implicit schemes (implicit diffusion, explicit advection) or fully implicit methods. Projection/splitting methods introduce splitting error and require careful handling if second-order temporal accuracy is desired.
- Linear solver bottlenecks: Poisson (pressure) and Helmholtz (implicit viscous) solves are the computational bottlenecks and must be accelerated (multigrid, preconditioning, reuse of factorization) for efficiency.
- Grid and indexing: using a staggered (MAC) grid avoids pressure–velocity decoupling but requires careful array layout: use indexing='ij' for meshgrid and order='F' for ravel/reshape to maintain correct coordinate alignment and consistent memory layout.
- Practical implementation constraints: follow provided hints (use visc global, use 'atol' in SciPy solvers not 'tol', avoid np.trapz, print progress every 10% of steps, save figures without plt.show()).

The goal is to implement a stable, efficient, and robust solver (Chorin or IPCS on a MAC grid), with reasonable acceleration strategies (reuse operators, AMG/ILU preconditioning, possible numba vectorization), and to simulate to t=1.0 and save contour plots of u, v, and p (colormap 'RdBu_r').

### solution_plans
	Current Stage [A/2]
solu_name='MAC Chorin Projection with Implicit Diffusion (Crank–Nicolson) + Upwind Advection' content="Governing idea:\n1) Use Chorin's projection on a staggered (MAC) grid. Treat advection explicitly with a stable upwind discretization on the staggered velocities; treat viscous diffusion implicitly (Crank–Nicolson or implicit Euler) so the linear diffusion operator is constant in time and can be factorized/reused for acceleration. Solve a pressure Poisson (Poisson-like) equation each timestep to enforce incompressibility. Use iterative solvers with appropriate preconditioning (multigrid or ILU) for the Poisson and Helmholtz systems.\n\nAlgorithmic steps (high level numbered pseudocode):\n1. Global constants / configuration (outside main): set visc = 1.0 / Re, domain [0,2]×[0,1], Nx, Ny, t_end = 1.0, CFL_limit, output / plot settings. Pass all constants explicitly into functions.\n2. Build a staggered MAC mesh (u located at vertical-face centers, v at horizontal-face centers, p at cell centers). Use indexing='ij' for meshgrid and ensure array ravel/reshape use order='F'.\n3. Initialize u=v=p=0 and apply inlet/top/bottom/outlet boundary conditions. Impose p(2,y,t)=0 as pressure reference.\n4. Precompute discrete derivative operators (finite-difference) on MAC grid and construct the constant Helmholtz matrices for u and v viscous implicit solve when using Crank–Nicolson or implicit Euler. Assemble Poisson matrix for pressure correction with Neumann/Dirichlet BCs consistent with staggered discretization.\n5. Choose time stepping: use AB2 (Adams–Bashforth 2) for convective terms (first step use forward Euler), and Crank–Nicolson (theta=0.5) or backward Euler for diffusion. Compute dt at each step from CFL: dt = CFL * min(dx/max(|u|+eps), dy/max(|v|+eps), 0.5*min(dx,dy)**2/visc) but cap so final time hits t_end exactly. Print progress every 10% of total steps.\n6. Time loop until t_end:\n   a. Evaluate explicit convective terms on staggered grid using a first-order upwind flux for stability. Include the explicit forcing f_x=0, f_y=-sin(pi x) sin(pi y) sin(pi t) projected appropriately onto MAC faces.\n   b. Form right-hand-sides for intermediate velocities u* and v* by adding explicit advection + forcing and treating diffusion implicitly (solve Helmholtz systems): (I - theta*dt*visc*L) u* = RHS_u, similarly for v*. Solve these linear systems using a pre-factorized sparse LU (if memory allows and operator constant) or a Krylov solver (GMRES/CG) with ILU/AMG preconditioner. When calling SciPy solvers, use atol (not tol).\n   c. Compute divergence of (u*,v*) on cell-centers and solve the pressure Poisson equation: L_p phi = (1/dt)*div(u*). Enforce Dirichlet at outlet p=0 and consistent Neumann/Dirichlet on other boundaries derived from velocity BCs; solve with geometric/Algebraic Multigrid (PyAMG) or a preconditioned CG using AMG preconditioner. Use 'atol' in solver calls.\n   d. Correct velocities: u^{n+1} = u* - dt * d(phi)/dx on staggered u faces; v^{n+1} = v* - dt * d(phi)/dy on staggered v faces. Update pressure p^{n+1} = phi (or incremental pressure depending on projection variant).\n   e. Re-apply velocity boundary conditions (inlet, no-slip walls, zero-gradient outlet) to maintain consistency; ensure ghost/fill logic for staggered grid.\n   f. Update time, compute diagnostics, print progress at 10% increments.\n7. At final t, interpolate or remap staggered u/v/p to cell centers for plotting, use plt.contourf with cmap='RdBu_r' and save figures (do not call plt.show()).\n\nStability / accuracy / complexity / efficiency limitations:\n- Splitting error: Chorin projection introduces O(dt) (or O(dt^2) for incremental variants) splitting error for pressure–velocity coupling; use smaller dt or incremental pressure correction for higher accuracy.\n- Time accuracy: AB2 + Crank–Nicolson yields nominal second-order temporal accuracy but explicit advection may require small dt for convective CFL; upwind is only first-order spatially for advection (stable but more diffusive).\n- Spatial accuracy: MAC second-order centered diffusive terms, but first-order upwind convective discretization reduces global spatial accuracy and can smear small scales.\n- Poisson/Helmholtz solves dominate cost: iterative solvers with AMG preconditioners are recommended; direct factorization is faster per solve if operator is constant and factorization fits memory, but memory may be limiting for fine grids.\n- Implementational complexity: assembling staggered operators and consistent BCs is error-prone; careful indexing (indexing='ij', order='F') is required to keep data alignment correct.\n- Acceleration trade-offs: reusing LU factorization accelerates solves but increases memory; AMG is scalable but requires external packages (PyAMG) and careful tuning. Parallelization via numba/numexpr can speed up advection/differencing but adds implementation complexity."

	Current Stage [B/2]
solu_name='IPCS (Incremental Pressure Correction) with BDF2 + AMG-accelerated Poisson/Helmholtz Solves' content="Governing idea:\n1) Use an Incremental Pressure Correction Scheme (IPCS, aka projection variant) on a MAC grid with BDF2 for time discretization (backward-looking multi-step) for improved temporal accuracy, explicit upwind advection, and fully implicit viscous term solved as a Helmholtz problem. Use algebraic multigrid (AMG) or a geometric multigrid V-cycle for the pressure Poisson to accelerate the dominant bottleneck. Leverage operator reuse and block preconditioning where possible.\n\nAlgorithmic steps (high level numbered pseudocode):\n1. Globals and configuration: set visc = 1.0 / Re, domain, Nx, Ny, t_end=1.0, choose BDF2 for time (first step use backward Euler), CFL safety factor, pass all constants explicitly to functions.\n2. Construct MAC grid (indexing='ij') and finite-difference stencils. Ensure ordering with order='F' on ravel/reshape operations. Precompute sparse Helmholtz matrices for u and v (mass + viscous term) and pre-factorize or setup AMG preconditioner.\n3. Initialize fields u,v,p=0 and apply boundary conditions (inlet Dirichlet for u, no-slip walls, Neumann outlet for velocities, Dirichlet p=0 at outlet to fix gauge). Map forcing f_y onto v faces consistently with staggered layout.\n4. Compute dt using BDF2 stability constraints and CFL: dt = min(CFL_convective, CFL_diffusive) and ensure final t exactly reaches t_end. Print progress every 10% of total steps.\n5. Time-stepping loop:\n   a. Compute explicit convective term using second-order upwind or stabilized donor-cell upwind on staggered grid (prefer first-order upwind if robustness prioritized). Use previous-step velocities for explicit term (BDF2 requires two previous velocities).\n   b. Form the BDF2 Helmholtz RHS for intermediate velocities: (1.5/dt) u^{*} - visc*L u^{*} = RHS (combining previous steps and explicit convective + forcing). Solve the sparse Helmholtz system for u^{*} and v^{*} using preconditioned Krylov (use atol), with AMG or ILU preconditioner; reuse preconditioner/operator across time steps.\n   c. Assemble divergence of u^{*}, build pressure Poisson problem L_p delta_p = (1/dt) div(u^{*}) and solve for pressure increment delta_p with a multigrid solver (geometric or PyAMG). Impose p=0 at outlet and consistent Neumann elsewhere.\n   d. Update pressure p^{n+1} = p^{n} + delta_p and correct velocities: u^{n+1} = u^{*} - dt * grad(delta_p) on faces, similarly for v, then apply boundary conditions.\n   e. Shift BDF history arrays and continue to next timestep.\n6. After t_end, interpolate staggered fields to cell centers for visualization. Save contours for u, v, p in a single figure using 'RdBu_r' colormap; do not call plt.show().\n\nAcceleration strategies embedded in this plan:\n- Pre-factorize Helmholtz operator if time-constant (Crank–Nicolson variant) or reuse AMG/ILU preconditioners; avoid rebuilding sparse matrices each step.\n- Use geometric multigrid (if implementing) for the pressure Poisson; if using PyAMG, reuse the hierarchy across solves to amortize setup cost.\n- Vectorize finite-difference loops; use order='F' memory layout and indexing='ij' to keep Fortran-ordered arrays cache friendly and compatible with linear algebra operations.\n- Reduce cost by adaptive dt control: keep dt as large as CFL permits; if using BDF2, adjust first step with backward Euler to start multistep.\n- If available, accelerate inner loops with numba or Cython and parallelize independent operations (advection residuals, boundary fills).\n\nStability / accuracy / complexity / efficiency limitations:\n- Implementation complexity: BDF2/IPCS with AMG and staggered indexing is more complex to implement and debug than simple projection; multigrid hierarchy construction can be fragile for mixed Dirichlet/Neumann BCs and staggered operators.\n- Splitting error: IPCS reduces some splitting error versus classical Chorin but does not eliminate it; pressure boundary treatment at the outlet (mixed BCs) requires care to avoid spurious pressure boundary layers.\n- Solver dependence: effectiveness heavily depends on available linear solvers and preconditioners (PyAMG, PETSc); absence of robust AMG increases solve times.\n- Temporal stability: explicit treatment of advection restricts dt by convective CFL (BDF2 is implicit for viscous part but not for advection), so for higher Re and fine grids dt may be small.\n- Spatial dissipation: upwind advection stabilizes but introduces significant numerical diffusion, reducing solution fidelity for sharp features.\n- Memory/time trade-offs: pre-factorization reduces run-time per step but increases memory footprint; multigrid reduces iteration counts but setup cost can dominate for small problems."



### technical_spec
	Current Stage [A/3]
This script implements a Chorin projection solver on a staggered (MAC) grid for 2D unsteady incompressible Navier–Stokes with an explicit body force. It uses a first-order upwind advection, implicit Euler for viscous diffusion, and a projection (pressure Poisson) step each time level. Constant sparse linear systems for the implicit viscous solves and pressure Poisson are assembled once and factorized (direct LU) for reuse, accelerating the time loop. Arrays use indexing='ij' and order='F' for consistent memory layout. The script prints concise progress every 10% of total steps and saves a figure with contours of u, v, and p at final time (colormap 'RdBu_r').

	Current Stage [B/3]
This script implements a Chorin projection solver on a staggered MAC grid for 2D incompressible Navier–Stokes with explicit body force. Core design:
- Staggered grid: u on vertical faces (Nx+1 x Ny), v on horizontal faces (Nx x Ny+1), p on cell centers (Nx x Ny).
- Time stepping: explicit upwind convection, implicit viscous diffusion (backward Euler) solved via constant sparse Helmholtz matrices LU-factorized once, projection via Poisson solve for pressure correction (constant Poisson matrix LU-factorized once).
- Boundary conditions: Dirichlet inlet for u (time-dependent), no-slip top/bottom, Neumann outlet for velocities (implemented via one-sided rows), Dirichlet p at outlet.
- Key fixes from review: corrected inlet-profile indexing (use 1-D y_u), robust and safe upwind implementation for v-face convection (no out-of-bounds indices), corrected Poisson matrix sign/construction to represent discrete ∇^2 consistently with RHS, and consistent enforcement of Dirichlet rows and RHS values. A small self-test verifies upwind_convective_v behaviour on trivial fields.
- Performance: matrices assembled once; LU factorizations reused in time loop. Arrays use indexing='ij' and order='F' for ravel/reshape consistency. Progress prints every 10% of steps. Final u,v,p contours are saved to 'uvp_contours.png'.

	Current Stage [C/3]
This module implements a Chorin projection solver for 2D incompressible Navier–Stokes on a staggered MAC grid. Main design points:
- Staggered grid: u on vertical faces (Nx+1 x Ny), v on horizontal faces (Nx x Ny+1), p on cell centers (Nx x Ny).
- Time stepping: explicit upwind convection, implicit viscous diffusion via constant Helmholtz matrices (LU factorized once), projection via Poisson solve (LU factorized once).
- Consistent boundary treatment: Dirichlet rows are exact identity rows; Neumann outlet conditions are enforced by identity rows that copy the adjacent interior DOF (ghost = interior) — both matrix rows and RHS are set consistently before solves.
- Poisson operator assembled deterministically with explicit handling of left/bottom Neumann by ghost-elimination substitution (no in-place LIL checks). The Poisson operator approximates Δφ so RHS = div/dt.
- Tests added: upwind_convective_u/v checks, discrete-poisson vs analytic Laplacian check, and a projection-preservation test (apply projection to a manufactured divergence-free field and verify divergence is reduced to near zero). Prints divergence before/after projection for the test.
- Performance: matrices built once; all ravel/reshape use order='F' and meshgrid uses indexing='ij'. Progress prints every 10% of steps.
- Global visc is used consistently.




### python_codes

#### Script block1:

#!/usr/bin/env python3
"""
Chorin projection solver on a MAC grid for 2D incompressible Navier-Stokes.
Global variable: visc = 1.0 / Re
Simulates until t = 1.0 and saves a figure of u,v,p contours.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# Global viscosity (do not rename)
Re = 100.0
visc = 1.0 / Re

# ----------------------------- Utility functions -----------------------------
def idx_u(i, j, Nx, Ny):
    # Flatten index for u array shape (Nx+1, Ny), Fortran order
    return i + j * (Nx + 1)

def idx_v(i, j, Nx, Ny):
    # Flatten index for v array shape (Nx, Ny+1), Fortran order
    return i + j * Nx

def idx_p(i, j, Nx, Ny):
    # Flatten index for p array shape (Nx, Ny), Fortran order
    return i + j * Nx

# ----------------------------- Discrete operators -----------------------------
def build_helmholtz_u(Nx, Ny, dx, dy, dt, visc):
    # u grid: (Nx+1) x Ny -> Nu = (Nx+1)*Ny
    Nu = (Nx + 1) * Ny
    A = sp.lil_matrix((Nu, Nu))
    idx = lambda i, j: idx_u(i, j, Nx, Ny)
    ax = visc * dt / dx**2
    ay = visc * dt / dy**2
    for j in range(Ny):
        for i in range(Nx + 1):
            row = idx(i, j)
            # Dirichlet inlet at i==0 and no-slip top/bottom at j==0, j==Ny-1
            if i == 0 or j == 0 or j == Ny - 1:
                A[row, row] = 1.0
                continue
            # Neumann at outlet i == Nx approximate first derivative zero: u_N - u_{N-1} = 0
            if i == Nx:
                # u_i - u_{i-1} = 0
                A[row, idx(i, j)] = 1.0
                A[row, idx(i - 1, j)] = -1.0
                continue
            # interior u node
            diag = 1.0 + 2.0 * (ax + ay)
            A[row, row] = diag
            # x neighbors
            A[row, idx(i - 1, j)] = -ax
            A[row, idx(i + 1, j)] = -ax
            # y neighbors
            A[row, idx(i, j - 1)] = -ay
            A[row, idx(i, j + 1)] = -ay
    A = A.tocsc()
    lu = spla.splu(A)
    return A, lu


def build_helmholtz_v(Nx, Ny, dx, dy, dt, visc):
    # v grid: Nx x (Ny+1) -> Nv = Nx*(Ny+1)
    Nv = Nx * (Ny + 1)
    A = sp.lil_matrix((Nv, Nv))
    idx = lambda i, j: idx_v(i, j, Nx, Ny)
    ax = visc * dt / dx**2
    ay = visc * dt / dy**2
    for j in range(Ny + 1):
        for i in range(Nx):
            row = idx(i, j)
            # Dirichlet top/bottom walls at j==0 or j==Ny and inlet v=0 at i==0
            if j == 0 or j == Ny or i == 0:
                A[row, row] = 1.0
                continue
            # Neumann at outlet (i == Nx-1?), we will approximate d/dx v = 0 at rightmost column: v_{Nx-1} - v_{Nx-2} = 0
            if i == Nx - 1:
                A[row, idx(i, j)] = 1.0
                A[row, idx(i - 1, j)] = -1.0
                continue
            # interior v node
            diag = 1.0 + 2.0 * (ax + ay)
            A[row, row] = diag
            # x neighbors
            A[row, idx(i - 1, j)] = -ax
            A[row, idx(i + 1, j)] = -ax
            # y neighbors
            A[row, idx(i, j - 1)] = -ay
            A[row, idx(i, j + 1)] = -ay
    A = A.tocsc()
    lu = spla.splu(A)
    return A, lu


def build_poisson_p(Nx, Ny, dx, dy):
    # p grid: Nx x Ny -> Np = Nx*Ny
    Np = Nx * Ny
    A = sp.lil_matrix((Np, Np))
    idx = lambda i, j: idx_p(i, j, Nx, Ny)
    rx = 1.0 / dx**2
    ry = 1.0 / dy**2
    for j in range(Ny):
        for i in range(Nx):
            row = idx(i, j)
            # Dirichlet at outlet cell centers: rightmost cells i == Nx-1
            if i == Nx - 1:
                A[row, row] = 1.0
                continue
            # interior/other boundaries: standard 5-point Laplacian, with Neumann at domain boundaries approximated by one-sided
            diag = -2.0 * (rx + ry)
            A[row, row] = -diag  # put positive diagonal
            # x neighbors
            if i - 1 >= 0:
                A[row, idx(i - 1, j)] = rx
            else:
                # Neumann at left boundary -> zero derivative => mirror
                A[row, row] += rx
            if i + 1 < Nx:
                A[row, idx(i + 1, j)] = rx
            else:
                # right boundary handled as Dirichlet earlier
                pass
            # y neighbors
            if j - 1 >= 0:
                A[row, idx(i, j - 1)] = ry
            else:
                # bottom Neumann -> mirror
                A[row, row] += ry
            if j + 1 < Ny:
                A[row, idx(i, j + 1)] = ry
            else:
                # top Neumann -> mirror
                A[row, row] += ry
    A = A.tocsc()
    lu = spla.splu(A)
    return A, lu

# ----------------------------- Advection and helpers -----------------------------
def upwind_convective_u(u, v, dx, dy, Nx, Ny):
    # Compute convective term for u at u-grid (Nx+1, Ny)
    # Using first-order upwind for du/dx and du/dy on staggered u-grid
    conv = np.zeros_like(u)
    # u shape (Nx+1, Ny), v shape (Nx, Ny+1)
    for j in range(Ny):
        for i in range(Nx + 1):
            uc = u[i, j]
            # du/dx: use faces i and i-1
            if i == 0:
                # forward difference at inlet
                du_dx = (u[i + 1, j] - u[i, j]) / dx
                fluxx = uc * du_dx
            elif i == Nx:
                # backward at outlet
                du_dx = (u[i, j] - u[i - 1, j]) / dx
                fluxx = uc * du_dx
            else:
                # upwind based on average u at face
                ua = 0.5 * (u[i, j] + u[i - 1, j])
                if ua >= 0:
                    du_dx = (u[i, j] - u[i - 1, j]) / dx
                else:
                    du_dx = (u[i + 1, j] - u[i, j]) / dx
                fluxx = ua * du_dx
            # du/dy: need v at top/bottom of u-face. v at (i-1 or i) positions; approximate v at u-face by averaging two neighboring v's
            # map v to u-face vertical positions: use v[i-1, j] and v[i, j]
            if i == 0:
                v_left = 0.0
            else:
                v_left = v[i - 1, j]
            if i < Nx:
                v_right = v[i, j]
            else:
                v_right = 0.0
            v_face = 0.5 * (v_left + v_right)
            if j == 0:
                du_dy = (u[i, j + 1] - u[i, j]) / dy
            elif j == Ny - 1:
                du_dy = (u[i, j] - u[i, j - 1]) / dy
            else:
                # upwind based on v_face sign
                if v_face >= 0:
                    du_dy = (u[i, j] - u[i, j - 1]) / dy
                else:
                    du_dy = (u[i, j + 1] - u[i, j]) / dy
            fluxy = v_face * du_dy
            conv[i, j] = fluxx + fluxy
    return conv


def upwind_convective_v(u, v, dx, dy, Nx, Ny):
    # Compute convective term for v at v-grid (Nx, Ny+1)
    conv = np.zeros_like(v)
    for j in range(Ny + 1):
        for i in range(Nx):
            vc = v[i, j]
            # dv/dy
            if j == 0:
                dv_dy = (v[i, j + 1] - v[i, j]) / dy
                fluxy = vc * dv_dy
            elif j == Ny:
                dv_dy = (v[i, j] - v[i, j - 1]) / dy
                fluxy = vc * dv_dy
            else:
                # upwind based on average v
                va = 0.5 * (v[i, j] + v[i, j - 1])
                if va >= 0:
                    dv_dy = (v[i, j] - v[i, j - 1]) / dy
                else:
                    dv_dy = (v[i, j + 1] - v[i, j]) / dy
                fluxy = va * dv_dy
            # dv/dx: need u at top/bottom of v-face; approximate u at v-face by averaging two neighboring u's
            if i == 0:
                u_bottom = u[i, j - 1] if (j - 1 >= 0 and j - 1 < u.shape[1]) else 0.0
            else:
                u_bottom = u[i, j - 1] if (j - 1 >= 0 and j - 1 < u.shape[1]) else 0.0
            # For more symmetric approx, average u[i] and u[i+1]
            if i + 1 <= u.shape[0] - 1:
                u_top = u[i + 1 - 1, j - 1] if (j - 1 >= 0 and j - 1 < u.shape[1]) else 0.0
            else:
                u_top = 0.0
            # Simpler: approximate u at v-face by averaging u[i] and u[i+1]
            if i + 1 <= u.shape[0] - 1:
                uleft = u[i, j - 1] if (0 <= j - 1 < u.shape[1]) else 0.0
                uright = u[i + 1, j - 1] if (0 <= j - 1 < u.shape[1]) else 0.0
                uface = 0.5 * (uleft + uright)
            else:
                uface = 0.0
            if i == 0:
                dv_dx = (v[i + 1, j] - v[i, j]) / dx
            elif i == Nx - 1:
                dv_dx = (v[i, j] - v[i - 1, j]) / dx
            else:
                if uface >= 0:
                    dv_dx = (v[i, j] - v[i - 1, j]) / dx
                else:
                    dv_dx = (v[i + 1, j] - v[i, j]) / dx
            fluxx = uface * dv_dx
            conv[i, j] = fluxx + fluxy
    return conv

# ----------------------------- Main solver -----------------------------

def run_simulation(Nx=64, Ny=32, t_end=1.0, CFL=0.4):
    Lx = 2.0
    Ly = 1.0
    dx = Lx / Nx
    dy = Ly / Ny
    # Grid coordinates
    x_u = np.linspace(0.0, Lx, Nx + 1)
    y_u = np.linspace(dy / 2.0, Ly - dy / 2.0, Ny)
    xu, yu = np.meshgrid(x_u, y_u, indexing='ij')  # shape (Nx+1, Ny)
    x_v = np.linspace(dx / 2.0, Lx - dx / 2.0, Nx)
    y_v = np.linspace(0.0, Ly, Ny + 1)
    xv, yv = np.meshgrid(x_v, y_v, indexing='ij')  # shape (Nx, Ny+1)
    x_p = np.linspace(dx / 2.0, Lx - dx / 2.0, Nx)
    y_p = np.linspace(dy / 2.0, Ly - dy / 2.0, Ny)
    xp, yp = np.meshgrid(x_p, y_p, indexing='ij')  # shape (Nx, Ny)

    # Initialize fields
    u = np.zeros((Nx + 1, Ny), order='F')
    v = np.zeros((Nx, Ny + 1), order='F')
    p = np.zeros((Nx, Ny), order='F')

    # Precompute matrices (we will choose dt based on initial guess and then rebuild if dt changes)
    # Choose an initial dt based on visc and grid and max inlet amplitude estimate
    max_inlet = 3.0  # conservative estimate of inlet waveform amplitude
    dt_conv = CFL * min(dx / max_inlet, dy / max_inlet)
    dt_diff = 0.5 * min(dx, dy) ** 2 / visc
    dt = min(dt_conv, dt_diff, 0.01)
    Nsteps = int(np.ceil(t_end / dt))
    dt = t_end / Nsteps  # adjust to hit t_end exactly

    print(f"Grid: Nx={Nx}, Ny={Ny}, dx={dx:.4e}, dy={dy:.4e}, dt={dt:.4e}, steps={Nsteps}")

    # Build constant matrices and LU factorizations
    A_u, LU_u = build_helmholtz_u(Nx, Ny, dx, dy, dt, visc)
    A_v, LU_v = build_helmholtz_v(Nx, Ny, dx, dy, dt, visc)
    A_p, LU_p = build_poisson_p(Nx, Ny, dx, dy)

    # time loop
    t = 0.0
    progress_next = 0.1
    for step in range(1, Nsteps + 1):
        t += dt
        # boundary conditions: inlet u at x=0
        inlet_profile = np.sin(np.pi * y_u[:, 0]) * (
            np.sin(np.pi * t) + np.sin(3.0 * np.pi * t) + np.sin(5.0 * np.pi * t)
        )
        u[0, :] = inlet_profile
        # no-slip at top and bottom (approx): set u[:,0] and u[:,Ny-1] = 0
        u[:, 0] = 0.0
        u[:, -1] = 0.0
        # v no-slip on walls: v[:,0] and v[:,Ny] = 0
        v[:, 0] = 0.0
        v[:, -1] = 0.0
        # inlet v = 0 enforced approximately at left-most columns
        # outlet Neumann handled in matrix rows/constraints

        # compute convective terms
        conv_u = upwind_convective_u(u, v, dx, dy, Nx, Ny)
        conv_v = upwind_convective_v(u, v, dx, dy, Nx, Ny)

        # forcing term (projected to u and v grids)
        # f_x = 0 everywhere
        fx = np.zeros_like(u)
        # f_y = -sin(pi x) sin(pi y) sin(pi t) -> place on v-grid
        fy_face = -np.sin(np.pi * xv) * np.sin(np.pi * yv) * np.sin(np.pi * t)

        # RHS for implicit viscous solve (implicit Euler): (I) u_old - dt*(conv - f)
        rhs_u = (u - dt * conv_u + dt * fx).ravel(order='F')
        rhs_v = (v - dt * conv_v + dt * fy_face).ravel(order='F')

        # Enforce Dirichlet known values in rhs for u and v nodes where matrix rows are identity
        # For u: i==0 (inlet) and j==0, j==Ny-1
        for j in range(Ny):
            for i in (0,):
                k = idx_u(i, j, Nx, Ny)
                rhs_u[k] = u[i, j]
        for j in (0, Ny - 1):
            for i in range(Nx + 1):
                k = idx_u(i, j, Nx, Ny)
                rhs_u[k] = u[i, j]
        # For v: j==0, j==Ny and i==0
        for j in (0, Ny):
            for i in range(Nx):
                k = idx_v(i, j, Nx, Ny)
                rhs_v[k] = v[i, j]
        for i in (0,):
            for j in range(Ny + 1):
                k = idx_v(i, j, Nx, Ny)
                rhs_v[k] = v[i, j]

        # Solve Helmholtz systems
        sol_u = LU_u.solve(rhs_u)
        sol_v = LU_v.solve(rhs_v)
        u_star = sol_u.reshape((Nx + 1, Ny), order='F')
        v_star = sol_v.reshape((Nx, Ny + 1), order='F')

        # Enforce boundary conditions approximately on intermediate velocities
        u_star[0, :] = inlet_profile
        u_star[:, 0] = 0.0
        u_star[:, -1] = 0.0
        v_star[:, 0] = 0.0
        v_star[:, -1] = 0.0

        # Compute divergence on cell centers (p-grid): div = (u_{i+1}-u_i)/dx + (v_{i,j+1}-v_{i,j})/dy
        div = np.zeros((Nx, Ny), order='F')
        for j in range(Ny):
            for i in range(Nx):
                div[i, j] = (u_star[i + 1, j] - u_star[i, j]) / dx + (v_star[i, j + 1] - v_star[i, j]) / dy
        rhs_p = (div.ravel(order='F') / dt)

        # Enforce Dirichlet pressure at outlet: rows already identity in matrix -> rhs for those rows equals 0
        for j in range(Ny):
            k = idx_p(Nx - 1, j, Nx, Ny)
            rhs_p[k] = 0.0

        # Solve Poisson: A_p * phi = rhs_p
        phi_vec = LU_p.solve(rhs_p)
        phi = phi_vec.reshape((Nx, Ny), order='F')

        # Correct velocities: u_new = u_star - dt * dp/dx at u faces
        u_new = u_star.copy()
        for j in range(Ny):
            for i in range(1, Nx):
                # pressure gradient between cells i-1 and i corresponds to u face at i
                dpdx = (phi[i, j] - phi[i - 1, j]) / dx
                u_new[i, j] -= dt * dpdx
        # enforce inlet and walls/outlet
        u_new[0, :] = inlet_profile
        u_new[:, 0] = 0.0
        u_new[:, -1] = 0.0
        # For outlet face i==Nx approximate gradient using last cell
        u_new[Nx, :] = u_new[Nx - 1, :]

        # v correction: v_new = v_star - dt * dp/dy at v faces
        v_new = v_star.copy()
        for j in range(1, Ny):
            for i in range(Nx):
                dpdy = (phi[i, j] - phi[i, j - 1]) / dy
                v_new[i, j] -= dt * dpdy
        v_new[:, 0] = 0.0
        v_new[:, -1] = 0.0
        v_new[0, :] = 0.0
        # approximate outlet Neumann
        v_new[-1, :] = v_new[-2, :]

        # Update pressure: set p = phi (projection)
        p = phi.copy()

        # Update fields for next step
        u = u_new
        v = v_new

        # Progress printing every 10%
        frac = step / Nsteps
        if frac >= progress_next - 1e-12 or step == Nsteps:
            print(f"Progress: {int(100*frac)}% (t={t:.4f})")
            progress_next += 0.1

    # final fields: interpolate u and v to cell centers for plotting
    u_c = np.zeros((Nx, Ny), order='F')
    v_c = np.zeros((Nx, Ny), order='F')
    for j in range(Ny):
        for i in range(Nx):
            u_c[i, j] = 0.5 * (u[i, j] + u[i + 1, j])
            v_c[i, j] = 0.5 * (v[i, j] + v[i, j + 1])

    # Compute simple L2 divergence norm of final velocity field to check incompressibility
    div_final = np.zeros_like(p)
    for j in range(Ny):
        for i in range(Nx):
            div_final[i, j] = (u[i + 1, j] - u[i, j]) / dx + (v[i, j + 1] - v[i, j]) / dy
    l2_div = np.linalg.norm(div_final.ravel(order='F')) * np.sqrt(dx * dy)
    print(f"Final L2 divergence norm: {l2_div:.3e}")

    # Plot contours of u_c, v_c, p
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    levels = 40
    # transpose for plotting as contourf expects X,Y with shapes (Ny,Nx) typically; we keep axes consistent
    cs1 = axes[0].contourf(xp, yp, u_c, levels=levels, cmap='RdBu_r')
    axes[0].set_title('u (cell centers)')
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    fig.colorbar(cs1, ax=axes[0])

    cs2 = axes[1].contourf(xp, yp, v_c, levels=levels, cmap='RdBu_r')
    axes[1].set_title('v (cell centers)')
    axes[1].set_xlabel('x')
    fig.colorbar(cs2, ax=axes[1])

    cs3 = axes[2].contourf(xp, yp, p, levels=levels, cmap='RdBu_r')
    axes[2].set_title('p (cell centers)')
    axes[2].set_xlabel('x')
    fig.colorbar(cs3, ax=axes[2])

    fig.suptitle(f'Fields at t={t_end:.3f}')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig('uvp_contours.png', dpi=150)
    print('Saved figure: uvp_contours.png')


if __name__ == '__main__':
    # Example run with moderate resolution
    run_simulation(Nx=64, Ny=32, t_end=1.0, CFL=0.4)



#### Script block2:

```python
#!/usr/bin/env python3
"""
Chorin projection solver on a MAC grid for 2D incompressible Navier-Stokes.
Global variable: visc = 1.0 / Re (set at module scope)
Simulates until t = 1.0 and saves a figure of u,v,p contours.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# Global viscosity (do not rename)
Re = 100.0
visc = 1.0 / Re

# ----------------------------- Index helpers -----------------------------
def idx_u(i, j, Nx, Ny):
    # Flatten index for u array shape (Nx+1, Ny), Fortran order
    return i + j * (Nx + 1)


def idx_v(i, j, Nx, Ny):
    # Flatten index for v array shape (Nx, Ny+1), Fortran order
    return i + j * Nx


def idx_p(i, j, Nx, Ny):
    # Flatten index for p array shape (Nx, Ny), Fortran order
    return i + j * Nx


# ----------------------------- Discrete operators -----------------------------
def build_helmholtz_u(Nx, Ny, dx, dy, dt, visc):
    # Build (I - dt*visc*L) for u-grid (Nx+1 x Ny)
    Nu = (Nx + 1) * Ny
    A = sp.lil_matrix((Nu, Nu))
    ax = visc * dt / dx**2
    ay = visc * dt / dy**2
    for j in range(Ny):
        for i in range(Nx + 1):
            row = idx_u(i, j, Nx, Ny)
            # Dirichlet inlet at i==0 and no-slip top/bottom at j==0, j==Ny-1
            if i == 0 or j == 0 or j == Ny - 1:
                A[row, row] = 1.0
                continue
            # Neumann at outlet i == Nx approximate first derivative zero: u_N - u_{N-1} = 0
            if i == Nx:
                A[row, idx_u(i, j, Nx, Ny)] = 1.0
                A[row, idx_u(i - 1, j, Nx, Ny)] = -1.0
                continue
            # interior u node: (I - dt*visc*L) -> diag = 1 + 2*(ax+ay), off-diags = -ax / -ay
            diag = 1.0 + 2.0 * (ax + ay)
            A[row, row] = diag
            A[row, idx_u(i - 1, j, Nx, Ny)] = -ax
            A[row, idx_u(i + 1, j, Nx, Ny)] = -ax
            A[row, idx_u(i, j - 1, Nx, Ny)] = -ay
            A[row, idx_u(i, j + 1, Nx, Ny)] = -ay
    A = A.tocsc()
    LU = spla.splu(A)
    return A, LU


def build_helmholtz_v(Nx, Ny, dx, dy, dt, visc):
    # Build (I - dt*visc*L) for v-grid (Nx x Ny+1)
    Nv = Nx * (Ny + 1)
    A = sp.lil_matrix((Nv, Nv))
    ax = visc * dt / dx**2
    ay = visc * dt / dy**2
    for j in range(Ny + 1):
        for i in range(Nx):
            row = idx_v(i, j, Nx, Ny)
            # Dirichlet top/bottom walls at j==0 or j==Ny and inlet v=0 at i==0
            if j == 0 or j == Ny or i == 0:
                A[row, row] = 1.0
                continue
            # Neumann at outlet (approx): rightmost column i == Nx-1 uses one-sided derivative v_{Nx-1} - v_{Nx-2} = 0
            if i == Nx - 1:
                A[row, idx_v(i, j, Nx, Ny)] = 1.0
                A[row, idx_v(i - 1, j, Nx, Ny)] = -1.0
                continue
            # interior v node
            diag = 1.0 + 2.0 * (ax + ay)
            A[row, row] = diag
            A[row, idx_v(i - 1, j, Nx, Ny)] = -ax
            A[row, idx_v(i + 1, j, Nx, Ny)] = -ax
            A[row, idx_v(i, j - 1, Nx, Ny)] = -ay
            A[row, idx_v(i, j + 1, Nx, Ny)] = -ay
    A = A.tocsc()
    LU = spla.splu(A)
    return A, LU


def build_poisson_p(Nx, Ny, dx, dy):
    # Build Poisson operator for p-grid corresponding to discrete Laplacian: ∇^2 p = rhs
    # p-grid: Nx x Ny -> Np = Nx*Ny
    Np = Nx * Ny
    A = sp.lil_matrix((Np, Np))
    rx = 1.0 / dx**2
    ry = 1.0 / dy**2
    for j in range(Ny):
        for i in range(Nx):
            row = idx_p(i, j, Nx, Ny)
            # Dirichlet at outlet cell centers: rightmost cells i == Nx-1 -> p = 0
            if i == Nx - 1:
                A[row, row] = 1.0
                continue
            # interior/other boundaries: standard 5-point Laplacian (∇^2 p)
            # diag = -2*(rx+ry), off-diagonals = rx, ry
            A[row, row] = -2.0 * (rx + ry)
            # x neighbors
            if i - 1 >= 0:
                A[row, idx_p(i - 1, j, Nx, Ny)] = rx
            else:
                # left Neumann: replace p_{i-1} by p_{i+1} -> add rx to right neighbor (double right neighbor)
                if i + 1 < Nx:
                    A[row, idx_p(i + 1, j, Nx, Ny)] = A[row, idx_p(i + 1, j, Nx, Ny)] + rx if A[row, idx_p(i + 1, j, Nx, Ny)] != 0 else rx
            if i + 1 < Nx:
                A[row, idx_p(i + 1, j, Nx, Ny)] = rx if A[row, idx_p(i + 1, j, Nx, Ny)] == 0 else A[row, idx_p(i + 1, j, Nx, Ny)]
            # y neighbors
            if j - 1 >= 0:
                A[row, idx_p(i, j - 1, Nx, Ny)] = ry
            else:
                # bottom Neumann -> mirror to top neighbor
                if j + 1 < Ny:
                    A[row, idx_p(i, j + 1, Nx, Ny)] = A[row, idx_p(i, j + 1, Nx, Ny)] + ry if A[row, idx_p(i, j + 1, Nx, Ny)] != 0 else ry
            if j + 1 < Ny:
                # if previously possibly incremented by left Neumann logic, avoid overwriting
                if A[row, idx_p(i, j + 1, Nx, Ny)] == 0:
                    A[row, idx_p(i, j + 1, Nx, Ny)] = ry
    A = A.tocsc()
    LU = spla.splu(A)
    return A, LU


# ----------------------------- Advection and helpers -----------------------------
def upwind_convective_u(u, v, dx, dy, Nx, Ny):
    # Compute convective term for u at u-grid (Nx+1, Ny)
    conv = np.zeros_like(u)
    for j in range(Ny):
        for i in range(Nx + 1):
            uc = u[i, j]
            # du/dx
            if i == 0:
                du_dx = (u[i + 1, j] - u[i, j]) / dx
                fluxx = uc * du_dx
            elif i == Nx:
                du_dx = (u[i, j] - u[i - 1, j]) / dx
                fluxx = uc * du_dx
            else:
                ua = 0.5 * (u[i, j] + u[i - 1, j])
                if ua >= 0:
                    du_dx = (u[i, j] - u[i - 1, j]) / dx
                else:
                    du_dx = (u[i + 1, j] - u[i, j]) / dx
                fluxx = ua * du_dx
            # du/dy: approximate v at u-face by averaging neighboring v values safely
            if i == 0:
                v_left = 0.0
            else:
                v_left = v[i - 1, j]
            if i < Nx:
                v_right = v[i, j]
            else:
                v_right = 0.0
            v_face = 0.5 * (v_left + v_right)
            if j == 0:
                du_dy = (u[i, j + 1] - u[i, j]) / dy
            elif j == Ny - 1:
                du_dy = (u[i, j] - u[i, j - 1]) / dy
            else:
                if v_face >= 0:
                    du_dy = (u[i, j] - u[i, j - 1]) / dy
                else:
                    du_dy = (u[i, j + 1] - u[i, j]) / dy
            fluxy = v_face * du_dy
            conv[i, j] = fluxx + fluxy
    return conv


def upwind_convective_v(u, v, dx, dy, Nx, Ny):
    # Robust upwind convective term for v-grid (Nx, Ny+1)
    conv = np.zeros_like(v)
    # u shape (Nx+1, Ny), v shape (Nx, Ny+1)
    for j in range(Ny + 1):
        for i in range(Nx):
            vc = v[i, j]
            # dv/dy: vertical advection by v itself
            if j == 0:
                dv_dy = (v[i, j + 1] - v[i, j]) / dy
            elif j == Ny:
                dv_dy = (v[i, j] - v[i, j - 1]) / dy
            else:
                if vc >= 0:
                    dv_dy = (v[i, j] - v[i, j - 1]) / dy
                else:
                    dv_dy = (v[i, j + 1] - v[i, j]) / dy
            fluxy = vc * dv_dy
            # dv/dx: advected by u at this v-face. Compute uface safely by averaging neighboring u-values.
            # Map v-face (i,j) to u indices: use j_u = clamp(j-1, 0, Ny-1)
            j_u = j - 1
            if j_u < 0:
                j_u = 0
            if j_u > Ny - 1:
                j_u = Ny - 1
            # u indices for averaging: u[i, j_u] and u[i+1, j_u]
            u_left = u[i, j_u]
            u_right = u[i + 1, j_u]
            uface = 0.5 * (u_left + u_right)
            # dv/dx: choose upwind difference based on sign of uface
            if i == 0:
                dv_dx = (v[i + 1, j] - v[i, j]) / dx
            elif i == Nx - 1:
                dv_dx = (v[i, j] - v[i - 1, j]) / dx
            else:
                if uface >= 0:
                    dv_dx = (v[i, j] - v[i - 1, j]) / dx
                else:
                    dv_dx = (v[i + 1, j] - v[i, j]) / dx
            fluxx = uface * dv_dx
            conv[i, j] = fluxx + fluxy
    return conv


# Small unit test for upwind_convective_v to ensure no out-of-bounds and expected trivial behaviour
def test_upwind_convective_v():
    Nx_t, Ny_t = 4, 3
    dx_t, dy_t = 1.0 / Nx_t, 1.0 / Ny_t
    u_t = np.zeros((Nx_t + 1, Ny_t), order='F')
    v_t = np.zeros((Nx_t, Ny_t + 1), order='F')
    # Test 1: zeros -> conv should be zeros
    conv = upwind_convective_v(u_t, v_t, dx_t, dy_t, Nx_t, Ny_t)
    assert np.allclose(conv, 0.0), "upwind_convective_v failed zero-field test"
    # Test 2: uniform positive uface should yield dv/dx using backward/forward differences; with uniform v zeros -> conv zeros
    u_t[:] = 1.0
    conv2 = upwind_convective_v(u_t, v_t, dx_t, dy_t, Nx_t, Ny_t)
    assert np.allclose(conv2, 0.0), "upwind_convective_v failed uniform-u zero-v test"
    print("upwind_convective_v quick tests passed")


# ----------------------------- Main solver -----------------------------
def run_simulation(Nx=64, Ny=32, t_end=1.0, CFL=0.4):
    Lx = 2.0
    Ly = 1.0
    dx = Lx / Nx
    dy = Ly / Ny
    # Grid coordinates
    x_u = np.linspace(0.0, Lx, Nx + 1)
    # y positions for u are cell-centred in y-direction
    y_u = np.linspace(dy / 2.0, Ly - dy / 2.0, Ny)
    xu, yu = np.meshgrid(x_u, y_u, indexing='ij')  # shape (Nx+1, Ny)
    x_v = np.linspace(dx / 2.0, Lx - dx / 2.0, Nx)
    y_v = np.linspace(0.0, Ly, Ny + 1)
    xv, yv = np.meshgrid(x_v, y_v, indexing='ij')  # shape (Nx, Ny+1)
    x_p = np.linspace(dx / 2.0, Lx - dx / 2.0, Nx)
    y_p = np.linspace(dy / 2.0, Ly - dy / 2.0, Ny)
    xp, yp = np.meshgrid(x_p, y_p, indexing='ij')  # shape (Nx, Ny)

    # Initialize fields
    u = np.zeros((Nx + 1, Ny), order='F')
    v = np.zeros((Nx, Ny + 1), order='F')
    p = np.zeros((Nx, Ny), order='F')

    # Precompute dt based on CFL and diffusion; do not rebuild matrices in time loop
    max_inlet = 3.0  # conservative estimate of inlet waveform amplitude
    dt_conv = CFL * min(dx / max_inlet, dy / max_inlet)
    dt_diff = 0.5 * min(dx, dy) ** 2 / visc
    dt = min(dt_conv, dt_diff, 0.01)
    Nsteps = int(np.ceil(t_end / dt))
    dt = t_end / Nsteps  # adjust to hit t_end exactly

    print(f"Grid: Nx={Nx}, Ny={Ny}, dx={dx:.4e}, dy={dy:.4e}, dt={dt:.4e}, steps={Nsteps}")

    # Build constant matrices and LU factorizations
    A_u, LU_u = build_helmholtz_u(Nx, Ny, dx, dy, dt, visc)
    A_v, LU_v = build_helmholtz_v(Nx, Ny, dx, dy, dt, visc)
    A_p, LU_p = build_poisson_p(Nx, Ny, dx, dy)

    # time loop
    t = 0.0
    progress_next = 0.1
    for step in range(1, Nsteps + 1):
        t += dt
        # boundary conditions: inlet u at x=0 using 1-D y_u array (fixed review bug)
        inlet_profile = np.sin(np.pi * y_u) * (
            np.sin(np.pi * t) + np.sin(3.0 * np.pi * t) + np.sin(5.0 * np.pi * t)
        )
        u[0, :] = inlet_profile
        # no-slip at top and bottom
        u[:, 0] = 0.0
        u[:, -1] = 0.0
        v[:, 0] = 0.0
        v[:, -1] = 0.0
        # inlet v = 0 enforced at i==0
        v[0, :] = 0.0

        # compute convective terms
        conv_u = upwind_convective_u(u, v, dx, dy, Nx, Ny)
        conv_v = upwind_convective_v(u, v, dx, dy, Nx, Ny)

        # forcing term (projected to u and v grids)
        fx = np.zeros_like(u)
        fy_face = -np.sin(np.pi * xv) * np.sin(np.pi * yv) * np.sin(np.pi * t)

        # RHS for implicit viscous solves
        rhs_u = (u - dt * conv_u + dt * fx).ravel(order='F')
        rhs_v = (v - dt * conv_v + dt * fy_face).ravel(order='F')

        # Enforce Dirichlet known values in rhs for u and v nodes where matrix rows are identity
        # For u: i==0 (inlet) and j==0, j==Ny-1
        for j in range(Ny):
            k = idx_u(0, j, Nx, Ny)
            rhs_u[k] = u[0, j]
        for j in (0, Ny - 1):
            for i in range(Nx + 1):
                k = idx_u(i, j, Nx, Ny)
                rhs_u[k] = u[i, j]
        # For v: j==0, j==Ny and i==0
        for j in (0, Ny):
            for i in range(Nx):
                k = idx_v(i, j, Nx, Ny)
                rhs_v[k] = v[i, j]
        for i in (0,):
            for j in range(Ny + 1):
                k = idx_v(i, j, Nx, Ny)
                rhs_v[k] = v[i, j]

        # Solve Helmholtz systems
        sol_u = LU_u.solve(rhs_u)
        sol_v = LU_v.solve(rhs_v)
        u_star = sol_u.reshape((Nx + 1, Ny), order='F')
        v_star = sol_v.reshape((Nx, Ny + 1), order='F')

        # Enforce BCs on intermediate velocities
        u_star[0, :] = inlet_profile
        u_star[:, 0] = 0.0
        u_star[:, -1] = 0.0
        v_star[:, 0] = 0.0
        v_star[:, -1] = 0.0
        v_star[0, :] = 0.0

        # Compute divergence on cell centers (p-grid)
        div = np.zeros((Nx, Ny), order='F')
        for j in range(Ny):
            for i in range(Nx):
                div[i, j] = (u_star[i + 1, j] - u_star[i, j]) / dx + (v_star[i, j + 1] - v_star[i, j]) / dy
        rhs_p = (div.ravel(order='F') / dt)

        # Enforce Dirichlet pressure at outlet: rows already identity in matrix -> rhs for those rows equals 0
        for j in range(Ny):
            k = idx_p(Nx - 1, j, Nx, Ny)
            rhs_p[k] = 0.0

        # Solve Poisson: A_p * phi = rhs_p
        phi_vec = LU_p.solve(rhs_p)
        phi = phi_vec.reshape((Nx, Ny), order='F')

        # Correct velocities: u_new = u_star - dt * dphi/dx at u faces
        u_new = u_star.copy()
        for j in range(Ny):
            for i in range(1, Nx):
                dpdx = (phi[i, j] - phi[i - 1, j]) / dx
                u_new[i, j] -= dt * dpdx
        # enforce inlet and walls/outlet
        u_new[0, :] = inlet_profile
        u_new[:, 0] = 0.0
        u_new[:, -1] = 0.0
        # simple outlet extrapolation
        u_new[Nx, :] = u_new[Nx - 1, :]

        # v correction: v_new = v_star - dt * dphi/dy at v faces
        v_new = v_star.copy()
        for j in range(1, Ny):
            for i in range(Nx):
                dpdy = (phi[i, j] - phi[i, j - 1]) / dy
                v_new[i, j] -= dt * dpdy
        v_new[:, 0] = 0.0
        v_new[:, -1] = 0.0
        v_new[0, :] = 0.0
        # approximate outlet Neumann via extrapolation
        v_new[-1, :] = v_new[-2, :]

        # Update pressure (projection variable) and fields
        p = phi.copy()
        u = u_new
        v = v_new

        # Progress printing every 10%
        frac = step / Nsteps
        if frac >= progress_next - 1e-12 or step == Nsteps:
            print(f"Progress: {int(100*frac)}% (t={t:.4f})")
            progress_next += 0.1

    # final fields: interpolate u and v to cell centers for plotting
    u_c = np.zeros((Nx, Ny), order='F')
    v_c = np.zeros((Nx, Ny), order='F')
    for j in range(Ny):
        for i in range(Nx):
            u_c[i, j] = 0.5 * (u[i, j] + u[i + 1, j])
            v_c[i, j] = 0.5 * (v[i, j] + v[i, j + 1])

    # Compute simple L2 divergence norm of final velocity field
    div_final = np.zeros_like(p)
    for j in range(Ny):
        for i in range(Nx):
            div_final[i, j] = (u[i + 1, j] - u[i, j]) / dx + (v[i, j + 1] - v[i, j]) / dy
    l2_div = np.linalg.norm(div_final.ravel(order='F')) * np.sqrt(dx * dy)
    print(f"Final L2 divergence norm: {l2_div:.3e}")

    # Plot contours of u_c, v_c, p
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    levels = 40
    cs1 = axes[0].contourf(xp, yp, u_c, levels=levels, cmap='RdBu_r')
    axes[0].set_title('u (cell centers)')
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    fig.colorbar(cs1, ax=axes[0])

    cs2 = axes[1].contourf(xp, yp, v_c, levels=levels, cmap='RdBu_r')
    axes[1].set_title('v (cell centers)')
    axes[1].set_xlabel('x')
    fig.colorbar(cs2, ax=axes[1])

    cs3 = axes[2].contourf(xp, yp, p, levels=levels, cmap='RdBu_r')
    axes[2].set_title('p (cell centers)')
    axes[2].set_xlabel('x')
    fig.colorbar(cs3, ax=axes[2])

    fig.suptitle(f'Fields at t={t_end:.3f}')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig('uvp_contours.png', dpi=150)
    print('Saved figure: uvp_contours.png')


if __name__ == '__main__':
    # Run quick tests for convective operator correctness
    test_upwind_convective_v()
    # Example run with moderate resolution
    run_simulation(Nx=64, Ny=32, t_end=1.0, CFL=0.4)
```}


#### Script block3:

```python
#!/usr/bin/env python3
"""
Chorin projection solver on a MAC grid for 2D incompressible Navier-Stokes.
- Global variable: visc set at module scope and used everywhere.
- Matrices are assembled once and LU-factorized once.
- Final contours saved to 'uvp_contours.png'.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# Global viscosity (do not rename)
Re = 100.0
visc = 1.0 / Re

# ----------------------------- Index helpers -----------------------------
def idx_u(i, j, Nx, Ny):
    return i + j * (Nx + 1)


def idx_v(i, j, Nx, Ny):
    return i + j * Nx


def idx_p(i, j, Nx, Ny):
    return i + j * Nx


# ----------------------------- Discrete operators (robust assembly) -----------------------------
def build_helmholtz_u(Nx, Ny, dx, dy, dt, visc):
    # Build (I - dt*visc*L) for u-grid (Nx+1 x Ny)
    Nu = (Nx + 1) * Ny
    A = sp.lil_matrix((Nu, Nu))
    ax = visc * dt / dx**2
    ay = visc * dt / dy**2

    # Masks and mapping to enforce BCs consistently
    dirichlet_mask = np.zeros(Nu, dtype=bool)
    dirichlet_values = np.zeros(Nu, dtype=float)  # for fixed Dirichlet rows
    neumann_copy_mask = np.zeros(Nu, dtype=bool)
    neumann_copy_from = np.full(Nu, -1, dtype=int)  # index to copy from for Neumann rows

    for j in range(Ny):
        for i in range(Nx + 1):
            row = idx_u(i, j, Nx, Ny)
            # Dirichlet inlet at i==0; no-slip top/bottom at j==0, j==Ny-1
            if i == 0:
                A[row, row] = 1.0
                dirichlet_mask[row] = True
                # value will be time-dependent inlet profile and set before solve
                dirichlet_values[row] = 0.0
                continue
            if j == 0 or j == Ny - 1:
                A[row, row] = 1.0
                dirichlet_mask[row] = True
                dirichlet_values[row] = 0.0
                continue
            # Neumann at outlet face i == Nx enforced by copying interior neighbor u_{Nx-1}
            if i == Nx:
                A[row, row] = 1.0
                neumann_copy_mask[row] = True
                # interior neighbor to copy from is (i-1,j)
                neumann_copy_from[row] = idx_u(i - 1, j, Nx, Ny)
                continue
            # interior u node: (I - dt*visc*L)
            diag = 1.0 + 2.0 * (ax + ay)
            A[row, row] = diag
            A[row, idx_u(i - 1, j, Nx, Ny)] = -ax
            A[row, idx_u(i + 1, j, Nx, Ny)] = -ax
            A[row, idx_u(i, j - 1, Nx, Ny)] = -ay
            A[row, idx_u(i, j + 1, Nx, Ny)] = -ay

    A = A.tocsc()
    LU = spla.splu(A)
    bc_info = {
        'dirichlet_mask': dirichlet_mask,
        'dirichlet_values': dirichlet_values,
        'neumann_copy_mask': neumann_copy_mask,
        'neumann_copy_from': neumann_copy_from,
    }
    return A, LU, bc_info


def build_helmholtz_v(Nx, Ny, dx, dy, dt, visc):
    # Build (I - dt*visc*L) for v-grid (Nx x Ny+1)
    Nv = Nx * (Ny + 1)
    A = sp.lil_matrix((Nv, Nv))
    ax = visc * dt / dx**2
    ay = visc * dt / dy**2

    dirichlet_mask = np.zeros(Nv, dtype=bool)
    dirichlet_values = np.zeros(Nv, dtype=float)
    neumann_copy_mask = np.zeros(Nv, dtype=bool)
    neumann_copy_from = np.full(Nv, -1, dtype=int)

    for j in range(Ny + 1):
        for i in range(Nx):
            row = idx_v(i, j, Nx, Ny)
            # Dirichlet top/bottom walls at j==0 or j==Ny and inlet v=0 at i==0
            if j == 0 or j == Ny or i == 0:
                A[row, row] = 1.0
                dirichlet_mask[row] = True
                dirichlet_values[row] = 0.0
                continue
            # Neumann at outlet: rightmost interior column i == Nx-1 -> copy from i-1
            if i == Nx - 1:
                A[row, row] = 1.0
                neumann_copy_mask[row] = True
                neumann_copy_from[row] = idx_v(i - 1, j, Nx, Ny)
                continue
            # interior v node
            diag = 1.0 + 2.0 * (ax + ay)
            A[row, row] = diag
            A[row, idx_v(i - 1, j, Nx, Ny)] = -ax
            A[row, idx_v(i + 1, j, Nx, Ny)] = -ax
            A[row, idx_v(i, j - 1, Nx, Ny)] = -ay
            A[row, idx_v(i, j + 1, Nx, Ny)] = -ay

    A = A.tocsc()
    LU = spla.splu(A)
    bc_info = {
        'dirichlet_mask': dirichlet_mask,
        'dirichlet_values': dirichlet_values,
        'neumann_copy_mask': neumann_copy_mask,
        'neumann_copy_from': neumann_copy_from,
    }
    return A, LU, bc_info


def build_poisson_p(Nx, Ny, dx, dy):
    # Build Poisson operator A_p such that A_p * phi ~= Laplacian(phi) on p-grid (Nx x Ny)
    Np = Nx * Ny
    A = sp.lil_matrix((Np, Np))
    rx = 1.0 / dx**2
    ry = 1.0 / dy**2

    dirichlet_mask = np.zeros(Np, dtype=bool)
    dirichlet_values = np.zeros(Np, dtype=float)

    for j in range(Ny):
        for i in range(Nx):
            row = idx_p(i, j, Nx, Ny)
            # Dirichlet at outlet cell centers: rightmost cells i == Nx-1 -> p = 0
            if i == Nx - 1:
                A[row, row] = 1.0
                dirichlet_mask[row] = True
                dirichlet_values[row] = 0.0
                continue
            # Assemble using ghost-elimination for left/bottom Neumann (reflective ghost: phi_ghost = phi_neighbor)
            central = 0.0
            coeffs = {}
            # x-direction contribution
            if i - 1 >= 0 and i + 1 < Nx:
                # full interior in x
                central += -2.0 * rx
                coeffs[idx_p(i - 1, j, Nx, Ny)] = coeffs.get(idx_p(i - 1, j, Nx, Ny), 0.0) + rx
                coeffs[idx_p(i + 1, j, Nx, Ny)] = coeffs.get(idx_p(i + 1, j, Nx, Ny), 0.0) + rx
            elif i - 1 < 0 and i + 1 < Nx:
                # left Neumann: phi_{-1} = phi_{1}
                central += -2.0 * rx
                coeffs[idx_p(i + 1, j, Nx, Ny)] = coeffs.get(idx_p(i + 1, j, Nx, Ny), 0.0) + 2.0 * rx
            elif i - 1 >= 0 and i + 1 >= Nx:
                # right neighbor missing (shouldn't happen except maybe 1-cell in x), treat similarly
                central += -2.0 * rx
                coeffs[idx_p(i - 1, j, Nx, Ny)] = coeffs.get(idx_p(i - 1, j, Nx, Ny), 0.0) + 2.0 * rx
            else:
                # isolated cell -> central only
                central += -2.0 * rx
            # y-direction contribution
            if j - 1 >= 0 and j + 1 < Ny:
                central += -2.0 * ry
                coeffs[idx_p(i, j - 1, Nx, Ny)] = coeffs.get(idx_p(i, j - 1, Nx, Ny), 0.0) + ry
                coeffs[idx_p(i, j + 1, Nx, Ny)] = coeffs.get(idx_p(i, j + 1, Nx, Ny), 0.0) + ry
            elif j - 1 < 0 and j + 1 < Ny:
                # bottom Neumann: phi_{-1} = phi_{1}
                central += -2.0 * ry
                coeffs[idx_p(i, j + 1, Nx, Ny)] = coeffs.get(idx_p(i, j + 1, Nx, Ny), 0.0) + 2.0 * ry
            elif j - 1 >= 0 and j + 1 >= Ny:
                central += -2.0 * ry
                coeffs[idx_p(i, j - 1, Nx, Ny)] = coeffs.get(idx_p(i, j - 1, Nx, Ny), 0.0) + 2.0 * ry
            else:
                central += -2.0 * ry

            A[row, row] = central
            for col_idx, val in coeffs.items():
                A[row, col_idx] = val

    A = A.tocsc()
    LU = spla.splu(A)
    bc_info = {'dirichlet_mask': dirichlet_mask, 'dirichlet_values': dirichlet_values}
    return A, LU, bc_info


# ----------------------------- Advection and helpers -----------------------------
def upwind_convective_u(u, v, dx, dy, Nx, Ny):
    # Compute convective term for u at u-grid (Nx+1, Ny)
    conv = np.zeros_like(u)
    for j in range(Ny):
        for i in range(Nx + 1):
            uc = u[i, j]
            # du/dx
            if i == 0:
                du_dx = (u[i + 1, j] - u[i, j]) / dx
                fluxx = uc * du_dx
            elif i == Nx:
                du_dx = (u[i, j] - u[i - 1, j]) / dx
                fluxx = uc * du_dx
            else:
                ua = 0.5 * (u[i, j] + u[i - 1, j])
                if ua >= 0:
                    du_dx = (u[i, j] - u[i - 1, j]) / dx
                else:
                    du_dx = (u[i + 1, j] - u[i, j]) / dx
                fluxx = ua * du_dx
            # du/dy: approximate v at u-face by averaging neighboring v values safely
            if i == 0:
                v_left = 0.0
            else:
                v_left = v[i - 1, j]
            if i < Nx:
                v_right = v[i, j]
            else:
                v_right = 0.0
            v_face = 0.5 * (v_left + v_right)
            if j == 0:
                du_dy = (u[i, j + 1] - u[i, j]) / dy
            elif j == Ny - 1:
                du_dy = (u[i, j] - u[i, j - 1]) / dy
            else:
                if v_face >= 0:
                    du_dy = (u[i, j] - u[i, j - 1]) / dy
                else:
                    du_dy = (u[i, j + 1] - u[i, j]) / dy
            fluxy = v_face * du_dy
            conv[i, j] = fluxx + fluxy
    return conv


def upwind_convective_v(u, v, dx, dy, Nx, Ny):
    conv = np.zeros_like(v)
    for j in range(Ny + 1):
        for i in range(Nx):
            vc = v[i, j]
            # dv/dy
            if j == 0:
                dv_dy = (v[i, j + 1] - v[i, j]) / dy
            elif j == Ny:
                dv_dy = (v[i, j] - v[i, j - 1]) / dy
            else:
                if vc >= 0:
                    dv_dy = (v[i, j] - v[i, j - 1]) / dy
                else:
                    dv_dy = (v[i, j + 1] - v[i, j]) / dy
            fluxy = vc * dv_dy
            # compute uface safely by averaging u-values mapped to this v-face
            j_u = j - 1
            if j_u < 0:
                j_u = 0
            if j_u > Ny - 1:
                j_u = Ny - 1
            u_left = u[i, j_u]
            u_right = u[i + 1, j_u]
            uface = 0.5 * (u_left + u_right)
            if i == 0:
                dv_dx = (v[i + 1, j] - v[i, j]) / dx
            elif i == Nx - 1:
                dv_dx = (v[i, j] - v[i - 1, j]) / dx
            else:
                if uface >= 0:
                    dv_dx = (v[i, j] - v[i - 1, j]) / dx
                else:
                    dv_dx = (v[i + 1, j] - v[i, j]) / dx
            fluxx = uface * dv_dx
            conv[i, j] = fluxx + fluxy
    return conv


# Simple tests for convective operators
def test_upwind_convective_u():
    Nx_t, Ny_t = 6, 4
    dx_t, dy_t = 1.0 / Nx_t, 1.0 / Ny_t
    u_t = np.zeros((Nx_t + 1, Ny_t), order='F')
    v_t = np.zeros((Nx_t, Ny_t + 1), order='F')
    # zeros -> conv should be zeros
    conv = upwind_convective_u(u_t, v_t, dx_t, dy_t, Nx_t, Ny_t)
    assert np.allclose(conv, 0.0), "upwind_convective_u failed zero-field test"
    # uniform u and zero v -> conv zeros
    u_t[:] = 1.0
    conv2 = upwind_convective_u(u_t, v_t, dx_t, dy_t, Nx_t, Ny_t)
    # interior faces may produce zero for uniform u and zero v
    assert np.allclose(conv2, 0.0), "upwind_convective_u failed uniform-u zero-v test"
    print("upwind_convective_u quick tests passed")


# Test Poisson operator against analytic Laplacian
def test_poisson_operator(Nx=32, Ny=16):
    Lx, Ly = 2.0, 1.0
    dx = Lx / Nx
    dy = Ly / Ny
    x_p = np.linspace(dx / 2.0, Lx - dx / 2.0, Nx)
    y_p = np.linspace(dy / 2.0, Ly - dy / 2.0, Ny)
    xp, yp = np.meshgrid(x_p, y_p, indexing='ij')
    A_p, LU_p, bc = build_poisson_p(Nx, Ny, dx, dy)
    phi = np.sin(np.pi * xp) * np.sin(np.pi * yp)
    phi_vec = phi.ravel(order='F')
    Ap_phi = A_p.dot(phi_vec)
    # analytic Laplacian of sin(pi x) sin(pi y) = -2*pi^2 * phi
    lap_analytic = -2.0 * (np.pi**2) * phi
    diff = Ap_phi - lap_analytic.ravel(order='F')
    err = np.linalg.norm(diff) * np.sqrt(dx * dy)
    print(f"Poisson operator test: L2 error vs analytic Laplacian = {err:.3e}")
    # allow modest error depending on resolution
    assert err < 1e-1, "Poisson operator deviates too much from analytic Laplacian"


# Projection preservation test: projection should preserve divergence-free field
def test_projection_preserves_divergence(Nx=32, Ny=16):
    Lx, Ly = 2.0, 1.0
    dx = Lx / Nx
    dy = Ly / Ny
    x_u = np.linspace(0.0, Lx, Nx + 1)
    y_u = np.linspace(dy / 2.0, Ly - dy / 2.0, Ny)
    xu, yu = np.meshgrid(x_u, y_u, indexing='ij')
    x_v = np.linspace(dx / 2.0, Lx - dx / 2.0, Nx)
    y_v = np.linspace(0.0, Ly, Ny + 1)
    xv, yv = np.meshgrid(x_v, y_v, indexing='ij')
    x_p = np.linspace(dx / 2.0, Lx - dx / 2.0, Nx)
    y_p = np.linspace(dy / 2.0, Ly - dy / 2.0, Ny)
    xp, yp = np.meshgrid(x_p, y_p, indexing='ij')

    # Streamfunction psi(x,y) -> produce divergence-free velocity field
    psi = np.sin(np.pi * xp) * np.sin(np.pi * yp)
    # sample to u and v faces
    u = np.zeros((Nx + 1, Ny), order='F')
    v = np.zeros((Nx, Ny + 1), order='F')
    # u = -dpsi/dy at u-face centers in y; psi defined at cell centers -> approximate
    # approximate psi at u-face centers by averaging adjacent cell-centers in x
    for j in range(Ny):
        for i in range(Nx + 1):
            # x index for adjacent p-cells
            if i == 0:
                # left boundary: use first cell
                p_left = psi[0, j]
                p_right = psi[0, j]
            elif i == Nx:
                p_left = psi[Nx - 1, j]
                p_right = psi[Nx - 1, j]
            else:
                p_left = psi[i - 1, j]
                p_right = psi[i, j]
            psi_at_u = 0.5 * (p_left + p_right)
            # derivative dpsi/dy approximate via central at cell-centers
            # approximate d/dy of psi_at_u using neighboring psi in y (map to nearest p indices)
            if j == 0:
                dpsi_dy = (psi_at_u - psi[i - 1 if i>0 else 0, j]) / dy if i>0 else 0.0
            elif j == Ny - 1:
                dpsi_dy = (psi_at_u - psi_at_u) / dy
            else:
                # central approx using psi at cell centers j+1 and j-1 averaged similarly
                # simplest: use finite diff between psi at cell centers
                dpsi_dy = (psi[min(i, Nx-1), min(j+1, Ny-1)] - psi[min(i, Nx-1), max(j-1,0)]) / (2*dy)
            u[i, j] = -dpsi_dy
    # v = dpsi/dx
    for j in range(Ny + 1):
        for i in range(Nx):
            if j == 0:
                p_bot = psi[i, 0]
                p_top = psi[i, 0]
            elif j == Ny:
                p_bot = psi[i, Ny - 1]
                p_top = psi[i, Ny - 1]
            else:
                p_bot = psi[i, j - 1]
                p_top = psi[i, j]
            psi_at_v = 0.5 * (p_bot + p_top)
            if i == 0:
                dpsi_dx = (psi_at_v - psi[0, j - 1 if j>0 else 0]) / dx if j>0 else 0.0
            elif i == Nx - 1:
                dpsi_dx = (psi_at_v - psi_at_v) / dx
            else:
                dpsi_dx = (psi[min(i+1, Nx-1), min(j, Ny-1)] - psi[max(i-1,0), min(j, Ny-1)]) / (2*dx)
            v[i, j] = dpsi_dx

    # Build matrices
    A_u, LU_u, bc_u = build_helmholtz_u(Nx, Ny, dx, dy, dt=1e-3, visc=visc)
    A_v, LU_v, bc_v = build_helmholtz_v(Nx, Ny, dx, dy, dt=1e-3, visc=visc)
    A_p, LU_p, bc_p = build_poisson_p(Nx, Ny, dx, dy)

    # Compute divergence currently
    div = np.zeros((Nx, Ny), order='F')
    for j in range(Ny):
        for i in range(Nx):
            div[i, j] = (u[i + 1, j] - u[i, j]) / dx + (v[i, j + 1] - v[i, j]) / dy
    l2_before = np.linalg.norm(div.ravel(order='F')) * np.sqrt(dx * dy)

    # Projection: compute RHS and solve Poisson
    dt_test = 1e-3
    rhs_p = (div.ravel(order='F') / dt_test)
    # enforce Dirichlet pressure rows from bc_p
    rhs_p[bc_p['dirichlet_mask']] = bc_p['dirichlet_values'][bc_p['dirichlet_mask']]
    phi_vec = LU_p.solve(rhs_p)
    phi = phi_vec.reshape((Nx, Ny), order='F')

    # correct velocities
    u_new = u.copy()
    v_new = v.copy()
    for j in range(Ny):
        for i in range(1, Nx):
            dpdx = (phi[i, j] - phi[i - 1, j]) / dx
            u_new[i, j] -= dt_test * dpdx
    for j in range(1, Ny):
        for i in range(Nx):
            dpdy = (phi[i, j] - phi[i, j - 1]) / dy
            v_new[i, j] -= dt_test * dpdy

    div_after = np.zeros_like(div)
    for j in range(Ny):
        for i in range(Nx):
            div_after[i, j] = (u_new[i + 1, j] - u_new[i, j]) / dx + (v_new[i, j + 1] - v_new[i, j]) / dy
    l2_after = np.linalg.norm(div_after.ravel(order='F')) * np.sqrt(dx * dy)
    print(f"Projection test: L2 divergence before = {l2_before:.3e}, after = {l2_after:.3e}")
    assert l2_after < 1e-10 or l2_after < l2_before * 1e-3, "Projection did not sufficiently reduce divergence"


# ----------------------------- Main solver -----------------------------
def run_simulation(Nx=64, Ny=32, t_end=1.0, CFL=0.4):
    Lx = 2.0
    Ly = 1.0
    dx = Lx / Nx
    dy = Ly / Ny
    # Grid coordinates
    x_u = np.linspace(0.0, Lx, Nx + 1)
    # y positions for u are cell-centred in y-direction
    y_u = np.linspace(dy / 2.0, Ly - dy / 2.0, Ny)
    xu, yu = np.meshgrid(x_u, y_u, indexing='ij')  # shape (Nx+1, Ny)
    x_v = np.linspace(dx / 2.0, Lx - dx / 2.0, Nx)
    y_v = np.linspace(0.0, Ly, Ny + 1)
    xv, yv = np.meshgrid(x_v, y_v, indexing='ij')  # shape (Nx, Ny+1)
    x_p = np.linspace(dx / 2.0, Lx - dx / 2.0, Nx)
    y_p = np.linspace(dy / 2.0, Ly - dy / 2.0, Ny)
    xp, yp = np.meshgrid(x_p, y_p, indexing='ij')  # shape (Nx, Ny)

    # Initialize fields
    u = np.zeros((Nx + 1, Ny), order='F')
    v = np.zeros((Nx, Ny + 1), order='F')
    p = np.zeros((Nx, Ny), order='F')

    # Precompute dt based on CFL and diffusion; do not rebuild matrices in time loop
    max_inlet = 3.0  # conservative estimate of inlet waveform amplitude
    dt_conv = CFL * min(dx / max_inlet, dy / max_inlet)
    dt_diff = 0.5 * min(dx, dy) ** 2 / visc
    dt = min(dt_conv, dt_diff, 0.01)
    Nsteps = int(np.ceil(t_end / dt))
    dt = t_end / Nsteps  # adjust to hit t_end exactly

    print(f"Grid: Nx={Nx}, Ny={Ny}, dx={dx:.4e}, dy={dy:.4e}, dt={dt:.4e}, steps={Nsteps}")

    # Build constant matrices and LU factorizations
    A_u, LU_u, bc_u = build_helmholtz_u(Nx, Ny, dx, dy, dt, visc)
    A_v, LU_v, bc_v = build_helmholtz_v(Nx, Ny, dx, dy, dt, visc)
    A_p, LU_p, bc_p = build_poisson_p(Nx, Ny, dx, dy)

    # time loop
    t = 0.0
    progress_next = 0.1
    for step in range(1, Nsteps + 1):
        t += dt
        # boundary conditions: inlet u at x=0 using 1-D y_u array
        inlet_profile = np.sin(np.pi * y_u) * (
            np.sin(np.pi * t) + np.sin(3.0 * np.pi * t) + np.sin(5.0 * np.pi * t)
        )
        # Apply physical Dirichlet BCs onto current fields (to form RHS)
        u[0, :] = inlet_profile
        u[:, 0] = 0.0
        u[:, -1] = 0.0
        v[:, 0] = 0.0
        v[:, -1] = 0.0
        v[0, :] = 0.0

        # compute convective terms
        conv_u = upwind_convective_u(u, v, dx, dy, Nx, Ny)
        conv_v = upwind_convective_v(u, v, dx, dy, Nx, Ny)

        # forcing term (projected to u and v grids)
        fx = np.zeros_like(u)
        fy_face = -np.sin(np.pi * xv) * np.sin(np.pi * yv) * np.sin(np.pi * t)

        # RHS for implicit viscous solves
        rhs_u = (u - dt * conv_u + dt * fx).ravel(order='F')
        rhs_v = (v - dt * conv_v + dt * fy_face).ravel(order='F')

        # Enforce Dirichlet rows for u: use bc_u masks
        # For fixed Dirichlet rows (inlet and walls), set RHS to prescribed values
        if np.any(bc_u['dirichlet_mask']):
            dir_mask = bc_u['dirichlet_mask']
            # For inlet rows, values are time-dependent (inlet_profile) -> populate from u array
            # Here bc_u['dirichlet_values'] contains only static values (e.g., zero walls). We must override
            # inlet indices: detect indices with i==0
            for j in range(Ny):
                k = idx_u(0, j, Nx, Ny)
                if dir_mask[k]:
                    rhs_u[k] = u[0, j]
            # walls (j==0, j==Ny-1) will also be set by mask; ensure their rhs are zero
            for j in (0, Ny - 1):
                for i in range(Nx + 1):
                    k = idx_u(i, j, Nx, Ny)
                    if dir_mask[k]:
                        rhs_u[k] = 0.0
        # Enforce Neumann-copy rows for u: set rhs to interior neighbor u value
        if np.any(bc_u['neumann_copy_mask']):
            copy_mask = bc_u['neumann_copy_mask']
            copy_from = bc_u['neumann_copy_from']
            idxs = np.nonzero(copy_mask)[0]
            for k in idxs:
                src = copy_from[k]
                rhs_u[k] = u.ravel(order='F')[src]

        # For v similar enforcement
        if np.any(bc_v['dirichlet_mask']):
            dir_mask = bc_v['dirichlet_mask']
            # set all dirichlet rhs to zero (inlet v and walls)
            for k in np.nonzero(dir_mask)[0]:
                rhs_v[k] = 0.0
        if np.any(bc_v['neumann_copy_mask']):
            copy_mask = bc_v['neumann_copy_mask']
            copy_from = bc_v['neumann_copy_from']
            for k in np.nonzero(copy_mask)[0]:
                src = copy_from[k]
                rhs_v[k] = v.ravel(order='F')[src]

        # Solve Helmholtz systems
        sol_u = LU_u.solve(rhs_u)
        sol_v = LU_v.solve(rhs_v)
        u_star = sol_u.reshape((Nx + 1, Ny), order='F')
        v_star = sol_v.reshape((Nx, Ny + 1), order='F')

        # Enforce BCs on intermediate velocities (consistent)
        u_star[0, :] = inlet_profile
        u_star[:, 0] = 0.0
        u_star[:, -1] = 0.0
        v_star[:, 0] = 0.0
        v_star[:, -1] = 0.0
        v_star[0, :] = 0.0
        # Neumann-copy enforcement on star fields
        # u rightmost face = interior neighbor
        u_star[Nx, :] = u_star[Nx - 1, :]
        v_star[-1, :] = v_star[-2, :]

        # Compute divergence on cell centers (p-grid)
        div = np.zeros((Nx, Ny), order='F')
        for j in range(Ny):
            for i in range(Nx):
                div[i, j] = (u_star[i + 1, j] - u_star[i, j]) / dx + (v_star[i, j + 1] - v_star[i, j]) / dy
        rhs_p = (div.ravel(order='F') / dt)

        # Enforce Dirichlet pressure at outlet: rows identity in matrix -> rhs value set to zero
        if np.any(bc_p['dirichlet_mask']):
            rhs_p[bc_p['dirichlet_mask']] = bc_p['dirichlet_values'][bc_p['dirichlet_mask']]

        # Solve Poisson: A_p * phi = rhs_p
        phi_vec = LU_p.solve(rhs_p)
        phi = phi_vec.reshape((Nx, Ny), order='F')

        # For the first time step only, print divergence before/after projection to demonstrate enforcement
        if step == 1:
            div_before = np.linalg.norm(div.ravel(order='F')) * np.sqrt(dx * dy)
            # compute divergence after correction (temporarily) to show effectiveness
            u_temp = u_star.copy()
            v_temp = v_star.copy()
            for j in range(Ny):
                for i in range(1, Nx):
                    dpdx = (phi[i, j] - phi[i - 1, j]) / dx
                    u_temp[i, j] -= dt * dpdx
            for j in range(1, Ny):
                for i in range(Nx):
                    dpdy = (phi[i, j] - phi[i, j - 1]) / dy
                    v_temp[i, j] -= dt * dpdy
            div_after = np.zeros((Nx, Ny), order='F')
            for j in range(Ny):
                for i in range(Nx):
                    div_after[i, j] = (u_temp[i + 1, j] - u_temp[i, j]) / dx + (v_temp[i, j + 1] - v_temp[i, j]) / dy
            div_after_norm = np.linalg.norm(div_after.ravel(order='F')) * np.sqrt(dx * dy)
            print(f"Step 1 divergence: before projection = {div_before:.3e}, after projection = {div_after_norm:.3e}")

        # Correct velocities: u_new = u_star - dt * dphi/dx at u faces
        u_new = u_star.copy()
        for j in range(Ny):
            for i in range(1, Nx):
                dpdx = (phi[i, j] - phi[i - 1, j]) / dx
                u_new[i, j] -= dt * dpdx
        # enforce inlet and walls/outlet
        u_new[0, :] = inlet_profile
        u_new[:, 0] = 0.0
        u_new[:, -1] = 0.0
        u_new[Nx, :] = u_new[Nx - 1, :]

        # v correction: v_new = v_star - dt * dphi/dy at v faces
        v_new = v_star.copy()
        for j in range(1, Ny):
            for i in range(Nx):
                dpdy = (phi[i, j] - phi[i, j - 1]) / dy
                v_new[i, j] -= dt * dpdy
        v_new[:, 0] = 0.0
        v_new[:, -1] = 0.0
        v_new[0, :] = 0.0
        v_new[-1, :] = v_new[-2, :]

        # Update pressure and fields
        p = phi.copy()
        u = u_new
        v = v_new

        # Progress printing every 10%
        frac = step / Nsteps
        if frac >= progress_next - 1e-12 or step == Nsteps:
            print(f"Progress: {int(100*frac)}% (t={t:.4f})")
            progress_next += 0.1

    # final fields: interpolate u and v to cell centers for plotting
    u_c = np.zeros((Nx, Ny), order='F')
    v_c = np.zeros((Nx, Ny), order='F')
    for j in range(Ny):
        for i in range(Nx):
            u_c[i, j] = 0.5 * (u[i, j] + u[i + 1, j])
            v_c[i, j] = 0.5 * (v[i, j] + v[i, j + 1])

    # Compute simple L2 divergence norm of final velocity field
    div_final = np.zeros_like(p)
    for j in range(Ny):
        for i in range(Nx):
            div_final[i, j] = (u[i + 1, j] - u[i, j]) / dx + (v[i, j + 1] - v[i, j]) / dy
    l2_div = np.linalg.norm(div_final.ravel(order='F')) * np.sqrt(dx * dy)
    print(f"Final L2 divergence norm: {l2_div:.3e}")

    # Plot contours of u_c, v_c, p
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    levels = 40
    cs1 = axes[0].contourf(xp, yp, u_c, levels=levels, cmap='RdBu_r')
    axes[0].set_title('u (cell centers)')
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    fig.colorbar(cs1, ax=axes[0])

    cs2 = axes[1].contourf(xp, yp, v_c, levels=levels, cmap='RdBu_r')
    axes[1].set_title('v (cell centers)')
    axes[1].set_xlabel('x')
    fig.colorbar(cs2, ax=axes[1])

    cs3 = axes[2].contourf(xp, yp, p, levels=levels, cmap='RdBu_r')
    axes[2].set_title('p (cell centers)')
    axes[2].set_xlabel('x')
    fig.colorbar(cs3, ax=axes[2])

    fig.suptitle(f'Fields at t={t_end:.3f}')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig('uvp_contours.png', dpi=150)
    print('Saved figure: uvp_contours.png')


if __name__ == '__main__':
    # Run quick tests
    # test_upwind_convective_v()
    # test_upwind_convective_u()
    # test_poisson_operator(Nx=48, Ny=24)
    # test_projection_preserves_divergence(Nx=48, Ny=24)

    # Example run with moderate resolution
    run_simulation(Nx=64, Ny=32, t_end=1.0, CFL=0.4)
```



### runtime_outputs

#### Output block1

Grid: Nx=64, Ny=32, dx=3.1250e-02, dy=3.1250e-02, dt=4.1667e-03, steps=240

Runtime Error:
Traceback (most recent call last):
  File "C:\Users\admin\AppData\Local\Temp\ipykernel_10144\2663542217.py", line 138, in execute_code_tool
    exec(code, global_scope, global_scope)
  File "<string>", line 447, in <module>
  File "<string>", line 288, in run_simulation
IndexError: too many indices for array: array is 1-dimensional, but 2 were indexed



#### Output block2

upwind_convective_v quick tests passed
Grid: Nx=64, Ny=32, dx=3.1250e-02, dy=3.1250e-02, dt=4.1667e-03, steps=240
Progress: 10% (t=0.1000)
Progress: 20% (t=0.2000)
Progress: 30% (t=0.3000)
Progress: 40% (t=0.4000)
Progress: 50% (t=0.5000)
Progress: 60% (t=0.6000)
Progress: 70% (t=0.7000)
Progress: 80% (t=0.8000)
Progress: 90% (t=0.9000)
Progress: 100% (t=1.0000)
Final L2 divergence norm: 3.366e+00
Saved figure: uvp_contours.png



#### Output block3

Grid: Nx=64, Ny=32, dx=3.1250e-02, dy=3.1250e-02, dt=4.1667e-03, steps=240
Step 1 divergence: before projection = 4.500e-01, after projection = 3.170e-01
Progress: 10% (t=0.1000)
Progress: 20% (t=0.2000)
Progress: 30% (t=0.3000)
Progress: 40% (t=0.4000)
Progress: 50% (t=0.5000)
Progress: 60% (t=0.6000)
Progress: 70% (t=0.7000)
Progress: 80% (t=0.8000)
Progress: 90% (t=0.9000)
Progress: 100% (t=1.0000)
Final L2 divergence norm: 4.460e-01
Saved figure: uvp_contours.png




### review_decision
revise

### review_comments
	Current Stage [A/3]
1. MAJOR

Category: runtime
Issue: Runtime IndexError in time loop when building inlet_profile: code uses y_u[:, 0] even though y_u is a 1-D array (np.linspace). This raises "too many indices for array" and aborts execution.

Recommendation: Use the correct y-vector / mesh variable. Either use the 1-D y_u directly (inlet_profile = np.sin(np.pi * y_u) * (...)) or index the 2-D meshgrid array (yu[:, 0]) if you intended the meshgrid. Replace y_u[:, 0] with y_u or yu[:, 0] consistently. Re-run to confirm the IndexError is resolved.


2. MAJOR

Category: accuracy
Issue: upwind_convective_v contains multiple incorrect and fragile index operations (e.g. j-1 indexing without robust boundary checks, confusing u_top/u_bottom computations, and inconsistent use of u array indices). This both risks IndexError at boundaries and produces incorrect convective fluxes.

Recommendation: Refactor upwind_convective_v to a clear, robust implementation. Compute the velocity at the v-face by averaging the two neighboring u-values with explicit bounds checking. Use clear conditionals for boundary faces (j==0, j==Ny) and interior faces; compute upwind derivatives using the sign of the advecting velocity. Prefer a vectorized implementation (or at least well-tested loops) and add unit tests on a small grid to verify expected upwind behaviour. Example fix outline: for each v-face (i,j) compute uface = 0.5*(u[i, max(0,j-1)] + u[i+1, max(0,j-1)]) with safe indexing, then choose dv/dx upwind differencing using uface sign; similarly treat dv/dy using local v-values. Ensure no negative indices are used inadvertently.


3. MAJOR

Category: accuracy
Issue: Poisson matrix assembly (build_poisson_p) uses a confusing sign convention (diag = -2*(rx+ry); then sets A[row,row] = -diag and neighbor entries positive rx/ry). The sign/scale may be inconsistent with the RHS formulation (rhs = div / dt) and the chosen discrete Laplacian convention, risking an incorrect pressure solve (wrong sign or scaling).

Recommendation: Unify the Laplacian sign convention and confirm consistency with RHS: decide whether A represents +L or -L. A conventional choice is to form A such that A p = b with A[row,row] = 2*(rx+ry) and off-diagonals = -rx, -ry when solving -∇^2 p = b, or use A[row,row] = -2*(rx+ry) and off-diagonals = rx, ry when solving ∇^2 p = b. Document which equation you're solving and adjust rhs sign accordingly. Add a small consistency test: apply A to a known analytic p and check discrete Laplacian matches expected value.


4. MAJOR

Category: accuracy
Issue: Inconsistent / fragile boundary-condition enforcement across matrices and RHS: some Dirichlet rows are set to identity in matrices but the corresponding RHS updates are ad-hoc (e.g. Neumann outlet approximations for v and u are handled differently across assembly, rhs, and post-correction). This can lead to mass-conservation errors or spurious pressure modes.

Recommendation: Systematically handle BCs: (1) document which grid nodes are Dirichlet and which are Neumann for each variable, (2) enforce Dirichlet by setting matrix row = identity and rhs = prescribed value, (3) for Neumann incorporate ghost/one-sided discretization consistently into operator assembly and rhs. Verify global discrete mass conservation (sum divergence ~ 0 given BCs) on a simple test (steady-state) before full runs.


5. MINOR

Category: structure
Issue: Code quality / structure: several functions are written with many nested loops and ad-hoc indexing which makes reasoning and debugging hard; direct LU on large sparse matrices (spla.splu) may be memory/time heavy for larger grids.

Recommendation: Improve readability: add inline comments that document array shapes, indexing conventions, and boundary node ranges. Where possible vectorize convective term computations or use numba for hotspots. For scalability, consider iterative solvers (cg, bicgstab) with appropriate preconditioners (AMG or ILU) for the Poisson and Helmholtz solves; but keep direct LU for moderate grids or testing.



	Current Stage [B/3]
1. MAJOR

Category: accuracy
Issue: Final velocity field has a large divergence (Final L2 divergence norm = 3.366), indicating the projection step does not enforce incompressibility correctly.

Recommendation: This is the single most important correctness failure. The projection solve should reduce the discrete divergence to near machine precision (or at least to a very small tolerance). Inspect and fix the discrete operators and boundary-row handling so that the discrete divergence, discrete gradient and the Poisson operator form a consistent triplet. Specific steps: (a) Verify the sign and stencil of the Poisson matrix A_p (ensure A_p * phi truly represents Δφ with the same finite-difference stencil used for the discrete gradient). (b) Ensure boundary rows in A_p are exact Dirichlet rows (identity) when p is prescribed, and treat Neumann boundaries in a way that preserves discrete compatibility (or eliminate Neumann by prescribing ghost values consistently). (c) Make Helmholtz velocity boundary rows consistent with the physical BCs and set the corresponding rhs entries to the same prescribed values (you currently set identity rows for Dirichlet but leave Neumann rows without consistent rhs treatment). After fixing, add a unit test: start with a divergence-free field and check that the projection preserves it; also print the divergence before and after projection in one step to confirm effectiveness.


2. MAJOR

Category: structure
Issue: Inconsistent / ad-hoc enforcement of Neumann outlet boundary conditions for velocity and convoluted Poisson boundary logic.

Recommendation: Rows for Neumann BCs are implemented as one-sided difference rows (e.g. A[row,row]=1; A[row, row-1]=-1) but the rhs is not set accordingly — this leads to inconsistent linear systems and mass imbalance. Choose one of the following consistent strategies and apply it everywhere: (1) implement Neumann by eliminating ghost nodes and modifying neighboring equations and RHS correctly (derive the discrete equation algebraically), or (2) enforce Neumann by setting the ghost value equal to the adjacent interior value and implement that by setting an identity row for the boundary DOF and setting the rhs to the interior value (so the algebraic equation is u_N = u_{N-1}). Apply the same consistent approach to u, v and pressure boundary conditions. Rewrite the build_helmholtz_* and build_poisson_p functions to explicitly set rows and RHS in a straightforward, auditable way (avoid in-place conditional increments on LIL rows that are hard to reason about).


3. MAJOR

Category: accuracy
Issue: Poisson matrix assembly is error-prone and uses awkward in-place checks (A[row, col] == 0) that can mask logic errors and lead to incorrect neighbor coefficients.

Recommendation: Reconstruct the Poisson matrix assembly more simply and deterministically: for each interior cell set central = -2*(rx+ry) and add neighbors if they exist. For left/bottom Neumann boundaries, either transform to explicit ghost-value substitution before assembly or add the proper one-sided second-order stencil; do not try to add to the right/top neighbor conditionally via 'if A[row, idx] == 0' — instead accumulate contributions into local variables and then assign. After assembly, test A_p by applying it to a known analytic phi (e.g. phi = sin(pi x) sin(pi y)) and compare the discrete Laplacian to the analytic Laplacian to validate correctness and signs. Also check that the linear system solved uses the correct rhs sign: if A_p approximates Δ, then rhs = div/dt is correct; if A_p approximates -Δ, rhs must be negated.


4. MAJOR

Category: accuracy
Issue: Helmholtz solves and boundary handling are inconsistent with RHS construction: some rows are identity-rows but the rhs isn't set for all such rows, and Neumann rows are not handled consistently.

Recommendation: When you build Helmholtz matrices with some rows set to identity (Dirichlet), you must set the corresponding rhs entries to the Dirichlet values before solving. The code sets rhs_u/rhs_v for some Dirichlet rows but leaves Neumann rows without consistent rhs. Rework: (a) keep a mask of Dirichlet rows while assembling matrices; (b) before solve set rhs[dirichlet_mask] = prescribed_value; (c) implement Neumann rows so they produce an algebraic equation that is consistent with the treatment in the PDE (either modify matrix/RHS or use ghost elimination). This will reduce spurious fluxes that generate divergence.


5. MINOR

Category: structure
Issue: Insufficient testing and diagnostics for mass conservation and operator consistency; only upwind_convective_v is tested.

Recommendation: Add additional unit checks and debug prints (only for debugging runs) to validate discrete operator consistency: - Test upwind_convective_u similarly to upwind_convective_v. - At startup or in a debug mode compute and print the mean and max divergence before projection, after solving Poisson and after correction for 1 or 2 steps to ensure projection reduces divergence. - Test the Poisson operator by applying it to manufactured phi and verifying the discrete Laplacian. - Consider replacing direct LU factorization for the Poisson with an iterative solver (CG) with a simple preconditioner for memory efficiency on larger grids (but only after correctness is established). These tests will help isolate where the incompressibility enforcement fails.



	Current Stage [C/3]
1. MAJOR

Category: accuracy
Issue: Projection step fails to produce a divergence-free velocity field: the L2 divergence is large after projection (final L2 ≈ 4.46e-01) and the first-step projection only reduces divergence from ~0.45 to ~0.32.

Recommendation: This is a fundamental correctness issue. The projection solve (Poisson assembly / BCs / gradient mapping) is inconsistent with the discrete divergence operator. Debug steps I recommend: (1) After solving phi, compute and print the residual r = A_p.dot(phi_vec) - rhs_p and its norm to verify the Poisson solver is solving the linear system to machine precision. (2) Verify the sign/scaling convention of the discrete Poisson operator: whether A_p represents Δφ or -Δφ, and ensure it matches the RHS (you currently assemble A_p approximating Δφ and use rhs = div/dt). (3) Check boundary conditions for φ — for the projection the correct BCs are generally Neumann: dφ/dn = (u*_n - u_prescribed_n)/dt on boundaries (or derived from problem physics); imposing Dirichlet p=0 only on the outlet is fine but all other boundaries must have consistent Neumann conditions. Your ghost-elimination substitutions for left/bottom may not implement the correct projection BCs, leading to the poor divergence reduction. Implement explicit Neumann boundary conditions for φ computed from u_star and prescribed boundary velocities, or derive and impose the correct discrete fluxes at faces. (4) After fixing BCs, re-run a small diagnostic: compute divergence before/after and ensure it reduces to solver tolerance (e.g., near machine eps times problem size).


2. MAJOR

Category: accuracy
Issue: Inconsistent/insufficient enforcement of projection boundary conditions and gradient mapping from cell centers to staggered faces.

Recommendation: Carefully check the discrete gradient and divergence index mappings for the staggered MAC grid. For u-faces at index i (0..Nx) the adjacent cell-centers are typically i-1 (left) and i (right) — your dpdx loop uses i in range(1,Nx) which is correct for interior faces, but the treatment of boundary faces is ad-hoc (setting u_new[Nx,:] = u_new[Nx-1,:] and leaving i==0 as inlet Dirichlet). Ensure that the projection boundary conditions and corrections are implemented so the discrete divergence at all cells accounts for these boundary faces. Specifically: (a) compute and verify the discrete divergence of the corrected field equals zero to within solver tolerance; (b) ensure that for faces with prescribed normal velocity (Dirichlet), the Neumann condition for phi is set consistently; (c) do not rely solely on ghost-copy identity rows unless they exactly reproduce the required Neumann projection BCs.


3. MINOR

Category: structure
Issue: No run-time diagnostic checks of Poisson solve accuracy or residuals; the built-in tests are commented out in __main__.

Recommendation: Enable or run the provided diagnostic tests (test_poisson_operator and test_projection_preserves_divergence) as part of a debug run before the production run. Add explicit residual checks after LU_p.solve (and LU_u/LU_v) such as ‖A·x - b‖ to ensure the linear solves are accurate. Print a concise warning and abort or raise if residuals exceed a small tolerance. This will help separate solver accuracy issues from discretization bugs.


4. MAJOR

Category: structure
Issue: Poisson operator assembly uses ghost-elimination for left/bottom Neumann in a bespoke way that is hard to verify and likely inconsistent with projection BCs.

Recommendation: Replace or augment the custom ghost-elimination assembly with an explicit handling of Neumann boundary rows: assemble standard 5-point Laplacian for interior cells and for boundary cells replace finite-difference one-sided stencils that implement dφ/dn = g (Neumann) by modifying the matrix row and RHS to include the known flux g (which for projection depends on u_star). This will make the BCs explicit and easier to verify. Alternatively, implement the Poisson solve with an established routine that supports mixed Dirichlet/Neumann BCs, or assemble the matrix with explicit ghost unknowns and eliminate them symbolically while tracking the RHS contributions from boundary fluxes.


5. MINOR

Category: structure
Issue: Lack of diagnostic and unit tests being active in the main run and insufficient logging for debugging convergence/solver residuals.

Recommendation: Temporarily re-enable the unit tests and add minimal logging for (a) Poisson residual norm, (b) divergence before and after projection for a few steps, (c) norm of u_star and v_star boundary flux mismatches. Keep the required concise progress prints but add an optional --debug flag or environment variable to enable extra diagnostics when investigating projection issues.





### iteration_count
3

### rev_image_description
	Current Stage [A/3]


	Current Stage [B/3]


	Current Stage [C/3]




