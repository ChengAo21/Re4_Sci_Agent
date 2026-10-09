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
We must solve the 2D incompressible Navier–Stokes equations around a circular cylinder (diameter d = 1.0) in a rectangular domain using a numerical method that: enforces incompressibility, imposes inflow/outflow/free-slip/no-slip boundaries, represents the cylinder via an immersed boundary with volume penalization, and runs stably up to t = 20.0 for Re = 100 (visc = 1/Re). Primary requirements and constraints from the statement: use a staggered (MAC) grid, treat penalization implicitly (incorporate into the momentum matrix diagonal), rely solely on penalization to enforce the cylinder no-slip, fix pressure at a single node to remove the nullspace, explicitly pass constants or define them globally (visc must be defined globally), print progress every 10% of total steps, and save the final velocity magnitude contour (RdBu_r) with the cylinder marked. Additional hints constrain solver call signatures (use atol not tol) and avoid removed NumPy functions. We are asked to implement a stable, efficient, and reasonably accelerated solver; the discretization and solver choices must balance accuracy, stability and computational cost. Primary mathematical and numerical challenges: - Nonlinearity & convection: convective acceleration terms (u·∇u) are nonlinear and require stable time-integration; explicit treatment yields strict CFL limits; implicit treatment requires nonlinear solves or linearization. - Incompressibility constraint: divergence-free velocity requires solving a pressure/Poisson problem or a saddle-point system; the pressure has a nullspace (constant offset), so the pressure must be pinned at one node. - Immersed boundary with penalization: the cylinder is represented via a mask χ(x,y); penalization adds a stiff term (χ/η) (u - u_b) that enforces no-slip inside the body; if treated explicitly this would force very small Δt; implicit treatment reduces the time-step constraint but requires modifying linear operators. - Staggered (MAC) grid coupling: velocity components and pressure are stored at different locations; discrete divergence/gradient operators must be consistent to maintain discrete conservation and pressure stability. - Boundary conditions: free-slip at top/bottom (v = 0, ∂u/∂y = 0), Dirichlet inlet (u=1), Neumann outlet for pressure (∂p/∂n = 0) — these must be discretized consistently in the projection algorithm. - Linear algebra cost: projection requires solving a Poisson or saddle-point system each time step; diffusion/penalization solves require solving large sparse linear systems. For tractability we must reuse factorizations or use efficient iterative solvers with good preconditioning (algebraic multigrid or incomplete factorizations). - Accuracy vs cost: second-order time-stepping (BDF2 or Crank–Nicolson/Adams–Bashforth) improves accuracy for long-time integration but complicates matrix reuse and solver preconditioning. - Implementation constraints: enforce explicit passing/definition of constants and use visc globally; adhere to specified SciPy solver arguments (atol) and avoid deprecated NumPy functions. Numerical parameters to choose in the implementation (examples): grid resolution (Nx, Ny) sufficiently fine around cylinder (recommendations below), penalization parameter η (small, e.g., 1e-3–1e-4), time-step selection via CFL for convective terms and viscous stability (with implicit diffusion the Δt can be larger), and tolerance/solver settings for iterative solvers (atol). Practical acceleration measures: reuse LU/Cholesky factors of constant matrices, precompute and cache sparse matrix factorizations, use algebraic multigrid for Poisson, choose semi-implicit schemes to keep convective CFL manageable, optionally perform local grid refinement near the cylinder. Output requirement: save a contour plot of velocity magnitude at final time with the cylinder marked (use RdBu_r colormap) and do not call plt.show().

### solution_plans
	Current Stage [A/2]
solu_name='MAC-projection with semi-implicit AB2/CN time stepping and implicit penalization' content='Governing idea:\n- Use a staggered (MAC) grid for discretization of u, v and p to maintain discrete conservation and avoid spurious pressure modes. Advance velocities with a semi-implicit scheme: explicit Adams–Bashforth-2 (AB2) for convection, Crank–Nicolson (CN) or implicit Euler for viscous diffusion, and treat the penalization term implicitly by adding χ/η to the diagonal of the velocity linear system. Enforce incompressibility with a projection (pressure-Poisson) step; fix pressure at one node to remove the nullspace. Reuse sparse factorizations (or AMG operators) for the constant linear solves to accelerate time stepping.\n\nAlgorithmic steps (numbered):\n1) Spatial discretization: build a uniform MAC grid (Nx × Ny pressure cells, u on vertical faces, v on horizontal faces). Construct discrete operators: gradient G (pressure→face), divergence D (face→cell), Laplacians L_u, L_v for faces using consistent second-order central differences, and interpolations for convective fluxes.\n2) Cylinder mask and penalization: build boolean mask χ defined at face locations (u- and v-nodes). Choose penalization parameter eta (η, e.g., 1e-3–1e-4). The penalization term for velocity becomes (χ/η)*(u - u_b) where u_b = 0 on the cylinder. Incorporate χ/η into diagonal of velocity operator to treat it implicitly.\n3) Time-stepping selection: choose Δt based on convective CFL (|u|max*Δt/Δx ≤ CFL_conv, CFL_conv≈0.4–0.6) and optionally viscous criterion; because diffusion and penalization are implicit, Δt can be moderately large but keep Δt small enough for convective accuracy. Use AB2 for nonlinear terms: N^{n+1/2} ≈ 1.5N^n - 0.5N^{n-1}. For the first step use RK2 or a single-step method.\n4) Build constant linear matrices once: A_u = I - (Δt/2) visc * L_u + Δt * diag(χ/η) (for CN) and similarly for A_v. These matrices remain constant if visc, Δt and χ are constant; factorize (LU or ILU) or prepare AMG preconditioners and cache them to use every time-step. For the Poisson matrix P = D * A_face^{-1} * G (approximate) or the discrete Laplacian for pressure use AMG or sparse direct factorization; pin one pressure DOF (set pressure at one cell to zero) to remove nullspace and make the Poisson unique.\n5) Time-step loop until t >= 20.0 with step counter and printing at every 10% of total steps: (a) compute convective terms N^{n+1/2} using AB2 extrapolation; (b) form right-hand sides for tentative velocity solves including explicit convective and implicit viscous/penalization contributions; (c) solve A_u * u* = RHS_u and A_v * v* = RHS_v using cached factorizations or iterative solves (use atol argument in SciPy solvers, not tol); (d) assemble right-hand side for the pressure Poisson ∇·(u*/Δt) (with consistent boundary conditions), solve Poisson for pressure correction p^{n+1} (fix one pressure node to zero); (e) update velocities: u^{n+1} = u* - Δt * G_x p^{n+1}, v^{n+1} = v* - Δt * G_y p^{n+1}; (f) enforce boundary conditions at faces (inlet Dirichlet, free-slip on top/bottom, penalization implicitly enforces cylinder no-slip); (g) update time and AB2 history.\n6) Diagnostics and output: compute global CFL, kinetic energy, drag/lift integrals if desired (use appropriate interpolation onto cylinder surface if needed), and print concise progress every 10% of steps. At final time compute velocity magnitude on cell centers, generate contour with RdBu_r colormap, overplot cylinder contour from mask and save the figure (no plt.show()).\n\nStability / accuracy / complexity / efficiency limitations:\n- Stability: explicit AB2 for convection imposes a convective CFL limit; although viscous and penalization terms are implicit, convective explicitness limits Δt. For Re=100 a moderate Δt is feasible, but highly under-resolved grids may require smaller Δt. The implicit penalization prevents the penalization term from forcing Δt → 0. Boundary condition discretization on MAC must be consistent or spurious oscillations can appear.\n- Accuracy: AB2/CN is formally second order in time (aside from the first-step startup). Spatial discretization is second-order; immersed boundary penalization is first- to second-order accurate depending on mask smoothing and the way face-centered mask is constructed — sharp masks produce O(Δx) boundary-layer errors inside the penalized region. Grid resolution near the cylinder controls drag/lift accuracy.\n- Complexity: the method requires assembling and factoring (or preconditioning) three major sparse systems per time-step: two for velocities and one for pressure. Reusing LU/AMG preconditioners reduces amortized cost but initial setups are nontrivial. Memory for factorization can be large if the grid is very fine.\n- Efficiency: reuse of factorizations/AMG hierarchies makes each time-step cheap, but changes in Δt, visc or mask (if moving body) invalidate the factorization. Algebraic multigrid is highly recommended for the Poisson solve; sparse direct solvers are robust but may be memory-limited. The penalization parameter η trades enforcement strength vs condition number: very small η enforces no-slip better but increases conditioning of the velocity system; implicit diagonal treatment alleviates but does not eliminate conditioning issues. Choosing moderate η (1e-3 to 1e-4) is a practical compromise.\n\nStop.'

	Current Stage [B/2]
solu_name='Implicit BDF2 projection with Krylov block solvers and AMG + local refinement acceleration' content='Governing idea:\n- Use a fully implicit second-order backward differentiation (BDF2) treatment of transient + viscous + penalization terms while keeping convective terms implicit via Newton or simplified Picard linearization. Enforce incompressibility using an incremental pressure-correction (rotational form) to reduce splitting errors. Solve the coupled linear systems arising at each (nonlinear) Newton/Picard iteration using Krylov solvers (GMRES) with physics-aware block preconditioning: block-diagonal approximate inverse of the viscous operator for velocity and AMG for the Schur complement/Poisson part. Use block-structured local mesh refinement (or clustering) near the cylinder to reduce DOF while maintaining resolution where needed. Implicit penalization remains in the diagonal of the velocity block to avoid restrictive time-steps.\n\nAlgorithmic steps (numbered):\n1) Spatial discretization: choose MAC grid with optional block-structured refinement near the cylinder (coarsen away from cylinder). Build consistent discrete operators and prolongation/restriction between refinement levels. Construct face-centered mask χ and mapping to refined grids.\n2) Time discretization: BDF2 for time derivative: (3u^{n+1} - 4u^n + u^{n-1})/(2Δt). Treat viscous and penalization terms fully implicitly. Linearize convective term at each time-step using a single Picard iteration (u·∇)_{n+1} ≈ u^{*}·∇ u^{n+1} with u^{*} from previous time-level to avoid costly Newton while retaining stability.\n3) Build coupled saddle-point linear system for (u^{n+1}, v^{n+1}, p^{n+1}) (or incremental form that reduces to solving velocity block then Schur complement). The velocity block includes mass, viscous Laplacian, and diagonal penalization (χ/η). The Schur complement approximates the pressure Laplacian.\n4) Solver strategy per time-step: (a) form the linear system or its block components; (b) precondition using a block preconditioner: approximate inverse of velocity block via AMG or pre-factored ILU; approximate Schur complement by pressure mass-lumped operator or AMG Poisson; (c) solve the linear system with GMRES until residual atol is met (use SciPy or PETSc interfaces but pass atol not tol). Pin one pressure DOF to remove nullspace when assembling the Schur complement solve.\n5) Nonlinear iteration: perform 1–3 Picard iterations per time-step (often a single Picard step is sufficient for moderate Re) to update linearization; update boundary conditions after each Newton/Picard step.\n6) Acceleration techniques: (a) reuse preconditioners across time-steps since viscous and penalization coefficients are constant; (b) if available, use algebraic multigrid (pyamg or AMG4Py) for both velocity-block approximate inverses and pressure Poisson; (c) use local refinement to reduce global DOF while keeping the cylinder boundary well resolved; (d) adaptive time stepping controlling Δt by convective CFL and nonlinear residual decrease; print progress every 10% of total steps.\n7) Final output: compute velocity magnitude on a chosen uniform or composite grid for visualization; plot contour using RdBu_r and overlay the cylinder boundary (circumference) as a circle; save figure only.\n\nStability / accuracy / complexity / efficiency limitations:\n- Stability: fully implicit BDF2 with Picard linearization is more stable than semi-implicit AB2 for stiff problems and allows larger Δt; however, linearization error makes implicit treatment less accurate unless enough inner iterations are performed. The choice of number of Picard/Newton iterations per time-step trades computational cost vs nonlinear accuracy.\n- Accuracy: BDF2 is second-order accurate in time; the Picard linearization may reduce overall order if not iterated sufficiently. Spatial accuracy remains second-order on uniform MAC grids; with local refinement special interpolation must preserve conservative properties to avoid spurious divergence.\n- Complexity: constructing and preconditioning the coupled saddle-point problem is more complex to implement than a simple projection; block preconditioners and AMG setup cost are higher but amortized over many time-steps. Local refinement increases implementation complexity (mesh management, prolongation operators) but reduces DOF for a target accuracy.\n- Efficiency: GMRES with good block preconditioning can be very efficient per unknown (especially combined with AMG), but convergence depends on the quality of the Schur complement approximation and penalization conditioning. Implicit penalization in the velocity block keeps Δt large, but very small η increases conditioning and may require more robust preconditioning. Reusing preconditioners between steps reduces runtime; however, changing Δt or adaptive refinement invalidates reuse.\n\nStop.'



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
This solver implements a semi-implicit MAC-grid Navier–Stokes solver with immersed-boundary volume penalization for flow around a cylinder. Core flow: convection is treated explicitly (AB2 extrapolation), viscous and penalization terms treated implicitly (single sparse factorization reused each timestep), and pressure correction via Poisson projection (pressure pinned at the cell nearest the cylinder center). Staggered-grid DOF ordering uses Fortran (column-major) flattening consistently with index-mapping functions, avoiding mismatches between sparse-matrix assembly and vector layouts. Cylinder masks are built with meshgrid(indexing='ij') so mask shapes match face-centered arrays. Convective interpolations are vectorized (no Python nested loops). Matrices Au, Av and Lp are assembled and factorized once before time-stepping; penalization enters matrices via diagonal terms (implicit). Diagnostics print at 10% intervals; final velocity magnitude contour is saved with the cylinder boundary overlay.



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


# Index mapping consistent with Fortran-order flattening (i fastest)
def idx_u(i, j, Nx, Ny):
    return j * (Nx + 1) + i


def idx_v(i, j, Nx, Ny):
    return j * Nx + i


def idx_p(i, j, Nx, Ny):
    return j * Nx + i


def build_cylinder_masks(x_u, y_u, x_v, y_v, radius=0.5):
    # Use 'ij' indexing so Xu.shape == (len(x_u), len(y_u)) == (Nx+1, Ny)
    Xu, Yu = np.meshgrid(x_u, y_u, indexing='ij')
    Xv, Yv = np.meshgrid(x_v, y_v, indexing='ij')

    chi_u = ((Xu**2 + Yu**2) <= radius**2).astype(float)
    chi_v = ((Xv**2 + Yv**2) <= radius**2).astype(float)

    return chi_u, chi_v


def build_laplacian_u(Nx, Ny, dx, dy, inlet_mask):
    Nu = (Nx + 1) * Ny
    data = []
    rows = []
    cols = []

    for j in range(Ny):
        for i in range(Nx + 1):
            row = idx_u(i, j, Nx, Ny)
            if inlet_mask[row]:
                rows.append(row); cols.append(row); data.append(1.0)
                continue

            center = 0.0
            # left
            if i - 1 >= 0:
                left = idx_u(i - 1, j, Nx, Ny)
                rows.append(row); cols.append(left); data.append(1.0 / dx**2)
                center += -1.0 / dx**2
            # right
            if i + 1 <= Nx:
                right = idx_u(i + 1, j, Nx, Ny)
                rows.append(row); cols.append(right); data.append(1.0 / dx**2)
                center += -1.0 / dx**2
            # down
            if j - 1 >= 0:
                down = idx_u(i, j - 1, Nx, Ny)
                rows.append(row); cols.append(down); data.append(1.0 / dy**2)
                center += -1.0 / dy**2
            # up
            if j + 1 <= Ny - 1:
                up = idx_u(i, j + 1, Nx, Ny)
                rows.append(row); cols.append(up); data.append(1.0 / dy**2)
                center += -1.0 / dy**2

            rows.append(row); cols.append(row); data.append(center)

    L = sp.csr_matrix((data, (rows, cols)), shape=(Nu, Nu))
    return L


def build_laplacian_v(Nx, Ny, dx, dy, v_dirichlet_mask):
    Nv = Nx * (Ny + 1)
    data = []
    rows = []
    cols = []

    for j in range(Ny + 1):
        for i in range(Nx):
            row = idx_v(i, j, Nx, Ny)
            if v_dirichlet_mask[row]:
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
                rows.append(row); cols.append(left); data.append(1.0 / dx**2)
                center += -1.0 / dx**2
            if i + 1 <= Nx - 1:
                right = idx_p(i + 1, j, Nx, Ny)
                rows.append(row); cols.append(right); data.append(1.0 / dx**2)
                center += -1.0 / dx**2
            if j - 1 >= 0:
                down = idx_p(i, j - 1, Nx, Ny)
                rows.append(row); cols.append(down); data.append(1.0 / dy**2)
                center += -1.0 / dy**2
            if j + 1 <= Ny - 1:
                up = idx_p(i, j + 1, Nx, Ny)
                rows.append(row); cols.append(up); data.append(1.0 / dy**2)
                center += -1.0 / dy**2

            rows.append(row); cols.append(row); data.append(center)

    L = sp.csr_matrix((data, (rows, cols)), shape=(Np, Np))
    return L


def divergence(u, v, Nx, Ny, dx, dy):
    # u shape: (Nx+1, Ny), v shape: (Nx, Ny+1)
    return (u[1:Nx+1, :] - u[0:Nx, :]) / dx + (v[:, 1:Ny+1] - v[:, 0:Ny]) / dy


def gradient_p_to_u(p, Nx, Ny, dx, dy):
    dpdx = np.zeros((Nx + 1, Ny))
    dpdx[1:Nx, :] = (p[1:Nx, :] - p[0:Nx-1, :]) / dx
    dpdx[0, :] = (p[0, :] - p[0, :]) / dx
    dpdx[Nx, :] = (p[Nx - 1, :] - p[Nx - 1, :]) / dx
    return dpdx


def gradient_p_to_v(p, Nx, Ny, dx, dy):
    dpdy = np.zeros((Nx, Ny + 1))
    if Ny >= 2:
        dpdy[:, 1:Ny] = (p[:, 1:Ny] - p[:, 0:Ny-1]) / dy
    dpdy[:, 0] = (p[:, 0] - p[:, 0]) / dy
    dpdy[:, Ny] = (p[:, Ny - 1] - p[:, Ny - 1]) / dy
    return dpdy


# Flatten/reshape consistently in Fortran order to match idx_* ordering
def flatten(arr):
    return arr.ravel(order='F')


def unflatten_u(vec, Nx, Ny):
    return vec.reshape((Nx + 1, Ny), order='F')


def unflatten_v(vec, Nx, Ny):
    return vec.reshape((Nx, Ny + 1), order='F')


def unflatten_p(vec, Nx, Ny):
    return vec.reshape((Nx, Ny), order='F')


def compute_convective(u, v, Nx, Ny, dx, dy):
    # u: (Nx+1,Ny), v: (Nx,Ny+1)
    # du/dx
    du_dx = np.zeros_like(u)
    du_dx[1:Nx, :] = (u[2:Nx+1, :] - u[0:Nx-1, :]) / (2 * dx)
    du_dx[0, :] = (u[1, :] - u[0, :]) / dx
    du_dx[Nx, :] = (u[Nx, :] - u[Nx - 1, :]) / dx

    # du/dy
    du_dy = np.zeros_like(u)
    if Ny >= 3:
        du_dy[:, 1:Ny-1] = (u[:, 2:Ny] - u[:, 0:Ny-2]) / (2 * dy)
    if Ny >= 2:
        du_dy[:, 0] = (u[:, 1] - u[:, 0]) / dy
        du_dy[:, Ny - 1] = (u[:, Ny - 1] - u[:, Ny - 2]) / dy

    # interpolate v to u locations: v_at_u shape (Nx+1,Ny)
    v_at_u = np.zeros_like(u)
    # interior faces
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
        dv_dy[:, Ny] = (v[:, Ny] - v[:, Ny - 1]) / dy
    else:
        dv_dy[:] = 0.0

    # dv/dx
    dv_dx = np.zeros_like(v)
    if Nx >= 3:
        dv_dx[1:Nx-1, :] = (v[2:Nx, :] - v[0:Nx-2, :]) / (2 * dx)
    if Nx >= 2:
        dv_dx[0, :] = (v[1, :] - v[0, :]) / dx
        dv_dx[Nx-1, :] = (v[Nx-1, :] - v[Nx-2, :]) / dx

    # interpolate u to v locations: u_at_v shape (Nx, Ny+1)
    u_at_v = np.zeros_like(v)
    # interior in y: j=1..Ny-1
    if Ny >= 2:
        # four-point average for interior faces
        u_at_v[:, 1:Ny] = 0.25 * (
            u[0:Nx, 0:Ny-1] + u[1:Nx+1, 0:Ny-1] + u[0:Nx, 1:Ny] + u[1:Nx+1, 1:Ny]
        )
    # bottom and top faces: use 2-point average
    u_at_v[:, 0] = 0.5 * (u[0:Nx, 0] + u[1:Nx+1, 0])
    u_at_v[:, Ny] = 0.5 * (u[0:Nx, Ny-1] + u[1:Nx+1, Ny-1])

    N_v = u_at_v * dv_dx + v * dv_dy
    return N_u, N_v


def main():
    x1 = 3.0; x2 = 8.0; y1 = 3.0; y2 = 3.0
    # moderate grid for runtime
    Nx = 220
    Ny = 120

    dx, dy, x_u, y_u, x_v, y_v, x_p, y_p = build_grid(x1, x2, y1, y2, Nx=Nx, Ny=Ny)

    # masks for cylinder
    chi_u, chi_v = build_cylinder_masks(x_u, y_u, x_v, y_v, radius=0.5)

    Nu = (Nx + 1) * Ny
    Nv = Nx * (Ny + 1)
    Np = Nx * Ny

    # inlet u-face flattened indices
    inlet_u_indices = [idx_u(0, j, Nx, Ny) for j in range(Ny)]
    inlet_u_mask = np.zeros(Nu, dtype=bool)
    inlet_u_mask[inlet_u_indices] = True

    # v Dirichlet (top/bottom)
    v_zero_mask = np.zeros(Nv, dtype=bool)
    for i in range(Nx):
        v_zero_mask[idx_v(i, 0, Nx, Ny)] = True
        v_zero_mask[idx_v(i, Ny, Nx, Ny)] = True

    print('Building operators...')

    L_u = build_laplacian_u(Nx, Ny, dx, dy, inlet_u_mask)
    L_v = build_laplacian_v(Nx, Ny, dx, dy, v_zero_mask)

    # choose pressure pin near (0,0)
    i0 = int(np.argmin(np.abs(x_p - 0.0)))
    j0 = int(np.argmin(np.abs(y_p - 0.0)))
    pinned_idx = idx_p(i0, j0, Nx, Ny)
    L_p = build_laplacian_p(Nx, Ny, dx, dy, pinned_idx=pinned_idx)

    # Time-stepping
    T = 20.0
    Umax = 1.0
    CFL = 0.45
    dt = CFL * min(dx, dy) / Umax
    nsteps = int(floor(T / dt))
    dt = T / nsteps

    print(f'Nx={Nx} Ny={Ny} dx={dx:.4f} dy={dy:.4f} dt={dt:.5f} steps={nsteps}')

    # Build velocity solve matrices once
    Iu = sp.eye(Nu, format='csr')
    Iv = sp.eye(Nv, format='csr')

    chi_u_flat = chi_u.ravel(order='F')
    chi_v_flat = chi_v.ravel(order='F')

    Pen_u = sp.diags(dt * (chi_u_flat / eta), format='csr')
    Pen_v = sp.diags(dt * (chi_v_flat / eta), format='csr')

    Au = Iu - dt * visc * L_u + Pen_u
    Av = Iv - dt * visc * L_v + Pen_v

    # enforce Dirichlet rows in Au/Av
    Au = Au.tolil()
    for idx in inlet_u_indices:
        Au.rows[idx] = [idx]
        Au.data[idx] = [1.0]
    Au = Au.tocsc()

    Av = Av.tolil()
    v_zero_indices = np.nonzero(v_zero_mask)[0].tolist()
    for idx in v_zero_indices:
        Av.rows[idx] = [idx]
        Av.data[idx] = [1.0]
    Av = Av.tocsc()

    print('Factorizing velocity matrices...')
    Au_fac = spla.splu(Au)
    Av_fac = spla.splu(Av)

    print('Factorizing pressure Poisson...')
    Lp_csc = L_p.tocsc()
    Lp_fac = spla.splu(Lp_csc)

    # initialize
    u = np.ones((Nx + 1, Ny)) * 1.0
    v = np.zeros((Nx, Ny + 1))
    p = np.zeros((Nx, Ny))

    # zero inside cylinder initially via penalization mask
    u = u * (1.0 - chi_u)
    v = v * (1.0 - chi_v)

    N_u, N_v = compute_convective(u, v, Nx, Ny, dx, dy)
    N_u_prev = N_u.copy()
    N_v_prev = N_v.copy()

    u_flat = flatten(u)
    v_flat = flatten(v)

    inlet_u_flat_indices = inlet_u_indices
    v_zero_flat_indices = v_zero_indices

    print('Starting time-stepping...')
    progress_interval = max(1, nsteps // 10)
    for step in range(1, nsteps + 1):
        N_u, N_v = compute_convective(u, v, Nx, Ny, dx, dy)
        N_u_ex = 1.5 * N_u - 0.5 * N_u_prev
        N_v_ex = 1.5 * N_v - 0.5 * N_v_prev

        rhs_u = flatten(u) - dt * flatten(N_u_ex)
        rhs_v = flatten(v) - dt * flatten(N_v_ex)

        for idx in inlet_u_flat_indices:
            rhs_u[idx] = 1.0
        for idx in v_zero_flat_indices:
            rhs_v[idx] = 0.0

        u_star_flat = Au_fac.solve(rhs_u)
        v_star_flat = Av_fac.solve(rhs_v)

        u_star = unflatten_u(u_star_flat, Nx, Ny)
        v_star = unflatten_v(v_star_flat, Nx, Ny)

        div_u_star = divergence(u_star, v_star, Nx, Ny, dx, dy)
        rhs_p = (1.0 / dt) * flatten(div_u_star)
        rhs_p[pinned_idx] = 0.0

        p_corr_flat = Lp_fac.solve(rhs_p)
        p_corr = unflatten_p(p_corr_flat, Nx, Ny)

        dpdx = gradient_p_to_u(p_corr, Nx, Ny, dx, dy)
        dpdy = gradient_p_to_v(p_corr, Nx, Ny, dx, dy)

        u_new = u_star - dt * dpdx
        v_new = v_star - dt * dpdy

        # enforce BCs
        u_new[0, :] = 1.0
        u_new[:, 0] = u_new[:, 1]
        u_new[:, Ny - 1] = u_new[:, Ny - 2]
        u_new[Nx, :] = u_new[Nx - 1, :]
        v_new[:, 0] = 0.0
        v_new[:, Ny] = 0.0

        u = u_new.copy()
        v = v_new.copy()
        p = p_corr.copy()

        N_u_prev = N_u.copy()
        N_v_prev = N_v.copy()

        if (step % progress_interval) == 0 or step == nsteps:
            perc = int(100 * step / nsteps)
            uc = 0.5 * (u[0:Nx, :] + u[1:Nx + 1, :])
            vc = 0.5 * (v[:, 0:Ny] + v[:, 1:Ny + 1])
            KE = 0.5 * np.sum(uc**2 + vc**2) * dx * dy
            umax = max(np.abs(u).max(), np.abs(v).max())
            print(f'Step {step}/{nsteps} ({perc}%) t={step*dt:.3f} umax={umax:.3f} KE={KE:.5f}')

    # final field and plotting
    u_center = 0.5 * (u[0:Nx, :] + u[1:Nx + 1, :])
    v_center = 0.5 * (v[:, 0:Ny] + v[:, 1:Ny + 1])
    vel_mag = np.sqrt(u_center**2 + v_center**2)

    Xp, Yp = np.meshgrid(x_p, y_p, indexing='xy')

    fig, ax = plt.subplots(figsize=(10, 4.5))
    cp = ax.contourf(Xp, Yp, vel_mag.T, levels=50, cmap='RdBu_r')
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
  File "C:\Users\admin\AppData\Local\Temp\ipykernel_47416\2663542217.py", line 138, in execute_code_tool
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
Step 264/888 (29%) t=5.946 umax=1.435 KE=35.95108
Step 352/888 (39%) t=7.928 umax=1.479 KE=36.52780
Step 440/888 (49%) t=9.910 umax=1.486 KE=37.12368
Step 528/888 (59%) t=11.892 umax=1.491 KE=37.85590
Step 616/888 (69%) t=13.874 umax=1.539 KE=38.46993
Step 704/888 (79%) t=15.856 umax=1.696 KE=40.01608
Step 792/888 (89%) t=17.838 umax=1.660 KE=40.11593
Step 880/888 (99%) t=19.820 umax=1.512 KE=38.53156
Step 888/888 (100%) t=20.000 umax=1.517 KE=38.83548
Saved figure to velocity_magnitude.png




### review_decision
accept

### review_comments
	Current Stage [A/2]
1. MAJOR

Category: runtime
Issue: Runtime crash when applying cylinder masks: shapes (Nx+1,Ny) and (Ny,Nx+1) are incompatible, causing a broadcasting error at u = u * (1 - chi_u).

Recommendation: In build_cylinder_masks the use of np.meshgrid(..., indexing='xy') yields arrays with shape (len(y), len(x)). This produces chi_u with shape (Ny, Nx+1) while the code expects (Nx+1, Ny). Fix by generating masks with indexing='ij' (or transpose the result) so chi_u has shape (Nx+1, Ny) and chi_v has shape (Nx, Ny+1). Example change: Xu, Yu = np.meshgrid(x_u, y_u, indexing='ij') (and similarly for Xv, Yv). After this change the mask multiply u*(1-chi_u) will broadcast correctly.


2. MAJOR

Category: structure
Issue: Inconsistent degree-of-freedom ordering between the sparse-matrix indexing (idx_u/idx_v/idx_p) and the Python flattening/unflattening used for RHS/solution vectors. This will silently scramble data mapping between arrays and sparse matrices if not corrected.

Recommendation: Choose and use a single consistent linear ordering for flattening all staggered arrays. Two practical fixes: (a) use Fortran (column-major) flatten/reshape everywhere to match idx_* which uses index = j*(Nx+1)+i (i varies fastest), e.g. replace flatten(arr) with arr.ravel(order='F') and replace unflatten_u(vec,...) with vec.reshape((Nx+1, Ny), order='F'). (b) or change idx_* mapping to match C-order flattening. Also ensure chi_u_flat and chi_v_flat (penalization diagonals) are created with the same order (ravel(order='F')). Make these changes consistently for u, v and p vectors and for construction of Pen_u/Pen_v so matrix rows/cols align with vector entries and boundary indices.


3. MINOR

Category: accuracy
Issue: Pressure pinning uses pinned_idx=0 arbitrarily; that node may be located at a corner rather than a physically sensible reference location (e.g. near the body center). Also the Neumann/outlet approximations in the Laplacians are simplistic and may produce non-physical pressure offsets or poor accuracy at the outlet.

Recommendation: Select a pinned pressure node close to a well-defined location (for example the pressure cell nearest x= x_center and y=0). Compute pinned_idx programmatically by finding the i,j that minimize distance to (0,0). Example: find i0 = np.argmin(np.abs(x_p - 0.0)), j0 = np.argmin(np.abs(y_p - 0.0)); pinned_idx = j0*Nx + i0 (or consistent ordering). Revisit the Neumann/outlet treatments: if possible impose the outlet as a Dirichlet/Robin on pressure or use ghost-cell consistent stencils to better represent zero pressure-gradient boundary.


4. MINOR

Category: structure
Issue: Convective term implementation uses explicit Python loops over i,j which will be very slow at the selected grid sizes, and many sparse-matrix-building loops also check Python set membership repeatedly (costly).

Recommendation: Vectorize the convective/interpolation computations using NumPy array operations (central differences and averaged interpolations can be expressed without explicit Python loops). When assembling sparse operators, precompute boolean arrays for Dirichlet locations (instead of set membership inside nested loops) or build matrices using kronecker/sparse diagonal operators to speed assembly. If memory is limited, consider lowering Nx,Ny for prototyping or using iterative solvers + AMG rather than direct LU (see next item).


5. MINOR

Category: structure
Issue: Memory/time cost: factorizing full LU (spla.splu) for the velocity and pressure matrices at Nx=220, Ny=120 is expensive in time and memory and may not scale to finer meshes.

Recommendation: For larger problems prefer iterative solvers (e.g., spla.cg or spla.gmres) with appropriate preconditioners (amg via pyamg or SciPy incomplete factorization) for the Poisson and velocity solves. If you keep splu, optionally lower the resolution for testing and provide a runtime option to switch between direct and iterative solves. Also consider reusing factorization only when matrices are constant (you already do so) but ensure matrices are small enough for LU to be feasible.



	Current Stage [B/2]
1. MINOR

Category: structure
Issue: Dirichlet BCs are imposed by replacing rows only; columns are left unchanged.

Recommendation: Also zero the corresponding columns (or move boundary contributions to the RHS) when imposing Dirichlet rows. Currently you set the Dirichlet rows to identity which correctly pins the DOF, but you leave the column entries intact. That is acceptable for a direct solver and will enforce the boundary strongly, yet zeroing the columns (and keeping the diagonal 1) makes the modification explicit and avoids subtle coupling/readability concerns. Example: after setting row i to identity, set the i-th column to zero and the diagonal to 1, or subtract the known boundary contributions from the RHS of other rows.


2. MINOR

Category: accuracy
Issue: Boundary treatment for pressure gradient at domain faces is implemented implicitly (set to zero) without explicit documentation.

Recommendation: Make the outlet pressure Neumann BC explicit and document the ghost/one-sided discretization used. Right now gradient_p_to_u/gradient_p_to_v set dpdx/dpdy to zero at boundary faces which effectively imposes zero pressure gradient, but the treatment is terse. Consider adding comments and/or one-sided difference formulas consistent with the Poisson stencil to make the implementation and assumptions explicit.


3. MINOR

Category: structure
Issue: Direct LU factorization (spla.splu) is used for all linear solves which is robust but memory- and time-intensive for larger grids.

Recommendation: For larger problems consider switching the Poisson and velocity solves to iterative solvers (e.g. conjugate-gradient/BiCGStab with AMG/ILU preconditioning) instead of direct LU to reduce memory and runtime. You're already reusing factorizations which is good; mention in comments that splu was chosen for robustness on moderate grids but does not scale well. If portability is a concern, provide an option to use scipy.sparse.linalg.cg or pyamg for L_p.


4. MINOR

Category: accuracy
Issue: Explicit AB2 + central differencing for convection may be marginally non-dissipative; energy behaviour could be improved.

Recommendation: Convective discretization and time-integration are standard AB2+central differences; this works but can be slightly non-dissipative and may produce small high-frequency noise. If you observe stability/noise issues at other Re or coarser grids consider using a skew-symmetric form for the convective term or a low-dissipation flux limiter/upwinding for the convective fluxes to improve energy behaviour while keeping AB2.


5. MINOR

Category: structure
Issue: Plotting uses mixed indexing conventions and a transpose to match array orientations, which is slightly confusing.

Recommendation: The plotting uses meshgrid(..., indexing='xy') combined with vel_mag.T which is correct but easy to confuse. Consider using consistent 'ij' indexing throughout (and avoid transposing) or add a short comment to clarify why the transpose is needed so future readers are not surprised by the orientation.





### iteration_count
2

### rev_image_description
	Current Stage [A/2]


	Current Stage [B/2]




