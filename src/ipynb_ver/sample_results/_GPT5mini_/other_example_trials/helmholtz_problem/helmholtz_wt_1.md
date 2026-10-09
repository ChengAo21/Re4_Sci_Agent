### 🛌🏻 Re4gent with Multi-Modal Review

#### Configs 🏷️:
MaxCount: 2, MultiModalReview: False, AlternativeNum: 2, MaxComments: 5, ConsModule_ON: True, CUT_OUTEXT: 2000

### prob_todo

The PDE of Poisson-Boltzmann(Helmholtz) 2D equation is given by:
\begin{equation}
-\Delta u + k^2 u = f(x,y)
\end{equation}
The function f(x) is defined as:
\begin{equation}
f(x)=A\cdot\left(\sum_{i}\mu_{i}^{2}+x_{i}^{2}\right)\sin(\mu_{1}\pi x_{1})\sin(\mu_{2}\pi x_{2})
\end{equation}
The parameter values are:
\begin{equation}
\mu_1 = 1,\quad \mu_2 = 4,\quad k = 8,\quad A = 10
\end{equation}

The domain is \(\Omega = [- 1,1]^{2}\) with several circles removed.
The circles \(\Omega_{\text{circle}} = \bigcup_{i = 1}^{4} R_i\) are defined as:
\begin{cases}
R_{1}=\left\{(x,y):(x - 0.5)^{2}+(y - 0.5)^{2}\leq0.2^{2}\right\}\\
R_{2}=\left\{(x,y):(x - 0.4)^{2}+(y + 0.4)^{2}\leq0.4^{2}\right\}\\
R_{3}=\left\{(x,y):(x + 0.2)^{2}+(y + 0.7)^{2}\leq0.1^{2}\right\}\\
R_{4}=\left\{(x,y):(x + 0.6)^{2}+(y - 0.5)^{2}\leq0.3^{2}\right\}\\
\end{cases}
The boundary conditions are:
\begin{cases}
u = 0.2, & x\in\partial\Omega_{\text{rec}}\\
u = 1, & x\in\partial\Omega_{\text{circle}}
\end{cases}

Implement a stable and efficient method to solve this problem.
Define visc = 1.0 / Re at the GLOBAL scope (outside main), and strictly use visc instead of nu throughout the code to avoid confusion with u.
Ensure all constants are passed explicitly to functions or defined globally to prevent NameError.
Plot the contour of u in one figure using 'RdBu_r' colormap, and mark the circles in the plot.
Just save figs do not use plt.show() in the code.

[HINTS]:
Address the geometric challenges of curved boundaries using appropriate techniques (e.g., Immersed Boundary Method).
Print necessary solving information to facilitate the reliability check of the solution. 


### expanded_prob
We must solve a 2D linear elliptic PDE with a positive mass term: -Δu + k^2 u = f(x,y) on the bounded rectangular domain Ω = [-1,1]^2 with four interior circular holes removed. Dirichlet boundary conditions are prescribed: u = 0.2 on the outer rectangle boundary and u = 1 on each circle boundary. f(x,y) is smooth and separable (product of sines) with known parameters μ1=1, μ2=4, k=8, A=10. The practical implementation must be stable and efficient, define and use a global visc = 1.0/Re (outside main) and pass constants explicitly to avoid NameError, and produce a saved contour plot using colormap 'RdBu_r' with the circles highlighted.

Primary mathematical and numerical challenges:
- Geometry: multiple curved internal boundaries (circles) create a multiply-connected domain. Accurately imposing Dirichlet data on these curves (and resolving curvature) is nontrivial for simple Cartesian discretizations.
- Boundary enforcement: internal Dirichlet data on holes requires either a mesh that conforms to the circles (so boundary conditions can be applied strongly) or an embedded/immersed technique (ghost points, cut elements, Nitsche/CutFEM) to impose BCs accurately and stably.
- Resolution and accuracy: small circles (e.g., radius 0.1, 0.2) demand local mesh resolution to control geometric error; insufficient resolution near curved boundaries degrades global error and can pollute the solution interior.
- Conditioning and solver choice: the operator L = -Δ + k^2 I is elliptic and symmetric positive definite for k^2>0, but discretizations can still be ill-conditioned for fine meshes. Efficient linear solvers (multigrid, algebraic multigrid (AMG), or preconditioned Krylov) are essential for performance.
- Implementation robustness: constants must be passed explicitly and visc defined globally to avoid runtime NameError; printing solver diagnostics (residual history, iterations, time, estimated condition number) is needed for confidence.
- Visualization: overlaying circle boundaries on the contour plot and saving the figure (without plt.show()) are required.

These considerations determine trade-offs between ease-of-implementation, geometric fidelity, accuracy, and computational cost. Below are two rigorous, alternative solution strategies (one conforming finite-element approach and one Cartesian immersed/ghost technique).

### solution_plans
	Current Stage [A/2]
solu_name='Conforming FEM with AMG' content="Solu Name: Conforming Unstructured Finite Element Method (FEM) with Adaptive Refinement and AMG\n\nGoverning idea:\nUse an unstructured triangular mesh that exactly resolves the rectangle and the four circular holes (mesh conforms to all boundaries). Solve the variational form of -Δu + k^2 u = f with strong Dirichlet enforcement on outer and inner boundaries using a standard FEM library (FEniCS, Firedrake, or a custom FEM using gmsh + PETSc). Use adaptive or graded refinement near circles to control geometric and interpolation error. Solve the resulting SPD linear system with a robust preconditioned Krylov method (Conjugate Gradient) preconditioned by AMG.\n\nAlgorithmic steps (numbered):\n1. Problem setup and constants: declare global visc = 1.0 / Re at module scope; pass in k, μ1, μ2, A, visc explicitly to all solver functions. Compute f(x,y) analytically.\n2. Geometry & mesh: define rectangle [-1,1]^2 and subtract circular holes R1..R4. Use a mesh generator (gmsh) to create a boundary-fitted triangular mesh. Specify finer target size near circle boundaries (mesh size function or physical fields) and optionally use boundary curvature-based sizing.\n3. Function spaces & variational form: choose continuous Lagrange finite elements (P2 for improved geometry/solution accuracy; or P1 if speed needed). Write weak form: ∫(∇u·∇v + k^2 u v) dx = ∫ f v dx. Impose Dirichlet BCs strongly on outer and circle boundaries.\n4. Assembly: assemble stiffness matrix and load vector using quadrature consistent with chosen polynomial degree. For integrals near curved boundaries, use the conforming mesh so no special treatment required.\n5. Linear solver & preconditioning: use Conjugate Gradient (CG) since matrix symmetric positive definite. Precondition with AMG (PETSc+Hypre/BoomerAMG or pyamg). Provide solver monitor callbacks to print residual reduction, iteration counts, and timing. Optionally compute an estimate of condition number or power method diagnostic.\n6. Adaptive refinement (optional): compute a posteriori error indicator (residual-based or jump indicators), refine mesh locally around circle boundaries and re-solve until error tolerance or limit reached.\n7. Postprocessing: evaluate u on a uniform Cartesian grid for plotting (interpolate FEM solution onto grid), create contour plot with cmap='RdBu_r', overlay circle boundaries (parametric plotting of circumference). Save figure to file (plt.savefig) and do not call plt.show().\n8. Diagnostics: print mesh statistics (number of nodes/els), DOFs, solver iterations, final residual, solve time, and memory footprint.\n\nStability / accuracy / complexity / efficiency limitations:\n- Accuracy: conforming FEM with P2 elements gives good geometric fidelity and solution accuracy; errors scale with element size h^p (p polynomial degree). Geometric error is negligible if mesh resolves curvature.\n- Complexity: assembling cost O(N) and solving with AMG-preconditioned CG typically near-optimal (O(N) or O(N log N)), where N is DOFs. AMG setup cost can be significant for very large problems.\n- Stability: discrete operator is SPD; CG+AMG is stable. However, if k is large relative to mesh resolution, mass term dominates and solution gradients may require finer mesh.\n- Implementation effort: requires mesh generation tools (gmsh) and FEM library; adaptive refinement increases coding complexity.\n- Memory: unstructured meshes and AMG can be memory-intensive for very fine discretizations.\n\n"

	Current Stage [B/2]
solu_name='Cartesian Immersed/Ghost Method with Multigrid' content="Solu Name: Cartesian Grid with Immersed Boundary / Ghost-Fluid (Sharp Interface) and Geometric Multigrid\n\nGoverning idea:\nUse a uniform Cartesian grid on the bounding rectangle and embed the circular holes inside. Enforce Dirichlet BCs on internal circles via a sharp-interface ghost-point or ghost-fluid approach (second-order accurate). Discretize -Δu + k^2 u with standard central differences away from the interface and modify finite difference stencils that cross the circle boundary using either a ghost-value extrapolation (to enforce u=1 on the circle) or a short-stencil interpolation that incorporates the boundary position. Solve the resulting sparse linear system with a geometric multigrid tailored to Cartesian grids (or use AMG) for efficiency.\n\nAlgorithmic steps (numbered):\n1. Constants and globals: define visc = 1.0 / Re globally; pass k, μ1, μ2, A, visc, circle parameters (centers, radii) to all functions explicitly.\n2. Cartesian discretization: choose grid resolution Nx×Ny so that smallest circle is resolved by ≥8–12 points across radius (refine until convergence). Construct cell-centered (or node-centered) arrays and index mapping.\n3. Regular/irregular classification: for each grid point determine whether it lies inside domain (outside all circles) or inside a hole. Identify “regular” points whose 5-point Laplacian stencil lies entirely inside domain, and “irregular/near-interface” points whose stencil crosses a circle boundary.\n4. Interface treatment (ghost/ghost-fluid): for irregular points, compute intersection(s) of stencil segments with the circle boundary. Impose Dirichlet u=1 at the exact intersection by constructing ghost values via linear/quadratic extrapolation or by solving a small local interpolation problem that enforces the boundary condition strongly. This yields modified finite-difference equations that remain second-order accurate in L∞ provided the interface is smooth and the interpolation is consistent.\n   - Option: use the Immersed Interface Method (IIM) or CutFEM/Nitsche-like discrete weak enforcement for higher robustness (requires slightly more complex linear algebra).\n5. Assemble sparse matrix: fill sparse matrix entries for all interior (regular + modified irregular) equations; enforce outer rectangle Dirichlet nodes by replacing rows with identity and RHS=0.2.\n6. Solver: use a geometric multigrid (GMG) tailored to the Cartesian layout with modified coarse-grid operators near interfaces (or use AMG if GMG incorporation of interface is complex). Use preconditioned conjugate gradient (for SPD) or direct solver for moderate sizes. Print residual history, iteration counts, and timing.\n7. Convergence checks: grid-convergence study (refine grid and compute L2 differences) to verify second-order behavior and ensure that ghost interpolation preserves accuracy near circles.\n8. Plotting: interpolate/reshape solution to 2D arrays; contour with cmap='RdBu_r'; overlay circles by plotting their parametric curves; save figure using plt.savefig without plt.show().\n9. Diagnostics: print numbers of grid points, irregular points count, solver iterations, residual norm, and estimated error from grid-refinement.\n\nStability / accuracy / complexity / efficiency limitations:\n- Accuracy: properly implemented ghost/immersed-interface methods achieve second-order accuracy in smooth regions. However, if ghost extrapolation is low order or interface geometry is under-resolved, local accuracy near circles degrades and pollutes global error.\n- Stability: naive ghost-value extrapolation may cause stencil imbalance and loss of SPD property; use symmetric modification or IIM/CutFEM variants to preserve SPD and solver compatibility.\n- Complexity: Cartesian codes are simpler to implement and exploit fast geometric multigrid (optimal complexity O(N)). However, handling many irregular stencils and robust coarse-grid transfer near interfaces complicates a fully optimal GMG implementation.\n- Implementation effort: more work to implement robust ghost/interpolation and correct coarse-grid operators; CutFEM/Nitsche provides a principled alternative but is more advanced.\n- Memory/time: Cartesian-grid + multigrid is memory- and time-efficient for large problems, but achieving strict second-order accuracy near many small holes requires sufficiently fine grids, increasing N.\n\n"



### technical_spec
	Current Stage [A/3]
We implement a Cartesian finite-difference solver with immersed/ghost treatment for interior circular Dirichlet holes. A uniform node-centered grid covers the rectangle [-1,1]^2; nodes that fall inside any circle (holes) are excluded as unknowns. For 5-point Laplacian stencils that cross a circle boundary, we compute the exact intersection parameter along the edge and eliminate the ghost value by enforcing the Dirichlet condition at that intersection. The code assembles a sparse linear system once (no time loop), solves it using conjugate gradients (scipy.sparse.linalg.cg) while tracking residual history, and interpolates the solution back to the full grid for plotting. The plot uses 'RdBu_r' colormap and overlays the circle boundaries. Diagnostics printed: grid size, DOFs, number of irregular (ghost) edges handled, solver iterations, residual history and solve time. All constants (k, mu1, mu2, A, visc, circle params) are passed explicitly to functions.

	Current Stage [B/3]
Cartesian finite-difference solver with immersed/ghost treatment for multiple circular Dirichlet holes. Key steps: (1) build node classification distinguishing hole-interior, circle-boundary (|dist-r|<=eps) and domain nodes; (2) create index map for unknowns (domain + boundary-dirichlet nodes kept as known rows) and assemble sparse matrix for -Delta u + k^2 u with ghost-edge elimination for edges that cross hole interiors; (3) symmetrize the assembled matrix and attempt to build an ILU preconditioner; (4) solve with scipy.sparse.linalg.cg using a robust wrapper that tries different keyword combinations for compatibility across SciPy versions; fallback to spsolve if iterative solves fail; the callback for residual tracking is made robust to different callback signatures; (5) interpolate solution to full grid (assign circle Dirichlet values inside holes for plotting) and save a labeled contour plot with the circle boundaries overlaid. Diagnostics printed: grid info, DOFs, irregular edges handled, solver method/fallback info, residual history length and final residual norm, solve time, and sample solution values. All physical constants are passed explicitly; visc is defined at global scope and passed around but not applied to the PDE (printed for clarity).

	Current Stage [C/3]
This module solves the 2D Poisson-Helmholtz problem on a square domain with multiple circular Dirichlet holes using a Cartesian finite-difference discretization and an immersed/ghost-edge elimination when edges cross hole interiors. Key steps:
- Grid and node classification (domain, outer Dirichlet, circle-Dirichlet, hole-interior). Classification tolerance is configurable; counts are printed for diagnostics.
- Sparse matrix assembly of -Delta u + k^2 u using a 5-point stencil. For neighbor edges that cross into hole interior a ghost-elimination algebraic substitution is applied (linear interpolation to the Dirichlet boundary) giving consistent modifications to diagonal and RHS. Per-circle counts of irregular edges are recorded.
- Detailed matrix diagnostics: min/max/mean of stored entries, symmetry measure, zero-diagonal count, max row-sum, and inexpensive eigen estimates via eigsh when available.
- Solver strategy: first attempt a direct solve (spsolve) and report residual and timing to determine if assembly is correct. Then attempt iterative solve (CG) with ILU preconditioner and robust fallbacks; fall back to GMRES/BiCGSTAB if matrix is non-symmetric/indefinite.
- Unit tests: lightweight checks for edge-circle intersection and ghost-coefficient algebra to exercise irregular elimination calculations.
- Interpolation to full grid and a saved contour plot with circle boundaries overlaid.
All physical constants (including visc defined globally) are passed or printed; no matrices are rebuilt in time loops. Extensive printed diagnostics help diagnose conditioning and solver behavior.



### python_codes

#### Script block1:

```python
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
import time

# Global visc definition as requested (use visc instead of nu)
Re = 1.0
visc = 1.0 / Re

# --- Problem constants (passed explicitly to functions) ---
k = 8.0
mu1 = 1.0
mu2 = 4.0
A = 10.0

# Circle definitions: list of (cx, cy, r)
circles = [
    (0.5, 0.5, 0.2),
    (0.4, -0.4, 0.4),
    (-0.2, -0.7, 0.1),
    (-0.6, 0.5, 0.3),
]

# Boundary Dirichlet values
outer_bc = 0.2
circle_bc = 1.0

# Grid resolution (node-centered)
Nx = 201  # choose sufficiently fine to resolve smallest circle
Ny = 201
xmin, xmax = -1.0, 1.0
ymin, ymax = -1.0, 1.0

# --- Helper functions ---

def f_rhs(x, y, mu1, mu2, A):
    """Right-hand side f(x,y) as given in problem statement.
    We interpret f = A * ( (mu1^2 + x^2) + (mu2^2 + y^2) ) * sin(mu1*pi*x) * sin(mu2*pi*y)
    """
    return A * ((mu1 ** 2 + x ** 2) + (mu2 ** 2 + y ** 2)) * np.sin(mu1 * np.pi * x) * np.sin(mu2 * np.pi * y)


def point_in_any_circle(x, y, circles):
    for (cx, cy, r) in circles:
        if (x - cx) ** 2 + (y - cy) ** 2 <= r ** 2:
            return True
    return False


def edge_circle_intersection(p, q, circle):
    """Return t in (0,1) where segment p -> q intersects circle, or None.
    Solve |p + t d - c|^2 = r^2 for t.
    """
    cx, cy, r = circle
    px, py = p
    qx, qy = q
    dx = qx - px
    dy = qy - py
    a = dx * dx + dy * dy
    b = 2.0 * (dx * (px - cx) + dy * (py - cy))
    c0 = (px - cx) ** 2 + (py - cy) ** 2 - r * r
    disc = b * b - 4 * a * c0
    if disc < 0:
        return None
    sqrt_disc = np.sqrt(disc)
    t1 = (-b - sqrt_disc) / (2 * a)
    t2 = (-b + sqrt_disc) / (2 * a)
    for t in (t1, t2):
        if 1e-12 < t < 1 - 1e-12:
            return t
    return None


def assemble_system(Nx, Ny, xmin, xmax, ymin, ymax, circles, k, mu1, mu2, A, outer_bc, circle_bc, visc):
    hx = (xmax - xmin) / (Nx - 1)
    hy = (ymax - ymin) / (Ny - 1)
    assert abs(hx - hy) < 1e-12, "Non-uniform spacing not supported in this implementation"
    h = hx

    xs = np.linspace(xmin, xmax, Nx)
    ys = np.linspace(ymin, ymax, Ny)

    # Node flags and indexing: include nodes that are NOT inside holes
    inside_hole = np.zeros((Ny, Nx), dtype=bool)
    for j in range(Ny):
        for i in range(Nx):
            inside_hole[j, i] = point_in_any_circle(xs[i], ys[j], circles)

    # Map (i,j) to unknown index if not inside hole
    idx_map = -np.ones((Ny, Nx), dtype=int)
    counter = 0
    for j in range(Ny):
        for i in range(Nx):
            if not inside_hole[j, i]:
                idx_map[j, i] = counter
                counter += 1
    N = counter

    # Preallocate lists for COO assembly
    rows = []
    cols = []
    vals = []
    rhs = np.zeros(N, dtype=float)

    irregular_count = 0

    # Precompute f at nodes (for nodes not in holes)
    fvals = np.zeros((Ny, Nx), dtype=float)
    for j in range(Ny):
        for i in range(Nx):
            fvals[j, i] = f_rhs(xs[i], ys[j], mu1, mu2, A)

    # Stencil neighbor offsets
    neigh = [(-1, 0), (1, 0), (0, -1), (0, 1)]

    for j in range(Ny):
        for i in range(Nx):
            if inside_hole[j, i]:
                continue
            idx = idx_map[j, i]
            x = xs[i]
            y = ys[j]

            # Outer rectangle boundary nodes get Dirichlet = outer_bc
            if i == 0 or i == Nx - 1 or j == 0 or j == Ny - 1:
                rows.append(idx)
                cols.append(idx)
                vals.append(1.0)
                rhs[idx] = outer_bc
                continue

            # Interior domain node: assemble -Delta u + k^2 u = f
            diag = 4.0 / (h * h) + k * k
            rhs_val = fvals[j, i]

            for d in neigh:
                ni = i + d[0]
                nj = j + d[1]
                # neighbor should be within grid; since this is interior node, it is
                if inside_hole[nj, ni]:
                    # edge crosses into a hole; find which circle
                    p = (x, y)
                    q = (xs[ni], ys[nj])
                    found = False
                    for circ in circles:
                        t = edge_circle_intersection(p, q, circ)
                        if t is not None:
                            # ghost value g = (u_b - (1 - t) * u0) / t
                            # neighbor contributes -1/h^2 * g
                            # substitute: contributes -1/(t h^2) * u_b + (1 - t)/(t h^2) * u0
                            diag += (1.0 - t) / (t * h * h)
                            rhs_val += (1.0 / (t * h * h)) * circle_bc
                            irregular_count += 1
                            found = True
                            break
                    if not found:
                        # Should not happen; but handle robustly by treating neighbor as Dirichlet=1
                        diag += 1.0 / (h * h)
                        rhs_val += (1.0 / (h * h)) * circle_bc
                        irregular_count += 1
                else:
                    # neighbor is part of domain (could be outer boundary or interior)
                    nidx = idx_map[nj, ni]
                    # if neighbor is outer boundary interior or interior node
                    # If neighbor is outer boundary Dirichlet, its row will be identity but it's included here as known.
                    if (ni == 0 or ni == Nx - 1 or nj == 0 or nj == Ny - 1):
                        # known outer boundary value
                        rhs_val += (1.0 / (h * h)) * outer_bc
                    else:
                        # Unknown neighbor contributes off-diagonal
                        rows.append(idx)
                        cols.append(nidx)
                        vals.append(-1.0 / (h * h))

            # insert diagonal and rhs
            rows.append(idx)
            cols.append(idx)
            vals.append(diag)
            rhs[idx] = rhs_val

    A = sp.coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsr()
    return A, rhs, idx_map, xs, ys, inside_hole, irregular_count


def solve_system(A, rhs, tol=1e-8, maxiter=2000):
    # CG with residual history via callback
    residuals = []

    def make_callback(A, rhs, residuals):
        def cb(xk):
            rk = rhs - A.dot(xk)
            residuals.append(np.linalg.norm(rk))
        return cb

    callback = make_callback(A, rhs, residuals)
    start = time.time()
    x, info = spla.cg(A, rhs, tol=tol, maxiter=maxiter, callback=callback)
    elapsed = time.time() - start
    if info == 0:
        status = 'converged'
    elif info > 0:
        status = f'no convergence after {info} iterations'
    else:
        status = f'illegal input or breakdown (info={info})'
    return x, residuals, elapsed, status


def interpolate_solution_to_grid(x_vec, idx_map, xs, ys, inside_hole, outer_bc, circle_bc, circles):
    Ny, Nx = idx_map.shape
    sol = np.full((Ny, Nx), np.nan, dtype=float)
    # fill domain unknowns
    for j in range(Ny):
        for i in range(Nx):
            if inside_hole[j, i]:
                sol[j, i] = np.nan
            else:
                sol[j, i] = x_vec[idx_map[j, i]]
    # fill hole boundaries by Dirichlet value for plotting mask
    # For visualization, set inside holes to circle_bc
    for j in range(Ny):
        for i in range(Nx):
            if inside_hole[j, i]:
                sol[j, i] = circle_bc
    return sol


def plot_solution(sol, xs, ys, circles, fname='solution_contour.png'):
    X, Y = np.meshgrid(xs, ys)
    plt.figure(figsize=(8, 6))
    levels = 60
    cf = plt.contourf(X, Y, sol, levels=levels, cmap='RdBu_r')
    plt.colorbar(cf, label='u')
    # overlay circle boundaries
    theta = np.linspace(0, 2 * np.pi, 400)
    for (cx, cy, r) in circles:
        plt.plot(cx + r * np.cos(theta), cy + r * np.sin(theta), 'k-', linewidth=1.5)
    plt.title('Solution u (contour) with circle boundaries')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.close()


# --- Main execution ---
if __name__ == '__main__':
    print('Assembling system...')
    A, rhs, idx_map, xs, ys, inside_hole, irregular_count = assemble_system(
        Nx, Ny, xmin, xmax, ymin, ymax, circles, k, mu1, mu2, A, outer_bc, circle_bc, visc)

    N = A.shape[0]
    nnz = A.nnz
    print(f'Grid: Nx={Nx}, Ny={Ny}, total nodes={Nx * Ny}')
    print(f'Unknowns (nodes not inside holes) = {N}')
    print(f'Sparse matrix: shape={A.shape}, nnz={nnz}')
    print(f'Number of irregular (ghost) edge treatments = {irregular_count}')
    print('Starting solver (Conjugate Gradient)...')

    x, residuals, elapsed, status = solve_system(A, rhs, tol=1e-8, maxiter=5000)
    print(f'Solver status: {status}')
    print(f'Solve time: {elapsed:.3f} s')
    if len(residuals) > 0:
        print(f'Initial residual norm (first callback) = {residuals[0]:.3e}, final tracked residual = {residuals[-1]:.3e}')
        print(f'Number of tracked iterations = {len(residuals)}')
    else:
        print('No residual history recorded (callback not invoked)')

    # Compute final residual norm explicitly
    final_res = np.linalg.norm(rhs - A.dot(x))
    print(f'Explicit final residual norm ||b - A x||_2 = {final_res:.3e}')

    # Interpolate back to full grid and plot
    sol_grid = interpolate_solution_to_grid(x, idx_map, xs, ys, inside_hole, outer_bc, circle_bc, circles)
    plot_solution(sol_grid, xs, ys, circles, fname='solution_contour.png')
    print('Saved contour plot to solution_contour.png')

    # Print a few sample values for verification
    sample_points = [(-0.9, -0.9), (0.0, 0.0), (0.5, 0.5), (0.8, -0.8)]
    print('Sample solution values (x,y) -> u :')
    for (sx, sy) in sample_points:
        # find closest grid node
        i = int(round((sx - xmin) / (xmax - xmin) * (Nx - 1)))
        j = int(round((sy - ymin) / (ymax - ymin) * (Ny - 1)))
        val = sol_grid[j, i]
        print(f'  ({sx:.2f}, {sy:.2f}) -> u = {val:.6f}')

    # Summary including visc (to demonstrate usage of global visc)
    print(f'Viscosity (global visc) = {visc} (Re= {Re})')
```


#### Script block2:

```python
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse.linalg import LinearOperator
import matplotlib.pyplot as plt
import time

# Global visc definition as requested (use visc instead of nu)
Re = 1.0
visc = 1.0 / Re

# --- Problem constants (passed explicitly to functions) ---
k = 8.0
mu1 = 1.0
mu2 = 4.0
A_amp = 10.0  # amplitude (renamed from A to avoid name clash with matrix A)

# Circle definitions: list of (cx, cy, r)
circles = [
    (0.5, 0.5, 0.2),
    (0.4, -0.4, 0.4),
    (-0.2, -0.7, 0.1),
    (-0.6, 0.5, 0.3),
]

# Boundary Dirichlet values
outer_bc = 0.2
circle_bc = 1.0

# Grid resolution (node-centered)
Nx = 201  # choose sufficiently fine to resolve smallest circle
Ny = 201
xmin, xmax = -1.0, 1.0
ymin, ymax = -1.0, 1.0


# --- Helper functions ---

def f_rhs(x, y, mu1, mu2, A_amp):
    """Right-hand side f(x,y) used for forcing.
    Interpreted as: f = A * (mu1^2 + mu2^2 + x^2 + y^2) * sin(mu1*pi*x) * sin(mu2*pi*y)
    """
    return A_amp * ((mu1 ** 2 + mu2 ** 2) + x ** 2 + y ** 2) * np.sin(mu1 * np.pi * x) * np.sin(mu2 * np.pi * y)


def edge_circle_intersection(p, q, circle):
    """Return t in (0,1) where segment p->q intersects circle, or None.
    Solve |p + t d - c|^2 = r^2 for t.
    """
    cx, cy, r = circle
    px, py = p
    qx, qy = q
    dx = qx - px
    dy = qy - py
    a = dx * dx + dy * dy
    if a == 0:
        return None
    b = 2.0 * (dx * (px - cx) + dy * (py - cy))
    c0 = (px - cx) ** 2 + (py - cy) ** 2 - r * r
    disc = b * b - 4 * a * c0
    if disc < 0:
        return None
    sqrt_disc = np.sqrt(disc)
    t1 = (-b - sqrt_disc) / (2 * a)
    t2 = (-b + sqrt_disc) / (2 * a)
    for t in (t1, t2):
        if 1e-12 < t < 1 - 1e-12:
            return t
    return None


def classify_nodes(xs, ys, circles, h):
    """Classify each grid node relative to the circles.
    Returns an integer array node_type with values:
      0 -> hole interior (exclude)
      1 -> domain node (unknown)
      2 -> outer rectangle boundary (Dirichlet known)
      3 -> circle boundary node (Dirichlet known)
    The tolerance eps for boundary classification is eps = 0.5*h
    """
    Ny = len(ys)
    Nx = len(xs)
    node_type = np.ones((Ny, Nx), dtype=np.int8)  # default domain
    eps = 0.5 * h
    for j in range(Ny):
        for i in range(Nx):
            x = xs[i]
            y = ys[j]
            # outer rectangle boundary
            if i == 0 or i == Nx - 1 or j == 0 or j == Ny - 1:
                node_type[j, i] = 2
                continue
            # compute signed distance to nearest circle: signed = dist - r
            min_signed = np.inf
            for (cx, cy, r) in circles:
                dist = np.hypot(x - cx, y - cy)
                signed = dist - r
                if signed < min_signed:
                    min_signed = signed
            if min_signed < -eps:
                node_type[j, i] = 0  # strictly inside hole
            elif abs(min_signed) <= eps:
                node_type[j, i] = 3  # circle boundary node (Dirichlet enforced strongly)
            else:
                node_type[j, i] = 1  # domain unknown
    return node_type


def assemble_system(Nx, Ny, xmin, xmax, ymin, ymax, circles, k, mu1, mu2, A_amp, outer_bc, circle_bc, visc):
    """Assemble sparse linear system once. Returns matrix A (CSR), rhs, idx_map, xs, ys, node_type, irregular_count.
    Note: visc is passed but not used in PDE (printed later for diagnostics).
    """
    hx = (xmax - xmin) / (Nx - 1)
    hy = (ymax - ymin) / (Ny - 1)
    assert abs(hx - hy) < 1e-12, "Non-uniform spacing not supported in this implementation"
    h = hx

    xs = np.linspace(xmin, xmax, Nx)
    ys = np.linspace(ymin, ymax, Ny)

    # classify nodes
    node_type = classify_nodes(xs, ys, circles, h)

    # Map (i,j) to unknown index for all nodes that are NOT hole interior (types 0)
    idx_map = -np.ones((Ny, Nx), dtype=int)
    counter = 0
    for j in range(Ny):
        for i in range(Nx):
            if node_type[j, i] != 0:
                idx_map[j, i] = counter
                counter += 1
    N = counter

    # Preallocate lists for COO assembly
    rows = []
    cols = []
    vals = []
    rhs = np.zeros(N, dtype=float)

    irregular_count = 0

    # Precompute f at nodes (for nodes that are domain unknowns)
    fvals = np.zeros((Ny, Nx), dtype=float)
    for j in range(Ny):
        for i in range(Nx):
            fvals[j, i] = f_rhs(xs[i], ys[j], mu1, mu2, A_amp)

    # Stencil neighbor offsets (i,j offsets correspond to x,y grid indexing)
    neigh = [(-1, 0), (1, 0), (0, -1), (0, 1)]

    for j in range(Ny):
        for i in range(Nx):
            ntype = node_type[j, i]
            if ntype == 0:
                continue  # hole interior not represented
            idx = idx_map[j, i]
            x = xs[i]
            y = ys[j]

            # Outer rectangle boundary nodes get Dirichlet = outer_bc
            if ntype == 2:
                rows.append(idx); cols.append(idx); vals.append(1.0)
                rhs[idx] = outer_bc
                continue

            # Circle boundary nodes get Dirichlet = circle_bc (enforce strongly)
            if ntype == 3:
                rows.append(idx); cols.append(idx); vals.append(1.0)
                rhs[idx] = circle_bc
                continue

            # Interior domain node: assemble -Delta u + k^2 u = f
            diag = 4.0 / (h * h) + k * k
            rhs_val = fvals[j, i]

            for d in neigh:
                ni = i + d[0]
                nj = j + d[1]
                # neighbor should be within grid; since this is interior node, it is
                if node_type[nj, ni] == 0:
                    # edge crosses into a hole interior; find which circle intersects
                    p = (x, y)
                    q = (xs[ni], ys[nj])
                    found = False
                    for circ in circles:
                        t = edge_circle_intersection(p, q, circ)
                        if t is not None:
                            # eliminate ghost by enforcing Dirichlet at intersection
                            # contribution modifies diagonal and rhs
                            diag += (1.0 - t) / (t * h * h)
                            rhs_val += (1.0 / (t * h * h)) * circle_bc
                            irregular_count += 1
                            found = True
                            break
                    if not found:
                        # Robust fallback: treat as if neighbor is known circle boundary
                        diag += 1.0 / (h * h)
                        rhs_val += (1.0 / (h * h)) * circle_bc
                        irregular_count += 1
                else:
                    # neighbor is represented in system (could be outer boundary or circle boundary or domain unknown)
                    nidx = idx_map[nj, ni]
                    # If neighbor is Dirichlet known (outer or circle), add contribution to RHS
                    if node_type[nj, ni] in (2, 3):
                        known_val = outer_bc if node_type[nj, ni] == 2 else circle_bc
                        rhs_val += (1.0 / (h * h)) * known_val
                    else:
                        # Unknown neighbor contributes off-diagonal
                        rows.append(idx)
                        cols.append(nidx)
                        vals.append(-1.0 / (h * h))

            # insert diagonal and rhs
            rows.append(idx); cols.append(idx); vals.append(diag)
            rhs[idx] = rhs_val

    A_mat = sp.coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsr()

    return A_mat, rhs, idx_map, xs, ys, node_type, irregular_count


def make_callback_factory(A_mat, rhs, residuals):
    """Create a robust callback that handles different scipy versions:
    callback may receive solution vector xk, residual vector rk, or iteration count int.
    We try to detect and compute norm accordingly.
    """
    def callback(arg):
        try:
            # If arg is an integer (iteration count), ignore
            if isinstance(arg, (int, np.integer)):
                return
            arr = np.asarray(arg)
            # If arr length matches rhs, it could be solution or residual
            if arr.shape[0] == rhs.shape[0]:
                # Prefer interpreting arg as solution vector: try computing residual rhs - A x
                try:
                    rk = rhs - A_mat.dot(arr)
                    res = np.linalg.norm(rk)
                except Exception:
                    # fallback: treat arg as residual vector
                    res = np.linalg.norm(arr)
                residuals.append(res)
        except Exception:
            # Do not let callback exceptions break the solver
            return

    return callback


def try_build_preconditioner(A_mat):
    """Attempt to build an ILU preconditioner; return LinearOperator or None.
    Be permissive: if spilu fails, return None.
    """
    try:
        # spilu requires csc or csr; use csc for factorization
        A_csc = A_mat.tocsc()
        # Adjust drop tolerance and fill factor to be robust
        ilu = spla.spilu(A_csc, drop_tol=1e-4, fill_factor=10)
        Mx = lambda x: ilu.solve(x)
        M = LinearOperator(dtype=A_mat.dtype, shape=A_mat.shape, matvec=Mx)
        print('Built ILU preconditioner (spilu)')
        return M
    except Exception as e:
        print(f'Could not build ILU preconditioner: {e}')
        return None


def run_cg_robust(A_mat, rhs, tol=1e-8, maxiter=2000, callback=None, M=None):
    """Run scipy.sparse.linalg.cg with a sequence of compatibility fallbacks.
    Returns (x, info, used_method_str).
    """
    attempts = []
    # Try with tol
    attempts.append({'tol': tol, 'maxiter': maxiter, 'callback': callback, 'M': M})
    # Try with atol (some scipy use atol instead of tol)
    attempts.append({'atol': tol, 'maxiter': maxiter, 'callback': callback, 'M': M})
    # Try with only maxiter and callback and M
    attempts.append({'maxiter': maxiter, 'callback': callback, 'M': M})
    # Try with maxiter and callback (no M)
    attempts.append({'maxiter': maxiter, 'callback': callback})

    last_exc = None
    for kwargs in attempts:
        try:
            # Filter None values and ensure we don't pass unsupported M=None explicitly
            kw = {k: v for k, v in kwargs.items() if v is not None}
            x, info = spla.cg(A_mat, rhs, **kw)
            return x, info, f'cg with kwargs={list(kw.keys())}'
        except TypeError as e:
            last_exc = e
            # try next combination
            continue
        except Exception as e:
            # other failures (e.g. breakdown), propagate up to allow fallback
            last_exc = e
            break

    # If we reach here, iterative solver attempts failed in compatibility or execution
    print(f'CG iterative attempts failed (last exception: {last_exc}). Falling back to direct solver spsolve.')
    x = spla.spsolve(A_mat.tocsc(), rhs)
    return x, 0, 'spsolve (fallback)'


def solve_system(A_mat, rhs, tol=1e-8, maxiter=2000):
    residuals = []
    callback = make_callback_factory(A_mat, rhs, residuals)

    # Try to symmetrize small asymmetries introduced by elimination
    try:
        # Only do this in sparse form by averaging with its transpose
        A_sym = 0.5 * (A_mat + A_mat.T)
    except Exception:
        A_sym = A_mat

    # Attempt a preconditioner for robustness
    M = try_build_preconditioner(A_sym)

    start = time.time()
    x, info, method = run_cg_robust(A_sym, rhs, tol=tol, maxiter=maxiter, callback=callback, M=M)
    elapsed = time.time() - start

    if isinstance(info, int):
        if info == 0:
            status = 'converged'
        elif info > 0:
            status = f'no convergence after {info} iterations'
        else:
            status = f'illegal input or breakdown (info={info})'
    else:
        status = f'info={info}'

    return x, residuals, elapsed, status, method


def interpolate_solution_to_grid(x_vec, idx_map, xs, ys, node_type, outer_bc, circle_bc):
    Ny, Nx = idx_map.shape
    sol = np.full((Ny, Nx), np.nan, dtype=float)
    # fill nodes that are represented in the system
    for j in range(Ny):
        for i in range(Nx):
            if node_type[j, i] == 0:
                # strictly inside hole: for plotting set to circle_bc so holes are visible
                sol[j, i] = circle_bc
            else:
                sol[j, i] = x_vec[idx_map[j, i]]
    return sol


def plot_solution(sol, xs, ys, circles, fname='solution_contour.png'):
    X, Y = np.meshgrid(xs, ys)
    plt.figure(figsize=(8, 6))
    levels = 60
    cf = plt.contourf(X, Y, sol, levels=levels, cmap='RdBu_r')
    cbar = plt.colorbar(cf, label='u')
    # overlay circle boundaries
    theta = np.linspace(0, 2 * np.pi, 400)
    for (cx, cy, r) in circles:
        plt.plot(cx + r * np.cos(theta), cy + r * np.sin(theta), 'k-', linewidth=1.5)
    plt.title('Solution u (contour) with circle boundaries')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.close()


# --- Main execution ---
if __name__ == '__main__':
    print('Assembling system...')
    A_mat, rhs, idx_map, xs, ys, node_type, irregular_count = assemble_system(
        Nx, Ny, xmin, xmax, ymin, ymax, circles, k, mu1, mu2, A_amp, outer_bc, circle_bc, visc)

    N = A_mat.shape[0]
    nnz = A_mat.nnz
    print(f'Grid: Nx={Nx}, Ny={Ny}, total nodes={Nx * Ny}')
    print(f'Unknowns (represented nodes not strictly inside holes) = {N}')
    print(f'Sparse matrix: shape={A_mat.shape}, nnz={nnz}')
    print(f'Number of irregular (ghost) edge treatments = {irregular_count}')
    print('Starting solver (robust Conjugate Gradient with fallbacks)...')

    x, residuals, elapsed, status, method = solve_system(A_mat, rhs, tol=1e-8, maxiter=5000)
    print(f'Solver method/result: {method}')
    print(f'Solver status: {status}')
    print(f'Solve time: {elapsed:.3f} s')
    if len(residuals) > 0:
        print(f'Tracked residuals: count = {len(residuals)}, initial ~ {residuals[0]:.3e}, final ~ {residuals[-1]:.3e}')
    else:
        print('No residual history recorded (callback may not have been invoked or iterative fallback used)')

    # Compute final residual norm explicitly
    final_res = np.linalg.norm(rhs - A_mat.dot(x))
    print(f'Explicit final residual norm ||b - A x||_2 = {final_res:.3e}')

    # Interpolate back to full grid and plot
    sol_grid = interpolate_solution_to_grid(x, idx_map, xs, ys, node_type, outer_bc, circle_bc)
    plot_solution(sol_grid, xs, ys, circles, fname='solution_contour.png')
    print('Saved contour plot to solution_contour.png')

    # Print a few sample values for verification
    sample_points = [(-0.9, -0.9), (0.0, 0.0), (0.5, 0.5), (0.8, -0.8)]
    print('Sample solution values (x,y) -> u :')
    for (sx, sy) in sample_points:
        # find closest grid node
        i = int(round((sx - xmin) / (xmax - xmin) * (Nx - 1)))
        j = int(round((sy - ymin) / (ymax - ymin) * (Ny - 1)))
        i = max(0, min(Nx - 1, i))
        j = max(0, min(Ny - 1, j))
        val = sol_grid[j, i]
        print(f'  ({sx:.2f}, {sy:.2f}) -> u = {val:.6f}')

    # Summary including visc (to demonstrate usage of global visc)
    print(f'Viscosity (global visc) = {visc} (Re= {Re})')
```


#### Script block3:

```python
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse.linalg import LinearOperator
import matplotlib.pyplot as plt
import time

# Global visc definition as requested (use visc instead of nu)
Re = 1.0
visc = 1.0 / Re

# --- Problem constants (passed explicitly to functions) ---
k = 8.0
mu1 = 1.0
mu2 = 4.0
A_amp = 10.0  # amplitude (renamed from A to avoid name clash with matrix A)

# Circle definitions: list of (cx, cy, r)
circles = [
    (0.5, 0.5, 0.2),
    (0.4, -0.4, 0.4),
    (-0.2, -0.7, 0.1),
    (-0.6, 0.5, 0.3),
]

# Boundary Dirichlet values
outer_bc = 0.2
circle_bc = 1.0

# Grid resolution (node-centered)
Nx = 201  # choose sufficiently fine to resolve smallest circle
Ny = 201
xmin, xmax = -1.0, 1.0
ymin, ymax = -1.0, 1.0


# --- Helper functions ---

def f_rhs(x, y, mu1, mu2, A_amp):
    """Right-hand side f(x,y).
    f = A * (mu1^2 + mu2^2 + x^2 + y^2) * sin(mu1*pi*x) * sin(mu2*pi*y)
    """
    return A_amp * ((mu1 ** 2 + mu2 ** 2) + x ** 2 + y ** 2) * np.sin(mu1 * np.pi * x) * np.sin(mu2 * np.pi * y)


def edge_circle_intersection(p, q, circle):
    """Return t in (0,1) where segment p->q intersects circle (first interior intersection), or None.
    Solves |p + t*(q-p) - c|^2 = r^2.
    """
    cx, cy, r = circle
    px, py = p
    qx, qy = q
    dx = qx - px
    dy = qy - py
    a = dx * dx + dy * dy
    if a == 0:
        return None
    b = 2.0 * (dx * (px - cx) + dy * (py - cy))
    c0 = (px - cx) ** 2 + (py - cy) ** 2 - r * r
    disc = b * b - 4 * a * c0
    if disc < 0:
        return None
    sqrt_disc = np.sqrt(disc)
    t1 = (-b - sqrt_disc) / (2 * a)
    t2 = (-b + sqrt_disc) / (2 * a)
    # return the first intersection inside the open segment
    for t in (t1, t2):
        if 1e-12 < t < 1 - 1e-12:
            return float(t)
    return None


def classify_nodes(xs, ys, circles, h, eps_factor=1.0, use_nearest_dirichlet=False):
    """Classify each grid node relative to the circles.
    node_type values:
      0 -> hole interior (exclude)
      1 -> domain node (unknown)
      2 -> outer rectangle boundary (Dirichlet known)
      3 -> circle boundary node (Dirichlet known)
    eps = eps_factor * h defines tolerance for boundary classification.
    If use_nearest_dirichlet=True, any node with dist <= r + 0.5*h is treated as circle Dirichlet
    (simple baseline verification mode).
    Returns node_type and counts per type.
    """
    Ny = len(ys)
    Nx = len(xs)
    node_type = np.ones((Ny, Nx), dtype=np.int8)  # default domain
    eps = eps_factor * h
    for j in range(Ny):
        for i in range(Nx):
            x = xs[i]
            y = ys[j]
            # outer rectangle boundary
            if i == 0 or i == Nx - 1 or j == 0 or j == Ny - 1:
                node_type[j, i] = 2
                continue
            # compute signed distance to nearest circle: signed = dist - r
            min_signed = np.inf
            for (cx, cy, r) in circles:
                dist = np.hypot(x - cx, y - cy)
                signed = dist - r
                if signed < min_signed:
                    min_signed = signed
            if use_nearest_dirichlet:
                # baseline: mark nodes within r + 0.5*h as Dirichlet
                if min_signed <= 0.5 * h:
                    node_type[j, i] = 3
                elif min_signed < 0:
                    node_type[j, i] = 0
                else:
                    node_type[j, i] = 1
            else:
                if min_signed < -eps:
                    node_type[j, i] = 0  # strictly inside hole
                elif abs(min_signed) <= eps:
                    node_type[j, i] = 3  # circle boundary node (Dirichlet enforced strongly)
                else:
                    node_type[j, i] = 1  # domain unknown
    # counts
    unique, counts = np.unique(node_type, return_counts=True)
    counts_dict = {int(u): int(c) for u, c in zip(unique, counts)}
    return node_type, counts_dict


def compute_ghost_coeffs(t, h, u_b):
    """Given intersection fraction t (distance from center to boundary divided by h),
    return (diag_increase, rhs_increase) produced by eliminating ghost neighbor
    using linear interpolation to the Dirichlet boundary u_b.
    Algebraic derivation (1D linear interpolation):
      u_g = (1/t)*u_b - ((1-t)/t)*u_i
    substitution into off-diagonal -1/h^2 * u_g yields:
      contributes to diag: (1-t)/(t*h^2)
      contributes to rhs: (1/(t*h^2))*u_b
    """
    diag_inc = (1.0 - t) / (t * h * h)
    rhs_inc = (1.0 / (t * h * h)) * u_b
    return diag_inc, rhs_inc


def assemble_system(Nx, Ny, xmin, xmax, ymin, ymax, circles, k, mu1, mu2, A_amp,
                    outer_bc, circle_bc, visc, eps_factor=1.0, use_nearest_dirichlet=False):
    """Assemble sparse linear system once. Returns matrix A (CSR), rhs, idx_map, xs, ys, node_type,
    irregular_count, irregular_per_circle, and node_type_counts.
    """
    hx = (xmax - xmin) / (Nx - 1)
    hy = (ymax - ymin) / (Ny - 1)
    assert abs(hx - hy) < 1e-12, "Non-uniform spacing not supported in this implementation"
    h = hx

    xs = np.linspace(xmin, xmax, Nx)
    ys = np.linspace(ymin, ymax, Ny)

    # classify nodes
    node_type, node_type_counts = classify_nodes(xs, ys, circles, h, eps_factor, use_nearest_dirichlet)

    # Map (i,j) to unknown index for all nodes that are NOT hole interior (types 0)
    idx_map = -np.ones((Ny, Nx), dtype=int)
    counter = 0
    for j in range(Ny):
        for i in range(Nx):
            if node_type[j, i] != 0:
                idx_map[j, i] = counter
                counter += 1
    N = counter

    # Preallocate lists for COO assembly
    rows = []
    cols = []
    vals = []
    rhs = np.zeros(N, dtype=float)

    irregular_count = 0
    irregular_per_circle = {ci: 0 for ci in range(len(circles))}

    # Precompute f at nodes
    fvals = np.zeros((Ny, Nx), dtype=float)
    for j in range(Ny):
        for i in range(Nx):
            fvals[j, i] = f_rhs(xs[i], ys[j], mu1, mu2, A_amp)

    # Stencil neighbor offsets (i,j offsets correspond to x,y grid indexing)
    neigh = [(-1, 0), (1, 0), (0, -1), (0, 1)]

    for j in range(Ny):
        for i in range(Nx):
            ntype = node_type[j, i]
            if ntype == 0:
                continue  # hole interior not represented
            idx = idx_map[j, i]
            x = xs[i]
            y = ys[j]

            # Outer rectangle boundary nodes get Dirichlet = outer_bc
            if ntype == 2:
                rows.append(idx); cols.append(idx); vals.append(1.0)
                rhs[idx] = outer_bc
                continue

            # Circle boundary nodes get Dirichlet = circle_bc (enforce strongly)
            if ntype == 3:
                rows.append(idx); cols.append(idx); vals.append(1.0)
                rhs[idx] = circle_bc
                continue

            # Interior domain node: assemble -Delta u + k^2 u = f
            diag = 4.0 / (h * h) + k * k
            rhs_val = fvals[j, i]

            for d in neigh:
                ni = i + d[0]
                nj = j + d[1]
                # neighbor should be within grid; since this is interior node, it is
                if node_type[nj, ni] == 0:
                    # edge crosses into a hole interior; find which circle intersects
                    p = (x, y)
                    q = (xs[ni], ys[nj])
                    found = False
                    for ci, circ in enumerate(circles):
                        t = edge_circle_intersection(p, q, circ)
                        if t is not None:
                            diag_inc, rhs_inc = compute_ghost_coeffs(t, h, circle_bc)
                            diag += diag_inc
                            rhs_val += rhs_inc
                            irregular_count += 1
                            irregular_per_circle[ci] += 1
                            found = True
                            break
                    if not found:
                        # Robust fallback: treat as if neighbor is known circle boundary
                        diag += 1.0 / (h * h)
                        rhs_val += (1.0 / (h * h)) * circle_bc
                        irregular_count += 1
                else:
                    # neighbor is represented in system (could be outer boundary or circle boundary or domain unknown)
                    nidx = idx_map[nj, ni]
                    if node_type[nj, ni] in (2, 3):
                        known_val = outer_bc if node_type[nj, ni] == 2 else circle_bc
                        rhs_val += (1.0 / (h * h)) * known_val
                    else:
                        # Unknown neighbor contributes off-diagonal
                        rows.append(idx)
                        cols.append(nidx)
                        vals.append(-1.0 / (h * h))

            # insert diagonal and rhs
            rows.append(idx); cols.append(idx); vals.append(diag)
            rhs[idx] = rhs_val

    A_mat = sp.coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsr()

    return A_mat, rhs, idx_map, xs, ys, node_type, irregular_count, irregular_per_circle, node_type_counts


def make_callback_factory(A_mat, rhs, residuals):
    """Create a robust callback that handles different scipy versions."""
    def callback(arg):
        try:
            if isinstance(arg, (int, np.integer)):
                return
            arr = np.asarray(arg)
            if arr.shape[0] == rhs.shape[0]:
                try:
                    rk = rhs - A_mat.dot(arr)
                    res = np.linalg.norm(rk)
                except Exception:
                    res = np.linalg.norm(arr)
                residuals.append(res)
        except Exception:
            return

    return callback


def try_build_preconditioner(A_mat, drop_tol=1e-4, fill_factor=10):
    """Attempt to build an ILU preconditioner; return LinearOperator or None."""
    try:
        A_csc = A_mat.tocsc()
        ilu = spla.spilu(A_csc, drop_tol=drop_tol, fill_factor=fill_factor)
        Mx = lambda x: ilu.solve(x)
        M = LinearOperator(dtype=A_mat.dtype, shape=A_mat.shape, matvec=Mx)
        print(f'Built ILU preconditioner (spilu) drop_tol={drop_tol} fill_factor={fill_factor}')
        return M
    except Exception as e:
        print(f'Could not build ILU preconditioner: {e}')
        return None


def matrix_diagnostics(A_mat):
    """Print diagnostics for sparse matrix A_mat."""
    data = A_mat.data
    print('Matrix entry stats: min={:.3e}, max={:.3e}, mean_abs={:.3e}'.format(data.min(), data.max(), np.mean(np.abs(data))))
    # symmetry measure
    try:
        A_diff = A_mat - A_mat.T
        frob_diff = spla.norm(A_diff)
        frob_A = spla.norm(A_mat)
        sym_rel = frob_diff / (frob_A + 1e-30)
        print(f'Symmetry measure ||A - A^T||_2 / ||A||_2 = {sym_rel:.3e}')
    except Exception as e:
        print('Could not compute symmetry measure:', e)
    # zero diagonal count
    diag = A_mat.diagonal()
    zero_diag = np.sum(np.isclose(diag, 0.0))
    print(f'Zero diagonal entries: {zero_diag} / {A_mat.shape[0]}')
    # max absolute row sum
    abs_row_sum = np.max(np.sum(np.abs(A_mat), axis=1))
    print(f'Max absolute row-sum = {abs_row_sum:.3e}')
    # eigen estimates (symmetric part)
    try:
        A_sym = 0.5 * (A_mat + A_mat.T)
        # largest eigenvalue
        vals_large = spla.eigsh(A_sym, k=1, which='LA', return_eigenvectors=False, tol=1e-2)
        lam_max = float(vals_large[0])
        # smallest algebraic eigenvalue
        vals_small = spla.eigsh(A_sym, k=1, which='SA', return_eigenvectors=False, tol=1e-2)
        lam_min = float(vals_small[0])
        print(f'Estimated eigenvalues of sym(A): min={lam_min:.3e}, max={lam_max:.3e}, cond_est={abs(lam_max/lam_min) if lam_min!=0 else np.inf:.3e}')
    except Exception as e:
        print('Eigenvalue estimates unavailable:', e)


def try_direct_solve(A_mat, rhs):
    """Attempt direct solve and report residual and timing."""
    start = time.time()
    x_direct = spla.spsolve(A_mat.tocsc(), rhs)
    elapsed = time.time() - start
    res = np.linalg.norm(rhs - A_mat.dot(x_direct))
    print(f'Direct solve (spsolve) time: {elapsed:.3f} s, residual ||b - A x|| = {res:.3e}')
    return x_direct, res, elapsed


def run_iterative_solve(A_mat, rhs, tol=1e-8, maxiter=2000):
    """Run iterative solves with fallbacks and diagnostics."""
    residuals = []
    callback = make_callback_factory(A_mat, rhs, residuals)

    try:
        A_sym = 0.5 * (A_mat + A_mat.T)
    except Exception:
        A_sym = A_mat

    M = try_build_preconditioner(A_sym)

    # Try CG if matrix is symmetric enough, otherwise try GMRES
    symmetry_ok = True
    try:
        A_diff_norm = spla.norm(A_sym - A_sym.T)
        symmetry_ok = (A_diff_norm < 1e-10)
    except Exception:
        symmetry_ok = False

    start = time.time()
    method_used = None
    x = None
    info = None

    if symmetry_ok:
        try:
            x, info = spla.cg(A_sym, rhs, tol=tol, maxiter=maxiter, callback=callback, M=M)
            method_used = 'cg'
        except Exception as e:
            print('CG failed:', e)
            x = None
    else:
        try:
            x, info = spla.gmres(A_mat, rhs, tol=tol, restart=50, maxiter=maxiter, callback=callback, M=M)
            method_used = 'gmres'
        except Exception as e:
            print('GMRES failed, trying bicgstab:', e)
            try:
                x, info = spla.bicgstab(A_mat, rhs, tol=tol, maxiter=maxiter, callback=callback)
                method_used = 'bicgstab'
            except Exception as e2:
                print('bicgstab failed:', e2)
                x = None
    elapsed = time.time() - start

    if x is None:
        print('Iterative methods failed; returning None')
        return None, residuals, elapsed, 'failed', method_used

    status = 'converged' if getattr(info, '__int__', lambda: info)() == 0 else f'info={info}'
    return x, residuals, elapsed, status, method_used


def interpolate_solution_to_grid(x_vec, idx_map, xs, ys, node_type, outer_bc, circle_bc):
    Ny, Nx = idx_map.shape
    sol = np.full((Ny, Nx), np.nan, dtype=float)
    for j in range(Ny):
        for i in range(Nx):
            if node_type[j, i] == 0:
                sol[j, i] = circle_bc
            else:
                sol[j, i] = x_vec[idx_map[j, i]]
    return sol


def plot_solution(sol, xs, ys, circles, fname='solution_contour.png'):
    X, Y = np.meshgrid(xs, ys)
    plt.figure(figsize=(8, 6))
    levels = 60
    cf = plt.contourf(X, Y, sol, levels=levels, cmap='RdBu_r')
    cbar = plt.colorbar(cf, label='u')
    theta = np.linspace(0, 2 * np.pi, 400)
    for (cx, cy, r) in circles:
        plt.plot(cx + r * np.cos(theta), cy + r * np.sin(theta), 'k-', linewidth=1.5)
    plt.title('Solution u (contour) with circle boundaries')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.close()


# --- Unit tests for edge elimination and intersection ---
def test_edge_elimination_and_intersection():
    print('\nRunning unit tests for edge-circle intersection and ghost-coefficient algebra...')
    # Simple synthetic test for compute_ghost_coeffs
    h = 0.1
    t = 0.3
    u_b = 1.0
    diag_inc, rhs_inc = compute_ghost_coeffs(t, h, u_b)
    diag_expected = (1.0 - t) / (t * h * h)
    rhs_expected = (1.0 / (t * h * h)) * u_b
    assert np.allclose(diag_inc, diag_expected)
    assert np.allclose(rhs_inc, rhs_expected)
    print('  Ghost coefficient algebra: PASS')

    # Test edge_circle_intersection produces an interior t when a segment intersects
    p = (0.0, 0.0)
    q = (1.0, 0.0)
    # circle centered at 0.4, radius 0.1 intersects segment in (0,1)
    circ = (0.4, 0.0, 0.1)
    t_int = edge_circle_intersection(p, q, circ)
    assert t_int is not None and 0.0 < t_int < 1.0
    print('  Edge-circle intersection detection: PASS (t = {:.3f})'.format(t_int))
    print('Unit tests passed.\n')


# --- Main execution ---
if __name__ == '__main__':
    # Run unit tests first
    test_edge_elimination_and_intersection()

    print('Assembling system (immersed geometry) with improved diagnostics...')
    # Use eps_factor=1.0 (recommended) to make boundary classification more robust
    A_mat, rhs, idx_map, xs, ys, node_type, irregular_count, irregular_per_circle, node_type_counts = assemble_system(
        Nx, Ny, xmin, xmax, ymin, ymax, circles, k, mu1, mu2, A_amp, outer_bc, circle_bc, visc,
        eps_factor=1.0, use_nearest_dirichlet=False)

    N = A_mat.shape[0]
    nnz = A_mat.nnz
    print(f'Grid: Nx={Nx}, Ny={Ny}, total nodes={Nx * Ny}')
    print(f'Unknowns (represented nodes not strictly inside holes) = {N}')
    print(f'Sparse matrix: shape={A_mat.shape}, nnz={nnz}')
    print('Node type counts:', node_type_counts)
    print(f'Number of irregular (ghost) edge treatments = {irregular_count}')
    if irregular_count > 0:
        print('Irregular edges per circle:', irregular_per_circle)
    else:
        print('No ghost-edge eliminations recorded (interfaces may be captured as Dirichlet nodes).')

    # Matrix diagnostics
    matrix_diagnostics(A_mat)

    # Try direct solve first to verify assembly correctness
    print('\nAttempting direct solve (spsolve) as baseline reference...')
    try:
        x_direct, res_direct, t_direct = try_direct_solve(A_mat, rhs)
    except Exception as e:
        print('Direct solve failed:', e)
        x_direct = None
        res_direct = np.inf

    # If direct solve residual is large, try baseline nearest-grid-Dirichlet assembly to isolate
    if res_direct > 1e-6:
        print('\nDirect solve residual is large -> trying baseline nearest-grid-Dirichlet assembly to get a reference solution')
        A_base, rhs_base, idx_map_base, xs_b, ys_b, node_type_b, irr_count_b, irr_per_circle_b, node_type_counts_b = assemble_system(
            Nx, Ny, xmin, xmax, ymin, ymax, circles, k, mu1, mu2, A_amp, outer_bc, circle_bc, visc,
            eps_factor=1.0, use_nearest_dirichlet=True)
        print('Baseline (nearest-grid-Dirichlet) node counts:', node_type_counts_b)
        matrix_diagnostics(A_base)
        try:
            x_base, res_base, t_base = try_direct_solve(A_base, rhs_base)
            print('Baseline direct solve residual =', res_base)
        except Exception as e:
            print('Baseline direct solve failed:', e)

    # Now attempt iterative solves (CG/GMRES with ILU)
    print('\nStarting iterative solver attempts (CG/GMRES with ILU preconditioner where appropriate)...')
    x_iter, residuals, elapsed_iter, status_iter, method_iter = run_iterative_solve(A_mat, rhs, tol=1e-8, maxiter=5000)
    print(f'Iterative method attempted: {method_iter}, status: {status_iter}, time: {elapsed_iter:.3f} s')
    if len(residuals) > 0:
        print(f'Tracked residuals: count = {len(residuals)}, initial ~ {residuals[0]:.3e}, final ~ {residuals[-1]:.3e}')
    else:
        print('No residual history recorded from iterative solver')

    # Choose solution to use for plotting: prefer direct if available and residual small, else iterative if converged, else fallback to direct
    sol_vec = None
    if x_direct is not None and res_direct < 1e-6:
        print('Using direct spsolve solution (residual small)')
        sol_vec = x_direct
    elif x_iter is not None and status_iter == 'converged':
        print('Using iterative solver solution (converged)')
        sol_vec = x_iter
    elif x_direct is not None:
        print('Iterative did not converge but direct solve available; using direct solution')
        sol_vec = x_direct
    else:
        raise RuntimeError('No valid solution available from either direct or iterative methods')

    final_res = np.linalg.norm(rhs - A_mat.dot(sol_vec))
    print(f'Explicit final residual norm ||b - A x||_2 = {final_res:.3e}')

    # Interpolate back to full grid and plot
    sol_grid = interpolate_solution_to_grid(sol_vec, idx_map, xs, ys, node_type, outer_bc, circle_bc)
    plot_solution(sol_grid, xs, ys, circles, fname='solution_contour.png')
    print('Saved contour plot to solution_contour.png')

    # Print a few sample values for verification
    sample_points = [(-0.9, -0.9), (0.0, 0.0), (0.5, 0.5), (0.8, -0.8)]
    print('Sample solution values (x,y) -> u :')
    for (sx, sy) in sample_points:
        i = int(round((sx - xmin) / (xmax - xmin) * (Nx - 1)))
        j = int(round((sy - ymin) / (ymax - ymin) * (Ny - 1)))
        i = max(0, min(Nx - 1, i))
        j = max(0, min(Ny - 1, j))
        val = sol_grid[j, i]
        print(f'  ({sx:.2f}, {sy:.2f}) -> u = {val:.6f}')

    print(f'Viscosity (global visc) = {visc} (Re= {Re})')
```



### runtime_outputs

#### Output block1

Assembling system...
Grid: Nx=201, Ny=201, total nodes=40401
Unknowns (nodes not inside holes) = 30999
Sparse matrix: shape=(30999, 30999), nnz=150195
Number of irregular (ghost) edge treatments = 804
Starting solver (Conjugate Gradient)...

Runtime Error:
Traceback (most recent call last):
  File "C:\Users\admin\AppData\Local\Temp\ipykernel_48940\2663542217.py", line 138, in execute_code_tool
    exec(code, global_scope, global_scope)
  File "<string>", line 261, in <module>
  File "<string>", line 198, in solve_system
TypeError: cg() got an unexpected keyword argument 'tol'



#### Output block2

Assembling system...
Grid: Nx=201, Ny=201, total nodes=40401
Unknowns (represented nodes not strictly inside holes) = 31269
Sparse matrix: shape=(31269, 31269), nnz=149005
Number of irregular (ghost) edge treatments = 0
Starting solver (robust Conjugate Gradient with fallbacks)...
Built ILU preconditioner (spilu)
Solver method/result: cg with kwargs=['atol', 'maxiter', 'callback', 'M']
Solver status: no convergence after 5000 iterations
Solve time: 15.832 s
Tracked residuals: count = 5000, initial ~ 8.334e+04, final ~ 1.131e+04
Explicit final residual norm ||b - A x||_2 = 1.131e+04
Saved contour plot to solution_contour.png
Sample solution values (x,y) -> u :
  (-0.90, -0.90) -> u = -0.104956
  (0.00, 0.00) -> u = 0.296282
  (0.50, 0.50) -> u = 1.000000
  (0.80, -0.80) -> u = 0.615584
Viscosity (global visc) = 1.0 (Re= 1.0)



#### Output block3


Running unit tests for edge-circle intersection and ghost-coefficient algebra...
  Ghost coefficient algebra: PASS
  Edge-circle intersection detection: PASS (t = 0.300)
Unit tests passed.

Assembling system (immersed geometry) with improved diagnostics...
Grid: Nx=201, Ny=201, total nodes=40401
Unknowns (represented nodes not strictly inside holes) = 31624
Sparse matrix: shape=(31624, 31624), nnz=148422
Node type counts: {0: 8777, 1: 29606, 2: 800, 3: 1218}
Number of irregular (ghost) edge treatments = 0
No ghost-edge eliminations recorded (interfaces may be captured as Dirichlet nodes).
Matrix entry stats: min=-1.000e+04, max=4.006e+04, mean_abs=1.586e+04
Symmetry measure ||A - A^T||_2 / ||A||_2 = 0.000e+00
Zero diagonal entries: 0 / 31624
Max absolute row-sum = 8.006e+04
Estimated eigenvalues of sym(A): min=1.000e+00, max=7.976e+04, cond_est=7.976e+04

Attempting direct solve (spsolve) as baseline reference...
Direct solve (spsolve) time: 0.067 s, residual ||b - A x|| = 6.325e-10

Starting iterative solver attempts (CG/GMRES with ILU preconditioner where appropriate)...
Built ILU preconditioner (spilu) drop_tol=0.0001 fill_factor=10
CG failed: cg() got an unexpected keyword argument 'tol'
Iterative methods failed; returning None
Iterative method attempted: None, status: failed, time: 0.000 s
No residual history recorded from iterative solver
Using direct spsolve solution (residual small)
Explicit final residual norm ||b - A x||_2 = 6.325e-10
Saved contour plot to solution_contour.png
Sample solution values (x,y) -> u :
  (-0.90, -0.90) -> u = -0.104924
  (0.00, 0.00) -> u = 0.302314
  (0.50, 0.50) -> u = 1.000000
  (0.80, -0.80) -> u = 0.618568
Viscosity (global visc) = 1.0 (Re= 1.0)




### review_decision
accept

### review_comments
	Current Stage [A/3]
1. MAJOR

Category: runtime
Issue: Runtime failure: scipy.sparse.linalg.cg() raised TypeError: cg() got an unexpected keyword argument 'tol'.

Recommendation: Make the cg() call compatible with different SciPy versions. Wrap the call in a try/except and fall back to alternate keyword names or a direct solver. For example:
- Primary attempt: x, info = spla.cg(A, rhs, tol=tol, maxiter=maxiter, callback=callback)
- except TypeError: try x, info = spla.cg(A, rhs, maxiter=maxiter, callback=callback) (omit tol) or use spla.cg(A, rhs, atol=tol, maxiter=maxiter, callback=callback)
As a robust fallback, if iterative solver fails, solve directly with spla.spsolve(A, rhs) or use an alternative Krylov with consistent keyword names (gmres/minres). Also print a clear fallback message. Consider adding a preconditioner (spilu + LinearOperator) for performance/robustness.


2. MINOR

Category: runtime
Issue: Callback compatibility and residual tracking: the code assumes the CG callback receives the current solution vector xk and computes rk = b - A xk. Different SciPy versions pass either the current residual vector, the current solution, or iteration count to the callback, which can cause the callback to fail or produce misleading residual history.

Recommendation: Make the callback robust: implement a callback that accepts a single argument and try to detect its meaning. Example strategy:
- In the wrapper, attempt to compute norm assuming callback receives xk; if that raises a shape/TypeError, assume callback receives residual vector and compute its norm directly. If callback receives an int, stop tracking.
- Alternatively, avoid relying on callback and instead compute residuals after solve by running a short power iteration or by tracking with a custom solver loop. For diagnostics, printing the explicit final residual norm (already done) is essential.


3. MAJOR

Category: structure
Issue: Boundary node classification: nodes exactly on the circle boundaries are treated as 'inside_hole' (point_in_any_circle uses <= r^2) and therefore excluded from the linear system. The Dirichlet condition at the circle boundary is only enforced via ghost-edge interpolation on adjacent cells and by forcing plotting values inside holes to circle_bc. This can leave boundary grid nodes not strongly enforced in the linear system and creates a consistency / accuracy issue.

Recommendation: Classify and handle circle boundary nodes explicitly. Distinguish three node types: interior domain nodes, circle-boundary nodes (distance approximately equal to r), and strictly interior-to-hole nodes. For nodes close to the circle (|dist - r| < eps, e.g. eps = 0.5*h), mark them as Dirichlet boundary nodes and impose u = circle_bc strongly (append identity rows). This ensures strong enforcement and reduces geometric inconsistency. Alternatively refine the grid locally near circles or use higher-order immersed techniques.


4. MINOR

Category: accuracy
Issue: Accuracy and convergence concerns: the immersed/ghost first-order elimination along edges is reasonable but only first-order accurate at curved boundaries. No checks for symmetry/positive-definiteness or condition-number estimation are provided; CG without preconditioning may converge slowly or fail on larger grids.

Recommendation: Improve robustness and accuracy:
- Optionally symmetrize the assembled matrix (A = 0.5*(A + A.T)) to guard against tiny assembly asymmetries before running CG.
- Add a simple preconditioner (e.g. incomplete LU via spilu + LinearOperator) or use AMG (pyamg) for faster convergence.
- Report the number of CG iterations, final residual, and wall time (you already print residuals and time; keep/extend this). For accuracy: do a grid-refinement study (coarse->fine) and monitor L2 or max norm changes to confirm convergence order.


5. MINOR

Category: structure
Issue: Minor code-quality issues: the global variable visc is declared and passed to assemble_system but not used in the PDE; f_rhs docstring mentions a different interpretation wording which could be confusing.

Recommendation: Either use visc somewhere meaningful or remove it from the function signature to avoid confusion. Clarify f_rhs docstring to match the problem text exactly (show full expression A*(mu1^2 + mu2^2 + x^2 + y^2)*sin(mu1*pi*x)*sin(mu2*pi*y) if that is intended). Add small unit tests or asserts for input shapes and small sanity checks (e.g. check that irregular_count>0 when expected).



	Current Stage [B/3]
1. MAJOR

Category: runtime
Issue: Solver non-convergence and large final residual

Recommendation: The solver did not converge: final explicit residual ||b - A x||_2 = 1.13e+04 and CG reported no convergence after 5000 iterations. First reproduce with a direct solve (scipy.sparse.linalg.spsolve) and check the residual; if direct solve also gives a large residual then the matrix assembly is incorrect. Compute and print basic diagnostics before solving: min/max/mean(A.data), norm(A - A.T), number of zero diagonal entries, and an estimate of condition number (e.g. via power method or sparse eigsh for a few extreme eigenvalues). If A is non-symmetric or indefinite, use appropriate solvers (minres for symmetric indefinite, gmres/bicgstab for general). As an immediate robustness step, try solving once with spsolve and report residual and timing to isolate whether iterative solver or assembly is the problem.


2. MAJOR

Category: accuracy
Issue: Ghost-point elimination coefficients are likely incorrect or not validated

Recommendation: Investigate and fix the immersed/ghost-edge elimination formulas. The current ghost treatment (diag += (1-t)/(t*h*h); rhs += (1/(t*h*h))*circle_bc) is unusual and likely incorrect or inconsistent with the finite-difference stencil. Derive the elimination algebraically: for a neighbor ghost point g at fractional distance alpha (distance from grid node to boundary divided by h), express u_g in terms of u_i and boundary u_b, substitute into the 5-point Laplacian and collect terms to get the correct modifications to the diagonal and RHS. Implement and test the algebraic formula on a simple 1D/2D case with a single circular hole where an analytic/benchmark solution is known or where you can enforce Dirichlet on a cut with a manufactured solution. Add unit tests for the elimination routine (edge_circle_intersection + coefficient computation) so the irregular_count > 0 is exercised.


3. MAJOR

Category: structure
Issue: Too small tolerance for classifying circle boundary nodes; no diagnostics printed for classification counts

Recommendation: Classification of boundary vs interior nodes is fragile: eps = 0.5*h is very small and may under-represent the circle boundary and mis-locate the interface. Use a more robust threshold (e.g. eps = h or 1.5*h) or, better, compute fractional intersection distance per edge and treat the interface as subcell (level-set) geometry. Print counts of nodes in each node_type (counts of 0,1,2,3) to verify that circle boundary nodes and hole-interior nodes exist in expected numbers. Without sufficient representation of boundary nodes the ghost elimination paths or Dirichlet enforcement will be ineffective and can produce large errors/ill-conditioning.


4. MINOR

Category: structure
Issue: Insufficient solver and matrix diagnostics to diagnose non-convergence

Recommendation: Add and print robust matrix and solver diagnostics before/after solve: symmetricity measure ||A - A^T||_F / ||A||_F, smallest and largest diagonal entries, max absolute row-sum, and a cheaply computed estimate of the smallest positive eigenvalue (using scipy.sparse.linalg.eigsh for a few eigenpairs). Also provide counts of how many edges triggered the ghost elimination (irregular_count) and which circles generated them. These diagnostics are essential to debug conditioning and to decide on preconditioner strategy (ILU tuning, AMG, or direct solve). If ILU is used, print fill_factor and drop_tol and verify ilu.solve produces reasonable approximations; consider using pyamg.smoothed_aggregation_solver or PETSc/pyAMG for robust preconditioning.


5. MAJOR

Category: accuracy
Issue: No baseline/reference solution or fallback validated; immersed approach not validated against simple cases

Recommendation: Numerical and implementation simplifications to get a reliable baseline: (1) For verification, solve a manufactured-problem variant where circle boundaries coincide with grid lines (so no immersed geometry) and verify convergence and correctness. (2) Temporarily enforce Dirichlet on the nearest grid nodes to the circles (i.e. mark any grid node with dist <= r + h/2 as Dirichlet circle boundary) to avoid ghost elimination and confirm solver performance. (3) If iterative methods remain problematic, use spsolve for final runs (possibly with sparse Cholesky) to get a reference solution, then refine the immersed-boundary implementation to match that reference. (4) Review f_rhs vs problem statement (you included x^2+y^2 term — that is fine but document it explicitly).



	Current Stage [C/3]
1. MINOR

Category: runtime
Issue: Iterative solver attempt failed with an unexpected keyword argument 'tol' for cg(), causing the iterative path to abort even though a valid direct solution was obtained.

Recommendation: Make iterative-solver calls robust to SciPy API differences. Wrap solver calls in try/except and detect available kwarg names (e.g., tolerate 'tol' vs 'atol'), or call via a small wrapper that inspects the function signature and adapts kwargs. Also provide a clear fallback sequence (e.g., try cg, then gmres, then bicgstab) and print which signature branch was used. This will avoid spurious 'iterative failed' reports in environments with different SciPy versions.


2. MINOR

Category: accuracy
Issue: The immersed/ghost-edge elimination machinery exists but was not exercised on this run: irregular_count = 0. The code therefore enforced circle Dirichlet values by classifying nearby grid nodes as boundary nodes (abs(dist - r) <= eps) which produces a staircasing approximation of curved boundaries.

Recommendation: Either (a) reduce eps (use a fraction of h) so that at least some edges are detected as crossing into holes and the ghost-elimination logic is tested, or (b) explicitly detect cell-edge intersections during classification so the ghost-elimination path is triggered. If the goal is higher geometric fidelity, enable and validate the ghost-edge elimination path (unit tests that force ghost treatment are useful). Document the tradeoff between capturing boundary nodes as Dirichlet (simple, but staircased) vs. using ghost elimination (more accurate if implemented robustly).


3. MINOR

Category: structure
Issue: Global visc is defined and passed around but not actually used in any operator or computation; this may be confusing and departs from the instruction to 'strictly use visc instead of nu' (there is no nu in the code).

Recommendation: If visc is not needed for this elliptic problem, either remove it from function signatures or document explicitly why visc is defined but unused. If visc was intended for future extensions, keep it but use it consistently or annotate with a comment to avoid confusion.


4. MINOR

Category: structure
Issue: interpolate_solution_to_grid fills grid points inside the removed circles with the circle Dirichlet value (circle_bc) rather than masking them (e.g., setting to NaN). This can be misleading when visualizing the removed holes (plot will appear as if hole interiors are filled with a solution value).

Recommendation: For plotting and postprocessing, mask hole-interior nodes (e.g., set them to np.nan) or use a matplotlib mask so the holes are visibly empty. If you prefer to display the Dirichlet value inside holes, document that behavior clearly. Alternatively, use pcolormesh with a masked array so the true multiply-connected nature of the domain is preserved in visuals.


5. MINOR

Category: accuracy
Issue: Conditioning and preconditioning: the matrix diagnostics show a moderate condition estimate (~8e4). The code attempts ILU with specific drop_tol and fill_factor; however iterative solves may be sensitive to these parameters and to the chosen solver (CG requires symmetry/PD).

Recommendation: Document or expose ILU parameters and consider adding an AMG/pyamg fallback for larger problems where direct solve becomes expensive. For CG use, ensure the matrix is symmetric positive definite (current assembly yields a symmetric SPD matrix for the chosen discretization). Record and print preconditioner statistics (setup time, apply time) when available. If iterative solves are important, add an option to run an inexpensive diagonal scaling (Jacobi) as a baseline preconditioner and compare iteration counts.





### iteration_count
3

### rev_image_description
	Current Stage [A/3]


	Current Stage [B/3]


	Current Stage [C/3]




