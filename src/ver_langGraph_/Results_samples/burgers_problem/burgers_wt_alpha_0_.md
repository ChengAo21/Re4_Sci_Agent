### 🛌🏻 Re4gent with Multi-Modal Review

#### Configs 🏷️:
MaxCount: 2, MultiModalReview: False, AlternativeNum: 2, MaxComments: 5, ConsModule_ON: True, CUT_OUTEXT: 2000

### prob_todo

The PDE of 1-D Burgers equation is given by:
\begin{cases}
u_t + u u_x - \nu u_{xx} = 0, & (x,t) \in \Omega \\
\end{cases}

The parameter is \(\nu = \frac{0.01}{\pi}\).
The domain is \(\Omega = [-1, 1] \times [0, 1]\).
The initial and boundary conditions are:
\begin{cases}
u(x, 0) = -\sin(\pi x), & (x,t) \in \Omega \text{ with } t = 0 \\
u(-1, t) = u(1, t) = 0, & (x,t) \in \Omega \\
\end{cases}

Implement a stable and efficient method to solve this problem.
Define visc = 1.0 / Re at the GLOBAL scope (outside main), and strictly use visc instead of nu throughout the code to avoid confusion with u.
Ensure all constants are passed explicitly to functions or defined globally to prevent NameError.
Plot the contour of the velocity magnitude using 'RdBu_r' colormap) and velocity profile at timesteps t=0.2,0.4,0.6,0.8,1.0 in one figure.
Just save figs do not use plt.show() in the code.


### expanded_prob
We must solve the 1D viscous Burgers equation

    u_t + u u_x - \nu u_{xx} = 0,

on the space-time domain x ∈ [-1,1], t ∈ [0,1] with viscosity parameter given numerically by \nu = 0.01/\pi. The initial condition is u(x,0) = -sin(\pi x) and homogeneous Dirichlet boundary conditions u(-1,t)=u(1,t)=0. Implementation constraints and plotting requirements that must be respected by the final code: define visc = 1.0/Re at global scope (outside main) and use visc everywhere in place of \nu; pass all constants explicitly to functions or declare them globally so no NameError arises; produce (a) a contour plot in the x–t plane of the velocity magnitude using the 'RdBu_r' colormap and (b) a single figure that overlays velocity profiles at t = 0.2, 0.4, 0.6, 0.8, 1.0; save the figures to files and do not call plt.show().

Primary mathematical and numerical challenges:
- Nonlinearity: the convection term u u_x makes the problem nonlinear and can generate steep gradients (quasi-shocks) even though viscosity smooths them; the scheme must be able to capture steep gradients without introducing spurious oscillations.
- Stiffness from diffusion when viscosity is small: visc ≈ 0.003183 (0.01/π), so the diffusive term may introduce a restrictive time-step if treated explicitly. A stable method must either use implicit time stepping for diffusion or a stiff-aware time integrator.
- Boundary conditions: Dirichlet zero values at endpoints require spatial discretizations that respect non-periodic boundaries (finite differences, Chebyshev/sine spectral, or appropriate transforms). A naive Fourier periodic pseudospectral approach is not directly applicable unless an odd/sine extension is used.
- Conservation and monotonicity: conservative discretizations for the convective flux f = u^2/2 help preserve integral properties; monotonicity or TVD-style limiters may be necessary to avoid nonphysical oscillations near steep gradients.
- Accuracy vs efficiency tradeoffs: high-order methods (spectral or WENO) reduce number of points needed but are more complex to implement and may require dealiasing when using pseudospectral evaluation of nonlinearities. Low-order schemes are simpler but may require finer grids and smaller time steps.
- Time integration coupling: choosing IMEX (implicit-explicit) or integrating-factor/ETD approaches to combine stable handling of diffusion with accurate explicit treatment of nonlinearity.

Quantities to output/visualize: contour |u(x,t)| on x–t grid with colormap 'RdBu_r' and time-slice velocity profiles at t = {0.2,0.4,0.6,0.8,1.0} in one figure. Ensure saved files only (no plt.show()).

### solution_plans
	Current Stage [A/2]
solu_name='Finite-difference IMEX (SSP-RK3 convective + Crank–Nicolson diffusion) with conservative flux and tridiagonal solver' content="Governing idea:\nTreat the nonlinear convection explicitly with a strong-stability-preserving Runge–Kutta (SSP-RK3) time integrator using a conservative discretization of the flux f(u)=u^2/2 (with a high-resolution option for the spatial derivative), and treat diffusion implicitly with a second-order-in-time Crank–Nicolson / theta-step so the viscous term does not restrict the time step severely. The diffusion linear operator yields a tridiagonal system that can be solved efficiently with the Thomas algorithm.\n\nAlgorithmic steps (numbered):\n1) Grid and constants: choose Nx uniform grid on x ∈ [-1,1] with spacing dx, time-step dt chosen from a CFL condition for convection: dt ≤ CFL*dx/max(|u|) (use CFL ≈ 0.4–0.6). Define global visc = 1.0/Re externally and pass it to all routines.\n2) Spatial discretization:\n   a) Represent the PDE in conservative form: u_t + (f(u))_x = visc u_{xx}, where f(u)=u^2/2.\n   b) Discretize (f(u))_x using a conservative numerical flux. Two options:\n      - Option A (robust and simple): 2nd-order upwind-biased (MUSCL) flux with slope limiters (minmod or van Leer) for TVD properties.\n      - Option B (higher accuracy): 5th-order WENO reconstruction for the flux values at cell faces.\n   c) Discretize u_{xx} with a centered second-order finite-difference; assemble the implicit diffusion operator as a tridiagonal matrix L (dependent only on dx and visc).\n3) Time stepping — IMEX splitting using SSP-RK3 for the explicit convective part and Crank–Nicolson for diffusion:\n   For each full time step from t^n to t^{n+1}:\n   a) Compute convective right-hand side R(u) ≈ −(f(u))_x using the chosen flux discretization.\n   b) Use SSP-RK3 stages: for each RK stage k, form u_stage = combination of previous stages and explicit R(u_stage). After each explicit stage, apply an implicit Crank–Nicolson solve for diffusion:\n      (I − (dt*theta) visc Dxx) u_stage_new = RHS_stage where Dxx is the second-difference operator and theta = 0.5 for CN.\n      Solve the tridiagonal system by the Thomas algorithm (O(N) cost).\n   c) Proceed through the 3 SSP stages to produce u^{n+1}.\n4) Boundary conditions: enforce u=0 at x = ±1 explicitly before/after each solve; ensure tridiagonal solver accounts for Dirichlet rows.\n5) Diagnostics and adaptivity: monitor max(|u|) to update dt via CFL; optionally implement adaptive dt if desired.\n6) Postprocessing for plotting:\n   a) Save u(x,t) at time slices needed (t = 0.2,0.4,0.6,0.8,1.0) and accumulate u on an x–t mesh to build the contour of |u|.\n   b) Create contour plot of |u(x,t)| with colormap 'RdBu_r' and a separate figure overlaying the specified time-slice profiles; save figures to files (no plt.show()).\n\nStability / accuracy / complexity / efficiency limitations:\n- Stability: explicit treatment of convection enforces a CFL constraint dt ≲ CFL * dx / max(|u|). The implicit diffusion removes the parabolic timestep restriction, but if visc is extremely small and convection dominates, dt must still be small for stability and accuracy of the explicit RK3 stage.\n- Accuracy: spatial accuracy depends on flux reconstructor: MUSCL gives second-order, WENO up to fifth-order. Time accuracy is overall second-to-third order depending on the IMEX coupling. If Crank–Nicolson (theta=0.5) is used at each RK stage, the combined scheme is typically second-order in time unless a formal IMEX RK framework elevates global order.\n- Complexity: WENO reconstructions are more expensive (O(N) per flux evaluation with larger constants) but still linear in N. Each time step requires one or more tridiagonal solves (O(N) each) per RK stage.\n- Efficiency: overall cost per time step is moderate (flux reconstruction + O(N) linear solves). For smooth solutions a high-order WENO reduces required N; for sharp gradients MUSCL may be cheaper but require finer grid.\n- Implementation caveats: careful implementation of flux limiters or WENO weights is necessary to avoid spurious oscillations; ensure correct handling of Dirichlet boundaries in the flux routines. All physical constants (visc, dx, dt, CFL) must be passed explicitly to functions or declared globally to avoid NameError and to follow the coding constraint."

	Current Stage [B/2]
solu_name='Sine-transform pseudospectral (collocation) with integrating-factor / ETDRK4 time integration' content="Governing idea:\nUse a spectral spatial discretization tailored to the homogeneous Dirichlet endpoints by expanding u in an odd sine basis (discrete sine transform, DST), which diagonalizes the Laplacian with Dirichlet BCs. Treat the linear diffusion term exactly via an integrating factor or an exponential time-differencing Runge–Kutta (ETDRK4) integrator and evaluate the nonlinear convective term pseudospectrally (transform to physical space, form u u_x, then transform back). This gives very high spatial accuracy for smooth solutions and allows larger time steps than fully explicit diffusive schemes.\n\nAlgorithmic steps (numbered):\n1) Grid and transforms:\n   a) Choose Nx interior collocation points consistent with sine transform (e.g., x_j = -1 + j*dx, j=1..Nx, with x=±1 as boundaries; use DST type that enforces zero at endpoints). dx = 2/(Nx+1).\n   b) Precompute sine-wave wavenumbers k and the eigenvalues of Laplacian: λ_k = −(π k / L)^2 where L = 1 (after scaling) or determined by node spacing.\n   c) Define visc globally and pass to routines.\n2) Initialization: compute u(x,0) = −sin(π x) on collocation nodes and obtain spectral coefficients û via DST of the odd extension (i.e., DST-I or DST-II depending on indexing).\n3) Nonlinear term evaluation (pseudospectral): to compute N(u) = −(u u_x),\n   a) Transform û (spectral) → u (physical) with inverse DST.\n   b) Compute u_x in physical space: either by transforming spectral û multiplied by i k factors or by differentiating in spectral space (preferred: multiply spectral coefficients by (i k) and inverse transform) to get u_x at nodes.\n   c) Form nonlinear product u * u_x in physical space and then transform back to spectral space via DST.\n   d) To avoid aliasing use 2/3-rule padding: zero-pad to 3/2*Nx in spectral space or perform explicit dealiasing by zeroing high modes after multiplication.\n4) Time integration: use integrating-factor ETDRK4 (or IFRK4) adapted to the sine basis where the linear operator is diagonal in spectral space.\n   a) Precompute the integrating-factor exponentials exp(L*dt) and ETDRK4 coefficients that depend on L*dt (diagonal operations in spectral space).\n   b) At each time step advance û using ETDRK4 where nonlinear evaluations are done pseudospectrally as in step 3.\n5) Boundaries: sine basis enforces u(±1)=0 automatically; no special ghost points needed.\n6) Postprocessing & plotting:\n   a) Save physical-space u(x,t) at desired time slices and accumulate snapshot matrix for x–t contour.\n   b) Produce contour of |u| on the x–t grid with 'RdBu_r' colormap and a single figure overlaying profiles at t = 0.2,0.4,0.6,0.8,1.0; save figures (no plt.show()).\n\nStability / accuracy / complexity / efficiency limitations:\n- Stability: ETDRK4 / integrating-factor methods treat the linear diffusive stiff part exactly and permit larger dt than fully explicit diffusion methods. The remaining timestep restriction arises from resolving the nonlinear advection and accuracy rather than linear stiffness; practical dt must still be chosen to control nonlinear aliasing and temporal truncation error.\n- Accuracy: spectral in space gives exponential (spectral) convergence for smooth solutions; ETDRK4 is fourth-order in time (assuming nonlinear evaluation and de-aliasing are handled consistently). However, accuracy deteriorates if solution develops very steep gradients unless Nx is increased.\n- Complexity: each nonlinear evaluation requires forward and inverse DSTs (O(N log N) cost). With dealiasing (2/3 rule) the transforms operate on a padded array of size ≈ 3/2 N, increasing cost. Precomputations for ETDRK4 coefficients involve evaluating functions of diagonal matrix entries but are cheap and done once.\n- Efficiency: for smooth solutions this method is extremely efficient (fewer grid points for high accuracy). Implementation complexity is higher: correct DST conventions, padding/dealiasing, and indexing must be handled carefully. The method assumes solution regularity; if shocks or underresolved steep gradients appear, spectral methods produce Gibbs oscillations unless appropriate filtering or limiting is applied.\n- Implementation caveats: enforce global constant usage (visc) and explicit passing of parameters to subroutines; ensure consistent DST type to enforce exact Dirichlet BCs; include dealiasing (2/3 rule) to avoid energy pile-up at high modes. Ensure saved plots meet the colormap and time-slice requirements and that plt.show() is not called."



### technical_spec
	Current Stage [A/1]
Finite-difference IMEX splitting: explicit SSP-RK3 for convection (conservative Rusanov flux) followed by a single Crank–Nicolson implicit solve for diffusion per timestep. Domain discretized on a uniform nodal grid including Dirichlet endpoints; the implicit tridiagonal matrix for CN is assembled and LU-like factors precomputed once before time-stepping (no matrix rebuild in time loop). The code stores snapshots across time for an x–t contour and selected time slices, prints diagnostic norms, and saves two figures: (1) contour of |u(x,t)| with 'RdBu_r' colormap, (2) overlay of u(x) at t=0.2,0.4,0.6,0.8,1.0. All physical constants are explicit and visc is defined at global scope.



### python_codes

#### Script block1:

```python
import numpy as np
import matplotlib
matplotlib.use('Agg')  # ensure no interactive display
import matplotlib.pyplot as plt
import time

# Global physical parameter: define Re so visc = 1.0/Re equals 0.01/pi
Re = np.pi / 0.01
visc = 1.0 / Re  # equals 0.01/pi as required

# -- Numerical utilities -------------------------------------------------

def rusanov_convective_rhs(u, dx):
    """Compute convective RHS R = -d/dx f(u) using the Rusanov flux (local Lax-Friedrichs).
    u: full nodal array including Dirichlet boundary nodes (u[0], u[-1])
    dx: grid spacing
    Returns R array same size as u (zeros at boundaries).
    """
    f = 0.5 * u * u
    # flux at interfaces F_{i+1/2} for i=0..N-2
    # local wave speed alpha = max(|u_i|, |u_{i+1}|)
    umax = np.maximum(np.abs(u[:-1]), np.abs(u[1:]))
    F = 0.5 * (f[:-1] + f[1:]) - 0.5 * umax * (u[1:] - u[:-1])
    # derivative at interior nodes i = 1..N-2: (F_{i+1/2} - F_{i-1/2})/dx
    R = np.zeros_like(u)
    R[1:-1] = - (F[1:] - F[:-1]) / dx
    return R


def precompute_tridiagonal_factor(lower_val, diag_val, upper_val, n):
    """Precompute multipliers and modified diagonal for solving same tridiagonal A x = d
    where A has constant lower, diag, upper entries (except the first lower omitted).
    lower_val, diag_val, upper_val: scalar values for off-diagonals and diagonal
    n: size of system (number of interior points)
    Returns arrays bb (modified diagonal) and m (multipliers) and arrays a,c for convenience.
    """
    a = np.zeros(n)
    c = np.zeros(n)
    b = np.zeros(n)
    # fill arrays: a[0] unused (0), a[1..] = lower_val; c[0..n-2] = upper_val; b[:] = diag_val
    a[1:] = lower_val
    c[:-1] = upper_val
    b[:] = diag_val
    bb = np.empty(n)
    m = np.empty(n)
    bb[0] = b[0]
    m[0] = 0.0
    for i in range(1, n):
        m[i] = a[i] / bb[i - 1]
        bb[i] = b[i] - m[i] * c[i - 1]
    return a, b, c, bb, m


def solve_tridiagonal_precomputed(rhs, a, c, bb, m):
    """Solve tridiagonal system given precomputed factorization arrays.
    rhs: right-hand side vector length n
    a,c: original sub- and super-diagonal arrays (a[0]=0)
    bb,m: as from precompute_tridiagonal_factor
    Returns solution x.
    """
    n = rhs.size
    d = rhs.copy()
    # forward elimination using multipliers m
    for i in range(1, n):
        d[i] -= m[i] * d[i - 1]
    # back substitution
    x = np.empty_like(d)
    x[-1] = d[-1] / bb[-1]
    for i in range(n - 2, -1, -1):
        x[i] = (d[i] - c[i] * x[i + 1]) / bb[i]
    return x


# -- Main solver ---------------------------------------------------------

def burgers_1d_imex(nx=401, CFL=0.4, t_final=1.0, save_every=1):
    """Solve u_t + u u_x - visc u_xx = 0 on x in [-1,1], t in [0,t_final]
    using explicit SSP-RK3 for convection (Rusanov flux) and Crank-Nicolson for diffusion.
    Returns x, times, solution snapshots (nt+1 x nx)
    """
    # Spatial grid including endpoints
    x = np.linspace(-1.0, 1.0, nx)
    dx = x[1] - x[0]
    # initial condition u(x,0) = -sin(pi x)
    u0 = -np.sin(np.pi * x)
    # enforce Dirichlet BCs (should be zero at endpoints analytically)
    u0[0] = 0.0
    u0[-1] = 0.0

    # CFL-based guess for dt; keep dt constant throughout to avoid rebuilding matrices
    umax0 = np.max(np.abs(u0))
    if umax0 <= 0:
        umax0 = 1e-6
    dt_guess = CFL * dx / umax0
    # choose integer number of steps so final time reached exactly
    nsteps = int(np.ceil(t_final / dt_guess))
    dt = t_final / nsteps

    print(f"nx={nx}, dx={dx:.3e}, initial_max|u|={umax0:.6f}, dt={dt:.6e}, nsteps={nsteps}")
    # Precompute tridiagonal matrix for Crank-Nicolson implicit solve on interior nodes
    # System for interior M = nx-2 unknowns: A = I - 0.5*dt*visc*D where D is discrete Laplacian
    M = nx - 2
    lam = visc * dt / (dx * dx)
    # A diag = 1 + 2*0.5*lam = 1 + lam, offdiag = -0.5*lam
    diag_val = 1.0 + lam
    off_val = -0.5 * lam
    # lower and upper entries for a and c arrays used in precompute
    lower_val = off_val
    upper_val = off_val
    a, b, c, bb, m = precompute_tridiagonal_factor(lower_val, diag_val, upper_val, M)

    # time-stepping storage
    times = np.linspace(0.0, t_final, nsteps + 1)
    U = np.zeros((nsteps + 1, nx))
    U[0, :] = u0.copy()

    u = u0.copy()
    t = 0.0
    t0 = time.time()
    # time loop
    for n in range(1, nsteps + 1):
        # --- Explicit SSP-RK3 for convection only ---
        R0 = rusanov_convective_rhs(u, dx)
        u1 = u + dt * R0
        u1[0] = 0.0; u1[-1] = 0.0

        R1 = rusanov_convective_rhs(u1, dx)
        u2 = 0.75 * u + 0.25 * (u1 + dt * R1)
        u2[0] = 0.0; u2[-1] = 0.0

        R2 = rusanov_convective_rhs(u2, dx)
        u_explicit = (1.0 / 3.0) * u + (2.0 / 3.0) * (u2 + dt * R2)
        u_explicit[0] = 0.0; u_explicit[-1] = 0.0

        # --- Implicit Crank-Nicolson diffusion solve for full dt ---
        # RHS_interior = u_explicit_interior + 0.5*dt*visc*D u_explicit_interior
        ue = u_explicit
        # compute D*ue on interior nodes
        Du = (ue[:-2] - 2.0 * ue[1:-1] + ue[2:]) / (dx * dx)
        RHS = ue[1:-1] + 0.5 * dt * visc * Du
        # solve tridiagonal A * u_interior_new = RHS
        u_interior_new = solve_tridiagonal_precomputed(RHS, a, c, bb, m)
        # update u and enforce BCs
        u[1:-1] = u_interior_new
        u[0] = 0.0; u[-1] = 0.0

        U[n, :] = u.copy()
        t += dt
        # optionally print progress
        if n % max(1, nsteps // 10) == 0 or n == nsteps:
            print(f"step {n}/{nsteps}, t={t:.4f}, max|u|={np.max(np.abs(u)):.6f}")

    dt_wall = time.time() - t0
    print(f"Completed in {dt_wall:.2f} s. Final max|u|={np.max(np.abs(u)):.6f}, L2_norm={np.linalg.norm(u):.6f}")
    return x, times, U


# -- Plotting utilities --------------------------------------------------

def plot_xt_contour(x, times, U, fname='u_xt_contour.png'):
    """Contour plot of |u(x,t)| on x-t plane using 'RdBu_r' colormap and save to file."""
    T, X = np.meshgrid(times, x, indexing='xy')  # times along columns, x along rows
    # U has shape (nt+1, nx) where rows correspond to time; transpose for mesh
    UU = np.abs(U)  # magnitude
    fig, ax = plt.subplots(figsize=(8, 4.5))
    # pcolormesh expects X and T shapes matching UU.T or UU depending on arrangement
    pcm = ax.pcolormesh(x, times, UU, shading='auto', cmap='RdBu_r')
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_title('Velocity magnitude |u(x,t)|')
    fig.colorbar(pcm, ax=ax, label='|u|')
    fig.tight_layout()
    fig.savefig(fname, dpi=200)
    plt.close(fig)


def plot_profiles_at_times(x, times, U, targets=(0.2, 0.4, 0.6, 0.8, 1.0), fname='u_profiles.png'):
    """Overlay velocity profiles at specified times in a single figure and save."""
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for tt in targets:
        # find nearest time index
        idx = np.argmin(np.abs(times - tt))
        ax.plot(x, U[idx, :], label=f't={times[idx]:.3f}')
    ax.set_xlabel('x')
    ax.set_ylabel('u(x,t)')
    ax.set_title('Velocity profiles at selected times')
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(fname, dpi=200)
    plt.close(fig)


# -- Entrypoint ----------------------------------------------------------
if __name__ == '__main__':
    # solver configuration
    NX = 401  # number of grid points including endpoints
    CFL = 0.4
    T_FINAL = 1.0

    x, times, U = burgers_1d_imex(nx=NX, CFL=CFL, t_final=T_FINAL)

    # Save plots
    plot_xt_contour(x, times, U, fname='u_xt_contour.png')
    plot_profiles_at_times(x, times, U, targets=(0.2, 0.4, 0.6, 0.8, 1.0), fname='u_profiles.png')

    # Print final quantitative diagnostics
    final_u = U[-1, :]
    print('\nDiagnostics:')
    print(f'visc = {visc:.6e} (global)')
    print(f'Final time t = {times[-1]:.6f}')
    print(f'max(u) = {final_u.max():.6e}, min(u) = {final_u.min():.6e}')
    print(f'L2 norm of final u = {np.linalg.norm(final_u):.6e}')
```



### runtime_outputs

#### Output block1

nx=401, dx=5.000e-03, initial_max|u|=1.000000, dt=2.000000e-03, nsteps=500
step 50/500, t=0.1000, max|u|=0.994430
step 100/500, t=0.2000, max|u|=0.988944
step 150/500, t=0.3000, max|u|=0.983477
step 200/500, t=0.4000, max|u|=0.977911
step 250/500, t=0.5000, max|u|=0.970372
step 300/500, t=0.6000, max|u|=0.932213
step 350/500, t=0.7000, max|u|=0.873955
step 400/500, t=0.8000, max|u|=0.813537
step 450/500, t=0.9000, max|u|=0.756709
step 500/500, t=1.0000, max|u|=0.704984
Completed in 0.13 s. Final max|u|=0.704984, L2_norm=8.470062

Diagnostics:
visc = 3.183099e-03 (global)
Final time t = 1.000000
max(u) = 7.049842e-01, min(u) = -7.049842e-01
L2 norm of final u = 8.470062e+00




### review_decision
accept

### review_comments
	Current Stage [A/1]
1. MINOR

Category: accuracy
Issue: Time-splitting approach (explicit RK3 for convection followed by a single Crank–Nicolson diffusion step) is a Lie (sequential) splitting and is only first-order accurate in time overall.

Recommendation: If you need higher temporal accuracy, replace the Lie splitting with a true IMEX Runge–Kutta scheme (additive RK) or use Strang splitting (convect half-step → diffusion full-step → convect half-step) to achieve second-order accuracy. If first-order splitting is acceptable for your purposes, document this choice in the code/comments so users are aware.


2. MINOR

Category: structure
Issue: A few small code cleanliness issues: precompute_tridiagonal_factor returns b (diag) which is never used later, and plot_xt_contour creates meshgrid(T,X) but then does not use T,X.

Recommendation: Remove or document unused returns/variables (either stop returning b or use it) and delete the unused meshgrid call to avoid confusion. Consider adding a short docstring note in precompute_tridiagonal_factor that the function assumes constant tri-diagonal coefficients.


3. MINOR

Category: runtime
Issue: Time step dt is computed from the initial maximum speed only and then held constant for the whole run. If the solution were to develop steeper features that increase |u|, the CFL condition for the explicit convection step could be violated.

Recommendation: Either recompute dt dynamically from the current max(|u|) (and rebuild the tridiagonal factorization when dt changes) or choose a conservative CFL factor built from a proven upper bound on |u|. Alternatively, monitor CFL during the time loop and assert or adapt if it becomes dangerously large.


4. MINOR

Category: accuracy
Issue: No systematic verification or convergence checks are included (grid refinement, order verification, or comparison to an available reference).

Recommendation: Add short regression/validation tests: compute solutions with two or three grid resolutions and check norms (L2, max) converge at the expected rates; where possible compare with an analytical or highly resolved reference solution (e.g., Cole–Hopf transform result or very fine-grid solution) to quantify discretization error.


5. MINOR

Category: structure
Issue: Plotting uses absolute value with a diverging colormap (RdBu_r) — visually OK but possibly misleading since sign information is lost in the contour plot.

Recommendation: Either plot u (signed) with a diverging map (so positive/negative regions are visible) or explicitly state in the figure caption/filename that the contour is of |u|. If you intend to emphasize magnitude only, consider a sequential colormap (e.g., viridis) which is more conventional for magnitudes.





### iteration_count
1

### rev_image_description
	Current Stage [A/1]




