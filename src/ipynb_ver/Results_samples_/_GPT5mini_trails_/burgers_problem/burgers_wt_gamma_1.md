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
Plot the contour of the velocity magnitude using 'RdBu_r' colormap) and velocity profile at timesteps t=0.2,0.4,0.6,0.8,1.0 in only one figure.
Just save figs do not use plt.show() in the code.
Do not use 'np.trapz' as it is removed in NumPy 2.0, implement the integration manually.


### expanded_prob
We must solve the 1D viscous Burgers PDE u_t + u u_x - visc u_{xx} = 0 on x in [-1,1], t in [0,1], with visc = 0.01/pi (but the codebase should use a global variable named visc defined as 1.0/Re and then used everywhere). Initial condition u(x,0) = -sin(pi x) and homogeneous Dirichlet boundary conditions u(-1,t) = u(1,t) = 0 are prescribed. The equation combines nonlinear advection (u u_x) and linear diffusion (visc u_{xx}), so solution behavior depends strongly on visc; for small visc steep gradients form and require careful treatment to avoid nonphysical oscillations or excessive numerical dissipation. Numerical and implementation constraints given in the problem statement: (i) use the global variable visc (avoid naming conflicts with u), (ii) pass constants explicitly to functions or declare them globally to prevent NameError, (iii) produce a single figure that contains (a) a contour/colormap of velocity magnitude using colormap 'RdBu_r' and (b) velocity profiles at t = 0.2, 0.4, 0.6, 0.8, 1.0, (iv) save figures only (no plt.show()), and (v) do not use np.trapz — implement integration manually (composite trapezoid or Simpson constructed by explicit loops or vector ops). Primary mathematical and numerical challenges to address: 
- Nonlinear convection can produce steep fronts and near-discontinuities; low-dissipation schemes risk Gibbs oscillations while high-dissipation schemes smear features. 
- Stability: explicit time-stepping for convection and diffusion is constrained by CFL-like conditions: dt <= C_dx/ max|u| (for convection) and dt <= C_diff * dx^2/visc (for diffusion). For small visc the diffusion constraint can be very strict. 
- Accuracy order: we want at least second-order accuracy in space and time for moderate resolutions; higher-order methods (WENO, spectral) give better resolution but at higher complexity and care for aliasing and boundary handling. 
- Boundary conditions are Dirichlet at both ends — if using spectral or transform-based methods choose bases that naturally satisfy Dirichlet (e.g., sine series) or apply collocation methods with boundary row modifications. 
- Efficient linear solvers: implicit treatment of diffusion leads to tridiagonal systems solvable efficiently with Thomas algorithm; implicit treatment of nonlinear terms requires Newton/fixed-point iterations and is more expensive. 
- Implementation hygiene: use a single global visc, pass domain size, dx, dt, final time, arrays, and plotting parameters explicitly to routines; implement integration formula manually (explicit composite rule). 
- Plotting requirement: combine contour of u(x,t) (or |u| if desired) with overlaid profiles at requested times in one saved figure using 'RdBu_r'.

### solution_plans
	Current Stage [A/2]
solu_name='Finite-Difference IMEX with Conservative Upwind/WENO' content="Governing idea:\nUse a conservative finite-difference spatial discretization with an IMEX (implicit-explicit) time integrator: treat nonlinear convective flux explicitly (using a stable high-resolution flux like second-order upwind or WENO) and treat viscous diffusion implicitly (Crank–Nicolson or implicit second-order backward differentiation) to remove the diffusion stability constraint. This yields a tridiagonal linear system at each step (cheap Thomas solver). Ensure global variable 'visc' is used and all constants are passed into functions.\nAlgorithmic steps (numbered):\n1) Grid and globals: define Nx, x-grid in [-1,1], dx; set global visc = 1.0/Re outside main and refer to it in routines; set times t0=0, tf=1.0 and choose dt according to convection CFL (dt <= CFLc*dx/max|u_init|) but allow larger dt because diffusion is implicit; pass nx, dx, dt, visc, tspan to time integrator functions explicitly.\n2) Initial & BC: initialize u(x,0) = -sin(pi x) on interior nodes; enforce Dirichlet boundaries u(-1)=u(1)=0 at every step (use ghost nodes or simply keep boundary indices fixed).\n3) Spatial discretization of convection: implement conservative flux f(u)=0.5*u^2. For moderate accuracy use a second-order central + flux limiter (MUSCL) or a compact 3rd-order ENO/WENO-3 scheme. Compute numerical fluxes at interfaces (use local Lax–Friedrichs flux-splitting for stability) so the convective term is approximated by -(F_{i+1/2}-F_{i-1/2})/dx.\n4) Spatial discretization of diffusion: use second-order central finite difference for u_{xx} on interior nodes; the discretization yields a tridiagonal matrix A representing visc * u_{xx}.\n5) Time stepping (IMEX Crank–Nicolson/explicit RK): use a two-stage IMEX or explicit RK3 for convective term and CN for diffusion. Example two-step IMEX: (I - 0.5*dt*visc*L) u^{n+1} = (I + 0.5*dt*visc*L) u^n - dt * conv_explicit(u^n), where L is the second-difference operator; implement boundary conditions by modifying first/last equations. Solve the tridiagonal system with Thomas algorithm (O(Nx)).\n6) Nonlinear stability measures: impose a CFL condition for the explicit convective part: dt <= CFLc * dx / max_i |u_i|; enforce a maximum dt if computed dt violates safety limits. Recompute dt adaptively if desired.\n7) Diagnostics & manual integration: implement manual composite trapezoid integration routine (do not use np.trapz) for any integral diagnostics (e.g., L2 norm). Example: integral = sum_{i=0}^{N-2} 0.5*(f_i+f_{i+1})*dx using explicit vector operations or loop.\n8) Output & plotting: collect u(x,t) snapshots for t in {0.2,0.4,0.6,0.8,1.0} and the full time history to prepare a contour plot in the (x,t) plane. Create a single figure with a contour (imshow or contourf) of u(x,t) over the x-t grid using colormap 'RdBu_r', and in the same figure add an inset or lower subplot that overplots the five profiles (different line styles/legend). Save the figure to PNG using fig.savefig(...). Do not call plt.show().\nStability / accuracy / complexity / efficiency limitations:\n- Accuracy depends on the convective flux discretization; second-order with limiters gives robust non-oscillatory results but smears steep gradients; WENO gives higher resolution but costs more CPU and is more complex to implement. \n- IMEX eliminates the diffusion timestep constraint but the explicit convective step still enforces a CFL dt ~ O(dx/|u|). For very small visc and strong gradients, dt may still be small. \n- The implicit solve each step is O(N) (Thomas) and memory O(N) — very efficient. However, if the convective term were treated implicitly (to allow larger dt), one would need nonlinear solves per step (Newton or fixed point), raising cost significantly. \n- The method is conservative if fluxes are computed consistently; boundary handling and limiter choice can affect global conservation. \n- Implementing WENO/ENO increases implementation complexity (reconstruction, nonlinear weights), and requires careful boundary stencils or ghost-cell treatment. \n- If adaptive dt is used, ensure consistent snapshot times for the requested profile times (interpolate in time if necessary)."

	Current Stage [B/2]
solu_name='Sine-Transform Spectral + ETDRK4 (semi-analytic linear step)' content="Governing idea:\nExploit spectral accuracy in space for smooth parts of the solution by expanding in a sine series (Dirichlet homogeneous BCs) and integrate in time with an exponential time-differencing RK4 (ETDRK4) or an IMEX RK method: treat the linear viscous operator exactly in Fourier/sine space (diagonal) and compute the nonlinear convective term in physical space with dealiasing. Use sine discrete transforms (DST) so boundary conditions are enforced exactly. Precompute ETDRK4 coefficients to speed time stepping and use 2/3 rule for dealiasing.\nAlgorithmic steps (numbered):\n1) Grid and transform choice: choose Nx interior collocation points compatible with sine transform (e.g., x_j = j*dx, j=1..N with dx=2/(N+1) mapping [-1,1]); set global visc explicitly and pass it to routines. Use DST-I/II consistent with chosen collocation. \n2) Represent u(x,t) as series u(x,t)=sum_k a_k(t) sin(k*pi*(x+1)/2) so that u(-1)=u(1)=0 automatically. The Laplacian eigenvalues are simple: -k^2*(pi/2)^2 in this basis scaled by visc.\n3) Linear operator diagonalization: in transform (k) space the diffusion term is visc * lambda_k * a_k where lambda_k is known. Precompute exp(dt*L) and ETDRK4 coefficients (or use an IMEX-RK where diffusion is handled implicitly and requires only diagonal inversion).\n4) Nonlinear term evaluation: at each time step, transform a_k -> u(x) via inverse DST, compute nonlinear convective flux f(u) = 0.5*u^2 in physical space, compute its derivative by spectral differentiation (or transform f(u) back and differentiate in spectral space by multiplying by ik_k) and then transform the nonlinear term back. Apply 2/3 dealiasing: zero out high modes above 2/3*N before inverse transform to avoid aliasing errors.\n5) Time integration: use ETDRK4 (recommended for stiff linear + nonlinearity) with precomputed phi-functions; alternatively use IMEX-RK4 treat linear part exactly and nonlinear explicitly in physical space. ETDRK4 requires only multiplications and transforms — no linear solves. \n6) Adaptive dt consideration: ETDRK4 is stable for larger dt than explicit schemes, but the nonlinear convective CFL still imposes practical dt limits for accuracy; monitor max|u| and reduce dt if nonlinear error grows. Pass all parameters explicitly.\n7) Diagnostics & manual integration: compute integrals (L2 norm or mass) using manual composite trapezoid implemented directly on the collocation grid (explicit sum of element contributions). Keep all constants explicit; do not import trapz.\n8) Output & plotting: store time series at times t=0.2,0.4,..., and create one figure that has a contourf of u(x,t) over the x-t plane with colormap 'RdBu_r' and a second subplot or an overlaid axes plotting u(x,t_i) for the five requested times (labels, legend). Save the figure to disk (no plt.show()).\nStability / accuracy / complexity / efficiency limitations:\n- Spectral (sine) gives exponential convergence for smooth solutions; however, Burgers can develop sharp gradients which produce Gibbs oscillations and require careful dealiasing and possibly spectral filtering. \n- ETDRK4 avoids a linear CFL for diffusion and allows larger dt, but accuracy is still limited by the nonlinear convective dynamics; very steep gradients require small dt to resolve advective timescales. \n- Implementation complexity is higher: correct sine transform conventions, dealiasing, and precomputation of ETDRK4 phi-coefficients are needed and can be a source of bugs. \n- Each time step requires several transforms (DST/iDST) — cost is O(N log N) per transform; with multiple transforms per stage the per-step cost is higher than simple finite-difference + Thomas O(N) implicit solves, but for a given accuracy spectral methods can require fewer grid points. \n- Boundary conditions are enforced naturally, but non-smoothness near boundaries (or mismatch with initial data) can reduce convergence rate. \n- If solution becomes effectively discontinuous in the inviscid limit, spectral methods perform poorly unless supplemented with shock-capturing filters or hybrid approaches."



### technical_spec
	Current Stage [A/2]
We solve 1D viscous Burgers using an IMEX Crank–Nicolson scheme: the nonlinear convection is discretized explicitly with a Rusanov (local Lax-Friedrichs) conservative flux, while diffusion is handled implicitly by Crank–Nicolson leading to a constant tridiagonal linear system solved via the Thomas algorithm. The diffusion tridiagonal matrices are assembled once. The code enforces Dirichlet BCs, stores time-history for contour plotting, and extracts snapshots (with linear interpolation if timestep doesn't coincide). Manual trapezoid integration is provided for diagnostics. The figure saved contains a contour (x vs t) of u with 'RdBu_r' colormap and an overlaid subplot of velocity profiles at requested times.

	Current Stage [B/2]
IMEX solver overview: The code implements a 1D viscous Burgers equation solver using an IMEX Crank–Nicolson (CN) discretization for diffusion and an explicit Rusanov (local Lax–Friedrichs) conservative flux for convection. The diffusion CN linear system is assembled once into tridiagonal coefficients (a, b, c) and solved each time-step with a Thomas algorithm. The convective flux divergence is computed explicitly via numerical fluxes at interfaces and applied as an explicit source in the RHS term; a previous sign bug (applying -dt * conv when conv already represented -d/dx f) is fixed by applying +dt * conv. Time step (dt) is chosen conservatively up front to avoid rebuilding diffusion matrices inside the time-loop: dt is limited by convection (using a safe estimate of max velocity) and by a diffusion guideline dt <= C*dx^2/visc. Defensive checks are added to the Thomas solver and before solving to detect NaNs/Infs. The code stores time-history for contour plotting and linearly interpolates snapshots if timesteps don't coincide. The final figure contains a pcolormesh of u(x,t) with colormap 'RdBu_r' (axes t vs x) and an overlaid subplot of velocity profiles at requested times; the figure is saved (no plt.show()).



### python_codes

#### Script block1:

```python
import numpy as np
import matplotlib.pyplot as plt

# Global viscosity defined as 1.0 / Re (do not rename)
Re = np.pi / 0.01  # so visc = 0.01/pi
visc = 1.0 / Re

# ---------------------- Utility functions ----------------------

def trapz_manual(x, y):
    """Manual composite trapezoid integration over nonuniform grid x."""
    dx = np.diff(x)
    return np.sum(0.5 * (y[:-1] + y[1:]) * dx)


def rusanov_convective_term(u, dx):
    """Compute convective term d/dx (0.5*u^2) using Rusanov flux.
    Returns an array of same size as u with convective term at interior nodes,
    with zeros at boundaries (Dirichlet enforced separately).
    """
    N = u.size
    f = 0.5 * u * u
    # Numerical fluxes at interfaces i+1/2 for i=0..N-2
    # alpha = max(|u_i|, |u_{i+1}|)
    alpha = np.maximum(np.abs(u[:-1]), np.abs(u[1:]))
    F = 0.5 * (f[:-1] + f[1:]) - 0.5 * alpha * (u[1:] - u[:-1])
    conv = np.zeros_like(u)
    # interior nodes i = 1..N-2: conv[i] = -(F[i] - F[i-1]) / dx
    conv[1:-1] = -(F[1:] - F[:-1]) / dx
    # boundaries remain zero (Dirichlet), so conv[0]=conv[-1]=0
    return conv


def build_tridiagonal_CN(Ni, dx, dt, visc):
    """Build tridiagonal coefficients for Crank-Nicolson diffusion system
    for interior unknowns (size Ni = Nx-2). Returns arrays for Thomas solver
    (a lower diag, b main diag, c upper diag) for the left matrix M_left = I - 0.5*dt*A,
    and arrays (bR, offR) to multiply M_right @ u_interior efficiently.
    """
    # A = visc/dx^2 * T , T = tridiag(1,-2,1) (size Ni)
    coef = visc / (dx * dx)
    off = coef  # off-diagonal entries of A
    mainA = -2.0 * coef
    # M_left = I - 0.5*dt*A => mainL = 1 - 0.5*dt*mainA, offL = -0.5*dt*off
    # But mainA is negative, so mainL = 1 + dt*visc/dx^2
    mainL = 1.0 - 0.5 * dt * mainA
    offL = -0.5 * dt * off
    # M_right = I + 0.5*dt*A => mainR = 1 + 0.5*dt*mainA , offR = +0.5*dt*off
    mainR = 1.0 + 0.5 * dt * mainA
    offR = 0.5 * dt * off
    # Create tridiagonal arrays for Thomas solver (a lower diag, b main, c upper)
    b_left = np.full(Ni, mainL)
    a_left = np.full(Ni - 1, offL)  # lower diag
    c_left = np.full(Ni - 1, offL)  # upper diag
    # For RHS multiplication we only need mainR and offR (symmetric)
    b_right = np.full(Ni, mainR)
    off_right = offR
    return a_left, b_left, c_left, b_right, off_right


def thomas_solve(a, b, c, d):
    """Solve tridiagonal system with lower diag a, main diag b, upper diag c.
    a, b, c are 1D arrays with lengths n-1, n, n-1 respectively. d is RHS (n,).
    Returns solution x of length n. In-place-safe implementation.
    """
    n = b.size
    # make copies to avoid modifying inputs
    cp = np.empty(n - 1, dtype=float)
    dp = np.empty(n, dtype=float)
    x = np.empty(n, dtype=float)

    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]
    for i in range(1, n - 1):
        denom = b[i] - a[i - 1] * cp[i - 1]
        cp[i] = c[i] / denom
        dp[i] = (d[i] - a[i - 1] * dp[i - 1]) / denom
    denom = b[-1] - a[-1] * cp[-1]
    dp[-1] = (d[-1] - a[-1] * dp[-2]) / denom

    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


# ---------------------- Main solver ----------------------

def solve_burgers_imex(nx=401, tf=1.0, cfl=0.4):
    # Domain
    x0, x1 = -1.0, 1.0
    x = np.linspace(x0, x1, nx)
    dx = x[1] - x[0]

    # Initial condition u(x,0) = -sin(pi x)
    u0 = -np.sin(np.pi * x)
    # Apply Dirichlet BCs explicitly
    u0[0] = 0.0
    u0[-1] = 0.0

    # Fixed time step computed from initial condition to avoid rebuilding matrices
    umax = np.max(np.abs(u0))
    if umax <= 0:
        umax = 1.0
    dt = cfl * dx / umax
    Nt = int(np.ceil(tf / dt))
    dt = tf / Nt  # adjust dt so that Nt*dt = tf (keeps matrix constant)

    # Prebuild tridiagonal matrices for interior unknowns (size Ni = nx-2)
    Ni = nx - 2
    a_left, b_left, c_left, b_right, off_right = build_tridiagonal_CN(Ni, dx, dt, visc)

    # Time integration arrays
    u = u0.copy()
    u_prev = u.copy()

    # Store history for contour plotting (Nt+1 times including initial)
    history = np.zeros((Nt + 1, nx))
    t_history = np.zeros(Nt + 1)
    history[0, :] = u.copy()
    t_history[0] = 0.0

    # Snapshot times we need profiles for (may not align with timesteps)
    snapshot_times = [0.2, 0.4, 0.6, 0.8, 1.0]
    snapshots = {t_s: None for t_s in snapshot_times}

    t = 0.0
    for n in range(1, Nt + 1):
        t_prev = t
        u_prev[:] = u
        # Compute explicit convective term
        conv = rusanov_convective_term(u_prev, dx)
        # RHS for interior unknowns: M_right @ u_interior - dt * conv_interior
        u_in = u_prev[1:-1]
        # Multiply M_right @ u_interior (tridiagonal multiply)
        # (b_right * u_in) + off_right * (u_in shifted left/right)
        rhs = b_right * u_in.copy()
        rhs[:-1] += off_right * u_in[1:]
        rhs[1:] += off_right * u_in[:-1]
        # subtract dt * convective term on interior
        rhs -= dt * conv[1:-1]

        # Solve for new interior u using Thomas
        u_in_new = thomas_solve(a_left, b_left, c_left, rhs)
        # Compose full u with Dirichlet BCs
        u[0] = 0.0
        u[-1] = 0.0
        u[1:-1] = u_in_new

        t += dt
        history[n, :] = u.copy()
        t_history[n] = t

        # If any snapshot times are crossed, interpolate linearly between u_prev and u
        for t_snap in snapshot_times:
            if snapshots[t_snap] is None and t_prev < t_snap <= t:
                theta = (t_snap - t_prev) / (t - t_prev)
                snapshots[t_snap] = (1.0 - theta) * u_prev + theta * u

    # Compute diagnostics at final time
    u_final = history[-1, :]
    mass = trapz_manual(x, u_final)
    energy = trapz_manual(x, u_final * u_final)
    print(f"Final time t = {t:.6f}")
    print(f"u_min = {u_final.min():.6e}, u_max = {u_final.max():.6e}")
    print(f"Mass (integral u) = {mass:.6e}")
    print(f"L2 energy (integral u^2) = {energy:.6e}")

    return x, t_history, history, snapshots


# ---------------------- Plotting ----------------------

def plot_and_save(x, t_history, history, snapshots, fname='burgers_solution.png'):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), gridspec_kw={'height_ratios': [3, 1]})

    # Top: contour/pcolormesh of u(x,t)
    T, X = np.meshgrid(t_history, x)
    U_plot = history.T  # shape (Nx, Nt+1)
    pcm = ax1.pcolormesh(X, T, U_plot, cmap='RdBu_r', shading='auto')
    cbar = fig.colorbar(pcm, ax=ax1)
    cbar.set_label('u(x,t)')
    ax1.set_xlabel('x')
    ax1.set_ylabel('t')
    ax1.set_title('Burgers solution u(x,t) (RdBu_r)')

    # Bottom: snapshots at requested times
    colors = ['k', 'r', 'b', 'g', 'm']
    for (t_s, col) in zip(sorted(snapshots.keys()), colors):
        u_snap = snapshots[t_s]
        if u_snap is None:
            continue
        ax2.plot(x, u_snap, label=f't={t_s:.1f}', color=col)
    ax2.set_xlabel('x')
    ax2.set_ylabel('u')
    ax2.set_title('Profiles at selected times')
    ax2.legend()
    ax2.grid(True)

    plt.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)


# ---------------------- Main entry ----------------------

if __name__ == '__main__':
    # Solver parameters
    nx = 401  # spatial points (including boundaries)
    tf = 1.0
    cfl = 0.4

    x, t_history, history, snapshots = solve_burgers_imex(nx=nx, tf=tf, cfl=cfl)
    plot_and_save(x, t_history, history, snapshots, fname='burgers_solution.png')
```


#### Script block2:

```python
import numpy as np
import matplotlib.pyplot as plt

# Global viscosity defined as 1.0 / Re (do not rename)
Re = np.pi / 0.01  # so visc = 0.01/pi
visc = 1.0 / Re

# ---------------------- Utility functions ----------------------

def trapz_manual(x, y):
    """Manual composite trapezoid integration over nonuniform grid x."""
    dx = np.diff(x)
    return np.sum(0.5 * (y[:-1] + y[1:]) * dx)


def rusanov_convective_term(u, dx):
    """Compute convective term using Rusanov flux.
    This routine returns conv where conv[i] approximates -d/dx(f) at node i
    (i.e. conv = -d/dx(0.5*u^2)). Boundary entries are zero to respect Dirichlet BCs.
    """
    N = u.size
    f = 0.5 * u * u
    # Numerical fluxes at interfaces i+1/2 for i=0..N-2
    alpha = np.maximum(np.abs(u[:-1]), np.abs(u[1:]))
    F = 0.5 * (f[:-1] + f[1:]) - 0.5 * alpha * (u[1:] - u[:-1])
    conv = np.zeros_like(u)
    # interior nodes i = 1..N-2: conv[i] = -(F[i] - F[i-1]) / dx  => conv ~ -d/dx(f)
    conv[1:-1] = -(F[1:] - F[:-1]) / dx
    return conv


def build_tridiagonal_CN(Ni, dx, dt, visc):
    """Build tridiagonal coefficients for Crank-Nicolson diffusion system
    for interior unknowns (size Ni = Nx-2). Returns arrays for Thomas solver
    (a lower diag, b main diag, c upper diag) for the left matrix M_left = I - 0.5*dt*A,
    and arrays (bR, offR) to multiply M_right @ u_interior efficiently.
    """
    coef = visc / (dx * dx)
    off = coef
    mainA = -2.0 * coef
    # M_left = I - 0.5*dt*A
    mainL = 1.0 - 0.5 * dt * mainA
    offL = -0.5 * dt * off
    # M_right = I + 0.5*dt*A
    mainR = 1.0 + 0.5 * dt * mainA
    offR = 0.5 * dt * off
    b_left = np.full(Ni, mainL)
    a_left = np.full(Ni - 1, offL)  # lower diag
    c_left = np.full(Ni - 1, offL)  # upper diag (symmetric here)
    b_right = np.full(Ni, mainR)
    off_right = offR
    return a_left, b_left, c_left, b_right, off_right


def thomas_solve(a, b, c, d, eps=1e-14):
    """Solve tridiagonal system with lower diag a, main diag b, upper diag c.
    Defensive: guard small pivots and check finite-ness of inputs.
    """
    if not (np.all(np.isfinite(b)) and np.all(np.isfinite(d))):
        raise ValueError("Non-finite values detected in main diagonal or RHS before Thomas solve.")
    n = b.size
    if n == 0:
        return np.array([])
    # copies for safe in-place forward elimination
    cp = np.empty(n - 1, dtype=float)
    dp = np.empty(n, dtype=float)
    x = np.empty(n, dtype=float)

    # first row
    denom = b[0]
    if abs(denom) < eps:
        denom = eps if denom >= 0 else -eps
    cp[0] = c[0] / denom if n > 1 else 0.0
    dp[0] = d[0] / denom

    for i in range(1, n - 1):
        denom = b[i] - a[i - 1] * cp[i - 1]
        if abs(denom) < eps:
            denom = eps if denom >= 0 else -eps
        cp[i] = c[i] / denom
        dp[i] = (d[i] - a[i - 1] * dp[i - 1]) / denom

    if n > 1:
        denom = b[-1] - a[-1] * cp[-1]
        if abs(denom) < eps:
            denom = eps if denom >= 0 else -eps
        dp[-1] = (d[-1] - a[-1] * dp[-2]) / denom
    else:
        # single equation, dp already set
        pass

    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - (cp[i] * x[i + 1])
    if not np.all(np.isfinite(x)):
        raise np.linalg.LinAlgError('Thomas solver produced non-finite solution')
    return x


# ---------------------- Main solver ----------------------

def solve_burgers_imex(nx=401, tf=1.0, cfl=0.4):
    # Domain
    x0, x1 = -1.0, 1.0
    x = np.linspace(x0, x1, nx)
    dx = x[1] - x[0]

    # Initial condition u(x,0) = -sin(pi x)
    u0 = -np.sin(np.pi * x)
    # Apply Dirichlet BCs explicitly
    u0[0] = 0.0
    u0[-1] = 0.0

    # Conservative estimate for max velocity to choose a constant dt
    umax0 = np.max(np.abs(u0))
    umax_bound = max(umax0, 1.0)  # safety floor to avoid tiny denom and overly large dt

    # Diffusion-guided timestep limit (CN is unconditionally stable but may oscillate if dt too large)
    dt_diff = 0.5 * dx * dx / max(visc, 1e-16)

    # Convection-guided timestep
    dt_conv = cfl * dx / umax_bound

    dt_pre = min(dt_conv, dt_diff, tf)
    # adjust dt to divide tf into integer Nt and keep matrices constant
    Nt = int(np.ceil(tf / dt_pre))
    dt = tf / Nt

    # Prebuild tridiagonal matrices for interior unknowns (size Ni = nx-2)
    Ni = nx - 2
    a_left, b_left, c_left, b_right, off_right = build_tridiagonal_CN(Ni, dx, dt, visc)

    # Time integration arrays
    u = u0.copy()
    u_prev = u.copy()

    # Store history for contour plotting (Nt+1 times including initial)
    history = np.zeros((Nt + 1, nx))
    t_history = np.zeros(Nt + 1)
    history[0, :] = u.copy()
    t_history[0] = 0.0

    # Snapshot times we need profiles for (may not align with timesteps)
    snapshot_times = [0.2, 0.4, 0.6, 0.8, 1.0]
    snapshots = {t_s: None for t_s in snapshot_times}

    t = 0.0
    for n in range(1, Nt + 1):
        t_prev = t
        u_prev[:] = u
        # Compute explicit convective term (conv ~ -d/dx(f))
        conv = rusanov_convective_term(u_prev, dx)
        # RHS for interior unknowns: M_right @ u_interior + dt * conv_interior
        u_in = u_prev[1:-1]
        # Multiply M_right @ u_interior (tridiagonal multiply)
        rhs = b_right * u_in.copy()
        if Ni > 1:
            rhs[:-1] += off_right * u_in[1:]
            rhs[1:] += off_right * u_in[:-1]
        # Apply convective contribution. Note: conv already equals -d/dx(f), so add dt*conv.
        rhs += dt * conv[1:-1]

        # Defensive check before Thomas
        if not np.all(np.isfinite(rhs)):
            raise ValueError(f"Non-finite RHS detected at step {n}, t={t_prev:.6f}")

        # Solve for new interior u using Thomas
        u_in_new = thomas_solve(a_left, b_left, c_left, rhs)
        # Compose full u with Dirichlet BCs
        u[0] = 0.0
        u[-1] = 0.0
        u[1:-1] = u_in_new

        t += dt
        history[n, :] = u.copy()
        t_history[n] = t

        # If any snapshot times are crossed, interpolate linearly between u_prev and u
        for t_snap in snapshot_times:
            if snapshots[t_snap] is None and t_prev < t_snap <= t:
                theta = (t_snap - t_prev) / (t - t_prev)
                snapshots[t_snap] = (1.0 - theta) * u_prev + theta * u

    # Compute diagnostics at final time
    u_final = history[-1, :]
    mass = trapz_manual(x, u_final)
    energy = trapz_manual(x, u_final * u_final)
    print(f"Final time t = {t:.6f}")
    print(f"u_min = {u_final.min():.6e}, u_max = {u_final.max():.6e}")
    print(f"Mass (integral u) = {mass:.6e}")
    print(f"L2 energy (integral u^2) = {energy:.6e}")

    return x, t_history, history, snapshots


# ---------------------- Plotting ----------------------

def plot_and_save(x, t_history, history, snapshots, fname='burgers_solution.png'):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), gridspec_kw={'height_ratios': [3, 1]})

    # Top: pcolormesh of u(x,t). Use t on x-axis and x on y-axis to match arrays.
    pcm = ax1.pcolormesh(t_history, x, history.T, cmap='RdBu_r', shading='auto')
    cbar = fig.colorbar(pcm, ax=ax1)
    cbar.set_label('u(x,t)')
    ax1.set_xlabel('t')
    ax1.set_ylabel('x')
    ax1.set_title('Burgers solution u(x,t) (RdBu_r)')

    # Bottom: snapshots at requested times
    colors = ['k', 'r', 'b', 'g', 'm']
    for (t_s, col) in zip(sorted(snapshots.keys()), colors):
        u_snap = snapshots[t_s]
        if u_snap is None:
            continue
        ax2.plot(x, u_snap, label=f't={t_s:.1f}', color=col)
    ax2.set_xlabel('x')
    ax2.set_ylabel('u')
    ax2.set_title('Profiles at selected times')
    ax2.legend()
    ax2.grid(True)

    plt.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)


# ---------------------- Main entry ----------------------

if __name__ == '__main__':
    # Solver parameters
    nx = 401  # spatial points (including boundaries)
    tf = 1.0
    cfl = 0.4

    x, t_history, history, snapshots = solve_burgers_imex(nx=nx, tf=tf, cfl=cfl)
    plot_and_save(x, t_history, history, snapshots, fname='burgers_solution.png')
```



### runtime_outputs

#### Output block1

<string>:22: RuntimeWarning: overflow encountered in multiply
<string>:26: RuntimeWarning: overflow encountered in multiply
<string>:26: RuntimeWarning: invalid value encountered in subtract
<string>:29: RuntimeWarning: invalid value encountered in subtract
<string>:29: RuntimeWarning: overflow encountered in divide
Final time t = 1.000000
u_min = nan, u_max = nan
Mass (integral u) = nan
L2 energy (integral u^2) = nan



#### Output block2

Final time t = 1.000000
u_min = -7.047128e-01, u_max = 7.047128e-01
Mass (integral u) = 1.110223e-16
L2 energy (integral u^2) = 3.580332e-01




### review_decision
accept

### review_comments
	Current Stage [A/2]
1. MAJOR

Category: runtime
Issue: Numerical blow-up produces overflows and NaNs (runtime warnings) and final solution is NaN.

Recommendation: Root cause: a sign error in how the convective term is applied. In rusanov_convective_term() conv is computed as conv = -(F[i] - F[i-1]) / dx, i.e. conv ≈ -d/dx(f). In the IMEX update the RHS requires -dt * (d/dx f). With conv = -d/dx f the correct update is rhs += dt * conv (or equivalently define conv as d/dx f and keep rhs -= dt*conv). Change the line rhs -= dt * conv[1:-1] to rhs += dt * conv[1:-1] (or flip the sign in rusanov_convective_term and keep rhs subtract). After fixing this, rerun and confirm that overflows disappear and that u remains finite.


2. MINOR

Category: structure
Issue: CFL time-step fixed using only the initial maximum |u|; not updated during time integration.

Recommendation: For explicit convective stability recompute umax each time step (or at least periodically) and adapt dt = cfl*dx/umax. If you want to keep a constant matrix for CN diffusion, you can recompute dt only when it changes significantly, or choose dt based on a conservative bound (e.g. umax_bound = max(|u0|, expected growth) ). Add a safety floor to umax to avoid division by zero.


3. MINOR

Category: structure
Issue: Thomas solver can encounter division by small/zero pivots if diagonals become tiny or sign errors propagate NaNs, leading to divide-by-zero or large cp/dp values.

Recommendation: Defensive fixes: (1) check and assert that b (main diagonal) entries are not too small before division (e.g. np.all(np.abs(b) > eps)); (2) use a small epsilon to guard denominators; (3) add nan/inf checks on RHS before calling Thomas and fail early with informative error if values are invalid. After fixing the convective sign bug this is unlikely to trigger, but defensive checks help debugging.


4. MINOR

Category: structure
Issue: Plotting: pcolormesh coordinate arrays are created with meshgrid(t_history, x) as (T,X) but passed to pcolormesh as (X,T,U_plot) which is confusing and the axis labels appear swapped.

Recommendation: Prefer calling X, T = np.meshgrid(x, t_history, indexing='xy') or simply call ax1.pcolormesh(t_history, x, history.T, cmap=...) and set xlabel='t', ylabel='x' to make axes and data correspondence explicit. Verify the orientation visually or via extent. This is not causing the NaNs but improves clarity and avoids transposed plots.


5. MINOR

Category: accuracy
Issue: Crank–Nicolson for diffusion with explicit Rusanov convection may produce oscillations if dt is large relative to visc/dx^2; no explicit control of diffusion-related timestep is used.

Recommendation: Although diffusion is implicit and CN is unconditionally stable, CN is only neutrally stable and can produce non-physical oscillations for large dt. Consider either (a) using backward-Euler (more damping) for the diffusion solve if you want robustness for larger dt, or (b) limit dt also using a diffusion-motivated guideline (dt <= C * dx^2 / visc) when selecting dt. At minimum test convergence and monotonicity by halving dt and/or refining dx.



	Current Stage [B/2]
1. MINOR

Category: accuracy
Issue: Time-step selection is conservative but static: umax is taken only from the initial condition and a floor of 1.0 is enforced. That is arbitrary and may be over- or under-conservative for other initial data.

Recommendation: Either compute umax dynamically and adapt dt (rebuild CN matrices when dt changes), or choose a documented safety floor (explain why 1.0). For robustness consider an adaptive-step strategy (e.g. recompute dt from max|u| every few steps and rebuild matrices if needed) or use sub-stepping for the convective term so the CN matrix remains constant.


2. MINOR

Category: accuracy
Issue: Rusanov (local Lax–Friedrichs) flux used for the convective term is robust but relatively diffusive; the IMEX time-stepping is first-order accurate in time overall because convection is advanced explicitly with forward Euler.

Recommendation: If sharper resolution of steep gradients is desired, consider replacing the explicit time advance of convection with a second-order strong-stability-preserving RK2 (SSP-RK2) or a higher-order spatial reconstructor (MUSCL, WENO). If you adopt RK2 you will need to either sub-step the CN solve or assemble + solve the linear CN system twice per time-step.


3. MINOR

Category: structure
Issue: Boundary handling for the convective flux sets conv[0] and conv[-1] = 0 implicitly by zeroing the conv vector. While this is consistent with the Dirichlet BCs used, it is slightly ad-hoc and may hide boundary-flux inconsistencies for other problems.

Recommendation: For clarity and to avoid surprises, compute interface fluxes using ghost states set to the exact Dirichlet BC values (u_ghost = 0 here). That makes the numerical-flux assembly explicit and easier to generalize to other BCs. Also document this choice in the function docstring.





### iteration_count
2

### rev_image_description
	Current Stage [A/2]


	Current Stage [B/2]




