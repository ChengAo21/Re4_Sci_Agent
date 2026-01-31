# problem statements
"""
Note:
Compared to the "clean" and "physics-centric" prompts presented in the paper for evaluating 
autonomous reasoning, the prompts in this repository include additional engineering-level "Hints." 
These hints (e.g., library versioning constraints, and visualization specifics) 
are designed to ensure operational robustness and stability across diverse local runtime environments 
without compromising the core scientific computing challenges.
"""


burgers_problem = r"""
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
Plot the contour of the velocity magnitude using 'RdBu_r' colormap and velocity profile at timesteps t=0.2,0.4,0.6,0.8,1.0 in only one figure.

[HINTS]
Define visc = 1.0 / Re at the GLOBAL scope (outside main), and strictly use visc instead of nu throughout the code to avoid confusion with u.
Ensure all constants are passed explicitly to functions or defined globally to prevent NameError.
Just save figs do not use plt.show() in the code.
Do not use 'np.trapz' as it is removed in NumPy 2.0, implement the integration manually.
"""


cavity_flow_problem = r"""
The PDEs of 2-D steady incompressible Navier-Stokes equations is given by:
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
where $\alpha$ is 2. Implement a stable and efficient method to solve this problem.
Plot the contour of the velocity magnitude using the 'RdBu_r' colormap, overlaid with streamlines, and the convergence history in one figure.
Implement reasonable acceleration strategies to reduce computational cost.
Just save figs do not use plt.show() in the code.

[HINTS]:
Define visc = 1.0 / Re at the GLOBAL scope (outside main), and strictly use visc instead of nu throughout the code to avoid confusion with u.
Ensure all constants are passed explicitly to functions or defined globally to prevent NameError.
Never use the keyword argument 'tol' in any SciPy solver. Always use 'atol'. 
Print concise progress information ONLY every 10% of total steps.
Focus on generating runnable code with reasonable physical results, do not over-optimize.
Do NOT conclude [not converged] solely because max_iters is reached, assess convergence based on whether residuals are still decreasing. 
"""


sod_shock_problem = r"""
The Euler equations for compressible flow in 1-D are given by:
\begin{cases}
\frac{\partial \rho}{\partial t} + \frac{\partial (\rho u)}{\partial x} = 0, \\
\frac{\partial (\rho u)}{\partial t} + \frac{\partial (\rho u^2 + p)}{\partial x} = 0, \\
\frac{\partial (\rho E)}{\partial t} + \frac{\partial (\rho E u + p u)}{\partial x} = 0,
\end{cases}
in the interval \( x \in [0,1] \), where \( E = \frac{1}{2} u^2 + \frac{p}{(\gamma - 1)\rho} \), \( \gamma = 1.4 \).
At \( t = 0 \), the initial condition is
(\rho, u, p) = 
\begin{cases}
(1.0, 0.0, 1.0), & 0 < x \leq 0.5, \\
(0.125, 0.0, 0.1), & 0.5 \leq x < 1,
\end{cases}

Implement a stable and efficient method to solve this problem.
Plot the density, velocity, and pressure at t=0.1,0.2,0.3 in three figures, respectively.
Just save figs do not use plt.show() in the code.

[HINTS]:
Strictly use visc instead of nu throughout the code to avoid confusion with u.
Ensure all constants are passed explicitly to functions or defined globally to prevent NameError.
Print diagnostic information that helps assess the physical correctness of the solution, such as shock position and the absence of spurious oscillations.
Global conservation outputs are intended for qualitative reference only and should not be treated as strict acceptance criteria.
Do not use 'np.trapz' as it is removed in NumPy 2.0, implement the integration manually.
"""


longTm_NS_problem = r"""
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
The initial condition is:
\(u(x, y, 0) = v(x, y, 0) = 0, (x,y) \in \Omega\)\\

Implement a stable and efficient method to solve this problem.
Implement reasonable acceleration strategies to reduce computational cost.
Simulate until t=1.0. Plot contours of u, v, and p at the final step in one figure using 'RdBu_r' colormap.
Just save figs do not use plt.show() in the code.

[HINTS]:
Define visc = 1.0 / Re at the GLOBAL scope (outside main), and strictly use visc instead of nu throughout the code to avoid confusion with u.
Ensure all constants are passed explicitly to functions or defined globally to prevent NameError.
Never use the keyword argument 'tol' in any SciPy solver. Always use 'atol'. 
Do not use 'np.trapz' as it is removed in NumPy 2.0.
Use indexing='ij' for meshgrid and order='F' for all ravel/reshape operations to ensure consistent shapes and correct coordinate alignment.
Use Chorin's projection method on a Staggered (MAC) Grid using the upwind scheme to ensure stability.
Focus on generating runnable code with reasonable physical results, do not over-optimize.
Print concise progress information ONLY every 10% of total steps.
"""


poison_problem = r"""
The PDE of Poisson 2D equation is given by:
\begin{cases}
-\Delta u = 0 \\
\end{cases}
The domain is a rectangle minus four circles \(\Omega = \Omega_{\text{rec}} \setminus R_i\) where \(\Omega_{\text{rec}} = [-0.5, 0.5]^2\) is the rectangle and \(R_i\) denotes four circle areas:
\begin{cases}
R_1 = \left\{ (x,y) : (x - 0.3)^2 + (y - 0.3)^2 \leq 0.1^2 \right\} \\
R_2 = \left\{ (x,y) : (x + 0.3)^2 + (y - 0.3)^2 \leq 0.1^2 \right\} \\
R_3 = \left\{ (x,y) : (x - 0.3)^2 + (y + 0.3)^2 \leq 0.1^2 \right\} \\
R_4 = \left\{ (x,y) : (x + 0.3)^2 + (y + 0.3)^2 \leq 0.1^2 \right\} \\
\end{cases}

The boundary condition is
\begin{cases}
u = 0,\, x \in \partial R_i \\
u = 1,\, x \in \partial \Omega_{\text{rec}} \\
\end{cases}

Implement a stable and efficient method to solve this problem.
Plot the contours of the velocity magnitude in one figure using 'RdBu_r' colormap, and mark the circles in the plot.
Just save figs do not use plt.show() in the code.

[HINTS]:
Strictly use visc instead of nu throughout the code to avoid confusion with u.
Ensure all constants are passed explicitly to functions or defined globally to prevent NameError.
Address the geometric challenges of curved boundaries using appropriate techniques (e.g., Immersed Boundary Method).
Print necessary solving information to facilitate the reliability check of the solution. 
"""


helmholtz_problem = r"""
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
Plot the contour of u in one figure using 'RdBu_r' colormap, and mark the circles in the plot.
Just save figs do not use plt.show() in the code.

[HINTS]:
Strictly use visc instead of nu throughout the code to avoid confusion with u.
Ensure all constants are passed explicitly to functions or defined globally to prevent NameError.
Use 'atol' to avoid TypeError: gmres()/cg() got an unexpected keyword argument 'tol'.
Address the geometric challenges of curved boundaries using appropriate techniques (e.g., Immersed Boundary Method).
Print necessary solving information to facilitate the reliability check of the solution.
"""


cylinder_problem = r"""
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
"""


hilbert_mat_problem = r"""
Implement various appropriate methods from scratch to solve the system of linear algebra equations $H_n \mathbf{x} = \mathbf{b}$ accurately and efficiently, 
where $H_n$ is an $n \times n$ matrix defined by $H_n(i,j) = \frac{1}{i + j - 1}$, and the vector $\mathbf{b}$ is taken to ensure that the system admits an exact solution $\mathbf{x}^* = (x_i)_{n \times 1} = (1)_{n \times 1}$.

Compare the $L_\infty$ error of the numerical results with the exact solution $\mathbf{x}^*$ for $n = 5, 10, 15, 20, 25$.
Plot the $L_\infty$ error versus $n$ in one figure with linear scale for both axes. Just save figs do not use plt.show() in the code.

[HINTS]:
Ensure the error is below $1e-2$ for all tested $n$ values.
Print necessary solving information to facilitate the reliability check of the solution.
"""


Keyhole_problem = r"""
Please read the CSV data file from path: ./dataset_keyhole.csv and consider the data in columns 3, 4, 5, 6, 7, 8 and 11.
These columns correspond to seven physical quantities respectively:
the effective laser power (eta P), the laser scan speed (V_s), the laser beam radius (r_0), the thermal diffusivity (alpha), the material density (rho), the heat capacity (C_p), and the difference between melting and ambient temperatures (T_l-T_0).

Based on dimensional analysis and using the data of these physical quantities in the file, please identify the optimal dimensionless quantity formed by combining these parameters, which exhibits the highest coefficient of determination (R-squared) and thus dominates the variation in the keyhole aspect ratio e*.

Implement a roubust and reliable method from scratch for this data-driven dimensional analysis. 
Ensure that the the resulting dimensionless exponents are normalized by V_s and that the exponents of physical quantities be integers or rational fractions with absolute values not exceeding 3.

[HINTS]:
Plot one figure with linear scale for both axes to visualize the correlation between the optimal dimensionless quantity and e*.
Just save figs do not use plt.show() in the code.
To read the data file, just use df = pd.read_csv(file_path) without any additional parameters.
"""
