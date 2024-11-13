import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from Config import L, ny, tol, maxiter, H, tf, u_in, rho, Re

class Fluid:

    def __init__(self, config_file='Config.py'):

        # Domain
        self.H = H
        self.L = L  # Physical length of the domain.
        # Mesh
        self.ny = ny # Ensure cells are square
        self.nx = self.L*ny
        self.dy = H/self.ny
        self.dx = H*L/self.nx # They have to be equal.
        if self.dx != self.dy:
            raise SystemExit("Error: Cells are not square. Simulation terminated.")
        # Physical properties
        self.rho = rho
        self.Re = Re # Kinematic viscosity.
        self.u_in = u_in
        self.nu = (self.u_in * self.H) / self.Re
        self.tf = tf # final simulation time
        print (f'Running simulation with nu={self.nu}, Re={self.Re} and inlet vel={self.u_in}')
        # Numerical
        self.cfl = 0.1 # CFL number
        self.tol = tol  # Poisson solver tolerance threshold value.
        self.maxiter = maxiter # Max number of iterations on the Poisson solver.

        # Turbine parameters
        ny_turb = 6
        # D_turb = ny_turb * self.dy  # Cells span by the turbine's diamater.
        # A = 3.14 * 0.25 * D_turb**2 # Swept area.
        Ct = 0.75 # Thrust coefficient.

        # Need to add force per unit mass (N/kg) to the velocity transport equation (du/dt)
        # F = - 0.5 * rho * Ct * A * u_in**2 in Newtons (N)
        # Specific force per grid cell is Force / Volume of Action / fluid density
        # 'specific F' = F / (A * dx) / rho
        # 'specific F' = -0.5 * Ct * u_in**2 / dx
        self.F = -0.5 * Ct * self.u_in**2 / self.dx

        # Initialize arrays.
        self.u = np.ones((self.ny+2, self.nx+1))*self.u_in  # Set initial velocity to u_in everywhere.
        self.v = np.zeros((self.ny+1, self.nx+2))  # First, fill the entire array with zeros
        self.p = np.zeros((self.ny+2, self.nx+2)) + 1e-20 # Pressure field. Zero initialization creates trouble.

        # Add a turbine.
        self.AI = np.zeros_like(self.u) # Init force array.
        mid_y, mid_x = self.ny // 2, self.nx // 2 # In the middle of the domain.
        # Assign force to corresponding cells.
        self.AI[(mid_y - ny_turb//2):(mid_y + ny_turb//2),mid_x] = self.F

    def interpolate(self, u, v):
        u_avg = (u[:-1, :-1] + u[1:, :-1] + u[1:, 1:] + u[:-1, 1:])/4
        v_avg = (v[:-1, :-1] + v[1:, :-1] + v[1:, 1:] + v[:-1, 1:])/4

        return u_avg, v_avg

    def centers_and_corners(self, u, v):
        u_mid = 0.5*(u[:, 1:] + u[:, :-1])
        v_mid = 0.5*(v[1:, :] + v[:-1, :])

        u_cor = 0.5*(u[1:, :] + u[:-1, :])
        v_cor = 0.5*(v[:, 1:] + v[:, :-1])

        # # Stacking arrays is very very slow. If you really want to make
        # # arrays have a different shape, do this instead...
        # u_avg = np.zeros_like(v) # shape [n+1, n+2]
        # v_avg = np.zeros_like(u) # shape [n+2, n+1]
        # # the u_avg above is shape [n+1, n], 'missing' two columns, so...
        # u_avg[:, 1:-1] = (u[:-1, :-1] + u[1:, :-1] + u[1:, 1:] + u[:-1, 1:])/4
        # # the v_avg above is shape [n, n+1], 'missing' two rows, so...
        # v_avg[1:-1, :] = (v[:-1, :-1] + v[1:, :-1] + v[1:, 1:] + v[:-1, 1:])/4

        return u_mid, v_mid, u_cor, v_cor

    def dudx (self, u):
        # (u_i+1 - u_i-1)/(2*dx).
        return (u[1:-1, 2:] - u[1:-1, :-2])/(2*self.dx)

    def dudy (self, u):
        # (u_j+1 - u_j-1)/(2*dy).
        return (u[2:, 1:-1] - u[:-2, 1:-1])/(2*self.dy)

    def dpdx(self, p):
        # (p_i+1 - p_i)/dx.
        return (p[1:-1, 2:-1] - p[1:-1, 1:-2])/self.dx

    def dpdy(self, p):
        # (p_j+1 - p_j)/dy.
        return (p[2:-1, 1:-1] - p[1:-2, 1:-1])/self.dy

    def d2udx2 (self, u):
        # (u_i+1 - 2u_i + u_i-1)/(dx**2)
        return (u[1:-1, 2:] - 2 * u[1:-1, 1:-1] + u[1:-1, :-2])/(self.dx**2)

    def d2udy2 (self, u):
        # (u_j+1 - 2u_j + u_j-1)/(dy**2)
        return (u[2:, 1:-1] - 2 * u[1:-1, 1:-1] + u[0:-2, 1:-1])/(self.dy**2)

    def set_BCs(self, u, v):
        # Top wall. ZeroGradient BCs (neumann)
        v[-1, :] = v[-2, :] # dv/dy=0 at half cell from boundary. Backward difference..
        u[-1, :] = u[-2, :] # du/dy=0 at boundary.
        # Bottom wall. ZeroGradient BC (neumann)
        u[0, :] = u[1, :] # du/dy=0 at boundary.
        v[0, :] = v[1, :] # dv/dy=0 at half cell from boundary. Forward difference..
        # Inlet.
        u[:, 0] = u_in # To set u=u_in.
        v[:, 0] = -v[:, 1] # v=0
        # Outlet. Zero Gradient
        u[:, -1] = u[:, -2] # du/dx = 0 at half cell from boundary. Backward difference..
        v[:, -1] = v[:, -2] # dv/dx = 0 at outlet.

        return u, v

    def advect_u_alt(self, u, v, dt):
        u_mid = 0.5*(u[:, 1:] + u[:, :-1]) # u at cell centers
        u_cor = 0.5*(u[1:, :] + u[:-1, :]) # u at cell corners
        v_cor = 0.5*(v[:, 1:] + v[:, :-1]) # v at cell corners

        duudx = (u_mid[1:-1, 1:]**2 - u_mid[1:-1, :-1]**2)/self.dx
        duvdy = (u_cor[1:, 1:-1]*v_cor[1:, 1:-1] - u_cor[:-1, 1:-1]*v_cor[:-1, 1:-1])/self.dy
        Laplu = self.nu * (self.d2udx2(u) + self.d2udy2(u))

        u[1:-1, 1:-1] += dt * (Laplu - duudx - duvdy)

        return u

    def advect_u_alt_AI(self, u, v, dt):
        'This function adds an extra term to the momentum equation'
        u_mid = 0.5*(u[:, 1:] + u[:, :-1]) # u at cell centers
        u_cor = 0.5*(u[1:, :] + u[:-1, :]) # u at cell corners
        v_cor = 0.5*(v[:, 1:] + v[:, :-1]) # v at cell corners

        duudx = (u_mid[1:-1, 1:]**2 - u_mid[1:-1, :-1]**2)/self.dx
        duvdy = (u_cor[1:, 1:-1]*v_cor[1:, 1:-1] - u_cor[:-1, 1:-1]*v_cor[:-1, 1:-1])/self.dy
        Laplu = self.nu * (self.d2udx2(u) + self.d2udy2(u))
        axial_induction = self.AI[1:-1, 1:-1]  # Add body force term

        u[1:-1, 1:-1] += dt * (Laplu - duudx - duvdy + axial_induction)

        return u

    def correct_u(self, u, p, dt):
        u[1:-1, 1:-1] -= self.dpdx(p) * dt / self.rho
        return u

    def advect_v_alt(self, v, u, dt):
        v_mid = 0.5*(v[1:, :] + v[:-1, :]) # v at cell centers
        u_cor = 0.5*(u[1:, :] + u[:-1, :]) # u at cell corners
        v_cor = 0.5*(v[:, 1:] + v[:, :-1]) # v at cell corners

        dvvdy = (v_mid[1:, 1:-1]**2 - v_mid[:-1, 1:-1]**2)/self.dy
        duvdx = (u_cor[1:-1, 1:]*v_cor[1:-1, 1:] - u_cor[1:-1, :-1]*v_cor[1:-1, :-1])/self.dx
        Laplv = self.nu * (self.d2udx2(v) + self.d2udy2(v))

        v[1:-1, 1:-1] +=  dt * (Laplv - dvvdy - duvdx)
        return v

    def correct_v(self, v, p, dt):
        v[1:-1, 1:-1] -= self.dpdy(p) * dt / self.rho
        return v

    def unpslit_euler_b(self, u, v, u_avg, v_avg, dt):
        # Central difference (1 cell).
        du_dx = (u[1:-1, 1:] - u[1:-1, :-1]) / self.dx # (u_j+1,i+1-u_j+1,i)/dx, central difference (1cell)
        dv_dy = (v[1:, 1:-1] - v[:-1, 1:-1]) / self.dy #  v_j+1,i+1-v_j,i+1/dy, central difference (1cell)
        du_dy = (u_avg[1:, :] - u_avg[:-1, :]) / self.dy # (u_avg_j+1,i-u_avg_j,i)/dy, central difference (1cell)
        dv_dx = (v_avg[:, 1:] - v_avg[:, :-1]) / self.dx # (v_avg_j,i+1-v_avg_j,i)/dx, central difference (1cell)
        b = (du_dx + dv_dy)/dt - du_dx**2 - 2*du_dy*dv_dx - dv_dy**2

        # Colin introduced an error by having b get multiplied by rho twice,
        # once here, once in pressure_poison. He deleted the extra rho here.
        # Since rho = 1, this never actually showed up as a mistake.

        return b

    def split_chorin_b(self, u, v, dt):
        # Central difference (1 cell).
        du_dx = (u[1:-1, 1:] - u[1:-1, :-1]) / self.dx
        dv_dy = (v[1:, 1:-1] - v[:-1, 1:-1]) / self.dy

        # return divergence divided by dt
        return (du_dx + dv_dy) / dt

    def pressure_poisson(self, p, b):
        """Solve the Poisson equation using Jabobi's method.
        """
        err = np.inf # Initialize huge error.
        nit = 0 # Reset num iterations.
        pcoef = 0.5 / (self.dx**2 + self.dy**2) # Simplifies code
        b *= self.rho * self.dx**2 * self.dy**2 / (2*(self.dx**2 + self.dy**2))

        while err > tol and nit < maxiter:
            pn = p.copy()

            p[1:-1, 1:-1] = (pcoef * ((pn[1:-1, 2:] + pn[1:-1, :-2])*self.dy**2
                            + (pn[2:, 1:-1] + pn[:-2, 1:-1])*self.dx**2) - b)

            # BCs. Openfield.
            p[:, 0] = p[:, 1] # dp/dx=0 at x=0.
            p[:, -1] = -p[:, -2] # p = 0 at x = L.
            p[0, :] = p[1, :]   # dp/dy = 0 at y = 0.
            p[-1, :] = p[-2, :] # dp/dx = 0 at y = 2.

            err = np.mean((p[1:-1, 1:-1] - pn[1:-1, 1:-1])**2)**0.5
            nit += 1

        #print(f'p iters: {nit}')

        return p

    def chorin_projection_step(self, u, v, p, dt):
        """1) Advection-Diffusion.
            2) Solve Poisson.
            3) Correct velocities (add p term).
        """
        # update BC
        u, v = self.set_BCs(u, v)

        # use the convective-form u * du/dx...
        # u_avg, v_avg = interpolate(u, v)
        # u = advect_u(u, v_avg, p, dx, dy, dt, nu, rho)
        # v = advect_v(v, u_avg, p, dx, dy, dt, nu, rho)

        # ...or the conservative-form d(u*u)dx
        u = self.advect_u_alt_AI(u, v, dt,) # Include wind turbine.
        #u = self.advect_u_alt(u, v, dt,)
        v = self.advect_v_alt(v, u, dt)

        # update BC with new intermediate velocity
        u, v = self.set_BCs(u, v)

        # then compute b and p with correct u/v BCs
        b = self.split_chorin_b(u, v, dt) # Calculate b term.
        p = self.pressure_poisson(p, b) # Solve Poisson eq.

        # then correct the velocity with pressure
        u = self.correct_u(u, p, dt)
        v = self.correct_v(v, p, dt)

        return u, v, p

    def solve(self, u, v, p):
        t = 0.0
        n = 0
        Co = 0.5 # Diffusive courant number
        dt_diff = Co * (1/self.dx**2 + 1/self.dy**2)**-1 / (2 * self.nu)
        #next_save_time = 0.2  # Store u array every 5 seconds.
        centerline_data = []  # To store centerline values over time.
        u_point_story = [] # Initialize array of u values at a certain point at different moments.

        while t < tf:
            # Compute time step.
            # Advective Courant number.
            u_max = np.max(np.abs(u[1:-1, 1:-1])) + 1.0e-20
            v_max = np.max(np.abs(v[1:-1, 1:-1])) + 1.0e-20
            dt_adv = self.cfl / (u_max/self.dx + v_max/self.dy)
            dt = min(dt_adv, dt_diff) # Choose conservative.
            u, v, p = self.chorin_projection_step(u, v, p, dt)
            # Extract horizontal centerline profiles.
            # centerline_p = p[p.shape[0] // 2, :]
            # centerline_u = u[u.shape[0] // 2, :]

            if n % 100 == 0:
                print (f'Step: {n}, t = {t:0.3e}, dt = {dt:0.3e}')

            t += dt
            n += 1

        return u, v, p
        #return u, v, p, centerline_p, centerline_u