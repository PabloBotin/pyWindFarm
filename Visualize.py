import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import cm
import numpy as np
from matplotlib.colors import ListedColormap

class Visualize:
    def __init__(self, fluid_instance):
        self.fluid_instance = fluid_instance

    def plot_divergence(u, v, dx, dy, fig, ax):
        """
        Plots the 2D divergence field given the velocity components u and v.
        dudx and dvdy are calculated using central difference.
        """
        # Calculate du/dx and dv/dy using central difference
        du_dx = (u[1:-1, 1:] - u[1:-1, :-1]) / dx  # interior x-derivative
        dv_dy = (v[1:, 1:-1] - v[:-1, 1:-1]) / dy  # interior y-derivative
        
        # Compute divergence for interior points only
        divergence = du_dx + dv_dy

        # Plot divergence for interior points
        cp = ax.imshow(divergence[:, :-1], cmap='seismic', origin='lower', aspect='equal')
        fig.colorbar(cp, ax=ax).set_label('Divergence')
        ax.set_title('2D Divergence Field (Interior Points)')

    def plot_pressure(pressure, fig, ax):
        """
        Plots the 2D pressure field.
        """
        cp = ax.imshow(pressure[1:-1, 1:-1], cmap='plasma', origin='lower',
                    aspect='equal')
        fig.colorbar(cp, ax=ax).set_label('Pressure')
        ax.set_title('2D Pressure Field')

    def plot_u(u, fig, ax):
        """
        Plots the 2D x-velocity field.
        """
        # Plot the velocity field with specified vmin and vmax
        cp = ax.imshow(u[1:-1, 1:-1], cmap='plasma', origin='lower',
                    aspect='equal')
        fig.colorbar(cp, ax=ax).set_label('U-velocity')
        ax.set_title('2D U-velocity Field')

    def plot_u_centerlines(centerline_data, fig, ax):
        """
        Plots the centerline velocity profiles of u at different times.
        
        Parameters:
        - centerline_data: List of tuples (time, centerline array).
        - fig: Figure object.
        - ax: Axes object to plot on.
        """
        # Loop through each time-centerline pair and plot
        for time, centerline in centerline_data:
            ax.plot(centerline, label=f't = {time:.1f} s')
        
        # Customize the plot
        ax.set_xlabel('Horizontal Position (grid index)')
        ax.set_ylabel('U-velocity at Centerline')
        ax.set_title('Centerline Velocity Profiles Over Time')
        #ax.legend(loc="upper right", fontsize="small")
        ax.grid(True)

    def plot_u_point_story(u_values, tf, fig, ax):
        """
        Plots the evolution of u at a specific point as a function of time.

        Parameters:
        - u_values: List or array of recorded u values at each timestep.
        - tf: Final simulation time.
        - fig: Matplotlib figure object.
        - ax: Matplotlib axes object.
        """
        # Create a time array based on the length of u_values and tf
        time_array = np.linspace(0, tf, len(u_values))
        
        # Plot u values in function of time
        ax.plot(time_array, u_values, color='blue', marker='o', linestyle='-', markersize=4)
        
        # Label the plot
        ax.set_xlabel('Time')
        ax.set_ylabel('U-velocity at specified point')
        #ax.set_ylim(0.85, 0.87)  # Focus on centerline to check oscillations. 
        ax.set_title('U-velocity evolution over time at a single point behind turbine')
        ax.grid(True)
        #fig.suptitle('U-velocity evolution over time at a single point behind turbine')
        #fig.suptitle(f'Final time: {sim.tf}s, Simulation time: {duration:.2f}s, ')