import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import cm
import numpy as np
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
        #fig.colorbar(cp, ax=ax).set_label('Divergence')
        cbar = fig.colorbar(cp,orientation='horizontal')
        ax.set_title('2D Divergence Field (Interior Points)')

    def plot_pressure(pressure, fig, ax):
        """
        Plots the 2D pressure field.
        """
        cp = ax.imshow(pressure[1:-1, 1:-1], cmap='plasma', origin='lower',
                    aspect='equal')
        #fig.colorbar(cp, ax=ax).set_label('Pressure')
        cbar = fig.colorbar(cp,orientation='horizontal')
        ax.set_title('2D Pressure Field')

    def plot_u(u, fig, ax):
        """
        Plots the 2D x-velocity field with a color bar that matches the plot height and includes tick labels.
        """
        # Plot the velocity field with specified vmin and vmax
        cp = ax.imshow(u[1:-1, 1:-1], cmap='plasma', origin='lower', aspect='equal')

        # Manually add the colorbar as a separate axis at the bottom
        #cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.03])  # [left, bottom, width, height] for positioning
        #cbar = fig.colorbar(cp, cax=cbar_ax, orientation='horizontal')
        cbar = fig.colorbar(cp,orientation='horizontal')
        #cbar.set_label('U-velocity')

        # Set plot title
        ax.set_title('2D U-velocity Field', fontsize=14)

    def plot_u_OG(u, fig, ax):
        """
        Plots the 2D x-velocity field.
        """
        # Plot the velocity field with specified vmin and vmax
        cp = ax.imshow(u[1:-1, 1:-1], cmap='plasma', origin='lower',
                    aspect='equal')
        fig.colorbar(cp, ax=ax).set_label('U-velocity')
        ax.set_title('2D U-velocity Field')

    def plot_p_centerline(centerline_data, fig, ax):
        """
        Plots the centerline pressure profile.
        
        Parameters:
        - centerline_data: array. 
        - fig: Figure object.
        - ax: Axes object to plot on.
        """
        # Loop through each time-centerline pair and plot
        ax.plot(centerline_data)
        ax.set_xlabel('Horizontal Position (grid index)')
        ax.set_ylabel('Pressure at Centerline')
        ax.set_title('Centerline Pressure Profile')
        ax.grid(True)

    def plot_u_centerline(centerline_data, fig, ax):
        """
        Plots the centerline velocity profile.
        
        Parameters:
        - centerline_data: array.
        - fig: Figure object.
        - ax: Axes object to plot on.
        """
        # Loop through each time-centerline pair and plot
        ax.plot(centerline_data)
        ax.set_xlabel('Horizontal Position (grid index)')
        ax.set_ylabel('Velocity at Centerline')
        ax.set_title('Centerline Velocity Profile')
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

    def plot_u_with_streamlines(u, v, fig, ax, dx, dy):
        """
        Plots the 2D x-velocity field with streamlines.
        Need to be debugged. 
        Parameters:
            u (2D array): The x-velocity component of the field.
            v (2D array): The y-velocity component of the field.
            fig (matplotlib.figure.Figure): The figure object for the plot.
            ax (matplotlib.axes._subplots.AxesSubplot): The axis object for the plot.
            dx (float): Grid spacing in the x-direction.
            dy (float): Grid spacing in the y-direction.
        """
        # Create mesh grid for streamline plotting
        ny, nx = u.shape
        x = np.linspace(0, (nx - 2) * dx, nx - 2)
        y = np.linspace(0, (ny - 2) * dy, ny - 2)
        X, Y = np.meshgrid(x, y)
        
        # Plot the velocity magnitude as a background color map
        velocity_magnitude = np.sqrt(u[1:-1, 1:-1]**2 + v[1:-1, 1:-1]**2)
        cp = ax.imshow(velocity_magnitude, cmap='plasma', origin='lower',
                    extent=(0, (nx - 2) * dx, 0, (ny - 2) * dy), aspect='equal')
        fig.colorbar(cp, ax=ax).set_label('Velocity Magnitude')
        
        # Plot streamlines on top of the velocity field
        ax.streamplot(X, Y, u[1:-1, 1:-1], v[1:-1, 1:-1], color='white', density=1.5)
        
        # Set plot title and labels
        ax.set_title('2D U-velocity Field with Streamlines')
        ax.set_xlabel('X Position')
        ax.set_ylabel('Y Position')