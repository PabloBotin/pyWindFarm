from Fluid import Fluid
from Visualize import Visualize
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import time
import os
import shutil
from matplotlib import rcParams
import zipfile

# ------------------------------------------------------------------------------
# Preamble - copies this script to unique directory for each time it is run.
# ------------------------------------------------------------------------------
today = time.strftime('%Y_%m_%d')
run_id = f"run{str(int(time.time()))[-5:]}" # 5 digits of current time in milliseconds
root = os.getcwd()
mdir = f"{root}/figures/{today}/{run_id}"
os.makedirs(mdir, exist_ok=True)
this_file = os.path.basename(__file__)
shutil.copy(this_file, f"{mdir}/{os.path.basename(this_file)}") # main.py 
script_dir = os.path.dirname(os.path.abspath(__file__))
Fluid_source = os.path.join(script_dir, "Fluid.py")
Fluid_destination = f"{mdir}/Fluid.py"
shutil.copy(Fluid_source, Fluid_destination) # Fluid.py
Visualize_source = os.path.join(script_dir, "Visualize.py")
Visualize_destination = f"{mdir}/Visualize.py"
shutil.copy(Visualize_source, Visualize_destination) # Visualize.py
Visualize_source = os.path.join(script_dir, "Config.py")
Visualize_destination = f"{mdir}/Config.py"
shutil.copy(Visualize_source, Visualize_destination) # Visualize.py

rcParams.update({
    ## Constrained Layout and tiny margins are necessary
    'figure.constrained_layout.use': True,
    'axes.xmargin': 0.005,
    'axes.ymargin': 0.005,
    ## These are merely stylistic, feel free to adjust
    'figure.subplot.left': 0.05,
    'figure.subplot.right': 0.95,
    'figure.subplot.bottom': 0.05,
    'figure.subplot.top': 0.95,
    'figure.subplot.wspace': 0.1,
    'figure.subplot.hspace': 0.1,
    'figure.constrained_layout.h_pad': 0.05,
    'figure.constrained_layout.w_pad': 0.05,
    'figure.constrained_layout.hspace': 0.05,
    'figure.constrained_layout.wspace': 0.05
    })

def main():

    # Start the timer
    start_time = time.time()
    sim = Fluid() # Create instance of the class.

    # Run solver
    # u, v, p, centerline_p, centerline_u = sim.solve(sim.u, sim.v, sim.p)
    u, v, p = sim.solve(sim.u, sim.v, sim.p)

    # End the timer
    end_time = time.time()
    # Calculate the duration and print it
    duration = end_time - start_time
    print(f"Simulation took {duration:.2f} seconds to run.")

    # Plot. 
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300, constrained_layout=True)  # Disable layout engine
    Visualize.plot_u(sim.u, fig, ax)
    ax.set_title(f'U-velocity,tf={sim.tf},t_run={duration:.2f}s,{sim.ny}x{sim.nx},tol={sim.tol},Re={sim.Re},maxiter={sim.maxiter}',fontsize=12, fontweight='bold')
    fig.savefig(f'{mdir}/u.png')
    plt.close(fig)  # Close the figure to avoid displaying it

    # # Plot centerline u 
    # fig, ax = plt.subplots(figsize=(7, 6), dpi=300)  # Create a single plot
    # Visualize.plot_u_centerline(centerline_u, fig, ax)
    # fig.suptitle(f'tf:{sim.tf}s, t_run:{duration:.2f}s,Re:{sim.Re}')
    # fig.savefig(f'{mdir}/centerline_u.png')
    # plt.close(fig)  # Close the figure to avoid displaying it

    # # Plot centerline.
    # fig, ax = plt.subplots(figsize=(7, 6), dpi=300)  # Create a single plot
    # Visualize.plot_p_centerline(centerline_p, fig, ax)
    # fig.suptitle(f'tf:{sim.tf}s, t_run:{duration:.2f}s,Re:{sim.Re}')
    # fig.savefig(f'{mdir}/centerline_p.png')
    # plt.close(fig)  # Close the figure to avoid displaying it



if __name__ == "__main__": # Only run the functions defined in this code. 
    main()
