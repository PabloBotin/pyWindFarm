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
# this_file = os.path.basename(__file__)
# shutil.copy(this_file, f"{mdir}/{os.path.basename(this_file)}") # main.py 
# script_dir = os.path.dirname(os.path.abspath(__file__))
# Fluid_source = os.path.join(script_dir, "Fluid.py")
# Fluid_destination = f"{mdir}/Fluid.py"
# shutil.copy(Fluid_source, Fluid_destination) # Fluid.py
# Visualize_source = os.path.join(script_dir, "Visualize.py")
# Visualize_destination = f"{mdir}/Visualize.py"
# shutil.copy(Visualize_source, Visualize_destination) # Visualize.py


# # Create folder to save figures and scripts. 
# today = time.strftime('%Y_%m_%d') # Get date.
# run_id = f"run{str(int(time.time()))[-5:]}" # 5 digits of current time in milliseconds
# root = os.getcwd() # Get current folder path. 
# mdir = f"{root}/figures/{today}/{run_id}" # create folder
# scripts = os.path.join(root, "python_files_backup.zip") # Create zip folder.
# with zipfile.ZipFile(scripts, 'w') as zipf:
#     for root, _, files in os.walk(root):
#         for file in files:
#             if file.endswith(".py"):  # Only add Python files
#                 full_path = os.path.join(root, file)
#                 # Ensure unique paths within the zip file by preserving directory structure
#                 archive_name = os.path.relpath(full_path, root)
#                 zipf.write(full_path, archive_name)
# print(f"Created a folder for figures and compressed all Python files into a zip!")


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
    u, v, p, u_point_story = sim.solve(sim.u, sim.v, sim.p)

    # End the timer
    end_time = time.time()
    # Calculate the duration and print it
    duration = end_time - start_time
    print(f"Simulation took {duration:.2f} seconds to run.")
    
    fig, ax = plt.subplots(1, 3, figsize=(14, 4),
                        sharey=True, sharex=True, dpi=300)
    Visualize.plot_divergence(sim.u, sim.v, sim.dx, sim.dy, fig, ax[0])
    Visualize.plot_pressure(sim.p, fig, ax[1])
    Visualize.plot_u(sim.u, fig, ax[2])

    #fig.suptitle(f'cells={sim.ny}x{sim.nx}, tol={sim.tol}, maxiter={sim.maxiter}, final time: {sim.tf}s, simulation time:{duration:.2f}')
    fig.suptitle(f'tf={sim.tf},t_run={duration:.2f}s,{sim.ny}x{sim.nx},tol={sim.tol},nu={sim.nu},F={sim.F},u_in={sim.u_in},maxiter={sim.maxiter}')
    fig.savefig(f'{mdir}/triple_panel.png')
    plt.close(fig)  # Close the figure to avoid displaying it

    # Plot centerline u progression. 
    fig, ax = plt.subplots(figsize=(7, 6), dpi=300)  # Create a single plot
    #Visualize.plot_u_centerlines(centerline_data, fig, ax)
    Visualize.plot_u_point_story(u_point_story, sim.tf, fig, ax)
    fig.suptitle(f'tf={sim.tf},t_run={duration:.2f}s,{sim.ny}x{sim.nx},tol={sim.tol},nu={sim.nu},F={sim.F},u_in={sim.u_in},maxiter={sim.maxiter}')
    fig.savefig(f'{mdir}/u_point_story.png')
    plt.close(fig)  # Close the figure to avoid displaying it

if __name__ == "__main__": # Only run the functions defined in this code. 
    main()
