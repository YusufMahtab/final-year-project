# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc

'Create index values for time and space to iterate over'
N_z = 500 # Number of nodes spaced along the z-axis
N_t = 1000 # Number of iterations the algorithm will perform

times = np.arange(N_t)
z = np.arange(N_z)

'Predetermined constants'
epsilon = np.ones(N_z)*sc.epsilon_0 # Permittivity of free space
mu = np.ones(N_z)*sc.mu_0 # Permeability of free space
c = sc.speed_of_light # Speed of light
sc.physical_constants

# Set a relative permittivity value within the middle 20% of the z-axis
bound1 = int(0.4*N_z)
bound2 = N_z
epsilon[bound1:bound2] *= 4

# Collect above physical constants into two multiplying factors
update_constant = 1/(c*epsilon)
update_constant2 = 1/(c*mu)

'Reset values for all fields to base'
Ex = np.zeros(N_z)
Ey = np.zeros(N_z)
Hx = np.zeros(N_z)
Hy = np.zeros(N_z)

'Plotting the fields'
fig, (ax1, ax2) = plt.subplots(nrows = 2, ncols = 1)
fig.tight_layout()
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized() # Maximise the figure window

# Set axis labels for all graphs
ax1.set_xlabel("z")
ax2.set_xlabel("z")
ax1.set_ylabel("$E_x$ mode")
ax2.set_ylabel("$E_y$ mode")

# Plot all field values
line1, = ax1.plot(z, Ex, 'r-', label = '$E_x$')
line2, = ax1.plot(z, Hy, '#008080', label = '$H_y$')
line3, = ax2.plot(z, Ey, 'r-', label = '$E_y$')
line4, = ax2.plot(z, Hx, '#008080', label = '$H_x$')

# Separate legends for each ax
ax1.legend(loc=4)
ax2.legend(loc=4)

# Colour the dielectric zone
ax1.axvspan(bound1, bound2, facecolor='b', alpha=0.2)
ax2.axvspan(bound1, bound2, facecolor='b', alpha=0.2)
################################################################################
def source_function(t):
    'Set a source function for each mode'
    Ex_source = np.exp(-((t - 6*tau)/tau)**2)*np.cos(t/tau-6*tau)
    Ey_source = np.exp(-((t - 6*tau)/tau)**2)*np.cos(t/tau-6*tau)
    return Ex_source, Ey_source
################################################################################
def E_update():
    'Imposing absorbing boundary conditions for E'
    Ex[0] = Ex[1]
    Ey[0] = Ey[1]

    'Update rest of E as the algorithm specifies'
    Ex[1:] -= (Hy[1:] - Hy[0:-1])*update_constant[1:]
    Ey[1:] += (Hx[1:] - Hx[0:-1])*update_constant[1:]

    'Additive source for E'
    Ex[source_node] += source_function(t)[0]
    Ey[source_node] += source_function(t)[1]
################################################################################
def H_update():
    'Imposing absorbing boundary conditions for H'
    Hy[-1] = Hy[-2]
    Hx[-1] = Hx[-2]

    'Creating a TF/SF boundary at the source'
    Hy[source_node - 1] += source_function(t)[0]*update_constant2[source_node-1]
    Hx[source_node - 1] -= source_function(t)[1]*update_constant2[source_node-1]

    'Update all but the last H node'
    Hy[:-1] -= (Ex[1:] - Ex[:-1])*update_constant2[:-1]
    Hx[:-1] += (Ey[1:] - Ey[:-1])*update_constant2[:-1]
################################################################################

for t in times:
    source_node = int(0.2*N_z)
    tau = 20
    E_update()
    H_update()

    # Drawing a source boundary line
    ax1.axvline(source_node, linewidth=0.5, color='k', linestyle='--')
    ax2.axvline(source_node, linewidth=0.5, color='k', linestyle='--')

    # Update the title with the current time-step
    ax1.set_title(label=("Step %.2f of " %t + "%.2f" %(N_t-1)))

    # Update all lines with refreshed values
    line1.set_ydata(Ex)
    line2.set_ydata(377*Hy) # Scaled up to compare with Ex
    line3.set_ydata(Ey)
    line4.set_ydata(377*Hx) # Scaled up to compare with Ey

    ax1.set_ylim(-2, 2)
    ax2.set_ylim(-2, 2)

    # Refresh the entire figure
    fig.canvas.draw()
    fig.canvas.flush_events()
    plt.show()
