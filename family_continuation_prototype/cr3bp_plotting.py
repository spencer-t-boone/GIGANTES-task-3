# Plotting functions for CR3BP orbits
import matplotlib.pyplot as plt
from cr3bp_propagate import propagate
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from cr3bp_functions import jacobi_constant_family
from matplotlib import colormaps


# Plot orbit
def plot_orbit(X_i, t, mu, fig = None, N_points = 1000, orbit_color = 'blue'):
    
    # Propagate orbit and save state
    t_eval = np.linspace(0, t, N_points)
    y_trajectory = propagate(X_i, mu, t, t_eval = t_eval)
    X_trajectory = y_trajectory.y
    
    if fig is None:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
    else:
        ax = fig.gca()
        
    
    ax.scatter(1-mu, 0, 0, s=10,c='grey', marker='o')
    ax.set_xlabel('$x$ [nd]')
    ax.set_ylabel('$y$ [nd]')
    ax.set_zlabel('$z$ [nd]')
    plt.axis('Equal')
    plt.grid()
    plt.show()
    ax.plot(X_trajectory[0,:], X_trajectory[1,:], X_trajectory[2,:], color = orbit_color)
    

# Plot a family of orbits
def plot_family(states_family, periods_family, mu, fig = None, N_points = 1000, spacing = 1):
    
    num_orbits = int(len(periods_family)/spacing)
    
    jc_vec = jacobi_constant_family(states_family, mu)
    
    colormap = colormaps['plasma'].resampled(num_orbits)
    
    if fig is None:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        
    for i in range(0,len(periods_family),spacing):
        plot_orbit(states_family[i,:], periods_family[i], mu, fig, N_points, orbit_color = colormap(i/len(periods_family)) )
        
        
    
