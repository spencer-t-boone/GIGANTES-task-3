# Plotting functions for CR3BP orbits
import matplotlib.pyplot as plt
from cr3bp_propagate import propagate
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from cr3bp_functions import jacobi_constant_family
from matplotlib import colormaps
from copy import deepcopy


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
        
    
    ax.scatter(1-mu, 0, 0, s=10,c='teal', marker='o')
    ax.set_xlabel('$x$ [nd]')
    ax.set_ylabel('$y$ [nd]')
    ax.set_zlabel('$z$ [nd]')
    plt.axis('Equal')
    plt.grid()
    plt.show()
    ax.plot(X_trajectory[0,:], X_trajectory[1,:], X_trajectory[2,:], color = orbit_color)
    
    
# Plot orbit
def plot_orbit_km(X_i, t, mu, R_sec, r_sec, fig = None, N_points = 1000, orbit_color = 'blue'):
    
    # Propagate orbit and save state
    t_eval = np.linspace(0, t, N_points)
    y_trajectory = propagate(X_i, mu, t, t_eval = t_eval)
    X_trajectory = y_trajectory.y
    X_trajectory_km = X_trajectory[:3,:]
    X_trajectory_km[0,:] += (mu - 1)
    X_trajectory_km *= r_sec
    
    if fig is None:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
    else:
        ax = fig.gca()
        
    plot_sphere(0, 0, 0, R_sec, fig, color = 'teal')
    
    ax.scatter(0, 0, 0, s=10,c='teal', marker='o')
    ax.set_xlabel('$x$ [km]')
    ax.set_ylabel('$y$ [km]')
    ax.set_zlabel('$z$ [km]')
    plt.axis('Equal')
    plt.grid()
    plt.show()
    ax.plot(X_trajectory[0,:], X_trajectory[1,:], X_trajectory[2,:], color = orbit_color)
    
    
    

# Plot a family of orbits
def plot_family(states_family, periods_family, mu, fig = None, N_points = 1000, spacing = 1, frame = 'synodic', R_sec = 252.1, r_sec = 238400,
                variable_mu = 0):
    
    if variable_mu == 1:
        mu_vector = deepcopy(mu)
    
    num_orbits = int(len(periods_family)/spacing)
    
    jc_vec = jacobi_constant_family(states_family, mu)
    
    colormap = colormaps['plasma'].resampled(num_orbits)
    
    if fig is None:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        
    for i in range(0,len(periods_family),spacing):
        if variable_mu == 1:
            mu = mu_vector[i]
            
        if frame == 'synodic':
            plot_orbit(states_family[i,:], periods_family[i], mu, fig, N_points, orbit_color = colormap(i/len(periods_family)) )
        elif frame == 'sec-centric':
            plot_orbit_km(states_family[i,:], periods_family[i], mu, R_sec, r_sec, fig, N_points, orbit_color = colormap(i/len(periods_family)))
        
        
# Plot sphere
def plot_sphere(x_cen, y_cen, z_cen, R, fig, color = 'grey'):
    ax = fig.gca()
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)
    ax.plot_surface(x*R + x_cen, y*R + y_cen, z*R + z_cen, color=color, alpha = 0.1)    