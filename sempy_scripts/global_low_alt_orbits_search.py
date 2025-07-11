# This script will perform a global search over two-body Enceladus-centric
# orbital parameters in the Saturn-Enceladus CR3BP
# This corresponds to Task 1.B.



# %% Import modules
import numpy as np

from src.init.primary import Primary
from src.dynamics.cr3bp_environment import CR3BP

from src.propagation.propagator import Propagator

import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.colors as co


import matplotlib.pyplot as plt
import copy
from mpl_toolkits.mplot3d import Axes3D

# Set up SEMPy 
Saturn = Primary.SATURN_NO_MOONS
Enceladus = Primary.ENCELADUS
env = CR3BP(Saturn, Enceladus)

# -----------------------------------------------------------------------------
# INPUT PARAMETERS
# -----------------------------------------------------------------------------
# Set up simulation parameters
# Max integration time (in non-dim)
tau_days = 150 # days
tau = tau_days/env.adim_t*3600*24


# Ranges of orbital elements
# Vary radius of periapse, radius of apoapse, inclination, or RAAN
r_p_range = np.linspace(20+252.1,100+252.1,9)
r_a_range = np.linspace(50+252.1,200+252.1,16)
i_range = 65.0 # Currently script runs for fixed inclination value
raan_range = np.linspace(0,360,9)
arg_p_range = np.linspace(0,360,16)

# Note: the script currently will varry radius of periapse, radius of apoapse, RAAN
# and argument of periapse. One could also set up the script to vary inclination as well.


# -----------------------------------------------------------------------------
# MAIN SCRIPT
# -----------------------------------------------------------------------------
# Create Propagator
prop = Propagator(method='DOP853')



fig = go.FigureWidget()

def ms(x, y, z, radius, resolution=20):
    """Return the coordinates for plotting a sphere centered at (x,y,z)"""
    u, v = np.mgrid[0:2*np.pi:resolution*2j, 0:np.pi:resolution*1j]
    X = radius * np.cos(u)*np.sin(v) + x
    Y = radius * np.sin(u)*np.sin(v) + y
    Z = radius * np.cos(v) + z
    return (X, Y, Z)


# Function to compute distance from secondary in CR3BP units
def compute_distance_from_moon(x, mu_cr3bp):
    num_states = x.shape[0]
    moon_pos = np.array([1-mu_cr3bp, 0, 0])
    dist = np.zeros(num_states)
    for i in range(num_states):
        dist[i] = np.linalg.norm(x[i,:3] - moon_pos)
        
    return dist


# Convert state (Cartesian) to orbital elements
# Note: mu refers to the gravitational parameter of the central body, not 
# the non-dimensional CR3BP mu parameter
def X2oe(X, mu):
    r_vec = X[:3]
    r = np.linalg.norm(r_vec)
    v_vec = X[3:6]
    v = np.linalg.norm(v_vec)
    h_vec = np.cross(r_vec, v_vec)
    h = np.linalg.norm(h_vec)
    n_vec = np.cross(np.array([0, 0, 1]), h_vec)
    n = np.linalg.norm(n_vec)
    e_vec = ((v**2 - mu/r)*r_vec - np.dot(r_vec, v_vec)*v_vec)/mu
    e = np.linalg.norm(e_vec)
    
    spec_en = v**2/2 - mu/r
    a = -mu/(2*spec_en)
    p = a*(1 - e**2)
    
    i = np.arccos(h_vec[2]/h)
    RAAN = np.arccos(n_vec[0]/n)
    if n_vec[1] < 0:
        RAAN = 2*np.pi - RAAN
        
    omega = np.arccos(np.dot(n_vec,e_vec)/(n*e))
    if e_vec[2] < 0:
        omega = 2*np.pi - omega
        
    nu = np.arccos(np.dot(e_vec,r_vec)/(e*r))
    if np.dot(r_vec, v_vec) < 0:
        nu = 2*np.pi - nu
        
    return np.array([a, e, i, RAAN, omega, nu])


# Convert orbital elements to Cartesian state vector
# Note: mu refers to the gravitational parameter of the central body, not 
# the non-dimensional CR3BP mu parameter
def oe2X(oe, mu):
    a = oe[0]
    e = oe[1]
    i = oe[2]
    RAAN = oe[3]
    omega = oe[4]
    nu = oe[5]
    
    p = a*(1-e**2)
    
    r_pqw = np.array([p*np.cos(nu)/(1 + e*np.cos(nu)),
                      p*np.sin(nu)/(1 + e*np.cos(nu)),
                      0.0])
    
    v_pqw = np.array([-np.sqrt(mu/p)*np.sin(nu),
                      np.sqrt(mu/p)*(e + np.cos(nu)),
                      0.0])
    
    rot_ijk_pqw = np.array([[np.cos(RAAN)*np.cos(omega) - np.sin(RAAN)*np.sin(omega)*np.cos(i),    -np.cos(RAAN)*np.sin(omega) - np.sin(RAAN)*np.cos(omega)*np.cos(i),    np.sin(RAAN)*np.sin(i)],
                            [np.sin(RAAN)*np.cos(omega) + np.cos(RAAN)*np.sin(omega)*np.cos(i),    -np.sin(RAAN)*np.sin(omega) + np.cos(RAAN)*np.cos(omega)*np.cos(i),    -np.cos(RAAN)*np.sin(i)],
                            [np.sin(omega)*np.sin(i),    np.cos(omega)*np.sin(i),    np.cos(i)]])
    
    
    r_ijk = rot_ijk_pqw@r_pqw
    v_ijk = rot_ijk_pqw@v_pqw
    
    return np.hstack([r_ijk, v_ijk])


# Convert non-dimensional CR3BP state vector to state vector in the secondary-body
# inertial frame (in km)
def nd2km(X, env):
    X_km = np.zeros(6)
    X_km[:3] = copy.deepcopy(X[:3])
    X_km[0] += -(1-env.mu)
    X_km[:3] = X_km[:3] * env.adim_l
    X_km[3:6] = X[3:6] * env.adim_v
    
    return X_km

# Convert state vector in the secondary-body inertial frame (in km) to 
# non-dimensional CR3BP state vector
def km2nd(X, env):
    X_nd = np.zeros(6)
    X_nd[:3] = X[:3] / env.adim_l
    X_nd[3:6] = X[3:6] / env.adim_v
    
    X_nd[0] += 1 - env.mu
    
    return X_nd




# Event function to trigger in case of impact with surface of Enceladus
def event_impact(t,X):
    x = X[0]
    y = X[1]
    z = X[2]
    
    r = np.sqrt( (x - 1 + env.mu)**2 + y**2 + z**2)
    return (r - Enceladus.Rm/env.adim_l)


# Event function to trigger in case of escape from vicinity of Enceladus
def event_escape(t,X):
    r_escape = 800 # km
    x = X[0]
    y = X[1]
    z = X[2]
    
    r = np.sqrt( (x - 1 + env.mu)**2 + y**2 + z**2)
    return (r - r_escape/env.adim_l)

event_impact.terminal = True
event_escape.terminal = True


# Time vector to output states at various sizes
t_int = 0
t_vec = np.linspace(t_int, tau, 5000)


# Initialize arrays to save values
t_end_save = np.zeros((r_p_range.size, r_a_range.size, raan_range.size, arg_p_range.size))
init_oe_save = np.zeros((r_p_range.size, r_a_range.size, raan_range.size, arg_p_range.size, 6)) 
init_cond_save = np.zeros((r_p_range.size, r_a_range.size, raan_range.size, arg_p_range.size, 6)) 
prop_end_flag = np.zeros((r_p_range.size, r_a_range.size, raan_range.size, arg_p_range.size))


# Loop over initial conditions to perform global search
# Loop over radius of periapse
for i in range(r_p_range.size):
    r_p = r_p_range[i]
    print('r_p: ', r_p, ' km')
    
    # Loop over radius of apoapse
    for j in range(r_a_range.size):
        r_a = r_a_range[j]
        
        # Loop over RAAN
        for k in range(raan_range.size):
            init_raan = raan_range[k]
            
            for m in range(arg_p_range.size):
                init_arg_p = arg_p_range[m]
            
                # Set up initial orbital elements vector (with argument of periapse at 90 degrees)
                init_i = i_range
                init_a = (r_p + r_a)/2 
                init_e = 1 - np.min([r_p,r_a])/init_a
                
                init_oe = np.array([init_a, init_e, init_i/180*np.pi, init_raan/180*np.pi, init_arg_p/180*np.pi, np.pi])
                
                # Convert git pushorbital elements to non-dimensional CR3BP state vector
                init_cond_km = oe2X(init_oe, Enceladus.GM)
                init_cond_new = km2nd(init_cond_km, env)
                
                # Propagate initial conditions until either impact or escape occurs
                prop_int_iter = prop(env, t_vec, init_cond_new, events=[event_impact,event_escape])
                
                # Flag if impact, escape, or time limit ended propagation
                if prop_int_iter.t_events[0].size > 0:
                    prop_end_flag[i,j,k,m] = 1 # impact
                    t_end = prop_int_iter.t_events[0]
                    
                elif prop_int_iter.t_events[1].size > 0:
                    prop_end_flag[i,j,k,m] = 2 # escape
                    t_end = prop_int_iter.t_events[1]
                else:
                    prop_end_flag[i,j,k,m] = 3 # remained in orbit around Enceladus
                    t_end = tau
                    
                t_end_save[i,j,k,m] = copy.deepcopy(t_end)
                    
                
                # Store initial conditions and orbital elements
                states_save = copy.deepcopy(prop_int_iter.state_vec)          
                init_oe_save[i,j,k,m,:] = copy.deepcopy(init_oe)           
                init_cond_save[i,j,k,m,:] = copy.deepcopy(init_cond_new)
                
                
                # Plot orbits and characteristics for select interesting orbits
                if t_end > tau+1:
                    fig2 = plt.figure()
                    ax3 = fig2.add_subplot(projection='3d')
                    ax3.set_xlabel('x [km]')
                    ax3.set_ylabel('y [km]')
                    ax3.set_zlabel('z [km]')
                    
                    # Plot orbits 
                    ax3.plot((states_save[:,0]-1+env.mu)*env.adim_l,states_save[:,1]*env.adim_l,states_save[:,2]*env.adim_l,
                              linewidth = 1.25, color = 'blue')
                    
                    states_save_oe_vec = np.zeros((states_save.shape[0],6))
                    for s in range(states_save.shape[0]):
                        states_save_km = nd2km(states_save[s,:], env)
                        states_save_oe_vec[s,:] = X2oe(states_save_km, Enceladus.GM)
                        
                    # Plot argument of periapse vs eccentricity (frozen orbit parameters)
                    fig3 = plt.figure()
                    plt.scatter(states_save_oe_vec[:,4]*180/np.pi,states_save_oe_vec[:,1])
                    plt.xlabel('w [deg]')
                    plt.ylabel('e')
                                                                                                                                               
                        
                        
                    
    
# Identify orbital parameters corresponding to maximum lifetime orbit for 
# each r_p / r_a pair                   
max_ind = np.unravel_index(np.argmax(t_end_save, axis=None), t_end_save.shape)
t_end_save_max = np.max(t_end_save,axis=(2,3))


# Plot contour graph of r_p vs. r_a  vs. max lifetime orbit
plt.figure()
plt.contourf(r_p_range-252.1, r_a_range-252.1,t_end_save_max.T * env.adim_t/3600/24, 100)
plt.colorbar(label='Days')
plt.xlabel('$alt_p$ [km]')
plt.ylabel('$alt_a$ [km]')