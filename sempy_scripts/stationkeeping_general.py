# Function to perform preliminary stationkeeping analysis using an apse-targeting stationkeeping strategy

# %% Import modules
import numpy as np

from src.init.primary import Primary
from src.dynamics.cr3bp_environment import CR3BP

from src.propagation.propagator import Propagator


import matplotlib.pyplot as plt
import plotly.graph_objects as go

import copy
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import minimize, NonlinearConstraint, Bounds, approx_fprime

from astro_functions import X2oe, oe2X


# Select orbit type to perform analysis
orbit_type_options = ['L2 NRHO', 'L1 NRHO', 'butterfly', 'period-3 halo']
orbit_type = orbit_type_options[2]

# Import selected orbit and set up scenario parameters
# L2 NRHO
if orbit_type == 'L2 NRHO':
    X_i = np.array([ 1.0025797154342282, 0,  0.004882278399970316, 0, -0.005445984590084007,  0])
    tau = 2.287781993612279
    
    num_ref_points = 2 
    noise_freq = 2
    N_rev = 100
    
        

# L1 NRHO
elif orbit_type == 'L1 NRHO':
    X_i = np.array([ 0.9974105021225586,  0,  0.004875390378705563, 0, 0.00546419495753333,  0]) # apoapse
    tau = 2.282108846114326
    
    num_ref_points = 2 
    noise_freq = 2
    N_rev = 100
    
        
   
# Butterfly orbit 
elif orbit_type == 'butterfly':
    X_i = np.array([0.9979394855172609, 0, 0.003723544484837956, 0, -0.0008522149920766561, 0]) # Apoapse
    tau = 3.593492349524928
    
    num_ref_points = 4 
    noise_freq = 2
    N_rev = 50
    


# Period-3 halo   
elif orbit_type == 'period-3 halo':
    X_i = np.array([0.9973675664626385, 0, 0.004768373628104226, 0, 0.00570366409390525, 0]) # Apoapse
    tau = 6.856114820633476
    
    num_ref_points = 6 
    noise_freq = 2
    N_rev = 10


# Set up SEMPy environment
Saturn = Primary.SATURN_NO_MOONS
Enceladus = Primary.ENCELADUS
env = CR3BP(Saturn, Enceladus)

# Create Propagator
prop = Propagator(method = 'DOP853')



# Event function to be triggered when r_dot (derivative of distance from s/c to
# Enceladus) = 0. Computed analytically for CR3BP equations
def event_apse(t, X):
    x = X[0]
    y = X[1]
    z = X[2]
    xdot = X[3]
    ydot = X[4]
    zdot = X[5]
    
    r = np.sqrt( (x - 1 + env.mu)**2 + y**2 + z**2)
    rdot = 2*((x - 1 + env.mu)*xdot + y*ydot + z*zdot)/(2*r)
    
    return rdot


# Event function to detect crossing at 600 km from surface of Enceladus
def event_alt_crossing(t, X):
    x = X[0]
    y = X[1]
    z = X[2]
    xdot = X[3]
    ydot = X[4]
    zdot = X[5]
    
    r = np.sqrt( (x - 1 + env.mu)**2 + y**2 + z**2)
    
    alt_crossing = 600
    r_crossing = (alt_crossing + 252.1)/env.adim_l
    
    return r - r_crossing


# Function to compute maneuver delta-V
def compute_delta_v(u, X, args):
    
    return np.linalg.norm(u)

    
# Function to compute distance from reference at apse point
def compute_distance_from_ref_at_apse(u, X_apo_in, args):
    
    # Extract reference target
    X_target = args['X_ref']

    target_range = args['target_range']
    
    # Apply maneuver
    X_apo = copy.deepcopy(X_apo_in)
    X_apo[3:] += u
    
    prop_to_rp = prop(env, np.linspace(0,2.5*tau,100), X_apo, events=[event_apse])
    
    if prop_to_rp.t_events[0].size >= 1:
        X_peri = prop_to_rp.state_events[0][0]
        
        pos_diff = np.linalg.norm(X_target[:3] - X_peri[:3])
        #print(pos_diff-target_range)
        return -(pos_diff - target_range)
    
    else:
        return 1
    
    


event_alt_crossing.direction = 0
event_alt_crossing.terminal = False


t_vec = np.linspace(0, 2*tau, 5000)

X_apo = copy.deepcopy(X_i)


nav_error_pos_low = (0.250/3)/env.adim_l
nav_error_vel_low = ((0.01/3)/1e3)/env.adim_v
nav_error_pos_high = (0.750/3)/env.adim_l
nav_error_vel_high = ((0.05/3)/1e3)/env.adim_v

sigma_nav_error_low = np.array([nav_error_pos_low, nav_error_pos_low, nav_error_pos_low,
                                nav_error_vel_low, nav_error_vel_low, nav_error_vel_low])

sigma_nav_error_high = np.array([nav_error_pos_high, nav_error_pos_high, nav_error_pos_high,
                                nav_error_vel_high, nav_error_vel_high, nav_error_vel_high])

error_model = 'nominal'
if error_model == 'high':
    sigma_nav_error = sigma_nav_error_high
    sigma_man_error = 0.02
elif error_model == 'nominal':
    sigma_nav_error = sigma_nav_error_low
    sigma_man_error = 0.007
elif error_model == 'none':
    sigma_nav_error = np.zeros(6)
    sigma_man_error = 0.00

np.random.seed(3)


args = {'target_range' : 1.0 / env.adim_l}


X_traj_full = None
u_save = None
u_save_noisy = None


# Propagate initial state to 600 km crossing
X_crossing_next = prop(env, t_vec, X_i, events=[event_alt_crossing])

X_crossing_in = copy.deepcopy(X_crossing_next.state_events[0][0])
X_crossing_in_pert = copy.deepcopy(X_crossing_in)

# Compute reference 
prop_ref_points = prop(env, np.linspace(0, tau, 2), X_i, events=[event_apse])
ref_states = prop_ref_points.state_events[0][1:,:]

# Compute number of maneuvers
N_man = num_ref_points*N_rev


# Loop over number of maneuvers
for N in range(N_man):
     
    # Add navigation error at desired frequency
    if np.mod(N, noise_freq) == 0:
        
        # Add random error
        for i in range(6):
            
            # Add random error in each component
            X_crossing_in_pert[i] += np.random.normal(0, sigma_nav_error[i])
        
        
    # Optimize maneuver to target reference at apse
    ref_ind = np.mod(N, num_ref_points)
    args.update({'X_ref': ref_states[ref_ind, :]})
    args_scipy=(X_crossing_in_pert, args)
    cons = ({'type': 'ineq',
           'fun': compute_distance_from_ref_at_apse,
           'args': args_scipy
           })
    
    bnds = ((-0.5, 0.5),(-0.5, 0.5),(-0.5, 0.5))
    
    u_0 = np.zeros(3)
    
    # Run optimization
    u_opt = minimize(compute_delta_v, u_0, args_scipy, method='SLSQP', constraints = cons,
                     tol=1e-8, options = {'maxiter': 100}, bounds=bnds)
    
    
    
    u_rev = copy.deepcopy(u_opt.x )
    
    # Add maneuver execution errors
    u_rev_mag = np.linalg.norm(u_rev)
    u_rev_noise_std_mag = u_rev_mag*sigma_man_error
    u_rev_noise_std = np.array([u_rev_noise_std_mag]*3)
    
    du_rev_noise = np.zeros(3)
    
    for i in range(3):
        du_rev_noise[i] = np.random.normal(0,u_rev_noise_std[i])
        
    u_rev_noisy = u_rev + du_rev_noise
    
    print('maneuver ', N, ' computed with Delta v: ', np.linalg.norm(u_rev)*env.adim_v, ' km/s')
    
        
    # Apply maneuver and propagate to next 600 km crossing
    X_crossing_in_pert[3:] += u_rev_noisy
    t_vec = np.linspace(0, 2*tau, 5000)
    X_crossing_out_next = prop(env, t_vec, X_crossing_in_pert, events=[event_alt_crossing])
    
    # Extract state at next crossing, make sure not to stop propagation at current crossing
    for j in range(X_crossing_out_next.t_events[0].size):

        if X_crossing_out_next.t_events[0][j] > tau/30:
            X_crossing_out = copy.deepcopy(X_crossing_out_next.state_events[0][j])
            break
       
    # Integrate and save full trajectory
    X_crossing_out_for_save = prop(env, np.linspace(0,X_crossing_out_next.t_events[0][j],500),
                                   X_crossing_in_pert)
    
    # Save trajectory
    if X_traj_full is not None:
        X_traj_full = np.vstack([X_traj_full, X_crossing_out_for_save.state_vec])
        u_save = np.vstack([u_save, u_rev])
        u_save_noisy = np.vstack([u_save_noisy, u_rev_noisy])
    else:
        X_traj_full = copy.deepcopy(X_crossing_out_for_save.state_vec)
        u_save = copy.deepcopy(u_rev)
        u_save_noisy = copy.deepcopy(u_rev_noisy)
    
    
    # Update initial state for next stationkeeping maneuver
    X_crossing_in_pert = copy.deepcopy(X_crossing_out)
        
    
        
        
# Plot resulting trajectory        
fig2 = plt.figure()
ax3 = fig2.add_subplot(projection='3d')
ax3.set_xlabel('x [km]')
ax3.set_ylabel('y [km]')
ax3.set_zlabel('z [km]')

ax3.plot((X_traj_full[:,0] - (1-env.mu))*env.adim_l, X_traj_full[:,1]*env.adim_l, X_traj_full[:,2]*env.adim_l,
          label = 'Trajectory', linewidth = 1.25)        
        
        




# Conservative navigation error estimates
# position: 750 m 3 sigma, velocity 0.05 m/s 3 sigma
# Nominal navigation error estimates
# position: 250 m 3 sigma, velocity: 0.01 m/s 3 sigma
# Source: Enceladus Orbilander, page C-4

# Maneuver execution error estimates
# Model 1 RCS: 0.7% magnitude, 9 mrad direction (cutoff 0.3 m/s)
# Model 1 MEA: 0.02% magnitude, 0.6 mrad direction (cutoff > 0.3 m/s)
# Model 2 RCS: 2.0% magnitude, 12 mrad direction (cutoff 0.4 m/s)
# Model 2 MEA: 0.2% magnitude, 3.5 mrad direction (cutoff > 0.4 m/s)