# %% Import modules
import numpy as np

from src.init.primary import Primary
from src.dynamics.cr3bp_environment import CR3BP

from src.propagation.propagator import Propagator

import matplotlib.pyplot as plt

from scipy.optimize import fsolve
from copy import deepcopy

import matplotlib.pyplot as plt
import copy
from mpl_toolkits.mplot3d import Axes3D
import labellines

from astro_functions import *

# Select orbit type to perform analysis
orbit_type_options = ['L2 NRHO', 'L1 NRHO', 'butterfly', 'period-3 halo']
orbit_type = orbit_type_options[0]

# Selection location for maneuver 
# Options for nrho and period-3 halo maneuver locations: 'periapse', 'apoapse'
# Options for butterfly maneuver locations: 'apoapse L2', 'apoapose L1'
maneuver_location = 'periapse'


# Set up SEMPy environment
Saturn = Primary.SATURN_NO_MOONS
Enceladus = Primary.ENCELADUS
env = CR3BP(Saturn, Enceladus)

# Create Propagator
prop = Propagator(method = 'DOP853')

# Import selected orbit and set up scenario parameters
# L2 NRHO
if orbit_type == 'L2 NRHO':
    X_peri = np.array([1.000064575406844, 0, -0.001183305369127517, 9.288141002844474e-13, 0.01685176088135196, 1.251253788081552e-11])
    tau = 2.287781993612279
    
    if maneuver_location == 'apoapse':
        X_apo = prop(env, [0,tau/2], X_peri).state_vec[-1] + 0     
        X_i = X_apo
        min_delta_v = 0.08/env.adim_v # min Delta-V 
        max_delta_v = 0.3/env.adim_v  # Max Delta-V
        
    elif maneuver_location == 'periapse':
        X_i = X_peri
        min_delta_v = 0.03/env.adim_v # min Delta-V 
        max_delta_v = 0.3/env.adim_v  # Max Delta-V
        

# L1 NRHO
elif orbit_type == 'L1 NRHO':
    X_peri = np.array([0.9999391020169318, 0, -0.001183305369127517, 0, -0.01686127728398247, 0]) # periapse
    tau = 2.282108846114326
    
    if maneuver_location == 'apoapse':
        X_apo = prop(env, [0,tau/2], X_peri).state_vec[-1] + 0     
        X_i = X_apo
        min_delta_v = 0.08/env.adim_v # min Delta-V 
        max_delta_v = 0.3/env.adim_v  # Max Delta-V
        
    elif maneuver_location == 'periapse':
        X_i = X_peri
        min_delta_v = 0.03/env.adim_v # min Delta-V 
        max_delta_v = 0.3/env.adim_v  # Max Delta-V
        
   
# Butterfly orbit 
elif orbit_type == 'butterfly':
    X_apo_L2 = np.array([0.9979394855172609, 0, 0.003723544484837956, 0, -0.0008522149920766561, 0]) # Apoapse
    tau = 3.593492349524928
    
    if maneuver_location == 'apoapse L2':
        X_i = X_apo_L2
        min_delta_v = 0.05/env.adim_v # min Delta-V 
        max_delta_v = 0.3/env.adim_v  # Max Delta-V
        
    elif maneuver_location == 'apoapse L1':
        X_apo_L1 = prop(env, [0,tau/2], X_apo_L2).state_vec[-1] + 0
        X_i = X_apo_L1
        min_delta_v = 0.05/env.adim_v # min Delta-V 
        max_delta_v = 0.3/env.adim_v  # Max Delta-V


# Period-3 halo   
elif orbit_type == 'period-3 halo':
    X_apo = np.array([0.9973675664626385, 0, 0.004768373628104226, 0, 0.00570366409390525, 0]) # Apoapse
    tau = 6.856114820633476
    
    if maneuver_location == 'apoapse':
        X_i = X_apo
        min_delta_v = 0.08/env.adim_v # min Delta-V 
        max_delta_v = 0.3/env.adim_v  # Max Delta-V
        
    elif maneuver_location == 'periapse':
        X_peri = prop(env, [0,tau/2], X_apo).state_vec[-1] + 0
        X_i = X_peri
        min_delta_v = 0.03/env.adim_v # min Delta-V 
        max_delta_v = 0.3/env.adim_v  # Max Delta-V



# Set up conditions for simulation
num_dvs = 301  # Number of delta-v magnitudes to simulare
delta_v_mags = np.linspace(min_delta_v, max_delta_v, num_dvs)


# Set up Poincare section event
event1 = event_poincare_y_axis 
event1.terminal = False


# Set up event checking impact with Enceladus
# stop propagation if trajectory enters within 20 km of surface of Enceladus (or
# other pre-defined 'danger zone')
event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, env.mu, (Enceladus.Rm + 20)/env.adim_l)
event_impact_enceladus.terminal = True

# Set up time vector (note time is negative because we are propagating backwards)
t_vec = np.linspace(0, -50*tau, 5000) # Make sure time is sufficient to reach negative x-axis crossing


# Set up plot on which to plot manifolds
fig2 = plt.figure()
ax3 = fig2.add_subplot(projection='3d')
ax3.set_xlabel('$x$ [nd]')
ax3.set_ylabel('$y$ [nd]')
ax3.set_zlabel('$z$ [nd]')
ax3.set_aspect('equal')

# Initialize arrays to store data
X_i_dvs = np.zeros((num_dvs,6))
X_saturn_centric = np.zeros((num_dvs, 6))
oe_saturn_centric = np.zeros((num_dvs,6))
ra_saturn_centric_save = np.zeros(num_dvs)
rp_saturn_centric_save = np.zeros(num_dvs)
t_tp = np.zeros(num_dvs)

# Loop over delta-V magnitudes
for i in range(num_dvs):
    
    # Apply Delta V in in-track velocity direction
    X_i_dv = copy.deepcopy(X_i)
    X_i_dv[3:] += delta_v_mags[i]*X_i[3:]/np.linalg.norm(X_i[3:])
    X_i_dvs[i,:] = copy.deepcopy(X_i_dv)
    
    # Integrate state with delta-v applied (backwards)
    X_sol = prop(env, t_vec, X_i_dv, events = [event1, event_impact_enceladus])
    
    
    # If no crossings of y=0 occurred, continue to next iteration
    num_crossings = X_sol.t_events[0].size
    if num_crossings == 0:
        print('No crossing of negative x-axis, continuing to next delta-V magnitude')
        continue
    
    # Loop through event triggers to find first instance of negative x-axis crossing
    crossing_found = 0
    for j in range(num_crossings):
        if X_sol.state_events[0][j,0] < 0:
            crossing_found = 1
            break
        
    # Exit loop if trajectory impacted surface of Enceladus before first crossing
    if X_sol.t_events[1].size > 0:
        if X_sol.t_events[1][0] > X_sol.t_events[0][j]:
            print('Impacted surface of Enceladus, continuing to next delta-V magnitude')
            continue
    
    # Exit loop if no crossings of negative x-axis occurred
    if crossing_found == 0:
        print('No crossing of negative x-axis, continuing to next delta-V magnitude')
        continue
    
    else: # Crossing detected
        # Extract state and time data for first instance of negative x-axis crossing
        t_tp[i] = X_sol.t_events[0][j]
        X_int = X_sol.state_vec[X_sol.t_vec > t_tp[i]]
        X_tp = X_sol.state_events[0][j]
                                     
    
    # Compute radius of periapse and apoapse for Tisserand-Poincare maps
    X_saturn_centric[i,:] = convert_cr3bp_frame_to_saturn_centric(X_tp, env)
    oe_saturn_centric[i,:] = X2oe(X_saturn_centric[i,:], Saturn.GM)
    a = oe_saturn_centric[i,0]
    e = oe_saturn_centric[i,1]
    rp_saturn_centric = a*(1-e)
    ra_saturn_centric = a*(1+e)
            
    ra_saturn_centric_save[i] = ra_saturn_centric
    rp_saturn_centric_save[i] = rp_saturn_centric
    
    
    # Plot select manifolds in 3d
    if np.mod(i,100) == 0:
        ax3.plot(X_int[:,0], X_int[:,1], X_int[:,2])
        print('Iteration ', i)
    
   
    
# Function for creating Tisserand-Poincare graph lines    
def tp_graph_test(v_inf, r_p_range, a = Enceladus.a):
    
    # Compute normalized v_infinity
    vEnc = np.sqrt(Saturn.GM / a)   
    v_inf_norm = v_inf/vEnc
    
    # Loop over provided radii of periapse
    r_a_list = np.zeros(r_p_range.size)
    
    # Compute isoline values following Eqn. 11 from Campagnola and Russell
    # "Endgame Problem Part 2: Multibody Technique and the Tisserand–Poincaré Graph"
    for i in range(r_p_range.size):
        fun_to_solve = lambda r_a: isoline_fun(r_a, r_p_range[i], v_inf_norm, a)
        
        if i == 0:
            a_guess = deepcopy(a)
        else:
            a_guess = deepcopy(r_a_list[i-1])
        r_a_list[i] = fsolve(fun_to_solve, a_guess)
        
    return r_a_list
        
    
# Function for Eqn. 11 from Campagnola and Russell
def isoline_fun(r_a, r_p, v_inf, a):
    
    fun = 2*a/(r_a + r_p) + 2*np.sqrt(2*r_p*r_a/((r_a + r_p) * a)) - (3 - v_inf**2)
    
    return fun


# Generate Tisserand-Poincare map figure
fig_tp = plt.figure()
ax_tp = plt.gca()
plt.xlabel('$r_a$ [km]')
plt.ylabel('$r_p$ [km]')
plt.grid()


# Plot orbit data on TP-map, colored by Delta-V magnitude
p_tp = ax_tp.scatter(ra_saturn_centric_save[ra_saturn_centric_save > 0],
                     rp_saturn_centric_save[ra_saturn_centric_save > 0], c=delta_v_mags[ra_saturn_centric_save > 0]*env.adim_v, s=2)  
fig_tp.colorbar(p_tp,ax=ax_tp,orientation='vertical',label='$\Delta V$ [km/s]')


# Plot Tisserand-Poincare lines at desired V_infinity values
v_infs = [0.2,0.3,0.4, 0.5, 0.6, 0.75]

for v_inf in v_infs:
    rp_range = np.linspace(170000,240000,101)
    ra_list = tp_graph_test(v_inf,rp_range)
    ax_tp.plot(ra_list,rp_range,alpha=0.5, label=v_inf, color = 'teal')

labellines.labelLines(ax_tp.get_lines())

ax_tp.set_ylim(180000, 250000)
plt.axis('equal')



# Compute closest point to desired v_infinity curve, and compute cost to transfer
# from endgame conditions to approach trajectory
v_inf_to_target = 0.2

# Compute v_infinity curve points
rp_range = np.linspace(170000,240000,101)
ra_list = tp_graph_test(v_inf_to_target,rp_range)
import scipy.interpolate


# Here we will compute the Delta-V for two different scenarios for each 
# insertion maneuver delta-V magnitude
# First, we will compute the Delta-V required to modify the r_a to transfer to 
# the v_infinity curve
# Second, we will compute the Delta-V required to modify the r_p to transfer to
# the v_infinity curve
# The minimum total delta-V for each insertion maneuver magnitude will correspond
# to the minimum of the two methods
    

# Compute DV to modify r_a at periapse (Saturn-centric orbit)
# Create function to interpolate r_a value along v_infinity curve given r_p
interp_vinf = scipy.interpolate.CubicSpline(rp_range,ra_list)
delta_v_tot_peri = np.zeros(rp_saturn_centric_save.size)
delta_v_insert = np.zeros(rp_saturn_centric_save.size)
delta_v_transfer = np.zeros(rp_saturn_centric_save.size)
dist_from_curve = np.zeros(rp_saturn_centric_save.size)

for i in range(delta_v_mags.size):
    # Interpolate value of r_a along v_infinity curve at associated r_p
    ra_pred = interp_vinf(rp_saturn_centric_save[i])
    
    # Compute distance from v_infinity curve along r_a axis
    dist_from_curve[i] = np.abs(ra_pred - ra_saturn_centric_save[i])

    # Compute delta-V to transfer from v_infinity curve to transfer orbit
    ra_depart = ra_pred
    ra_target = ra_saturn_centric_save[i]
    rp_maneuver = rp_saturn_centric_save[i]   
    a_target = (ra_target + rp_maneuver)/2
    a_depart = (ra_depart + rp_maneuver)/2
    v_target = np.sqrt(Saturn.GM*(2/rp_maneuver - 1/a_target))
    v_depart = np.sqrt(Saturn.GM*(2/rp_maneuver - 1/a_depart))
    delta_v_hoh = (v_target - v_depart)
    
    # Convert insertion maneuver delta-V to Cartesian units
    delta_v_min = delta_v_mags[i]*env.adim_v
    
    # Compute total delta-V (transfer + insertion maneuvers)
    delta_v_insert[i] = np.abs(delta_v_min)
    delta_v_transfer[i] = np.abs(delta_v_hoh)
    delta_v_tot_peri[i] = delta_v_insert[i] + delta_v_transfer[i]
    
    
    
# Compute DV to modify r_p at apoapse (Saturn-centric orbit)
# Create function to interpolate r_p value along v_infinity curve given r_a
interp_vinf = scipy.interpolate.CubicSpline(ra_list,rp_range)
delta_v_tot_apo = np.zeros(rp_saturn_centric_save.size)
delta_v_insert_apo = np.zeros(rp_saturn_centric_save.size)
delta_v_transfer_apo = np.zeros(rp_saturn_centric_save.size)
dist_from_curve_apo = np.zeros(rp_saturn_centric_save.size)

for i in range(delta_v_mags.size):
    # Interpolate value of r_p along v_infinity curve at associated r_a
    rp_pred = interp_vinf(ra_saturn_centric_save[i])
    
    # Compute distance from v_infinity curve along r_p axis
    dist_from_curve_apo[i] = np.abs(rp_pred - rp_saturn_centric_save[i])
        
    # Compute delta-V to transfer from v_infinity curve to transfer orbit
    rp_depart = rp_pred
    rp_target = rp_saturn_centric_save[i]
    ra_maneuver = ra_saturn_centric_save[i]
    a_target = (rp_target + ra_maneuver)/2
    a_depart = (rp_depart + ra_maneuver)/2
    v_target = np.sqrt(Saturn.GM*(2/ra_maneuver - 1/a_target))
    v_depart = np.sqrt(Saturn.GM*(2/ra_maneuver - 1/a_depart))
    delta_v_hoh = (v_target - v_depart)
    
    # Convert insertion maneuver delta-V to Cartesian units
    delta_v_min = delta_v_mags[i]*env.adim_v
    
    # Compute total delta-V (transfer + insertion maneuvers)
    delta_v_insert_apo[i] = np.abs(delta_v_min)
    delta_v_transfer_apo[i] = np.abs(delta_v_hoh)
    delta_v_tot_apo[i] = np.abs(delta_v_min) + np.abs(delta_v_hoh)
    

    
# Compute absolute minimum delta_V to transfer from desired v_infinity value and
# insert into science orbit  
delta_v_tot_min = np.min([delta_v_tot_apo[delta_v_tot_apo>0], delta_v_tot_peri[delta_v_tot_peri>0]])

# Plot total transfer delta-V as a function of the insertion maneuver delta-V  
plt.figure()
plt.plot(delta_v_mags*env.adim_v, delta_v_tot_peri, label = 'Total Delta V (lower apoapse)')
plt.plot(delta_v_mags*env.adim_v, delta_v_tot_apo, label = 'Total Delta V (raise periapse)')
plt.legend()
plt.grid()

plt.xlabel('$\Delta V$ applied at insertion [km/s]')
plt.ylabel('Total transfer $\Delta V$ [km/s]')

