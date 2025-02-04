# Test script for generating families of periodic orbits

import numpy as np
from cr3bp_propagate import propagate, event_impact_secondary
from single_shooter import correct_ics
from cr3bp_functions import add_stm_i
from continuation import continue_family_npc, continue_family_palc
from cr3bp_plotting import plot_family


# Saturn-Enceladus CR3BP constants
mu = 1.9011497893288988e-7
R_enc_km = 252.1
r_enc = 238400
R_enc = R_enc_km/r_enc


# Tests: 
# L1 halo orbits using pseudo arclength continuation
# L2 halo orbits
# Butterfly orbits using pseudo arclength continuation
# Butterfly orbits using natural parameter continuation along z parameter

# L1 halo orbit

test = "L1 halo NPC"


if test == "L1 halo":
    
    X_i = np.array([ 9.96657814e-01, 0, -3.99088311e-05, 0, -3.84161723e-03, 0 ])
    t_f_guess =  3.07299480e+00
    
    
    free_vars = ["x", "z", "ydot", "t"]
    constraints = ["y", "xdot", "zdot"]
    
    continuation_var = []
    step = -2e-5
    event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
    
    orbit_family_states, orbit_family_periods, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                     constraints, step, N_orbits_max=2000, half_period = 1, 
                                                                     max_step = np.abs(step)*20, event_stop = event_impact_enceladus)

    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 100, frame = 'sec-centric', 
                    R_sec = R_enc_km, r_sec = r_enc)
        

if test == "L1 halo NPC":
    
    X_i = np.array([ 9.96657814e-01, 0, -3.99088311e-05, 0, -3.84161723e-03, 0 ])
    t_f_guess =  3.07299480e+00
    
    free_vars = ["z", "ydot", "t"]
    constraints = ["y", "xdot", "zdot"]
    
    
    continuation_var = ["x"]
    step = 1e-6
    event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
    
    orbit_family_states, orbit_family_periods, flag = continue_family_npc(X_i, mu, t_f_guess, continuation_var, free_vars, 
                                                                     constraints, step, N_orbits_max=3300, half_period = 1,
                                                                     event_stop = event_impact_enceladus)

    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 5)
        
        

elif test == "L2 halo":  
    X_i = np.array([ 1.00446305e+00, 0,  2.59982886e-05, 0, -3.53884567e-03, 0])
    t_f_guess = 3.08985431e+00
    
    free_vars = ["x", "z", "ydot", "t"]
    constraints = ["y", "xdot", "zdot"]
    
    continuation_var = []
    step = -5e-4
    event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
    
    orbit_family_states, orbit_family_periods, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                     constraints, step, N_orbits_max=2000, half_period = 1,
                                                                     event_stop = event_impact_enceladus)
    
    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 5)
        
        
elif test == 'butterfly':
    X_i = np.array([ 9.98431607e-01,  0,  3.79041307e-03,  0, -3.43371947e-03,0])
    t_f_guess = 4.69397956e+00
    
    free_vars = ["x", "z", "ydot", "t"]
    constraints = ["y", "xdot", "zdot"]
    
    continuation_var = []
    step = -0.01
    event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
    
    orbit_family_states, orbit_family_periods, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                     constraints, step, N_orbits_max=100, half_period = 1,
                                                                     event_stop = event_impact_enceladus)
    
    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 10, frame = 'sec-centric', 
                    R_sec = R_enc_km, r_sec = r_enc)
        
        
elif test == 'butterfly NPC':
    X_i = np.array([ 9.98431607e-01, 0, 3.79041307e-03, 0, -3.43371947e-03, 0])
    t_f_guess = 4.69397956e+00
    
    continuation_var = ["z"]
    free_vars = ["x", "ydot", "t"]
    constraints = ["y", "xdot", "zdot"]
    
    step = 1e-5
    event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
    
    orbit_family_states, orbit_family_periods, flag = continue_family_npc(X_i, mu, t_f_guess, continuation_var, free_vars, 
                                                                     constraints, step, N_orbits_max=80, half_period = 1,
                                                                     event_stop = event_impact_enceladus)
    
    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 2)
        
        


