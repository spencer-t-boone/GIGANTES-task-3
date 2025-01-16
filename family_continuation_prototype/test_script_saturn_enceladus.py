# Test script for generating families of periodic orbits

import numpy as np
from cr3bp_propagate import propagate
from single_shooter import correct_ics
from cr3bp_functions import add_stm_i
from continuation import continue_family_npc, continue_family_palc
from cr3bp_plotting import plot_family

mu = 1.9011497893288988e-7


# Tests: 
# L1 halo orbits using pseudo arclength continuation
# L2 halo orbits
# Butterfly orbits using pseudo arclength continuation
# Butterfly orbits using natural parameter continuation along z parameter

# L1 halo orbit

test = "butterfly"

if test == "L1 halo":

    X_i = np.array([ 9.97234212e-01, -7.28454572e-26, -1.37102221e-03, -9.09779652e-17,
           -5.67257640e-03,  5.23832922e-16])
    t_f_guess = 3.05133468e+00
    
    free_vars = ["x", "z", "ydot", "t"]
    constraints = ["y", "xdot", "zdot"]
    
    # y_int = propagate(X_i, mu, t_f_guess)

    # y_int_stm = propagate(add_stm_i(X_i), mu, t_f_guess, with_stm = 1)
    
    continuation_var = []
    step = -0.001
    orbit_family_states, orbit_family_periods, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                     constraints, step, N_orbits_max=1000, half_period = 0, use_xz_symmetry = 0)

    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 50)
        
        

elif test == "L2 halo":
    X_i = np.array([ 1.00290351e+00, -7.49987773e-14, -1.21076489e-03,
      1.95365532e-13,  5.27054580e-03, -8.69803361e-14])
    t_f_guess = 3.07287146e+00
    
    free_vars = ["x", "z", "ydot", "t"]
    constraints = ["y", "xdot", "zdot"]
    
    y_int = propagate(X_i, mu, t_f_guess)

    y_int_stm = propagate(add_stm_i(X_i), mu, t_f_guess, with_stm = 1)
    
    continuation_var = []
    step = -0.001
    orbit_family_states, orbit_family_periods, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                     constraints, step, N_orbits_max=1000, half_period = 0, use_xz_symmetry = 0)
    
    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 50)
        
        
elif test == 'butterfly':
    X_i = np.array([ 9.98431607e-01,  8.18082304e-21,  3.79041307e-03,  3.51117145e-16,
           -3.43371947e-03, -6.09114193e-15])
    t_f_guess = 4.69397956e+00
    
    free_vars = ["x", "z", "ydot", "t"]
    constraints = ["y", "xdot", "zdot"]
    
    y_int = propagate(X_i, mu, t_f_guess)

    y_int_stm = propagate(add_stm_i(X_i), mu, t_f_guess, with_stm = 1)
    
    continuation_var = []
    step = -0.01
    orbit_family_states, orbit_family_periods, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                     constraints, step, N_orbits_max=150, half_period = 0, use_xz_symmetry = 0)
    
    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 10)
        
        
elif test == 'butterfly NPC':
    X_i = np.array([ 9.98431607e-01,  8.18082304e-21,  3.79041307e-03,  3.51117145e-16,
           -3.43371947e-03, -6.09114193e-15])
    t_f_guess = 4.69397956e+00
    
    continuation_var = ["z"]
    free_vars = ["x", "ydot", "t"]
    constraints = ["y", "xdot", "zdot"]
    
    y_int = propagate(X_i, mu, t_f_guess)

    y_int_stm = propagate(add_stm_i(X_i), mu, t_f_guess, with_stm = 1)
    
    #continuation_var = []
    step = 0.00001
    orbit_family_states, orbit_family_periods, flag = continue_family_npc(X_i, mu, t_f_guess, continuation_var, free_vars, 
                                                                     constraints, step, N_orbits_max=80, half_period = 0, use_xz_symmetry = 0)
    
    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 2)
        
        


