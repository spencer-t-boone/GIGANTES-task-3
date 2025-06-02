# Test script for generating families of periodic orbits

import numpy as np
from cr3bp_propagate import propagate, event_impact_secondary
from single_shooter import correct_ics
from cr3bp_functions import add_stm_i
from continuation import continue_family_npc, continue_family_palc
from cr3bp_plotting import plot_family
from family_interpolation import interpolate_orbit
from bifurcation_functions import plot_broucke_diagram, detect_bifurcations_broucke
from copy import deepcopy


# Saturn-Enceladus CR3BP constants
mu = 1.9011497893288988e-7
R_enc_km = 252.1
r_enc = 238400
R_enc = R_enc_km/r_enc


# Validation test cases: 
# L2 halo orbits using pseudo arclength continuation
# L1 halo orbits using natural parameter continuation
# Butterfly orbits using natural parameter continuation along period
# Continue butterfly orbit from Earth-Moon system into Saturn-Enceladus system using NPC along mu
# Continue L2 Lyapunov family, detect bifurcation w/ halo orbits and continue halo orbit family from bifurcation
# Continue L2 halo orbit family, detect first period-doubling and tripling bifurcations, and continue these families
# L1 halo orbit family using PALC in Earth-Moon system

# Additional test cases:
# L1 halo orbits using pseudo arclength continuation
# Butterfly orbits using pseudo arclength continuation
# Butterfly orbits using natural parameter continuation along z parameter
# L1 Period-3 halo orbits using PALC
# Continue L2 halo orbit family, detect first period-doubling bifurcation, and continue this family
# Continue L2 halo orbit family, detect first period-tripling bifurcation, and continue this family

# L1 halo orbit
validation_case_options = ["L2 halo",
                     "L1 halo NPC",
                     "butterfly NPC period",
                     "butterfly NPC mu",
                     "L2 lyapunov",
                     "L2 bifurcations",
                     "L1 halo PALC earthmoon"]

additional_test_case_options = ["L1 halo",
                     "butterfly",
                     "butterfly NPC",
                     "L1 p3-halo",                    
                     "L2 period doubling",
                     "L2 period tripling",
                     "lyapunov to butterfly"]


test = validation_case_options[3]
test = additional_test_case_options[6]


# -----------------------------------------------------------------------------
# L1 halo orbit test case using PALC
# -----------------------------------------------------------------------------
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
        
        # Interpolate orbit with specific z_i value (30 km above surface of Enceladus at periapse)
        target_var = ["z"]
        free_vars = ["x", "ydot", "t"]
        target_value = -(R_enc + 30/r_enc)
        orbit_interp_state, orbit_interp_period, flag = interpolate_orbit(mu, orbit_family_states, orbit_family_periods,
                                                                          target_var, target_value,
                                                                          free_vars, constraints, half_period = 1)
        

        
        
        
# -----------------------------------------------------------------------------
# L1 halo orbit test case using NPC
# -----------------------------------------------------------------------------
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
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 5, frame = 'sec-centric', 
                    R_sec = R_enc_km, r_sec = r_enc)
        
        
# -----------------------------------------------------------------------------
# L2 halo orbit test case using PALC
# -----------------------------------------------------------------------------
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
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 5, frame = 'sec-centric', 
                    R_sec = R_enc_km, r_sec = r_enc)
        
        # Interpolate orbit with specific z_i value (30 km above surface of Enceladus at periapse)
        orbit_peri_states = np.zeros(orbit_family_states.shape)
        for i in range(orbit_family_states.shape[0]):
            y_prop = propagate(orbit_family_states[i,:], mu, orbit_family_periods[i]/2)
            orbit_peri_states[i,:] = y_prop.y[:6,-1]
        
        target_var = ["z"]
        free_vars = ["x", "ydot", "t"]
        target_value = -(R_enc + 30/r_enc)
        orbit_interp_state, orbit_interp_period, flag = interpolate_orbit(mu, orbit_peri_states, orbit_family_periods,
                                                                          target_var, target_value,
                                                                          free_vars, constraints, half_period = 1)
        
 
# -----------------------------------------------------------------------------
# Butterfly orbit test case using PALC
# -----------------------------------------------------------------------------
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
        
        # Interpolate orbit with specific x_i value
        target_var = ["x"]
        free_vars = ["z", "ydot", "t"]
        target_value = 0.998 # target x = 0.998
        orbit_interp_state, orbit_interp_period, flag = interpolate_orbit(mu, orbit_family_states, orbit_family_periods,
                                                                          target_var, target_value,
                                                                          free_vars, constraints, half_period = 1)
        
        # Interpolate orbit with specific period
        target_var = ["t"]
        free_vars = ["x", "z", "ydot"]
        target_value = 4.0 # target period = 4.0 nondim units
        orbit_interp_state_2, orbit_interp_period_2, flag = interpolate_orbit(mu, orbit_family_states, orbit_family_periods,
                                                                          target_var, target_value,
                                                                          free_vars, constraints, half_period = 1)


# -----------------------------------------------------------------------------
# Butterfly orbit test case using NPC along z parameter
# -----------------------------------------------------------------------------        
        
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
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 10, frame = 'sec-centric', 
                    R_sec = R_enc_km, r_sec = r_enc)
        
         

# -----------------------------------------------------------------------------
# Butterfly orbit test case using NPC along period
# -----------------------------------------------------------------------------        
elif test == 'butterfly NPC period':
    X_i = np.array([ 9.98431607e-01, 0, 3.79041307e-03, 0, -3.43371947e-03, 0])
    t_f_guess = 4.69397956e+00
    
    continuation_var = ["t"]
    free_vars = ["x", "z", "ydot"]
    constraints = ["y", "xdot", "zdot"]
    
    step = -0.02
    event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
    
    orbit_family_states, orbit_family_periods, flag = continue_family_npc(X_i, mu, t_f_guess, continuation_var, free_vars, 
                                                                     constraints, step, N_orbits_max=80, half_period = 1,
                                                                     event_stop = event_impact_enceladus)
    
    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 2, frame = 'sec-centric', 
                    R_sec = R_enc_km, r_sec = r_enc)
        
        
        
        
# -----------------------------------------------------------------------------
# L1 period-3 halo orbit test case using PALC
# ----------------------------------------------------------------------------- 
elif test == 'L1 p3-halo':
    X_i = np.array([0.99562262, 0, 0.00344389, 0, 0.01031047, 0])
    t_f_guess = 8.860077536892092

    
    free_vars = ["x", "z", "ydot", "t"]
    constraints = ["y", "xdot", "zdot"]
    
    continuation_var = []
    step = -0.002
    event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
    
    orbit_family_states, orbit_family_periods, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                     constraints, step, N_orbits_max=560, half_period = 1)
    
    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 40, frame = 'sec-centric', 
                    R_sec = R_enc_km, r_sec = r_enc)        


# ---------------------------------------------------------------------------------------------------
# L2 lyapunov orbit test case using PALC, detect halo orbit bifurcation and continue resulting family
# ---------------------------------------------------------------------------------------------------
elif test == "L2 lyapunov":  
    X_i = np.array([ 1.00401614e+00, -3.12950469e-20,  7.58158917e-29,  2.24329841e-14,
           -1.60598059e-04,  4.31816369e-28])
    t_f_guess = 3.04165111e+00
    
    free_vars = ["x", "z", "ydot", "t"]
    constraints = ["y", "xdot", "zdot"]
    
    continuation_var = []
    step = 1e-3
    event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
    
    orbit_family_states, orbit_family_periods, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                     constraints, step, N_orbits_max=1000, half_period = 1,
                                                                     event_stop = event_impact_enceladus)
    
    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 5, frame = 'sec-centric', R_sec = R_enc_km, r_sec = r_enc)
        
        # Plot Broucke diagram
        fig_broucke = plot_broucke_diagram(orbit_family_states, orbit_family_periods, mu)
        
        # Detect bifurcations
        bif_types = ['tangent']
        bif_states, bif_periods, bif_cont_states, bif_cont_periods, bif_types = detect_bifurcations_broucke(orbit_family_states, orbit_family_periods, mu,
                                                                                                            free_vars, constraints, bif_types = bif_types)

    # Continue L2 halo orbit family from resulting tangent bifurcation
    X_i = bif_cont_states[0]
    t_f_guess = bif_cont_periods[0]
    step = -5e-4
    event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
    
    orbit_family_states, orbit_family_periods, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                     constraints, step, N_orbits_max=2000, half_period = 1,
                                                                     event_stop = event_impact_enceladus)



# ---------------------------------------------------------------------------------------------------------------------
# L2 period-doubling halo orbit test case using PALC, detect bifurcation from halo orbits and continue resulting family
# ---------------------------------------------------------------------------------------------------------------------
elif test == 'L2 period doubling':
    
    X_i = np.array([ 1.00446305e+00, 0,  2.59982886e-05, 0, -3.53884567e-03, 0])
    t_f_guess = 3.08985431e+00
    
    free_vars = ["x", "z", "ydot", "t"]
    constraints = ["y", "xdot", "zdot"]
    
    continuation_var = []
    step = -5e-4
    event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
    
    orbit_family_states, orbit_family_periods, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                     constraints, step, N_orbits_max=1600, half_period = 1,
                                                                     event_stop = event_impact_enceladus)
    
    
    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 50, frame = 'sec-centric', 
                    R_sec = R_enc_km, r_sec = r_enc)
        
        # Plot Broucke diagram
        fig_broucke = plot_broucke_diagram(orbit_family_states, orbit_family_periods, mu)
        
        # Detect bifurcations
        bif_types = ['period_doubling']
        bif_states, bif_periods, bif_cont_states, bif_cont_periods, bif_types = detect_bifurcations_broucke(orbit_family_states, orbit_family_periods, mu,
                                                                                                            free_vars, constraints, bif_types = bif_types, bif_dir_step = 1e-5)
        
        # Continue first period-doubling family
        X_i = bif_cont_states[0]
        t_f_guess = bif_cont_periods[0]
        step = 5e-4
        event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
        
        orbit_family_states_1, orbit_family_periods_1, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                         constraints, step, N_orbits_max=800, half_period = 1,
                                                                         event_stop = event_impact_enceladus)
        
        if flag == 1:
            plot_family(orbit_family_states_1, orbit_family_periods_1, mu, spacing = 50, frame = 'sec-centric', 
                        R_sec = R_enc_km, r_sec = r_enc)
            

# ---------------------------------------------------------------------------------------------------------------------
# L2 period-tripling halo orbit test case using PALC, detect bifurcation from halo orbits and continue resulting family
# Needs a denser computation of the L2 halo orbit family in order to properly detect the sensitive bifurcation
# ---------------------------------------------------------------------------------------------------------------------
elif test == 'L2 period tripling':

    X_i = np.array([ 1.00398819,  0.        ,  0.00327934,  0.        , -0.00561346,
            0.        ])
    t_f_guess = 2.990128883855193
    
    free_vars = ["x", "z", "ydot", "t"]
    constraints = ["y", "xdot", "zdot"]
    
    continuation_var = []
    step = -2e-5
    event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
    
    orbit_family_states, orbit_family_periods, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                     constraints, step, N_orbits_max=2500, half_period = 1,
                                                                     event_stop = event_impact_enceladus)
    
    
    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 50, frame = 'sec-centric', 
                    R_sec = R_enc_km, r_sec = r_enc)
        
        # Plot Broucke diagram
        fig_broucke = plot_broucke_diagram(orbit_family_states, orbit_family_periods, mu)
        
        # Detect bifurcations
        bif_types = ['period_tripling']
        bif_states, bif_periods, bif_cont_states, bif_cont_periods, bif_types = detect_bifurcations_broucke(orbit_family_states, orbit_family_periods, mu,
                                                                                                            free_vars, constraints, bif_types = bif_types, bif_dir_step = 1e-6,
                                                                                                            half_period = 1)
        
        # Continue first period-tripling family
        X_i = bif_cont_states[0]
        t_f_guess = bif_cont_periods[0]
        step = 2e-5
        event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
        
        orbit_family_states_1, orbit_family_periods_1, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                         constraints, step, N_orbits_max=2000, half_period = 1, max_step = step*5,
                                                                         event_stop = event_impact_enceladus)
        
        if flag == 1:
            plot_family(orbit_family_states_1, orbit_family_periods_1, mu, spacing = 100, frame = 'sec-centric', 
                        R_sec = R_enc_km, r_sec = r_enc)
            
# ---------------------------------------------------------------------------------------------------------------------
# Continue L2 halo orbit family, detect first period-doubling and period-tripling bifurcations, and continue 
# both resulting families
# ---------------------------------------------------------------------------------------------------------------------
elif test == 'L2 bifurcations':

    X_i = np.array([ 1.00398819,  0.        ,  0.00327934,  0.        , -0.00561346,
            0.        ])
    t_f_guess = 2.990128883855193
    
    free_vars = ["x", "z", "ydot", "t"]
    constraints = ["y", "xdot", "zdot"]
    
    continuation_var = []
    step = -2e-5
    event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
    
    orbit_family_states, orbit_family_periods, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                     constraints, step, N_orbits_max=10000, half_period = 1,
                                                                     event_stop = event_impact_enceladus)
    
    
    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 50, frame = 'sec-centric', 
                    R_sec = R_enc_km, r_sec = r_enc)
        
        # Plot Broucke diagram
        fig_broucke = plot_broucke_diagram(orbit_family_states, orbit_family_periods, mu)
        
        # Detect bifurcations
        print('Detecting bifurcations')
        bif_types = ['period_doubling','period_tripling']
        bif_states, bif_periods, bif_cont_states, bif_cont_periods, bif_types = detect_bifurcations_broucke(orbit_family_states, orbit_family_periods, mu,
                                                                                                            free_vars, constraints, bif_types = bif_types, bif_dir_step = 1e-6,
                                                                                                            half_period = 1)
        
        # Continue first period-doubling family
        X_i = bif_cont_states[0]
        t_f_guess = bif_cont_periods[0]
        step = 5e-4
        event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
        
        orbit_family_states_1, orbit_family_periods_1, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                         constraints, step, N_orbits_max=400, half_period = 1, max_step = step*5,
                                                                         event_stop = event_impact_enceladus)
        
        if flag == 1:
            plot_family(orbit_family_states_1, orbit_family_periods_1, mu, spacing = 20, frame = 'sec-centric', 
                        R_sec = R_enc_km, r_sec = r_enc)
        

        # Continue first period-tripling family
        X_i = bif_cont_states[2]
        t_f_guess = bif_cont_periods[2]
        step = 2e-5
        event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
        
        orbit_family_states_2, orbit_family_periods_2, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                         constraints, step, N_orbits_max=2000, half_period = 1, max_step = step*5,
                                                                         event_stop = event_impact_enceladus)
        
        if flag == 1:
            plot_family(orbit_family_states_2, orbit_family_periods_2, mu, spacing = 100, frame = 'sec-centric', 
                        R_sec = R_enc_km, r_sec = r_enc)
            
            
            
# -----------------------------------------------------------------------------
# Compute Saturn-Enceladus butterfly orbit family by continuing along mu from 
# Earth-Moon butterfly orbit family
# ----------------------------------------------------------------------------- 
elif test == 'butterfly NPC mu':
    X_i = np.array([1.0669, 0, 0.1619, 0, -0.0189, 0]) # Initial conditions in Earth-moon system
    t_f_guess = 3.2255
    mu_i =  0.0121505856
    
    # Compute continuation step
    target_mu = mu
    num_step = 500
    step = (mu_i/target_mu)**(1/num_step)
    
    mu_vector = np.zeros(num_step)
    for i in range(num_step):
        mu_vector[i] = mu_i/step**i

    continuation_var = ["mu"]
    free_vars = ["x", "z", "ydot", "t"]
    constraints = [ "y", "xdot", "zdot"]
    
    
    orbit_family_states, orbit_family_periods, flag = continue_family_npc(X_i, mu_i, t_f_guess, continuation_var, free_vars, 
                                                                     constraints, step, N_orbits_max = num_step, half_period = 1)

    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu_vector, spacing = 50, variable_mu = 1)
        
        
# -----------------------------------------------------------------------------
# Compute L1 halo orbit family in Eath-Moon system
# ----------------------------------------------------------------------------- 
elif test == 'L1 halo PALC earthmoon':
    X_i = np.array([0.8234, 0, 0.0224, 0, 0.1343, 0]) # Initial conditions in Earth-moon system
    t_f_guess = 2.7430
    mu =  0.0121505856
    R_moon = 1737.4
    
    event_impact_moon = lambda t, X: event_impact_secondary(t, X, mu, R_moon)
    
    
    free_vars = ["x", "z", "ydot", "t"]
    constraints = ["y", "xdot", "zdot"]
    
    continuation_var = []
    step = -10e-4
    event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
    
    orbit_family_states, orbit_family_periods, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                     constraints, step, N_orbits_max=1000, half_period = 1, 
                                                                     max_step = np.abs(step), event_stop = event_impact_moon)
    
    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 20, frame = 'sec-centric', 
                    R_sec = R_moon, r_sec = 384400)
        
        
# ---------------------------------------------------------------------------------------------------------------------
# Continue L2 Lyapunov orbit family in Earth-Moon system, detect bifurcation with halo orbits, 
# continue L2 halo orbit family, detect bifurcation with butterfly orbits (4th period-doubling), continue butterfly orbit family
# ---------------------------------------------------------------------------------------------------------------------
elif test == 'lyapunov to butterfly':
    mu =  0.0121505856
    R_moon = 1737.4
    r_moon = 384400
    
    X_i = np.array([ 1.12894, 0, 0, 0, 0.135192, 0])
    X_i = np.array([1.1762, 0, 0, 0, -0.1231, 0])
    t_f_guess = 3.3981
    
    free_vars = ["x", "z", "ydot", "t"]
    constraints = ["y", "xdot", "zdot"]
    
    continuation_var = []
    
    # Generate L2 Lyapunov orbit family
    step = 5e-3
    event_impact_moon = lambda t, X: event_impact_secondary(t, X, mu, R_moon)
    
    orbit_family_states, orbit_family_periods, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                     constraints, step, N_orbits_max=300, half_period = 1,
                                                                     event_stop = event_impact_moon)
    
    
    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 50, frame = 'sec-centric', 
                    R_sec = R_moon, r_sec = r_moon)
        
        # Plot Broucke diagram
        fig_broucke = plot_broucke_diagram(orbit_family_states, orbit_family_periods, mu)
        
        # Detect bifurcations
        bif_types = ['tangent']
        bif_states, bif_periods, bif_cont_states, bif_cont_periods, bif_types = detect_bifurcations_broucke(orbit_family_states, orbit_family_periods, mu,
                                                                                                            free_vars, constraints, bif_types = bif_types, bif_dir_step = 1e-5)
        
        # Continue L2 halo orbit family
        X_i = bif_cont_states[0]
        t_f_guess = bif_cont_periods[0]
        step = -5e-4
        event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
        
        orbit_family_states_1, orbit_family_periods_1, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                          constraints, step, N_orbits_max=2400, half_period = 1)
        
        if flag == 1:
            plot_family(orbit_family_states_1, orbit_family_periods_1, mu, spacing = 50, frame = 'sec-centric', 
                        R_sec = R_moon, r_sec = r_moon)
            
            # Plot Broucke diagram
            fig_broucke = plot_broucke_diagram(orbit_family_states_1, orbit_family_periods_1, mu)
            
            # Detect bifurcations
            bif_types = ['period_doubling']
            bif_states_1, bif_periods_1, bif_cont_states_1, bif_cont_periods_1, bif_types = detect_bifurcations_broucke(orbit_family_states_1, orbit_family_periods_1, mu,
                                                                                                                free_vars, constraints, bif_types = bif_types, bif_dir_step = 1e-4)
               
        
            # Continue butterfly orbit family (fourth period-doubling bifurcation)
            X_i = bif_cont_states_1[3]
            t_f_guess = bif_cont_periods_1[3]
            step = -5e-4
            event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
            
            orbit_family_states_2, orbit_family_periods_2, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                              constraints, step, N_orbits_max=1200, half_period = 1)
            
            
            if flag == 1:
                plot_family(orbit_family_states_2, orbit_family_periods_2, mu, spacing = 50, frame = 'sec-centric', 
                            R_sec = R_moon, r_sec = r_moon)
            
            
            
        
        