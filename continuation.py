# This file contains a number of functions to perform numerical continuation
# of a family of orbits in the CR3BP (or any dynamic system)
# Both natural parameter and pseudo arc-length continuation methods are 
# implemented
# Author: Spencer Boone
import numpy as np
from single_shooter import correct_ics, map_vars_to_index
import copy
import scipy.linalg as sp



# Continue a family of orbits using the natural parameter continuation method
def continue_family_npc(X_initial, mu, period_initial, continuation_var, free_vars, constraint_vars, step, N_orbits_max = 50, half_period = 1,
                        use_xz_symmetry = 1, line_search_param = 0.75, min_step = 1e-6):
    
    # Convert variable along which to perform continuation to associated integer value
    continuation_var_index = map_vars_to_index(continuation_var)
    
    print('Continuting family along', continuation_var[0], 'parameter starting at ', continuation_var[0], '=', str(X_initial[continuation_var_index]))
    
    
    
    # Correct initial conditions in case they are not initially periodic
    X_initial_corrected, period_initial_corrected, flag = correct_ics(X_initial, mu, period_initial, free_vars, constraint_vars,
                                                                half_period = half_period, use_xz_symmetry = use_xz_symmetry)
    
    if flag == 0:
        print('ERROR: Initial conditions did not converge to periodic orbit')
    
    
    # Set up continuation loop
    orbit_family_states = np.zeros((N_orbits_max, 6))
    orbit_family_states[0,:] = X_initial_corrected
    orbit_family_periods = np.zeros(N_orbits_max)
    orbit_family_periods[0] = period_initial_corrected
    
    # Set up line search 
    max_step = np.abs(step)
    
    cont_success_flag = 1
    
    # Begin natural parameter continuation for N_orbits_max orbits
    for i in range(1,N_orbits_max):
        converged = 0
        
        while converged == 0:
            
            # New initial guess, Move next guess by 'step' along the desired parameter
            X_cont_guess = copy.deepcopy(orbit_family_states[i-1,:])
            X_cont_guess[continuation_var_index] += step 
            print('continuing orbit ', i, ', ', continuation_var[0], '=', str(X_cont_guess[continuation_var_index]))
            period_cont_guess = orbit_family_periods[i-1]
            
            
            # Correct new initial guess
            X_cont_corrected, period_cont_corrected, flag = correct_ics(X_cont_guess, mu, period_cont_guess, free_vars, constraint_vars,
                                                                  half_period = half_period, use_xz_symmetry = use_xz_symmetry)
            
            # Line search logic to reduce step size if failed (still in development)
            if flag == 0:
                step = step*line_search_param
                print('Continuation iteration failed, reducing step size')
                
                if np.abs(step) < min_step:
                    print('Continuation failed, try with different parameters')
                    cont_success_flag = 0
                    break
            else:
                
                # If next guess converged, save orbit and move on to next orbit in family
                orbit_family_states[i,:] = X_cont_corrected
                orbit_family_periods[i] = period_cont_corrected
                converged = 1
                
                # Increase line search parameter if it is lower than max_step
                if np.abs(step) < max_step:
                    step = np.sign(step)*np.min((np.abs(step)/line_search_param, max_step))
                    
                print('Converged')
                
        if cont_success_flag == 0:
            break
        
        
    return orbit_family_states, orbit_family_periods, cont_success_flag




# Continue a family of orbits using the psuedo arclength continuation method
def continue_family_palc(X_initial, mu, period_initial, free_vars, constraint_vars, step, N_orbits_max = 50, half_period = 1,
                        use_xz_symmetry = 1, line_search_param = 0.75, min_step = 1e-6):
    
   
    #continuation_var_index = map_vars_to_index(continuation_var)
    free_vars_index = map_vars_to_index(free_vars)
    
    print('Continuting family using Pseudo arc-length continuation')
    
    # Correct initial conditions in case they are not initially periodic
    X_initial_corrected, period_initial_corrected, flag, DF = correct_ics(X_initial, mu, period_initial, free_vars, constraint_vars,
                                                                half_period = half_period, use_xz_symmetry = use_xz_symmetry,
                                                                output_DF = 1)
    
    if flag == 0:
        print('ERROR: Initial conditions did not converge to periodic orbit')
    
    
    # Continuation loop
    orbit_family_states = np.zeros((N_orbits_max, 6))
    orbit_family_states[0,:] = X_initial_corrected
    orbit_family_periods = np.zeros(N_orbits_max)
    orbit_family_periods[0] = period_initial_corrected
    
    # Set up line search 
    max_step = np.abs(step)
    
    # Set up continuation loop
    cont_success_flag = 1
    null_vec = np.zeros(4)
    
    
    # Begin family continuation loop
    for i in range(1,N_orbits_max):
        converged = 0
        
        # logic to make sure that null vector is always in the correct direction (and not moving backwards along the family)
        if i ==1:
            prev_null_vec = np.zeros(4)
        else:
            prev_null_vec = copy.deepcopy(null_vec)
        
        
        while converged == 0:
            print('continuing orbit ', i)
            
            # New initial guess
            null_vec = (sp.null_space(DF)).flatten()
            
            # Check if null vector is in same direction as previous null vector, if not, reverse the sign
            if i >= 2:               
                null_vecs_dot = np.dot(null_vec, prev_null_vec)
                null_vec = null_vec * np.sign(null_vecs_dot)
                
            
            # Compute next guess for family member along direction of null vector
            X_cont_guess = copy.deepcopy(orbit_family_states[i-1,:])
            next_guess = step*null_vec
            X_cont_guess[free_vars_index[:-1]] += next_guess[:-1]
            period_cont_guess = orbit_family_periods[i-1] + next_guess[-1]
            
        
            # Arguments to be passed to corrector for PALC constraint
            palc_args = [step, null_vec, orbit_family_states[i-1,:], orbit_family_periods[i-1]]
            
            # Correct new guess
            X_cont_corrected, period_cont_corrected, flag, DF = correct_ics(X_cont_guess, mu, period_cont_guess, free_vars, constraint_vars,
                                                                  half_period = half_period, use_xz_symmetry = use_xz_symmetry,
                                                                  palc_continuation = 1, palc_args = palc_args)
            
            # Line search if guess did not converge (still needs testing)
            if flag == 0:
                step = step*line_search_param
                print('Continuation iteration failed, reducing step size')
                
                if np.abs(step) < min_step:
                    print('Continuation failed, try with different parameters')
                    cont_success_flag = 0
                    break
            else:
                
                # Save converged orbit
                orbit_family_states[i,:] = X_cont_corrected
                orbit_family_periods[i] = period_cont_corrected
                converged = 1
                
                # Increase line search parameter if it is lower than max_step
                if np.abs(step) < max_step:
                    step = np.sign(step)*np.min((np.abs(step)/line_search_param, max_step))
                    
                print('Converged')
                
        if cont_success_flag == 0:
            break
        
        
        
        
    return orbit_family_states, orbit_family_periods, cont_success_flag
        