# This file contains a number of functions to perform numerical continuation
# of a family of orbits in the CR3BP (or any dynamic system)
# Both natural parameter and pseudo arc-length continuation methods are 
# implemented
# Author: Spencer Boone
import numpy as np
from single_shooter import correct_ics, map_vars_to_index
from cr3bp_propagate import propagate
import copy
import scipy.linalg as sp



# Continue a family of orbits using the natural parameter continuation method
def continue_family_npc(X_initial, mu, period_initial, continuation_var, free_vars, constraint_vars, step, N_orbits_max = 50, half_period = 1,
                         line_search_param = 0.75, max_step = None, min_step = 1e-8, event_stop = None):
    
    # Convert variable along which to perform continuation to associated integer value
    continuation_var_index = map_vars_to_index(continuation_var)
    if continuation_var_index < 6:
        print('Continuting family along', continuation_var[0], 'parameter starting at ', continuation_var[0], '=', str(X_initial[continuation_var_index]))
    elif continuation_var_index == 6: # continuing along period
        print('Continuting family along', continuation_var[0], 'parameter starting at ', continuation_var[0], '=', str(period_initial))
        
    
    
    if half_period == 1:
        period_initial = period_initial/2
    
    
    # Correct initial conditions in case they are not initially periodic
    X_initial_corrected, period_initial_corrected, flag = correct_ics(X_initial, mu, period_initial, free_vars, constraint_vars)
    
    if flag == 0:
        print('ERROR: Initial conditions did not converge to periodic orbit')
    
    
    # Set up continuation loop
    orbit_family_states = np.zeros((N_orbits_max, 6))
    orbit_family_states[0,:] = X_initial_corrected
    orbit_family_periods = np.zeros(N_orbits_max)
    orbit_family_periods[0] = period_initial_corrected
    
    # Set up line search 
    if max_step == None:
        max_step = np.abs(step)
    
    cont_success_flag = 1
    exit_continuation = 0
    
    # Begin natural parameter continuation for N_orbits_max orbits
    for i in range(1,N_orbits_max):
        converged = 0
        
        while converged == 0:
            
            # New initial guess, Move next guess by 'step' along the desired parameter
            X_cont_guess = copy.deepcopy(orbit_family_states[i-1,:])
            period_cont_guess = copy.deepcopy(orbit_family_periods[i-1])
            
            if continuation_var_index < 6:
                X_cont_guess[continuation_var_index] += step
                print('continuing orbit ', i, ', ', continuation_var[0], '=', str(X_cont_guess[continuation_var_index]))
                            
            elif continuation_var_index == 6:
                
                if half_period == 1:
                    period_cont_guess *= 2
                    
                period_cont_guess += step
                print('continuing orbit ', i, ', ', continuation_var[0], '=', str(period_cont_guess))
                    
                if half_period == 1:
                    period_cont_guess = period_cont_guess/2
                    
                
                
            
            
            
            # Correct new initial guess
            X_cont_corrected, period_cont_corrected, flag = correct_ics(X_cont_guess, mu, period_cont_guess, free_vars, constraint_vars)
            
            # Line search logic to reduce step size if failed
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
                    print('Success, increasing step size to ', step)
                    
                print('Converged')
                
                # Check stop condition by integrating orbit with stopping event  
                y_prop = propagate(X_cont_corrected, mu, period_cont_corrected, with_stm = 0, events_fun = event_stop)
                
                # Check if stop criterion has been met
                if len(y_prop.t_events[0]) > 0:
                    
                    # If so, exit continuation loop and resize output arrays accordingly
                    print('Stopping criterion has been met, ending continuation')
                    exit_continuation = 1                    
                    orbit_family_states = orbit_family_states[:i,:]
                    orbit_family_periods = orbit_family_periods[:i]
                
        if cont_success_flag == 0:
            break
        
        if exit_continuation == 1:
            break
        
    if half_period == 1:
        orbit_family_periods = orbit_family_periods*2
        
        
    return orbit_family_states, orbit_family_periods, cont_success_flag




# Continue a family of orbits using the psuedo arclength continuation method
def continue_family_palc(X_initial, mu, period_initial, free_vars, constraint_vars, step, N_orbits_max = 50, half_period = 1,
                         line_search_param = 0.75, max_step = None, min_step = 1e-6, event_stop = None):
    
   
    #continuation_var_index = map_vars_to_index(continuation_var)
    free_vars_index = map_vars_to_index(free_vars)
    
    print('Continuting family using Pseudo arc-length continuation')
    
    # Use half period to enforce constraints if orbit is symmetric along x-z plane
    if half_period == 1:
        period_initial = period_initial/2
    
    
    # Correct initial conditions in case they are not initially periodic
    X_initial_corrected, period_initial_corrected, flag, DF = correct_ics(X_initial, mu, period_initial, free_vars, 
                                                                          constraint_vars, output_DF = 1)

        
    
    
    if flag == 0:
        print('ERROR: Initial conditions did not converge to periodic orbit')
    
    
    # Continuation loop
    orbit_family_states = np.zeros((N_orbits_max, 6))
    orbit_family_states[0,:] = X_initial_corrected
    orbit_family_periods = np.zeros(N_orbits_max)
    orbit_family_periods[0] = period_initial_corrected
    
    # Set up line search 
    if max_step == None:
        max_step = np.abs(step)
    
    # Set up continuation loop
    cont_success_flag = 1
    null_vec = np.zeros(4)
    exit_continuation = 0
    
    
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
                                                                  palc_continuation = 1, palc_args = palc_args)
            
            # Line search if guess did not converge
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
                    print('Success, increasing step size to ', step)
                    
                print('Converged')
                
                # Check stop condition by integrating orbit with stopping event  
                y_prop = propagate(X_cont_corrected, mu, period_cont_corrected, with_stm = 0, events_fun = event_stop)
                
                # Check if stop criterion has been met
                if event_stop != None:
                    if len(y_prop.t_events[0]) > 0:
                        
                        # If so, exit continuation loop and resize output arrays accordingly
                        print('Stopping criterion has been met, ending continuation')
                        exit_continuation = 1                    
                        orbit_family_states = orbit_family_states[:i,:]
                        orbit_family_periods = orbit_family_periods[:i]
                    
                
        if cont_success_flag == 0:
            break
        
        if exit_continuation == 1:
            break
            
    
    # Output full period if using only half-period for constraints
    if half_period == 1:
        orbit_family_periods = orbit_family_periods*2
        
        
    return orbit_family_states, orbit_family_periods, cont_success_flag




        