# This file contains functions to interpolate specific physical properties
# within a family of periodic orbits
# Author: Spencer Boone

import numpy as np
from single_shooter import correct_ics, map_vars_to_index


def interpolate_orbit(mu, orbit_family_states, orbit_family_periods, target_var,
                      target_value, free_vars, constraint_vars, half_period = 1):
    
    # Compute bounds for specific variable
    target_var_index = map_vars_to_index(target_var)
    
    if target_var_index < 6: # State component
        target_var_bounds = [np.min(orbit_family_states[:,target_var_index]), np.max(orbit_family_states[:,target_var_index])]
    elif target_var_index == 6: # Target periods
        target_var_bounds = [np.min(orbit_family_periods), np.max(orbit_family_periods)]
        
    if target_value < target_var_bounds[0] or target_value > target_var_bounds[1]:
        print("Target value not within family bounds, cannot interpolate")
        
    else:
        
        # Find closest value in family to target value and set as initial guess
        if target_var_index < 6:
            closest_ind = np.abs(orbit_family_states[:,target_var_index] - target_value).argmin()
        elif target_var_index == 6:
            closest_ind = np.abs(orbit_family_periods[target_var_index] - target_value).argmin()
            
        X_i_guess = orbit_family_states[closest_ind,:]
        period_guess = orbit_family_periods[closest_ind]
        
        # Fix target value
        if target_var_index < 6:
            X_i_guess[target_var_index] = target_value
        elif target_var_index == 6:
            period_guess = target_value
        
        # Use half period to enforce constraints if orbit is symmetric along x-z plane
        if half_period == 1:
            period_guess = period_guess/2
        
        # Correct initial conditions with fixed target parameter
        X_initial_corrected, period_corrected, flag = correct_ics(X_i_guess, mu, period_guess, free_vars, constraint_vars)
        
        if half_period == 1:
            period_corrected *= 2
            
        # Output orbit with fixed target parameter
        return X_initial_corrected, period_corrected, flag
    
    
    
    
    