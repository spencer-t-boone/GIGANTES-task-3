# This file contains the functions required to perform a single shooting 
# algorithm with a variety of possible inputs
# Author: Spencer Boone


import numpy as np
from cr3bp_propagate import propagate
from cr3bp_functions import add_stm_i, state_derivs_cr3bp
from copy import deepcopy

"""
# Main differential corrector function
Inputs:
X_guess ---------- initial state guess to be corrected
mu --------------- mass parameter
t_guess ---------- initial guess for orbit period to be corrected
free_vars -------- list of strings containing free variables for correction
constraint_vars -- list of strings containing constraint variables for correction
N_max ------------ maximum correction iterations
corr_tol --------- tolerance for correction procedure
half_period ------ target conditions at 'half-period' (useful for symmetric orbits such as halo and Lyapunov)
palc_continuation  flag indicating that pseudo arc-length constraint should be added
palc_args -------- list containing terms needed for pseudo arc-length conrinuation
output_DF -------- flag indicating whether the DF matrix should be outputted


Outputs:
X_guess ---------- converged state
period_guess ----- converged period
converged_flag --- flag indicating whether tolerances are met
DF --------------- derivative of free variable vector with respect to constraint vector
"""
    

def correct_ics(X_guess, mu, t_guess, free_vars, constraint_vars, N_max = 50, 
                corr_tol = 1e-12, palc_continuation = 0, palc_args = None, output_DF = 0):
    
    
        
    # Convert free variable and constraint variable strings to vectors
    free_vars_index = map_vars_to_index(free_vars)
    n_free = len(free_vars_index)
    if "t" in free_vars:
        n_free_in_stm = n_free - 1
    else:
        n_free_in_stm = n_free
        
    constraints_index = map_vars_to_index(constraint_vars)
    n_constraint = len(constraints_index)
    X_constraints_desired = X_guess[constraints_index]
    
    
        
    
    # Propagate initial conditions and STM to t_guess       
    y_prop_stm = propagate(add_stm_i(X_guess), mu, t_guess, with_stm = 1)
    
    # Extract state and STM from final time
    X_f_guess = y_prop_stm.y[:6,-1]
    stm_f = y_prop_stm.y[6:,-1].reshape((6,6))
    period_guess = y_prop_stm.t[-1]
        
        
    # Extract pseudo arc-length continuation parameters if using PALC
    if palc_continuation == 1:
        step = palc_args[0]
        null_vec = palc_args[1]
        X_prev_converged = palc_args[2]
        t_prev = palc_args[3]
        
        
        
    
    # Compute dXf for initial guess (achieved - desired)    
    dXf_desired = X_f_guess[constraints_index] - X_constraints_desired
    X_diff = np.linalg.norm(dXf_desired)
        
    
    # Indices of STM to be used in targeter
    stm_col_index = free_vars_index[free_vars_index < 6]
    stm_row_index = constraints_index[constraints_index < 6]
    
    # Initialize constraints gradient matrix
    DF = np.zeros((n_constraint, n_free))
    
    # Start correction count
    N_count = 0
    
    # Compute DF before loop in case constraints are already satisfied
    stm_temp = stm_f[stm_row_index,:]
    DF[:n_constraint, :n_free_in_stm] = stm_temp[:,stm_col_index]
    
    # Logic if time is a free variable
    if "t" in free_vars:
        DF[:,-1] = state_derivs_cr3bp(X_f_guess, mu)[constraints_index]

    # Run correction
    while X_diff > corr_tol and N_count < N_max:
        
        #print(X_diff)

        # Compute first order correction with Newton's method
        
        stm_temp = stm_f[stm_row_index,:]
        DF[:n_constraint, :n_free_in_stm] = stm_temp[:,stm_col_index]
        
        # Logic if time is a free variable
        if "t" in free_vars:
            DF[:,-1] = state_derivs_cr3bp(X_f_guess, mu)[constraints_index]
        
            
        # Compute additional constraint for PALC
        if palc_continuation == 1:
            g_add = np.dot(np.hstack([X_guess[stm_col_index] - X_prev_converged[stm_col_index], t_guess - t_prev]), null_vec) - step
            g_vec = np.hstack([dXf_desired, g_add])        
            DG = np.vstack([DF, null_vec])
           
            
        # Compute update
        if palc_continuation == 0:
            if n_free == n_constraint: 
                dX0_update = np.linalg.inv(DF)@dXf_desired
            elif n_free > n_constraint:
                dX0_update = np.linalg.pinv(DF)@dXf_desired
            else:
                dX0_update = np.linalg.lstsq(DF, dXf_desired)[0]
                
        elif palc_continuation == 1:
            if np.shape(DG[0]) == np.shape(DG[1]):
                dX0_update = np.linalg.inv(DG)@g_vec
            elif np.shape(DG[0]) > np.shape(DG[1]):
                dX0_update = np.linalg.pinv(DG)@g_vec
            else:
                dX0_update = np.linalg.lstsq(DG, g_vec)[0]
                
                
        #print(dX0_update)
            
        # Apply update       
        if "t" in free_vars:
            dX0_update_apply = dX0_update[:-1]
            t_guess += -dX0_update[-1]
        else:
            dX0_update_apply = dX0_update
        
        
        X_guess[stm_col_index] = X_guess[stm_col_index] - dX0_update_apply
        
        
        # Propagate updated guess and STM to new t_f
        y_prop_stm = propagate(add_stm_i(X_guess), mu, t_guess, with_stm = 1)
        
        # Propagate state and STM to final time
        X_f_guess = y_prop_stm.y[:6,-1]
        stm_f = y_prop_stm.y[6:,-1].reshape((6,6))
        period_guess = y_prop_stm.t[-1]
        
            
        # Compute dXf for updated guess  
        dXf_desired = X_f_guess[constraints_index] - X_guess[constraints_index]
        X_diff = np.linalg.norm(dXf_desired) # Will exit while loop if below tol
        
        N_count += 1
    
        
    # Set flag to indicate convergence
    if X_diff > corr_tol:
        converged_flag = 0
    else:
        converged_flag = 1
        
      
    
    
    if palc_continuation == 0 and output_DF == 0:
        return X_guess, period_guess, converged_flag
    else:
        return X_guess, period_guess, converged_flag, DF



    

    
    
# Map variable strings to specific integer values
def map_vars_to_index(vars_list = None):
    variable_dict = {"x": 0, "y": 1, "z": 2, "xdot": 3, "ydot": 4, "zdot": 5, "t": 6, "jc": 7}
    vars_index = []

    if vars_list is not None:
        for i in range(len(vars_list)):
            # Check and handle KeyError
            try:
                vars_index.append(variable_dict[vars_list[i]])
            except KeyError:
                print(vars_index[i], "is not a valid free variable")
                vars_index = 0
                break

    # Sort the indices to handle variables passed in any order
    vars_index = sorted(vars_index)

    return np.array(vars_index)
