# This file contains a number of functions for propagating a state in the CR3BP
# Author: Spencer Boone


import numpy as np
from scipy.integrate import solve_ivp
from cr3bp_functions import state_derivs_cr3bp, A_cr3bp, add_stm_i


# Propagate a state in the CR3BP
def propagate(X_i, mu, t_prop, int_method = 'DOP853', abs_tol = 1e-12, rel_tol = 1e-12, 
              with_stm = 0, events_fun = None, t_eval = None):
    
    # Only save initial and final states if not specified by user
    if t_eval is None:
        t_eval = [0,t_prop]
    
    # Propagate state using SciPy's default integrator
    y_obj = solve_ivp(lambda t,X: cr3bp_ode(t, X, mu, with_stm), 
                           [0,t_prop], X_i, atol=abs_tol, rtol=rel_tol, method = int_method,
                           events = events_fun, t_eval = t_eval)
    
    # Extract STM if STM was integrated
    if with_stm == 1:
        y_obj.stm = y_obj.y[6:,:].reshape(6,6,y_obj.y.shape[1])
    
    return y_obj
    
    

# This function defines the dynamics ODE
def cr3bp_ode(t, X, mu, with_stm):
    
    # Compute state derivatives
    dX = state_derivs_cr3bp(X[0:6], mu)
    
    # Compute STM derivatives if desired
    if with_stm == 1:
        phi = X[6:].reshape((6,6))
        A = A_cr3bp(X, mu)
        
        phidot = A@phi
        
        dX = np.hstack([dX, phidot.reshape(-1)])
    
    return dX



# EVENT FUNCTIONS
# Event function to be triggered when trajectory impacts surface of secondary body
def event_impact_secondary(t, X, mu, R_sec):
    x = X[0]
    y = X[1]
    z = X[2]
    
    r_2 = np.sqrt(y**2 + z**2 + (x - 1 + mu)**2)
    return r_2 - R_sec

# Event function to be triggered when specified state = 0
def event_state_crossing(t, X, state_ind):
    return X[state_ind]


# Compute stability index
def compute_stability_index(state, period, mu):
    
    y_prop_stm = propagate(add_stm_i(state), mu, period, with_stm = 1)
    stm_f = y_prop_stm.y[6:,-1].reshape((6,6))
    evals,evecs = np.linalg.eig(stm_f)    
    stab_ind = (np.abs(evals[0]) + np.abs(1/evals[0]))/2
        
    return stab_ind
    