# This file contains a number of functions to detect bifurcations within a family
# Author: Spencer Boone

import numpy as np
from cr3bp_propagate import propagate
from cr3bp_functions import add_stm_i
from copy import deepcopy
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from single_shooter import correct_ics
import scipy.linalg as sp
from single_shooter import map_vars_to_index


# Use Broucke stability diagram to detect bifurcations in a family of periodic orbits
def detect_bifurcations_broucke(orbit_family_states, orbit_family_periods, mu, free_vars, constraint_vars, 
                                bif_types = ["tangent", "period_doubling"], bif_dir_step = 1e-5, half_period = 1):
    
    
    # Map free variables to indices
    free_vars_index = map_vars_to_index(free_vars)


    # Generate Broucke parameters for each state vector in states
    states_size = orbit_family_states.shape[0] 
    alpha_vec = np.zeros(states_size)
    beta_vec = np.zeros(states_size)
    
    # Extract initial state and period
    print('Computing bifurcation parameters')
    for i in range(states_size):
        
        X_0 = orbit_family_states[i,:]
        t_f = orbit_family_periods[i]
        
        alpha, beta = get_bifurcation_params(X_0, t_f, mu)
        
        alpha_vec[i] = deepcopy(alpha)
        beta_vec[i] = deepcopy(beta)
        
    # Create interpolation of states and periods
    ind_vec = np.linspace(0,states_size-1,states_size)
    states_interp_fun = CubicSpline(ind_vec, orbit_family_states)
    period_interp_fun = CubicSpline(ind_vec, orbit_family_periods)
    
    
    
    bif_points = np.array([])
    bif_types_detected = []
    
    # Detect tangent bifurcations    
    if "tangent" in bif_types:
        tangent_bif_data = 2 + beta_vec + 2*alpha_vec
        tangent_bif_fun = CubicSpline(ind_vec, tangent_bif_data)
        tangent_bif_points = tangent_bif_fun.roots()
        bif_points = np.hstack([bif_points, tangent_bif_points])
        for i in range(tangent_bif_points.size):
            bif_types_detected.append("tangent")
            
        
    # Detect period-doubling bifurcations   
    if "period_doubling" in bif_types:
        pd_bif_data = 2 + beta_vec - 2*alpha_vec
        pd_bif_fun = CubicSpline(ind_vec, pd_bif_data)
        pd_bif_points = pd_bif_fun.roots()
        bif_points = np.hstack([bif_points, pd_bif_points])
        for i in range(pd_bif_points.size):
            bif_types_detected.append("period_doubling")
            
    # Detect period-tripling bifurcations   
    if "period_tripling" in bif_types:
        pd_bif_data = -1 + beta_vec - alpha_vec
        pd_bif_fun = CubicSpline(ind_vec, pd_bif_data)
        pd_bif_points = pd_bif_fun.roots()
        bif_points = np.hstack([bif_points, pd_bif_points])
        for i in range(pd_bif_points.size):
            bif_types_detected.append("period_tripling")
        
    if bif_points.size == 0:
        print('No bifurcations detected')
        return [],[]
        
    num_bifs = bif_points.size
    print(str(num_bifs), ' bifurcation points detected')
    
    bif_states_list = []
    bif_periods_list = []
    bif_states_guess_list = []
    bif_periods_guess_list = []
    bif_types_computed = []
    for i in range(num_bifs):
        if bif_points[i] < 0 or bif_points[i] > float(states_size-1):
            print("Skipping bifurcation at ", str(bif_points[i]))
        else:
            print("Interpolating bifurcation at ", str(bif_points[i]))
            bif_types_computed.append(bif_types_detected[i])
            X_i_guess = states_interp_fun(bif_points[i])
            period_guess = period_interp_fun(bif_points[i])
            
            # Correct interpolated state
            # Use half period to enforce constraints if orbit is symmetric along x-z plane
            if half_period == 1:
                period_guess = period_guess/2
                
            if bif_types_detected[i] == "period_doubling":
                period_guess *= 2
            elif bif_types_detected[i] == "period_tripling":
                period_guess *= 3
            
            # Correct initial conditions with fixed target parameter
            X_initial_corrected, period_corrected, flag, DF = correct_ics(X_i_guess, mu, period_guess, free_vars, constraint_vars, output_DF = 1)
            
            if half_period == 1:
                period_corrected *= 2
            
            if flag == 1:
                bif_states_list.append(X_initial_corrected)
                bif_periods_list.append(period_corrected)
                
                # Compute bifurcation direction - corresponds to second smallest
                # singular value of SVD, see Zimovan dissertation page 78
                sigma, V = np.linalg.svd(DF)[1:]
                bif_dir = V[-2,:]/np.linalg.norm(V[-2,:])
                
                bif_dir_mapped = np.zeros(7)
                bif_dir_mapped[free_vars_index] = bif_dir
                
                
                X_next = X_initial_corrected + bif_dir_mapped[:6]*bif_dir_step
                period_next = period_corrected + bif_dir_mapped[6]*bif_dir_step
                
                bif_states_guess_list.append(X_next)
                bif_periods_guess_list.append(period_next)
    
    
    return np.array(bif_states_list), np.array(bif_periods_list), np.array(bif_states_guess_list), np.array(bif_periods_guess_list), bif_types_computed



# Plot Broucke diagram
def plot_broucke_diagram(orbit_family_states, orbit_family_periods, mu):
    
    # Generate Broucke parameters for each state vector in states
    states_size = orbit_family_states.shape[0]  
    alpha_vec = np.zeros(states_size)
    beta_vec = np.zeros(states_size)
    
    # Extract initial state and period
    for i in range(states_size):
        X_0 = orbit_family_states[i,:]
        t_f = orbit_family_periods[i]
        
        alpha, beta = get_bifurcation_params(X_0, t_f, mu)
        
        alpha_vec[i] = deepcopy(alpha)
        beta_vec[i] = deepcopy(beta)

        
    # Plot Broucke diagram
    fig_broucke = plt.figure()
    plt.plot(alpha_vec,beta_vec, label = 'Family')
    print('alpha_vec: ', alpha_vec)
    print('beta_vec: ', beta_vec)
    
    # Plot bifurcation lines
    # Plot tangent bifurcation line
    alphas = np.linspace(-np.max(np.abs(alpha_vec)),np.max(np.abs(alpha_vec)),100)
    betas_tangent = -2 - 2*alphas
    plt.plot(alphas,betas_tangent,'--',label = 'Tangent')
    
    # Plot period doubling bifuraction line
    betas_period_doubling = -2 + 2*alphas
    plt.plot(alphas,betas_period_doubling,'--',label = 'Period doubling')
    
    # Plot period tripling bifurcation line
    betas_period_tripling = alphas + 1
    plt.plot(alphas,betas_period_tripling,'--',label = 'Period tripling')
    
    # Plot period quadrupling bifurcation line
    betas_period_quadrupling = (np.ones(alphas.size))*2
    plt.plot(alphas,betas_period_quadrupling,'--',label = 'Period quadrupling')
    
    plt.legend()
    
    plt.xlabel(r'$\alpha$')
    plt.ylabel(r'$\beta$')
    
    return fig_broucke
    
    
    
def get_bifurcation_params(X, t, mu):
    
    # Propagate state and STM for one orbit
    y_prop_stm = propagate(add_stm_i(X), mu, t, with_stm = 1)
    
    # Extract state and STM
    X_f = y_prop_stm.y[:6,-1]
    stm = y_prop_stm.y[6:,-1].reshape((6,6))
    
    # Compute Broucke parameters
    alpha = 2 - np.trace(stm)
    beta = 0.5*(alpha**2 + 2 - np.trace(stm @ stm))
    
    return alpha, beta
        
