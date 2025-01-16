# This file contains a number of functions for calculating various quantities
# in the circular restricted three-body problem (CR3BP)
# Author: Spencer Boone

import numpy as np


# Compute state derivatives in CR3BP
def state_derivs_cr3bp(X, mu):

    x = X[0]
    y = X[1]
    z = X[2]
    xdot = X[3]
    ydot = X[4]
    zdot = X[5]
    
    rho1 = (y**2 + z**2 + (-1 + x + mu)**2)
    rho2 = (y**2 + z**2 + (x + mu)**2)
    
    xdotdot = 2*ydot + x - (mu*(-1 + x + mu))/rho1**(3/2) - (1-mu)*(x+mu)/rho2**(3/2)
    ydotdot = -2*xdot + y - y*mu/rho1**(3/2) - y*(1-mu)/rho2**(3/2)
    zdotdot = -z*mu/rho1**(3/2) - z*(1-mu)/rho2**(3/2)
    
    
    dX = np.array([xdot, ydot, zdot, xdotdot, ydotdot, zdotdot])
    return dX
    

# Compute A matrix for propagating STM (phi_dot = A * phi)
def A_cr3bp(X, mu):
    
    x = X[0]
    y = X[1]
    z = X[2]
    xdot = X[3]
    ydot = X[4]
    zdot = X[5]
    
    rho1 = (y**2 + z**2 + (-1 + x + mu)**2)
    rho2 = (y**2 + z**2 + (x + mu)**2)
    r1 = rho2**(1/2)
    r2 = rho1**(1/2)            
    
    Uxx = -(1-mu)/r1**3 + 3*(x+mu)**2*(1-mu)/r1**5 - mu/r2**3 + 3*(x+mu-1)**2*mu/r2**5 + 1
    Uxy = 3*(x+mu)*(1-mu)*y/r1**5 + 3*mu*y*(x+mu-1)/r2**5
    Uxz = 3*z*(x+mu)*(1-mu)/r1**5 + 3*mu*z*(x+mu-1)/r2**5
    Uyy = -(1-mu)/r1**3 + 3*y**2*(1-mu)/r1**5 - mu/r2**3 + 3*mu*y**2/r2**5 + 1
    Uyz = 3*z*y*(1-mu)/r1**5 + 3*y*z*mu/r2**5
    Uzz = -(1-mu)/r1**3 + 3*z**2*(1-mu)/r1**5 - mu/r2**3 + 3*mu*z**2/r2**5
    
    Umat = np.array([[Uxx, Uxy, Uxz],[Uxy, Uyy, Uyz],[Uxz, Uyz, Uzz]])
    Kmat = np.array([[0,2,0],[-2,0,0],[0,0,0]])
    
    A = np.vstack((np.hstack((np.zeros((3,3)),np.identity(3))),np.hstack((Umat,Kmat))))
    
    return A


# Compute Jacobi constant in CR3BP for given state
def jacobi_constant(X,mu):
    x = X[0];
    y = X[1];
    z = X[2];
    xdot = X[3];
    ydot = X[4];
    zdot = X[5];
    
    r2 = np.sqrt(y**2 + z**2 + (-1 + x + mu)**2)
    r1 = np.sqrt(y**2 + z**2 + (x + mu)**2)
    v = np.linalg.norm(X[3:6])
    
    U = 1/2*(x**2 + y**2) + (1 - mu)/r1 + mu/r2
    
    C = 2*U - v**2
    
    return C


# Compute Jacobi constants for all members in a family of orbits/states
def jacobi_constant_family(states, mu):
    num_orbits = states.shape[0]
    jc_vec = np.zeros(num_orbits)
    
    for i in range(states.shape[0]):
        jc_vec[i] = jacobi_constant(states[i,:], mu)
        
    return jc_vec[i]
        
    
 # Compute distances between state and primary/secondary bodies   
def rel_dist_cr3bp(X, mu):

    
    r_1 = ((X[0] + mu) ** 2 + X[1] ** 2 + X[2] ** 2) ** 0.5
    r_2 = ((X[0] - 1 + mu) ** 2 + X[1] ** 2 + X[2] ** 2) ** 0.5

    return r_1, r_2


# Add phi_0 (initial STM) as a vector to a state vector (for integrating STM)
def add_stm_i(X):
    nX = np.size(X)
    phi0 = np.eye(nX)
    return np.hstack([X, phi0.reshape(-1)])
