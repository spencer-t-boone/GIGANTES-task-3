# Functions used for astrodynamics calculations
# %% Import modules
import numpy as np
from copy import deepcopy


# Convert state (Cartesian) to orbital elements
# Note: mu refers to the gravitational parameter of the central body, not 
# the non-dimensional CR3BP mu parameter
def X2oe(X, mu):
    r_vec = X[:3]
    r = np.linalg.norm(r_vec)
    v_vec = X[3:6]
    v = np.linalg.norm(v_vec)
    h_vec = np.cross(r_vec, v_vec)
    h = np.linalg.norm(h_vec)
    n_vec = np.cross(np.array([0, 0, 1]), h_vec)
    n = np.linalg.norm(n_vec)
    e_vec = ((v**2 - mu/r)*r_vec - np.dot(r_vec, v_vec)*v_vec)/mu
    e = np.linalg.norm(e_vec)
    
    spec_en = v**2/2 - mu/r
    a = -mu/(2*spec_en)
    
    i = np.arccos(h_vec[2]/h)
    RAAN = np.arccos(n_vec[0]/n)
    if n_vec[1] < 0:
        RAAN = 2*np.pi - RAAN
        
    omega = np.arccos(np.dot(n_vec,e_vec)/(n*e))
    if e_vec[2] < 0:
        omega = 2*np.pi - omega
        
    nu = np.arccos(np.dot(e_vec,r_vec)/(e*r))
    if np.dot(r_vec, v_vec) < 0:
        nu = 2*np.pi - nu
        
    return np.array([a, e, i, RAAN, omega, nu])


# Convert orbital elements to Cartesian state vector
# Note: mu refers to the gravitational parameter of the central body, not 
# the non-dimensional CR3BP mu parameter
def oe2X(oe, mu):
    a = oe[0]
    e = oe[1]
    i = oe[2]
    RAAN = oe[3]
    omega = oe[4]
    nu = oe[5]
    
    p = a*(1-e**2)
    
    r_pqw = np.array([p*np.cos(nu)/(1 + e*np.cos(nu)),
                      p*np.sin(nu)/(1 + e*np.cos(nu)),
                      0.0])
    
    v_pqw = np.array([-np.sqrt(mu/p)*np.sin(nu),
                      np.sqrt(mu/p)*(e + np.cos(nu)),
                      0.0])
    
    rot_ijk_pqw = np.array([[np.cos(RAAN)*np.cos(omega) - np.sin(RAAN)*np.sin(omega)*np.cos(i),    -np.cos(RAAN)*np.sin(omega) - np.sin(RAAN)*np.cos(omega)*np.cos(i),    np.sin(RAAN)*np.sin(i)],
                            [np.sin(RAAN)*np.cos(omega) + np.cos(RAAN)*np.sin(omega)*np.cos(i),    -np.sin(RAAN)*np.sin(omega) + np.cos(RAAN)*np.cos(omega)*np.cos(i),    -np.cos(RAAN)*np.sin(i)],
                            [np.sin(omega)*np.sin(i),    np.cos(omega)*np.sin(i),    np.cos(i)]])
    
    
    r_ijk = rot_ijk_pqw@r_pqw
    v_ijk = rot_ijk_pqw@v_pqw
    
    return np.hstack([r_ijk, v_ijk])


# Function for plotting 3D sphere in plotly
def ms(x, y, z, radius, resolution=20):
    """Return the coordinates for plotting a sphere centered at (x,y,z)"""
    u, v = np.mgrid[0:2*np.pi:resolution*2j, 0:np.pi:resolution*1j]
    X = radius * np.cos(u)*np.sin(v) + x
    Y = radius * np.sin(u)*np.sin(v) + y
    Z = radius * np.cos(v) + z
    return (X, Y, Z)


# Event function for crossing of x axis
def event_poincare_y_axis(t, X):

    return X[1]


# Event function to be triggered when trajectory impacts surface of secondary body
def event_impact_secondary(t, X, mu, R_sec):
    x = X[0]
    y = X[1]
    z = X[2]
    
    r_2 = np.sqrt(y**2 + z**2 + (x - 1 + mu)**2)
    return r_2 - R_sec


# Convert state in CR3BP frame to Saturn-centric frame (in km or nd units)
def convert_cr3bp_frame_to_saturn_centric(X_cr3bp, env, output_units = 'km'):
    X_inertial = np.array([X_cr3bp[0], X_cr3bp[1], X_cr3bp[2], 
                           X_cr3bp[3] - X_cr3bp[1], X_cr3bp[4] + X_cr3bp[0], X_cr3bp[5]])
    
    X_saturn_centric = deepcopy(X_inertial)
    X_saturn_centric[0] += env.mu
    X_saturn_centric[4] += env.mu
    
    
    if output_units == 'km':
        X_saturn_centric_km = deepcopy(X_saturn_centric)
        X_saturn_centric_km[0:3] *= env.adim_l
        X_saturn_centric_km[3:6] *= env.adim_v
        
        return X_saturn_centric_km
    else:
        return X_saturn_centric


