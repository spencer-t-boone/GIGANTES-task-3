# Package 1A: Generation of families of periodic orbits in the CR3BP
This folder contains the scripts for the Package 1A prototype written as part of the software prototype for Task 3 of the GIGANTES project. The objective of package 1A is the generation of families of periodic orbits in the circular restricted three-body problem (CR3BP). 

This package contains a number of functionalities including:
* Propagation in the circular restricted three-body problem (CR3BP)
* Single shooting method to correct an initial guess into a periodic orbit
* Two methods to continue families of periodic orbits:
  * Natural parameter continuation (NPC)
  * Pseudo-arclength continuation (PALC)
* A method to detect bifurcations from a family, and continue the new resulting family
* Interpolation method to generate orbits with specific parameters within a family of periodic orbits


## Installation and requirements

To work with Package 1A, one can simply clone the repository in the local machine:

```bash
git clone "https://github.com/spencer-t-boone/GIGANTES-task-3"
```

Package 1A is found in the *family_continuation_prototype* folder.

Only invited developers can contribute to the folder, and each should create a separate branch or fork, otherwise push requests will not be accepted on main branch modifications. This work is under [European Space Agency Public License (ESA-PL) Permissive (Type 3) - v2.4 licence](https://essr.esa.int/license/european-space-agency-public-license-v2-4-permissive-type-3). See [LICENSE](https://github.com/MacPau/FLYbyENCELADUS/blob/main/LICENSE.txt)) file for details.

This prototype requires the installation of Python, version 3.6 or newer is recommended.

To run the software, the user will require the following Python modules: NumPy, SciPy, copy, Matplotlib.

## Usage and test cases

The script *test_script_saturn_enceladus.py* contains a number of validation cases and additional test cases that can be used to test the functionalities of this package. The validation cases were selected to fully validate all the required functionalities in the MIDAS implementation of this prototype. The validation and test cases are listed at the top of the script. To select a specific script to run, modify the *test* variable accordingly:

```python
# Select test case to run
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
                     "lyapunov to butterfly"]


test = validation_case_options[0]
```

The validation cases involve the following (in the Saturn-Enceladus CR3BP unless otherwise specified):
* "L2 halo": Generate L2 halo orbit family using PALC 
* "L1 halo NPC": Generate L1 halo orbit family using NPC along *x* direction 
* "butterfly NPC period": Generate butterfly orbit family using NPC along period 
* "butterfly NPC mu": Continue butterfly orbit in Earth-Moon CR3BP ($\mu = 0.0121505856$) to Saturn-Enceladus system ($\mu = 1.9011497893 \times 10^{-7}$) using NPC along $\mu$
* "L2 lyapunov": Continue L2 Lyapunov orbit family, identify bifurcation with L2 halo orbits, and continue resulting L2 halo orbit family
* "L2 bifurcations": Continue L2 halo orbit family, identify period-doubling and period-tripling bifurcations, and continue families corresponding to first period-doubling bifurcation and first period-tripling bifurcation
* "L1 halo PALC earthmoon": Generate L1 halo orbit family in Earth-Moon CR3BP using PALC

Several additional test cases are included in the script to demonstrate specific functionalities, these include
* "L1 halo": Generate $L_1$ halo orbit family using PALC
* "butterfly": Generate butterfly orbit family using PALC
* "butterfly NPC": Generate butterfly orbit family using NPC along *z* parameter
* "L1 p3-halo": Generate $L_1$ period-3 halo orbit family using PALC
* "lyapunov to butterfly": Continue $L_2$ Lyapunov orbit family in Earth-Moon system, detect bifurcation with halo orbits, then continue $L_2$ halo orbit family, detect bifurcation with butterfly orbits (4th period-doubling), continue butterfly orbit family

Some of the validation and test cases are described in further detail to demonstrate the use of the scripts.

### Validation case 1: Generate $L_2$ halo orbit family using PALC
```python
# -----------------------------------------------------------------------------
# L2 halo orbit test case using PALC
# -----------------------------------------------------------------------------
elif test == "L2 halo":
```
Define an initial guess for a member of the orbit family (state and period):
```python
    X_i = np.array([ 1.00446305e+00, 0,  2.59982886e-05, 0, -3.53884567e-03, 0])
    t_f_guess = 3.08985431e+00
```

Define the free variables and constraints for the corrections procedure:
```python  
    free_vars = ["x", "z", "ydot", "t"]
    constraints = ["y", "xdot", "zdot"]
```

Define the initial step size for the continuation method:
```python
    step = -5e-4
```

Create event function to use as stop criterion for continuation procedure:
```python
    event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
```

Run PALC continuation function with desired inputs using the *continue_family_palc()* function from the *continuation.py* file:
```python
    orbit_family_states, orbit_family_periods, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                     constraints, step, N_orbits_max=2000, half_period = 1,
                                                                     event_stop = event_impact_enceladus)
```

If continuation returns a success, then plot resulting family using the *plot_family()* function from the *cr3bp_plotting.py* file:
```python
    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 5, frame = 'sec-centric', 
                    R_sec = R_enc_km, r_sec = r_enc)
```

This gives the following plot:

![val_l2_halo_palc](https://github.com/user-attachments/assets/327ab9f9-f774-4ab9-b0cb-fab6488972cd)

We can then interpolate an orbit with a specific *z* value at periapse, which corresponds to 30 km above the surface of Enceladus. First, since the family contains the states at apoapse, we need to propagate each of the states one half period to periapse:

```python
# Interpolate orbit with specific z_i value (30 km above surface of Enceladus at periapse)
orbit_peri_states = np.zeros(orbit_family_states.shape)
for i in range(orbit_family_states.shape[0]):
    y_prop = propagate(orbit_family_states[i,:], mu, orbit_family_periods[i]/2)
    orbit_peri_states[i,:] = y_prop.y[:6,-1]
```

Next, we can interpolate the specific orbit with the desired target *z* value using the *interpolate_orbit()* function from the *family_interpolation.py* file:

```python
target_var = ["z"]
free_vars = ["x", "ydot", "t"]
target_value = -(R_enc + 30/r_enc)
orbit_interp_state, orbit_interp_period, flag = interpolate_orbit(mu, orbit_peri_states, orbit_family_periods,
                                                                          target_var, target_value,
                                                                          free_vars, constraints, half_period = 1)
```

Verifying in the console that this successfully targets a periodic orbit with the desired value:

```python
orbit_interp_state
Out: 
array([ 1.00006458e+00,  1.45568510e-12, -1.18330537e-03,  9.28814100e-13,
        1.68517609e-02,  1.25125379e-11])

target_value
Out: -0.0011833053691275167

orbit_interp_state[2]
Out: -0.0011833053691275167
```


### Validation case 2: Generate $L_1$ halo orbit using NPC along *x* direction 
For the second validation case, the $L_1$ halo orbit family is continued using the natural parameter continuation method along the *x* direction. Most of the variables are defined similarly to Validation case 1, but a continuation variable must be defined as well. Options for the continuation variable include *x*, *y*, *z*, *xdot*, *ydot*, *zdot*, *t*, and *mu*. The continuation procedure is then run with the desired inputs using the *continue_family_npc()* function from the *continuation.py* file:

```python
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
        plot_family(orbit_family_states, orbit_family_periods, mu, spacing = 20, frame = 'sec-centric', 
                    R_sec = R_enc_km, r_sec = r_enc)
```

This gives the following plot of the $L_1$ halo orbit family:

![val_l1_halo_npc](https://github.com/user-attachments/assets/7b97b84a-1a2d-4c4b-8a80-7b8e24738230)


### Validation case 3: Generate butterfly orbit family using NPC along period 
For the third validation case, the butterfly orbit family is generated using the natural parameter continuation method and continuing along the period, starting from an initial guess with a period of $\tau = 3.2255$. The continuation procedure is run using the *continue_family_npc()* function from the *continuation.py* file:

```python
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
```

This gives the following plot of the butterfly orbit family:

![val_butterfly_npc_period](https://github.com/user-attachments/assets/0b36196b-fec9-4937-bf72-9c6be63e58b2)

### Validation case 4: Continue butterfly orbit in Earth-Moon CR3BP to Saturn-Enceladus CR3BP using NPC along $\mu$
For the fourth validation case, a member of the butterfly orbit family in the Earth-Moon system ($\mu = 0.0121505856$) is continued into the Saturn-Enceladus system ($\mu = 1.9011497893 \times 10^{-7}$) using the natural parameter continuation method, by varying the $\mu$ parameter. In order to do so, the value of $\mu$ is multiplied by a fixed value $\Delta s$ at each step. In this specific validation case, the number of steps is set to $n_{step}  = 500$. 

```python
# -----------------------------------------------------------------------------
# Compute Saturn-Enceladus butterfly orbit family by continuing along mu from 
# Earth-Moon butterfly orbit family
# ----------------------------------------------------------------------------- 
elif test == 'butterfly NPC mu':
    X_i = np.array([1.0669, 0, 0.1619, 0, -0.0189, 0]) # Initial conditions in Earth-moon system
    t_f_guess = 3.2255
    mu_i =  0.0121505856
```
Here the value of the step size is computed given the desired number of steps:
```python
    # Compute continuation step
    target_mu = mu
    num_step = 500
    step = (mu_i/target_mu)**(1/num_step)
    
    mu_vector = np.zeros(num_step)
    for i in range(num_step):
        mu_vector[i] = mu_i/step**i
```

The orbit can then be continued using the NPC method while varying the mass parameter until the desired system is reached:

```python

    continuation_var = ["mu"]
    free_vars = ["x", "z", "ydot", "t"]
    constraints = [ "y", "xdot", "zdot"]
    
    
    orbit_family_states, orbit_family_periods, flag = continue_family_npc(X_i, mu_i, t_f_guess, continuation_var, free_vars, 
                                                                     constraints, step, N_orbits_max = num_step, half_period = 1)

    if flag == 1:
        plot_family(orbit_family_states, orbit_family_periods, mu_vector, spacing = 50, variable_mu = 1)
```

This test gives the following plot:

![val_butterfly_npc_mu](https://github.com/user-attachments/assets/93fca609-f427-4803-8347-e3d512770e18)

This strategy is useful for computing the butterfly orbits in the Saturn-Enceladus CR3BP because they do not directly bifurcate from the halo orbit family in this sytem (while they do in the Earth-Moon system).

### Validation case 5: Continue $L_2$ Lyapunov orbit family, identify bifurcation with L2 halo orbits, and continue resulting $L_2$ halo orbit family

For the next validation case, the bifurcation detection methods are demonstrated. In this case, the $L_2$ Lyapunov orbit family is first generated using the PALC method, as in the other validation cases:
        
```python
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
```

This gives the following plot of the $L_2$ Lyapunov orbit family:

![val_l2_lyapunov](https://github.com/user-attachments/assets/aaa307be-82c4-4653-854d-6c4987466b7f)

Next, the Broucke diagram for the $L_2$ Lyapunov orbit family is generated using the *plot_broucke_diagram()* function from the *bifurcation_functions* file:

```python
        # Plot Broucke diagram
        fig_broucke = plot_broucke_diagram(orbit_family_states, orbit_family_periods, mu)
```

This generates the following Broucke diagram:

![broucke_lyapunov](https://github.com/user-attachments/assets/8826dcc7-0fd5-4586-b7b3-1c3fa1c2e2dd)

We can see that the family crosses the tangent bifurcation line twice. The first of these crossings corresponds to the bifurcation with the halo orbit family. Note that it is not necessary to generate the Broucke diagram to later compute the bifurcations, but it can help easily visualize the location and existence of any bifurcations within the particular subset of the family of interest.

Next, all tangent bifurcations are detected using the *detect_bifurcations_broucke()* function from the *bifurcation_functions* file:

```python
        # Detect bifurcations
        bif_types = ['tangent']
        bif_states, bif_periods, bif_cont_states, bif_cont_periods, bif_types = detect_bifurcations_broucke(orbit_family_states, orbit_family_periods, mu,
                                                                                                            free_vars, constraints, bif_types = bif_types)
```

We can verify in the console that the two tangent bifurcations are indeed detected:

```python
bif_states
Out[3]: 
array([[ 1.00446306e+00, -3.12950469e-20, -8.39949907e-20,
         2.24329841e-14, -3.53858139e-03,  4.31816369e-28],
       [ 1.00546538e+00, -3.12950469e-20,  3.12302556e-28,
         2.24329841e-14, -1.08772737e-02,  4.31816369e-28]])

bif_periods
Out[4]: array([3.08985815, 4.13121377])

bif_types
Out[5]: ['tangent', 'tangent']
```

Finally, the initial guess for the first member of the bifurcated family (in this case, the halo orbit family) is used to continue the halo orbit family using PALC:

```python
    # Continue L2 halo orbit family from resulting tangent bifurcation
    X_i = bif_cont_states[0]
    t_f_guess = bif_cont_periods[0]
    step = -5e-4
    event_impact_enceladus = lambda t, X: event_impact_secondary(t, X, mu, R_enc)
    
    orbit_family_states, orbit_family_periods, flag = continue_family_palc(X_i, mu, t_f_guess, free_vars, 
                                                                     constraints, step, N_orbits_max=2000, half_period = 1,
                                                                     event_stop = event_impact_enceladus)
```

### Validation case 6: $L_2$ halo bifurcations

For the next validation case, the $L_2$ halo orbit family is continued, and the period-doubling and period-tripling bifurcations are detected. The script then continues the families corresponding to first period-doubling bifurcation and first period-tripling bifurcation.

First, the $L_2$ halo orbit family is generated using the PALC method:

```python
# ---------------------------------------------------------------------------------------------------------------------
# Continue L2 halo orbit family, detect first period-doubling and period-tripling bifurcations, and continue 
# both resulting families
# ---------------------------------------------------------------------------------------------------------------------
elif test == 'L2 bifurcations':

    X_i = np.array([ 1.00398819,  0.        ,  0.00327934,  0.        , -0.00561346, 0.        ])
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
```

Next, the Broucke diagram is generated for a subset of the $L_2$ halo orbit family:

```python
        # Plot Broucke diagram
        fig_broucke = plot_broucke_diagram(orbit_family_states, orbit_family_periods, mu)
```

![broucke_halo](https://github.com/user-attachments/assets/8bdfed9c-c69e-4fed-b3ca-88b2eb7f727a)

One can see two crossings of the period-doubling bifurcation line and one crossing of the period-tripling bifurcation line (it can be difficult to discern, but the two crossings of the period-doubling line are quite close together). These bifurcations are detected using the *detect_bifurcations_broucke()* function from the *bifurcation_functions* file:

```python
        # Detect bifurcations
        print('Detecting bifurcations')
        bif_types = ['period_doubling','period_tripling']
        bif_states, bif_periods, bif_cont_states, bif_cont_periods, bif_types = detect_bifurcations_broucke(orbit_family_states, orbit_family_periods, mu,
                                                                                                            free_vars, constraints, bif_types = bif_types, bif_dir_step = 1e-6,
                                                                                                            half_period = 1)
```

Using the console, we can verify that the three bifurcations are indeed detected:

```python
bif_types
Out[7]: ['period_doubling', 'period_doubling', 'period_tripling']
```

The family corresponding to the first period-doubling bifurcation is then continued using the initial guess provided by the *detect_bifurcations_broucke()* function.

```python
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
```

This results in the following family:

![val_period_doubling](https://github.com/user-attachments/assets/cef9c786-5443-4503-a392-c2d1b56c5b2a)

Next, the family corresponding to the first period-tripling bifurcation is continued using the initial guess provided by the *detect_bifurcations_broucke()* function.     

```python
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
```

This gives the following family:

![val_period_tripling](https://github.com/user-attachments/assets/9eeff87a-1bdd-4ac7-b10d-e9f335b8c8af)

### Validation case 7: $L_1$ halo orbit family in Earth-Moon system
For the final validation case, the $L_1$ halo orbit family is generated in the Earth-Moon system, to demonstrate the applicability of the algorithms in the prototype to a different system. In this case, the $\mu$ parameter and radius of the moon must be specified for the new system.

```python
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
```

The following plot showing the $L_1$ halo orbit family in the Earth-Moon system is generated:

![val_l1_halo_em](https://github.com/user-attachments/assets/840a887c-8618-4552-a2fd-9d7c6ebc04e3)

### Additional test case 4: Generate $L_1$ period-3 halo orbit family using PALC

Additional test cases are provided to demonstrate the generation of specific periodic orbit families. Two of these (Additional test cases 4 and 5) are discussed in detail in the Readme file.

First, the $L_1$ period-3 halo orbit family can be generated using the PALC method. Note that in the Saturn-Enceladus system, this family does not directly bifurcate from the halo orbit family. However, it does bifurcate from the halo orbit family in the Earth-Moon system.

```python
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
```

This gives the following plot of select members of the $L_1$ period-3 halo orbit family:

![test_period3_halo](https://github.com/user-attachments/assets/bd468ae5-c49f-495f-88c0-59a5c7003485)

### Additional test case 5: Starting from $L_2$ Lyapunov orbit family in Earth-Moon system, detect bifurcation with halo orbit family and continue $L_2$ halo orbit family, detect bifurcation with butterfly orbits and continue butterfly orbit family
This additional test case demonstrates how the butterfly orbit family can be obtained starting from the $L_2$ Lyapunov orbit family using bifurcation theory. Note that the $L_2 Lyapunov orbit family begins at the $L_2$ Lagrange point.

First continue the $L_2$ Lyapunov orbit family from an initial guess using the PALC method:
```python
# ---------------------------------------------------------------------------------------------------------------------
# Continue L2 Lyapunov orbit family in Earth-Moon system, detect bifurcation with halo orbits, 
# continue L2 halo orbit family, detect bifurcation with butterfly orbits (4th period-doubling), continue butterfly orbit family
# ---------------------------------------------------------------------------------------------------------------------
elif test == 'lyapunov to butterfly':
    mu =  0.0121505856
    R_moon = 1737.4
    r_moon = 384400
    
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

```

Next, detect the tangent bifurcation with the $L_2$ halo orbit family:

```python
        # Plot Broucke diagram
        fig_broucke = plot_broucke_diagram(orbit_family_states, orbit_family_periods, mu)
        
        # Detect bifurcations
        bif_types = ['tangent']
        bif_states, bif_periods, bif_cont_states, bif_cont_periods, bif_types = detect_bifurcations_broucke(orbit_family_states, orbit_family_periods, mu,
                                                                                                            free_vars, constraints, bif_types = bif_types, bif_dir_step = 1e-5)
```

Continue the $L_2$ halo orbit family using the initial guess obtained from the bifurcation detection function:

```python        
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
```

Generate the Broucke diagram for the $L_2$ halo orbit family in the Earth-Moon system:

```python
            # Plot Broucke diagram
            fig_broucke = plot_broucke_diagram(orbit_family_states_1, orbit_family_periods_1, mu)
```

This gives the following Broucke diagram (zoomed in on the portion of interest):

![broucke_halo_em](https://github.com/user-attachments/assets/9d2055ce-63e3-4a89-a8f9-3bfee732d8b1)

Next, detect the period-doubling bifurcations in the $L_2$ halo orbit family:

```python
            # Detect bifurcations
            bif_types = ['period_doubling']
            bif_states_1, bif_periods_1, bif_cont_states_1, bif_cont_periods_1, bif_types = detect_bifurcations_broucke(orbit_family_states_1, orbit_family_periods_1, mu,
                                                                                                                free_vars, constraints, bif_types = bif_types, bif_dir_step = 1e-4)
```
                                                                                                                
In the Earth-Moon system, there are four period-doubling bifurcations from the $L_2$ halo orbit family. We can use the console to verify that these are all detected:                                                                                                              
```python
bif_types
Out[6]: ['period_doubling', 'period_doubling', 'period_doubling', 'period_doubling']
```              

The fourth period-doubling bifurcation in the Earth-Moon system corresponds to the butterfly orbit family. Using the initial guess obtained from the bifurcation detection algorithm, we can continue the butterfly orbit family:

```python
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
```

This gives the following plot of the Earth-Moon butterfly orbit family:

![test_butterfly_em](https://github.com/user-attachments/assets/9eeebd03-4225-448b-b9a7-623d977f10f9)


                            

