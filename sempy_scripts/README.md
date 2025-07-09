# SEMPy scripts (Packages 1B, 3 and 4)
This folder contains the scripts that depend on SEMpy, ISAE-SUPAERO's open source library dedicated to astrodynamics calculations in non-Keplerian environments (see https://gitlab.isae-supaero.fr/sempy). 

These scripts include:
* *global_low_alt_orbits_search.py* (Package 1B): Global analysis of initial conditions for low-altitude elliptical orbits
* *stationkeeping_general.py* (Package 3): Preliminary navigation analysis
* *direct_insertion_tp_general.py* (Package 4): Connecting science orbits with end-game conditions

Each of these scripts is independent. Packages 3 and 4 contain a number of test cases that can be selected within the script, corresponding to the different periodic science orbits identified in the CR3BP. These science orbits are:
* $L_1$ near rectilinear halo orbit (NRHO)
* $L_2$ NRHO
* Butterfly orbit
* $L_1$ period-3 halo orbit

Initial conditions for the science orbits were generated using Package 1A. Each science orbit corresponds to a periapse at South pole of around 30 km above the surface of Enceladus.


## Installation and requirements

To work with these Packages, one must first install the SEMpy environment. While a public version of SEMpy is available on the public GitLab repository, a newer (currently unreleased) version is required to run 
the scripts in this folder. To access this version, please contact Thibault Gateau (Thibault.GATEAU@isae-supaero.fr) to gain access to the private repository, and contact Thibault Gateau or Spencer Boone (spencer.boone@isae-supaero.fr) for installation instructions. Note that the branch should be switched to 'gigantes'.

Only invited developers can contribute to the folder, and each should create a separate branch or fork, otherwise push requests will not be accepted on main branch modifications. This work is under [European Space Agency Public License (ESA-PL) Permissive (Type 3) - v2.4 licence](https://essr.esa.int/license/european-space-agency-public-license-v2-4-permissive-type-3). See [LICENSE](https://github.com/MacPau/FLYbyENCELADUS/blob/main/LICENSE.txt)) file for details.

This prototype requires the installation of Python, version 3.6 or newer is recommended.

## Usage and test cases

### Package 1B: Global analysis of initial conditions for low-altitude elliptical orbits

Low-altitude elliptical orbits have been proposed as potential science orbits for an Enceladus orbiter. Following the strategies developed in the most recent work investigating this type of orbit [Ref. 1], an algorithm performing a global search over a variety of initial conditions was developed as part of the software prototype. In this algorithm, the initial orbital elements $r_p$ (radius of periapse), $r_a$ (radius of apoapse), $\Omega$ (right ascension of ascending node), and $\omega$ (argument of periapse) are varied, while the initial inclination $i$ is held fixed for each run of the algorithm. The resulting states are converted into the Saturn-Enceladus CR3BP frame, and propagated forward in the CR3BP until either the trajectory collides with the surface of Enceladus, escapes from the vicinity of Enceladus, or a specified maximum integration time $t_{max}$ is reached. Note that this strategy could also be employed with higher-fidelity dynamics to produce a more realistic analysis. Due to the effect of Saturn's gravity on the orbits, most orbits at high inclinations will impact the surface of Enceladus relatively quickly.

To run Package 1B, open the *global_low_alt_orbits_search.py* script. In the "INPUT PARAMETERS" section, set the desired ranges for each orbital parameter to be varied:

```python
# -----------------------------------------------------------------------------
# INPUT PARAMETERS
# -----------------------------------------------------------------------------
# Set up simulation parameters
```

First select the desired maximum integration time (time after which to stop the integration). At low inclinations ($i < 65^{\circ}), this is necessary to prevent the script from running too long, since stable orbits exist at these inclinations. At higher inclinations, long-term stable orbits do not exist, so generally no orbits will reach the maximum integration time if it is sufficiently long.

```python
# Max integration time (in non-dim)
tau_days = 150 # days
tau = tau_days/env.adim_t*3600*24
```

Next, select the ranges of orbital elements to vary. This script varies $r_p$, $r_a$, $\Omega$ and $\omega$. Select the (fixed) inclination to run the analysis at. In the MIDAS implementation, the inclination is added as an additional parameter that can be varied. 

```python
# Ranges of orbital elements
# Vary radius of periapse, radius of apoapse, inclination, or RAAN
r_p_range = np.linspace(20+252.1,100+252.1,9) # 20 km to 100 km
r_a_range = np.linspace(50+252.1,200+252.1,16) # 50 km to 200 km
i_range = 65.0 # Currently script runs for fixed inclination value
raan_range = np.linspace(0,360,9) # 0 to 360 degrees
arg_p_range = np.linspace(0,360,16) # 0 to 360 degrees
```

The output of the script will be an array containing all the different sets of initial orbital elements that were tested, and an array with a maximum integration time corresponding to each set of initial orbital elements. It is helpful to plot the results in 2 dimensions, so the maximum integration time for each $r_p$ - $r_a$ pair is saved, and plotted as output.

For an inclination of $90^{\circ}$ (fully polar initial orbit conditions), this produces the following plot:

![low_altitude_90deg](https://github.com/user-attachments/assets/04a06252-3316-487d-b198-642f3331540a)

As can be seen from this plot, the longest-lifetime orbits have a lifetime of around 9 days. 

When at lower inclinations, long-lifetime orbits emerge. For example, for $i = 55^{\circ}$, the following plot is obtained:

![low_altitude_55deg](https://github.com/user-attachments/assets/d9366f8c-f3cc-4540-ba4b-0bcec77d1f71)

At $i = 65^{\circ}$, the following plot is obtained:

![low_altitude_65deg](https://github.com/user-attachments/assets/37809ca1-fa7b-4e96-a161-761ec31f37ba)

At $i = 70^{\circ}$, the following plot is obtained, showing that no long-lifetime orbits exist at this inclination. This lines up with the findings from the literature [Refs. 1-3].

![low_altitude_70deg](https://github.com/user-attachments/assets/38b02e9d-a3ce-4c4b-bdb6-c7816b847eba)


### Package 3: Preliminary navigation analysis
In practice, any science orbit will not perfectly conform to its baseline orbit. Small uncertainties such as navigation errors, thruster performance errors, and errors in the dynamics model will cause the spacecraft to deviate from its nominal trajectory. If not accounted for, the spacecraft will potentially depart from the vicinity of Enceladus or impact its surface. Station-keeping maneuvers are necessary to maintain the overall geometry of the orbits.

As part of the software prototype, the preliminary station-keeping strategy developed in [Ref. 4] was adapted for all of the low-energy periodic orbits. In this strategy, a maneuver is performed at the crossing of a specific altitude (in this case, 600 km). This maneuver targets the position at the next apse location (either periapse or apoapse). Navigation errors are added at the locations of the maneuvers, and maneuver execution errors are considered. In the prototype, demonstrations of this station-keeping strategy are given for the different science orbits, for both nominal and high levels of navigation errors and maneuver execution errors, obtained from the analysis in [Ref. 4]. This station-keeping strategy was found to be sufficient to maintain stability of all orbit configurations.

Package 3 was ultimately not implemented in the final software delivery.

To run the script, open the *stationkeeping_general.py* script. Select the science orbit for which to run the script.

```python
# -----------------------------------------------------------------------------
# INPUT PARAMETERS
# -----------------------------------------------------------------------------
# Select orbit type to perform analysis
orbit_type_options = ['L2 NRHO', 'L1 NRHO', 'butterfly', 'period-3 halo']
orbit_type = orbit_type_options[0]
orbit_type = 'L2 halo'
```

Then select the desired error model:

```python
# Select error model
error_model_options = ['nominal', 'high', 'none']
error_model = error_model_options[1]
```

The script will then run the station-keeping algorithm for the number of revolutions specified in the scenario parameters. These can be modified in the following lines of code (only shown for $L_2$ NRHO):

```python
# Import selected orbit and set up scenario parameters
# L2 NRHO
if orbit_type == 'L2 NRHO':
    X_i = np.array([ 1.0025797154342282, 0,  0.004882278399970316, 0, -0.005445984590084007,  0])
    tau = 2.287781993612279
    
    num_ref_points = 2 
    noise_freq = 2
    N_rev = 100
```

The number of revolutions *N_rev* can be modified to run the simulation for the desired length.

For example, when run for the $L_2$ NRHO scenario with the 'high' error model, the script will generate the following plot:

![stationkeeping_results_l2nrho](https://github.com/user-attachments/assets/dbe5bed6-08fa-4923-989d-8d637670b8cb)

When run for the period-3 halo orbit scenario with the 'nominal' error model, the script will generate the following plot:

![stationkeeping_results_period3](https://github.com/user-attachments/assets/9e827079-85c0-400f-a64a-c99f6c15d6ef)



### Package 4: Connecting science orbits with end-game conditions
The various science orbits must be entered into after a series of flybys of Saturn’s many moons (also known as moon tour) is performed to lower the spacecraft’s Saturn-centric semi-major axis. This moon tour was the focus of the work completed as part of Task 2. At the end of the moon tour, one or more maneuvers are required to insert into the science orbits from a Saturn-centric approach trajectory. In order to fully characterize these orbits’ suitability for future missions, these insertion costs can be computed for each configuration in Package 4.

In Package 4, the direct-insertion strategy is implemented using a two-maneuver strategy. First, an intermediate maneuver is applied to transfer from the endgame conditions to the transfer trajectory (that will target the science orbit). Then, an insertion maneuver is performed to enter the science orbit. This procedure is shown below. 

![direct_insertion_diagram](https://github.com/user-attachments/assets/0fa9e720-e6c7-4f75-a523-059503134154)
 
A backwards strategy was employed to compute the estimated minimum $\Delta V$ values for these maneuvers. The insertion maneuver is applied at a point along the periodic orbit and the resulting state is propagated backwards in time until the negative $x$-axis crossing. At the $x$-axis crossing, the trajectories’ orbital parameters are plotted on a Tisserand-Poincaré graph. For each of these points, the cost to perform the intermediate maneuver to transfer from the end-game conditions (i.e. $V_{\infty}$ = 200 m/s) onto the associated trajectory can be computed.

To run the script, open the *direct_insertion_tp_general.py* script in the *sempy_scripts* folder. Select the science orbit for which to run the analysis:

```python
# -----------------------------------------------------------------------------
# INPUT PARAMETERS
# -----------------------------------------------------------------------------
# Select orbit type to perform analysis
orbit_type_options = ['L2 NRHO', 'L1 NRHO', 'butterfly', 'period-3 halo']
orbit_type = orbit_type_options[0]
```

Next, select the maneuver location. For the halo orbits, this can be set to apoapse or periapse. For the butterfly orbit, this can be set to the apoapse at either the $L_1$ or $L_2$ 'lobes' of the orbit. 

```python
# Selection location for maneuver 
# Options for nrho and period-3 halo maneuver locations: 'periapse', 'apoapse'
# Options for butterfly maneuver locations: 'apoapse L2', 'apoapose L1'
maneuver_location = 'apoapse'
```

When run, the script generates a Tisserand-Poincaré graph with the points corresponding to the various maneuver magnitudes. Two options for the intermediate maneuver are then computed, either to raise the periapse from the 200 m/s line on the Tisserand-Poincaré graph, or to lower the apoapse, and enter onto the transfer trajectory.

An example of how the intermediate maneuver options is computed is shown below:

![butterfly_tp_l2](https://github.com/user-attachments/assets/87cd61dc-aa45-40a3-9954-2940441ec907)

For the $L_2$ halo orbit scenario, which was used as the validation case for the software implemention, the following TP-graph is obtained when run with maneuvers applied at apoapse:

![l2_halo_tp_apoapse](https://github.com/user-attachments/assets/51dcd9a4-e637-4bdd-b996-d052fbb60593)

And the following TP-graph is obtained when run with the maneuvers applied at periapse:

![l2_halo_tp_periapse](https://github.com/user-attachments/assets/eaea9b48-2593-43e4-b1b4-b12a09411fd3)

The minimum total $\Delta V$ cost (intermediate maneuver + transfer maneuver) is stored in the *delta_v_tot_min* variable (in km/s). Using the console we can check this value:

```python
delta_v_tot_min
Out[4]: 0.14539735196369163
```

The value obtained for the $L_2$ halo insertion at apoapsis is thus found to be 145.4 m/s, while the value for the insertion at periapsis is 63.5 m/s.





## References

[Ref. 1] Parihar, T., Hussmann, H., Wickhusen, K., Oberst, J., Stark, A., Caritá, G., Benidikter, A., Matar, J., Filho, E. S. R., and Galas, R., “Numerical analysis of polar orbits for future Enceladus missions,” EGU General Assembly 2024, 2024.
[Ref. 2] Russell, R. P., and Lara, M., “On the design of an Enceladus science orbit,” Acta Astronautica, Vol. 65, No. 1, 2009, pp. 27–39.
[Ref. 3] Lara, M., Palacian, J., and Russell, R., “Mission design through averaging of perturbed Keplerian systems: the paradigm of an Enceladus orbiter,” Celestial Mechanics and Dynamical Astronomy, Vol. 108, 2010, pp. 1–22.
[Ref. 4] MacKenzie, S., Neveu, M., Davila, A., Lunine, J., Craft, K., Cable, M., Phillips-Lander, C., Hofgartner, J., Eigenbrode, J., Waite, J., Glein, C., Gold, R., Greenauer, P., Kirby, K., Bradburne, C., Kounaves, S., Malaska, M., Postberg, F., Patterson, G., and Spilker, L., “The Enceladus Orbilander Mission Concept: Balancing Return and Resources in the Search fo Life,” The Planetary Science Journal, Vol. 2, 2021, p. 77.
