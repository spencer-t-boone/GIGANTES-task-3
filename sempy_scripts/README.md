# SEMPy scripts
# Packages 1B, 3 and 4
This folder contains the scripts that depend on SEMpy, ISAE-SUPAERO's open source library dedicated to astrodynamics calculations in non-Keplerian environments (see https://gitlab.isae-supaero.fr/sempy). 

These scripts include:
* *global_low_alt_orbits_search.py* (Package 1B): Global analysis of initial conditions for low-altitude elliptical orbits
* *stationkeeping_general.py* (Package 3): Preliminary navigation analysis
* *direct_insertion_tp_general.py* (Package 4): Connecting science orbits with end-game conditions

Each of these scripts is independent and has a number of test cases that can be selected within the script, corresponding to the different periodic science orbits identified in the CR3BP. These science orbits are:
* $L_1$ near rectilinear halo orbit (NRHO)
* $L_2$ NRHO
* Butterfly orbit
* $L_1$ period-3 halo orbit

Initial conditions for the science orbits were generated using Package 1A. Each science orbit corresponds to a periapse at South pole of around 30 km above the surface of Enceladus.


## Installation and requirements

To work with these Packages, one must first install the SEMpy environment. While a public version of SEMpy is available on the public GitLab repository, a newer (currently unreleased) version is required to run 
the scripts in this folder.

, one can simply clone the repository in the local machine:

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
