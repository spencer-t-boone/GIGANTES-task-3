# Package 1A: Generation of families of periodic orbits in the CR3BP
This folder contains the scripts for the Package 1A prototype written as part of the software prototype for Task 3 of the GIGANTES project. The objective of package 1A is the generation of families of periodic orbits in the circular restricted three-body problem (CR3BP). 

This package contains a number of functionalities including:
* Propagation in the circular restricted three-body problem (CR3BP)
* Single shooting method to correct an initial guess into a periodic orbit
* Two methods to continue families of periodic orbits:
** Natural parameter continuation
** Pseudo-arclength continuation
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

The script *test_script_saturn_enceladus.py* contains a number of validation cases and additional test cases that can be used to test the functionalities of this package.
