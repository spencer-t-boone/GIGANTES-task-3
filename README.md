# Science orbit analysis for an Enceladus orbiter
This folder contains several packages written as part of the software prototype for Task 3 of the GIGANTES project. The purpose of Task 3 is the development of tools to enable the construction of families of orbits at Enceladus (and any planetary moon), and the analysis of their stability and station-keeping cost. The work for this task was subdivided into separate packages to organize the different functionalities. These packages are:

* Package 1A: Generation of families of periodic orbits in the CR3BP
* Package 1B: Global analysis of initial conditions for low-altitude elliptical orbits
* Package 2: Refinement of periodic orbits in ephemeris model
* Package 3: Preliminary navigation analysis
* Package 4: Connecting science orbits with end-game conditions

These packages are found in the various subfolders in this codebase. These are organized according to their implementation in the following folders:

* family_continuation_prototype: This folder contains the scripts for Package 1A. These scripts are written in Python and are independent of any software package, requiring only the standard Python NumPy and SciPy libraries.
* sempy_scripts: This folder contains the scripts for Packages 1B, 3 and 4. These scripts require the installation of SEMPy, an open-source multi-body dynamics library developed at ISAE-SUPAERO. These scripts are written in Python.
* transition_to_ephemeris model: This folder contains the scripts for Package 2. These scripts require the installation of GODOT and MIDAS, ESA's flight dynamics and mission analysis libraries. 
