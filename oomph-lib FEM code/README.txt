This folder contains the files regarding the continuum modelling.

Continuum description of the system is solved with oomph-lib.

Oomph-lib is a FEM package that you can obtain from
https://oomph-lib.github.io/oomph-lib/doc/html/index.html

For installation instructions, see
https://oomph-lib.github.io/oomph-lib/doc/the_distribution/html/index.html

The folder contains three subfolders. 

- Diffusivities_2D_Augustae contain the Maple script used to compute the diffusion and swimming
direction components in the oomph-lib model.
- user_driver contains the driver codes, the codes containing the actual model for the problem.
The folder also contains a Matlab script for postprocessing the data (postprocess_euler.m and drift_disp_euler.m) and
a code piece to compute velocity scaling factor (velsc_calculator.m). Move GYRO_OSC_Euler folder under
user_driver directory under oomph-lib to start.
- user_src folder contains the libraries used for running this model. Both the source files and 
the driver codes need to be installed into the oomph-lib distribution, see the directions in the 
second link in this README file