The file descriptions:

Folder Passive Particle Model:
	- lee_fig6_data.m contains the data from Lee et al.'s (2014) Fig. 6 used in validation.
	- PAS_OSC_D.m is a spatial-diffusion based model of individual particles.

- GYRO_NO_OSC_sim.m simulates gyrotactic particles in non-oscillating (Poiseuille) flow.
- GYRO_OSC_sim.m is the simulation code for gyrotactic particles in oscillatory flows.
- postprocess_lagrangian.m evaluates key parameters such as the drift, dispersion etc.
- PP_GYRO_OSC_AUGUSTAE.mat contains the key parameters for C. augustae evaluated with postprocess_lagrangian.m.
- PP_GYRO_OSC_SALINA.mat contains the key parameters for D. salina evaluated with postprocess_lagrangian.m.
- PP_GYRO_OSC_STRONG.mat contains the key parameters for strongly gyrotactic particles (see text)
 evaluated with postprocess_lagrangian.m.
 
 Note that GYRO_OSC_sim.m records large files containing position information of each particle simulated,
 thus they are not included here. Those are required to run postprocess_lagrangian.m.