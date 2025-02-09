# LB-2D-thrombolysis
Lattice-Boltzmann simulation of thrombolysis in 2D blood circuits

ABSTRACT

This study explores the use of the Lattice-Boltzmann (LB) method to simulate the dissolution of blood clots in a vascular system, introducing a novel computational approach to model thrombolysis. Blood clots, while essential for healing, can lead to severe medical conditions when formed inappropriately, necessitating effective treatment strategies. This work presents an innovative 2D LB framework that models clot interactions with blood flow and tissue plasminogen activator (tPA) through a force field representation of clot resistance. Unlike conventional methods that explicitly define clot permeability, our approach enforces clot resistance through a velocity-proportional force, effectively linking it to Darcy’s law. Combining this force-based clot model with advection-diffusion dynamics for tPA transport and binding, this work lays the foundation for a new modeling approach for studying clot degradation.

CODES
- lb_2D_fluid_with_clot.py : Runs an LB simulation in 2D geometry with a clot and converges the fluid which is then saved in a sperate file.
- lb_2D_thrombolysis : Lysis of the clot by tPA. Must be run after lb_2D_fluid_with_clot.py to use an already converged fluid.
- functionsLB2.py : Contains all functions necessary to execute the simulation.
- functionsMonitoring2.py : Contains functions to monitor the progress and values of the simulation.
