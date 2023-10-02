Supplementing material for the paper

## "Dispersion in Shallow Moment Equations"
    Ullika Scholz, Department of Mathematics, RWTH Aachen University, scholz@acom.rwth-aachen.de
    Julia Kowalski, Institute for Applied Geophysics and Geothermal Energy, RWTH Aachen University, julia.kowalski@rwth-aachen.de
    Manuel Torrilhon, Department of Mathematics, RWTH Aachen University, mt@acom.rwth-aachen.de


#### Abstract:  
Shallow moment models are extensions of the hyperbolic shallow water equations. They admit variations in the vertical profile of the horizontal velocity. This paper introduces a nonhydrostatic pressure to this framework and shows the systematic derivation of dimensionally reduced dispersive equation systems which still hold information on the vertical profiles of the flow variables. The derivation from a set of balance laws is based on a splitting of the pressure followed by a same-degree polynomial expansion of the velocity and pressure fields in vertical direction. Dimensional reduction is done via Galerkin projections with weak enforcement of the boundary conditions at the bottom and at the free surface. The resulting equation systems of order zero and one are presented in linear and nonlinear form for Legendre basis functions and an analysis of dispersive properties is given. A numerical experiment shows convergence towards the resolved reference model in the linear stationary case and demonstrates the reconstruction of vertical profiles.

#### How to use (requires [Wolfram Mathematica](https://www.wolfram.com/mathematica/) or the [Free Wolfram Engine for Developers](https://www.wolfram.com/engine/))

The code is provided twice: In Mathematica notebook format and in Wolfram Script (.wls) format.
The .wls-format is human-readable and can be executed with the Free Wolfram Engine available for download at the Wolfram website.
In case a Mathematica license is available, the notebook format may be more convenient.

The folder "notebooks" contains the code in notebook format. To get started mark all cells (ctrl+a) and run (shift+enter).

The folder "scripts" contains the equivalent Wolfram scripts.
They can be executed via "wolfram -script /path/to/script.wls". 
Plots produced by the scripts will be stored in a separate folder next to the scripts folder.

The folder "engine" contains the core packages for the numerical simulations, which
are used by both, the notebooks and the scripts.

The code was written for the Mathematica 13.0 Kernel. 
For help with the execution of the code or any further questions
feel free to contact the authors.
