# dual-fishery-management-model
A collection of code resources used to complete my MSc. in Mathematics. The focus was on the dual (non)management of a single fishery. This is an extension and partial re-write of existing code. From a programatic perspective, this is not good coding, but it works for what was needed. The usage of globals in this context was not intended, but is a byproduct of the original code.

## Simulation_Model [.m :: MATLAB]

Contains files coded in MATLAB. Use this to run a simulation for a single (vanilla) or dual managed fishery situation.

* EWS_Community          :: Simulation driver.
* EWS_dynamics_vanilla   :: Single fishery dynamics.
* EWS_dynamics_ross      :: Dual fishery dynamics.
* results_filter         :: Simple analysis of output from EWS_Community.

## Analysis [.m :: Mathematica] (Invalid -- Do Not use -- Archive Only)

Mostly for backup, however, this contains Mathematica code, programmed to assist myself in the analysis of my simulation model above, absent from any stochastic effects. As it stands, trying to read this on the repo will be next to impossible. Recommended you download to use, if needed.

* detj               :: Produces eigenvalues of a 3 ODE system of equations from jacobian determinant.
* detj-pre-reduction :: Produces eigenvalues of a 3 ODE system of equations from jacobian determinant, with constraints pre-baked.
* GXSXLX             :: Graphing files, uses the vanilla parameters from the simulation model (with any extended parameters outlined in EWS_Community).
