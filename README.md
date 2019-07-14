# Systems Biology

## Problem Statement
Two bacteria species are competing for space, food, etc. in a environment. One species (the controller) produces antibiotic; the antibiotic production rate is determined by a state variable, u, and the concentration of antibiotic. Presence of antibiotic in the environment slows the growth of the second species (the wild-type) more than the controller, so the goal is study the "controllability" of the system.

## Directories
#### Receding_horizon_controller_original_model
Contains files to simulate the system under the influence of a RHC (receding horizon control) algorithm. The algorithm is a flavor of MPC (model predictive control). The goal is to determine under what conditions the wild-type species concentration can be driven to a chosen value.

#### Receding_horizon_controller_ratio_model
Similar to the previous directory except the characterization of the underlying (mathematical) model is different.

#### Square_wave_controller_ratio_model
No RHC algorithm is used; the input, u, is a square wave normalized between 0 and 1.

#### No_controller_original_model
No control algorithm is used at all. u is defined at the start of the simulation.

#### Bifurcation_analysis
Plots diagrams to analyze changes in the system's behavior under parameter changes.