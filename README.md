# Systems Biology

## Problem Statement
Two bacteria species are competing for space, food, etc. in a environment. One species (the controller) produces antibiotic; the antibiotic production rate is determined by a state variable, u, and the concentration of antibiotic. Presence of antibiotic in the environment slows the growth of the second species (the wild-type) more than the controller.

## Directories
#### Bifurcation Analysis
(Deprecated) Generates a plot demonsrating for what combinations of paramter values the wild-type species' concentration remains bounded in time. Needs heavy improvement.

#### Models
Mathematical models for the wild-type, controller, antibiotic, and metabolite growth rates. Currently there are three models, each more advanced at slowing the growth rate of the controller:
1) Degradation term.
2) Degradation term and carrying capacity.
3) Degradation term, carrying capacity, and metabolite.

#### RHC
Runs a receding horizon control (RHC) algorithm for a model from Models. Various parameters for the algorithm can be adjusted from inside the script.