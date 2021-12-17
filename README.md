# NeuronGrowth_IGAcollocation
2D Neuron growth using Phase-field model implemented using isogeometric collocation method (IGA-collocation).

## User Guide
This code is the implementation of the phase field model using isogeometric collocation to simulate multiple stages of neuron growth by considering the effect of tubulin. The stages modeled include lamellipodia formation, initial neurite outgrowth, axon differentiation, and dendrite formation considering the effect of intracellular transport of tubulin on neurite outgrowth. The gradient computation of Φ is carried out using cubic B-splines to increase the smoothness of the solution.

## File structures
- *[Neuron_Growth_Code](https://github.com/CMU-CBML/NeuronGrowth_IGAcollocation/tree/main/Neuron_Growth_Code)*: contains the implementation of 2D Phase field model for neuron growth code using isogeometric collocation method.
- *[Simulation_Cases_in_paper](https://github.com/CMU-CBML/NeuronGrowth_IGAcollocation/tree/main/Simulation_Cases_in_paper)*: contains simulation cases used in paper (see below). Each case can be reproduced by using the same random seed saved in the folder. (Note that these simulations were restarted, so random seed needs to be set based on log file)
	- *[IGA_collocation_algorithm](https://github.com/CMU-CBML/NeuronGrowth_IGAcollocation/tree/main/IGA_collocation_algorithm)*: contains functions built and used in the model.
	- *[case_Xneuron_X_paper]*: contains code and saved data necessary to reproduce simulation results shown in the paper.

## How to run
1. Valid installation of Matlab (code written with Matlab 2021a)
2. Navigate to *[Neuron_Growth_Code](https://github.com/CMU-CBML/NeuronGrowth_IGAcollocation/tree/main/Neuron_Growth_Code)*
3. Run *main.m*. For simulation cases in paper, run *main.m* in each folder for that specific case.

## 
The matlab code implemented banded matrix calculation to speed up multiplications. To see the full equations, please refer to our paper:
*Modeling Neuron Growth Using Isogeometric Collocation Based Phase Field Method* (in preparation)