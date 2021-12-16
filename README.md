# NeuronGrowth_IGAcollocation
IGA-collocation Implementation for 2D neuron growth using Phase-field model

## User Guide
This code is the implementation of the phase field model with isogeometric collocation to simulate different stages of neuron growth by
considering the effect of tubulin. The stages modeled include lamellipodia formation, initial neurite outgrowth, axon differentiation, and
dendrite formation considering the effect of intracellular transport of tubulin on neurite outgrowth. The gradient computation of $\phi$ is carried out using cubic B-splines to increase the smoothness of the solution.

### File structures
*Neuron_Growth_Code folder* contains the implementation of 2D Phase field model for neuron growth code using isogeometric collocation method. 
*IGA_collocation_algorithm* folder contains functions built and used in the model.
*case_Xneuron_X_paper* folders contain code and saved data necessary to reproduce simulation results shown in the paper.

### How to run
1. Have a valid installation of Matlab
2. Navigate to *Neuron_Growth_Code folder*
3. Run *main.m*

### 
The matlab code implemented banded matrix calculation to speed up multiplications. To see the full equations, please refer to our paper:
*Modeling Neuron Growth Using Isogeometric Collocation Based Phase Field Method*