# Multiscale-PBPK-Nanoparticle-Biodistribution-Model
This repository contains all matlab code files relevant to the Physiologically Based Multiscale Pharmacokinetic Model for Determining the Temporal Biodistribution of Targeted Nanoparticles paper.

ModelA.m:
Matlab code for the ODE model with 17 linear ODEs that describes the biodistribution of NP in five organ compartments: lung, heart, kidney, liver, and spleen. 

ModelB.m:
Matlab code for the ODE model with 23 linear ODEs tht describes the biodistribution of NP in seven organ compartments: lung, heart, kidney, liver, spleen, gut, and other compartments. This model also has an updated architecture that is more physiologically relevant.

BranchedModel.m:
Matlab code for the ODE model with 457 linear ODEs that describe the biodistribution of NP in the same organ compartments as model b, but updated with an arterial branching architecture that is most physiologically relevant. 
