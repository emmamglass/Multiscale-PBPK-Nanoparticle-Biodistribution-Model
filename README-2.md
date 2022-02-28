# Multiscale-PBPK-Nanoparticle-Biodistribution-Model
This repository contains all matlab code files relevant to the Physiologically Based Multiscale Pharmacokinetic Model for Determining the Temporal Biodistribution of Targeted Nanoparticles paper.

ModelA.m:
Matlab code for the ODE model with 17 linear ODEs that describes the biodistribution of NP in five organ compartments: lung, heart, kidney, liver, and spleen. 

ModelB.m:
Matlab code for the ODE model with 23 linear ODEs tht describes the biodistribution of NP in seven organ compartments: lung, heart, kidney, liver, spleen, gut, and other compartments. This model also has an updated architecture that is more physiologically relevant.

BranchedModel.m:
Matlab code for the ODE model with 457 linear ODEs that describe the biodistribution of NP in the same organ compartments as model b, but updated with an arterial branching architecture that is most physiologically relevant. 

modelA_julia.ipynb:
Julia code for model A

modelB_julia.ipynb:
Julia code for model B

branched_model_julia.ipynb:
Julia code for branched model

NeuralNet_QSSA_julia.ipynb:
Neural Net solver incorporating quasi-steady-state approximation

# Running the Julia Code
There are two ways to run the available Julia code; 1) Using VS code, 2) Using Jupyter Notebook.
Running the code using VS code
i) Download and install Julia for you machine: https://julialang.org/downloads/
ii) Open Julia terminal on your machine and run the following commands.
    using Pkg
    Pkg.add(["DifferentialEquations", "Flux", "Plots", "PyPlot", "IJulia"])
ii) Download VS code for your machine from: https://code.visualstudio.com/download
iii) Launch VS code and select "extensions" from the icons list on the left.
iv) Type Julia in the search bar and select "Julia Language Support" from the search     results and install the extension on VS code
v) Open any of the .jl files in VS code and click on "Julia: Execute active File in REPL" icon to run the code

Running the code using Jupyter Notebook
i) Download and install Julia for you machine: https://julialang.org/downloads/
ii) Open Julia terminal on your machine and run the following commands.
    using Pkg
    Pkg.add(["DifferentialEquations", "Flux", "Plots", "PyPlot", "IJulia"])
iii) Download Anaconda for your machine from: https://www.anaconda.com/products/individual
iv) Launch an instance of Jupyter Notebook on your machine: https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/execute.html
v) Open any of the .ipynb file inside Jupyter Notebook session and run the code by clicking on "Cell -> Run All"