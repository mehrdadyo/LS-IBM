clc
clear all
% close all
% profile on

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

[StateVar,VARIABLES,DOMAIN,BC,IBM,LS, LSCase]= setUpVariablesNonDim;
%% ==============================
% load 'initial at any time '


%% object location and IBM coefficients
% [IBM_coeffU,IBM_coeffV,IBM_coeffP]=IBMcoeffs(IBM,DOMAIN);

Flux = getDiffFlux(VARIABLES, DOMAIN,BC);

ControlVar = setUpControlVar(VARIABLES, DOMAIN);


 [StateVar, LS, ControlVar, VARIABLES] = ...
    modelSimulation(StateVar, ControlVar, Flux, ...
    DOMAIN, VARIABLES, IBM,BC, LS);
    
% [StateVar.PSI,StateVar.VOR]=PlotFlowField(CurrentStateVar,IBM,DOMAIN);
% 
% profile viewer