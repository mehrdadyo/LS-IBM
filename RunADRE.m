clc
clear all
close all
% profile on

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

[StateVar,VARIABLES,DOMAIN,BC,IBM,LS, LSCase]= setUpVariablesNonDim3;
%% ==============================
% load 'initial at any time '


%% object location and IBM coefficients
% [IBM_coeffU,IBM_coeffV,IBM_coeffP]=IBMcoeffs(IBM,DOMAIN);

Flux = getDiffFlux(VARIABLES, DOMAIN,BC);

ControlVar = setUpControlVar(VARIABLES, DOMAIN);


%% Level set at time t = dt

[IBM_coeffU,IBM_coeffV,IBM_coeffP] = LSIBMcoeffs(IBM,DOMAIN,LS, StateVar.phi);



for iTime = 1:5
        
    ControlVar.resi=1;
    ControlVar.ii=0;
    ControlVar.time = ControlVar.time+VARIABLES.dt;
    
%% ========================================================================
%                      LEVEL SET EQUATION            
%% ========================================================================

    if ~mod(iTime, VARIABLES.nLSupdate)

        [LS] = LSeqSolve(LS,StateVar,VARIABLES,DOMAIN);
        [IBM_coeffU,IBM_coeffV,IBM_coeffP] = LSIBMcoeffs(IBM,DOMAIN,LS, StateVar.phi);
        StateVar.phi = (~IBM_coeffP.flag_p).* StateVar.phi;
    end
    

%% =======================================================================
%                       SOLVE FOR FLOW
% ========================================================================

    [StateVar,ControlVar, Soln] = ... 
        SolveUVP (ControlVar,Flux,DOMAIN,VARIABLES,StateVar,IBM,IBM_coeffU,...
        IBM_coeffV,IBM_coeffP,BC);    
    
    StateVar.U_old = StateVar.U;
    StateVar.V_old = StateVar.V;
    StateVar.P_old = StateVar.P;

%% ========================================================================
%                       SCALAR TRANSPORT            
%  ========================================================================
    
    if ControlVar.flow_steady == 0
        [StateVar] = SolveTransportADRE(ControlVar,DOMAIN,VARIABLES,... 
            StateVar,IBM,IBM_coeffP,BC,Flux, LS);
        
        [StateVar] = CheckforNegative(ControlVar,DOMAIN,VARIABLES,...
            StateVar,IBM,BC,Flux, LS);
        
    end
    
   
%% ========================================================================
%                       CurrentStateVar            
%  ========================================================================     
    
    CurrentStateVar =  struct('U',StateVar.U,'V',StateVar.V,...
            'P',StateVar.P,'phi',StateVar.phi, 'psi', LS.psi);
    
  
    
%% ========================================================================
%                      PLOTS and SAVE            
%% ========================================================================
    
    
    if ~mod(iTime, 100)
        flowfilename = strcat('dataRDE',num2str(iTime),'dt.mat');
        matfile = strcat('Output/', flowfilename);
 
        save(matfile,'-v7.3','-struct','CurrentStateVar');




    end


end

% [StateVar.PSI,StateVar.VOR]=PlotFlowField(CurrentStateVar,IBM,DOMAIN);
% 
% 
% profile viewer