function [StateVar, ControlVar] = flowSStransportSS(StateVar, ControlVar, Flux,...
    DOMAIN, VARIABLES, IBM,IBM_coeffU, IBM_coeffV,IBM_coeffP,BC, LS)

% model 3 :SS flow and transport and LS unsteady

        
%% =======================================================================
%                       SOLVE FOR FLOW
% ========================================================================
    
[StateVar,ControlVar, Soln] = ... 
    SolveUVP (ControlVar,Flux,DOMAIN,VARIABLES,StateVar,IBM,IBM_coeffU,...
    IBM_coeffV,IBM_coeffP,BC);  
% if ControlVar.reduceTime
%     return;
% end
disp(ControlVar.messageFlow)

%% ========================================================================
%                       SCALAR TRANSPORT            
%  ========================================================================

[StateVar] = SolveTransportADRE(ControlVar,DOMAIN,VARIABLES,... 
    StateVar,IBM,IBM_coeffP,BC,Flux, LS);



