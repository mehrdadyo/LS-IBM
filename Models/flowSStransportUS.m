function [StateVar] = flowSStransportUS(StateVar, ControlVar, Flux,...
    DOMAIN, VARIABLES, IBM,IBM_coeffU, IBM_coeffV,IBM_coeffP,BC, LS)

% model 4 : SS flow and solve transport to SS and then solve unsteady LS

%% =======================================================================
%                       SOLVE FOR FLOW
% ========================================================================
[StateVar,ControlVar, Soln] = ... 
    SolveUVP (ControlVar,Flux,DOMAIN,VARIABLES,StateVar,IBM,IBM_coeffU,...
    IBM_coeffV,IBM_coeffP,BC);  
disp(ControlVar.messageFlow)
time = ControlVar.time;        
for iTime = 1:ControlVar.noTime
    time = time + VARIABLES.dt;
    fprintf('\n ~~~~~~~~~~~~~~~~~~~~ time = %8.6f ~~~~~~~~~~~~~~~~~~~~~ \n'...
        ,time);
    
%% ========================================================================
%                       SCALAR TRANSPORT            
%  ========================================================================

    [StateVar] = SolveTransportADRE(ControlVar,DOMAIN,VARIABLES,... 
        StateVar,IBM,IBM_coeffP,BC,Flux, LS);


%% ========================================================================
%                       BREAK if ERROR            
%  ========================================================================
    if sum(sum(isnan(StateVar.U))) || max(max(StateVar.phi))>1e4
        break
    end

end

