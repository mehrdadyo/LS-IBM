function [StateVar, LS, ControlVar, VARIABLES] = ...
    modelSimulation(StateVar, ControlVar, Flux,...
        DOMAIN, VARIABLES, IBM,BC, LS)

% model 1 : march all 3 eqns US together in time
% model 2 : slove flow and transport untill SS and then march LS
% model 3 : SS flow and transport and LS unsteady
% model 4 : SS flow and solve transport to SS and then solve unsteady LS

model = ControlVar.model;

if model == 1
    message = 'Model 1 : Simulation starts with marching all 3 eqns US together in time:';
elseif model == 2
    message = 'Model 2 : Simulation starts with sloving flow and transport untill SS and then march LS:';
elseif model == 3
    message = 'Model 3 : Simulation starts with flow and transport being steady state and LS being unsteady:';
elseif model == 4
    message = 'Model 4 : Simulation starts with flow steady state and transport and LS being unsteady:';
end
fprintf('=========================================================================== \n');
fprintf(['       ', message,'           \n']);
fprintf('=========================================================================== \n');

% for iTime = 1:ControlVar.noLStime
ControlVar.endTime = ControlVar.noLStime * VARIABLES.dt;
iTime = 0;
while ControlVar.time < ControlVar.endTime         
    ControlVar.time = ControlVar.time+VARIABLES.dt;
    iTime = iTime + 1; %ceil(ControlVar.time/ VARIABLES.dt);
    fprintf('=========================================================================== \n');
    fprintf('                      updating geometry at time = %8.6f                      \n'...
        ,ControlVar.time);
    fprintf('=========================================================================== \n');
    
%% ========================================================================
%                      LEVEL SET EQUATION            
%% ========================================================================
%     LS_temp = LS;        

    [LS] = LSeqSolve(LS,StateVar,VARIABLES,DOMAIN);
    [IBM_coeffU,IBM_coeffV,IBM_coeffP] = ...
        LSIBMcoeffs(IBM,DOMAIN,LS, StateVar.phi);
    disp(IBM_coeffP.message)
    
    if model == 1   
        ControlVar.PISO = 1;
        ControlVar.flow_steady = 0;
        ControlVar.transport_steady = 0;


        ControlVar.noTime = VARIABLES.nLSupdate ; 
        VARIABLES.dt = VARIABLES.dt/VARIABLES.nLSupdate;
        
        [StateVar] = flowUStransportUS(StateVar, ControlVar, Flux,...
            DOMAIN, VARIABLES, IBM,IBM_coeffU, IBM_coeffV,IBM_coeffP,BC, LS);
        
        VARIABLES.dt = VARIABLES.dt * VARIABLES.nLSupdate;
        
    elseif model == 2
        
    elseif model == 3
        ControlVar.reduceTime = 0;
        ControlVar.PISO = 0;
        ControlVar.flow_steady = 1;
        ControlVar.transport_steady = 1;
%         StateVar_temp = StateVar;
        [StateVar, ControlVar] = flowSStransportSS(StateVar, ControlVar, Flux,...
            DOMAIN, VARIABLES, IBM,IBM_coeffU, IBM_coeffV,IBM_coeffP,BC, LS);
%         if ControlVar.reduceTime
%             StateVar = StateVar_temp;
%             LS = LS_temp;
%             ControlVar.time = ControlVar.time - VARIABLES.dt;
%             VARIABLES.dt = 0.1*VARIABLES.dt;
%         end
            
        
    elseif model == 4
        
        ControlVar.PISO = 0;
        ControlVar.flow_steady = 1;
        ControlVar.transport_steady = 0;
        
        ControlVar.noTime = VARIABLES.nLSupdate ; 
        VARIABLES.dt = VARIABLES.dt/VARIABLES.nLSupdate;
        
        [StateVar] = flowSStransportUS(StateVar, ControlVar, Flux,...
            DOMAIN, VARIABLES, IBM,IBM_coeffU, IBM_coeffV,IBM_coeffP,BC, LS);
        
        VARIABLES.dt = VARIABLES.dt * VARIABLES.nLSupdate;
        
        
    end
%% ========================================================================
%                       CurrentStateVar            
%  ========================================================================     
    if isfield(LS, 'LS1')
        CurrentStateVar =  struct('U',StateVar.U,'V',StateVar.V,...
        'P',StateVar.P,'phi',StateVar.phi, 'psi1', LS.LS1.psi, ...
        'psi2', LS.LS2.psi, 'time', ControlVar.time);
    else
        CurrentStateVar =  struct('U',StateVar.U,'V',StateVar.V,...
                'P',StateVar.P,'phi',StateVar.phi, 'psi', LS.psi, 'time', ...
                ControlVar.time);
    end
  
    
%% ========================================================================
%                      PLOTS and SAVE            
%% ========================================================================
    
    
    if ~mod(iTime, 10)
        flowfilename = strcat('dataRDE',num2str(iTime),'dt.mat');
        matfile = strcat('Output/', flowfilename);
 
        save(matfile,'-v7.3','-struct','CurrentStateVar');

    end
    

end
