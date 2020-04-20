function [LS] = LSeqSolve(LS,StateVar,VARIABLES,DOMAIN)

%% ======================================================================
% this function calculates the LS function based on different temporal 
% methods. The schemes implemented are RK1 (Euler step methods), RK2, RK3
% ONE thing I'm not quite sure about is whether I have to use velocities at
% the n+1 and n+1/2 in RK stages (to be checked later).
%% =======================================================================
% extrapolate the velocity
if isfield(LS, 'LS1')
    % extrapolate the velocity
    tempVARIABLES.LSgamma = VARIABLES.LSgamma;
    tempVARIABLES.LSbeta = VARIABLES.LSbeta;
    tempVARIABLES.Pe = VARIABLES.Pe;
    tempVARIABLES.dt = VARIABLES.dt;

    tempVARIABLES.dissolution = VARIABLES.dissolution.dissolution1;
    tempVARIABLES.Da = VARIABLES.Da.Da1;
    [LS.LS1] = LSVelocityExtrapolation(tempVARIABLES, StateVar, DOMAIN, ...
        LS.LS1);
    

    % save phi_o for getting the N_o the tube for reinitialization eq
    psi_o = LS.LS1.psi;
    % solve LS eq to get tilda psi
    LS.LS1.psi = solveHJEq(LS.LS1, VARIABLES, DOMAIN, psi_o);
    % reinitialize LS to a signed distance function
    LS.LS1.psi = LSreinitialization(VARIABLES, DOMAIN, LS.LS1, psi_o);
    % compute the normals
    [LS.LS1] = LSnormals(LS.LS1,DOMAIN);    
    
    % extrapolate the velocity
    tempVARIABLES.dissolution = VARIABLES.dissolution.dissolution2;
    tempVARIABLES.Da = VARIABLES.Da.Da2;
    [LS.LS2] = LSVelocityExtrapolation(tempVARIABLES, StateVar, DOMAIN, ...
        LS.LS2);
    % save phi_o for getting the N_o the tube for reinitialization eq
    psi_o = LS.LS2.psi;
    % solve LS eq to get tilda psi
    LS.LS2.psi = solveHJEq(LS.LS2, VARIABLES, DOMAIN, psi_o);
    % reinitialize LS to a signed distance function
    LS.LS2.psi = LSreinitialization(VARIABLES, DOMAIN, LS.LS2, psi_o);
    % compute the normals
    [LS.LS2] = LSnormals(LS.LS2,DOMAIN);    
    
    
else
    % extrapolate the velocity
    [LS] = LSVelocityExtrapolation(VARIABLES, StateVar, DOMAIN, LS);
    % save phi_o for getting the N_o the tube for reinitialization eq
    psi_o = LS.psi;
    % solve LS eq to get tilda psi
    LS.psi = solveHJEq(LS, VARIABLES, DOMAIN, psi_o);
    % reinitialize LS to a signed distance function
    LS.psi = LSreinitialization(VARIABLES, DOMAIN, LS, psi_o);
    % compute the normals
    [LS] = LSnormals(LS,DOMAIN);
end
% % save phi_o for getting the N_o the tube for reinitialization eq
% psi_o = LS.psi;
% % solve LS eq to get tilda psi
% LS.psi = solveHJEq(LS, VARIABLES, DOMAIN, psi_o);
% % reinitialize LS to a signed distance function
% LS.psi = LSreinitialization(VARIABLES, DOMAIN, LS, psi_o);
% % compute the normals
% [LS] = LSnormals(LS,DOMAIN);



