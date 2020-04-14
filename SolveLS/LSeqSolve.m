function [LS] = LSeqSolve(LS,StateVar,VARIABLES,DOMAIN)

%% ======================================================================
% this function calculates the LS function based on different temporal 
% methods. The schemes implemented are RK1 (Euler step methods), RK2, RK3
% ONE thing I'm not quite sure about is whether I have to use velocities at
% the n+1 and n+1/2 in RK stages (to be checked later).
%% =======================================================================
% extrapolate the velocity
[LS.u,LS.v] = LSVelocityExtrapolation(VARIABLES, StateVar, DOMAIN, LS);
% save phi_o for getting the N_o the tube for reinitialization eq
psi_o = LS.psi;
% solve LS eq to get tilda psi
LS.psi = solveHJEq(LS, VARIABLES, DOMAIN, psi_o);
% reinitialize LS to a signed distance function
LS.psi = LSreinitialization(VARIABLES, DOMAIN, LS, psi_o);
% compute the normals
[LS] = LSnormals(LS,DOMAIN);



