function [LS] = LSeqSolve(LS,StateVar,VARIABLES,DOMAIN)

%% ======================================================================
% this function calculates the LS function based on different temporal 
% methods. The schemes implemented are RK1 (Euler step methods), RK2, RK3
% ONE thing I'm not quite sure about is whether I have to use velocities at
% the n+1 and n+1/2 in RK stages (to be checked later).
%% =======================================================================

[LS.u,LS.v] = LSVelocityExtrapolation(VARIABLES, StateVar, DOMAIN, LS);


LS.psi = solveHJEq(LS, VARIABLES, DOMAIN);

LS.psi = LSreinitialization(VARIABLES, DOMAIN, LS);
[LS] = LSnormals(LS,DOMAIN);



