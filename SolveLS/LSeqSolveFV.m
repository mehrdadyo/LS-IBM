function [psi] = LSeqSolveFV(LS,psinp,dt,u,v,DOMAIN,scheme)

%% ======================================================================
% this function calculates the LS function based on different temporal 
% methods. The schemes implemented are RK1 (Euler step methods), RK2, RK3
% ONE thing I'm not quite sure about is whether I have to use velocities at
% the n+1 and n+1/2 in RK stages (to be checked later).
%% =======================================================================

psi_n = LS.psi;

imax = DOMAIN.imax;
jmax = DOMAIN.jmax;
dy = DOMAIN.dyv;
dx = DOMAIN.dxu;
dV = DOMAIN.dV_p;

dtdV = dt./dV;

equation = "LevelSetEqn";

% convective fluxes
[Fe1,Fw1,Fn2,Fs3] = LSconvFluxe(u,v,dx,dy,imax,jmax);


[psinp] = LSFindDerivative(psi_n,DOMAIN,equation,StateVar);

% psie = psinp.psie;
% psiw = psinp.psiw;
% psin = psinp.psin;
% psis = psinp.psis;


psi_n1 = psi_n + dtdV.*( Fw1.*psinp.psiw - Fe1.*psinp.psie + ...
    Fs3.*psinp.psis - Fn2.*psinp.psin );

if scheme == 'RK1'

    psi = psi_n1; 

elseif scheme == 'RK2' || scheme == 'Rk3'
    
    [psinp] = LSFindDerivative(psi_n1,DOMAIN,equation,StateVar);


    psi_n2 = ...
        psi_n1 + dtdV.*( Fw1.*psinp.psiw - Fe1.*psinp.psie + ...
        Fs3.*psinp.psis - Fn2.*psinp.psin);

    psi = 0.5 * (psi_n + psi_n2);

elseif scheme == 'RK3'

    psi_n12 = 0.75 * psi_n + 0.25 * psi_n2;
    
    [psinp] = LSFindDerivative(psi_n12,DOMAIN,equation,StateVar);

    psi_n32 = psi_n12 + dtdV.*( Fw1.*psinp.psiw - Fe1.*psinp.psie + ...
        Fs3.*psinp.psis - Fn2.*psinp.psin);
    
    psi = (psi_n + 2* psi_n32)/3;

end






