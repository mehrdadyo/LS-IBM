function [LS] = LSeqSolve(LS,StateVar,VARIABLES,DOMAIN)

%% ======================================================================
% this function calculates the LS function based on different temporal 
% methods. The schemes implemented are RK1 (Euler step methods), RK2, RK3
% ONE thing I'm not quite sure about is whether I have to use velocities at
% the n+1 and n+1/2 in RK stages (to be checked later).
%% =======================================================================
phi = StateVar.phi;
psi_n = LS.psi;
Da = VARIABLES.Da;
Pe = VARIABLES.Pe;
dt = VARIABLES.dt;
n_iter = VARIABLES.n_iter_ReLS;
dtau = VARIABLES.dtau;

[u,v] = LSVelocityExtrapolation(phi,Da,Pe,DOMAIN,LS);

LS.u = u;
LS.v = v;
equation = "LevelSetEqn";
scheme = VARIABLES.TimeSchemeLS;
% convective fluxes
[psinp] = LSFindDerivative(u,v,psi_n,DOMAIN,equation);

% psie = psinp.psie;
% psiw = psinp.psiw;
% psin = psinp.psin;
% psis = psinp.psis;
epsSign = min(min(DOMAIN.dxp))*2;

psi_n1 = psi_n - dt*( u.*psinp.psi_x + v.*psinp.psi_y );

if scheme == "RK1"

    psi = psi_n1; 

elseif scheme == "RK2" 
    
    [psinp] = LSFindDerivative(u,v,psi_n1,DOMAIN,equation);



    psi_n2 = psi_n1 - dt*( u.*psinp.psi_x + v.*psinp.psi_y );

    psi = 0.5 * (psi_n + psi_n2);
    

elseif scheme == "RK3"
    
    [psinp] = LSFindDerivative(u,v,psi_n1,DOMAIN,equation);
    
    psi_n2 = psi_n1 - dt*( u.*psinp.psi_x + v.*psinp.psi_y );


    psi_n12 = 0.75 * psi_n + 0.25 * psi_n2;
    
    [psinp] = LSFindDerivative(u,v,psi_n1,DOMAIN,equation);
    
    psi_n32 = psi_n12 - dt * ( u.*psinp.psi_x + v.*psinp.psi_y );
    
    psi = (psi_n + 2* psi_n32)/3;

end

psi_n = psi;

%% Reinitialize the distance function
equation = "ReinitializationEqn";

scheme = VARIABLES.TimeSchemeRLS;

for i=1:n_iter
      
    [psinp] = LSFindDerivative(u,v,psi,DOMAIN,equation);
    
    [G]= LSreInitilization(psinp,psi_n);
    fl = double(G ~= 0);
    SignPsi = psi_n./sqrt(psi_n.^2+epsSign.^2);
%     psi_m1 = psi - dtau*sign(psi_n).*(G-1);
    psi_m1 = psi - dtau*SignPsi.*(G-1).*fl;

    
    if scheme == "RK1"

        psi = psi_m1; 

    elseif scheme == "RK2" 

        [psinp] = LSFindDerivative(u,v,psi_m1,DOMAIN,equation);
        [G]= LSreInitilization(psinp,psi_n);
        fl = double(G ~= 0);
%         psi_m2 = psi_m1 - dtau*sign(psi_n).*(G-1);
        psi_m2 = psi_m1 - dtau*SignPsi.*(G-1).*fl;

        psi = 0.5 * (psi + psi_m2);


    elseif scheme == "RK3"

        [psinp] = LSFindDerivative(u,v,psi_m1,DOMAIN,equation);
        [G]= LSreInitilization(psinp,psi_n);
        fl = double(G ~= 0);

%         psi_m2 = psi_m1 - dtau*sign(psi_n).*(G-1);
        psi_m2 = psi_m1 - dtau*SignPsi.*(G-1).*fl;

        psi_m12 = 0.75 * psi + 0.25 * psi_m2;

        [psinp] = LSFindDerivative(u,v,psi_m12,DOMAIN,equation);
        [G]= LSreInitilization(psinp,psi_n);
        fl = double(G ~= 0);

%         psi_m32 = psi_m12 - dtau*sign(psi_n).*(G-1);
        psi_m32 = psi_m12 - dtau*SignPsi.*(G-1).*fl;    

        psi = (psi + 2* psi_m32)/3;

    end
    
    
end
LS.psi = psi;
[LS] = LSnormals(LS,DOMAIN);



