%% =======================
% This function finds the negative or positive direction derivative 
% WENO 5th order. 
% The inputs are:
% - the calculated first order derivatives
% - the weights
% - the equation type: level set equation or reinitialization equation

% ( so far this is limited to uniform grids)

%%


function [psinp] = LSFindDerivativeFV(psi,DOMAIN,equation,StateVar)

dx = min(min(DOMAIN.dxp));

[dir] = LSdirDerivates(dx,psi,equation,StateVar);
[weight] = LSWeights(dir,equation,eps);

if equation == "LevelSetEqn"
    
    psi1e = dir.psi1e;
    psi2e = dir.psi2e;
    psi3e = dir.psi3e;
    psi4e = dir.psi4e;
    psi5e = dir.psi5e;
    
    psi1w = dir.psi1w;
    psi2w = dir.psi2w;
    psi3w = dir.psi3w;
    psi4w = dir.psi4w;
    psi5w = dir.psi5w;

    psi1n = dir.psi1n;
    psi2n = dir.psi2n;
    psi3n = dir.psi3n;
    psi4n = dir.psi4n;
    psi5n = dir.psi5n;
    
    psi1s = dir.psi1s;
    psi2s = dir.psi2s;
    psi3s = dir.psi3s;
    psi4s = dir.psi4s;
    psi5s = dir.psi5s;
    
    %% weights 
    
    we1 = weight.we1;
    we2 = weight.we2;
    we3 = weight.we3;

    ww1 = weight.ww1;
    ww2 = weight.ww2;
    ww3 = weight.ww3;
    
    wn1 = weight.wn1;
    wn2 = weight.wn2;
    wn3 = weight.wn3;
    
    ws1 = weight.ws1;
    ws2 = weight.ws2;
    ws3 = weight.ws3;
    
    psinp.psie = ...
        we1.*(1/3*psi1e-7/6*psi2e+11/6*psi3e)+ ...
        we2.*(-1/6*psi2e+5/6*psi3e+1/3*psi4e)+ ...
        we3.*(1/3*psi3e+5/6*psi4e-1/6*psi5e);
    
    psinp.psiw = ...
        ww1.*(1/3*psi1w-7/6*psi2w+11/6*psi3w)+ ...
        ww2.*(-1/6*psi2w+5/6*psi3w+1/3*psi4w)+ ...
        ww3.*(1/3*psi3w+5/6*psi4w-1/6*psi5w);
        
    psinp.psin = ...
        wn1.*(1/3*psi1n-7/6*psi2n+11/6*psi3n)+ ...
        wn2.*(-1/6*psi2n+5/6*psi3n+1/3*psi4n)+ ...
        wn3.*(1/3*psi3n+5/6*psi4n-1/6*psi5n);
    
    psinp.psis = ...
        ws1.*(1/3*psi1s-7/6*psi2s+11/6*psi3s)+ ...
        ws2.*(-1/6*psi2s+5/6*psi3s+1/3*psi4s)+ ...
        ws3.*(1/3*psi3s+5/6*psi4s-1/6*psi5s);
    
elseif equation == 'ReinitializationEqn'
    
    vxn1 = dir.vxn1;
    vxn2 = dir.vxn2;
    vxn3 = dir.vxn3;
    vxn4 = dir.vxn4;
    vxn5 = dir.vxn5;
    
    vxp1 = dir.vxp1;
    vxp2 = dir.vxp2;
    vxp3 = dir.vxp3;
    vxp4 = dir.vxp4;
    vxp5 = dir.vxp5;

    vyn1 = dir.vyn1;
    vyn2 = dir.vyn2;
    vyn3 = dir.vyn3;
    vyn4 = dir.vyn4;
    vyn5 = dir.vyn5;

    vyp1 = dir.vyp1;
    vyp2 = dir.vyp2;
    vyp3 = dir.vyp3;
    vyp4 = dir.vyp4;
    vyp5 = dir.vyp5;
    
    
    wxn1 = weight.wxn1;
    wxn2 = weight.wxn2;
    wxn3 = weight.wxn3;
    
    wxp1 = weight.wxp1;
    wxp2 = weight.wxp2;
    wxp3 = weight.wxp3;

    wyn1 = weight.wyn1;
    wyn2 = weight.wyn2;
    wyn3 = weight.wyn3;
    
    wyp1 = weight.wyp1;
    wyp2 = weight.wyp2;
    wyp3 = weight.wyp3;
    
    
    psinp.psi_xn = wxn1.*(1/3*vxn1-7/6*vxn2+11/6*vxn3)+ ...
        wxn2.*(-1/6*vxn2+5/6*vxn3+1/3*vxn4)+...
        wxn3.*(1/3*vxn3+5/6*vxn4-1/6*vxn5);
    
    psinp.psi_yn = wyn1.*(1/3*vyn1-7/6*vyn2+11/6*vyn3)+ ...
        wyn2.*(-1/6*vyn2+5/6*vyn3+1/3*vyn4)+...
        wyn3.*(1/3*vyn3+5/6*vyn4-1/6*vyn5);
    
    
    psinp.psi_xp = wxp1.*(1/3*vxp1-7/6*vxp2+11/6*vxp3)+ ...
        wxp2.*(-1/6*vxp2+5/6*vxp3+1/3*vxp4)+...
        wxp3.*(1/3*vxp3+5/6*vxp4-1/6*vxp5);
    
    psinp.psi_yp = wyp1.*(1/3*vyp1-7/6*vyp2+11/6*vyp3)+ ...
        wyp2.*(-1/6*vyp2+5/6*vyp3+1/3*vyp4)+...
        wyp3.*(1/3*vyp3+5/6*vyp4-1/6*vyp5);
    
end

    
    
    