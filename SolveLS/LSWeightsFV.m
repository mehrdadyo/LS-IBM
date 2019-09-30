%% =======================
% This function finds the weights required to calculate the directional 
% derivatives by WENO 5th order. 
% The inputs are:
% - the calculated first order derivatives
% - the equation type: level set equation or reinitialization equation
% - epsilon for avoiding division by zero. 


%%

function [Weight] = LSWeightsFV(dir,equation,eps)

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

    %% Smoothness
    Se1 = 13/12*(psi1e-2*psi2e+psi3e).^2 + 1/4*(psi1e-4*psi2e+3*psi3e).^2;
    Se2 = 13/12*(psi2e-2*psi3e+psi4e).^2 + 1/4*(psi2e-psi4e).^2;
    Se3 = 13/12*(psi3e-2*psi4e+psi5e).^2 + 1/4*(3*psi3e-4*psi4e+psi5e).^2;
    
    Sw1 = 13/12*(psi1w-2*psi2w+psi3w).^2 + 1/4*(psi1w-4*psi2w+3*psi3w).^2;
    Sw2 = 13/12*(psi2w-2*psi3w+psi4w).^2 + 1/4*(psi2w-psi4w).^2;
    Sw3 = 13/12*(psi3w-2*psi4w+psi5w).^2 + 1/4*(3*psi3w-4*psi4w+psi5w).^2;

    Sn1 = 13/12*(psi1n-2*psi2n+psi3n).^2 + 1/4*(psi1n-4*psi2n+3*psi3n).^2;
    Sn2 = 13/12*(psi2n-2*psi3n+psi4n).^2 + 1/4*(psi2n-psi4n).^2;
    Sn3 = 13/12*(psi3n-2*psi4n+psi5n).^2 + 1/4*(3*psi3n-4*psi4n+psi5n).^2;
    
    Ss1 = 13/12*(psi1s-2*psi2s+psi3s).^2 + 1/4*(psi1s-4*psi2s+3*psi3s).^2;
    Ss2 = 13/12*(psi2s-2*psi3s+psi4s).^2 + 1/4*(psi2s-psi4s).^2;
    Ss3 = 13/12*(psi3s-2*psi4s+psi5s).^2 + 1/4*(3*psi3s-4*psi4s+psi5s).^2;

    %% Weights
    ae1 = 1/10*(1./(eps+Se1).^2);
    ae2 = 6/10*(1./(eps+Se2).^2);
    ae3 = 3/10*(1./(eps+Se3).^2);
    
    aw1 = 1/10*(1./(eps+Sw1).^2);
    aw2 = 6/10*(1./(eps+Sw2).^2);
    aw3 = 3/10*(1./(eps+Sw3).^2);

    an1 = 1/10*(1./(eps+Sn1).^2);
    an2 = 6/10*(1./(eps+Sn2).^2);
    an3 = 3/10*(1./(eps+Sn3).^2);
    
    as1 = 1/10*(1./(eps+Ss1).^2);
    as2 = 6/10*(1./(eps+Ss2).^2);
    as3 = 3/10*(1./(eps+Ss3).^2);

    we1 = ae1./(ae1+ae2+ae3);
    we2 = ae2./(ae1+ae2+ae3);
    we3 = ae3./(ae1+ae2+ae3);
    
    ww1 = aw1./(aw1+aw2+aw3);
    ww2 = aw2./(aw1+aw2+aw3);
    ww3 = aw3./(aw1+aw2+aw3);

    wn1 = an1./(an1+an2+an3);
    wn2 = an2./(an1+an2+an3);
    wn3 = an3./(an1+an2+an3);
    
    ws1 = as1./(as1+as2+as3);
    ws2 = as2./(as1+as2+as3);
    ws3 = as3./(as1+as2+as3);


    %% update struct 
    Weight = struct('we1',we1,'we2',we2,'we3',we3, ...
        'ww1',ww1,'ww2',ww2,'ww3',ww3,'wn1',wn1,'wn2',wn2,'wn3',wn3, ...
        'ws1',ws1,'ws2',ws2,'ws3',ws3);
        
elseif equation == "ReinitializationEqn"
    
    % east negative
    psi1eN = dir.psi1eN;
    psi2eN = dir.psi2eN;
    psi3eN = dir.psi3eN;
    psi4eN = dir.psi4eN;
    psi5eN = dir.psi5eN;
    
    % east positive
    psi1eP = dir.psi1eP;
    psi2eP = dir.psi2eP;
    psi3eP = dir.psi3eP;
    psi4eP = dir.psi4eP;
    psi5eP = dir.psi5eP;
    
    % west negative
    psi1wN = dir.psi1wN;
    psi2wN = dir.psi2wN;
    psi3wN = dir.psi3wN;
    psi4wN = dir.psi4wN;
    psi5wN = dir.psi5wN;
    
    % west positive
    psi1wP = dir.psi1wP;
    psi2wP = dir.psi2wP;
    psi3wP = dir.psi3wP;
    psi4wP = dir.psi4wP;
    psi5wP = dir.psi5wP;

    % north negative
    psi1nN = dir.psi1nN;
    psi2nN = dir.psi2nN;
    psi3nN = dir.psi3nN;
    psi4nN = dir.psi4nN;
    psi5nN = dir.psi5nN;
    
    % north positive
    psi1nP = dir.psi1nP;
    psi2nP = dir.psi2nP;
    psi3nP = dir.psi3nP;
    psi4nP = dir.psi4nP;
    psi5nP = dir.psi5nP;
    
    % south negative
    psi1sN = dir.psi1sN;
    psi2sN = dir.psi2sN;
    psi3sN = dir.psi3sN;
    psi4sN = dir.psi4sN;
    psi5sN = dir.psi5sN;
    
    % south positive
    psi1sP = dir.psi1sP;
    psi2sP = dir.psi2sP;
    psi3sP = dir.psi3sP;
    psi4sP = dir.psi4sP;
    psi5sP = dir.psi5sP;
    %% Smoothness
    Sxn1 = 13/12*(vxn1-2*vxn2+vxn3).^2 + 1/4*(vxn1-4*vxn2+3*vxn3).^2;
    Sxn2 = 13/12*(vxn2-2*vxn3+vxn4).^2 + 1/4*(vxn2-vxn4).^2;
    Sxn3 = 13/12*(vxn3-2*vxn4+vxn5).^2 + 1/4*(3*vxn3-4*vxn4+vxn5).^2;
    
    Sxp1 = 13/12*(vxp1-2*vxp2+vxp3).^2 + 1/4*(vxp1-4*vxp2+3*vxp3).^2;
    Sxp2 = 13/12*(vxp2-2*vxp3+vxp4).^2 + 1/4*(vxp2-vxp4).^2;
    Sxp3 = 13/12*(vxp3-2*vxp4+vxp5).^2 + 1/4*(3*vxp3-4*vxp4+vxp5).^2;

    Syn1 = 13/12*(vyn1-2*vyn2+vyn3).^2 + 1/4*(vyn1-4*vyn2+3*vyn3).^2;
    Syn2 = 13/12*(vyn2-2*vyn3+vyn4).^2 + 1/4*(vyn2-vyn4).^2;
    Syn3 = 13/12*(vyn3-2*vyn4+vyn5).^2 + 1/4*(3*vyn3-4*vyn4+vyn5).^2;
    
    Syp1 = 13/12*(vyp1-2*vyp2+vyp3).^2 + 1/4*(vyp1-4*vyp2+3*vyp3).^2;
    Syp2 = 13/12*(vyp2-2*vyp3+vyp4).^2 + 1/4*(vyp2-vyp4).^2;
    Syp3 = 13/12*(vyp3-2*vyp4+vyp5).^2 + 1/4*(3*vyp3-4*vyp4+vyp5).^2;

    %% Weights
    axn1 = 1/10*(1./(eps+Sxn1).^2);
    axn2 = 6/10*(1./(eps+Sxn2).^2);
    axn3 = 3/10*(1./(eps+Sxn3).^2);
    
    axp1 = 1/10*(1./(eps+Sxp1).^2);
    axp2 = 6/10*(1./(eps+Sxp2).^2);
    axp3 = 3/10*(1./(eps+Sxp3).^2);

    ayn1 = 1/10*(1./(eps+Syn1).^2);
    ayn2 = 6/10*(1./(eps+Syn2).^2);
    ayn3 = 3/10*(1./(eps+Syn3).^2);
    
    ayp1 = 1/10*(1./(eps+Syp1).^2);
    ayp2 = 6/10*(1./(eps+Syp2).^2);
    ayp3 = 3/10*(1./(eps+Syp3).^2);

    wxn1 = axn1./(axn1+axn2+axn3);
    wxn2 = axn2./(axn1+axn2+axn3);
    wxn3 = axn3./(axn1+axn2+axn3);
    
    wxp1 = axp1./(axp1+axp2+axp3);
    wxp2 = axp2./(axp1+axp2+axp3);
    wxp3 = axp3./(axp1+axp2+axp3);

    wyn1 = ayn1./(ayn1+ayn2+ayn3);
    wyn2 = ayn2./(ayn1+ayn2+ayn3);
    wyn3 = ayn3./(ayn1+ayn2+ayn3);

    wyp1 = ayp1./(ayp1+ayp2+ayp3);
    wyp2 = ayp2./(ayp1+ayp2+ayp3);
    wyp3 = ayp3./(ayp1+ayp2+ayp3);
    %% update struct 
    Weight.wxn1 = wxn1;
    Weight.wxn2 = wxn2;
    Weight.wxn3 = wxn3;
    
    Weight.wxp1 = wxp1;
    Weight.wxp2 = wxp2;
    Weight.wxp3 = wxp3;

    Weight.wyn1 = wyn1;
    Weight.wyn2 = wyn2;
    Weight.wyn3 = wyn3;
    
    Weight.wyp1 = wyp1;
    Weight.wyp2 = wyp2;
    Weight.wyp3 = wyp3;
    
end
