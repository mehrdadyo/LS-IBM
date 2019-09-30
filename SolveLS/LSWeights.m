%% =======================
% This function finds the weights required to calculate the directional 
% derivatives by WENO 5th order. 
% The inputs are:
% - the calculated first order derivatives
% - the equation type: level set equation or reinitialization equation
% - epsilon for avoiding division by zero. 


%%

function [Weight] = LSWeights(dir,equation,eps)

if equation == "LevelSetEqn"

    vx1 = dir.vx1;
    vx2 = dir.vx2;
    vx3 = dir.vx3;
    vx4 = dir.vx4;
    vx5 = dir.vx5;

    vy1 = dir.vy1;
    vy2 = dir.vy2;
    vy3 = dir.vy3;
    vy4 = dir.vy4;
    vy5 = dir.vy5;

    %% Smoothness
    Sx1 = 13/12*(vx1-2*vx2+vx3).^2 + 1/4*(vx1-4*vx2+3*vx3).^2;
    Sx2 = 13/12*(vx2-2*vx3+vx4).^2 + 1/4*(vx2-vx4).^2;
    Sx3 = 13/12*(vx3-2*vx4+vx5).^2 + 1/4*(3*vx3-4*vx4+vx5).^2;

    Sy1 = 13/12*(vy1-2*vy2+vy3).^2 + 1/4*(vy1-4*vy2+3*vy3).^2;
    Sy2 = 13/12*(vy2-2*vy3+vy4).^2 + 1/4*(vy2-vy4).^2;
    Sy3 = 13/12*(vy3-2*vy4+vy5).^2 + 1/4*(3*vy3-4*vy4+vy5).^2;

    %% Weights
    ax1 = 1/10*(1./(eps+Sx1).^2);
    ax2 = 6/10*(1./(eps+Sx2).^2);
    ax3 = 3/10*(1./(eps+Sx3).^2);

    ay1 = 1/10*(1./(eps+Sy1).^2);
    ay2 = 6/10*(1./(eps+Sy2).^2);
    ay3 = 3/10*(1./(eps+Sy3).^2);

    wx1 = ax1./(ax1+ax2+ax3);
    wx2 = ax2./(ax1+ax2+ax3);
    wx3 = ax3./(ax1+ax2+ax3);

    wy1 = ay1./(ay1+ay2+ay3);
    wy2 = ay2./(ay1+ay2+ay3);
    wy3 = ay3./(ay1+ay2+ay3);


    %% update struct 
    Weight.wx1 = wx1;
    Weight.wx2 = wx2;
    Weight.wx3 = wx3;

    Weight.wy1 = wy1;
    Weight.wy2 = wy2;
    Weight.wy3 = wy3;
    
elseif equation == "ReinitializationEqn"
    
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
