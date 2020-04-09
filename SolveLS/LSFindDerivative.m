%% =======================
% This function finds the negative or positive direction derivative 
% WENO 5th order. 
% The inputs are:
% - the calculated first order derivatives
% - the weights
% - the equation type: level set equation or reinitialization equation

% ( so far this is limited to uniform grids)

%%


function [psinp] = LSFindDerivative(u,v,psi,DOMAIN,equation, h)



[dir] = LSdirDerivates(DOMAIN,psi,u,v,equation, h);
[weight] = LSWeights(dir,equation,eps);


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
    
    
    wx1 = weight.wx1;
    wx2 = weight.wx2;
    wx3 = weight.wx3;

    wy1 = weight.wy1;
    wy2 = weight.wy2;
    wy3 = weight.wy3;
    
    psinp.psi_x = wx1.*(1/3*vx1-7/6*vx2+11/6*vx3)+ ...
        wx2.*(-1/6*vx2+5/6*vx3+1/3*vx4)+wx3.*(1/3*vx3+5/6*vx4-1/6*vx5);
    
    psinp.psi_y = wy1.*(1/3*vy1-7/6*vy2+11/6*vy3)+ ...
        wy2.*(-1/6*vy2+5/6*vy3+1/3*vy4)+wy3.*(1/3*vy3+5/6*vy4-1/6*vy5);
        
    
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

    
    
    