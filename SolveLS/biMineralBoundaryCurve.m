function [BoundaryCurve] = biMineralBoundaryCurve(DOMAIN, fy, b, lambda, a, x_0,...
    refine_x, refine_y, xt_s, xt_e, xb_s, xb_e)
%% 
% This function creates the boundary curve for bi mineral fracture surfaces 

dx_ref = min(min(DOMAIN.dxp))/refine_x;
dy_ref = min(min(DOMAIN.dyp))/refine_y;
xt = [];
yt = [];

xb = [];
yb = [];



%%
for i = 1:length(xt_s)
    if xt_s(i) ~= 0  % start vertical line
        y_cons_value = b/2 + fy(xt_s(i), a, lambda, x_0);
        y_cons_len = floor((DOMAIN.ly/2 - y_cons_value)/dy_ref);
        x_cons = zeros(1, y_cons_len);
        x_cons(:) = xt_s(i);
        y_cons = y_cons_value+dy_ref:dy_ref:DOMAIN.ly/2;
        xt = [xt x_cons];
        yt = [yt y_cons];
    end
    x_sin = xt_s(i):dx_ref:xt_e(i);
    y_sin = b/2 + fy(x_sin, a, lambda, x_0);
    xt = [xt x_sin];
    yt = [yt y_sin];
    
    if xt_s(i) ~= DOMAIN.lx
        y_cons_value = b/2 + fy(xt_e(i), a, lambda, x_0);
        y_cons_len = floor((DOMAIN.ly/2 - y_cons_value)/dy_ref);
        x_cons = zeros(1, y_cons_len);
        x_cons(:) = xt_e(i);
        xt = [xt x_cons];

        y_cons = y_cons_value+dy_ref:dy_ref:DOMAIN.ly/2;

        yt = [yt y_cons];
    end
end
    
 %% Bottom fracture curve
for i = 1:length(xb_s)
    if xb_s(i) ~= 0  % start vertical line
        y_cons_value = -b/2 - fy(xb_s(i), a, lambda, x_0);
        y_cons_len = floor((y_cons_value + DOMAIN.ly/2)/dy_ref);
        x_cons = zeros(1, y_cons_len);
        x_cons(:) = xb_s(i);
        y_cons = -DOMAIN.ly/2:dy_ref:y_cons_value - dy_ref;
        xb = [xb x_cons];
        yb = [yb y_cons];
    end
    x_sin = xb_s(i):dx_ref:xb_e(i);
    y_sin = -b/2 - fy(x_sin, a, lambda, x_0);
    xb = [xb x_sin];
    yb = [yb y_sin];
    
    if xb_s(i) ~= DOMAIN.lx
        y_cons_value = -b/2 - fy(xb_e(i), a, lambda, x_0);
        y_cons_len = floor((y_cons_value + DOMAIN.ly/2)/dy_ref);
        x_cons = zeros(1, y_cons_len);
        x_cons(:) = xb_e(i);
        y_cons = -DOMAIN.ly/2:dy_ref:y_cons_value - dy_ref;
        xb = [xb x_cons];
        yb = [yb y_cons];
    end
end

%%
BoundaryCurve.xt = xt;
BoundaryCurve.xb = xb;
BoundaryCurve.yt = yt;
BoundaryCurve.yb = yb;

BoundaryCurve.a = a;
BoundaryCurve.lambda = lambda;
BoundaryCurve.b = b;
BoundaryCurve.x_0 = x_0;

BoundaryCurve.xt_s = xt_s;
BoundaryCurve.xt_e = xt_e;
BoundaryCurve.xb_s = xb_s;
BoundaryCurve.xb_e = xb_e;

BoundaryCurve.boundfunctionTop = ...
    @(x, a, lambda, x_0, b) b/2 + a*sin(2*pi/lambda * (x-x_0));
BoundaryCurve.bFunctionBottom = ...
    @(x, a, lambda, x_0, b) -b/2 - a*sin(2*pi/lambda * (x-x_0));
BoundaryCurve.endBotIndx = ceil(length(DOMAIN.yp)/2);
