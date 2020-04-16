function BoundaryCurve = uniMineralroughFracture(DOMAIN, ...
    b, fy, lambda, a, x_0, refine)
%% ======= rough fracture ===============
xt = 0:min(min(DOMAIN.dxp))/refine:DOMAIN.lx; % x-fracture

BoundaryCurve.yt = b/2 + fy(xt, a, lambda, x_0);
BoundaryCurve.yb = -b/2 - fy(xt, a, lambda, x_0);
%%%============================================================

%%
BoundaryCurve.a = a;
BoundaryCurve.lambda = lambda;
BoundaryCurve.b = b;
BoundaryCurve.x_0 = x_0;
BoundaryCurve.xt = xt;

BoundaryCurve.boundfunctionTop = ...
    @(x, a, lambda, x_0, b) b/2 + a*sin(2*pi/lambda * (x-x_0));
BoundaryCurve.bFunctionBottom = ...
    @(x, a, lambda, x_0, b) -b/2 - a*sin(2*pi/lambda * (x-x_0));

BoundaryCurve.endBotIndx = ceil(length(DOMAIN.yp)/2);
%%
