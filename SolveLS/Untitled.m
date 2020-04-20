
%% ======= rough fracture ===============
Boundary1.xt_s = [0 0.8 1.5 2.5];
Boundary1.xt_e = [0.5 1.2 2 4];

Boundary1.xb_s = [0 2 3.5];
Boundary1.xb_e = [0.8 3 4];


[xt, xb, yt, yb] = boundaryCurve(DOMAIN, fy, b, lambda, a, x_0,...
    refine_x, refine_y, xt_s, xt_e, xb_s, xb_e);

Boundary2.xt_s = [0.5 1.2 2];
Boundary2.xt_e = [0.8 1.5 2.5];

Boundary2.xb_s = [0.8 3];
Boundary2.xb_e = [2 3.5];

[xt2, xb2, yt2, yb2] = boundaryCurve(DOMAIN, fy, b, lambda, a, x_0,...
    refine_x, refine_y, xt_s, xt_e, xb_s, xb_e);
figure(1)
scatter(LSCase.BoundaryCurve1.xt, LSCase.BoundaryCurve1.yt)
hold on
scatter(LSCase.BoundaryCurve1.xb, LSCase.BoundaryCurve1.yb)
ylim([-DOMAIN.ly/2 DOMAIN.ly/2])
xlim([0 DOMAIN.lx])
figure(2)
hold on 
scatter(LSCase.BoundaryCurve2.xt, LSCase.BoundaryCurve2.yt)
hold on
scatter(LSCase.BoundaryCurve2.xb, LSCase.BoundaryCurve2.yb)
ylim([-DOMAIN.ly/2 DOMAIN.ly/2])
xlim([0 DOMAIN.lx])