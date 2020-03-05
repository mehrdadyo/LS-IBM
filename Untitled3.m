figure(3)
contourf(DOMAIN.Xp,DOMAIN.Yp,IBM_coeffP.flag_p, 50, 'LineStyle', 'none'); colormap jet
hold on
plot(LSCase.BoundaryCurve.xt, LSCase.BoundaryCurve.yt+0.7, 'k')
hold on
plot(LSCase.BoundaryCurve.xt, LSCase.BoundaryCurve.yb+0.7, 'k')

hold on
scatter(DOMAIN.xp(IBM_coeffP.I1_p), DOMAIN.yp(IBM_coeffP.J1_p))
hold on 
scatter(DOMAIN.xp(IBM_coeffP.I2_p), DOMAIN.yp(IBM_coeffP.J2_p))
hold on 
scatter(DOMAIN.xp(IBM_coeffP.I3_p), DOMAIN.yp(IBM_coeffP.J3_p))
hold on 
scatter(DOMAIN.xp(IBM_coeffP.I4_p), DOMAIN.yp(IBM_coeffP.J4_p))
hold on 
scatter(DOMAIN.xp(IBM_coeffP.I_g_p), DOMAIN.yp(IBM_coeffP.J_g_p))
hold on 

plot(LSCase.BoundaryCurve.xt,LSCase.BoundaryCurve.yb+DOMAIN.ly/2, 'k')


    [StateVar,ControlVar] = ... 
        SolveUVP (ControlVar,Flux,DOMAIN,VARIABLES,StateVar,IBM,IBM_coeffU,...
        IBM_coeffV,IBM_coeffP,BC);    
    
    StateVar.U_old = StateVar.U;
    StateVar.V_old = StateVar.V;
    StateVar.P_old = StateVar.P;
    
contourf(DOMAIN.Xp(11:end-10,:),DOMAIN.Yp(11:end-10,:),(LS.psi(11:end-10,:)>0).*StateVar.phi(11:end-10,:),20, 'LineStyle', 'none'); colormap jet

% contourf(DOMAIN.Xp,DOMAIN.Yp,LS.psi,50, 'LineStyle', 'none'); colormap jet
hold on
plot(LSCase.BoundaryCurve.xt, LSCase.BoundaryCurve.yt+0.7, 'k')
hold on
plot(LSCase.BoundaryCurve.xt, LSCase.BoundaryCurve.yb+0.7, 'k')
xlim([DOMAIN.xp(11), DOMAIN.xp(end-10)])    