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
    
contourf(DOMAIN.Xp(11:end-10,:),DOMAIN.Yp(11:end-10,:),(psi(11:end-10,:)>0).*phi(11:end-10,:),20, 'LineStyle', 'none'); colormap jet

figure(2);contourf(DOMAIN.Xp ,DOMAIN.Yp,phi.*(psi>0),50, 'LineStyle', 'none'); colormap jet
% contourf(DOMAIN.Xp,DOMAIN.Yp,LS.psi,50, 'LineStyle', 'none'); colormap jet
hold on
plot(LSCase.BoundaryCurve.xt, LSCase.BoundaryCurve.yt+0.7,'w','LineWidth', 0.5)
hold on
plot(LSCase.BoundaryCurve.xt, LSCase.BoundaryCurve.yb+0.7,'w','LineWidth', 0.5)
axis equal
xlim([DOMAIN.xp(25), DOMAIN.xp(end-10)])   

dyb = diff(LSCase.BoundaryCurve.yb);
dx = diff(LSCase.BoundaryCurve.xt);
dyt = diff(LSCase.BoundaryCurve.yt);
sum(sqrt(dx.^2+dyt.^2))*2


% load('G:\My Drive\8 - Stanford\Research\4th Paper\Simulations\NewCases\Case8\Output1\dataRDE1390dt.mat')
load('G:\My Drive\8 - Stanford\Research\4th Paper\Simulations\NewCases\Case33\Output\dataRDE8000dt.mat')
sum((psi(1,:)>0))
sum((psi(end-44,:)>0))
phi_f = 0.5*(phi(1:end-1,:) + phi(2:end,:));
q = U.*phi_f;
R = ( sum(q(1,:)) - sum(q(end-44,:)) ) * 0.004


StateVar.U = U;
StateVar.V = V;
StateVar.P = P;
LS.psi = psi;
StateVar.U_old = U;
StateVar.V_old = V;
StateVar.P_old = P;
StateVar.phi = phi;
StateVar.phi_old = phi;
ControlVar.time = time;
LS.psi = psi;
[LS] = LSnormals(LS,DOMAIN);
[IBM_coeffU,IBM_coeffV,IBM_coeffP] = LSIBMcoeffs(IBM,DOMAIN,LS, StateVar.phi);
VARIABLES.dt = VARIABLES.dt/ 10;


nLSupdate = VARIABLES.nLSupdate;
VARIABLES.nLSupdate = 1;

load('G:\My Drive\8 - Stanford\Research\4th Paper\Simulations\NewCases\Case30\Output\dataRDE8000dt.mat')
