L = sqrt(diff(LSCase.BoundaryCurve.xt).^2 + ...
    diff(LSCase.BoundaryCurve.yt).^2);

SFR = (sum(L) - 2*DOMAIN.x_0)/(DOMAIN.lx - 2*DOMAIN.x_0)


dpsi = (LS.psi<0) - (psi<0);
DR = nnz(dpsi)*DOMAIN.dxp(1,1)^2
phi_out = 0.5*(phi(end-10,:)+phi(end-11, :));
q_out = U(end-11,:).*phi_out;
% plot(q_out)
phi_in = 0.5*(phi(3,:)+phi(4,:));
q_in = U(3,:).*phi_in;
RR = sum(q_in)- sum(q_out)

figure(1);
contourf(DOMAIN.Xp,...
DOMAIN.Yp,phi.*(psi>0),50, 'LineStyle', 'none'); colormap jet
% contourf(DOMAIN.Xp,DOMAIN.Yp,phi,50, 'LineStyle', 'none'); colormap jet
% hold on
% scatter(LSCase.BoundaryCurve.xt, LSCase.BoundaryCurve.yt+0.4,2.5,'w','filled')
% hold on
% scatter(LSCase.BoundaryCurve.xt, LSCase.BoundaryCurve.yb+0.4,2.5,'w','filled')
axis equal



% figure(2);
% contourf(DOMAIN.Xp(:,3:end-4),DOMAIN.Yp(:,3:end-4),...
%     phi(:,3:end-4).*(psi(:,3:end-4)>0),50, 'LineStyle', 'none'); colormap jet
% % contourf(DOMAIN.Xp,DOMAIN.Yp,phi,50, 'LineStyle', 'none'); colormap jet
% % hold on
% % scatter(LSCase.BoundaryCurve.xt, LSCase.BoundaryCurve.yt+0.4,2.5,'w','filled')
% % hold on
% % scatter(LSCase.BoundaryCurve.xt, LSCase.BoundaryCurve.yb+0.4,2.5,'w','filled')
% axis equal

figure(3);
contourf(DOMAIN.Xp,...
DOMAIN.Yp,(LS.psi<0)-(psi<0),50, 'LineStyle', 'none'); colormap jet
% contourf(DOMAIN.Xp,DOMAIN.Yp,phi,50, 'LineStyle', 'none'); colormap jet
% hold on
% scatter(LSCase.BoundaryCurve.xt, LSCase.BoundaryCurve.yt+0.4,2.5,'w','filled')
% hold on
% scatter(LSCase.BoundaryCurve.xt, LSCase.BoundaryCurve.yb+0.4,2.5,'w','filled')
axis equal


figure(4);
contourf(DOMAIN.Xu,...
DOMAIN.Yu,U,50, 'LineStyle', 'none'); colormap jet
% contourf(DOMAIN.Xp,DOMAIN.Yp,phi,50, 'LineStyle', 'none'); colormap jet
% hold on
% scatter(LSCase.BoundaryCurve.xt, LSCase.BoundaryCurve.yt+0.4,2.5,'w','filled')
% hold on
% scatter(LSCase.BoundaryCurve.xt, LSCase.BoundaryCurve.yb+0.4,2.5,'w','filled')
axis equal

figure(5);
contourf(DOMAIN.Xp,...
DOMAIN.Yp,psi>0,50, 'LineStyle', 'none'); colormap jet
% contourf(DOMAIN.Xp,DOMAIN.Yp,phi,50, 'LineStyle', 'none'); colormap jet
% hold on
% scatter(LSCase.BoundaryCurve.xt, LSCase.BoundaryCurve.yt+0.4,2.5,'w','filled')
% hold on
% scatter(LSCase.BoundaryCurve.xt, LSCase.BoundaryCurve.yb+0.4,2.5,'w','filled')
axis equal

