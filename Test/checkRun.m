uvI = sqrt(uI.^2+vI.^2);



figure(1)
contourf(DOMAIN.xp(1:end-10),DOMAIN.yp,(phi(1:end-10,:))',200,...
    'LineStyle','none');
colormap jet;colorbar
axis equal
hold on
plot(xc+xcir,yc+ycir,'k');

figure(2)
contourf(DOMAIN.xp(1:end-10),DOMAIN.yp,(psi(1:end-10,:)<0)',200,...
    'LineStyle','none');
colormap jet;colorbar
axis equal
hold on
plot(xc+xcir,yc+ycir,'k');

figure(3)
contourf(DOMAIN.xu(1:end-10),DOMAIN.yu,(U(1:end-10,:))',200,...
    'LineStyle','none');
colormap jet;colorbar
axis equal
hold on
plot(xc+xcir,yc+ycir,'k');

figure(4)
contourf(DOMAIN.xp(1:end-10),DOMAIN.yp,uvI(1:end-10,:)',200,...
    'LineStyle','none');
colormap jet;colorbar
axis equal
hold on
plot(xc+xcir,yc+ycir,'k');

