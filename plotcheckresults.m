psi = LS.psi;
dx = min(min(DOMAIN.dxp));

psi = psi.*(abs(psi)<10*dx);

figure(2)
r= IBM.diamcyl/2;
ang=0:0.01:2*pi;
xcir=r*cos(ang);
ycir=r*sin(ang);
xc = IBM.xc;
yc = IBM.yc;
figure(13)
contourf(DOMAIN.xp(1:end-10),DOMAIN.yp,(psi(1:end-10,:))',200,...
    'LineStyle','none');
colormap jet
axis equal
hold on
plot(xc+xcir,yc+ycir,'k');
hold on
plot([xc, 0], [yc, 0])

for i =2:imax
    for j = 2:jmax
        %         if abs(psi(i,j))<10*dx
        nx(i,j) = (psi(i+1,j)-psi(i-1,j))/(2*dx);
        ny(i,j) = (psi(i,j+1)-psi(i,j-1))/(2*dy);
        grad = sqrt(nx(i,j)^2 + ny(i,j)^2);
        %             nx(i,j) = nx(i,j)/grad;
        %             ny(i,j) = ny(i,j)/grad;
        %         end
    end`
end

dpsi = sqrt(nx.^2+ny.^2);


[row,col,valu] = find(LS.u);

imax = DOMAIN.imax;
jmax = DOMAIN.jmax;

psi = LS.psi;

[nx,ny,grad] = deal(zeros(imax+1,jmax+1));

dx = min(min(DOMAIN.dxp));
dy = min(min(DOMAIN.dyp));

for i =2:imax
    for j = 2:jmax
        
        if abs(psi_n(i,j))<10*dx  

            nx(i,j) = (psi_n(i+1,j)-psi_n(i-1,j))/(2*dx);
            ny(i,j) = (psi_n(i,j+1)-psi_n(i,j-1))/(2*dy);
            grad(i,j) = sqrt(nx(i,j)^2 + ny(i,j)^2);

            
            
        end

    end
end

figure(5)
contourf(DOMAIN.xp(48:168),DOMAIN.yp(78:198),...
    phi(48:168,78:198)',200,'LineStyle','none');
colormap jet
colorbar
axis equal
hold on
plot(xc+xcir,yc+ycir,'w');


figure(5)
contourf(DOMAIN.xp,DOMAIN.yp,...
    nx',200,'LineStyle','none');
colormap jet
colorbar
axis equal
hold on
plot(xc+xcir,yc+ycir,'w');

[I,J,valu] = find(u);
plot(valu)

rhsu = reshape(Soln.RHS_U, DOMAIN.imax-2, DOMAIN.jmax-1);
figure(14)
contourf(DOMAIN.xu(1:end-2),DOMAIN.yu(1:end-2),rhsu',200,...
    'LineStyle','none');
colormap jet
axis equal
(imax-2)*(jmax-1)


scatter(DOMAIN.xu(IBM_coeffU.I_g_u), DOMAIN.yu(IBM_coeffU.J_g_u))