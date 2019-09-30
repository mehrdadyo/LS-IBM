function [PSI,VOR]=PlotFlowField(StateVar,IBM,DOMAIN)


%% Retrieve DATA
U = StateVar.U;
V = StateVar.V;
P = StateVar.P;
phi = StateVar.phi;

xp = DOMAIN.xp;
yp = DOMAIN.yp;
dxp = DOMAIN.dxp;
dyp = DOMAIN.dyp;
CoEWu = DOMAIN.CoEWu;
CoNSv = DOMAIN.CoNSv;

diamcyl = IBM.diamcyl;
nrgrainx = IBM.nrgrainx;
nrgrainy = IBM.nrgrainy;

xc = IBM.xc;
yc = IBM.yc;

% for i=1:length(xp)-1
%     for j=1:length(yp)-1
%         if IBMu.flag_u == 1
%             U(i,j)=0;
%         end
%     end
% end
% 
% for i=1:length(xp)-1
%     for j=1:length(yp)-1
%         if IBMv.flag_v == 1
%             V(i,j)=0;
%         end
%     end
% end


r=diamcyl/2;
ang=0:0.01:2*pi; 
xcir=r*cos(ang);
ycir=r*sin(ang);
UP= U(2:end,2:end-1).*CoEWu(:,2:end-1) +...
    U(1:end-1,2:end-1).*(1-CoEWu(:,2:end-1));

VP= V(2:end-1,2:end).*CoNSv(2:end-1,:) +...
    V(2:end-1,1:end-1).*(1-CoNSv(2:end-1,:));
% for j = 2:jmax
%   for i = 2:imax
%     UP(i,j) = 0.5*(U(i,j) + U(i-1,j));
%     VP(i,j) = 0.5*(V(i,j) + V(i,j-1));
%   end
% end

%% Streamfunction contours
[PSI] = STREAMFUNCTION(dxp,dyp,UP,VP);
figure(1)
% contour(xp(1:end-1),yp(1:end-1),PSI',200)

contour(xp(1:end-1),yp(1:end-1),PSI',200)
% contour(xp(25:68),yp(13:62),PSI(25:68,13:62)',200)

% contour(xp(ilow:ihigh),yp(jlow:jhigh),PSI(ilow:ihigh,jlow:jhigh)',200);
% contour(xp(1:end-1),yp(1:end-1),PSI',200);
colormap('jet')
hold on
for i=1:nrgrainx
    for j=1:nrgrainy
        pos = [xc(i,j)-r yc(i,j)-r diamcyl diamcyl];
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor','w',...
            'EdgeColor','w')
        hold on
        plot(xc(i,j)+xcir,yc(i,j)+ycir,'k');
    end

    
end

axis equal

%% U- velocity contour
figure(2)

contourf(xp(2:end-1),yp(2:end-1),UP',200,...
                'LineStyle','none');
colormap('jet')
            
hold on
for i=1:nrgrainx
    for j=1:nrgrainy
        pos = [xc(i,j)-r yc(i,j)-r diamcyl diamcyl];
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor','w',...
            'EdgeColor','w')
        hold on
        plot(xc(i,j)+xcir,yc(i,j)+ycir,'k');
    end

    
end
axis equal
%% V- velocity contours
figure(3)
contourf(xp(2:end-1),yp(2:end-1),VP',200,...
                'LineStyle','none');
colormap('jet')
            
hold on
for i=1:nrgrainx
    for j=1:nrgrainy
        pos = [xc(i,j)-r yc(i,j)-r diamcyl diamcyl];
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor','w',...
            'EdgeColor','w')
        hold on
        plot(xc(i,j)+xcir,yc(i,j)+ycir,'k');
    end

    
end
axis equal
%% Pressure Contours 
figure(4)

            
contourf(xp(2:end-1),yp(2:end-1),P(2:end-1,2:end-1)',200,...
                'LineStyle','none');
% contourf(xp(25:68),yp(13:62),P(25:68,13:62)',200,...
%                 'LineStyle','none');


colormap('jet')

hold on
for i=1:nrgrainx
    for j=1:nrgrainy
        pos = [xc(i,j)-r yc(i,j)-r diamcyl diamcyl];
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor','w',...
            'EdgeColor','w')
        hold on
        plot(xc(i,j)+xcir,yc(i,j)+ycir,'k');
    end

    
end
axis equal
%% Vorticity Contours 
figure(5)

[VOR] = VORTICITY(dxp,dyp,UP,VP);
contour(xp(1:end-1),yp(1:end-1),VOR',200)
% contour(xp(25:68),yp(13:62),VOR(25:68,13:62)',200)

colormap('jet')


hold on
for i=1:nrgrainx
    for j=1:nrgrainy
        pos = [xc(i,j)-r yc(i,j)-r diamcyl diamcyl];
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor','w',...
            'EdgeColor','w')
        hold on
        plot(xc(i,j)+xcir,yc(i,j)+ycir,'k');
    end

    
end
axis equal
%% Velocity magnitude Contours 
figure(6)

UV=sqrt(UP.^2+VP.^2);
contourf(xp(2:end-1),yp(2:end-1),UV',200,'LineStyle','none');
colormap('jet')


hold on
for i=1:nrgrainx
    for j=1:nrgrainy
        pos = [xc(i,j)-r yc(i,j)-r diamcyl diamcyl];
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor','w',...
            'EdgeColor','w')
        hold on
        plot(xc(i,j)+xcir,yc(i,j)+ycir,'k');
    end

    
end
axis equal

%% Concentration contours
figure(7)

            
contourf(xp(2:end-1),yp(2:end-1),phi(2:end-1,2:end-1)',200,...
                'LineStyle','none');
% contourf(xp(25:68),yp(13:62),P(25:68,13:62)',200,...
%                 'LineStyle','none');


colormap('jet')

hold on
for i=1:nrgrainx
    for j=1:nrgrainy
        pos = [xc(i,j)-r yc(i,j)-r diamcyl diamcyl];
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor','w',...
            'EdgeColor','w')
        hold on
        plot(xc(i,j)+xcir,yc(i,j)+ycir,'k');
    end

    
end
axis equal