function [X_ib_s,Y_ib_s,dS_s,X_forc,Y_forc,X_virt1,Y_virt1,...
    X_virt2,Y_virt2,X_virt1_x,Y_virt1_x,X_virt1_y,Y_virt1_y,nx,ny]=...
    Scalarpoints(X_ib_g_p,Y_ib_g_p,xc,yc,dx,dy)

% This function gives one the points needed to interpolate the pressure and
% velocity close to the boundary in order to calculate the forces as well
% as other needed parameters. 


I=zeros(length(X_ib_g_p),1);
[d,d_s]=deal(zeros(length(X_ib_g_p),length(X_ib_g_p)));
i_s=zeros(length(X_ib_g_p),2);
for i=1:length(X_ib_g_p)
    % find all the distances of lagrangian points from each other
    d(i,:)=sqrt((X_ib_g_p(i)-X_ib_g_p).^2+(Y_ib_g_p(i)-Y_ib_g_p).^2);  
    d_s(i,:)=sort(d(i,:));      % sort the distances and put them in rows             
    i_s(i,1)=find(d(i,:)==d_s(i,2));    % find first closest point
    i_s(i,2)=find(d(i,:)==d_s(i,3));    % find second closest point
end

for i=1:length(X_ib_g_p)
    if i==1
        I(i)=i;
    elseif i==2
        I(i)=i_s(i-1,1);
    else
         if i_s(I(i-1),1)==I(i-2)   % check the direction of the points
             I(i)=i_s(I(i-1),2);
         else 
             I(i)=i_s(I(i-1),1);
         end
    end
        
end



X_ib_s=X_ib_g_p(I);
Y_ib_s=Y_ib_g_p(I);
dS_s=[sqrt((X_ib_s(2:end)-X_ib_s(1:end-1)).^2+...
    (Y_ib_s(2:end)-Y_ib_s(1:end-1)).^2) sqrt((X_ib_s(end)-X_ib_s(1)).^2+...
    (Y_ib_s(end)-Y_ib_s(1)).^2)];


X_forc=1/2*[X_ib_s(2:end)+X_ib_s(1:end-1) X_ib_s(end)+X_ib_s(1)];
Y_forc=1/2*[Y_ib_s(2:end)+Y_ib_s(1:end-1) Y_ib_s(end)+Y_ib_s(1)];

%% find two exterior virtual  points
treshhold=1e-5;
Hypo=sqrt((X_forc-xc).^2+(Y_forc-yc).^2);
nx=(X_forc-xc)./Hypo;
ny=(Y_forc-yc)./Hypo;
nx(abs(nx)<treshhold)=0;
ny(abs(ny)<treshhold)=0;
X_virt1=sqrt(2)*dx*nx+X_forc;
Y_virt1=sqrt(2)*dy*ny+Y_forc;
X_virt2=1.5*sqrt(2)*dx*nx+X_forc;
Y_virt2=1.5*sqrt(2)*dy*ny+Y_forc;
X_virt1_x=sqrt(2)*dx*nx+X_forc;
Y_virt1_x=Y_forc;
X_virt1_y=X_forc;
Y_virt1_y=sqrt(2)*dy*ny+Y_forc;

[f]=find(ny==min(abs(ny)));
f=f(1);
I=[I(f:end); I(1:f-1)];

X_ib_s=[X_ib_s(f:end) X_ib_s(1:f-1)];
Y_ib_s=[Y_ib_s(f:end) Y_ib_s(1:f-1)];
dS_s=[dS_s(f:end) dS_s(1:f-1)];
X_forc=[X_forc(f:end) X_forc(1:f-1)];
Y_forc=[Y_forc(f:end) Y_forc(1:f-1)];
X_virt1=[X_virt1(f:end) X_virt1(1:f-1)];
Y_virt1=[Y_virt1(f:end) Y_virt1(1:f-1)];
X_virt2=[X_virt2(f:end) X_virt2(1:f-1)];
Y_virt2=[Y_virt2(f:end) Y_virt2(1:f-1)];
X_virt1_x=[X_virt1_x(f:end) X_virt1_x(1:f-1)];
Y_virt1_x=[Y_virt1_x(f:end) Y_virt1_x(1:f-1)];
X_virt1_y=[X_virt1_y(f:end) X_virt1_x(1:f-1)];
Y_virt1_y=[Y_virt1_y(f:end) Y_virt1_y(1:f-1)];
nx=[nx(f:end) nx(1:f-1)];
ny=[ny(f:end) ny(1:f-1)];



% I=find(xu<X_virt1_y(1));  %% Find the i- index of the points around virtual Point
% I=I(end);
% J=find(yu<Y_virt1_y(1));  %% Find the j- index of the points around virtual Point
% J=J(end);
% 
% 
% scatter(X_virt1(1:5),Y_virt1(1:5))
% hold on
% scatter(X_virt1_x(1:2),Y_virt1_x(1:2))
% hold on
% scatter(X_virt1_y(1:4),Y_virt1_y(1:4))
% hold on
% scatter(X_forc(1:3),Y_forc(1:3))
% hold on
% scatter(xu(I),yu(J),'x')
% hold on
% scatter(xu(I+1),yu(J),'x')
% hold on
% scatter(xu(I),yu(J+1),'x')
% hold on
% scatter(xu(I+1),yu(J+1),'x')
% 
% axis equal




%% find bilinear interpolation coefficients for exterior virtual points

