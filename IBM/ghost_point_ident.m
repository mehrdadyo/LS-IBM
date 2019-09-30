function [landa,landa_g_1,landa_g_2,landa_g_3,landa_g_4,A1_g,numg,flag,...
    I_g,J_g,I_solid,J_solid,I_e,J_e,I1,J1,I2,J2,I3,J3,I4,J4,...
    landa_g_5,landa_g_6,I5,J5,I6,J6,nx_g,ny_g]=...
                 ghost_point_ident(xc,yc,diamcyl,dx,dy,x,y,alpha,beta,q,BQ)
% X=XYp;
% x=xu;
% y=yu;
% -alpha dphi - beta phi = q 
%size of coordinate matrix n=size[x], m=size[y]
Nx=size(x,2);
Ny=size(y,2);
% flag for all points
flag=zeros(Nx,Ny);
r=diamcyl/2;
nrgrain=size(xc,2);


%% reducing size of search domain
xmin=xc-3*r;                  %left
xmax=xc+3*r;                  %right 
ymin=yc-3*r;                  %bottom
ymax=yc+3*r;                  %up
%it's corresponding i and j
imin=find(x<(xmin+dx)& x>(xmin-dx));
imin=imin(1);
imax=find(x<(xmax+dx)& x>(xmax-dx));
imax=imax(1);
jmin=find(y<(ymin+dy)& y>(ymin-dy));
jmin=jmin(1);
jmax=find(y<(ymax+dy)& y>(ymax-dy));
jmax=jmax(1);
% flag=-1 ==>extrapolation cells (fluid)
% flag=1  ==>Mirror or ghost cells (soild)
% flag=0  ==> all others (fluid)
% flag=2  ==> body cells (solid)
% first everywhere fluid (flag=-1), second inside points solid (flag=1),
% third ib points (flag=-1), fourth all other points zero.
% numg=0;  
% numg=0;
% shows number of ib, and lagrangian points
%% Solid points
l=0;
for k=1:nrgrain
    for i=imin(k):imax(k)
        for j=jmin(k):jmax(k)
            d=(x(i)-xc(k))^2+(y(j)-yc(k))^2-r^2;
            if d<0
                l=l+1;
                flag(i,j)=2;
                
            end
        end
    end
end

%% Ghost points and its corresponding boundary points
l=0;
jjj=0;
for k=1:nrgrain
    for j=jmin(k):jmax(k)
        for i=imin(k):imax(k)
%% IB points and normal lines
            if (flag(i,j)==2 && (flag(i,j+1)==0 || flag(i,j-1)==0 || ...
                flag(i+1,j)==0 || flag(i-1,j)==0)) %point outside and at 
                                                  %least one solid neighbor  
                flag(i,j)=1;
                l=l+1;                      %number of ib points
                I_g(l)=i;
                J_g(l)=j;
                            
                X_g(l)=x(i);
                Y_g(l)=y(j);
                %find lines parameters between ib points and center (normal
                %line)
                a1=(yc(k)-Y_g(l))/(xc(k)-X_g(l));
                b=-a1*xc(k)+yc(k);
                    
            
            
                %% Ghost Boundary points: lines and circle eqn intersection 
                A=a1^2+1;
                B=2*a1*(b-yc(k))-2*xc(k);
                C=xc(k)^2+(b-yc(k))^2-r^2;
                delta=B^2-4*A*C;
                %which point is desired
                xla(1)=(-B+sqrt(delta))/(2*A);
                yla(1)=a1*xla(1)+b;
                xla(2)=(-B-sqrt(delta))/(2*A);
                yla(2)=a1*xla(2)+b;
                h2(1)=sqrt((xla(1)-X_g(l))^2+(yla(1)-Y_g(l))^2);
                h2(2)=sqrt((xla(2)-X_g(l))^2+(yla(2)-Y_g(l))^2);
                if (h2(1)>h2(2))
                    X_ib_g(l)=xla(2);
                    Y_ib_g(l)=yla(2);
                else %(d2(1)<d2(2))
                    X_ib_g(l)=xla(1);
                    Y_ib_g(l)=yla(1);
                end
                if (X_g(l)==xc(k)) && (Y_g(l)>yc(k))
                    X_ib_g(l)=xc(k);
                    Y_ib_g(l)=yc(k)+r;
                elseif (X_g(l)==xc(k)) && (Y_g(l)<yc(k))
                    X_ib_g(l)=xc(k);
                    Y_ib_g(l)=yc(k)-r;
                end
            end
        end
    end
end



%% Virtual points
% [landa_g_1,landa_g_2,landa_g_3,landa_g_4,...
%     A1_g,I_e,J_e,numg]=mirror_points(xc,yc,dx,dy,x,y,alpha,beta,q,...
%     flag,X_g,Y_g,X_ib_g,Y_ib_g,r);

[landa,landa_g_1,landa_g_2,landa_g_3,landa_g_4,A1_g,I_e,J_e,...
    I1,J1,I2,J2,I3,J3,I4,J4,X_e_g,Y_e_g,numg,nx_g,ny_g,...
    landa_g_5,landa_g_6,I5,J5,I6,J6]=...
    mirror_pointsBQ(xc,yc,x,y,alpha,beta,q,X_g,Y_g,X_ib_g,Y_ib_g,r,BQ);



% if i==47
%             X_e_g(47)
%             Y_e_g(4)
%             
%             X_ib_g(i)
%             Y_ib_g(i)
% end



l=0;
for k=1:nrgrain
    for i=imin(k):imax(k)
        for j=jmin(k):jmax(k)
            d=(x(i)-xc(k))^2+(y(j)-yc(k))^2-r^2;
            if d<0 && flag(i,j)==2
                l=l+1;
                I_solid(l)=i;
                J_solid(l)=j;
            end
        end
    end
end










