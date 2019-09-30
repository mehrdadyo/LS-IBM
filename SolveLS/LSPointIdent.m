function [landa,landa_g_1,landa_g_2,landa_g_3,landa_g_4,A1_g,numg,flag,...
    I_g,J_g,I_solid,J_solid,I_e,J_e,I1,J1,I2,J2,I3,J3,I4,J4,...
    landa_g_5,landa_g_6,I5,J5,I6,J6,nx,ny]=...
                 LSPointIdent(DOMAIN,alpha,beta,q,BQ,LS,UVP)

% -alpha dphi - beta phi = q 
%size of coordinate matrix n=size[x], m=size[y]






PSI = LS.psi;




% flag=-1 ==>extrapolation cells (fluid)
% flag=1  ==>Mirror or ghost cells (soild)
% flag=0  ==> all others (fluid)
% flag=2  ==> body cells (solid)
% first everywhere fluid (flag=-1), second inside points solid (flag=1),
% third ib points (flag=-1), fourth all other points zero.
% numg=0;  
% numg=0;

X = DOMAIN.Xp;
Y = DOMAIN.Yp;

%% find out which variable is to be reconstructed by looking at UVP

if UVP == -1  % u-velocity
    
    x = DOMAIN.xu;
    y = DOMAIN.yu;
    Xq = DOMAIN.Xu;
    Yq = DOMAIN.Yu;
    
%     psi = interp2(X',Y',PSI',Xq',Yq')';  

    psi = zeros(length(DOMAIN.xu),length(DOMAIN.yu));
    psi(:,:) = ( PSI(1:end-1,:)+ PSI(2:end,:) )/2;

    
    dx=min(min(DOMAIN.dxu));
    dy=min(min(DOMAIN.dyu));

    
elseif UVP == 0 % v-velocity
    
    x = DOMAIN.xv;
    y = DOMAIN.yv;
    Xq = DOMAIN.Xv;
    Yq = DOMAIN.Yv;
    
%     psi = interp2(X',Y',PSI',Xq',Yq')';
    psi = zeros(length(DOMAIN.xv),length(DOMAIN.yv));
    psi(:,:) = ( PSI(:,1:end-1)+PSI(:,2:end) )/2;
    
    dx=min(min(DOMAIN.dxv));
    dy=min(min(DOMAIN.dyv));
    
    
elseif UVP == 1 % scalar variables
    
    x = DOMAIN.xp;
    y = DOMAIN.yp;
    
    psi = PSI;
    
    dx=min(min(DOMAIN.dxp));
    dy=min(min(DOMAIN.dyp));
    
    
end

Nx = size(x,2);
Ny = size(y,2);
% flag for all points
flag = zeros(Nx,Ny);

%% solid points
flag = 2*double(psi<0) + flag;

%% Ghost points and its corresponding boundary points

[NBabsSolid,NBabs,absNB] = deal(zeros(Nx,Ny));

% NBabsSolid (2:end-1,2:end-1) = flag(3:end,2:end-1) + ...
%     flag(1:end-2,2:end-1) + flag(2:end-1,3:end) + ...
%     flag(2:end-1,1:end-2);
% 
% SolidCorrection = double((flag == 0) .* (NBabsSolid == 8));
% FluidCorrection = double((flag == 2) .* (NBabsSolid == 0));
% 
% 
% flag = flag + 2*SolidCorrection -2*FluidCorrection;
% 

NBabs (2:end-1,2:end-1) =  flag(3:end,2:end-1) + ...
    flag(1:end-2,2:end-1) + flag(2:end-1,3:end) + ...
    flag(2:end-1,1:end-2);   % neighbors 


% absNB (2:end-1,2:end-1) =  abs(psi(3:end,2:end-1)+psi(1:end-2,2:end-1) +...
%     psi(2:end-1,3:end) + psi(2:end-1,1:end-2));

% flag = flag - double(NBabs ~= absNB) .* double(flag == 2);
flag (2:end-1,2:end-1) = flag (2:end-1,2:end-1)...
    - double(NBabs (2:end-1,2:end-1) ~= 8) .* double(flag (2:end-1,2:end-1) == 2);



[I_g,J_g,~] = find(flag == 1);
I_g = I_g';
J_g = J_g';
X_g = x(I_g);
Y_g = y(J_g);




%% Virtual points

[landa,landa_g_1,landa_g_2,landa_g_3,landa_g_4,A1_g,I_e,J_e,...
    I1,J1,I2,J2,I3,J3,I4,J4,numg,landa_g_5,landa_g_6,I5,J5,I6,J6,nx,ny] = ...
    LSmirPointsBQ(x,y,alpha,beta,q,X_g,Y_g,BQ,dx,dy,I_g,J_g,LS);

%% Solid Points indices
[J_solid,I_solid] = find(flag' == 2);

J_solid = J_solid';
I_solid = I_solid';









