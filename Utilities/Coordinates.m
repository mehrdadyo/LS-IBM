function DOMAIN=Coordinates(lx,ly,diamcyl,Grid)
% this function sets up the staggered grid points for u, v, p. we have a 
% refined area in the middle of the domain.



%% ---------- define first parameters --------------
% dvdxdy=31.25;

r = Grid.r;
expon = Grid.expon;
A = Grid.A;
dvdxdy = Grid.dvdxdy;
uniform = Grid.uniform;
lengthUnit = Grid.lengthUnit;

zoomedAreax = Grid.zoomedAreax;
zoomedAreay = Grid.zoomedAreay;

Lx_l = Grid.Lx_l;
Ly_b = Grid.Ly_b;


imax_fine = ceil(dvdxdy*lx/lengthUnit)+1;  % Grid size in x direction (imax=3 0*8+2-1)
jmax_fine = ceil(dvdxdy*ly/lengthUnit)+1;  % Grid size in y direction (jmax=18*8+2-1)
dx = lx/(imax_fine-1);    % Grid spacing in x direction
dy = ly/(jmax_fine-1);    % Grid spacing in y direction

% exponential growth
 n = ceil(log(A)/log(r));
 indx = 1:n;
 
 if ~uniform
    %% refined area  
    dx_r=dx/A;                      % refined grid size
    xuRef = Lx_l + (dx_r : dx_r : zoomedAreax * diamcyl);  
    %% refined area  
    dy_r=dy/A;                      % refined grid size
    yvRef = Ly_b + (dy_r : dy_r : zoomedAreay * diamcyl);

    if expon
        %% transional zone
        x_trans = ((1 - r.^indx)/(1 - r)) * dx_r;
        y_trans = ((1 - r.^indx)/(1 - r)) * dy_r;
        resX = Lx_l - x_trans(end);
        resY = Ly_b - y_trans(end);
        if resY<dx || resX<dx
            error('The transition grid exceeds the coarse grid area length. Please reduce the corase grid cell or A.')
        end
    elseif ~expon
        
        x_trans = [];
        y_trans = [];
        
        resX = Lx_l;
        resY = Ly_b;
    end
    
    %% coarse grid area
    %==================  Lx_l ====================================
    % ++++++++ xu +++++++++++++
    m = floor(resX/dx);
    indxCoarse = 0:1:m-1;

    dx_mod = resX/m;
    x_coarse_L = indxCoarse .*dx_mod;
    x_transL = Lx_l - [fliplr(x_trans) 0];
    xu = [x_coarse_L x_transL xuRef];


    %==================  Lx_r ====================================
    xu = [xu xu(end)+x_trans];
    res = lx - xu(end);
    m = floor(res/dx);
    dx_mod = res/m;
    indxCoarse = 1:m;
    x_coarse_R = xu(end) + indxCoarse .*dx_mod;
    xu = [xu x_coarse_R];

    % ++++++++ xv +++++++++++++
    xv = (xu(1:end-1)+xu(2:end))/2;
    xv = [-dx_mod/2 xv xv(end)+dx_mod];
    
    
    %==================  Ly_b ====================================
    % ++++++++ xu +++++++++++++
    
    m = floor(resY/dy);
    indxCoarse = 0:1:m-1;

    dy_mod = resY/m;
    y_coarse_B = indxCoarse .*dy_mod;
    y_transB = Ly_b - [fliplr(y_trans) 0];
    yv = [y_coarse_B y_transB yvRef];


    %==================  Ly_u ====================================
    yv = [yv yv(end)+y_trans];
    res = ly - yv(end);
    m = floor(res/dy);
    dy_mod = res/m;
    indxCoarse = 1:m;
    y_coarse_U = yv(end) + indxCoarse .*dy_mod;
    yv = [yv y_coarse_U];

    % ++++++++ xv +++++++++++++
    yu = (yv(1:end-1)+yv(2:end))/2;
    yu = [-dy_mod/2 yu yu(end)+dy_mod];



 elseif uniform 
    xu = 0:dx:lx;
    xv = -dx/2:dx:lx+dx/2;

    yu = -dy/2:dy:ly+dy/2;
    yv = 0:dy:ly;

     
     
     
 end
 
 


%% scalar nodes --------------------------
xp=xv;
yp=yu;

imax=length(xp)-1;
jmax=length(yp)-1;

% size dx_u=Nx-1 x Ny+1
dxu=diff(xu);
dxv=diff(xv);
dxp=diff(xp);



dyu=diff(yu);
dyv=diff(yv);
dyp=diff(yp);


CoEWu=(xp(2:end-1)-xu(1:end-1))./dxu;
CoNSu=(xu-xv(1:end-1))./dxv;
dV_u=[0 dxv(2:end)]'*[0 dyv 0];

CoEWv=(yv(1:end)-yp(1:end-1))./dyu;
CoNSv=(yp(2:end-1)-yv(1:end-1))./dyv;
dV_v=[0 dxu 0]'*[0 dyu(2:end-1) 0];

CoEWp=(xu(1:end)-xp(1:end-1))./dxp;
CoNSp=(yv(1:end)-yp(1:end-1))./dyp;
dV_p=[0 dxu 0]'*[0 dyv 0];



%% %%%%%%%%%%%%%%%%% QUICK non-uniform coefficients: %%%%%%%%%%%%%%%%%%%%%%
%% =========================== U-velocity =================================


%% east 

% positive flow
x=xp(3:end-2);
xU=xu(2:end-2);
xUU=xu(1:end-3);
xD=xu(3:end-1);
[g1u_e_p,g2u_e_p]=QuickInterp(x,xU,xUU,xD);
g1u_e_p=[0 g1u_e_p 0 0];
g2u_e_p=[0 g2u_e_p 0 0];
g1u_e_p=g1u_e_p'*ones(1,jmax+1);
g2u_e_p=g2u_e_p'*ones(1,jmax+1);

% negative flow
x=xp(3:end-2);
xU=xu(3:end-1);
xUU=xu(4:end);
xD=xu(2:end-2);
[g1u_e_n,g2u_e_n]=QuickInterp(x,xU,xUU,xD);
g1u_e_n=[0 g1u_e_n 0 0];
g2u_e_n=[0 g2u_e_n 0 0];
g1u_e_n=g1u_e_n'*ones(1,jmax+1);
g2u_e_n=g2u_e_n'*ones(1,jmax+1);

%% west 

% positive flow
x=xp(3:end-2);
xU=xu(2:end-2);
xUU=xu(1:end-3);
xD=xu(3:end-1);
[g1u_w_p,g2u_w_p]=QuickInterp(x,xU,xUU,xD);
g1u_w_p=[0 0 g1u_w_p 0];
g2u_w_p=[0 0 g2u_w_p 0];
g1u_w_p=g1u_w_p'*ones(1,jmax+1);
g2u_w_p=g2u_w_p'*ones(1,jmax+1);

% negative flow
x=xp(3:end-2);
xU=xu(3:end-1);
xUU=xu(4:end);
xD=xu(2:end-2);
[g1u_w_n,g2u_w_n]=QuickInterp(x,xU,xUU,xD);
g1u_w_n=[0 0 g1u_w_n 0];
g2u_w_n=[0 0 g2u_w_n 0];
g1u_w_n=g1u_w_n'*ones(1,jmax+1);
g2u_w_n=g2u_w_n'*ones(1,jmax+1);

%% north 

% positive flow
y=yv(2:end-1);
yU=yu(2:end-2);
yUU=yu(1:end-3);
yD=yu(3:end-1);
[g1u_n_p,g2u_n_p]=QuickInterp(y,yU,yUU,yD);
g1u_n_p=[0 g1u_n_p 0 0];
g2u_n_p=[0 g2u_n_p 0 0];
g1u_n_p=ones(imax,1)*g1u_n_p;
g2u_n_p=ones(imax,1)*g2u_n_p;

% negative flow
y=yv(2:end-1);
yU=yu(3:end-1);
yUU=yu(4:end);
yD=yu(2:end-2);
[g1u_n_n,g2u_n_n]=QuickInterp(y,yU,yUU,yD);
g1u_n_n=[0 g1u_n_n 0 0];
g2u_n_n=[0 g2u_n_n 0 0];
g1u_n_n=ones(imax,1)*g1u_n_n;
g2u_n_n=ones(imax,1)*g2u_n_n;

%% south 

% positive flow
y=yv(2:end-1);
yU=yu(2:end-2);
yUU=yu(1:end-3);
yD=yu(3:end-1);
[g1u_s_p,g2u_s_p]=QuickInterp(y,yU,yUU,yD);
g1u_s_p=[0 0 g1u_s_p 0];
g2u_s_p=[0 0 g2u_s_p 0];
g1u_s_p=ones(imax,1)*g1u_s_p;
g2u_s_p=ones(imax,1)*g2u_s_p;

% negative flow
y=yv(2:end-1);
yU=yu(3:end-1);
yUU=yu(4:end);
yD=yu(2:end-2);
[g1u_s_n,g2u_s_n]=QuickInterp(y,yU,yUU,yD);
g1u_s_n=[0 0 g1u_s_n 0];
g2u_s_n=[0 0 g2u_s_n 0];
g1u_s_n=ones(imax,1)*g1u_s_n;
g2u_s_n=ones(imax,1)*g2u_s_n;

%% =========================== V-velocity =================================


%% east 

% positive flow
x=xu(2:end-1);
xU=xv(2:end-2);
xUU=xv(1:end-3);
xD=xv(3:end-1);
[g1v_e_p,g2v_e_p]=QuickInterp(x,xU,xUU,xD);
g1v_e_p=[0 g1v_e_p 0 0];
g2v_e_p=[0 g2v_e_p 0 0];
g1v_e_p=g1v_e_p'*ones(1,jmax);
g2v_e_p=g2v_e_p'*ones(1,jmax);

% negative flow
x=xu(2:end-1);
xU=xv(3:end-1);
xUU=xv(4:end);
xD=xv(2:end-2);
[g1v_e_n,g2v_e_n]=QuickInterp(x,xU,xUU,xD);
g1v_e_n=[0 g1v_e_n 0 0];
g2v_e_n=[0 g2v_e_n 0 0];
g1v_e_n=g1v_e_n'*ones(1,jmax);
g2v_e_n=g2v_e_n'*ones(1,jmax);

%% west 

% positive flow
x=xu(2:end-1);
xU=xv(2:end-2);
xUU=xv(1:end-3);
xD=xv(3:end-1);
[g1v_w_p,g2v_w_p]=QuickInterp(x,xU,xUU,xD);
g1v_w_p=[0 0 g1v_w_p 0];
g2v_w_p=[0 0 g2v_w_p 0];
g1v_w_p=g1v_w_p'*ones(1,jmax);
g2v_w_p=g2v_w_p'*ones(1,jmax);

% negative flow
x=xu(2:end-1);
xU=xv(3:end-1);
xUU=xv(4:end);
xD=xv(2:end-2);
[g1v_w_n,g2v_w_n]=QuickInterp(x,xU,xUU,xD);
g1v_w_n=[0 0 g1v_w_n 0];
g2v_w_n=[0 0 g2v_w_n 0];
g1v_w_n=g1v_w_n'*ones(1,jmax);
g2v_w_n=g2v_w_n'*ones(1,jmax);

%% north 

% positive flow
y=yp(3:end-2);
yU=yv(2:end-2);
yUU=yv(1:end-3);
yD=yv(3:end-1);
[g1v_n_p,g2v_n_p]=QuickInterp(y,yU,yUU,yD);
g1v_n_p=[0 g1v_n_p 0 0];
g2v_n_p=[0 g2v_n_p 0 0];
g1v_n_p=ones(imax+1,1)*g1v_n_p;
g2v_n_p=ones(imax+1,1)*g2v_n_p;

% negative flow
y=yp(3:end-2);
yU=yv(3:end-1);
yUU=yv(4:end);
yD=yv(2:end-2);
[g1v_n_n,g2v_n_n]=QuickInterp(y,yU,yUU,yD);
g1v_n_n=[0 g1v_n_n 0 0];
g2v_n_n=[0 g2v_n_n 0 0];
g1v_n_n=ones(imax+1,1)*g1v_n_n;
g2v_n_n=ones(imax+1,1)*g2v_n_n;

%% south 

% positive flow
y=yp(3:end-2);
yU=yv(2:end-2);
yUU=yv(1:end-3);
yD=yv(3:end-1);
[g1v_s_p,g2v_s_p]=QuickInterp(y,yU,yUU,yD);
g1v_s_p=[0 0 g1v_s_p 0];
g2v_s_p=[0 0 g2v_s_p 0];
g1v_s_p=ones(imax+1,1)*g1v_s_p;
g2v_s_p=ones(imax+1,1)*g2v_s_p;

% negative flow
y=yp(3:end-2);
yU=yv(3:end-1);
yUU=yv(4:end);
yD=yv(2:end-2);
[g1v_s_n,g2v_s_n]=QuickInterp(y,yU,yUU,yD);
g1v_s_n=[0 0 g1v_s_n 0];
g2v_s_n=[0 0 g2v_s_n 0];
g1v_s_n=ones(imax+1,1)*g1v_s_n;
g2v_s_n=ones(imax+1,1)*g2v_s_n;


%% =========================== Transport ==================================


%% east 

% positive flow
x=xu(2:end-1);
xU=xp(2:end-2);
xUU=xp(1:end-3);
xD=xp(3:end-1);
[g1c_e_p,g2c_e_p]=QuickInterp(x,xU,xUU,xD);
g1c_e_p=[0 g1c_e_p 0 0];
g2c_e_p=[0 g2c_e_p 0 0];
g1c_e_p=g1c_e_p'*ones(1,jmax+1);
g2c_e_p=g2c_e_p'*ones(1,jmax+1);

% negative flow
x=xu(2:end-1);
xU=xp(3:end-1);
xUU=xp(4:end);
xD=xp(2:end-2);
[g1c_e_n,g2c_e_n]=QuickInterp(x,xU,xUU,xD);
g1c_e_n=[0 g1c_e_n 0 0];
g2c_e_n=[0 g2c_e_n 0 0];
g1c_e_n=g1c_e_n'*ones(1,jmax+1);
g2c_e_n=g2c_e_n'*ones(1,jmax+1);

%% west 

% positive flow
x=xu(2:end-1);
xU=xp(2:end-2);
xUU=xp(1:end-3);
xD=xp(3:end-1);
[g1c_w_p,g2c_w_p]=QuickInterp(x,xU,xUU,xD);
g1c_w_p=[0 0 g1c_w_p 0];
g2c_w_p=[0 0 g2c_w_p 0];
g1c_w_p=g1c_w_p'*ones(1,jmax+1);
g2c_w_p=g2c_w_p'*ones(1,jmax+1);

% negative flow
x=xu(2:end-1);
xU=xp(3:end-1);
xUU=xp(4:end);
xD=xp(2:end-2);
[g1c_w_n,g2c_w_n]=QuickInterp(x,xU,xUU,xD);
g1c_w_n=[0 0 g1c_w_n 0];
g2c_w_n=[0 0 g2c_w_n 0];
g1c_w_n=g1c_w_n'*ones(1,jmax+1);
g2c_w_n=g2c_w_n'*ones(1,jmax+1);

%% north 

% positive flow
y=yv(2:end-1);
yU=yp(2:end-2);
yUU=yp(1:end-3);
yD=yp(3:end-1);
[g1c_n_p,g2c_n_p]=QuickInterp(y,yU,yUU,yD);
g1c_n_p=[0 g1c_n_p 0 0];
g2c_n_p=[0 g2c_n_p 0 0];
g1c_n_p=ones(imax+1,1)*g1c_n_p;
g2c_n_p=ones(imax+1,1)*g2c_n_p;

% negative flow
y=yv(2:end-1);
yU=yp(3:end-1);
yUU=yp(4:end);
yD=yp(2:end-2);
[g1c_n_n,g2c_n_n]=QuickInterp(y,yU,yUU,yD);
g1c_n_n=[0 g1c_n_n 0 0];
g2c_n_n=[0 g2c_n_n 0 0];
g1c_n_n=ones(imax+1,1)*g1c_n_n;
g2c_n_n=ones(imax+1,1)*g2c_n_n;

%% south 

% positive flow
y=yv(2:end-1);
yU=yp(2:end-2);
yUU=yp(1:end-3);
yD=yp(3:end-1);
[g1c_s_p,g2c_s_p]=QuickInterp(y,yU,yUU,yD);
g1c_s_p=[0 0 g1c_s_p 0];
g2c_s_p=[0 0 g2c_s_p 0];
g1c_s_p=ones(imax+1,1)*g1c_s_p;
g2c_s_p=ones(imax+1,1)*g2c_s_p;

% negative flow
y=yv(2:end-1);
yU=yp(3:end-1);
yUU=yp(4:end);
yD=yp(2:end-2);
[g1c_s_n,g2c_s_n]=QuickInterp(y,yU,yUU,yD);
g1c_s_n=[0 0 g1c_s_n 0];
g2c_s_n=[0 0 g2c_s_n 0];
g1c_s_n=ones(imax+1,1)*g1c_s_n;
g2c_s_n=ones(imax+1,1)*g2c_s_n;

[Xu,Yu] = meshgrid(xu,yu);
Xu = Xu';
Yu = Yu';


[Xv,Yv] = meshgrid(xv,yv);
Xv = Xv';
Yv = Yv';

[Xp,Yp] = meshgrid(xp,yp);
Xp = Xp';
Yp = Yp';



%% ========================================================================

%% %%%%%%%%%%%%%% CENTRAL DIFFERENCE INTERPOLATION COEFFICIENTS %%%%%%%%%%%

% size dx_u=Nx-1 x Ny+1
dxu=dxu'*ones(1,jmax+1);
CoEWu=CoEWu'*ones(1,jmax+1);
% size dx_v=Nx x Ny
dxv=dxv'*ones(1,jmax);
CoNSu=CoNSu'*ones(1,jmax);
% size dx_u=Nx x Ny+1
dxp=dxp'*ones(1,jmax+1);
CoEWp=CoEWp'*ones(1,jmax+1);

% size dy_u=Nx x Ny
dyu=ones(imax,1)*dyu;
CoEWv=ones(imax,1)*CoEWv;
% size dy_v=Nx+1 x Ny-1
dyv=ones(imax+1,1)*dyv;
CoNSv=ones(imax+1,1)*CoNSv;
% size dy_p=Nx+1 x Ny
dyp=ones(imax+1,1)*dyp;
CoNSp=ones(imax+1,1)*CoNSp;

DOMAIN = struct('xu',xu,'yu',yu,'xv',xv,'yv',yv,'xp',xp,'yp',yp,...
    'imax',imax,'jmax',jmax,'dxu',dxu,'dyu',dyu,'dxv',dxv,'dyv',dyv,...
    'dxp',dxp,'dyp',dyp,'lx',lx,'ly',ly,...
    'dV_u',dV_u,'dV_v',dV_v,'dV_p',dV_p,...
    'CoEWu',CoEWu,'CoNSu',CoNSu,'CoEWv',CoEWv,'CoNSv',CoNSv,...
    'CoEWp',CoEWp,'CoNSp',CoNSp,...
    'g1c_e_p',g1c_e_p,'g2c_e_p',g2c_e_p,'g1c_w_p',g1c_w_p,...
    'g2c_w_p',g2c_w_p,'g1c_n_p',g1c_n_p,'g2c_n_p',g2c_n_p,...
    'g1c_s_p',g1c_s_p,'g2c_s_p',g2c_s_p,'g1c_e_n',g1c_e_n,...
    'g2c_e_n',g2c_e_n,'g1c_w_n',g1c_w_n,'g2c_w_n',g2c_w_n,...
    'g1c_n_n',g1c_n_n,'g2c_n_n',g2c_n_n,'g1c_s_n',g1c_s_n,...
    'g2c_s_n',g2c_s_n,...
    'g1u_e_p',g1u_e_p,'g2u_e_p',g2u_e_p,'g1u_w_p',g1u_w_p,...
    'g2u_w_p',g2u_w_p,'g1u_n_p',g1u_n_p,'g2u_n_p',g2u_n_p,...
    'g1u_s_p',g1u_s_p,'g2u_s_p',g2u_s_p,'g1u_e_n',g1u_e_n,...
    'g2u_e_n',g2u_e_n,'g1u_w_n',g1u_w_n,'g2u_w_n',g2u_w_n,...
    'g1u_n_n',g1u_n_n,'g2u_n_n',g2u_n_n,'g1u_s_n',g1u_s_n,...
    'g2u_s_n',g2u_s_n,...
    'g1v_e_p',g1v_e_p,'g2v_e_p',g2v_e_p,'g1v_w_p',g1v_w_p,...
    'g2v_w_p',g2v_w_p,'g1v_n_p',g1v_n_p,'g2v_n_p',g2v_n_p,...
    'g1v_s_p',g1v_s_p,'g2v_s_p',g2v_s_p,'g1v_e_n',g1v_e_n,...
    'g2v_e_n',g2v_e_n,'g1v_w_n',g1v_w_n,'g2v_w_n',g2v_w_n,...
    'g1v_n_n',g1v_n_n,'g2v_n_n',g2v_n_n,'g1v_s_n',g1v_s_n,...
    'g2v_s_n',g2v_s_n,'Xu',Xu,'Yu',Yu,'Xv',Xv,'Yv',Yv,'Xp',Xp,'Yp',Yp);





