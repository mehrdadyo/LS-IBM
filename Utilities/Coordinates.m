function DOMAIN=Coordinates(xc,yc,lx,ly,A,diamcyl,dvdxdy,expon,r,uniform)
% this function sets up the staggered grid points for u, v, p. we have a 
% refined area in the middle of the domain.



%% ---------- define first parameters --------------
% dvdxdy=31.25;
imax_fine = dvdxdy*lx+1;  % Grid size in x direction (imax=3 0*8+2-1)
jmax_fine = dvdxdy*ly+1;  % Grid size in y direction (jmax=18*8+2-1)
dx = lx/(imax_fine-1);    % Grid spacing in x direction
dy = ly/(jmax_fine-1);    % Grid spacing in y direction

% exponential growth
 if expon==1
    ind1=0:1:250;
    dxy=1/125*r.^ind1;
    dxyrev=fliplr(dxy);
 end


%% ---------- refined area--------------
if uniform == 0
    Lx_r=4*diamcyl;                 % horizontal length of fine area
    dx_r=dx/A;                      % refined grid size
    x_r=(dx_r:dx_r:Lx_r);           % refined grid size

    Ly_r=4*diamcyl;                 % vertical length of fine area 
    dy_r=dy/A;                      % refined grid size 
    y_r=(dy_r/2:dy_r:Ly_r-dy_r/2);  % refined grid size     

    % Left
    Lx_l=xc-2*diamcyl;
    % Right starting
    Lx_r=xc+2*diamcyl;
    % Bottom
    Ly_b=yc-2*diamcyl;             % vertical length of coarse area at bottom
    % Top starting
    Ly_t=yc+2*diamcyl;

    xu_2=x_r+Lx_l;
    xv_2=x_r+Lx_l-dx_r/2;


    yv_2=y_r+Ly_b+dy_r/2;
    yu_2=y_r+Ly_b;

    %% ---------- coarsed area--------------

    if expon==0

    %% U nodes ---------------------------
        xu_1=(0:dx:Lx_l);
        xu_3=(Lx_r+dx:dx:lx);
        xu=[xu_1 xu_2 xu_3];         

        yu_1=(-dy/2:dy:Ly_b-dy/2);
        yu_3=(Ly_t+dy/2:dy:ly+dy/2);
        yu=[yu_1 yu_2 yu_3];

    %% V nodes -----------------------------
        xv_1=(-dx/2:dx:Lx_l-dx/2);
        xv_3=(Lx_r+dx/2:dx:lx+dx/2);
        xv=[xv_1 xv_2 xv_3];         

        yv_1=(0:dy:Ly_b);
        yv_3=(Ly_t+dy:dy:ly);
        yv=[yv_1 yv_2 yv_3];
    elseif expon==1

    %% xU nodes ---------------------------
        xu_1=zeros(1,length(dxy));
        xu_3=zeros(1,length(dxy)-1);

        for i=2:length(dxyrev)
            xu_1(i)=xu_1(i-1)+dxyrev(i);
        end

        for i=2:length(dxy)
            xu_3(i)=xu_3(i-1)+dxy(i-1);
        end
        xu_3=xu_3(2:end)+xu_2(end);
        lx_rem=lx-xu_3(end);
        i_rem=floor(lx_rem/(xu_3(end)-xu_3(end-1)));
        xu_4=xu_3(end)+(1:1:i_rem)*(xu_3(end)-xu_3(end-1));
        xu=[xu_1 xu_2 xu_3 xu_4];




    %% yV nodes -----------------------------
        yv_1=zeros(1,length(dxy));
        yv_3=zeros(1,length(dxy)-1);

        for i=2:length(dxyrev)
            yv_1(i)=yv_1(i-1)+dxyrev(i);
        end

        for i=2:length(dxy)
            yv_3(i)=yv_3(i-1)+dxy(i-1);
        end
        yv_3=yv_3(2:end)+yv_2(end);
        yv=[yv_1 yv_2 yv_3];

    %% yU and xV
        yu=1/2*(yv(1:end-1)+yv(2:end));
        yu=[2*yv(1)-yu(1) yu 2*yv(end)-yu(end)];

        xv=1/2*(xu(1:end-1)+xu(2:end));
        xv=[2*xu(1)-xv(1) xv 2*xu(end)-xv(end)];

    end    
elseif uniform == 1
    xu = 0:dx:lx;
    xv = -dx/2:dx:lx+dx/2;
    
    yu = -dy/2:dy:ly+dy/2;
    yv = 0:dy:ly;
end
%% P nodes --------------------------
xp=xv;
yp=yu;

imax=length(xp)-1;
jmax=length(yp)-1;

% size dx_u=Nx-1 x Ny+1
dxu=xu(2:end)-xu(1:end-1);
dxv=xv(2:end)-xv(1:end-1);
dxp=xp(2:end)-xp(1:end-1);



dyu=yu(2:end)-yu(1:end-1);
dyv=yv(2:end)-yv(1:end-1);
dyp=yp(2:end)-yp(1:end-1);


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





