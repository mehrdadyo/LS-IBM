function [StateVar,VARIABLES,DOMAIN,BC,IBM,LS]=SetUpVariables
%% =========================== Domain size ================================

diamcyl = 1;                  % Diameter of cylinders

nrgrainx = 1;                 % Number of Cylinders in one row
nrgrainy = 1;

S = 2*diamcyl;                % Space between centers of objects

freeEast = 5*diamcyl;        %exit length after cylinders
freeWest = 5*diamcyl;      % entrance length after
freeNorth = 5*diamcyl;
freeSouth = 5*diamcyl;

lx = freeWest + nrgrainx*diamcyl + (nrgrainx-1)*S + ...
    freeEast; % Length in x direction
ly = freeSouth + nrgrainy*diamcyl + (nrgrainy-1)*S + ...
    freeNorth; % Length in x direction

% ly = 7*diamcyl;             % Length of computational domain in y dir





%% ======================== Additional variables ==========================

Re = 1;
uinflow = 1;
Pe = 1;
D = 1/Pe;
phi_inlet = 1;
phi_init = 1;
%% =============================== BCs ====================================

%% velocities and flow



% -alpha(dphi/dn_w)-beta phi=q
 
q = 0;
alpha = 0;
beta = 1;

BQu = 0;
BQv = 0;

BC_e_u = 3;
BC_w_u = 1;
BC_n_u = 3;         %BC_n=1 ==> Dirichlet,   BC_n=3 ==> Neumann
BC_s_u = 3;         %BC_s=1 ==> Dirichlet,   BC_s=3 ==> Neumann


BC_n_v = 3;         %BC_n=1 ==> Dirichlet,   BC_n=3 ==> Neumann
BC_s_v = 3;         %BC_s=1 ==> Dirichlet,   BC_s=3 ==> Neumann
BC_e_v = 3;
BC_w_v = 1;

BC_n_p = 3;         %BC_n=1 ==> Dirichlet,   BC_n=3 ==> Neumann
BC_s_p = 3;         %BC_s=1 ==> Dirichlet,   BC_s=3 ==> Neumann
BC_e_p = 3;
BC_w_p = 3;

P0_e = 0;
%% Transport
% -alpha(dphi/dn_w)-beta phi=q
q_phi = 0;
alpha_phi = D;
beta_phi = 1;
BQp = 0;

BC_e_phi=1;
BC_w_phi=3;
BC_n_phi=3;         %BC_n=1 ==> Dirichlet,   BC_n=3 ==> Neumann
BC_s_phi=3;         %BC_s=1 ==> Dirichlet,   BC_s=3 ==> Neumann



%% ========================= Objects ======================================

% yc=ly/2;
xcent = freeWest+diamcyl/2:2+diamcyl:lx-freeEast-diamcyl/2 ;
ycent = freeSouth+diamcyl/2:2+diamcyl:ly-freeNorth-diamcyl/2;
% for i=1:nrgrainx
%     xc(i)=entlength+(i-1)*S;
% end
phi_inside_u=0;                           % concen inside of objects
phi_inside_phi=q_phi;
[xc,yc] =deal(zeros(nrgrainx,nrgrainy));
for i=1:nrgrainx
    for j=1:nrgrainy
        xc(i,j) = xcent(i);
        yc(i,j) = ycent(j);
    end
end
        




%% ======================== Coordinates ===================================
uniform = 1;
dvdxdy= 81;
A=4;
expon=0;
r=1.00941835135602;
[DOMAIN]=Coordinates(xc,yc,lx,ly,A,diamcyl,dvdxdy,expon,r,uniform);

%% ========================= Time Steps ===================================

dt_man=5e-3;                             %% Arbitrary time step
dt_diff =  Re/((1/min(min(DOMAIN.dxu)))^2 + (1/min(min(DOMAIN.dyu)))^2); %% Diffisiun time step                  
dt_cour = 1/sqrt(2*(uinflow^2)*...
    ((1/min(min(DOMAIN.dxu)))^2 + (1/min(min((DOMAIN.dyu)))^2))); % Courant time
Dt=[dt_man,dt_diff,dt_cour];
% dt=min(Dt);
dt=dt_man;
% dt=0.008;
%% =============== Define marix variables / storage =======================
imax=DOMAIN.imax;
jmax=DOMAIN.jmax;
%U- velocity ---------------------
U=zeros(imax,jmax+1);

U_a=zeros(1,jmax+1);
U_b=zeros(1,jmax+1);
U_c=zeros(imax,1);
U_d=zeros(imax,1);

%V- velocity ----------------------
V=zeros(imax+1,jmax);

V_a=zeros(1,jmax);
V_b=zeros(1,jmax);
V_c=zeros(1,imax+1);
V_d=zeros(1,imax+1);


%Scalars -------------------------
P=zeros(imax+1,jmax+1);
phi=zeros(imax+1,jmax+1);

phi_a=zeros(1,jmax+1);
phi_b=zeros(1,jmax+1);
phi_c=zeros(1,jmax+1);
phi_d=zeros(1,jmax+1);

%% ================= Initialization of matrix fields ======================


U(:,:)=uinflow;
U_a(:)=U(1,:);
U_c(:)=uinflow;
U_d(:)=uinflow;





% V(:,:)=0.1*uinflow;



%% =================== initialize LS
h = 4.5;

LSCase.case = 1;  % 1 ==> circles obstacles, 2 ==> fracture
LSCase.xc = xc;
LSCase.yc = yc;
LSCase.diamcyl = diamcyl;
LSCase.h = h;


psi = LSInitialize(DOMAIN, LSCase);

% psi = zeros(imax+1,jmax+1);
% for i=1:imax+1
%     for j=1:jmax+1
%         d = sqrt( (DOMAIN.xp(i)-xc).^2 + (DOMAIN.yp(j)-yc).^2 )-diamcyl/2;
%         center = min(min(abs(d)));
%         [I,J,~] = find(abs(d) == center);
%         psi(i,j) = d(I(1),J(1));
%         
%     end
% end
psiU = interp2(DOMAIN.yp.*ones(imax+1,1), ones(jmax+1,1)'.*DOMAIN.xp', ...
    psi ,DOMAIN.yu.*ones(imax,1), ones(jmax+1,1)'.*DOMAIN.xu');

psiV = interp2(DOMAIN.yp.*ones(imax+1,1), ones(jmax+1,1)'.*DOMAIN.xp', ...
    psi ,DOMAIN.yv.*ones(imax+1,1), ones(jmax,1)'.*DOMAIN.xv');

U = double(psiU>0).*U;
V = double(psiV>0).*V;

%% initialize phi

phi_a(1,:)=phi_inlet;
phi(:,:)=phi_init;
phi = double(psi>0).*phi;

phi_old=phi;



%% LS solver variables

fac = 0.5;    % dta factor 
dtau = fac*min(min(DOMAIN.dxp));  % psedue time step

n_iter_ReLS = 4;

TimeSchemeLS = "RK3";
TimeSchemeRLS = "RK3";

% psi = sqrt( (DOMAIN.Xp-xc).^2 + (DOMAIN.Yp-yc).^2 )-diamcyl/2;
% LS = double(d>0)*d + double(d<0)*(d-r);

%% ========================== Time integration ============================

alpha_u=1;
alpha_v=alpha_u;
alpha_p=1;

VARIABLES = struct('D',D,'Re',Re,...
    'alpha_u',alpha_u,'alpha_v',alpha_v,...
    'alpha_p',alpha_p,'dt',dt,'dtVec',Dt,'Da',beta_phi,'Pe',Pe,...
    'n_iter_ReLS',n_iter_ReLS,'dtau',dtau,'TimeSchemeLS',TimeSchemeLS,...
    'TimeSchemeRLS',TimeSchemeRLS);


BC = struct('U_a',U_a,'U_b',U_b,'U_c',U_c,'U_d',U_d,...
    'V_a',V_a,'V_b',V_b,'V_c',V_c,'V_d',V_d,...
    'phi_a',phi_a,'phi_b',phi_b,'phi_c',phi_c,'phi_d',phi_d,...
    'BC_e_u',BC_e_u,'BC_w_u',BC_w_u,'BC_n_u',BC_n_u,'BC_s_u',BC_s_u,...
    'BC_e_v',BC_e_v,'BC_w_v',BC_w_v,'BC_n_v',BC_n_v,'BC_s_v',BC_s_v,...
    'BC_e_phi',BC_e_phi,'BC_w_phi',BC_w_phi,'BC_n_phi',BC_n_phi,...
    'BC_s_phi',BC_s_phi, ...
    'BC_n_p', BC_n_p, 'BC_s_p',BC_s_p, 'BC_w_p',BC_w_p, 'BC_e_p',BC_e_p, ...
    'P0_e', P0_e);

IBM = struct('q',q,'alpha',alpha,'beta',beta,'q_phi',q_phi,'alpha_phi',...
    alpha_phi,'beta_phi',beta_phi,...
    'phi_inside_u',phi_inside_u,'phi_inside_phi',phi_inside_phi,...
    'xc',xc,'yc',yc,'diamcyl',diamcyl,'nrgrainx',nrgrainx,...
    'nrgrainy',nrgrainy,'BQu',BQu,'BQv',BQv,'BQp',BQp);

StateVar = struct('U',U,'V',V,'P',P,'phi_old',phi_old);

LS = struct('psi',psi);
[LS] = LSnormals(LS,DOMAIN);
