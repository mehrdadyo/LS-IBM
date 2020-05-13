function [StateVar,VARIABLES,DOMAIN,BC,IBM,LS, LSCase]= setUpVariablesNonDim
%% ===================== geometry type ====================================
% possible options : 'grain', 'fracture'

geometry = 'grain';

lengthUnit = 1;
%% =========================== Domain size ================================
if strcmp(geometry, 'grain')

%%% ======= circular grain ===============
    LSCase.case = 1;   % 1 ==> circles 
                    
    
    diamcyl = 2*lengthUnit;                  % Diameter of cylinders in cm

    nrgrainx = 1;                 % Number of Cylinders in one row
    nrgrainy = 1;

    S = 2*diamcyl;                % Space between centers of objects

    freeEast = 2*diamcyl;        %exit length after cylinders
    freeWest = 2*diamcyl;      % entrance length after
    freeNorth = 0.75*diamcyl;
    freeSouth = 0.75*diamcyl;

    lx = freeWest + nrgrainx*diamcyl + (nrgrainx-1)*S + ...
        freeEast; % Length in x direction
    ly = freeSouth + nrgrainy*diamcyl + (nrgrainy-1)*S + ...
        freeNorth; % Length in x direction

    
elseif strcmp(geometry, 'fracture')
    
%%% ======= rough fracture ===============
    LSCase.case = 3;    % 2 ==> flat fracture
                        % 3 ==> uni-mineral rough fracture
                        % 4 ==>  bi-mineral rough fracture
                   
    b = 0.2;        % half of average fracture width
    lambda = b;     % wave length
    a = 0.2*b;

    lx = 12 * b + 2 * b;
    x_0 = b;
    
    ly = 4 * b;
    
    diamcyl = [];
    xc = [];
    yc = [];
    nrgrainx = [];
    nrgrainy = [];
    
    
end




%% ======================== Additional variables ==========================
% - Re is in fact kinematic viscosity
% - Pe is 1/(diffusion coeffcient)
Re = 0.6;      % Re is the Kinematic viscosity in dimensional case
uinflow = 1;
Pe = 600;           
D = 1/Pe;
phi_inlet = 1;
phi_init = 0;
density = 1;
dimensional = false;
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
BC_n_u = 1;         %BC_n=1 ==> Dirichlet,   BC_n=3 ==> Neumann
BC_s_u = 1;         %BC_s=1 ==> Dirichlet,   BC_s=3 ==> Neumann


BC_n_v = 1;         %BC_n=1 ==> Dirichlet,   BC_n=3 ==> Neumann
BC_s_v = 1;         %BC_s=1 ==> Dirichlet,   BC_s=3 ==> Neumann
BC_e_v = 3;
BC_w_v = 1;

BC_n_p = 3;         %BC_n=1 ==> Dirichlet,   BC_n=3 ==> Neumann
BC_s_p = 3;         %BC_s=1 ==> Dirichlet,   BC_s=3 ==> Neumann
BC_e_p = 3;
BC_w_p = 3;

P0_e = 0;
%% Transport
% -alpha(dphi/dn_w)-beta phi=q
uniMineral = true;
BQp = 0;

if uniMineral 
    q_phi = 0;
    alpha_phi = -1;
    beta_phi = 178;     % beta is K in the reaction
    dissolution = true;
else
    q_phi.q_phi1 = 0;
    alpha_phi.alpha_phi1 = 1;
    beta_phi.beta_phi1 = -0.244;     % beta is K in the reaction
    dissolution.dissolution1 = true;  

    q_phi.q_phi2 = 0;
    alpha_phi.alpha_phi2 = 1;
    beta_phi.beta_phi2 = 0;     % beta is K in the reaction
    dissolution.dissolution2 = true;  
end
    

BC_e_phi=1;
BC_w_phi=3;
BC_n_phi=3;         %BC_n=1 ==> Dirichlet,   BC_n=3 ==> Neumann
BC_s_phi=3;         %BC_s=1 ==> Dirichlet,   BC_s=3 ==> Neumann



%% ========================= Objects ======================================

% ==== values of velocity and scalar inside the solids 
phi_inside_u=0;                           % concen inside of objects
phi_inside_phi=0;


% ==== define grain centers
if strcmp(geometry, 'grain')

    xcent = freeWest+diamcyl/2:2+diamcyl:lx-freeEast-diamcyl/2 ;
    ycent = freeSouth+diamcyl/2:2+diamcyl:ly-freeNorth-diamcyl/2;
    % for i=1:nrgrainx
    %     xc(i)=entlength+(i-1)*S;
    % end
    [xc,yc] =deal(zeros(nrgrainx,nrgrainy));
    for i=1:nrgrainx
        for j=1:nrgrainy
            xc(i,j) = xcent(i);
            yc(i,j) = ycent(j);
        end
    end
        
end
% xc = [];
% yc = [];
% 

%% ======================== Coordinates ===================================
Grid.r = 1.1;
Grid.expon = 1;
Grid.A = 2;
Grid.dvdxdy = 32;
Grid.uniform = true;
Grid.zoomedAreax = 1.5;
Grid.zoomedAreay = 1.5;
Grid.Lx_l = 1.5/50;
Grid.Ly_b = 0.5/50;
Grid.lengthUnit = lengthUnit;

[DOMAIN]=Coordinates(lx,ly,diamcyl,Grid);
%% ========================= Time Steps ===================================

dt_man=1e-2;                             %% Arbitrary time step
dt_diff =  Re/((1/min(min(DOMAIN.dxu)))^2 + (1/min(min(DOMAIN.dyu)))^2); %% Diffisiun time step                  
dt_cour = 1/sqrt(2*(uinflow^2)*...
    ((1/min(min(DOMAIN.dxu)))^2 + (1/min(min((DOMAIN.dyu)))^2))); % Courant time
Dt=[dt_man,dt_diff,dt_cour];
% dt=min(Dt);
dt=dt_man;
% dt=0.008;
if dt > 0.25*min(min(DOMAIN.dxp))/(beta_phi/Pe)
    disp('Warning!!! The time step seems too high for the the reaction rate')
    disp('LS will move too fast, please consider reducing time step')
    pause(1)
end
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
V_c=zeros(imax+1,1);
V_d=zeros(imax+1,1);


%Scalars -------------------------
P=zeros(imax+1,jmax+1);
phi=zeros(imax+1,jmax+1);

phi_a=zeros(1,jmax+1);
phi_b=zeros(1,jmax+1);
phi_c=zeros(1,jmax+1);
phi_d=zeros(1,jmax+1);

%% =================== LS parameters
symmetry = false;  % symmetry BC not currently working, 
                   % set FALSE until fixed
                   
% LSCase.case = 3;   % 1 ==> circles 
%                    % 2 ==> flat fracture
%                    % 3 ==> uni-mineral rough fracture
%                    % 4 ==>  bi-mineral rough fracture

if ~uniMineral && strcmp(geometry, 'fracture')
    LSCase.case = 4;
end
%%% ======= cylinder ===============
if LSCase.case == 1
    LSCase.xc = xc;
    LSCase.yc = yc;
    LSCase.diamcyl = diamcyl;
end
%%% ======= simple fracture ===============

if LSCase.case == 2
    LSCase.h = b;
end


%%

%%% ======= rough fracture ===============
if LSCase.case  == 3 || LSCase.case == 4
%%%=========boundary generating function for rough fractures ==============
    fy = @(x, a, lambda, x_0)a.*sin(2*pi/lambda * (x-x_0));
    DOMAIN.x_0 = x_0;

end
%%%========= uni mineral fracture case
if LSCase.case == 3
    refine = 10;
    LSCase.BoundaryCurve = uniMineralBoundaryCurve(DOMAIN, ...
        b, fy, lambda, a, x_0, refine);
    LSCase.BoundaryCurve.symmetry = symmetry;
end

%%%============= biMineral fracture case
if LSCase.case == 4
    refine_x = 10;
    refine_y = 10;

    % the start and end of top and bottom for mineral 1
    xt_s = [0 0.8 1.5 3];
    xt_e = [0.5 1.2 2 3.64];

    xb_s = [0 2 3.4];
    xb_e = [0.8 3 3.64]; 

    
%     xt_s = [0 0.8 1.5 2.5];
%     xt_e = [0.5 1.2 2 4];
% 
%     xb_s = [0 2 3.5];
%     xb_e = [0.8 3 4]; 

    [LSCase.BoundaryCurve1] = biMineralBoundaryCurve(DOMAIN, fy, b, ...
        lambda, a, x_0, refine_x, refine_y, xt_s, xt_e, xb_s, xb_e);

    % the start and end of top and bottom for mineral 2
    xt_s = [0.5 1.2 2];
    xt_e = [0.8 1.5 3];

    xb_s = [0.8 3];
    xb_e = [2 3.4];

%     xt_s = [0.5 1.2 2];
%     xt_e = [0.8 1.5 2.5];
% 
%     xb_s = [0.8 3];
%     xb_e = [2 3.5];

    [LSCase.BoundaryCurve2] = biMineralBoundaryCurve(DOMAIN, fy, b, ...
        lambda, a, x_0, refine_x, refine_y, xt_s, xt_e, xb_s, xb_e);
    
    
    % the whole solid matrix
    LSCase.BoundaryCurve_s = uniMineralBoundaryCurve(DOMAIN, ...
    b, fy, lambda, a, x_0, refine_x);
    LSCase.BoundaryCurve_s.symmetry = symmetry;


end

%% Initialize level set 
psi = LSInitialize(DOMAIN, LSCase);
%%

psi_f = psi;
if LSCase.case == 4
%     psi_f = min(psi.psi1, psi.psi2);
%     psi_f_sign = (psi.psi1>0).*(psi.psi2>0);
%     psi_f = psi_f.*psi_f_sign;
    psi_f = psi.psi_s;
end
    
psiU = interp2(DOMAIN.yp.*ones(imax+1,1), ones(jmax+1,1)'.*DOMAIN.xp', ...
    psi_f ,DOMAIN.yu.*ones(imax,1), ones(jmax+1,1)'.*DOMAIN.xu');

psiV = interp2(DOMAIN.yp.*ones(imax+1,1), ones(jmax+1,1)'.*DOMAIN.xp', ...
    psi_f ,DOMAIN.yv.*ones(imax+1,1), ones(jmax,1)'.*DOMAIN.xv');


%% initialize U, V, phi

U(:,2:end-1)=uinflow;
V(:,:) = 0;

U = double(psiU>0).*U;
V = double(psiV>0).*V;

U(1,:) = uinflow * (U(1,:)>0);
U_a(:)=U(1,:);
U_c(:)=0;
U_d(:)=0;




if LSCase.case ~= 4

    phi(:,:)=phi_init;
    phi = double(psi>0).*phi;
    phi_a(1,:) = ones(size(phi(1,:)))*phi_inlet.*(psi(1,:)>0);
    phi(1,:) = phi_a;
    phi_old=phi;
    LS.psi = psi;
else
    phi(:,:)=phi_init;
    phi = double(psi_f>0).*phi;
    phi_a(1,:) = ones(size(phi(1,:)))*phi_inlet.*(psi_f(1,:)>0);
    phi(1,:) = phi_a;
    phi_old=phi;
end

%% interface velocity parameters
methodInterfaceVelocity = 2;    % chooses the denominator of the 

if LSCase.case ~= 4
    molarVol = 36.9;
    interfaceVelocityCoeff = beta_phi * molarVol/D;
    %                   velocity, 1 for conditionally stable, 2 for stable
    Da = beta_phi;
else 
    molarVol.molarVol1 = 36.9;
    molarVol.molarVol2 = 36.9;

    interfaceVelocityCoeff.interfaceVelocityCoeff1 = ...
        beta_phi.beta_phi1 * molarVol.molarVol1/D;
     interfaceVelocityCoeff.interfaceVelocityCoeff2 = ...
        beta_phi.beta_phi2 * molarVol.molarVol2/D;   
    
    Da.Da1 = beta_phi.beta_phi1;
    Da.Da2 = beta_phi.beta_phi2;
end
%% LS solver variables

fac = 0.5;    % dta factor 
dtau = fac*min(min(DOMAIN.dxp));  % psedue time step

n_iter_ReLS = 4;

TimeSchemeLS = "RK3";
TimeSchemeRLS = "RK3";

% psi = sqrt( (DOMAIN.Xp-xc).^2 + (DOMAIN.Yp-yc).^2 )-diamcyl/2;
% LS = double(d>0)*d + double(d<0)*(d-r);

%% ========================== Time integration ============================

alpha_u=0.7;
alpha_v=alpha_u;
alpha_p=1;
alpha_q = 1;
VARIABLES = struct('D',D,'Re',Re,'density', density,'molarVol',molarVol,...
    'dimensional',dimensional,...
    'intVelCoeff',interfaceVelocityCoeff,'intVelMethod',methodInterfaceVelocity,...
    'alpha_u',alpha_u,'alpha_v',alpha_v,'alpha_p',alpha_p,'alpha_q',alpha_q,...
    'dt',dt,'dtVec',Dt,'Da',Da,'Pe',Pe,...
    'n_iter_ReLS',n_iter_ReLS,'dtau',dtau,'TimeSchemeLS',TimeSchemeLS,...
    'TimeSchemeRLS',TimeSchemeRLS, 'dissolution', dissolution);


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
    'nrgrainy',nrgrainy,'BQu',BQu,'BQv',BQv,'BQp',BQp, 'treshold', 100 *eps);

StateVar = struct('U',U,'V',V,'P',P,'phi_old',phi_old);
StateVar.P_cor_vec = zeros(1,(DOMAIN.imax-1)*(DOMAIN.jmax-1))';
if LSCase.case == 4
    LS.LS1 = struct('psi',psi.psi1);
    [LS.LS1] = LSnormals(LS.LS1,DOMAIN);
    
    LS.LS2 = struct('psi',psi.psi2);
    [LS.LS2] = LSnormals(LS.LS2,DOMAIN);    
    
    LS.LS_s = struct('psi',psi.psi_s);
    [LS.LS_s] = LSnormals(LS.LS_s,DOMAIN);    
    
else 
    [LS] = LSnormals(LS,DOMAIN);    

end

StateVar.U_old = StateVar.U;
StateVar.V_old = StateVar.V;
StateVar.P_old = StateVar.P;

% StateVar.U_star = StateVar.U;
% StateVar.V_star = StateVar.V;
% StateVar.U_star_old = StateVar.U_star;
% StateVar.V_star_old = StateVar.V_star;
        
%% ==============================

StateVar.phi = StateVar.phi_old;
VARIABLES.LSband = 10*min(min(DOMAIN.dxp));
VARIABLES.nLSupdate = 10;
VARIABLES.dtLS =  VARIABLES.dt;
VARIABLES.LSgamma = 6 * min(min(DOMAIN.dxp));
VARIABLES.LSbeta = 3 * min(min(DOMAIN.dxp));
VARIABLES.Pe_vel = VARIABLES.Pe;
LS.case = LSCase.case;
% end