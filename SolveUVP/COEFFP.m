function [CM,CM2,a_p,a_w,a_e,a_s,a_n, S_p] = COEFFP(d_u,d_v,BC,DOMAIN)
% dy==>dyv(j-1)
% dx==>dxu(i-1)
% flag=flag_c;
% landa_ext=landa_ext_c(:,:,:,:);
% X_e_ext=X_e_ext_c(:,:,:,:);
% X=XYp;
                      
% D=0.3;
% Desciption:
% 
% This function defines the coefficient matrix for the pressure solver.
% The pressure equation is a 2D poisson equation with:
% 1. dp/dx = 0 (hom. Neumann condition) on left boundary
% 2. p = 0 (Dirichlet condition) on the right boundary
% 3. dp/dy = 0 on bottom boundary
% 4. dp/dy = 0 on top boundary

%% Define parameters
imax = DOMAIN.imax;
jmax = DOMAIN.jmax;
dxu = DOMAIN.dxu;
dyv = DOMAIN.dyv;
BC_e_p = BC.BC_e_p;
BC_w_p = BC.BC_w_p;
BC_n_p = BC.BC_n_p;
BC_s_p = BC.BC_s_p;


L    = (imax-1)*(jmax-1); % Length of the matrix
jump = imax-1;
% I_g=0;
% I_ext=zeros(1,nrgrain);

%% Initialize Storage

%% sparce diagonal indexes


AE_I=1:L-1;
AE_J=2:L;

AW_I=2:L;
AW_J=1:L-1;

AN_I=1:L-jump;
AN_J=jump+1:L;

AS_I=jump+1:L;
AS_J=1:L-jump;

AP_I=1:L;
AP_J=1:L;

[a_p,a_w,a_e,a_s,a_n, S_p]=deal(zeros(imax+1,jmax+1));


AE = -d_u(3:imax-1,2:jmax).*dyv(3:imax-1,1:jmax-1);
      
AW = -d_u(2:imax-2,2:jmax).*dyv(2:imax-2,1:jmax-1);
      
AN = -d_v(2:imax,3:jmax-1).*dxu(1:imax-1,3:jmax-1);
      
AS = -d_v(2:imax,2:jmax-2).*dxu(1:imax-1,2:jmax-2);


%% West Boundary i=2
AW_w = zeros(1,size(AW,2));


if BC_w_p == 1
    AE_w = -d_u(2,2:jmax).*dyv(2,1:jmax-1);
    S_p(2,2:jmax) = d_u(1,2:jmax).*dyv(1,1:jmax-1);
elseif BC_w_p == 3 
    AE_w = -d_u(2,2:jmax).*dyv(2,1:jmax-1);
end

%% East Boundary i=imax
[AE_e]=deal(zeros(1,size(AE,2)));

if BC_e_p == 1
    AW_e = -d_u(imax-1,2:jmax).*dyv(imax-1,1:jmax-1);
    S_p(imax,2:jmax) =  d_u(imax,2:jmax).*dyv(imax,1:jmax-1);
elseif BC_e_p == 3 
    AW_e = -d_u(imax-1,2:jmax).*dyv(imax-1,1:jmax-1);
end
%% South Boundary j=2

if BC_s_p == 1
    AN_s = -d_v(2:imax,2).*dxu(1:imax-1,2);
elseif BC_s_p == 3 
    AN_s = -d_v(2:imax,2).*dxu(1:imax-1,2);
end

%% North Boundary j=jmax
if BC_n_p == 1
    AS_n = -d_v(2:imax,jmax-1).*dxu(1:imax-1,jmax-1);
elseif BC_n_p == 3 
    AS_n = -d_v(2:imax,jmax-1).*dxu(1:imax-1,jmax-1);
end


AW=[AW_w; AW; AW_e];
AE=[AE_w; AE; AE_e];
AS=[AS AS_n];
AN=[AN_s AN];


a_w(2:imax,2:jmax) = AW;
a_e(2:imax,2:jmax) = AE;
a_s(2:imax,3:jmax) = AS;
a_n(2:imax,2:jmax-1) = AN;

a_p(2:imax,2:jmax) = -(a_w(2:imax,2:jmax) + a_e(2:imax,2:jmax) +...
                        a_s(2:imax,2:jmax) + a_n(2:imax,2:jmax)) + S_p(2:imax,2:jmax);

%% Form Vectors
AW = reshape(AW,[1,numel(AW)]);
AW = AW(2:end);


AE = reshape(AE,[1,numel(AE)]);
AE = AE(1:end-1);

AS = reshape(AS,[1,numel(AS)]);

AN = reshape(AN,[1,numel(AS)]);

AP = a_p(2:imax,2:jmax);
AP = reshape(AP,[1,numel(AP)]);

AP ( AP==0 ) = 1;

CM = sparse(AE_I,AE_J,AE,L,L)+sparse(AW_I,AW_J,AW,L,L)+...
  sparse(AN_I,AN_J,AN,L,L)+sparse(AS_I,AS_J,AS,L,L)+...
  sparse(AP_I,AP_J,AP,L,L);

CM2=CM;
J = ceil((jmax-1)/2);
I = (J-2) * (DOMAIN.imax-1) + imax-1;
CM2(I,:) = [];
CM2(:, I) = [];
% CM2(end,:) = [];
% CM2(:, end) = [];

