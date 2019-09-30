function [RHSP,RHSP2] = RHSP(u_star,v_star,DOMAIN,...
    IBM_coeffU,IBM_coeffV)
% Desciption:
% FV
%% Define parameters
%% Define parameters
imax = DOMAIN.imax;
jmax = DOMAIN.jmax;
dxu = DOMAIN.dxu;
dyv = DOMAIN.dyv;

flag_u = IBM_coeffU.flag_u;
flag_v = IBM_coeffV.flag_v;
% dx==>dxu(i-1)
% dy==>dyv(j-1)
% I_g=0;
% I_ext=zeros(1,nrgrain);

%% Vectors
%AP
% RHSSCALAR = zeros(L,1);
phi_inside=0;

for i=1:imax
    for j=1:jmax+1
        if flag_u(i,j)==1
            u_star(i,j)=phi_inside;
        end
    end
end


for i=1:imax+1
    for j=1:jmax
        if flag_v(i,j)==1
            v_star(i,j)=0;
        end
    end
end
%% Fill in sparse matrix


S0= u_star(1:imax-1,2:jmax).*dyv(2:imax,1:jmax-1) - ...
    u_star(2:imax,2:jmax).*dyv(2:imax,1:jmax-1) + ...
    v_star(2:imax,1:jmax-1).*dxu(1:imax-1,2:jmax) - ...
    v_star(2:imax,2:jmax).*dxu(1:imax-1,2:jmax);

RHSP=reshape(S0,[numel(S0),1]);

% S0(dirich)=0;
RHSP2=RHSP;
% J = ceil((jmax-1)/2);
% I = (J-2) * (DOMAIN.imax-1) + imax-1;
% RHSP2(I,:) = []; 
RHSP2(end,:) = []; 

