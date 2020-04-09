function [RHSP,RHSP2] = RHSP_PISO(Soln,u_star,v_star,DOMAIN,IBM_coeffU,...
    IBM_coeffV)
% Desciption:
% Soln.RHSP,Soln.RHSP2] = RHSP_PISO(Soln,StateVar.U_star,...
%             StateVar.V_star,DOMAIN,IBM_coeffU,IBM_coeffV);
% % FV
%% Define parameters
%% Define parameters
imax = DOMAIN.imax;
jmax = DOMAIN.jmax;
dxu = DOMAIN.dxu;
dyv = DOMAIN.dyv;

[Su1, Su2, Sv1, Sv2] = deal(zeros(imax-1,jmax-1));

dU = Soln.dU_star;
dV = Soln.dV_star;

a_u_w = Soln.a_u_w;
a_u_e = Soln.a_u_e;
a_u_s = Soln.a_u_s;
a_u_n = Soln.a_u_n;

a_v_w = Soln.a_v_w;
a_v_e = Soln.a_v_e;
a_v_s = Soln.a_v_s;
a_v_n = Soln.a_v_n;

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
% Su1(1,:) = u_star(1,2:jmax).*dyv(1,1:jmax-1);
Su1(2:imax-1,:) = -Soln.d_u_sec(2:imax-1,2:jmax).* ...
    ( a_u_w(2:imax-1,2:jmax).*dU(1:imax-2,2:jmax) + ... 
    a_u_e(2:imax-1,2:jmax).*dU(3:imax,2:jmax) + ...
    a_u_s(2:imax-1,2:jmax).*dU(2:imax-1,1:jmax-1) + ...
    a_u_n(2:imax-1,2:jmax).*dU(2:imax-1,3:jmax+1) );

% Su2(imax-1,:) = u_star(imax,2:jmax).*dyv(imax,1:jmax-1);
Su2(1:imax-2,:) = -Soln.d_u_sec(2:imax-1,2:jmax).* ...
    ( a_u_w(2:imax-1,2:jmax).*dU(1:imax-2,2:jmax) + ... 
    a_u_e(2:imax-1,2:jmax).*dU(3:imax,2:jmax) + ...
    a_u_s(2:imax-1,2:jmax).*dU(2:imax-1,1:jmax-1) + ...
    a_u_n(2:imax-1,2:jmax).*dU(2:imax-1,3:jmax+1) );

% Sv1(:,1) = v_star(2:imax,1).*dxu(1:imax-1,1);
Sv1(:,2:jmax-1) = -Soln.d_v_sec(2:imax,2:jmax-1).* ...
    ( a_v_w(2:imax,2:jmax-1).*dV(2:imax,1:jmax-2) + ... 
    a_v_e(2:imax,2:jmax-1).*dV(2:imax,3:jmax) + ...
    a_v_s(2:imax,2:jmax-1).*dV(1:imax-1,2:jmax-1) + ...
    a_v_n(2:imax,2:jmax-1).*dV(3:imax+1,2:jmax-1) );

% Sv2(:,jmax-1) = v_star(2:imax,jmax).*dxu(1:imax-1,jmax);
Sv2(:,1:jmax-2) = -Soln.d_v_sec(2:imax,2:jmax-1).* ...
    ( a_v_w(2:imax,2:jmax-1).*dV(2:imax,1:jmax-2) + ... 
    a_v_e(2:imax,2:jmax-1).*dV(2:imax,3:jmax) + ...
    a_v_s(2:imax,2:jmax-1).*dV(1:imax-1,2:jmax-1) + ...
    a_v_n(2:imax,2:jmax-1).*dV(3:imax+1,2:jmax-1) );

S0 = Su1 - Su2 + Sv1 - Sv2;



RHSP=reshape(S0,[numel(S0),1]);

% S0(dirich)=0;
% RHSP2=RHSP(1:end-1,1);
RHSP2=RHSP;
J = ceil((jmax-1)/2);
I = (J-2) * (DOMAIN.imax-1) + imax-1;
RHSP2(I,:) = []; 
