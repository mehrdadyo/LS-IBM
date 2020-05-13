function [U,V] = FORMUV(U_vec,V_vec,BC,iter,DOMAIN,IBM_coeffU,IBM_coeffV)
% Description:
%
% This function solves the pressure Poisson equation, CM*PVECTOR = RHSP,
% for PVECTOR (correction pressure field in vector form). CM is the 
% coefficient matrix and RHSP is the RHS of the Poisson equation in 
% vector from.

imax = DOMAIN.imax;
jmax = DOMAIN.jmax;

U_a = BC.U_a;
U_b = BC.U_b;
U_c = BC.U_c;
U_d = BC.U_d;


V_a = BC.V_a;
V_b = BC.V_b;
V_c = BC.V_c;
V_d = BC.V_d;

BC_e_p = BC.BC_e_p;

%% Initialize Storage
U       = zeros(imax,jmax+1);
V       = zeros(imax+1,jmax);
%% Put solution for U into matrix form
if BC_e_p ~= 1
    
    U(2:imax-1,2:jmax) = reshape(U_vec, [imax-2, jmax-1]);

elseif BC_e_p == 1
    
    U(2:imax,2:jmax) = reshape(U_vec, [imax-1, jmax-1]);

end

% ind = 0;
% for j = 2:jmax
%   for i = 2:imax-1
%     ind = ind + 1;
%     U(i,j) = U_vec(ind);
%   end
% end

%% Impose boundary conditions on U and V
% % U 
U(1,2:jmax)     = U_a(1,2:jmax).* double(IBM_coeffU.flag_u(1,2:jmax)~=2);                   % Left  (inlet, u=U_inlet(1,j))

if BC_e_p ~= 1

    U(imax,2:jmax-1)= U(imax-1,2:jmax-1);              % Right  (outlet, du/dx=0)
    
    if iter>1
        M_in=sum(U(1,2:end-1));
        M_out=sum(U(end,2:end-1));
        U(imax,2:jmax-1)=M_in/M_out*U(imax-1,2:jmax-1); 
    end 

elseif BC_e_p == 1
    
    if iter>1
        M_in=sum(U(1,2:end-1));
        M_out=sum(U(end,2:end-1));
        U(imax,2:jmax-1)=M_in/M_out*U(imax,2:jmax-1); 
    end 
    
end

%% North BC for U
BC_n=BC.BC_n_u;
if BC_n==1 
    U(1:imax,jmax+1)  = 2*U_d -U(1:imax,jmax);       % Top    (wall, u=0)
elseif BC_n==3
    U(1:imax,jmax+1)  = U(1:imax,jmax);       % Top    (wall, u=0)
end

%% South BC for U
BC_s=BC.BC_s_u;
if BC_s==1
    U(1:imax,1)       = 2*U_c-U(1:imax,2);          % Bottom (wall, u=0)
elseif BC_s==3
    U(1:imax,1)       = U(1:imax,2);          % Bottom (wall, u=0)
end
%% Extra B.C.'s needed for calculation of vorticity in corners
% U(1,1)         = U_a(1,1).* double(IBM_coeffU.flag_u(1,1)~=2);                  % Left- Bottom
% U(1,jmax+1)    = U_a(1,end).* double(IBM_coeffU.flag_u(1,end)~=2);                    % Left- Top
U(imax,1)      = U(imax,2);                 % Right- Bottom
U(imax,jmax+1) = U(imax,jmax);              % Right- Top



%% Put solution for V into matrix form
ind = 0;
for j = 2:jmax-1
  for i = 2:imax
    ind = ind + 1;
    V(i,j) = V_vec(ind);
  end
end
%% North BC for V
BC_n=BC.BC_n_v;
if BC_n==1
    V(2:imax,jmax)  = V_d(2:imax,1);           % Top    (wall, v=0)
elseif BC_n==3
    V(2:imax,jmax)  = V(2:imax,jmax-1);            % Top    (wall, v=0)
end
%% South BC for V
BC_s=BC.BC_s_v;
if BC_s==1
    V(2:imax,1)     =  V_c(2:imax,1);             % Bottom (wall, v=0)
elseif BC_s==3
    V(2:imax,1)     = V(2:imax,2);                 % Bottom (wall, v=0)
end


V(1,2:jmax-1)      = -V(2,2:jmax-1);                % Left  (inlet, v=0)
V(imax+1,2:jmax-1)= V(imax,2:jmax-1);          % Right  (outlet, dv/dx=0)

%% Extra B.C.'s for calculation of vorticity in corner points
% V(1,1)         = 0;%-V(2,1);                    % Bottom- Left
% V(1,jmax)      = 0;%-V(1,jmax-1);               % Top- Left
% V(imax+1,1)    = 0;%V(imax,1);                  % Bottom- Right
% V(imax+1,jmax) = 0;%V(imax+1,jmax-1);           % Top- Right

return
end
