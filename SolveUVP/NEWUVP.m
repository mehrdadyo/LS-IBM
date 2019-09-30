% ====================================================================
function [U,V,P] = NEWUVP (USTAR,VSTAR,P_star,Soln,...
    ControlVar,DOMAIN,BC,IBM_coeffU,IBM_coeffV,IBM_coeffP,VARIABLES,PISO)
%% Description:
%
% This function performs the velocity corrector step using the fractional
% time-step velocity and the correction pressure. The corrector step is 
% given by:
%
%    u^(n+1) = u** - dt*grad(PCOR)
%
% where n+1 denotes the new time level.

%% Retrieve Variables
iter = ControlVar.ii;
imax = DOMAIN.imax;
jmax = DOMAIN.jmax;



d_u = Soln.d_u;
d_v = Soln.d_v;

a_u_w = Soln.a_u_w;
a_u_e = Soln.a_u_e;
a_u_s = Soln.a_u_s;
a_u_n = Soln.a_u_n;
a_u_p = Soln.a_u_p;

a_v_w = Soln.a_v_w;
a_v_e = Soln.a_v_e;
a_v_s = Soln.a_v_s;
a_v_n = Soln.a_v_n;
a_v_p = Soln.a_v_p;

U_a = BC.U_a;
U_b = BC.U_b;
U_c = BC.U_c;
U_d = BC.U_d;
V_a = BC.V_a;
V_b = BC.V_b;
V_c = BC.V_c;
V_d = BC.V_d;
BC_n_u = BC.BC_n_u;
BC_s_u = BC.BC_s_u;
BC_n_v = BC.BC_n_v;
BC_s_v = BC.BC_s_v;

BC_e_p = BC.BC_e_p;


alpha_p = VARIABLES.alpha_p;

%% Initialize Storage
U = zeros(imax,jmax+1);
V = zeros(imax+1,jmax);
P = zeros(imax+1,jmax+1);
%% Define parameters
% dtdxi = dt/dx;
% dtdyi = dt/dy;
% alpha_p=0.78;

%% Compute U velocity
if PISO == 1
    dU = Soln.dU_star;
    dV = Soln.dV_star;
end
                                                                                                                                                                                                                                
for j=2:jmax
    for i = 2:imax-1
        if PISO == 1
            
            U(i,j) = USTAR(i,j) + d_u(i,j)* (Soln.PCOR(i,j)- Soln.PCOR(i+1,j))-...
                ( a_u_w(i,j)*dU(i-1,j) + a_u_e(i,j)*dU(i+1,j) + ...
                a_u_s(i,j)*dU(i,j-1) + a_u_n(i,j)*dU(i,j+1) )...
                /a_u_p(i,j);

        else     
            U(i,j) = USTAR(i,j) + d_u(i,j)* (Soln.PCOR(i,j)-Soln.PCOR(i+1,j));
        end


    end
    
end

if BC_e_p == 1     
    if PISO == 1
            U(imax,2:jmax) = USTAR(imax,2:jmax) + ...
                d_u(imax,2:jmax).* (Soln.PCOR(imax,2:jmax)- ...
                Soln.PCOR(imax+1,2:jmax))-...
                ( a_u_w(imax,2:jmax).*dU(imax-1,2:jmax) + ...
                a_u_s(imax,2:jmax).*dU(imax,1:jmax-1) + ...
                a_u_n(imax,2:jmax).*dU(imax,3:jmax+1) )...
                ./a_u_p(imax,2:jmax);

        else     
            U(imax,2:jmax) = USTAR(imax,2:jmax) + ...
                d_u(imax,2:jmax).* ...
                (Soln.PCOR(imax,2:jmax)-Soln.PCOR(imax+1,2:jmax));
    end
end
    
    
    



% for j=2:jmax
%   for i = 2:imax-1
%             U(i,j) = USTAR(i,j) - (PCOR(i+1,j)-PCOR(i,j))*dtdxi;
%   end
% end
%% Apply BC's to U
U(1,:)     = U_a(1,:).* double(IBM_coeffU.flag_u(1,:)~=2);               % Left  (inlet, u=U_inlet(1,j))
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
if BC_n_u==1
    U(2:imax-1,jmax+1)  = 2*U_d(1,2:imax-1)'-U(2:imax-1,jmax);% Top    (wall, u=0)
elseif BC_n_u==3
    U(2:imax-1,jmax+1)  = U(2:imax-1,jmax);       % Top    (wall, u=0)
end

%% South BC for U
if BC_s_u==1
    U(2:imax-1,1)       = 2*U_c(1,2:imax-1)'-U(2:imax-1,2); % Bottom (wall, u=0)
elseif BC_s_u==3
     U(1:imax,1)       = U(1:imax,2);          % Bottom (wall, u=0)
end
% %% Extra B.C.'s needed for calculation of vorticity in corners
U(1,1)         = U_a(1,1);                  % Left- Bottom
U(1,jmax)    = U(1,end);                    % Left- Top
U(imax,1)      = U(imax,2);                 % Right- Bottom
U(imax,jmax+1) = U(imax,jmax);              % Right- Top

%% Compute V velocity
for j = 2:jmax-1
  for i = 2:imax
      if PISO == 1
                V(i,j) = VSTAR(i,j) + d_v(i,j)* (Soln.PCOR(i,j)-Soln.PCOR(i,j+1))-...
                    ( a_v_w(i,j)*dV(i-1,j) + a_v_e(i,j)*dV(i+1,j) + ...
                    a_v_s(i,j)*dV(i,j-1) + a_v_n(i,j)*dV(i,j+1) )...
                    /a_v_p(i,j);
                
       else     
                V(i,j) = VSTAR(i,j) + d_v(i,j)* (Soln.PCOR(i,j)-Soln.PCOR(i,j+1));
       end
      
  end
end

% for j = 2:jmax-1
%   for i = 2:imax
%         V(i,j) = VSTAR(i,j) - (PCOR(i,j+1)-PCOR(i,j))*dtdyi;
%       
%   end
% end

%% Apply BC's to V
%% North BC for V
if BC_n_v==1
    V(2:imax,jmax)  = V_d(1,2:imax)';           % Top    (wall, v=0)
elseif BC_n_v==3
    V(2:imax,jmax)  = V(2:imax,jmax-1);            % Top    (wall, v=0)
end

%% South BC for V
if BC_s_v==1
    V(2:imax,1)     =  V_c(1,2:imax)';             % Bottom (wall, v=0)
elseif BC_s_v==3
    V(2:imax,1)     = V(2:imax,2);                 % Bottom (wall, v=0)
end

% V(1,2:jmax-1)     = 2*V_a(1,2:jmax-1)-V(2,2:jmax-1); 
V(1,2:jmax-1)      = V(2,2:jmax-1);                % Left  (inlet, v=0)
V(imax+1,2:jmax-1)= V(imax,2:jmax-1);          % Right  (outlet, dv/dx=0)
%% Extra B.C.'s for calculation of vorticity in corner points
V(1,1)         = 0;%-V(2,1);                    % Bottom- Left
V(1,jmax)      = 0;%-V(1,jmax-1);               % Top- Left
V(imax+1,1)    = 0;%V(imax,1);                  % Bottom- Right
V(imax+1,jmax) = 0;%V(imax+1,jmax-1);           % Top- Right

%% Compute new pressure
for j = 2:jmax
  for i = 2:imax
    P(i,j) = P_star(i,j) +alpha_p* Soln.PCOR(i,j);
  end
end

% % %% No BC's needed for pressure, but computed here for completeness
P(2:imax,1)       = P(2:imax,2);     % Bottom (dp/dy=0)
P(2:imax,jmax+1)  = P(2:imax,jmax);  % Top    (dp/dy=0)
P(1,2:jmax)       = P(2,2:jmax);     % Left   (dp/dx=0)
P(imax+1,2:jmax)  = -P(imax,2:jmax); % Right  (p=0)
return
end
% ==============================================================
