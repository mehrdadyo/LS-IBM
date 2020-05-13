function [phi] = FORMPHI(DOMAIN,phi_vec,BC,phi_A,flag,phi_inside, LS)

% Description:
%
% This function solves the pressure Poisson equation, CM*PVECTOR = RHSP,
% for PVECTOR (correction pressure field in vector form). CM is the 
% coefficient matrix and RHSP is the RHS of the Poisson equation in 
% vector from.

imax = DOMAIN.imax;
jmax = DOMAIN.jmax;

%% Initialize Storage
phi       = zeros(imax+1,jmax+1);

%% Put solution for pressure into matrix form
ind = 0;

for j = 2:jmax
  for i = 2:imax
      
    ind = ind + 1;
    phi(i,j) = phi_vec(ind);
  end
end

%% Impose boundary conditions on correction phi
phi(2:imax,1)       = phi(2:imax,2);     % Bottom (dpcor/dy=0)
phi(2:imax,jmax+1)  = phi(2:imax,jmax);  % Top    (dpcor/dy=0)

if BC.BC_w_phi == 1
    phi(1,1:jmax+1)       = 2*phi_A(1,:)-phi(2,1:jmax+1);     % Left   (dpcor/dx=0)
elseif BC.BC_w_phi == 3
    phi(1,1:jmax+1)       =  phi(2,1:jmax+1);     % Left   (dpcor/dx=0)    
end

phi(imax+1,1:jmax+1)  = phi(imax,1:jmax+1); % Right  (pcor=0)
for j=1:jmax+1
    for i=1:imax+1
        if flag(i,j)== 2 
            phi(i,j)=phi_inside;
        end
    end
end
% % get indices of all cells which are in fluids and have negative phi
% negIndx = phi<0 & LS.psi>0;
% % get number of cells that are inside fluid and have positive phi
% insideVol = sum(sum((phi>0 & LS.psi>0)));
% % get the sum of all the negative phis
% totNegflux = sum(sum(phi(negIndx)));
% % add the negative mass to the fluid cells
% phi(phi>0 & LS.psi>0) = phi(phi>0 & LS.psi>0) + totNegflux/insideVol;
% % reset the negative cells to zero
% phi(phi>0 & LS.psi>0) = 0;

return
end
