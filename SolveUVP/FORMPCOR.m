function [PCOR] = FORMPCOR(PCORVEC,DOMAIN)

% Description:
%
% This function solves the pressure Poisson equation, CM*PVECTOR = RHSP,
% for PVECTOR (correction pressure field in vector form). CM is the 
% coefficient matrix and RHSP is the RHS of the Poisson equation in 
% vector from.


imax = DOMAIN.imax;
jmax = DOMAIN.jmax;
%% Initialize Storage
PCOR       = zeros(imax+1,jmax+1);
% PCORVEC1 = [PCORVEC ; 0];

PCOR       = zeros(imax+1,jmax+1);
J = ceil((jmax-1)/2);
I = (J-2) * (DOMAIN.imax-1) + imax-1;
PCORVEC1 = [PCORVEC(1:I-1,:) ; 0; PCORVEC(I:end,:)];

%% Put solution for pressure into matrix form
PCOR(2:imax,2:jmax) = reshape(PCORVEC1, [imax-1, jmax-1]);

% for j = 2:jmax
%   for i = 2:imax
%     ind = ind + 1;
%     if (i==imax) && (j==jmax)
%         PCOR(i,j) = 0;
%     else 
%         PCOR(i,j) =PCORVEC(ind);
%     end
%   end
% end
% PCOR(imax,jmax)=0;
% %% Impose boundary conditions on correction pressure

PCOR(2:imax,1)       = PCOR(2:imax,2);     % Bottom (dpcor/dy=0)
PCOR(2:imax,jmax+1)  = PCOR(2:imax,jmax);  % Top    (dpcor/dy=0)
PCOR(1,2:jmax)       = PCOR(2,2:jmax);     % Left   (dpcor/dx=0)
PCOR(imax+1,2:jmax)  = -PCOR(imax,2:jmax); % Right  (pcor=0)

return
end
