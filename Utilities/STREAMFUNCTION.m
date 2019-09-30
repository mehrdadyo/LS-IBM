function [PSI] = STREAMFUNCTION(dx,dy,U,V)

% Description:
%
% This function solves for the stream function using the V velocity. The
% streamfunction is defined as follows:
%
%    u =  d(psi)/dy
%    v = -d(psi)/dx
%
% psi is defined here at the top-right corner of a grid cell (whereas
% pressure is defined at cell center). psi(1,1)=0 lies exactly on the
% lower-left corner; psi(imax,jmax) lies exactly on the top-right corner.

%% Initialize psi
[Nx,Ny]=size(U);
PSI = zeros(Nx+1,Ny+1);

%% Compute stream function
for j=2:Ny
  PSI(1,j) = PSI(1,j-1) + dy(1,j+1)*U(1,j); % PSI(1,1)=0
end
for i=2:Nx
  PSI(i,2:end) = PSI(i-1,2:end) - dx(i+1,2:end-1).*V(i,:);
end

return
end
