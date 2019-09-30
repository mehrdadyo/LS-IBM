function [VOR] = VORTICITY(dx,dy,U,V)

% Description:
%
% This function computes the vorticity on the top-right corner of the
% grid cells.

%% Define parameters
dxi = 1./dx;
dyi = 1./dy;

%% Initialize Storage
[Nx,Ny]=size(U);
VOR = zeros(Nx+1,Ny+1);
%% Compute Vorticity
for j=1:Ny-1
  for i=1:Nx-1
    VOR(i,j) = (V(i+1,j)-V(i,j))*dxi(i+1,j) - (U(i,j+1)-U(i,j))*dyi(i,j+1);
  end
end

return
end
