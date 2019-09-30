function [phi_x] = biLinearInterpolation(x,y,phi,x_pr,y_pr)


I = find(x<x_pr);
J = find(y<y_pr);

I=I(end);
J=J(end);
                
I1 = I;
I2 = I1+1;
I3 = I1;
I4 = I1+1;
      
J1 = J;
J2 = J1;   
J3 = J1+1;
J4 = J1+1;
                
x1=x(I1)-x(I1);
y1=y(J1)-y(J1);
            
x2=x(I2)-x(I1);
y2=y(J2)-y(J1);

x3=x(I3)-x(I1);
y3=y(J3)-y(J1);
                   
x4=x(I4)-x(I1);
y4=y(J4)-y(J1);

%% the x_pr, y_pr location within the rectangular local coordinates 
x_dp = x_pr-x(I1);
y_dp = y_pr -y(J1);
            
A = [ 1 x1 y1 x1*y1; 1 x2 y2 x2*y2;...
    1 x3 y3 x3*y3; 1 x4 y4 x4*y4];
                
                %% phi_x_pr calculation
b = [phi(I1,J1) phi(I2,J2) phi(I3,J3) phi(I4,J4)]';
a = A\b;
phi_x = a(1) + a(2)*x_dp + a(3)*y_dp + a(4)*x_dp*y_dp ;