function [CM,AP,AW,AE,AS,AN] = COEFFPHIDRE(DOMAIN,Flux,IBM_coeffP)
            
% D=0.3;
% Desciption:
% 
% This function defines the coefficient matrix for the pressure solver.
% The pressure equation is a 2D poisson equation with:
% 1. dp/dx = 0 (hom. Neumann condition) on left boundary
% 2. p = 0 (Dirichlet condition) on the right boundary
% 3. dp/dy = 0 on bottom boundary
% 4. dp/dy = 0 on top boundary

%  BC=1 ==>     Dirichlet
%  BC=2 ==>     Constant total flux
%  BC=3 ==>     % zero diffusive

imax = DOMAIN.imax;
jmax = DOMAIN.jmax;
%% Define parameters
De = Flux.Diffu_P.De_p;
Dw = Flux.Diffu_P.Dw_p;
Dn = Flux.Diffu_P.Dn_p;
Ds = Flux.Diffu_P.Ds_p;

D1_a = Flux.Diffu_P.D1_a;
D2_a = Flux.Diffu_P.D2_a;

A0_p = Flux.Diffu_P.A0_p_p;

landa_g_1 = IBM_coeffP.landa_g_1_p;
landa_g_2 = IBM_coeffP.landa_g_2_p;
landa_g_3 = IBM_coeffP.landa_g_3_p;
landa_g_4 = IBM_coeffP.landa_g_4_p;

I_solid = IBM_coeffP.I_solid_p;
J_solid = IBM_coeffP.J_solid_p;
I_g = IBM_coeffP.I_g_p;
J_g = IBM_coeffP.J_g_p;
I_e = IBM_coeffP.I_e_p;
J_e = IBM_coeffP.J_e_p;

numg = IBM_coeffP.numg_p;
%% Define parameters

L    = (imax-1)*(jmax-1); % Length of the matrix
jump = imax-1;

%% Initialize Storage



%% Vectors
%AP
AP_I=zeros(1,L);
AP_J=zeros(1,L);
AP=zeros(1,L);
Su=zeros(1,L);

%AW
AW_I=zeros(1,L-1);
AW_J=zeros(1,L-1);
AW=zeros(1,L-1);

%AE
AE_I=zeros(1,L-1);
AE_J=zeros(1,L-1);
AE=zeros(1,L-1);

%AS
AS_I=zeros(1,L-jump);
AS_J=zeros(1,L-jump);
AS=zeros(1,L-jump);

%AN
AN_I=zeros(1,L-jump);
AN_J=zeros(1,L-jump);
AN=zeros(1,L-jump);

%A_gs
A1_I_g=zeros(1,numg);
A2_I_g=zeros(1,numg);
A3_I_g=zeros(1,numg);
A4_I_g=zeros(1,numg);

A1_J_g=zeros(1,numg);
A2_J_g=zeros(1,numg);
A3_J_g=zeros(1,numg);
A4_J_g=zeros(1,numg);

A1_g=zeros(1,numg);
A2_g=zeros(1,numg);
A3_g=zeros(1,numg);
A4_g=zeros(1,numg);



%% Fill in sparse matrix
ind = 0;
for j = 2:jmax      % Loop over all interior pressure points
  for i = 2:imax
    ind = ind+1; 

    if(i == 2 && j == 2)                    % Bottom-left corner
        AE_I(ind)=ind;
        AE_J(ind)=ind+1;
        AE(ind)=-(De(i,j)-D1_a(j));
             
        AN_I(ind)=ind;
        AN_J(ind)=ind+jump;
        AN(ind)=-(Dn(i,j));
      
        AP_I(ind)=ind;
        AP_J(ind)=ind;
        S_p=D2_a(j);
        AP(ind)=-(AE(ind)+AN(ind))+A0_p(i,j)+S_p;
          
    elseif(i == 2 && j == jmax)             % Top-left corner
      
        AE_I(ind)=ind;
        AE_J(ind)=ind+1;
        AE(ind)=-(De(i,j)-D1_a(j));
      
      
        AW_I(ind-1)=ind;
        AW_J(ind-1)=ind-1;
        AW(ind-1)=0;
      
        AS_I(ind-jump)=ind;
        AS_J(ind-jump)=ind-jump;
        AS(ind-jump)=-(Ds(i,j));
      
        AP_I(ind)=ind;
        AP_J(ind)=ind;
        S_p=D2_a(j);
        AP(ind)=-(AE(ind)+AS(ind-jump))+A0_p(i,j)+S_p;
      

    elseif(i == imax && j == 2)             % Bottom-right corner
        AE_I(ind)=ind;
        AE_J(ind)=ind+1;
        AE(ind)=0;
      
        AW_I(ind-1)=ind;
        AW_J(ind-1)=ind-1;
        AW(ind-1)=-(Dw(i,j));
      
        AN_I(ind)=ind;
        AN_J(ind)=ind+jump;
        AN(ind)=-(Dn(i,j));
        
        S_p=D2_a(j);
        AP_I(ind)=ind;
        AP_J(ind)=ind;
        AP(ind)=-(AW(ind-1)+AN(ind))+A0_p(i,j)+S_p;
        
    elseif(i == imax && j == jmax)          % Top-right corner
        AW_I(ind-1)=ind;
        AW_J(ind-1)=ind-1;
        AW(ind-1)=-(Dw(i,j));
      
        AS_I(ind-jump)=ind;
        AS_J(ind-jump)=ind-jump;
        AS(ind-jump)=-(Ds(i,j));
      
        S_p=0;
        AP_I(ind)=ind;
        AP_J(ind)=ind;
        AP(ind)=-(AW(ind-1)+AS(ind-jump))+A0_p(i,j)+S_p;
        
          
    elseif(i == 2)      % Left edge
        AE_I(ind)=ind;
        AE_J(ind)=ind+1;
        AE(ind)=-(De(i,j)-D1_a(j));
      
        AW_I(ind-1)=ind;
        AW_J(ind-1)=ind-1;
        AW(ind-1)=0;
      
        AN_I(ind)=ind;
        AN_J(ind)=ind+jump;
        AN(ind)=-(Dn(i,j));
      
        AS_I(ind-jump)=ind;
        AS_J(ind-jump)=ind-jump;
        AS(ind-jump)=-(Ds(i,j));
      
        S_p=D2_a(j);
        AP_I(ind)=ind;
        AP_J(ind)=ind;
        AP(ind)=-(AE(ind)+AW(ind-1)+AN(ind)+AS(ind-jump))+A0_p(i,j)+S_p;
      
    elseif(i == imax)                       % Right edge
        AE_I(ind)=ind;
        AE_J(ind)=ind+1;
        AE(ind)=0;
      
        AW_I(ind-1)=ind;
        AW_J(ind-1)=ind-1;
        AW(ind-1)=-(Dw(i,j));
      
      
        AN_I(ind)=ind;
        AN_J(ind)=ind+jump;
        AN(ind)=-(Dn(i,j));
      
        AS_I(ind-jump)=ind;
        AS_J(ind-jump)=ind-jump;
        AS(ind-jump)=-(Ds(i,j));
      
        S_p=0;
        AP_I(ind)=ind;
        AP_J(ind)=ind;
        AP(ind)=-(AE(ind)+AW(ind-1)+AN(ind)+AS(ind-jump))+A0_p(i,j)+S_p;
      
    elseif(j == 2)                          % Bottom edge
            
      AE_I(ind)=ind;
      AE_J(ind)=ind+1;
      AE(ind)=-(De(i,j));
      
      AW_I(ind-1)=ind;
      AW_J(ind-1)=ind-1;
      AW(ind-1)=-(Dw(i,j));
      
      AN_I(ind)=ind;
      AN_J(ind)=ind+jump;
      AN(ind)=-(Dn(i,j));
      
      S_p=0;
      AP_I(ind)=ind;
      AP_J(ind)=ind;
      AP(ind)=-(AE(ind)+AW(ind-1)+AN(ind))+A0_p(i,j)+S_p;
    elseif(j == jmax)                       % Top edge
            
      AE_I(ind)=ind;
      AE_J(ind)=ind+1;
      AE(ind)=-(De(i,j));
      
      AW_I(ind-1)=ind;
      AW_J(ind-1)=ind-1;
      AW(ind-1)=-(Dw(i,j));
      
      AS_I(ind-jump)=ind;
      AS_J(ind-jump)=ind-jump;
      AS(ind-jump)=-(Ds(i,j));
      
      S_p=0;
      AP_I(ind)=ind;
      AP_J(ind)=ind;
      AP(ind)=-(AE(ind)+AW(ind-1)+AS(ind-jump))+A0_p(i,j)+S_p;
    else                                    % All others
        AE_I(ind)=ind;
        AE_J(ind)=ind+1;
        AE(ind)=-(De(i,j));
      
        AW_I(ind-1)=ind;
        AW_J(ind-1)=ind-1;
        AW(ind-1)=-(Dw(i,j));
      
        AN_I(ind)=ind;
        AN_J(ind)=ind+jump;
        AN(ind)=-(Dn(i,j));
      
        AS_I(ind-jump)=ind;
        AS_J(ind-jump)=ind-jump;
        AS(ind-jump)=-(Ds(i,j));
      
        S_p=0;
        AP_I(ind)=ind;
        AP_J(ind)=ind;
        AP(ind)=-(AE(ind)+AW(ind-1)+AN(ind)+AS(ind-jump))+A0_p(i,j)+S_p;
        
            
    end
  end
end

%% Body Cells

%****************************
for i=1:size(I_solid,2)
%     change values for An As Ae Aw Ap and RHS.
    ind=(J_solid(i)-2)*(imax-1)+I_solid(i)-1;
    
    AE_I(ind)=ind;
    AE_J(ind)=ind+1;
    AE(ind)=0;
      
    AW_I(ind-1)=ind;
    AW_J(ind-1)=ind-1;
    AW(ind-1)=0;
      
    AN_I(ind)=ind;
    AN_J(ind)=ind+jump;
    AN(ind)=0;
      
    AS_I(ind-jump)=ind;
    AS_J(ind-jump)=ind-jump;
    AS(ind-jump)=0;
      
    AP_I(ind)=ind;
    AP_J(ind)=ind;
    AP(ind)=1;
   
end

%****************************
for i=1:size(I_g,2)
    
    ind=(J_g(i)-2)*(imax-1)+I_g(i)-1;
    % FIND Ghost VIRTUAL POINT LOCATION
    ii_e=I_e(i);
    jj_e=J_e(i);
        
    i_e=ii_e;
    j_e=jj_e;
        
    A1_I_g(1,i)=ind;
    A1_J_g(1,i)=(j_e-2)*(imax-1)+i_e-1;
    A1_g(1,i)=-landa_g_1(1,i);
        
    i_e=ii_e+1;
    j_e=jj_e;
    A2_I_g(1,i)=ind;
    A2_J_g(1,i)=(j_e-2)*(imax-1)+i_e-1;
    A2_g(1,i)=-landa_g_2(1,i);
        

    i_e=ii_e;
    j_e=jj_e+1;
    A3_I_g(1,i)=ind;
    A3_J_g(1,i)=(j_e-2)*(imax-1)+i_e-1;
    A3_g(1,i)=-landa_g_3(1,i);
        
    i_e=ii_e+1;
    j_e=jj_e+1;
    A4_I_g(1,i)=ind;
    A4_J_g(1,i)=(j_e-2)*(imax-1)+i_e-1;
    A4_g(1,i)=-landa_g_4(1,i);
        


        
        
    AE_I(ind)=ind;
    AE_J(ind)=ind+1;
    AE(ind)=0;
      
    AW_I(ind-1)=ind;
    AW_J(ind-1)=ind-1;
    AW(ind-1)=0;
      
    AN_I(ind)=ind;
    AN_J(ind)=ind+jump;
    AN(ind)=0;
      
    AS_I(ind-jump)=ind;
    AS_J(ind-jump)=ind-jump;
    AS(ind-jump)=0;
      
    AP_I(ind)=ind;
    AP_J(ind)=ind;
    AP(ind)=1;
                
    
end

A_g_sparse= sparse(A1_I_g,A1_J_g,A1_g,L,L)+...
            sparse(A2_I_g,A2_J_g,A2_g,L,L)+...
            sparse(A3_I_g,A3_J_g,A3_g,L,L)+...
            sparse(A4_I_g,A4_J_g,A4_g,L,L);

CM=sparse(AE_I,AE_J,AE,L,L)+sparse(AW_I,AW_J,AW,L,L)+...
  sparse(AN_I,AN_J,AN,L,L)+sparse(AS_I,AS_J,AS,L,L)+...
  sparse(AP_I,AP_J,AP,L,L)+A_g_sparse;
