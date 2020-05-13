function [CM,ap_p,aw_p,ae_p,as_p,an_p] = ... 
    COEFFPHIADRE(DOMAIN,Flux,IBM_coeffP,IBM,VARIABLES)
            
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
alpha_q = VARIABLES.alpha_q;
imax = DOMAIN.imax;
jmax = DOMAIN.jmax;
%% Define parameters
De = Flux.Diffu_P.De_p;
Dw = Flux.Diffu_P.Dw_p;
Dn = Flux.Diffu_P.Dn_p;
Ds = Flux.Diffu_P.Ds_p;

D1_a = Flux.Diffu_P.D1_a;
D2_a = Flux.Diffu_P.D2_a;

Fe = Flux.ConvF_phi.Fe_p;
Fw = Flux.ConvF_phi.Fw_p;
Fn = Flux.ConvF_phi.Fn_p;
Fs = Flux.ConvF_phi.Fs_p;
A0_p = Flux.ConvF_phi.A0_p_p;
dF = Flux.ConvF_phi.dF_p;

alphae = Flux.ConvF_phi.alphae_p;
alphaw = Flux.ConvF_phi.alphaw_p;
alphan = Flux.ConvF_phi.alphan_p;
alphas = Flux.ConvF_phi.alphas_p;

landa_g_1 = IBM_coeffP.landa_g_1_p;
landa_g_2 = IBM_coeffP.landa_g_2_p;
landa_g_3 = IBM_coeffP.landa_g_3_p;
landa_g_4 = IBM_coeffP.landa_g_4_p;

I_solid = IBM_coeffP.I_solid_p;
J_solid = IBM_coeffP.J_solid_p;
I_g = IBM_coeffP.I_g_p;
J_g = IBM_coeffP.J_g_p;


I1=IBM_coeffP.I1_p;
I2=IBM_coeffP.I2_p;
I3=IBM_coeffP.I3_p;
I4=IBM_coeffP.I4_p;

J1=IBM_coeffP.J1_p;
J2=IBM_coeffP.J2_p;
J3=IBM_coeffP.J3_p;
J4=IBM_coeffP.J4_p;


numg = IBM_coeffP.numg_p;
%% Define parameters

L    = (imax-1)*(jmax-1); % Length of the matrix
jump = imax-1;

%% Initialize Storage

%% Vectors

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


if IBM.BQp == 1
    
    landa_g_5=IBM_coeffU.landa_g_5_u;
    landa_g_6=IBM_coeffU.landa_g_6_u;
    
    I5=IBM_coeffU.I5_u;
    I6=IBM_coeffU.I6_u;
    
    J5=IBM_coeffU.J5_u;
    J6=IBM_coeffU.J6_u;
    
    A5_I_g=zeros(1,numg);
    A6_I_g=zeros(1,numg);
    
    A5_J_g=zeros(1,numg);
    A6_J_g=zeros(1,numg);
    
    A5_g=zeros(1,numg);
    A6_g=zeros(1,numg);
end

%% Matrices
aw_p=zeros(imax+1,jmax+1);
ae_p=zeros(imax+1,jmax+1);
as_p=zeros(imax+1,jmax+1);
an_p=zeros(imax+1,jmax+1);
ap_p=zeros(imax+1,jmax+1);
S_p=zeros(imax+1,jmax+1);



S_p_w=zeros(1,jmax-1);
S_p_e=zeros(1,jmax-1);

S_p_s=zeros(imax-1,1);
S_p_n=zeros(imax-1,1);


%% sparse diagonal indexes

AE_I=1:L-1;
AE_J=2:L;

AW_I=2:L;
AW_J=1:L-1;

AN_I=1:L-jump;
AN_J=jump+1:L;

AS_I=jump+1:L;
AS_J=1:L-jump;

AP_I=1:L;
AP_J=1:L;


%% Inner Points: 

% We have (Nx+1,Ny+1) Scalar-cells. All other points on the margin are ghost
% points or we have their values from BC. Therefore they are not evaluated,
% through discretized equations. We have only (Nx-1, Ny-1) nodes 
%(2:Nx, 2:Ny). For this compuitational cells the marginal coefficients 
% are modified due to the BC. This leaves us with only U(3:Nx-1, 3:Ny-1). 
% Since AE and AW are not changed for top and bottom boundary cells 
% ==> AE and AW have the size of (3:Nx-1,2:Ny), 
% for AN and AS ==> (2:Nx, 3:Ny-1). 

AE= -( De(3:imax-1,2:jmax) - (1-alphae(3:imax-1,2:jmax)) ...
    .*Fe(3:imax-1,2:jmax) );

AW= -( Dw(3:imax-1,2:jmax) + alphaw(3:imax-1,2:jmax) ...
    .*Fw(3:imax-1,2:jmax) );

AN= -( Dn(2:imax,3:jmax-1) - (1-alphan(2:imax,3:jmax-1)) ...
    .*Fn(2:imax,3:jmax-1) );

AS= -( Ds(2:imax,3:jmax-1) + alphas(2:imax,3:jmax-1) ...
    .*Fs(2:imax,3:jmax-1) );

%% West Boundary i=2
AW_w=zeros(1,size(AW,2));

AE_w = -( De(2,2:jmax) - (1-alphae(2,2:jmax) ) .*Fe(2,2:jmax) - ...
    D1_a(1,2:jmax) );

S_p_w = Fw(2,2:jmax)+D2_a(1,2:jmax);


%% East Boundary i=Nx-1

AE_e=zeros(1,size(AE,2));

AW_e= -( Dw(imax,2:jmax) + alphaw(imax,2:jmax)...
    .*Fw(imax,2:jmax) );



%% South Boundary j=2
AN_s= -( Dn(2:imax,2) - (1- alphan(2:imax,2) ) .*Fn(2:imax,2) );



%% North Boundary j=Ny
AS_n= -( Ds(2:imax,jmax) + alphas(2:imax,jmax)...
    .*Fs(2:imax,jmax) );

    
S_p(2,2:jmax)   =   S_p(2,2:jmax)+S_p_w;
S_p(imax,2:jmax)  =   S_p(imax,2:jmax)+S_p_e;
S_p(2:imax,2) =   S_p(2:imax,2)+S_p_s;
S_p(2:imax,jmax)  =   S_p(2:imax,jmax)+S_p_n;


%% 
AW=[AW_w; AW; AW_e];
AE=[AE_w; AE; AE_e];
AS=[AS AS_n];
AN=[AN_s AN];

aw_p(2:imax,2:jmax)=AW;
ae_p(2:imax,2:jmax)=AE;
as_p(2:imax,3:jmax)=AS;
an_p(2:imax,2:jmax-1)=AN;


ap_p(2:imax,2:jmax)= (-( aw_p(2:imax,2:jmax) +ae_p(2:imax,2:jmax) ...
        +as_p(2:imax,2:jmax) +an_p(2:imax,2:jmax) )+ ... 
        A0_p(2:imax,2:jmax) + dF(2:imax,2:jmax) + ...
        S_p(2:imax,2:jmax))/alpha_q;

    
    
%% Immersed Boundary Treating

%****************************
% Body Cells

for i=1:size(I_solid,2)
%     change values for An As Ae Aw Ap and RHS.
    
 
    ae_p(I_solid(i),J_solid(i))=0;
      
 
    aw_p(I_solid(i),J_solid(i))=0;
      
    
    an_p(I_solid(i),J_solid(i))=0;
      
    
    as_p(I_solid(i),J_solid(i))=0;
      
    
    ap_p(I_solid(i),J_solid(i))=1;
end

%****************************    
% Ghost Cells
for i=1:length(I_g)
    
    ind=(J_g(i)-2)*(imax-1)+I_g(i)-1;
    % FIND Ghost VIRTUAL POINT LOCATION
    
        
    A1_I_g(1,i)=ind;
    A1_J_g(1,i)=(J1(i)-2)*(imax-1)+I1(i)-1;
    A1_g(1,i)=-landa_g_1(1,i);
        
    
    A2_I_g(1,i)=ind;
    A2_J_g(1,i)=(J2(i)-2)*(imax-1)+I2(i)-1;
    A2_g(1,i)=-landa_g_2(1,i);
        

    
    A3_I_g(1,i)=ind;
    A3_J_g(1,i)=(J3(i)-2)*(imax-1)+I3(i)-1;
    A3_g(1,i)=-landa_g_3(1,i);
        
    
    A4_I_g(1,i)=ind;
    A4_J_g(1,i)=(J4(i)-2)*(imax-1)+I4(i)-1;
    A4_g(1,i)=-landa_g_4(1,i);
    
    
    if IBM.BQv == 1
    
        A5_I_g(1,i)=ind;
        A5_J_g(1,i)=(J5(i)-2)*(imax-1)+I5(i)-1;
        A5_g(1,i)=-landa_g_5(1,i);

    
        A6_I_g(1,i)=ind;
        A6_J_g(1,i)=(J6(i)-2)*(imax-1)+I6(i)-1;
        A6_g(1,i)=-landa_g_6(1,i);
        
    end
        
        
   
    ae_p(I_g(i),J_g(i))=0;
      
    
    aw_p(I_g(i),J_g(i))=0;
      
    
    an_p(I_g(i),J_g(i))=0;
      
    
    as_p(I_g(i),J_g(i))=0;
      
    
    ap_p(I_g(i),J_g(i))=1;
                
               
    
end


%% Form Vectors
AW=aw_p(2:imax,2:jmax);
AE=ae_p(2:imax,2:jmax);
AS=as_p(2:imax,3:jmax);
AN=an_p(2:imax,2:jmax-1);
AP=ap_p(2:imax,2:jmax);


AW=reshape(AW,[1,numel(AW)]);
AW=AW(2:end);

AE=reshape(AE,[1,numel(AE)]);
AE=AE(1:end-1);

AS=reshape(AS,[1,numel(AS)]);

AN=reshape(AN,[1,numel(AS)]);

AP=reshape(AP,[1,numel(AP)]);

if IBM.BQu == 1

    A_g_sparse= sparse(A1_I_g,A1_J_g,A1_g,L,L)+...
                sparse(A2_I_g,A2_J_g,A2_g,L,L)+...
                sparse(A3_I_g,A3_J_g,A3_g,L,L)+...
                sparse(A4_I_g,A4_J_g,A4_g,L,L)+...
                sparse(A5_I_g,A5_J_g,A5_g,L,L)+...
                sparse(A6_I_g,A6_J_g,A6_g,L,L);
            
elseif IBM.BQu == 0
    
    A_g_sparse= sparse(A1_I_g,A1_J_g,A1_g,L,L)+...
                sparse(A2_I_g,A2_J_g,A2_g,L,L)+...
                sparse(A3_I_g,A3_J_g,A3_g,L,L)+...
                sparse(A4_I_g,A4_J_g,A4_g,L,L);
            
end



CM=sparse(AE_I,AE_J,AE,L,L)+sparse(AW_I,AW_J,AW,L,L)+...
  sparse(AN_I,AN_J,AN,L,L)+sparse(AS_I,AS_J,AS,L,L)+...
  sparse(AP_I,AP_J,AP,L,L)+A_g_sparse;
