function [CM_v,S,aw_v,ae_v,as_v,an_v,ap_v,d_v] = COEFFV1(V_old,P,...
        V_star_old,DOMAIN,Diffu_V,ConvF_V,IBM_coeffV,IBM,VARIABLES,BC,...
        disc_scheme)

        
%         
%         V_old,P,...
%             De,Dw,Dn,Ds,Fe,Fw,Fn,Fs,dF,A0_p,...
%             landa_g_1,landa_g_2,landa_g_3,landa_g_4,...
%             I_solid,J_solid,I_g,J_g,I_e,J_e,RHS_V_g,numg,...
%             alpha_v,phi_inside,V_star_old)
% CoNSu 
% CoNSv
% dx ==> dxu

%  BC=1 ==>     Dirichlet
%  BC=2 ==>     Constant total flux
%  BC=3 ==>     % zero diffusive
% global CoNSu V_a V_b V_c V_d BC_e_v BC_w_v BC_n_v BC_s_v dxu imax jmax

%% Retrieve variables from structs

imax = DOMAIN.imax;
jmax = DOMAIN.jmax;
dxu = DOMAIN.dxu;
CoNSu = DOMAIN.CoNSu;
CoNSv = DOMAIN.CoNSv;

De=Diffu_V.De_v;
Dw=Diffu_V.Dw_v;
Dn=Diffu_V.Dn_v;
Ds=Diffu_V.Ds_v;

Fe=ConvF_V.Fe_v;
Fw=ConvF_V.Fw_v;
Fn=ConvF_V.Fn_v;
Fs=ConvF_V.Fs_v;
A0_p=ConvF_V.A0_p_v;
dF=ConvF_V.dF_v;

V_a = BC.V_a;
V_b = BC.V_b;
V_c = BC.V_c;
V_d = BC.V_d;

BC_e_v = BC.BC_e_v;
BC_w_v = BC.BC_w_v;
BC_n_v = BC.BC_n_v;
BC_s_v = BC.BC_s_v;

landa_g_1=IBM_coeffV.landa_g_1_v;
landa_g_2=IBM_coeffV.landa_g_2_v;
landa_g_3=IBM_coeffV.landa_g_3_v;
landa_g_4=IBM_coeffV.landa_g_3_v;


I_solid = IBM_coeffV.I_solid_v;
J_solid = IBM_coeffV.J_solid_v;
I_g = IBM_coeffV.I_g_v;
J_g = IBM_coeffV.J_g_v;
% I_e = IBM_coeffV.I_e_v;
% J_e = IBM_coeffV.J_e_v;

I1=IBM_coeffV.I1_v;
I2=IBM_coeffV.I2_v;
I3=IBM_coeffV.I3_v;
I4=IBM_coeffV.I4_v;

J1=IBM_coeffV.J1_v;
J2=IBM_coeffV.J2_v;
J3=IBM_coeffV.J3_v;
J4=IBM_coeffV.J4_v;


RHS_V_g = IBM_coeffV.A1_g_v;
numg = IBM_coeffV.numg_v;
phi_inside = IBM.phi_inside_u;


alpha_v = VARIABLES.alpha_v;
%% Define parameters


L    = (imax-1)*(jmax-2); % Length of the matrix
jump = imax-1;
coef=1-alpha_v;
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



if IBM.BQv == 1 
    
    landa_g_5=IBM_coeffV.landa_g_5_v;
    landa_g_6=IBM_coeffV.landa_g_6_v;
    
    I5=IBM_coeffV.I5_v;
    I6=IBM_coeffV.I6_v;
 
    J5=IBM_coeffV.J5_v;
    J6=IBM_coeffV.J6_v;

    A5_I_g=zeros(1,numg);
    A6_I_g=zeros(1,numg);

    A5_J_g=zeros(1,numg);
    A6_J_g=zeros(1,numg);
    
    A5_g=zeros(1,numg);
    A6_g=zeros(1,numg);
    
end
%% Matrices
aw_v=zeros(imax+1,jmax);
ae_v=zeros(imax+1,jmax);
as_v=zeros(imax+1,jmax);
an_v=zeros(imax+1,jmax);
ap_v=zeros(imax+1,jmax);
d_v=zeros(imax+1,jmax);
S_p=zeros(imax+1,jmax);

S0=zeros(imax+1,jmax);
S_Pres=zeros(imax+1,jmax);
S_ur=zeros(imax+1,jmax);
S_u=zeros(imax+1,jmax);

S_p_w=zeros(1,jmax-2);
S_p_e=zeros(1,jmax-2);

S_p_s=zeros(imax-1,1);
S_p_n=zeros(imax-1,1);

S_u_w=zeros(1,jmax-2);
S_u_e=zeros(1,jmax-2);

S_u_s=zeros(imax-1,1);
S_u_n=zeros(imax-1,1);

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


if disc_scheme == 1 || disc_scheme == 3
    alphae = ConvF_V.alphae_v;
    alphaw = ConvF_V.alphaw_v;
    alphan = ConvF_V.alphan_v;
    alphas = ConvF_V.alphas_v;
%% Inner Points: 

% We have (Nx+1,Ny) V-cells. All other points on the margin are ghost
% points or we have their values from BC. Therefore they are not evaluated,
% through discretized equations. We have only (Nx-1, Ny-2) nodes 
%(2:Nx, 2:Ny-1). For this compuitational cells the marginal coefficients 
% are modified due to the BC. This leaves us with only U(3:Nx-1, 3:Ny-2). 
% Since AE and AW are not changed for top and bottom boundary cells 
% ==> AE and AW have the size of (3:Nx-1,2:Ny-1), 
% for AN and AS ==> (2:Nx, 3:Ny-2). 

    AE= -( De(3:imax-1,2:jmax-1) - (1-alphae(3:imax-1,2:jmax-1)) ...
        .*Fe(3:imax-1,2:jmax-1) );
    AW= -( Dw(3:imax-1,2:jmax-1) + alphaw(3:imax-1,2:jmax-1) ...
        .*Fw(3:imax-1,2:jmax-1) );

    AN= -( Dn(2:imax,3:jmax-2) - (1-alphan(2:imax,3:jmax-2)) ...
        .*Fn(2:imax,3:jmax-2) );
    AS= -( Ds(2:imax,3:jmax-2) + alphas(2:imax,3:jmax-2) ...
        .*Fs(2:imax,3:jmax-2) );


%% West Boundary i=2
    AW_w=zeros(1,size(AW,2));

    if BC_w_v==1          % 1 for Dirichlet, 3 for Neumann
        AE_w = -( De(2,2:jmax-1) - (1-alphae(2,2:jmax-1) ) .*Fe(2,2:jmax-1) );
%     S_p(2,2:jmax-1)= Dw(2,2:jmax-1)./(1- CoNSu(2,2:jmax-1))  + Fw(2,2:jmax-1);
%     S_p_w= Dw(2,2:jmax-1)./(1- CoNSu(2,2:jmax-1))  + Fw(2,2:jmax-1);
        S_p_w= Dw(2,2:jmax-1)./(1-CoNSu(1,2:jmax-1)) + Fw(2,2:jmax-1);
    elseif BC_w_v==3
        AE_w = -( De(2,2:jmax-1) - (1-alphae(2,2:jmax-1) ) .*Fe(2,2:jmax-1) );
    %S_p(2,2:Ny)= 0
    end


%% East Boundary i=Nx-1

    AE_e=zeros(1,size(AE,2));

    if BC_e_v==1          % 1 for Dirichlet, 3 for Neumann
        AW_e= -( Dw(imax,2:jmax-1) + alphaw(imax,2:jmax-1) .*Fw(imax,2:jmax-1) );
%     S_p(imax,2:jmax-1)=  De(imax,2:jmax-1)./ ...
%         CoNSu(imax,2:jmax-1) + Fe(imax,2:jmax-1) ;
        S_p_e=  De(imax,2:jmax-1)./...
            CoNSu(imax,2:jmax-1) + Fe(imax,2:jmax-1) ;
    
    elseif BC_e_v==3
        AW_e= -( Dw(imax,2:jmax-1) + alphaw(imax,2:jmax-1) .*Fw(imax,2:jmax-1) );
    % S_p(Nx-1,2:Ny)= 0;
    end


%% South Boundary j=2
    if BC_s_v==1          % 1 for Dirichlet, 3 for Neumann
        AN_s= -( Dn(2:imax,2) - (1- alphan(2:imax,2) ) .*Fn(2:imax,2) );
%     S_p(2:imax,2)= Ds(2:imax,2) + alphas(2:imax,2)  .*Fs(2:imax,2);
        S_p_s= Ds(2:imax,2) + alphas(2:imax,2)  .*Fs(2:imax,2);

    elseif BC_s_v==3
        AN_s= -( Dn(2:imax,2) - (1- alphan(2:imax,2) ) .*Fn(2:imax,2) );
    % S_p(2:Nx-1,2)= 0;
    end

%% North Boundary j=Ny
    if BC_n_v==1          % 1 for Dirichlet, 3 for Neumann
        AS_n= -( Ds(2:imax,jmax-1) + alphas(2:imax,jmax-1) .*Fs(2:imax,jmax-1) );
%     S_p(2:imax,jmax-1)= ( Dn(2:imax,jmax-1) ...
%         - (1-alphan(2:imax,jmax-1)) .*Fn(2:imax,jmax-1) );
        S_p_n= ( Dn(2:imax,jmax-1) ...
            - (1-alphan(2:imax,jmax-1)) .*Fn(2:imax,jmax-1) );
    elseif BC_n_v==3
        AS_n=  -( Ds(2:imax,jmax-1) + alphas(2:imax,jmax-1) .*Fs(2:imax,jmax-1) );
    % S_p(2:Nx-1,Ny)=0;
    end

    
S_p(2,2:jmax-1) = S_p(2,2:jmax-1)+S_p_w;
S_p(imax,2:jmax-1) = S_p(imax,2:jmax-1)+S_p_e;
S_p(2:imax,2) = S_p(2:imax,2)+S_p_s;
S_p(2:imax,jmax-1) = S_p(2:imax,jmax-1)+S_p_n;

elseif disc_scheme ==2
    %% Inner Points: 

% We have (Nx+1,Ny) V-cells. All other points on the margin are ghost
% points or we have their values from BC. Therefore they are not evaluated,
% through discretized equations. We have only (Nx-1, Ny-2) nodes 
%(2:Nx, 2:Ny-1). For this compuitational cells the marginal coefficients 
% are modified due to the BC. This leaves us with only U(3:Nx-1, 3:Ny-2). 
% Since AE and AW are not changed for top and bottom boundary cells 
% ==> AE and AW have the size of (3:Nx-1,2:Ny-1), 
% for AN and AS ==> (2:Nx, 3:Ny-2). 
    
    AE= -( De(3:imax-1,2:jmax-1) - ...
        CoNSu(3:imax-1,2:jmax-1) .*Fe(3:imax-1,2:jmax-1) );
    
    AW= -( Dw(3:imax-1,2:jmax-1) + ... 
    (1- CoNSu(2:imax-2,2:jmax-1)).*Fw(3:imax-1,2:jmax-1) );
    
    
    AN= -( Dn(2:imax,3:jmax-2) - ...
        (CoNSv(2:imax,3:jmax-2)).*Fn(2:imax,3:jmax-2) );
    
    AS= -( Ds(2:imax,3:jmax-2) + ...
        (1- CoNSv(2:imax,2:jmax-3)).*Fs(2:imax,3:jmax-2) );
    
    

%% West Boundary i=2
    AW_w=zeros(1,size(AW,2));

    if BC_w_v==1          % 1 for Dirichlet, 3 for Neumann
        AE_w = -( De(2,2:jmax-1) - CoNSu(2,2:jmax-1) .*Fe(2,2:jmax-1) );
        S_p_w= Dw(2,2:jmax-1)./(1-CoNSu(1,2:jmax-1)) + Fw(2,2:jmax-1);
    elseif BC_w_v==3
        AE_w = -( De(2,2:jmax-1) - CoNSu(2,2:jmax-1) .*Fe(2,2:jmax-1) );
    %S_p(2,2:Ny)= 0
    end


%% East Boundary i=Nx-1

    AE_e=zeros(1,size(AE,2));

    if BC_e_v==1          % 1 for Dirichlet, 3 for Neumann
        AW_e= -( Dw(imax,2:jmax-1) + ...
            (1- CoNSu(imax-1,2:jmax-1)) .*Fw(imax,2:jmax-1) );
%     S_p(imax,2:jmax-1)=  De(imax,2:jmax-1)./ ...
%         CoNSu(imax,2:jmax-1) + Fe(imax,2:jmax-1) ;
        S_p_e=  De(imax,2:jmax-1)./...
            CoNSu(imax,2:jmax-1) + Fe(imax,2:jmax-1) ;
    
    elseif BC_e_v==3
        AW_e= -( Dw(imax,2:jmax-1) + ...
            (1- CoNSu(imax-1,2:jmax-1)) .*Fw(imax,2:jmax-1) );
    % S_p(Nx-1,2:Ny)= 0;
    end

  
%% South Boundary j=2
    if BC_s_v==1          % 1 for Dirichlet, 3 for Neumann
        AN_s= -( Dn(2:imax,2) - (CoNSv(2:imax,2)).*Fn(2:imax,2) );
%     S_p(2:imax,2)= Ds(2:imax,2) + alphas(2:imax,2)  .*Fs(2:imax,2);
        S_p_s= -( Ds(2:imax,2) + ...
        (1- CoNSv(2:imax,1)).*Fs(2:imax,2) );

    elseif BC_s_v==3
        AN_s= -( Dn(2:imax,2) - (CoNSv(2:imax,2)).*Fn(2:imax,2) );
    % S_p(2:Nx-1,2)= 0;
    end

    
    
%% North Boundary j=Ny
    if BC_n_v==1          % 1 for Dirichlet, 3 for Neumann
        AS_n= -( Ds(2:imax,jmax-1) + ...
            (1- CoNSv(2:imax,jmax-2)).*Fs(2:imax,jmax-1) );
%     S_p(2:imax,jmax-1)= ( Dn(2:imax,jmax-1) ...
%         - (1-alphan(2:imax,jmax-1)) .*Fn(2:imax,jmax-1) );
        S_p_n= ( Dn(2:imax,jmax-1) ...
            - (CoNSv(2:imax,jmax-1)) .*Fn(2:imax,jmax-1) );
    elseif BC_n_v==3
        AS_n= -( Ds(2:imax,jmax-1) + ...
            (1- CoNSv(2:imax,jmax-2)).*Fs(2:imax,jmax-1) );
    % S_p(2:Nx-1,Ny)=0;
    end
    
    S_p(2,2:jmax-1) = S_p(2,2:jmax-1)+S_p_w;
    S_p(imax,2:jmax-1) = S_p(imax,2:jmax-1)+S_p_e;
    S_p(2:imax,2) = S_p(2:imax,2)+S_p_s;
    S_p(2:imax,jmax-1) = S_p(2:imax,jmax-1)+S_p_n;

end
    
S_p(2,2:jmax-1) = S_p(2,2:jmax-1)+S_p_w;
S_p(imax,2:jmax-1) = S_p(imax,2:jmax-1)+S_p_e;
S_p(2:imax,2) = S_p(2:imax,2)+S_p_s;
S_p(2:imax,jmax-1) = S_p(2:imax,jmax-1)+S_p_n;
%% 
AW=[AW_w; AW; AW_e];
AE=[AE_w; AE; AE_e];
AS=[AS AS_n];
AN=[AN_s AN];

aw_v(2:imax,2:jmax-1)=AW;
ae_v(2:imax,2:jmax-1)=AE;
as_v(2:imax,3:jmax-1)=AS;
an_v(2:imax,2:jmax-2)=AN;

ap_v(2:imax,2:jmax-1)= (-( aw_v(2:imax,2:jmax-1) +ae_v(2:imax,2:jmax-1) ...
                +as_v(2:imax,2:jmax-1) +an_v(2:imax,2:jmax-1) )+ ... 
                A0_p(2:imax,2:jmax-1) + dF(2:imax,2:jmax-1) + ...
                S_p(2:imax,2:jmax-1))/alpha_v;
            
d_v(2:imax,2:jmax-1)=dxu(1:imax-1,2:jmax-1) ./...
    ( ap_v(2:imax,2:jmax-1)+ aw_v(2:imax,2:jmax-1)+ ae_v(2:imax,2:jmax-1)+...
    as_v(2:imax,2:jmax-1) + an_v(2:imax,2:jmax-1) );


%% ============ Source Term Discretization ================================
% This terms only evolve in time loop not the QUICK defered-correction loop
% S0 ==> is the contribution from time value
% S_Pres ==> the pressure gradient discretization
% S_ur ==> under-relaxed value
% S_u ==> Contribution from boundary conditions
%==========================================================================
S0(2:imax,2:jmax-1) = A0_p(2:imax,2:jmax-1) .*V_old(2:imax,2:jmax-1);
S_Pres(2:imax,2:jmax-1) = ( P(2:imax,2:jmax-1)-P(2:imax,3:jmax) ) .*dxu(1:imax-1,2:jmax-1);
S_ur(2:imax,2:jmax-1) = ap_v(2:imax,2:jmax-1)* coef.* V_star_old(2:imax,2:jmax-1);

% ==== Boundary Condition Contribution
%WEST   i=2
if BC_w_v==1
    S_u_w= V_a(:,2:jmax-1).*...
        S_p_w;
elseif BC_w_v==3
    S_u_w=0;
end

%EAST   i=Nx
if BC_e_v==1
    S_u_e=V_b(:,2:jmax-1).*...
        S_p_e;
elseif BC_e_v==3
    S_u_e=0;
end

%SOUTH   j=2
if BC_s_v==1
    S_u_s= S_p_s.*V_c(:,2:imax);
elseif BC_s_v==3
    S_u_s=0;
end

%NORTH   j=Ny-1
if BC_s_v==1
    S_u_n=S_p_n.*V_d(:,2:imax);
elseif BC_s_v==3
    S_u_n=0;
end
%=====


S_u(2,2:jmax-1)   =   S_u(2,2:jmax-1)+S_u_w;
S_u(imax,2:jmax-1)  =   S_u(imax,2:jmax-1)+S_u_e;
S_u(2:imax,2) =   S_u(2:imax,2)+S_u_s;
S_u(2:imax,jmax-1)  =   S_u(2:imax,jmax-1)+S_u_n;


% AP=ap_v(2:imax,2:jmax-1);



%% Immersed Boundary Treating

% Body Cells

%****************************
for i=1:size(I_solid,2)
%     change values for An As Ae Aw Ap and RHS.
    
 
    ae_v(I_solid(i),J_solid(i))=0;
      
 
    aw_v(I_solid(i),J_solid(i))=0;
      
    
    an_v(I_solid(i),J_solid(i))=0;
      
    
    as_v(I_solid(i),J_solid(i))=0;
      
    
    ap_v(I_solid(i),J_solid(i))=1;
    d_v(I_solid(i),J_solid(i))=0;  %% I Changed Here
%     ap_v(I_solid(i),J_solid(i))=AP(I_solid(i),J_solid(i));
    
    S0(I_solid(i),J_solid(i))=phi_inside;
    S_Pres(I_solid(i),J_solid(i))=0;
    S_ur(I_solid(i),J_solid(i))=0;
end

%****************************
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
        
        
   
    ae_v(I_g(i),J_g(i))=0;
      
    
    aw_v(I_g(i),J_g(i))=0;
      
    
    an_v(I_g(i),J_g(i))=0;
      
    
    as_v(I_g(i),J_g(i))=0;
      
    
    ap_v(I_g(i),J_g(i))=1;
                
               
    S0(I_g(i),J_g(i))=RHS_V_g(1,i);
    S_Pres(I_g(i),J_g(i))=0;
    S_ur(I_g(i),J_g(i))=0;
    S_u(I_g(i),J_g(i))=0;
    
    d_v(I_g(i),J_g(i))=0;
%     ap_v(I_g(i),J_g(i))=AP(I_g(i),J_g(i));

%     change values for An As Ae Aw Ap A1_g A2_g A3_g A4_g and its RHS
end
%% Form Vectors
AW=aw_v(2:imax,2:jmax-1);
AE=ae_v(2:imax,2:jmax-1);
AS=as_v(2:imax,3:jmax-1);
AN=an_v(2:imax,2:jmax-2);
AP=ap_v(2:imax,2:jmax-1);


AW=reshape(AW,[1,numel(AW)]);
AW=AW(2:end);


AE=reshape(AE,[1,numel(AE)]);
AE=AE(1:end-1);

AS=reshape(AS,[1,numel(AS)]);

AN=reshape(AN,[1,numel(AS)]);

AP=reshape(AP,[1,numel(AP)]);


S=S0+S_Pres+S_ur+S_u;

S= S(2:imax,2:jmax-1);    

S=reshape(S,[numel(S),1]);

if IBM.BQv == 1

    A_g_sparse= sparse(A1_I_g,A1_J_g,A1_g,L,L)+...
            sparse(A2_I_g,A2_J_g,A2_g,L,L)+...
            sparse(A3_I_g,A3_J_g,A3_g,L,L)+...
            sparse(A4_I_g,A4_J_g,A4_g,L,L)+...
            sparse(A5_I_g,A5_J_g,A5_g,L,L)+...
            sparse(A6_I_g,A6_J_g,A6_g,L,L);
        
elseif IBM.BQv == 0
     A_g_sparse= sparse(A1_I_g,A1_J_g,A1_g,L,L)+...
            sparse(A2_I_g,A2_J_g,A2_g,L,L)+...
            sparse(A3_I_g,A3_J_g,A3_g,L,L)+...
            sparse(A4_I_g,A4_J_g,A4_g,L,L);
    
end
            
CM_v=sparse(AE_I,AE_J,AE,L,L)+sparse(AW_I,AW_J,AW,L,L)+...
  sparse(AN_I,AN_J,AN,L,L)+sparse(AS_I,AS_J,AS,L,L)+...
  sparse(AP_I,AP_J,AP,L,L)+A_g_sparse;



