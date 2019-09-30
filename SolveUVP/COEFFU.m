function [CM_u,S,aw_u,ae_u,as_u,an_u,ap_u,d_u,A_g_sparse,M,d_u_sec] = COEFFU(U_old,P,...
    U_star_old,DOMAIN,Diffu_U,ConvF_U,IBM_coeffU,IBM,VARIABLES,BC,...
    disc_scheme)
   

%% Retrieve variables from structs
 
imax = DOMAIN.imax;
jmax = DOMAIN.jmax;
dyv = DOMAIN.dyv;
CoEWv = DOMAIN.CoEWv;
CoEWu = DOMAIN.CoEWu;
 
De=Diffu_U.De_u;
Dw=Diffu_U.Dw_u;
Dn=Diffu_U.Dn_u;
Ds=Diffu_U.Ds_u;

Fe=ConvF_U.Fe_u;
Fw=ConvF_U.Fw_u;
Fn=ConvF_U.Fn_u;
Fs=ConvF_U.Fs_u;
A0_p=ConvF_U.A0_p_u;
dF=ConvF_U.dF_u;

U_a = BC.U_a;
U_b = BC.U_b;
U_c = BC.U_c;
U_d = BC.U_d;

BC_e_u = BC.BC_e_u;
BC_w_u = BC.BC_w_u;
BC_n_u = BC.BC_n_u;
BC_s_u = BC.BC_s_u;

BC_e_p = BC.BC_e_p;


landa_g_1=IBM_coeffU.landa_g_1_u;
landa_g_2=IBM_coeffU.landa_g_2_u;
landa_g_3=IBM_coeffU.landa_g_3_u;
landa_g_4=IBM_coeffU.landa_g_3_u;


I_solid = IBM_coeffU.I_solid_u;
J_solid = IBM_coeffU.J_solid_u;
I_g = IBM_coeffU.I_g_u;
J_g = IBM_coeffU.J_g_u;
% I_e = IBM_coeffU.I_e_u;
% J_e = IBM_coeffU.J_e_u;

I1=IBM_coeffU.I1_u;
I2=IBM_coeffU.I2_u;
I3=IBM_coeffU.I3_u;
I4=IBM_coeffU.I4_u;

J1=IBM_coeffU.J1_u;
J2=IBM_coeffU.J2_u;
J3=IBM_coeffU.J3_u;
J4=IBM_coeffU.J4_u;


RHS_U_g = IBM_coeffU.A1_g_u;
numg = IBM_coeffU.numg_u;
phi_inside = IBM.phi_inside_u;


alpha_u = VARIABLES.alpha_u;

        
% CoEWu 
% CoEWv
% dy ==> dyv

        
%  BC=1 ==>     Dirichlet
%  BC=2 ==>     Constant total flux
%  BC=3 ==>     % zero diffusive




%% Define parameters
% get no of grids in every direction


L    = (imax-2)*(jmax-1); % size of the system


% distance between three center bands and two farther bands 
jump = imax-2;
if BC_e_p == 1
    L    = (imax-1)*(jmax-1);
    jump = imax-1;

end

coef=1-alpha_u;
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


if IBM.BQu == 1
    
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
aw_u=zeros(imax,jmax+1);
ae_u=zeros(imax,jmax+1);
as_u=zeros(imax,jmax+1);
an_u=zeros(imax,jmax+1);
ap_u=zeros(imax,jmax+1);
d_u=zeros(imax,jmax+1);
d_u_sec = zeros(imax,jmax+1);
S_p=zeros(imax,jmax+1);

if BC_e_p == 1
    aw_u=zeros(imax+1,jmax+1);
    ae_u=zeros(imax+1,jmax+1);
    as_u=zeros(imax+1,jmax+1);
    an_u=zeros(imax+1,jmax+1);
    ap_u=zeros(imax+1,jmax+1);
end 

S0=zeros(imax,jmax+1);
S_Pres=zeros(imax,jmax+1);
S_ur=zeros(imax,jmax+1);
S_u=zeros(imax,jmax+1);


S_p_w=zeros(1,jmax-1);
S_p_e=zeros(1,jmax-1);

S_p_s=zeros(imax-2,1);
S_p_n=zeros(imax-2,1);

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
    
    alphae = ConvF_U.alphae_u;
    alphaw = ConvF_U.alphaw_u;
    alphan = ConvF_U.alphan_u;
    alphas = ConvF_U.alphas_u;
%% Inner Points: 

% We have (Nx,Ny+1) U-cells. All other points on the margin are ghost
% points or we have their values from BC. Therefore they are not evaluated,
% through discretized equations. We have only (Nx-2, Ny-1) nodes 
%(2:Nx-1, 2:Ny). For this compuitational cells the marginal coefficients 
% are modified due to the BC. This leaves us with only U(3:Nx-2, 3:Ny-1). 
% Since AE and AW are not changed for top and bottom boundary cells 
% ==> AE and AW have the size of (3:Nx-2,2:Ny), 
% for AN and AS ==> (2:Nx-1, 3:Ny-1). 

    
    
    if BC_e_p == 1
    
        AE= -( De(2:imax-1,2:jmax) - (1-alphae(2:imax-1,2:jmax)) ...
            .*Fe(2:imax-1,2:jmax) );

        AW= -( Dw(3:imax,2:jmax) + alphaw(3:imax,2:jmax) ...
            .*Fw(3:imax,2:jmax) );

        AN= -( Dn(2:imax,2:jmax-1) - (1-alphan(2:imax,2:jmax-1)) ...
            .*Fn(2:imax,2:jmax-1) );

        AS= -( Ds(2:imax,3:jmax) + alphas(2:imax,3:jmax) ...
            .*Fs(2:imax,3:jmax) );
        
    else
        
        AE= -( De(2:imax-2,2:jmax) - (1-alphae(2:imax-2,2:jmax)) ...
        .*Fe(2:imax-2,2:jmax) );
    
        AW= -( Dw(3:imax-1,2:jmax) + alphaw(3:imax-1,2:jmax) ...
            .*Fw(3:imax-1,2:jmax) );

        AN= -( Dn(2:imax-1,2:jmax-1) - (1-alphan(2:imax-1,2:jmax-1)) ...
            .*Fn(2:imax-1,2:jmax-1) );

        AS= -( Ds(2:imax-1,3:jmax) + alphas(2:imax-1,3:jmax) ...
            .*Fs(2:imax-1,3:jmax) );               
        
    end


%% West Boundary i=2
    AW_w=zeros(1,size(AW,2));

    S_p_w=(Dw(2,2:jmax) + alphaw(2,2:jmax)  .*Fw(2,2:jmax)) * (BC_w_u==1);
    


%% East Boundary i=Nx-1

    AE_e=zeros(1,size(AE,2));
    

    
    
    


    
    
    if BC_e_p == 1
        S_p_e= ( De(imax,2:jmax) ...
            - (1-alphae(imax,2:jmax) ) .*Fe(imax,2:jmax) )...
            * (BC_e_u==1);
        
        S_p_s= (Ds(2:imax,2)./(1-CoEWv(2:imax,2)) + Fs(2:imax,2)) * ...
        (BC_s_u==1);
    
        S_p_n= (Dn(2:imax,jmax)./CoEWv(2:imax,jmax) - Fn(2:imax,jmax))...
            * (BC_n_u==1);
        
        S_p(2,2:jmax)   =   S_p(2,2:jmax)+S_p_w;
        S_p(imax,2:jmax)  =   S_p(imax,2:jmax)+S_p_e;
        S_p(2:imax,2) =   S_p(2:imax,2)+S_p_s;
        S_p(2:imax,jmax)  =   S_p(2:imax,jmax)+S_p_n;

    else
        
        S_p_e= ( De(imax-1,2:jmax) ...
            - (1-alphae(imax-1,2:jmax) ) .*Fe(imax-1,2:jmax) )...
            * (BC_e_u==1);
        
        S_p_s= (Ds(2:imax-1,2)./(1-CoEWv(2:imax-1,2)) + Fs(2:imax-1,2)) * ...
            (BC_s_u==1);
        
        S_p_n= (Dn(2:imax-1,jmax)./CoEWv(2:imax-1,jmax) - Fn(2:imax-1,jmax))...
            * (BC_n_u==1);
        
        S_p(2,2:jmax)   =   S_p(2,2:jmax)+S_p_w;
        S_p(imax-1,2:jmax)  =   S_p(imax-1,2:jmax)+S_p_e;
        S_p(2:imax-1,2) =   S_p(2:imax-1,2)+S_p_s;
        S_p(2:imax-1,jmax)  =   S_p(2:imax-1,jmax)+S_p_n;


        
    end
    
%% 
elseif disc_scheme == 2  %% CENTRAL DIFFERENCE scheme
    %% Inner Points: 

% We have (Nx,Ny+1) U-cells. All other points on the margin are ghost
% points or we have their values from BC. Therefore they are not evaluated,
% through discretized equations. We have only (Nx-2, Ny-1) nodes 
%(2:Nx-1, 2:Ny). For this compuitational cells the marginal coefficients 
% are modified due to the BC. This leaves us with only U(3:Nx-2, 3:Ny-1). 
% Since AE and AW are not changed for top and bottom boundary cells 
% ==> AE and AW have the size of (3:Nx-2,2:Ny), 
% for AN and AS ==> (2:Nx-1, 3:Ny-1). 

    AE= -( De(2:imax-2,2:jmax) - ...
        CoEWu(2:imax-2,2:jmax) .*Fe(2:imax-2,2:jmax) );
    
    AW= -( Dw(3:imax-1,2:jmax) + ... 
    (1- CoEWu(2:imax-2,2:jmax)).*Fw(3:imax-1,2:jmax) );

    AN= -( Dn(2:imax-1,2:jmax-1) - ...
        (CoEWv(2:imax-1,2:jmax-1)).*Fn(2:imax-1,2:jmax-1) );
    
    AS= -( Ds(2:imax-1,3:jmax) + ...
        (1- CoEWv(2:imax-1,2:jmax-1)).*Fs(2:imax-1,3:jmax) );
    
    
    if BC_e_p == 1
    
        AE= -( De(2:imax-1,2:jmax) - ...
            CoEWu(2:imax-1,2:jmax) .*Fe(2:imax-1,2:jmax) );

        AW= -( Dw(3:imax,2:jmax) + ... 
            (1- CoEWu(2:imax-1,2:jmax)).*Fw(3:imax,2:jmax) );

        AN= -( Dn(2:imax,2:jmax-1) - ...
            (CoEWv(2:imax,2:jmax-1)).*Fn(2:imax,2:jmax-1) );

        AS= -( Ds(2:imax,3:jmax) + ...
            (1- CoEWv(2:imax,2:jmax-1)).*Fs(2:imax,3:jmax) );
        
    end


%% West Boundary i=2
    AW_w=zeros(1,size(AW,2));

    S_p_w = (Dw(2,2:jmax) + (1- CoEWu(1,2:jmax))  .*Fw(2,2:jmax))...
        *(BC_w_u==1);
    


%% East Boundary i=Nx-1

    AE_e=zeros(1,size(AE,2));




    



    S_p(2,2:jmax)   =   S_p(2,2:jmax)+S_p_w;
    
    
    if BC_e_p == 1
        S_p_e= ( De(imax,2:jmax) ...
        -  CoEWu(imax-1,2:jmax).*Fe(imax,2:jmax) )*(BC_e_u==1);
        
        S_p_s= (Ds(2:imax,2)./(1-CoEWv(2:imax,2)) + Fs(2:imax,2)) * ...
        (BC_s_u==1);
    
        S_p_n= (Dn(2:imax,jmax)./CoEWv(2:imax,jmax) - Fn(2:imax,jmax))...
            * (BC_n_u==1);
        
        S_p(imax,2:jmax)  =   S_p(imax,2:jmax)+S_p_e;
        S_p(2:imax,2) =   S_p(2:imax,2)+S_p_s;
        S_p(2:imax,jmax)  =   S_p(2:imax,jmax)+S_p_n;
        
    else 
        
        S_p_e= ( De(imax-1,2:jmax) ...
            -  CoEWu(imax-1,2:jmax).*Fe(imax-1,2:jmax) )*(BC_e_u==1);
        
        S_p_s= (Ds(2:imax-1,2)./(1-CoEWv(2:imax-1,1)) + Fs(2:imax-1,2))...
            * (BC_s_u==1);
        
        S_p_n= (Dn(2:imax-1,jmax)./CoEWv(2:imax-1,jmax) - Fn(2:imax-1,jmax))...
            * (BC_n_u==1); 
        
        S_p(imax-1,2:jmax)  =   S_p(imax-1,2:jmax)+S_p_e;
        S_p(2:imax-1,2) =   S_p(2:imax-1,2)+S_p_s;
        S_p(2:imax-1,jmax)  =   S_p(2:imax-1,jmax)+S_p_n;
        

    end
end

AW=[AW_w; AW];
AE=[AE; AE_e];

if BC_e_p == 1
    
    aw_u(2:imax,2:jmax)=AW;
    ae_u(2:imax,2:jmax)=AE;
    as_u(2:imax,3:jmax)=AS;
    an_u(2:imax,2:jmax-1)=AN;
    
    ap_u(2:imax,2:jmax)= (-( aw_u(2:imax,2:jmax) +ae_u(2:imax,2:jmax) ...
        +as_u(2:imax,2:jmax) +an_u(2:imax,2:jmax) )+ ... 
        A0_p(2:imax,2:jmax) + dF(2:imax,2:jmax) + ...
        S_p(2:imax,2:jmax))/alpha_u;
    
    d_u(2:imax,2:jmax)=dyv(2:imax,1:jmax-1) ./...
        ( ap_u(2:imax,2:jmax)+ aw_u(2:imax,2:jmax)+ ae_u(2:imax,2:jmax)+...
        as_u(2:imax,2:jmax) + an_u(2:imax,2:jmax) );

    d_u_sec(2:imax,2:jmax)=dyv(2:imax,1:jmax-1) ./...
         ap_u(2:imax,2:jmax) ;    


       
        
else 
    aw_u(2:imax-1,2:jmax)=AW;
    ae_u(2:imax-1,2:jmax)=AE;
    as_u(2:imax-1,3:jmax)=AS;
    an_u(2:imax-1,2:jmax-1)=AN;
    
    ap_u(2:imax-1,2:jmax)= (-( aw_u(2:imax-1,2:jmax) +ae_u(2:imax-1,2:jmax) ...
        +as_u(2:imax-1,2:jmax) +an_u(2:imax-1,2:jmax) )+ ... 
        A0_p(2:imax-1,2:jmax) + dF(2:imax-1,2:jmax) + ...
        S_p(2:imax-1,2:jmax))/alpha_u;
    
    d_u(2:imax-1,2:jmax)=dyv(2:imax-1,1:jmax-1) ./...
        ( ap_u(2:imax-1,2:jmax)+ aw_u(2:imax-1,2:jmax)+ ae_u(2:imax-1,2:jmax)+...
        as_u(2:imax-1,2:jmax) + an_u(2:imax-1,2:jmax) );

    d_u_sec(2:imax-1,2:jmax)=dyv(2:imax-1,1:jmax-1) ./...
         ap_u(2:imax-1,2:jmax) ;    



end



            
     

%% ============ Source Term Discretization ================================
% This terms only evolve in time loop not the QUICK defered-correction loop
% S0 ==> is the contribution from time value
% S_Pres ==> the pressure gradient discretization
% S_ur ==> under-relaxed value
% S_u ==> Contribution from boundary conditions
%==========================================================================
S0(2:imax,2:jmax) = A0_p(2:imax,2:jmax) .*U_old(2:imax,2:jmax);
S_Pres(2:imax,2:jmax) = ( P(2:imax,2:jmax)-P(3:imax+1,2:jmax) ) .*dyv(2:imax,1:jmax-1);
S_ur(2:imax,2:jmax) = ap_u(2:imax,2:jmax)* coef.* U_star_old(2:imax,2:jmax);

% ==== Boundary Condition Contribution
%WEST  i=2
S_u_w= U_a(:,2:jmax).*S_p_w * (BC_w_u==1);

%EAST  i=Nx-1
S_u_e=U_b(:,2:jmax).*S_p_e * (BC_e_u==1);

S_u(2,2:jmax)   =   S_u(2,2:jmax)+S_u_w;
S_u(imax-1,2:jmax)  =   S_u(imax-1,2:jmax)+S_u_e;


if BC_e_p ~= 1
    %SOUTH  j=2
    S_u_s= U_c(2:imax-1,:).*S_p_s * (BC_s_u==1);

    %NORTH  j=Ny
    S_u_n=U_d(2:imax-1,:).*S_p_n * (BC_n_u==1);
    
    S_u(2:imax-1,2) =   S_u(2:imax-1,2)+S_u_s;
    S_u(2:imax-1,jmax)  =   S_u(2:imax-1,jmax)+S_u_n;

elseif BC_e_p == 1
    %SOUTH  j=2
    S_u_s= U_c(2:imax,:).*S_p_s * (BC_s_u==1);

    %NORTH  j=Ny
    S_u_n=U_d(2:imax,:).*S_p_n * (BC_n_u==1);
    
    S_u(2:imax,2) =   S_u(2:imax,2)+S_u_s;
    S_u(2:imax,jmax)  =   S_u(2:imax,jmax)+S_u_n;

end
    
%=====
% AP=ap_u(2:imax-1,2:jmax);



%% Immersed Boundary Treating

for i=1:size(I_solid,2)
%     change values for An As Ae Aw Ap and RHS.
    
 
    ae_u(I_solid(i),J_solid(i))=0;
      
 
    aw_u(I_solid(i),J_solid(i))=0;
      
    
    an_u(I_solid(i),J_solid(i))=0;
      
    
    as_u(I_solid(i),J_solid(i))=0;
      
    
    ap_u(I_solid(i),J_solid(i))=1;
    d_u(I_solid(i),J_solid(i))=0;  %% I Changed Here
    d_u_sec(I_solid(i),J_solid(i))=0;
%     ap_u(I_solid(i),J_solid(i))=AP(I_solid(i),J_solid(i));
    
    S0(I_solid(i),J_solid(i))=phi_inside;
    S_Pres(I_solid(i),J_solid(i))=0;
    S_ur(I_solid(i),J_solid(i))=0;
    S_u(I_solid(i),J_solid(i))=0;

end

%****************************



for i=1:length(I_g)
    
    if BC_e_p ~= 1

    
        ind=(J_g(i)-2)*(imax-2)+I_g(i)-1;
        % FIND Ghost VIRTUAL POINT LOCATION


        A1_I_g(1,i)=ind;
        A1_J_g(1,i)=(J1(i)-2)*(imax-2)+I1(i)-1;
        A1_g(1,i)=-landa_g_1(1,i);


        A2_I_g(1,i)=ind;
        A2_J_g(1,i)=(J2(i)-2)*(imax-2)+I2(i)-1;
        A2_g(1,i)=-landa_g_2(1,i);



        A3_I_g(1,i)=ind;
        A3_J_g(1,i)=(J3(i)-2)*(imax-2)+I3(i)-1;
        A3_g(1,i)=-landa_g_3(1,i);


        A4_I_g(1,i)=ind;
        A4_J_g(1,i)=(J4(i)-2)*(imax-2)+I4(i)-1;
        A4_g(1,i)=-landa_g_4(1,i);

        if IBM.BQu ==1 
            A5_I_g(1,i)=ind;
            A5_J_g(1,i)=(J5(i)-2)*(imax-2)+I5(i)-1;
            A5_g(1,i)=-landa_g_5(1,i);

            A6_I_g(1,i)=ind;
            A6_J_g(1,i)=(J6(i)-2)*(imax-2)+I6(i)-1;
            A6_g(1,i)=-landa_g_6(1,i);

        end
        
    elseif BC_e_p == 1
        
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

        if IBM.BQu ==1 
            A5_I_g(1,i)=ind;
            A5_J_g(1,i)=(J5(i)-2)*(imax-1)+I5(i)-1;
            A5_g(1,i)=-landa_g_5(1,i);

            A6_I_g(1,i)=ind;
            A6_J_g(1,i)=(J6(i)-2)*(imax-1)+I6(i)-1;
            A6_g(1,i)=-landa_g_6(1,i);

        end
        
    end
        
        
    
    
        
        
   
    ae_u(I_g(i),J_g(i))=0;
      
    
    aw_u(I_g(i),J_g(i))=0;
      
    
    an_u(I_g(i),J_g(i))=0;
      
    
    as_u(I_g(i),J_g(i))=0;
      
    
    ap_u(I_g(i),J_g(i))=1;
                
               
    S0(I_g(i),J_g(i))=RHS_U_g(1,i);
    S_Pres(I_g(i),J_g(i))=0;
    S_ur(I_g(i),J_g(i))=0;
    S_u(I_g(i),J_g(i))=0;
    
    d_u(I_g(i),J_g(i))=0;
    d_u_sec(I_g(i),J_g(i))=0;
    
    
%     ap_u(I_g(i),J_g(i))=AP(I_g(i),J_g(i));

%     change values for An As Ae Aw Ap A1_g A2_g A3_g A4_g and its RHS
end

%% Form Vectors

if BC_e_p ~= 1

    AW=aw_u(2:imax-1,2:jmax);
    AE=ae_u(2:imax-1,2:jmax);
    AS=as_u(2:imax-1,3:jmax);
    AN=an_u(2:imax-1,2:jmax-1);
    AP=ap_u(2:imax-1,2:jmax);
    
elseif BC_e_p == 1
    
    AW=aw_u(2:imax,2:jmax);
    AE=ae_u(2:imax,2:jmax);
    AS=as_u(2:imax,3:jmax);
    AN=an_u(2:imax,2:jmax-1);
    AP=ap_u(2:imax,2:jmax);
    
end

AW=reshape(AW,[1,numel(AW)]);
AW=AW(2:end);

AE=reshape(AE,[1,numel(AE)]);
AE=AE(1:end-1);

AS=reshape(AS,[1,numel(AS)]);

AN=reshape(AN,[1,numel(AS)]);

AP=reshape(AP,[1,numel(AP)]);


S=S0+S_Pres+S_ur+S_u;

if BC_e_p ~= 1

    S= S(2:imax-1,2:jmax); 
    
elseif BC_e_p == 1
    
        S= S(2:imax,2:jmax); 
    
end


S=reshape(S,[numel(S),1]);

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


CM_u=   sparse(AE_I,AE_J,AE,L,L)+sparse(AW_I,AW_J,AW,L,L)+...
            sparse(AN_I,AN_J,AN,L,L)+sparse(AS_I,AS_J,AS,L,L)+...
            sparse(AP_I,AP_J,AP,L,L)+A_g_sparse;

M_g=A1_g;
M_g(:,:)=1-1/20;
M=sparse(A1_I_g,A1_I_g,M_g,L,L);
M=speye(L,L)-M;
