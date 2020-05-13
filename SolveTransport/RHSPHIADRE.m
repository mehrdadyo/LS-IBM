function [RHSSCALAR_c] = RHSPHIADRE(StateVar,BC,...
        Flux,DOMAIN,IBM_coeffP,IBM, VARIABLES, Soln)
% Desciption:
%
%% Define parameters
%% Define parameters
imax = DOMAIN.imax;
jmax = DOMAIN.jmax;
L    = (imax-1)*(jmax-1); % Length of the matrix

alpha_q = VARIABLES.alpha_q;
ap_p = Soln.ap_p;

Fe = Flux.ConvF_phi.Fe_p;
Fw = Flux.ConvF_phi.Fw_p;
Fn = Flux.ConvF_phi.Fn_p;
Fs = Flux.ConvF_phi.Fs_p;

A0_p = Flux.ConvF_phi.A0_p_p;

alphae = Flux.ConvF_phi.alphae_p;
alphaw = Flux.ConvF_phi.alphaw_p;
alphan = Flux.ConvF_phi.alphan_p;
alphas = Flux.ConvF_phi.alphas_p;

D2_a = Flux.Diffu_P.D2_a;

g1_e_p = DOMAIN.g1c_e_p;
g2_e_p = DOMAIN.g2c_e_p;
g1_w_p = DOMAIN.g1c_w_p;
g2_w_p = DOMAIN.g2c_w_p;
g1_n_p = DOMAIN.g1c_n_p;
g2_n_p = DOMAIN.g2c_n_p;
g1_s_p = DOMAIN.g1c_s_p;
g2_s_p = DOMAIN.g2c_s_p;

g1_e_n = DOMAIN.g1c_e_n;
g2_e_n = DOMAIN.g2c_e_n;
g1_w_n = DOMAIN.g1c_w_n;
g2_w_n = DOMAIN.g2c_w_n;
g1_n_n = DOMAIN.g1c_n_n;
g2_n_n = DOMAIN.g2c_n_n;
g1_s_n = DOMAIN.g1c_s_n;
g2_s_n = DOMAIN.g2c_s_n;

I_solid = IBM_coeffP.I_solid_p;
J_solid = IBM_coeffP.J_solid_p;
I_g = IBM_coeffP.I_g_p;
J_g = IBM_coeffP.J_g_p;

phi = StateVar.phi;
phi_old = StateVar.phi_old;

phi_inside = IBM.phi_inside_phi;
phi_a = BC.phi_a;

RHS_phi_g = IBM_coeffP.A1_g_p;
%% Vectors
%AP
S0=zeros(imax+1,jmax+1);
S_e=zeros(imax+1,jmax+1);
S_w=zeros(imax+1,jmax+1);
S_n=zeros(imax+1,jmax+1);
S_s=zeros(imax+1,jmax+1);
S_ur=zeros(imax+1,jmax+1);



%% all others 

S0(2:imax,2:jmax) = A0_p(2:imax,2:jmax) .*phi_old(2:imax,2:jmax);

S_e(2:imax-1,2:jmax) = Fe(2:imax-1,2:jmax).*( alphae(2:imax-1,2:jmax).* ...
    ( phi(2:imax-1,2:jmax).* ...
    ( g1_e_p(2:imax-1,2:jmax)-g2_e_p(2:imax-1,2:jmax) ) + ...
    phi(1:imax-2,2:jmax).*g2_e_p(2:imax-1,2:jmax) - ...
    phi(3:imax,2:jmax).*g1_e_p(2:imax-1,2:jmax) )...
    +(1-alphae(2:imax-1,2:jmax)).* ( phi(3:imax,2:jmax).* ...
    ( g1_e_n(2:imax-1,2:jmax)-g2_e_n(2:imax-1,2:jmax) ) +...
    phi(4:imax+1,2:jmax).*g2_e_n(2:imax-1,2:jmax)...
    - phi(2:imax-1,2:jmax).*g1_e_n(2:imax-1,2:jmax) ) );
       
S_w(3:imax,2:jmax)=Fw(3:imax,2:jmax).*( alphaw(3:imax,2:jmax).* ...
    ( phi(2:imax-1,2:jmax).* ...
    ( g2_w_p(3:imax,2:jmax)-g1_w_p(3:imax,2:jmax) ) + ...
    phi(3:imax,2:jmax).*g1_w_p(3:imax,2:jmax) - ...
    phi(1:imax-2,2:jmax).*g2_w_p(3:imax,2:jmax) )...
    +(1-alphaw(3:imax,2:jmax)).* ( phi(3:imax,2:jmax).*...
    ( g2_w_n(3:imax,2:jmax)-g1_w_n(3:imax,2:jmax) ) +...
    phi(2:imax-1,2:jmax).*g1_w_n(3:imax,2:jmax) - ...
    phi(4:imax+1,2:jmax).*g2_w_n(3:imax,2:jmax) ) );

S_n(2:imax,2:jmax-1)=Fn(2:imax,2:jmax-1).*( alphan(2:imax,2:jmax-1).* ...
    ( phi(2:imax,2:jmax-1).* ...
    ( g1_n_p(2:imax,2:jmax-1)-g2_n_p(2:imax,2:jmax-1) ) + ...
    phi(2:imax,1:jmax-2).*g2_n_p(2:imax,2:jmax-1) - ...
    phi(2:imax,3:jmax).*g1_n_p(2:imax,2:jmax-1) )...
    +(1-alphan(2:imax,2:jmax-1)).* ( phi(2:imax,3:jmax).* ...
    ( g1_n_n(2:imax,2:jmax-1)-g2_n_n(2:imax,2:jmax-1) ) +...
    phi(2:imax,4:jmax+1).*g2_n_n(2:imax,2:jmax-1) - ...
    phi(2:imax,2:jmax-1).*g1_n_n(2:imax,2:jmax-1) ) );

S_s(3:imax-1,3:jmax)=Fs(3:imax-1,3:jmax).*( alphas(3:imax-1,3:jmax).* ...
    ( phi(3:imax-1,2:jmax-1).* ...
    ( g2_s_p(3:imax-1,3:jmax)-g1_s_p(3:imax-1,3:jmax) ) + ...
    phi(3:imax-1,3:jmax).*g1_s_p(3:imax-1,3:jmax) - ...
    phi(3:imax-1,1:jmax-2).*g2_s_p(3:imax-1,3:jmax) )...
    +(1-alphas(3:imax-1,3:jmax)).* ( phi(3:imax-1,3:jmax).* ...
    ( g2_s_n(3:imax-1,3:jmax)-g1_s_n(3:imax-1,3:jmax) ) +...
    phi(3:imax-1,2:jmax-1).*g1_s_n(3:imax-1,3:jmax) - ...
    phi(3:imax-1,4:jmax+1).*g2_s_n(3:imax-1,3:jmax) ) );


S_w(2,2:jmax) = (Fw(2,2:jmax)+ D2_a(1,2:jmax).*phi_a(1,2:jmax) );

S_ur(2:imax,2:jmax) = ap_p(2:imax,2:jmax) .* (1-alpha_q) .* phi(2:imax,2:jmax);
%% Immersed Boundary Treating

%****************************
% Body Cells

for i=1:size(I_solid,2)
    
    S_e(I_solid(i),J_solid(i)) = 0;
    S_w(I_solid(i),J_solid(i)) = 0;
    S_n(I_solid(i),J_solid(i)) = 0;
    S_s(I_solid(i),J_solid(i)) = 0;
    S0(I_solid(i),J_solid(i)) = phi_inside;
    S_ur(I_solid(i),J_solid(i)) = 0;

        
end

%****************************    
% Ghost Cells
                        
for i=1:size(I_g,2)

    ind=(J_g(i)-2)*(imax-1)+I_g(i)-1;
    
    S_e(I_g(i),J_g(i))=0;
    S_w(I_g(i),J_g(i))=0;
    S_n(I_g(i),J_g(i))=0;
    S_s(I_g(i),J_g(i))=0;
    S0(I_g(i),J_g(i))=RHS_phi_g(1,i);
    S_ur(I_g(i),J_g(i))=0;
end
    
S = S_e + S_w +S_n +S_s + S0 + S_ur;

S = S(2:imax,2:jmax);    

RHSSCALAR_c = reshape(S,[numel(S),1]);




