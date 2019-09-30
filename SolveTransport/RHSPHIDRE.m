function [RHSSCALAR_c] = RHSPHIDRE(StateVar,BC,...
        Flux,DOMAIN,IBM_coeffP,IBM)
% Desciption:
%
%% Define parameters
%% Define parameters
imax = DOMAIN.imax;
jmax = DOMAIN.jmax;
L    = (imax-1)*(jmax-1); % Length of the matrix



A0_p = Flux.Diffu_P.A0_p_p;


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
% RHSSCALAR = zeros(L,1);
S_e=zeros(1,L);
S_w=zeros(1,L);
S_n=zeros(1,L);
S_s=zeros(1,L);
S0=zeros(1,L);


%% Fill in sparse matrix
ind = 0;
for j = 2:jmax      % Loop over all interior pressure points
  for i = 2:imax
    ind = ind+1;        
    if(i == 2 && j == 2)                    % Bottom-left corner
          
      S0(ind)=A0_p(i,j)*phi_old(i,j);

    elseif(i == 2 && j == jmax)             % Top-left corner

        S0(ind)=A0_p(i,j)*phi_old(i,j);


    elseif(i == imax && j == 2)             % Bottom-right corner
            
     S0(ind)=A0_p(i,j)*phi_old(i,j);
     

    elseif(i == imax && j == jmax)          % Top-right corner
                       
        S0(ind)=A0_p(i,j)*phi_old(i,j);

    elseif(i == 2)      % Left edge
                  
        
        S0(ind)=A0_p(i,j)*phi_old(i,j);

    elseif(i == imax)                       % Right edge
              
      S0(ind)=A0_p(i,j)*phi_old(i,j);
 
    elseif(j == 2)                          % Bottom edge
            
      
      S0(ind)=A0_p(i,j)*phi_old(i,j);
    elseif(j == jmax)                       % Top edge
            
      
      S0(ind)=A0_p(i,j)*phi_old(i,j);
    else                % All others
      
            S0(ind)=A0_p(i,j)*phi_old(i,j);
            
            
    end
  end
end

%% Body Cells

for i=1:size(I_solid,2)
    
    ind=(J_solid(i)-2)*(imax-1)+I_solid(i)-1;
    
    S_e(ind)=0;
    S_w(ind)=0;
    S_n(ind)=0;
    S_s(ind)=0;
    S0(ind)=phi_inside;
    
end

                        
for i=1:size(I_g,2)

    ind=(J_g(i)-2)*(imax-1)+I_g(i)-1;
    
    S_e(ind)=0;
    S_w(ind)=0;
    S_n(ind)=0;
    S_s(ind)=0;
    S0(ind)=RHS_phi_g(1,i);
end
    



RHSSCALAR_c=S_s'+S_n'+S_w'+S_e'+S0';


