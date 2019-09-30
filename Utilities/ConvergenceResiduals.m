function [resid,resi]=ConvergenceResiduals(U,V,...
    RHS_P,ii,time,PCOR,CM_u,CM_v,RHS_U,RHS_V,DOMAIN,disc_scheme,...
    tol,BC,iter_qq_u,iter_qq_v,err_q_u,err_q_v)
%% Retrieve variable
imax = DOMAIN.imax;
jmax = DOMAIN.jmax;

Lu    = (imax-2)*(jmax-1); % Length of the matrix
U_star_v=zeros(Lu,1);



Lv    = (imax-1)*(jmax-2);
V_star_v=zeros(Lv,1);

% if nargin==10
%     
%     resi1=abs(sum(sum(sqrt(U(2:end-1,2:end-1)-U_star_old(2:end-1,2:end-1)).^2)))/((imax-1)*(jmax));
%     resi2=abs(sum(sum(sqrt(V(2:end-1,2:end-1)-V_star_old(2:end-1,2:end-1)).^2)))/((imax)*(jmax-1));
%     resi3=abs(sum(sum(sqrt(RHS_P).^2)))/((imax+1)*(jmax+1));%max(RHS_P);%
%     resi4=abs(sum(sum(sqrt(PCOR).^2)))/((imax+1)*(jmax+1));%max(RHS_P);%
% 
% elseif nargin==8
% 
%     %% Maximum Residual
%     resi1=max(max(abs(U(2:end-1,2:end-1)-U_star_old(2:end-1,2:end-1))));
%     resi2=max(max(abs(V(2:end-1,2:end-1)-V_star_old(2:end-1,2:end-1))));
%     resi3=max(RHS_P);
%     resi4=max(max(abs(PCOR(2:end-1,2:end-1))));
% 
% elseif nargin>10
ind=0;
for j=2:jmax
    for i=2:imax-1
        ind=ind+1;
        U_star_v(ind)=U(i,j);
    end
end

if BC.BC_e_p == 1
    Lu    = (imax-1)*(jmax-1); % Length of the matrix
    U_star_v=zeros(Lu,1);
    
    ind=0;
    for j=2:jmax
        for i=2:imax
            ind=ind+1;
            U_star_v(ind)=U(i,j);
        end
    end

end


ind=0;
for j=2:jmax-1
    for i=2:imax
        ind=ind+1;
        V_star_v(ind)=V(i,j);            
    end
end
    
    %% Normalized Residual
%     % max-norm
    resi1=max(abs(CM_u*U_star_v-RHS_U))./max(abs(diag(CM_u).*U_star_v));
    resi2=max(abs(CM_v*V_star_v-RHS_V))./max(abs(diag(CM_v).*V_star_v));
    resi3=max(RHS_P);
    resi4=max(max(abs(PCOR(2:end-1,2:end-1))));
    
    % second-norm
%     resi1 = norm(CM_u*U_star_v-RHS_U)/length(CM_u*U_star_v-RHS_U);
%     resi2 = norm(CM_v*V_star_v-RHS_V)/length(CM_v*V_star_v-RHS_V);
%     resi3=norm(RHS_P)/length(RHS_P);
%     resi4=max(max(abs(PCOR(2:end-1,2:end-1))));

% end

%% 
%%
resid=[resi1,resi2,resi3 resi4];
resi=max(resid(1:3));
% c='k';
% s=3;
% hh=figure(1);
% hold on
% subplot(3,1,1)
% scatter(ii,resi1,s,c,'fill')    % U- res
% set(gca,'yscale','log')
% hold on
% subplot(3,1,2)
% scatter(ii,resi2,s,c,'fill')    % V- res
% set(gca,'yscale','log')
% hold on
% subplot(3,1,3)
% scatter(ii,resi3,s,c,'fill')    %P- res
% set(gca,'yscale','log')
% hold on
% resi4=max(max(abs(PCOR(2:end-1,2:end-1))));
% M_in=sum(U(1,2:end-1));
% M_out=sum(U(end,2:end-1));
% MM=M_out;
if disc_scheme == 3
    if resi>tol || ii<2
        displ = [ii iter_qq_u iter_qq_v err_q_u err_q_v resi];
        s1= 'At SIMPLE iter : %3.1f , QUICK iter for U and V is %3.1f';
        s2=', %3.1f the error is %8.6f and %8.6f max resi is %8.6f. \n';
        formatSpec = strcat(s1,s2);
        fprintf(formatSpec,displ);
    end



    if resi<tol && ii>1
        s1 = '\n ~~~~~~~~~~~~~~~~~~~~ time = %8.6f';
        s2 = '~~~~~~~~~~~~~~~~~~~~~ \n';
        formatSpec = strcat(s1,s2);
        fprintf(formatSpec,time);
    
        s1 = 'Residual for U, V and mass is %8.6f, %8.6f, %8.6f';
        s2 = 'at iter = %3.1f.\n\n';
        formatSpec = strcat(s1,s2);
        fprintf(formatSpec,[resid(1:3),ii]);
        formatSpec = '\n ~~~~~~~~~~~~~~~~~~~~~~%8.6f~~~~~~~~~~~~~~~~~~ \n';
        fprintf(formatSpec,[]);

    end
    
elseif disc_scheme == 1 || disc_scheme == 2
    if resi>tol || ii<2
        displ = [ii resid(1:3)];
        s1= 'At SIMPLE iter : %3.0f , Residual for U, V and mass is ';
        s2=' %8.8f, %8.8f, %8.8f.\n';
        formatSpec = strcat(s1,s2);
        fprintf(formatSpec,displ);
    end



    if resi<tol && ii>1
        s1 = '\n ~~~~~~~~~~~~~~~~~~~~ time = %8.6f';
        s2 = '~~~~~~~~~~~~~~~~~~~~~ \n';
        formatSpec = strcat(s1,s2);
        fprintf(formatSpec,time);
    
        s1 = 'Residual for U, V and mass is %8.8f, %8.8f, %8.8f';
        s2 = 'in %3.1f iterations.\n\n';
        formatSpec = strcat(s1,s2);
        fprintf(formatSpec,[resid(1:3),ii]);
        formatSpec = '\n ~~~~~~~~~~~~~~~~~~~~~~%8.6f~~~~~~~~~~~~~~~~~~ \n';
        fprintf(formatSpec,[]);

    end
end
% displ = [iter_qq_u iter_qq_v err_q_u err_q_v];
% formatSpec = 'The QUICK iteration for U and V is %8.3f, %8.3f and the error is %8.3f and %8.3f. \n';
% fprintf(formatSpec,displ)

% disp([' # of Iterataion is  = ',num2str(ii)...
%      ,' at time  = ',num2str(time)]); 
%  disp([' The iteration for QUICK U-def correction is = ',...
%      num2str(iter_qq_u),' for V is = ',num2str(iter_qq_v)...
%      ,', the error for U and V is respectively  = ',num2str(err_q_u), ...
%      'and   ',num2str(err_q_v)]);
%  
 
 