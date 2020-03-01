clear all

xt = 0:0.01:1000;

b = 50;
lambda = 40;
a = 0.2 * lambda;
x_0 = 0.5*lambda;
y = a*sin(2*pi/lambda * (xt-x_0));

yt = b/2 + y;
yb = -b/2 - y;


x = 0:0.4:1000;
y = -70:0.4:70;

[X,Y] = meshgrid(x,y);
X = X';
Y = Y';
psi = zeros(size(X));
% profile on
%% Bottom
endBotIndx = ceil(length(y)/2);
dispIndx = ceil(length(x)/20);
for i =1:length(x)
    if ~mod(i,dispIndx)
        display(strcat('...', num2str(i/length(x) *100), '.............'))
    end
    [~,J] = find(xt >= X(i,1));
    indx = abs(xt-X(i,1)) < lambda/2;
    xt_temp = xt(indx);
    yt_temp = yt(indx);
    yb_temp = yb(indx);
    %%%%% bottom
    dist = sqrt((xt_temp' - X(i,1:endBotIndx)).^2 + ...
        (yb_temp' - Y(i,1:endBotIndx)).^2 );
    y_curve = -b/2 - a*sin(2*pi/lambda * (X(i,1)-x_0));
    
    psi(i,1:endBotIndx) = min(dist) .* (Y(i,1:endBotIndx) >= y_curve) - ...
         min(dist) .* (Y(i,1:endBotIndx) <= y_curve) ;
    
    %%%%% top
    dist = sqrt((xt_temp' - X(i,endBotIndx+1:end)).^2 + ...
        (yt_temp' - Y(i,endBotIndx+1:end)).^2 );
    y_curve = b/2 + a*sin(2*pi/lambda * (X(i,1)-x_0));
    
    psi(i,endBotIndx+1:end) = -min(dist) .* (Y(i,endBotIndx+1:end) >= y_curve) + ...
         min(dist) .* (Y(i,endBotIndx+1:end) <= y_curve) ; 
     
end
% profile viewer
%%
contourf(DOMAIN.Xp,DOMAIN.Yp,psi, 50, 'LineStyle', 'none')
colormap jet
hold on
plot(xt,yt, 'k')
hold on 
plot(xt,yb, 'k')
% axis equal

