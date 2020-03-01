function psi = FractureInitialize(DOMAIN, BoundaryCurve)

xt = BoundaryCurve.xt;
yt = BoundaryCurve.yt;
yb = BoundaryCurve.yb;
fbTop = BoundaryCurve.boundfunctionTop;
fbBottom = BoundaryCurve.bFunctionBottom;

endBotIndx = BoundaryCurve.endBotIndx;

b = BoundaryCurve.b;
lambda = BoundaryCurve.lambda;
a = BoundaryCurve.a;
x_0 = BoundaryCurve.x_0;

X = DOMAIN.Xp;
Y = DOMAIN.Yp - DOMAIN.ly/2;
psi = zeros(size(X));
% profile on
% dispIndx = ceil(length(x)/20);
for i =1:length(DOMAIN.xp)
%     if ~mod(i,dispIndx)
%         display(strcat('...', num2str(i/length(x) *100), '.............'))
%     end
    indx = abs(xt-X(i,1)) < lambda/2;
    xt_temp = xt(indx);
    yt_temp = yt(indx);
    yb_temp = yb(indx);
    %%%%% bottom
    dist = sqrt((xt_temp' - X(i,1:endBotIndx)).^2 + ...
        (yb_temp' - Y(i,1:endBotIndx)).^2 );
    y_curve = fbBottom(X(i,1), a, lambda, x_0, b);

    psi(i,1:endBotIndx) = min(dist) .* (Y(i,1:endBotIndx) >= y_curve) - ...
         min(dist) .* (Y(i,1:endBotIndx) <= y_curve) ;

    %%%%% top
    dist = sqrt((xt_temp' - X(i,endBotIndx+1:end)).^2 + ...
        (yt_temp' - Y(i,endBotIndx+1:end)).^2 );
    y_curve = fbTop(X(i,1), a, lambda, x_0, b);

    psi(i,endBotIndx+1:end) = -min(dist) .* (Y(i,endBotIndx+1:end) >= y_curve) + ...
         min(dist) .* (Y(i,endBotIndx+1:end) <= y_curve) ; 

end

end
% profile viewer
%%
% contourf(X,Y,psi, 50, 'LineStyle', 'none')
% colormap jet
% hold on
% plot(xt,yt, 'k')
% hold on 
% plot(xt,yb, 'k')
% axis equal

