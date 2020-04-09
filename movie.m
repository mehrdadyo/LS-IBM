figure(1)
set(gca, 'nextplot', 'replacechildren');  
caxis manual;          % allow subsequent plots to use the same color limits
caxis([0 1]);         % set the color axis scaling to your min and max color limits

vidfile = VideoWriter('testmovie.mp4','MPEG-4');
open(vidfile);
file_name = 'D:\GitHub\LS-IBM\Output\dataRDE';
for i = 1:314
    file_no = num2str(i*100);
    file_name1 = strcat(file_name, file_no, 'dt.mat');
    load(file_name1);
    z = (psi>0).*phi;
    x = DOMAIN.Xp;
    y = DOMAIN.Yp;
    contourf(x,y,z,20, 'LineStyle', 'none'); colormap jet;colorbar
    hold on
    plot(LSCase.BoundaryCurve.xt, LSCase.BoundaryCurve.yt+0.7, 'k')
    hold on
    plot(LSCase.BoundaryCurve.xt, LSCase.BoundaryCurve.yb+0.7, 'k')
    axis equal
    xlim([DOMAIN.xp(11), DOMAIN.xp(end-410)])   

    drawnow
    F(i) = getframe(gcf); 
    writeVideo(vidfile,F(i));
end
close(vidfile)
