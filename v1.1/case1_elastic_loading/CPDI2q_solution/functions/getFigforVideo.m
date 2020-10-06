function [figVid] = getFigforVideo(meD,mpD,xB)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    figVid = figure(1);
    set(figVid,'Units','pixels','Position',[100 100 1127 283]);
    clf
    h = scatter(mpD.x(:,1),mpD.x(:,2));
    h.CData = mpD.epII;
    Cdata = h.CData;
    cmap = colormap('parula');
    % make it into a index image.
    cmin = min(Cdata(:));
    cmax = max(Cdata(:));
    m = length(cmap);
    index = fix((Cdata-cmin)/(cmax-cmin)*m)+1; %A
    % Then to RGB
    RGB = ind2rgb(index,cmap);
    
    xs = [mpD.xc mpD.xc(:,1)];
    ys = [mpD.yc mpD.yc(:,1)];
    xn = reshape(meD.x,meD.nNy,meD.nNx);
    yn = reshape(meD.y,meD.nNy,meD.nNx);
    
    plot(xn,yn,'k-',xn',yn','k-');
    hold on
    fill(xs',ys',RGB);
    hold off
    axis equal;
    box on;
    xlim(xB);
    ylim([0 11.8]);
    xlabel('$x$ (m)');ylabel('$y$ (m)');
    set(gca,'FontSize',15,'TickLabelInterpreter','latex');
    drawnow;
end

