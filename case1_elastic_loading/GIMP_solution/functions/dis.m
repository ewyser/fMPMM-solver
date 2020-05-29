function dis(data,Xp,Yp,time,plotprop)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
set(0,'defaulttextinterpreter','latex');
clf;

check = size(Xp);
da    = size(data);

if(da~=check)
   data = data'; 
end

Data = sortrows([Xp Yp data],3);
hold on;


scatter(Data(:,1),Data(:,2),plotprop.ps,(Data(:,3)),'filled');
plot([plotprop.BoundX(1,1) plotprop.BoundX(1,2)],[0 0],'r--',...
     [plotprop.BoundX(1,1) plotprop.BoundX(1,1)],[0 3*plotprop.BoundY(1,2)],'r--',...
     [plotprop.BoundX(1,2) plotprop.BoundX(1,2)],[0 3*plotprop.BoundY(1,2)],'r--',...
     'MarkerSize',2,'LineWidth',2);
hold off;
axis equal;
xlim(plotprop.LimX);ylim(plotprop.LimY);
xlabel('$x$ (m)');ylabel('$y$ (m)');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
box on;
if(length(plotprop.caxis)==1)
    if(min(data)<0)
        caxis([-plotprop.caxis plotprop.caxis]);
    else
        caxis([0 plotprop.caxis]);
    end
end
if(length(plotprop.caxis)==2)
    caxis([plotprop.caxis(1) plotprop.caxis(2)]);
end
%caxis([0 plotprop.caxis])
ax = gca;
ax.TickDir = 'out';

cb=colormap(flipud(jet(plotprop.cbclass)));
cb=colormap((jet(plotprop.cbclass)));
%cb=colormap(parula);
cb=colorbar('FontSize',14,'TickLabelInterpreter','latex','Fontsize',14,'Location','south');
cb.Position         =plotprop.cbpos;
cb.Label.String     =plotprop.cbchar;
cb.Label.FontSize   =14;
cb.Label.Interpreter='Latex';
cb.Label.Position   =plotprop.cblpos;
title(plotprop.tit);

drawnow;
end

