clear,close
set(0,'defaulttextinterpreter','latex')                                   ;%

CPDI   = load('.\data\LS_CPDI_data_1.mat');
CPDINH = load('.\data\NH_CPDI_data_1.mat');
GIMP   = load('.\data\LS_GIMPM_data_1.mat');
GIMPNH = load('.\data\NH_GIMPM_data_1.mat');

CPDI_Sad    = load('.\Sadeghirad_solution\data_CPDI_vertical_deflection_du.txt')
FEM_Sad    = load('.\Sadeghirad_solution\data_FEM_vertical_deflection_du.txt')

dt = CPDI.DT;
nit= CPDI.nit;

fig1=figure(1);
clf
set(fig1,'Units','pixels','Position',[100 100 541 277]);
hold on
ax6=plot(CPDI_Sad(1:3:end,1),CPDI_Sad(1:3:end,2),'ko','LineWidth',1)
ax5=plot(FEM_Sad(1:1:end,1),FEM_Sad(1:1:end,2),'ks','LineWidth',1)
ax1=plot(CPDI.DT(2:end),-CPDI.duy,'b-' ,'LineWidth',2)
ax2=plot(GIMP.DT(2:end),-GIMP.duy,'r-','LineWidth',2)
ax3=plot(CPDINH.DT(2:end),-CPDINH.duy,'b--','LineWidth',2)
ax4=plot(GIMPNH.DT(2:end),-GIMPNH.duy,'r--','LineWidth',2)

hold off
ylim([-3.5 0])
grid on
box on
tit = {'FEM','CPDI','CPDI','CPDI','cpGIMP','cpGIMP'};
h1=legend([ax5 ax1 ax3 ax6 ax2 ax4],tit);
set(h1,'Interpreter','latex','FontSize',12,'Position',[0.3353 0.5814 0.2979 0.3044],'NumColumns',2);
title(h1,'Vertical deflection');
xlabel('$t$ (s)')
ylabel('$\Delta u$ (m)')
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
print(fig1,'-depsc','figCantileverBeamProblem.eps');
print(fig1,'-dpng' ,'figCantileverBeamProblem.png')






CPDI   = load('.\CPDI_solution\data\CPDI_data_1.mat');
GIMP   = load('.\GIMPM_solution\data\GIMPM_data_1.mat');

ps = 10

fig3=figure(3);
set(fig3,'Units','pixels','Position',[33.6667 108.3333 613.3333 532.6667]);
subplot(221)
xs=[CPDI.mpD.xc CPDI.mpD.xc(:,1)];
ys=[CPDI.mpD.yc CPDI.mpD.yc(:,1)];
plot(xs',ys','k-');axis equal;axis tight;
ylabel('$y$ (m)');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
xlim([0 CPDI.meD.L(1)]);
ylim([0 CPDI.meD.L(2)]);
yticks([0 2.5 5]);
xticks([0 2.5 5]);
title('Finite deformation')
set(gca,'XTickLabel',[])
subplot(222)
du = CPDI.mpD.s(2,:)./1e3;%sqrt(mpD.u(:,1).^2+mpD.u(:,2).^2);
ax2=scatter(CPDI.mpD.x(:,1),CPDI.mpD.x(:,2),ps,du,'filled');
axis equal;
box on;
%xlabel('$x$ (m)');ylabel('$y$ (m)');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
colormap(gca,(jet));
cb1=colorbar('FontSize',14,'TickLabelInterpreter','latex','Fontsize',14,'Location','East');
cb1.Label.FontSize   =14;
cb1.Label.Interpreter='Latex';
xlim([0 CPDI.meD.L(1)]);
ylim([0 CPDI.meD.L(2)]);
yticks([0 2.5 5]);
xticks([0 2.5 5]);
title(['$\sigma_{yy}$ [kPa]']);
set(gca,'YTickLabel',[],'XTickLabel',[])


xc = repmat(GIMP.mpD.x(:,1),1,4)+GIMP.mpD.l(:,1).*[-1 1 1 -1];
yc = repmat(GIMP.mpD.x(:,2),1,4)+GIMP.mpD.l(:,2).*[-1 -1 1 1];
subplot(223)
xs=[xc xc(:,1)];
ys=[yc yc(:,1)];
plot(xs',ys','k-');axis equal;axis tight;
xlabel('$x$ (m)');
ylabel('$y$ (m)');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
xlim([0 GIMP.meD.L(1)]);
ylim([0 GIMP.meD.L(2)]);
yticks([0 2.5 5]);
xticks([0 2.5 5]);
%title('Finite deformation')
subplot(224)
du = GIMP.mpD.s(2,:)./1e3;%sqrt(mpD.u(:,1).^2+mpD.u(:,2).^2);
ax2=scatter(GIMP.mpD.x(:,1),GIMP.mpD.x(:,2),ps,du,'filled');
axis equal;
box on;
xlabel('$x$ (m)');
%ylabel('$y$ (m)');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
colormap(gca,(jet));
cb2=colorbar('FontSize',14,'TickLabelInterpreter','latex','Fontsize',14,'Location','East');
cb2.Label.FontSize   =14;
cb2.Label.Interpreter='Latex';
%%cb2.Label.Position   =plotprop.cblpos;
xlim([0 GIMP.meD.L(1)]);
ylim([0 GIMP.meD.L(2)]);
xticks([0 2.5 5]);
set(gca,'YTickLabel',[])
%title(['$\sigma_{yy}$ [kPa]']);

str = '\textbf{(a)}';
t1=text(-1.4858,1.4732,str,'Interpreter','latex','Units','normalized','FontSize',18);
str = '\textbf{(b)}';
t2=text(0.0013,1.4732,str,'Interpreter','latex','Units','normalized','FontSize',18);
str = '\textbf{(c)}';
t3=text(-1.4858,0.0835,str,'Interpreter','latex','Units','normalized','FontSize',18);
str = '\textbf{(d)}';
t4=text(0.0068,0.0854,str,'Interpreter','latex','Units','normalized','FontSize',18);

str = '$t=1.5$ s';
t5=text(-0.4853,1.1930,str,'Interpreter','latex','Units','normalized','FontSize',18);

set(gca,'FontSize',15,'TickLabelInterpreter','latex');
print(fig3,'-depsc','figCantileverBeamProblemStress.eps');
print(fig3,'-dpng' ,'figCantileverBeamProblemStress.png')





