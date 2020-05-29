clear,close
set(0,'defaulttextinterpreter','latex')                                   ;%

CPDI   = load('.\CPDI_solution\data\CPDI_data_1.mat');
GIMP   = load('.\GIMP_solution\data\GIMP_data_1.mat');

dt = CPDI.dt;
nit= CPDI.nit;

fig1=figure(1);
clf
set(fig1,'Units','pixels','Position',[100 100 541 277]);
hold on
ax1=plot(CPDI.dt:CPDI.dt:CPDI.dt*CPDI.nit,-CPDI.duy,'b-' ,'LineWidth',2)
ax2=plot(GIMP.dt:GIMP.dt:GIMP.dt*GIMP.nit,-GIMP.duy,'r-','LineWidth',2)

hold off
ylim([-3.5 0])
grid on
box on
tit = {'CPDI','cpGIMP'};
h1=legend([ax1 ax2],tit);
set(h1,'Interpreter','latex','FontSize',12,'Position',[0.3353 0.5814 0.2979 0.3044],'NumColumns',2);
title(h1,'Vertical deflection');
xlabel('$t$ (s)')
ylabel('$\Delta u$ (m)')
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
print(fig1,'-dpng','figCantileverBeamProblem.png')
