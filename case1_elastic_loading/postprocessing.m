clear,close
set(0,'defaulttextinterpreter','latex')                                   ;%
fslab  = 14; fsleg  = 14; fstit  = 14; fstick = 14;

l0 = 10;

for i = 1:12
name1 = '.\CPDI_solution\data_convergence\CPDI_data_';
name2 = num2str(i);
name3 = '.mat';

CPDI      = load([name1 name2 name3]);
CPDI_h(i,:)    = CPDI.meD.h;
H(i).h0   = l0;
H(i).yp   = H(i).h0-CPDI.mpD.x(:,2);
H(i).syyp = CPDI.mpD.s(2,:)./1e3;
H(i).syya = -CPDI.rho0.*CPDI.g.*(l0-CPDI.y0)./1e3;
H(i).V0   = (CPDI_h(i,1)*l0);
CPDI_error(i)= sum(abs(H(i).syyp-H(i).syya').*CPDI.mpD.V',2)./(CPDI.meD.h(1).*l0.*(CPDI.g*CPDI.rho0.*l0))
end
for i = 1:12
name1 = '.\iCPDI_solution\data_convergence\iCPDI_data_';
name2 = num2str(i);
name3 = '.mat';

iCPDI      = load([name1 name2 name3]);
iCPDI_h(i,:)    = iCPDI.meD.h;
H(i).h0   = l0;
H(i).yp   = H(i).h0-iCPDI.mpD.x(:,2);
H(i).syyp = iCPDI.mpD.s(2,:)./1e3;
H(i).syya = -iCPDI.rho0.*iCPDI.g.*(l0-iCPDI.y0)./1e3;
H(i).V0   = (iCPDI_h(i,1)*l0);
iCPDI_error(i)= sum(abs(H(i).syyp-H(i).syya').*iCPDI.mpD.V',2)./(iCPDI.meD.h(1).*l0.*(iCPDI.g*iCPDI.rho0.*l0))
end
for i = 1:12
name1 = '.\GIMP_solution\data_convergence\GIMP_data_';
name2 = num2str(i);
name3 = '.mat';

GIMP      = load([name1 name2 name3]);
GIMP_h(i,:)    = GIMP.meD.h;
H(i).h0   = l0;
H(i).yp   = H(i).h0-GIMP.mpD.x(:,2);
H(i).syyp = GIMP.mpD.s(2,:)./1e3;
H(i).syya = -GIMP.rho0.*GIMP.g.*(l0-GIMP.y0)./1e3;
H(i).V0   = (GIMP_h(i,1)*l0);
GIMP_error(i)= sum(abs(H(i).syyp-H(i).syya').*GIMP.mpD.V',2)./(GIMP.meD.h(1).*l0.*(GIMP.g*GIMP.rho0.*l0))



end

X = 1./CPDI_h(1:7,2)
Y = CPDI_error(1,1:7)

a = 2.316e-5
b = -1.0
x = 0.005:0.5:500
y = a.*x.^b

a = 1.25e-4
b = -2.0
x2 = 0.005:0.5:500
y2 = a.*x.^b


fig1=figure(1)
clf
sim = 3
set(fig1,'Units','pixels','Position',[100 100 541 277]);
hold on
ax1 = plot(H(sim).yp,H(sim).syyp,'o'                        ,'MarkerSize',10,'LineWidth',2)
ax2 = plot(H(sim+2).yp(5:5:end),H(sim+2).syyp(5:5:end),'d'  ,'MarkerSize',10,'LineWidth',2)
ax3 = plot(H(sim+4).yp(1:10:end),H(sim+4).syyp(1:10:end),'s','MarkerSize',10,'LineWidth',2)
ax4 = plot(H(sim+4).yp,H(sim+4).syya,'r--'                   ,'LineWidth',2)
hold off
box on
tit = {['$1/h=',num2str(1./GIMP_h(sim,2)),'$ m$^{-1}$'],['$1/h=',num2str(1./GIMP_h(sim+2,2),'%.1f'),'$ m$^{-1}$'],['$1/h=',num2str(1./GIMP_h(sim+4,2),'%.1f'),'$ m$^{-1}$'],'$\sigma_{yy}(y_0)$'};
h1=legend([ax1 ax2 ax3 ax4],tit);
set(h1,'Interpreter','latex','FontSize',12,'Position',[0.5745 0.6002 0.3171 0.2971],'NumColumns',1);
xlabel('$y_0$ (m)')
ylabel('$\sigma_{yy}$ (kPa)')
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
str = '\textbf{(a)}';
t=text(0.0197,0.1203,str,'Interpreter','latex','Units','normalized','FontSize',18);
print(fig1,'figElasticCompactionStress','-depsc')
print(fig1,'figElasticCompactionStress','-dpng')


fig2=figure(2);
clf
set(fig2,'Units','pixels','Position',[100 100 541 277]);
hold on
ax4=loglog(x,y,'k--','LineWidth',2)
ax5=loglog(x2,y2,'k:','LineWidth',2)
ax1=loglog(1./CPDI_h(:,1),CPDI_error,'bo','MarkerSize',10,'LineWidth',2)
ax2=loglog(1./iCPDI_h(:,1),iCPDI_error,'go','MarkerSize',6,'LineWidth',2)
ax3=loglog(1./GIMP_h(:,1),GIMP_error,'rx','MarkerSize',10,'LineWidth',2)
hold off
grid on
tit = {'CPDI','iCPDI','cpGIMP','$e \propto (1/h)^{-1}$','$e \propto (1/h)^{-2}$'};
h1=legend([ax1 ax2 ax3 ax4 ax5],tit);
set(h1,'Interpreter','latex','FontSize',12,'Position',[0.1642 0.2576 0.3171 0.3687],'NumColumns',1);
set(gca, 'YScale', 'log','XScale','log')
box on;
xlabel('$1/h$ (m$^{-1}$)')
ylabel('$\mathrm{error}$ (-)')
yticks([1e-8 1e-6 1e-4])
xlim([min(1./CPDI_h(:,1))-min(1./CPDI_h(:,1))/5 max(1./CPDI_h(:,1))+max(1./CPDI_h(:,1))/5])
ylim([min(iCPDI_error)-min(iCPDI_error)/5 max(iCPDI_error)+max(iCPDI_error)/5])
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
%str = '\textbf{(b)}';
%t=text(0.0197,0.1203,str,'Interpreter','latex','Units','normalized','FontSize',18);
print(fig2,'figElasticCompactionConvergence','-depsc')
print(fig2,'figElasticCompactionConvergence','-dpng')




