clear,close
set(0,'defaulttextinterpreter','latex')                                   ;%
fslab  = 14; fsleg  = 14; fstit  = 14; fstick = 14;

% % cd 'C:\Users\Manu\Desktop\case1_ElasticLoading_new\GIMPM_iterative_solution'
% % run UNIL_2DGIMPM_E__Loading_loop
% % cd 'C:\Users\Manu\Desktop\case1_ElasticLoading_new\GIMPM_solution'
% % run UNIL_2DGIMPM_E__Loading
% % cd 'C:\Users\Manu\Desktop\case1_ElasticLoading_new\CPDI2q_solution'
% % run UNIL_2DCPDI2q_E__Loading


cd 'C:\Users\agnes\Desktop\MPM\work\technical_note\case1_ElasticLoading_new'
l0 = 10;

for i = 1:11
name1 = 'C:\Users\agnes\Desktop\MPM\work\technical_note\case1_ElasticLoading_new\CPDI2q_solution\data\CPDI2q_data_';
name2 = num2str(i);
name3 = '.mat';

CPDI      = load([name1 name2 name3]);
CPDI_h(i,:)    = CPDI.meD.h;
H(i).h0   = l0;
H(i).yp   = H(i).h0-CPDI.mpD.x(:,2);
H(i).yp1  = CPDI.mpD.x(:,2);
H(i).syyp = CPDI.mpD.s(2,:)./1e3;
H(i).syya = -CPDI.rho0.*CPDI.g.*(l0-CPDI.y0)./1e3;
H(i).V0   = (CPDI_h(i,1)*l0);
CPDI_error(i)= sum(abs(H(i).syyp-H(i).syya').*CPDI.mpD.V',2)./(CPDI.meD.h(1).*l0.*(CPDI.g*CPDI.rho0.*l0))

CPDI2q_tsolve(i,1) = CPDI.meD.h(1);
CPDI2q_tsolve(i,2) = CPDI.mpD.n;
name1 = 'C:\Users\agnes\Desktop\MPM\work\technical_note\case1_ElasticLoading_new\CPDI2q_solution\data\CPDI2q_time_vectorized_';
name2 = num2str(i);
name3 = '.mat';
CPDI  = load([name1 name2 name3]);  
CPDI2q_tsolve(i,3) = CPDI.tsolve;
end
x = H(end).yp1;
y = H(end).syyp;
plot(y,x,'o-');



fig1=figure(1)
clf
sim = 3
set(fig1,'Units','pixels','Position',[100 100 541 277]);
hold on
% ax1 = plot(H(sim  ).syyp,H(sim).yp1,'o'                        ,'MarkerSize',10,'LineWidth',2)
% ax2 = plot(H(sim+2).syyp(5:5:end),H(sim+2).yp1(5:5:end),'d'  ,'MarkerSize',10,'LineWidth',2)
% ax3 = plot(H(sim+4).syyp(1:10:end),H(sim+4).yp1(1:10:end),'s','MarkerSize',10,'LineWidth',2)

ax1 = plot(sort(H(sim  ).syyp),sort(H(sim  ).yp1),'-o'                        ,'MarkerSize',4,'LineWidth',2)
ax2 = plot(sort(H(sim+2).syyp(5:5:end)),sort(H(sim+2).yp1(5:5:end)),'-d'  ,'MarkerSize',4,'LineWidth',2)
ax3 = plot(sort(H(sim+4).syyp(1:10:end)),sort(H(sim+4).yp1(1:10:end)),'-s','MarkerSize',4,'LineWidth',2)

ax4 = plot(sort(H(sim+4).syya),sort(H(sim+4).yp1),'r--'                   ,'LineWidth',2)
hold off
box on
tit = {['$1/h=',num2str(1./CPDI_h(sim,2)),'$ m$^{-1}$'],['$1/h=',num2str(1./CPDI_h(sim+2,2),'%.1f'),'$ m$^{-1}$'],['$1/h=',num2str(1./CPDI_h(sim+4,2),'%.1f'),'$ m$^{-1}$'],'analytical solution'};
h1=legend([ax1 ax2 ax3 ax4],tit);
set(h1,'Interpreter','latex','FontSize',12,'Position',[0.1402 0.6074 0.3319 0.2971],'NumColumns',1);
ylabel('$y_p$ (m)')
xlabel('$\sigma_{yy}$ (kPa)')
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
str = '\textbf{(b)}';
t=text(0.0197,0.1203,str,'Interpreter','latex','Units','normalized','FontSize',18);
print(fig1,'figElasticCompactionStress','-depsc')
print(fig1,'figElasticCompactionStress','-dpng')

for i = 1:11
name1 = '.\GIMPM_solution\data\GIMPM_data_';
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

GIMP_tsolve(i,1) = GIMP.meD.h(1);
GIMP_tsolve(i,2) = GIMP.mpD.n;
name1 = 'C:\Users\agnes\Desktop\MPM\work\technical_note\case1_ElasticLoading_new\GIMPM_solution\data\GIMPM_time_vectorized_';
name2 = num2str(i);
name3 = '.mat';
GIMP  = load([name1 name2 name3]);  
GIMP_tsolve(i,3) = GIMP.tsolve;
end
clear GIMP_h
for i = 1:10
name1 = '.\GIMPM_iterative_solution\data\GIMPM_data_';
name2 = num2str(i);
name3 = '.mat';

GIMP_it      = load([name1 name2 name3]);
GIMP_h(i,:)    = GIMP_it.meD.h;
H(i).h0   = l0;
H(i).yp   = H(i).h0-GIMP_it.mpD.x(:,2);
H(i).syyp = GIMP_it.mpD.s(2,:)./1e3;
H(i).syya = -GIMP_it.rho0.*GIMP_it.g.*(l0-GIMP_it.y0)./1e3;
H(i).V0   = (GIMP_h(i,1)*l0);
GIMP_it_error(i)= sum(abs(H(i).syyp-H(i).syya').*GIMP_it.mpD.V',2)./(GIMP_it.meD.h(1).*l0.*(GIMP_it.g*GIMP_it.rho0.*l0))

plot(H(i).yp,H(i).syya,'o',H(i).yp,H(i).syyp,'x')
drawnow
GIMP_it_tsolve(i,1) = GIMP_it.meD.h(1);
GIMP_it_tsolve(i,2) = GIMP_it.mpD.n;
name1 = 'C:\Users\agnes\Desktop\MPM\work\technical_note\case1_ElasticLoading_new\GIMPM_iterative_solution\data\GIMPM_time_iterative_';
name2 = num2str(i);
name3 = '.mat';
GIMP_it  = load([name1 name2 name3]);  
GIMP_it_tsolve(i,3) = GIMP_it.tsolve;
end

for i = 1:11
name1 = '.\iCPDI_solution\data\iCPDI_data_';
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


X = 1./CPDI_h(1:7,2)
Y = CPDI_error(1,1:7)

a = 2.316e-6/1.3
b = -2.0
x = 0.005:0.5:100
y = a.*x.^b
% a = 1.7e-5
% b = -1.0
% x = 0.005:0.5:100
% y = a.*x.^b

fig2=figure(2);
clf
set(fig2,'Units','pixels','Position',[100 100 541 277]);
hold on
ax4=loglog(x,y,'k--','LineWidth',2)
ax1=loglog(1./CPDI_h(:,1),CPDI_error,'b-o','MarkerSize',4,'LineWidth',2)
ax2=loglog(1./CPDI_h(:,1),GIMP_error,'r-x','MarkerSize',4,'LineWidth',2)
ax3=loglog(1./GIMP_h(:,1),GIMP_it_error,'g--d','MarkerSize',4,'LineWidth',2)
ax5=loglog(1./iCPDI_h(:,1),iCPDI_error,':o','MarkerSize',4,'LineWidth',2,'Color',[1 0 1])
hold off

tit = {'CPDI2q (vectorised)','cpGIMP (vectorised)','cpGIMP (iterative)','iCPDI2q'};
h1=legend([ax1 ax2 ax3 ax5],tit);
set(h1,'Interpreter','latex','FontSize',12,'Position',[0.5277 0.6074 0.3664 0.2971],'NumColumns',1);
set(gca, 'YScale', 'log','XScale','log')
box on;
xlabel('$1/h$ (m$^{-1}$)')
ylabel('$\mathrm{error}$ (-)')
yticks([1e-6 1e-4])
xlim([min(1./CPDI_h(:,1))-min(1./CPDI_h(:,1))/5 max(1./CPDI_h(:,1))+max(1./CPDI_h(:,1))/5])
ylim([min(iCPDI_error)-min(iCPDI_error)/5 max(CPDI_error)+max(CPDI_error)/5])
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
str = '\textbf{(a)}';
t=text(0.0197,0.1203,str,'Interpreter','latex','Units','normalized','FontSize',18);
%str = '\textbf{(b)}';
%t=text(0.0197,0.1203,str,'Interpreter','latex','Units','normalized','FontSize',18);
print(fig2,'figElasticCompactionConvergence','-depsc')
print(fig2,'figElasticCompactionConvergence','-dpng')

fig3=figure(3);
clf
set(fig3,'Units','pixels','Position',[100 100 541 277]);
hold on
ax1=loglog(GIMP_it_tsolve(:,2),GIMP_it_tsolve(:,3),'o-','MarkerSize',4,'LineWidth',2)
ax2=loglog(GIMP_tsolve(:,2),GIMP_tsolve(:,3),'s-','MarkerSize',4,'LineWidth',2)
ax3=loglog(CPDI2q_tsolve(:,2),CPDI2q_tsolve(:,3),'d--','MarkerSize',4,'LineWidth',2)
% ax4=loglog([min(GIMP_tsolve(:,2));max(GIMP_tsolve(:,2))],[1;1],'r:','LineWidth',1)
% ax4=loglog([min(GIMP_tsolve(:,2));max(GIMP_tsolve(:,2))],[60;60],'r:','LineWidth',1)
% ax4=loglog([min(GIMP_tsolve(:,2));max(GIMP_tsolve(:,2))],[3600;3600],'r:','LineWidth',1)

axis tight
hold off
tit = {'cpGIMP (iterative)','cpGIMP (vectorised)','CPDI2q (vectorised)'};
h1=legend([ax1 ax2 ax3],tit);
set(h1,'Interpreter','latex','FontSize',12,'Position',[0.1507 0.6146 0.3664 0.2971],'NumColumns',1);
set(gca, 'YScale', 'log','XScale','log')
box on;
yticks([1 60 3600])
% ylim([0 2*3600])
xlabel('$n_p$ (-)')
ylabel('Wall-clock time (s)')
set(gca,'FontSize',15,'TickLabelInterpreter','latex');

print(fig3,'figElasticCompactionWalltime','-depsc')
print(fig3,'figElasticCompactionWalltime','-dpng')