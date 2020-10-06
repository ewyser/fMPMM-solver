% fMPMM-solver: Fast and efficient MATLAB-based MPM solver
%
% Copyright (C) 2020  Emmanuel Wyser, Yury Alkhimenkov, Michel Jaboyedoff, Yury Podladchikov
% -------------------------------------------------------------------------%
% version    : v1.1
% date       : october, 2020
% description: explicit mpm solver based on an updated Lagrangian
% frame in plane strain condition for elasto-plastic problem discretized on
% a 4-noded quadrilateral mesh
%% FANCY CALLS
clear                                                                     ;%
addpath('functions')                                                      ;%   
version    = 'GIMPM_EP_'                                                    ;%
plasticity = false                                                         ;%         
typeD='double'                                                            ;%
disp('EXACT PLASTIC CORRECTION')                                          ;%
disp('------------------------')                                          ;%

set(0,'defaulttextinterpreter','latex')                                   ;%
fslab  = 14; fsleg  = 14; fstit  = 14; fstick = 14;

numel = repmat(20,1,1)                                                    ;%
for sim=1:length(numel)
    disp('------------------------')                                      ;%
    disp(['Run ',num2str(sim),': nel = ',num2str(numel(sim)),''])         ;%
    disp('------------------------')                                      ;%
    
    %% NON-DIMENSIONAL CONSTANT
    ar      = 0.8                                                         ;% aspect ratio thickness/width
    nu      = 0.3                                                         ;% poisson ratio
    ni      = 3                                                           ;% number of mp in h(1) direction
    nstr    = 3                                                           ;% number of stress components
    %---------------------------------------------------------------------%
    
    %% INDEPENDANT PHYSICAL CONSTANT
    g       = 10.0                                                        ;% gravitationnal acceleration [m/s^2]
    E       = 1.0e6                                                       ;% young's modulus             [Pa]
    Gc      = E/(2*(1+nu))                                                ;% shear modulus               [Pa]
    Kc      = E/(3*(1-2*nu))                                              ;% bulk modulus                [Pa]
    lame1   = 0.5*E/(1+nu);
    lame2   = nu*E/((1+nu)*(1-2*nu));
    rho0    = 1050                                                        ;% density                     [kg/m^3]
    n0      = 0.25                                                        ;% initial porosity            [-]
    yd      = sqrt((Kc+4/3*Gc)/rho0)                                      ;% elastic wave velocity       [m/s]
    t       = 1.5                                                         ;% simulation time             [s]
    te      = 0.0                                                         ;% elastic loading             [s]
    %---------------------------------------------------------------------%
    
    %% MESH & MP INITIALIZATION
    [meD,bc] = meSetup(numel(sim),typeD)                                  ;% - see function   
    ly      = meD.L(2)                                                    ;% layer thickness [m]
    [mpD]   = mpSetup(meD,ni,ly,0.0,0.0,0.0,0.0,n0,rho0,nstr,typeD)   ;% - see function
    mpq     = find(mpD.x(:,1)==max(mpD.x(:,1)) &  mpD.x(:,2)==min(mpD.x(:,2)));
    y0      = mpD.x(mpq,2);
    % ISOTROPIC ELASTIC MATRIX
    Del     = [ Kc+4/3*Gc,Kc-2/3*Gc,0.0;...
                Kc-2/3*Gc,Kc+4/3*Gc,0.0;...
                0.0      ,0.0      ,Gc]                                   ;%
    %---------------------------------------------------------------------%
    
    %% DISPLAY PARAMETERS AND RUNTIME INITIALIZATION
    fps  = 25                                                             ;% image per second
    %COURANT-FRIEDRICH-LEVY CONDITION
    C    = 0.1                                                            ;%                 
    dt   = C*meD.h(1)/yd                                                  ;% unconditionally stable timestep
    nit  = ceil(t/dt)                                                     ;% maximum number of interation
    nf   = max(2,ceil(round(1/dt)/fps))                                   ;% number of frame interval
    % GRAVITY INCREMENT
    dg   = linspace(0,g,round(((te)/1.5)/dt))                             ;% gravity increase
    % RUNTIME PARAMETERS
    nc   = 0                                                              ;% initialize iteration counter                                                         
    it   = 1                                                              ;% initialize iteration
    tw   = 0.0                                                            ;% initialize time while statement
    cycle_time = zeros(nit,7)                                             ;%
    % PLOT SETTING
    [pp] = plotSetup(meD.xB(1),meD.xB(2),meD.xB(3),ly,meD.y,bc.y)                     ;% - see function
    %---------------------------------------------------------------------%
    DT = 0
    %% MPM DM ALGORITHM EXPLICIT SOLVER FOR INFINITESIMAL STRAIN
    tsolve=tic                                                            ;%
    fprintf('MPM SOLVER ON: %.0f elements \n'      ,meD.nEx*meD.nEy)      ;%
    fprintf('               %.0f nodes \n'         ,meD.nN)               ;%
    fprintf('               %.0f material point \n',mpD.n)                ;%                                                           ;
    while((tw<t)||(sum(isnan(mpD.v(:,1))+isnan(mpD.v(:,2)))>0.0))% BEGIN WHILE LOOP
        time_it= tic                                                      ;%
        dpi    = time_it                                                  ;% CURRENT ITERATION TIMER BEGIN
        duy(it)= y0-mpD.x(mpq,2)                                          ;%
        %% CFL CONDITION & LINEAR INCREASE OF GRAVITY FOR EQUILIBRIUM
        c  = [max(yd+abs(mpD.v(:,1))),max(yd+abs(mpD.v(:,2)))]            ;%
        dt = C*min(meD.h./c)                                              ;%
        DT(it+1) = DT(it)+dt;
        %------------------------------------------------------------------%

        %% TRACK MATERIAL POINTS (p) IN ELEMENTS (E)
        p2.e = (floor((mpD.x(:,2)-min(meD.y))./meD.h(2))+1)+...
               (meD.nEy).*floor((mpD.x(:,1)-min(meD.x))./meD.h(1))        ;% index mp to element
        p2.N = meD.e2N(p2.e,:)                                            ;% index mp (to element) to node 
        l2g  = [meD.DoF*p2.N-1;meD.DoF*p2.N]                              ;% local to global node index list [x_I,y_I]
        neon = size(unique(p2.e),1)                                       ;% number of active element
        %------------------------------------------------------------------%
        
        %% BASIS FUNCTIONS
        tic;
        [mpD] = SdS(meD,mpD,p2.N)                                         ;% - see function 
        cycle_time(it,1)=toc;
        %------------------------------------------------------------------%
        %% PROJECTION FROM MATERIAL POINTS (p) TO NODES (N)
        tic
        [meD] = p2Nsolve(meD,mpD,g,dt,l2g,p2.N,bc)                        ;% - see function
        cycle_time(it,2)=toc; 
        %------------------------------------------------------------------%

        %% INTERPOLATION FROM NODAL SOLUTIONS (N) TO MATERIAL POINTS (p)
        tic
        [meD,mpD] = mapN2p(meD,mpD,dt,l2g,p2.N,bc)                        ;% - see function
        cycle_time(it,3)=toc;
        %------------------------------------------------------------------%
        
        %% UPDATE INCREMENTAL DEFORMATION & STRAIN 
        tic
        [mpD] = DefUpdate(meD,mpD,l2g)                                    ;% - see function
        cycle_time(it,4)=toc;
        %------------------------------------------------------------------%   
        
        %% ELASTO-PLASTIC RELATION: ELASTIC PREDICTOR - PLASTIC CORRECTOR
        tic
        % ELASTIC PREDICTOR: TRIAL ELASTIC STEP
        [mpD] = elastic(mpD,Del,lame1,lame2)                              ;% - see function
        cycle_time(it,5)=toc                                              ;%
        dpi=toc(dpi)                                                      ;% CURRENT ITERATION TIMER END
        cycle_time(it,6)=dpi                                              ;%
        %------------------------------------------------------------------%
        
        %% TERMINAL DISPLAY
        if(mod(it,nf)==1)
            rt    = ((nit-it)*toc(time_it))                               ;%
            dpi   = mean(1./cycle_time(1:it,6))                           ;%
            clocktimer(rt,'Remaining estimated time:',dpi)                ;%
        end
        if((it*dt<te)&&(it*dt>te-dt))
            figure(2)
            pp.cbchar='$p$ [kPa]';
            pp.cbpos =[0.42 0.5 0.2 0.05];
            pos            = pp.cbpos;
            pp.cblpos=[pos(1) pos(2)+2];
            pp.caxis =(rho0*g*max(mpD.x(:,2)))/1e3;
            pp.ps    =3.0;
            pp.tit   =['time: ',num2str(it*dt-te,'%.2f'),' (s), $g=',num2str(g,'%.2f'),'$ (m/s$^2$)'];
            pp.cbclass = 20;
            dd = -((mpD.s(1,:)+mpD.s(2,:)+mpD.s(3,:))./3)./1e3;
            dis(dd,mpD.x(:,1),mpD.x(:,2),it*dt,pp);
            print(figure(2),['./data/' version 'pressure_equilibrium' '_',num2str(sim),'' ],'-dpng');
        end
        %------------------------------------------------------------------%        
        
        %% ITERATION INCREMENT
        tw=tw+dt                                                          ;%
        it=it+1                                                           ;%
        %------------------------------------------------------------------%
        
    end% END WHILE LOOP
    tsolve = toc(tsolve)                                                  ;%
    clocktimer(tsolve,'Runtime MPM solver:',mean(1./cycle_time(:,6)))     ;%
    

    %% SAVE DATA
    name=['.\data\GIMPM_time_vectorized_' num2str(sim) '.mat'];
    save(name,'tsolve','cycle_time');
    
    name=['.\data\GIMPM_data_' num2str(sim) '.mat'];
    save(name,'mpD','meD','ni','rho0','E','nu','DT','nit','duy');
 
end

mpD.xc = repmat(mpD.x(:,1),1,4)+mpD.l(:,1).*[-1 1 1 -1];
mpD.yc = repmat(mpD.x(:,2),1,4)+mpD.l(:,2).*[-1 -1 1 1];

fig3=figure(3);
set(fig3,'Units','pixels','Position',[139 293 588.6667 264]);
subplot(121)
xs=[mpD.xc mpD.xc(:,1)];
ys=[mpD.yc mpD.yc(:,1)];
plot(xs',ys','k-');axis equal;axis tight;
ylabel('$y$ (m)');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
xlim([0 meD.L(1)]);
ylim([0 meD.L(2)]);
title('Finite deformation')
subplot(122)
du = sqrt(mpD.u(:,1).^2+mpD.u(:,2).^2);
ax2=scatter(mpD.x(:,1),mpD.x(:,2),pp.ps,du,'filled');
axis equal;
box on;
xlabel('$x$ (m)');ylabel('$y$ (m)');
set(gca,'FontSize',15,'TickLabelInterpreter','latex');
colormap(gca,(jet));
xlim([0 meD.L(1)]);
ylim([0 meD.L(2)]);
title(['$\Delta u$, max($\Delta u$) = ',num2str(max(du),'%.2f'),' [m]']);
print(fig3,['./data/' version 'summary' '_',num2str(sim),'' ],'-dpng');


%% REFERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
% Dunatunga, S., Kamrin, K. (2017) Continuum modeling of projectile impact
% and penetration in dry granular media. J. Mech. Phys. Solids. 100. 45-60
%-------------------------------------------------------------------------%
% Dunatunga, S., Kamrin, K. (2015) Continuum modelling and simulation of
% granular flows through their many phases. J. Fluid Mech. 779. 483-513.
%-------------------------------------------------------------------------%
% Stomakhin, A., Schroeder, C., Chai, L., Teran, C., Selle, A. (2013) A
% material point method for snow simulation. ACM Trans. Graph. 32.
%-------------------------------------------------------------------------%
% Gaume, J., Gast, T., Teran, J., van Herwijnen, A., Jiang, C. (2018)
% Dynamic anticrack propagation in snow. Nat. Com. 9
%-------------------------------------------------------------------------%
% Wang, B., Vardon, P.J., Hicks, M.A. (2016) Investigation of retrogressive
% and progressive slope failure mechanisms using the material point method.
% Comp. and Geot. 78. 88-98
%-------------------------------------------------------------------------%
% Liang, W., Zhao, J. (2018) Multiscale modeling of large deformation in
% geomechanics. Int. J. Anal. Methods Geomech. 43. 1080-1114
%-------------------------------------------------------------------------%
% Baumgarten, A.S., Kamrin, K. (2018) A general fluid-sediment mixture
% model and constitutive theory validated in many flow regimes. J. Fluid
% Mech. 861. 7211-764.
%-------------------------------------------------------------------------%
% Huang, P., Zhang, X., Ma, S., Huang, X. (2011) Contact algorithms for the
% material point method in impact and penetration simulation. Int. J. Numer.
% Meth. Engng. 85. 498-517.
%-------------------------------------------------------------------------%
% Homel, M., Herbold, E.B. (2017) Field-gradient partitioning for fracture
% and frictional contact in the material point method. Int. J. Numer. Meth.
% Engng. 109. 1013-1044
%-------------------------------------------------------------------------%
% Ma. S., Zhank, X., Qui, X.M. (2009) Comparison study of MPM and SPH in
% modeling hypervelocity impact problems. Int. J. Imp. Engng. 36. 272-282.
%-------------------------------------------------------------------------%
% Bardenhagen, S.G., Brackbill, J.U., Sulsky, D. (2000) The material-point
% method for granular materials. Comput. Methods Appl. Mech. Engrg. 187.
% 529-541.
%-------------------------------------------------------------------------%
% York, A.R., Sulsky, D., Schreyer, H.L. (1999) The material point method
% for simulation of thin membranes. Int. J. Numer. Meth. Engng. 44. 1429-1456
%-------------------------------------------------------------------------%
% Nairn, J.A., Bardenhagen, S.G., Smith, G.D. (2018) Generalized contact
% and improved frictional heating in the material point method. Comp. Part.
% Mech. 3. 285-296
%-------------------------------------------------------------------------%
% Nairn, J.A. (2013) Modeling imperfect interfaces in the material point
% method using multimaterial methods. Comp. Mod. Engrg. Sci. 92. 271-299.
%-------------------------------------------------------------------------%
% Hammerquist, C.C., Nairn, J.A. (2018) Modeling nanoindentation using the
% material point method. J. Mater. Res. 33. 1369-1381
%-------------------------------------------------------------------------%
% Hamad, F., Stolle, D., Moormann, C. (2016) Material point modelling of
% releasing geocontainers from a barge. Geotext. Geomembr. 44. 308-318.
%-------------------------------------------------------------------------%
% Bhandari, T., Hamad, F., Moormann, C., Sharma, K.G., Westrich, B. (2016)
% Comp. Geotech. 75. 126-134.
%-------------------------------------------------------------------------%
% Wang, B., Vardon, P.J., Hicks, M.A. (2016) Investigation of retrogressive
% and progressive slope failure mechanisms using the material point method.
% Comp. Geotech. 78. 88-98.
%-------------------------------------------------------------------------%
% Keller, T., May, D.A., Kaus, B.J.P. (2013) Numerical modelling of magma
% dynamics oupled to tectonic deformation of lithosphere and crust.
% Geophys. J. Int. 195. 1406-1442.
%
%
%
%
%
%
%
%
% CONTACT ALGORITHM: p.127 in Nguyen, An introduction to Material Point
% Method
% https://math.unm.edu/~sulsky/papers/CMAME.pdf Bardenhagen, 2000


%N-----N%(i  ,j)%(i  ,j+1)%1-----3
%|-p,E-|%%%%%%%%p%%%%%%%%%%|-p,E-|
%N-----N%(i+1,j)%(i+1,j+1)%2-----4


% % % %
% % % %     if(iterative==1)
% % % %         fid = fopen(['./data/' '',num2str(sim),'' '_fem_iterative_strain.txt'],'wt');
% % % %     elseif(exact==1)
% % % %         fid = fopen(['./data/' '',num2str(sim),'' '_fem_exact_strain.txt'],'wt');
% % % %     end
% % % %
% % % %     for ii = 1:length(d1)
% % % %         fprintf(fid,'%.2f \t %.2f \t %.2f \r \n',xp(ii),yp(ii),d1(ii));
% % % %     end
% % % %     fclose(fid);
% % % %
% % % %     if(iterative==1)
% % % %         output = fopen(['./data/' '',num2str(sim),'' '_fem_iterative_runtime.txt'],'w');
% % % %     elseif(exact==1)
% % % %         output = fopen(['./data/' '',num2str(sim),'' '_fem_exact_runtime.txt'],'w');
% % % %     end
% % % %     fprintf(output,'numelx [-]: %.2f  \r\n',nelx);
% % % %     fprintf(output,'numely [-]: %.2f  \r\n',nely);
% % % %     fprintf(output,'active elements [-]: %.2f  \r\n',num_e_on);
% % % %     fprintf(output,'runtime MPM solver [s]: %.2f  \r\n',time_solver);
% % % %     fclose(output);
% % % %

