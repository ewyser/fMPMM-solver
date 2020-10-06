% fMPMM-solver: Fast and efficient MATLAB-based MPM solver
%
% Copyright (C) 2020  Emmanuel Wyser, Yury Alkhimenkov, Michel Jaboyedoff, Yury Podladchikov
% -------------------------------------------------------------------------%
% version    : v1.1
% date       : october, 2020
% description: explicit mpm (CPDI2q) solver based on an updated Lagrangian
% frame in plane strain condition for elasto-plastic problem discretized on
% a 4-noded quadrilateral mesh
%% FANCY CALLS
clear                                                                     ;%
addpath('functions')                                                      ;%   
version    = 'GIMPM_EP_'                                                    ;%
plasticity = true                                                         ;%         
typeD='double'                                                            ;%
disp('EXACT PLASTIC CORRECTION')                                          ;%
disp('------------------------')                                          ;%

set(0,'defaulttextinterpreter','latex')                                   ;%
fslab  = 14; fsleg  = 14; fstit  = 14; fstick = 14;

numel = [1 2 5 10 20 40 80 160 320 640 1280]                                       ;%
run   = zeros(length(numel),6)                                            ;%

for sim=1:length(numel)
    disp('------------------------')                                      ;%
    disp(['Run ',num2str(sim),': nel = ',num2str(numel(sim)),''])         ;%
    disp('------------------------')                                      ;%
    
    %% NON-DIMENSIONAL CONSTANT
    ar      = 0.8                                                         ;% aspect ratio thickness/width
    nu      = 0.0                                                         ;% poisson ratio
    ni      = 2                                                           ;% number of mp in h(1) direction
    nstr    = 3                                                           ;% number of stress components
    %---------------------------------------------------------------------%
    
    %% INDEPENDANT PHYSICAL CONSTANT
    g       = 9.81                                                        ;% gravitationnal acceleration [m/s^2]
    E       = 10.0e3                                                      ;% young's modulus             [Pa]
    Gc      = E/(2*(1+nu))                                                ;% shear modulus               [Pa]
    Kc      = E/(3*(1-2*nu))                                              ;% bulk modulus                [Pa]
    rho0    = 80                                                          ;% density                     [kg/m^3]
    n0      = 0.25                                                        ;% initial porosity            [-]
    yd      = sqrt(E/rho0)                                                ;% elastic wave velocity       [m/s]
    coh0    = 20.0e3                                                      ;% cohesion                    [Pa]
    phi0    = 20.0*pi/180                                                 ;% friction angle              [Rad]
    H       = -60e3                                                       ;% softening modulus           [Pa]
    cohr    =  4.0e3                                                      ;% residual cohesion           [Pa]
    phir    = 7.5*pi/180                                                  ;% residual friction angle     [Rad]
    l0      = 10                                                          ;% initial column height
    t       = ceil((1/yd)*(2*l0)*40)                                      ;% simulation time             [s]
    te      = ceil((1/yd)*(2*l0)*40)                                      ;% elastic loading             [s]           [s]
    %---------------------------------------------------------------------%
    
    %% MESH & MP INITIALIZATION
    [meD,bc] = meSetup(numel(sim),typeD,l0)                               ;% - see function   
    yinf     = 0.0;
    [mpD,p2] = mpSetup(meD,ni,l0,yinf,coh0,cohr,phi0,phir,n0,rho0,nstr,typeD);% - see function
    y0       = mpD.x(:,2);
    % ISOTROPIC ELASTIC MATRIX
    Del     = [ Kc+4/3*Gc,Kc-2/3*Gc,0.0;...
                Kc-2/3*Gc,Kc+4/3*Gc,0.0;...
                0.0      ,0.0      ,Gc]                                   ;%
    Hp      = H*meD.h(1)                                                  ;%
    %---------------------------------------------------------------------%
    
    %% DISPLAY PARAMETERS AND RUNTIME INITIALIZATION
    fps  = 25                                                             ;% image per second
    %COURANT-FRIEDRICH-LEVY CONDITION
    C    = 0.5                                                            ;%                 
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
    [pp] = plotSetup(meD.xB(1),meD.xB(2),yinf,l0,meD.y,bc.y)                     ;% - see function
    %---------------------------------------------------------------------%

    %% MPM DM ALGORITHM EXPLICIT SOLVER FOR INFINITESIMAL STRAIN
    tsolve=tic                                                            ;%
    fprintf('MPM SOLVER ON: %.0f elements \n'      ,meD.nEx*meD.nEy)      ;%
    fprintf('               %.0f nodes \n'         ,meD.nN)               ;%
    fprintf('               %.0f material point \n',mpD.n)                ;%                                                           ;
    while(tw<t)% BEGIN WHILE LOOP
        time_it= tic                                                      ;%
        dpi    = time_it                                                  ;% CURRENT ITERATION TIMER BEGIN
        
        %% CFL CONDITION & LINEAR INCREASE OF GRAVITY FOR EQUILIBRIUM
        c  = [max(yd+abs(mpD.v(:,1))),max(yd+abs(mpD.v(:,2)))]            ;%
        dt = C*min(meD.h./c)                                              ;%
        if(tw<=te)
            g = 9.81/te*tw                                                ;% incremented gravity load
        else
            g = 9.81                                                      ;%
        end
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
        [mpD] = constitutive(mpD,Del,Hp,cohr,plasticity,dt,it,te)         ;%
        cycle_time(it,5)=toc                                              ;%
        dpi=toc(dpi)                                                      ;% CURRENT ITERATION TIMER END
        cycle_time(it,6)=dpi                                              ;%
        cycle_time(it,7)=length(unique(p2.N))                             ;%
        %------------------------------------------------------------------%
        
        %% TERMINAL DISPLAY
        if(mod(it,nf)==1)
            rt    = ((nit-it)*toc(time_it))                               ;%
            dpi   = mean(1./cycle_time(1:it,6))                           ;%
            clocktimer(rt,'Remaining estimated time:',dpi)                ;%
        end

        %------------------------------------------------------------------%        
        
        %% CONDITION: NO ACTIVE ELEMENT
        if(isempty(neon)==1 || sum(isnan(mpD.v(:,1))+isnan(mpD.v(:,2)))>0.0)
            break                                                         ;%
        end
        %------------------------------------------------------------------%
        
        %% ITERATION INCREMENT
        tw=tw+dt                                                          ;%
        it=it+1                                                           ;%
        %------------------------------------------------------------------%
        
    end% END WHILE LOOP
    tsolve = toc(tsolve)                                                  ;%
    clocktimer(tsolve,'Runtime MPM solver:',mean(1./cycle_time(:,6)))     ;%
    
    
    
    yp   = y0;
    syyp = mpD.s(2,:)./1e3;
    syya = -rho0.*g.*(l0-yp)./1e3;
    figure(1)
    plot(yp,syyp,'s',yp,syya,'r-')

    name=['.\data\GIMPM_time_vectorized_' num2str(sim) '.mat'];
    save(name,'cycle_time','tsolve');
    name=['.\data\GIMPM_data_' num2str(sim) '.mat'];
    save(name,'mpD','meD','y0','ni','rho0','E','nu','g');
    
    mpD.xc = repmat(mpD.x(:,1),1,4)+mpD.l(:,1).*[-1 1 1 -1];
    mpD.yc = repmat(mpD.x(:,2),1,4)+mpD.l(:,2).*[-1 -1 1 1];
    fig3=figure(3);
    set(fig3,'Units','pixels','Position',[139 98.3333 504 458.6667]);
    subplot(311)
    xs=[mpD.xc mpD.xc(:,1)];
    ys=[mpD.yc mpD.yc(:,1)];
    plot(xs',ys','k-');axis equal;axis tight;
    ylabel('$y$ (m)');
    set(gca,'FontSize',15,'TickLabelInterpreter','latex');
    xlim([0 meD.L(1)]);
    ylim([0 12]);
    title('Finite deformation')
    subplot(312)
    ax1=scatter(mpD.x(:,1),mpD.x(:,2),pp.ps,mpD.epII,'filled');
    axis equal;
    box on;
    ylabel('$y$ (m)');
    set(gca,'FontSize',15,'TickLabelInterpreter','latex');
    colormap(gca,jet);
    xlim([0 meD.L(1)]);
    ylim([0 12]);
    title('$\epsilon_{II}$')
    title(['$\epsilon_{II}$, max($\epsilon_{II}$) = ',num2str(max(mpD.epII),'%.2f'),' [-]']);
    subplot(313)
    du = sqrt(mpD.u(:,1).^2+mpD.u(:,2).^2);
    ax2=scatter(mpD.x(:,1),mpD.x(:,2),pp.ps,du,'filled');
    axis equal;
    box on;
    xlabel('$x$ (m)');ylabel('$y$ (m)');
    set(gca,'FontSize',15,'TickLabelInterpreter','latex');
    colormap(gca,(jet));
    xlim([0 meD.L(1)]);
    ylim([0 12]);
    title(['$\Delta u$, max($\Delta u$) = ',num2str(max(du),'%.2f'),' [m]']);
    print(fig3,['./data/' version 'summary' '_',num2str(sim),'' ],'-dpng');
end



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

