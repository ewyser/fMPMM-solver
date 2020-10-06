% fMPMM-solver: Fast and efficient MATLAB-based MPM solver
%
% Copyright (C) 2020  Emmanuel Wyser, Yury Alkhimenkov, Michel Jaboyedoff, Yury Podladchikov
% -------------------------------------------------------------------------%
% version    : v1.1
% date       : october, 2020
% description: implicit mpm (CPDI2q) solver based on an updated Lagrangian
% frame in plane strain condition for elasto-plastic problem discretized on
% a 4-noded quadrilateral mesh
%% FANCY CALLS
clear                                                                     ;%
addpath('functions')                                                      ;%
version    = 'iCPDI2q_EP_'                                                ;%
datapath   = 'C:\Users\agnes\Desktop\MPM\work\inprogress\implicit\data\'   ;%
plasticity = false                                                        ;%
typeD='double'                                                            ;%
disp('EXACT PLASTIC CORRECTION')                                          ;%
disp('------------------------')                                          ;%

set(0,'defaulttextinterpreter','latex')                                   ;%
fslab  = 14; fsleg  = 14; fstit  = 14; fstick = 14;

numel = repmat(1000,1,1)                                                    ;%
numel = [1 2 5 10 20 40 80 160 320];
numel = [1 2 5 10 20 40 80 160 320 640 1280]
run   = zeros(length(numel),6)                                            ;%

for sim=1:length(numel)
    disp('------------------------')                                      ;%
    disp(['Run ',num2str(sim),': nel = ',num2str(numel(sim)),''])         ;%
    disp('------------------------')                                      ;%
    
    %% NON-DIMENSIONAL CONSTANT
    nu      = 0.0                                                         ;% poisson ratio
    ni      = 2                                                           ;% number of mp in h(1) direction
    nstr    = 4                                                           ;% number of stress components
    %---------------------------------------------------------------------%
    
    %% INDEPENDANT PHYSICAL CONSTANT
    g       = 9.8                                                         ;% gravitationnal acceleration [m/s^2]
    E       = 10e3                                                         ;% young's modulus             [Pa]
    Gc      = E/(2*(1+nu))                                                ;% shear modulus               [Pa]
    Kc      = E/(3*(1-2*nu))                                              ;% bulk modulus                [Pa]
    rho0    = 80                                                        ;% density                     [kg/m^3]
    n0      = 0.25                                                        ;% initial porosity            [-]
    coh0    = 25.0e3                                                      ;% cohesion                    [Pa]
    phi0    = 30.0*pi/180                                                 ;% friction angle              [Rad]
    H       = -60e3                                                       ;% softening modulus           [Pa]
    cohr    = 4.0e3                                                       ;% residual cohesion           [Pa]
    phir    = 10*pi/180                                                  ;% residual friction angle     [Rad]           [s]
    %---------------------------------------------------------------------%
    
    %% MESH INITIALIZATION
    [meD] = meSetupCOMPRESSION(numel(sim),typeD)                                     ;% - see function
    % BOUNDARY CONDITIONS
    du      = 0.0;
    [bc,BC,xB,yB] = setBCs(meD,du)                                        ;%
    %---------------------------------------------------------------------%
    %% MPM DISCRETIZATION
    ly      = 10                                                          ;% layer thickness [m]
    yinf    = 0.0;
    [mpD,p2]= mpSetupCOMPRESSION(meD,ni,ly,xB(1),xB(2),yinf,coh0,cohr,phi0,phir,n0,rho0,nstr,typeD);% - see function
    
    y0 = mpD.x(:,2);
    figure(1),plot(meD.x,meD.y,'s',meD.x(BC),meD.y(BC),'gs',mpD.x(:,1),mpD.x(:,2),'x');axis equal;drawnow
    % ISOTROPIC ELASTIC MATRIX
    Del = [ Kc+4/3*Gc,Kc-2/3*Gc,Kc-2/3*Gc,0.0;...
            Kc-2/3*Gc,Kc+4/3*Gc,Kc-2/3*Gc,0.0;...
            Kc-2/3*Gc,Kc-2/3*Gc,Kc+4/3*Gc,0.0;...
            0.0      ,0.0      ,0.0      ,Gc]                             ;%                                                   ;%
    %---------------------------------------------------------------------%
    
    lstps = 50;
    lstp  = 1;
    fErr    = eps;
    tol     = 1e-9;
    NRitMax = 10;
    D       = zeros(size(mpD.epII));
    LPtimer = zeros(lstps,1);
    
%     y0 = mpD.x(:,2);
    %% MPM MUSL VARIANT EXPLICIT SOLVER FOR INFINITESIMAL STRAIN
    disp(['MPM SOLVER ON: ',num2str(meD.nN),' nodes, ',num2str(mpD.n),' material points']);
    tsolve=tic;
    while( lstp<=lstps )% BEGIN WHILE LOOP
        %% LOADSTEP INCREMENT
        fprintf('\n LOADSTEP %.0f / %.0f : previous NR error %8.2e',lstp,lstps,fErr);
        fprintf('\n ---------------------------------------------- \n');
        time_it= tic                                                      ;%
        dpi    = time_it                                                  ;% CURRENT ITERATION TIMER BEGIN
        %% TRACK MATERIAL POINTS (p) IN ELEMENTS (E)
        c2e  = reshape((permute((floor((mpD.yc-min(meD.y))./meD.h(2))+1)+...
               (meD.nEy).*floor((mpD.xc-min(meD.x))./meD.h(1)),[2 1])),mpD.n*4,1);%
        neon = length(unique(c2e))                                        ;%
        c2N  = reshape((meD.e2N(c2e,:)'),meD.nNe,mpD.n)'                  ;%
        l2g  = [meD.DoF*c2N-1;meD.DoF*c2N]                                ;% local to global node index list [x_I,y_I]

        %------------------------------------------------------------------%
        %% FREE DEGREES OF FREEDOM
        inC    = [meD.DoF*unique(c2N)-1;meD.DoF*unique(c2N)]              ;%
        fd     = 1:meD.nDoF(2)                                            ;%                                                    % zero fixed displacement BCs
        fd     = fd(inC)                                                  ;%  
        
        iDx    = meD.DoF*c2N-1                                            ;%
        iDy    = iDx+1                                                    ;%
        %% BASIS FUNCTIONS
        [mpD,N] = SdS(mpD,meD,c2N);
        %------------------------------------------------------------------%
        %% GET EXTERNAL FORCE
        [fext] = getFext(mpD,meD,g,l2g)                                   ;%
        fext   = fext.*(lstp/lstps) ;
        %% STORE PREVIOUSLY CONVERGED MP PROPERTIES
        vp     = mpD.V;
        sig0   = mpD.s;
        F      = mpD.F;
        %% INITIALIZATION FOR N-R ITERATIONS
        % SCALAR
        fErr   = 1                                                        ;% guess error
        NRit   = 0                                                        ;% zero NR counter
        % VECTOR
        uvw    = zeros(meD.nDoF(2),1)                                     ;% zero displacement over a loadstep
        %% NEWTON-RAPHSON ITERATIONS
        while((fErr > tol) && (NRit < NRitMax+1))              % global equilibrium loop
            %--------------------------------------------------------------%
            [mpD]=incrDef(mpD,uvw,F,meD.DoF,iDx,iDy,nstr)                 ;% update deformation
            [sig] = elasticTrial(mpD,sig0,Del)                            ;% update elastic predictor
            if(plasticity)
                [sig,mpD] = plasticCorrector(sig,mpD,Hp,cohr,Del)         ;% update plastic corrector
            end
            [K,fint] = getStiffness(meD,mpD,sig,vp,c2N,Del)               ;% get global internal force & global stiffness matrix
            oobf     = fext-fint                                          ;% out-of-balance force
            %--------------------------------------------------------------%
            duvw=zeros(meD.nDoF(2),1)                                     ;% zero displacement correction    
            [uvw,duvw,fErr]=solve(K,oobf,uvw,duvw,bc,fd,NRit)             ;% solve for incremental displacement and correction
            %--------------------------------------------------------------%
            NRit=NRit+1                                                   ;% increment the NR counter
            if(NRit>1)
               fprintf('  iteration %.0f , NR error %8.2e \n',NRit-1,fErr);% terminal display
            end
        end
        
        %% UPDATE MATERIAL POINT DOMAIN AND COORDINATE
        C=[1:4;5:8;9:12;13:16];
        % CORNERS COORDINATE UPDATE
        for c=1:4
            nc = C(c,:);
            mpD.xc(:,c) = mpD.xc(:,c) + sum(N(:,nc).*uvw(meD.DoF.*c2N(:,nc)-1),2);
            mpD.yc(:,c) = mpD.yc(:,c) + sum(N(:,nc).*uvw(meD.DoF.*c2N(:,nc)  ),2);
        end
        % MP'S COORDINATE UPDATE
        mpD.x = [mean(mpD.xc,2) mean(mpD.yc,2)];
        mpD.u = mpD.u+[sum(mpD.S.*uvw(iDx),2) sum(mpD.S.*uvw(iDy),2)];
        % MP'S VOLUME UPDATE
        mpD.V = ((mpD.xc(:,1).*mpD.yc(:,2)-mpD.xc(:,2).*mpD.yc(:,1))+...
                 (mpD.xc(:,2).*mpD.yc(:,3)-mpD.xc(:,3).*mpD.yc(:,2))+...
                 (mpD.xc(:,3).*mpD.yc(:,4)-mpD.xc(:,4).*mpD.yc(:,3))+...
                 (mpD.xc(:,4).*mpD.yc(:,1)-mpD.xc(:,1).*mpD.yc(:,4))).*0.5;
        % MATERIAL POINT UPDATE
        mpD.s    = sig                                                    ;% stress 
        % MP'S SECOND INVARIANT PLASTIC STRAIN
        mpD.epII = mpD.epII+mpD.depII                                     ;% second invariant of the plastic strain

        %------------------------------------------------------------------%
        %% TERMINAL DISPLAY
        LPtimer(lstp)=toc(dpi)                                            ;% CURRENT ITERATION TIMER END
        timer = mean(LPtimer(1:lstp))*(lstps-lstp);
        fprintf('-----------------------------------\n');
        fprintf('Remaining time : %.2f seconds \n',timer);
        %------------------------------------------------------------------%
        lstp = lstp+1                                                     ;%
    end% END WHILE LOOP                                                   ;%
    tsolve = toc(tsolve)                                                  ;%
%    print(gcf,[datapath version 'vertical_stress'],'-dpng');
    
    
    name=['.\data\iCPDI_data_' num2str(sim) '.mat'];
    save(name,'mpD','meD','y0','ni','rho0','E','nu','g','tsolve');
    
    fig3=figure(3);
    set(fig3,'Units','pixels','Position',[139 98.3333 504 458.6667]);
    pp.ps = 10;
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

