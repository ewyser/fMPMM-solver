% fMPMM-solver: Fast and efficient MATLAB-based MPM solver
%
% Copyright (C) 2020  Emmanuel Wyser, Michel Jaboyedoff, Yury Podladchikov
% -------------------------------------------------------------------------%
% version    : v1.0
% date       : february, 2020
% description: implicit mpm (CPDI2q) solver based on an updated Lagrangian
% frame in plane strain condition for elasto-plastic problem discretized on
% a 4-noded quadrilateral mesh
% -------------------------------------------------------------------------%
%% FANCY CALLS
clear                                                                     ;%
addpath('functions')                                                      ;%   
version    = 'GIMP_EP_'                                                    ;%
plasticity = false                                                         ;%         
typeD='double'                                                            ;%
disp('EXACT PLASTIC CORRECTION')                                          ;%
disp('------------------------')                                          ;%

set(0,'defaulttextinterpreter','latex')                                   ;%
fslab  = 14; fsleg  = 14; fstit  = 14; fstick = 14;

numel = [1 2 5 10 20 40 80 160 320 640];
run   = zeros(length(numel),6)                                            ;%                                ;%
for sim=1:length(numel)
    disp('------------------------')                                      ;%
    disp(['Run ',num2str(sim),': nel = ',num2str(numel(sim)),''])         ;%
    disp('------------------------')                                      ;%
    
    %% NON-DIMENSIONAL CONSTANT
    ar      = 0.8                                                         ;% aspect ratio thickness/width
    nu      = 0.0                                                         ;% poisson ratio
    ni      = 2                                                           ;% number of mp in h(1) direction
    nstr    = 4                                                           ;% number of stress components
    %---------------------------------------------------------------------%
    
    %% INDEPENDANT PHYSICAL CONSTANT
    g       = 9.81                                                        ;% gravitationnal acceleration [m/s^2]
    E       = 1.0e6                                                       ;% young's modulus             [Pa]
    Gc      = E/(2*(1+nu))                                                ;% shear modulus               [Pa]
    Kc      = E/(3*(1-2*nu))                                              ;% bulk modulus                [Pa]
    rho0    = 80                                                          ;% density                     [kg/m^3]
    n0      = 0.25                                                        ;% initial porosity            [-]
    yd      = sqrt(E/rho0)                                                ;% elastic wave velocity       [m/s]
    coh0    = 20.0e3                                                      ;% cohesion                    [Pa]
    phi0    = 30.0*pi/180                                                 ;% friction angle              [Rad]
    H       = -60e3                                                       ;% softening modulus           [Pa]
    cohr    =  4.0e3                                                      ;% residual cohesion           [Pa]
    phir    = 7.5*pi/180                                                  ;% residual friction angle     [Rad]
    l0      = 50                                                          ;% initial column height
    t       = ceil((1/yd)*(2*l0)*40)                                                        ;% simulation time             [s]
    te      = ceil((1/yd)*(2*l0)*40)                                                         ;% elastic loading             [s]    %---------------------------------------------------------------------%
    
    %% MESH INITIALIZATION
    [meD] = meSetupElasticLoading(numel(sim),typeD,l0)                                     ;% - see function   
    % BOUNDARY CONDITIONS
    yinf    = 0.0                                                         ;%
    xB      = [-meD.L(1)/2 meD.L(1)/2]                                    ;%
    [row]   = find(meD.y<=yinf & meD.y >yinf-meD.h(2))                    ;%
    BC.yi   = row                                                         ;%
    [row]   = find(meD.x<=xB(1))                                          ;%
    BC.xi   = row                                                         ;%   
    [row]   = find(meD.x>=xB(2))                                          ;%
    BC.xs   = row                                                         ;%   
    clear row                                                             ;%
    %---------------------------------------------------------------------%
    
    %% MPM DISCRETIZATION                                                 ;% layer thickness [m]
    [mpD,p2]= mpSetupElasticLoading(meD,ni,l0,xB(1),xB(2),yinf,coh0,cohr,phi0,phir,n0,rho0,nstr,typeD);% - see function
    % ISOTROPIC ELASTIC MATRIX
    Del     = [ Kc+4/3*Gc,Kc-2/3*Gc,Kc-2/3*Gc,0.0;...
                Kc-2/3*Gc,Kc+4/3*Gc,Kc-2/3*Gc,0.0;...
                Kc-2/3*Gc,Kc-2/3*Gc,Kc+4/3*Gc,0.0;...
                0.0      ,0.0      ,0.0      ,Gc]                         ;%
    Hp      = H*meD.h(1)                                                  ;%
    %---------------------------------------------------------------------%
    
    %% DISPLAY PARAMETERS AND RUNTIME INITIALIZATION
    fps  = 25                                                             ;% image per second
    %COURANT-FRIEDRICH-LEVY CONDITION
    dt   = 0.5*meD.h(1)/yd                                                ;% unconditionally stable timestep
    nit  = ceil(t/dt)                                                     ;% maximum number of interation
    nf   = max(2,ceil(round(1/dt)/fps))                                   ;% number of frame interval
    % GRAVITY INCREMENT
    dg   = linspace(0,g,round(((te))/dt))                             ;% gravity increase
    %STORAGE DATA ALLOCATION
    runt = zeros(ceil(nit/nf),2)                                         ;% computational activity
    % RUNTIME PARAMETERS
    nc   = 0                                                              ;% initialize iteration counter                                                         
    it   = 1                                                              ;% initialize iteration
    tw   = 0.0                                                            ;% initialize time while statement
    cycle_time = zeros(nit,7)                                             ;%
    % PLOT SETTING
    [pp] = plotSetup(xB(1),xB(2),yinf,l0,meD.y,BC.yi)                     ;% - see function
    %---------------------------------------------------------------------%
    y0 = mpD.x(:,2);
    %% MPM MUSL VARIANT EXPLICIT SOLVER FOR INFINITESIMAL STRAIN
    disp(['MPM SOLVER ON: ',num2str(meD.nN),' nodes, ',num2str(mpD.n),' material points']);
    tsolve=tic                                                            ;
    while(tw<t)% BEGIN WHILE LOOP
        time_it= tic                                                      ;%
        dpi    = time_it                                                  ;% CURRENT ITERATION TIMER BEGIN
        %% LINEAR INCREASE OF GRAVITY FOR EQUILIBRIUM
        if(it<=size(dg,2))
            g = dg(it)                                                    ;% incremented gravity load
        end
        %------------------------------------------------------------------%
% %          %% ADAPTATIVE MP DOMAIN
% %         if(it*dt>te)
% %             l = 1.25;
% %             [mpD] = f_mpDomain(l,mpD);
% %             mpD.n = length(mpD.x);
% %         end
% %         %------------------------------------------------------------------%

        %% TRACK MATERIAL POINTS (p) IN ELEMENTS (E)
        p2.e = (floor((mpD.x(:,2)-min(meD.y))./meD.h(2))+1)+...
               (meD.nEy).*floor((mpD.x(:,1)-min(meD.x))./meD.h(1))        ;% index mp to element
        p2.N = meD.e2N(p2.e,:)                                            ;% index mp (to element) to node 
        l2g  = [meD.DoF*p2.N-1;meD.DoF*p2.N]                              ;% local to global node index list [x_I,y_I]
        neon = size(unique(p2.e),1)                                       ;% number of active element
        %------------------------------------------------------------------%
        
        %% BASIS FUNCTIONS
        tic;
        [mpD] = NdN(meD,mpD,p2.N)                                         ;% - see function 
        cycle_time(it,1)=toc;
        %------------------------------------------------------------------%
        %% MAPPING FROM MATERIAL POINTS (p) TO NODES (N)
        tic
        [meD] = p2Nsolve(meD,mpD,g,dt,l2g,p2.N,BC)                        ;% - see function
        cycle_time(it,2)=toc; 
        %------------------------------------------------------------------%

        %% MAPPING FROM NODES (N) TO MATERIAL POINTS (p)
        tic
        [meD,mpD] = mapN2p(meD,mpD,dt,l2g,p2.N,BC)                        ;% - see function
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
        [mpD] = elastic(mpD,Del)                                          ;% - see function
        % PLASTIC CORRECTION: EXACT SOLUTION f<1E-11
        if((plasticity)&&(dt*it>te))
            [mpD,dat] = plastic(mpD,Hp,cohr,Del)                          ;% - see function
        end
        cycle_time(it,5)=toc                                              ;%
        dpi=toc(dpi)                                                      ;% CURRENT ITERATION TIMER END
        cycle_time(it,6)=dpi                                              ;%
        cycle_time(it,7)=length(unique(p2.N))                             ;%
        %------------------------------------------------------------------%
        
        %% TERMINAL DISPLAY
        if(mod(it,nf)==1)
            nc    = nc+1                                                  ;%
            clocktimer(((nit-it)*toc(time_it)),'Remaining estimated time:',1/dpi);%
            runt(nc,:) = [1/dpi mean(runt(1:nc,1))]                       ;%
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
    clocktimer(tsolve,'Runtime MPM solver:',mean(runt(:,2)))              ;%


    name=['.\data\GIMP_time_vectorized_' num2str(sim) '.mat'];
    save(name,'runt','cycle_time');
    
    name=['.\data\GIMP_data_' num2str(sim) '.mat'];
    save(name,'mpD','meD','y0','ni','rho0','E','nu','g');
end