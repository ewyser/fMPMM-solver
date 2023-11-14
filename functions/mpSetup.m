function [mpD] = mpSetup(meD,ni,ly,coh0,phi0,n0,rho0,nstr,typeD)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% MPM DISCRETIZATION
% wall features
wall_width       = 10.0
wall_location    = 60.0
wall_left_coord  = wall_location-wall_width
wall_right_coord = wall_location+wall_width

% layer initialization
xinf    = meD.xB(1);
xsup    = meD.xB(2);
yinf    = meD.xB(3)
xL      = xinf+(0.5*meD.h(1)/ni):meD.h(1)/ni:xsup   ;
yL      = yinf+(0.5*meD.h(2)/ni):meD.h(2)/ni:  ly-(0.5*meD.h(2)/ni)       ;
[xl,yl] = meshgrid(xL,yL)                                                 ;
wl      = 15                                                         ; % 0.15


xl = xl(:); 
yl = yl(:); 

x=linspace(min(xl),max(xl),200);
a= -1.0;
y= a.*x+ly/4;
x= x+meD.L(1)/2;

xlt = [];
ylt = [];
cp = [];
isplast = [];
for mp=1:length(xl)
    for p = 1:length(y)
        DX = xl(mp)-x(p);
        DY = yl(mp)-y(p);
        nx = a;
        ny = -1;
        s = DX*nx+DY*ny;
        if(s>0)
            pos = 1;
        else
            pos = 0;
        end
        if(yl(mp)<wl)
            pos = 1;
        end
    end
    if(pos==1)
        xlt     = [xlt xl(mp)];
        ylt     = [ylt yl(mp)];
        isplast = [isplast true];

    end
    if(xl(mp)>wall_left_coord && xl(mp)<wall_right_coord)
        xlt     = [xlt xl(mp)];
        ylt     = [ylt yl(mp)];
        isplast = [isplast false];
    end
end
xl = xlt;
yl = ylt;

%% MATERIAL POINT:
% SCALAR & VECTOR
mpD.n  = length(xl(:))                                                    ;% number of material point
mpD.n0 = ones(mpD.n,1,typeD).*n0                                          ;% porosity
mpD.l0 = ones(mpD.n,2,typeD).*(meD.h(1)/ni)./2.0                          ;% reference domain dimension

mpD.r01 = [2*mpD.l0(:,1)    zeros(mpD.n,1)];
mpD.r02 = [zeros(mpD.n,1) 2*mpD.l0(:,1)   ];

mpD.r1 = mpD.r01;
mpD.r2 = mpD.r02;

mpD.V0 = ones(mpD.n,1,typeD).*(2.*mpD.l0(:,1).*2.*mpD.l0(:,2))            ;% reference volume
mpD.V  = mpD.V0                                                           ;% current volume
mpD.m  = rho0.*mpD.V                                                      ;% mass
mpD.x  = [xl(:) yl(:)]                                                    ;% coordinate
mpD.u  = zeros(mpD.n,2,typeD)                                             ;% displacement
mpD.v  = zeros(mpD.n,2,typeD)                                             ;% velocity
mpD.p  = zeros(mpD.n,2,typeD)                                             ;% momentum
mpD.coh= coh0.*ones(1,1,typeD)                                            ;% cohesion
mpD.phi= ones(1,mpD.n,typeD).*phi0;                                       ;% friction
mpD.J  = ones(mpD.n,1,typeD)                                              ;% determinant of deformation gradient
mpD.P  = zeros(1,mpD.n,typeD)                                             ;% pressure
mpD.epV  = zeros(1,mpD.n,typeD)                                           ;% volumetric plastic strain
mpD.epII = zeros(1,mpD.n,typeD)                                           ;% second invariant of the deviatoric plastic strain
% TENSOR
mpD.w  = zeros(mpD.n,3,typeD)                                             ;% [x,y,z] axis spin
mpD.U  = zeros(mpD.n,meD.nDoF(1),typeD)                                    ;% strain
mpD.B  = zeros(nstr,meD.nDoF(1),mpD.n)                                        ;% strain-displacement matrix

mpD.dD = zeros(mpD.n,4,typeD)                                             ;% incremental deformation gradient
mpD.F  = zeros(mpD.n,4,typeD)+repmat([1 0 0 1],mpD.n,1)                   ;% deformation gradient
mpD.b  = zeros(mpD.n,4,typeD)+repmat([1 0 0 1],mpD.n,1)                   ;% left cauchy green deformation gradient

mpD.e  = zeros(nstr,mpD.n,typeD)                                          ;% strain
mpD.s  = zeros(nstr,mpD.n,typeD)                                          ;% stress tensor
mpD.sn = zeros(nstr,mpD.n,typeD)                                          ;% stress tensor
mpD.tau = zeros(nstr,mpD.n,typeD)                                         ;% deviatoric stess
mpD.ep = zeros(nstr,mpD.n,typeD)                                          ;% plastic strain
mpD.devep = zeros(nstr,mpD.n,typeD)                                       ;% deviatoric plastic strain
% ADDITIONAL QUANTITES
mpD.xc = repmat(mpD.x(:,1),1,4)+mpD.r1(:,1).*[-0.5 0.5 0.5 -0.5]+...
    mpD.r2(:,1).*[-0.5 -0.5 0.5 0.5];
mpD.yc = repmat(mpD.x(:,2),1,4)+mpD.r1(:,2).*[-0.5 0.5 0.5 -0.5]+...
    mpD.r2(:,2).*[-0.5 -0.5 0.5 0.5];
mpD.isPlast = isplast

end

