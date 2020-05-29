function [mpD,p2] = mpSetup(meD,ni,ly,xinf,xsup,yinf,coh0,cohr,phi0,phir,n0,rho0,nstr,typeD)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% MPM DISCRETIZATION
% layer initialization
xL      = xinf+(0.5*meD.h(1)/ni):meD.h(1)/ni:xsup   ;
yL      = yinf+(0.5*meD.h(2)/ni):meD.h(2)/ni:  ly-(0.5*meD.h(2)/ni)       ;
[xl,yl] = meshgrid(xL,yL)                                                 ;
dx      = meD.h(1)/0.5;
dy      = meD.h(2)/0.5;
%% MATERIAL POINT:
% SCALAR & VECTOR
mpD.n  = length(xl(:))                                                    ;% number of material point
PF     = 0.99                                                             ;% PIC-FLIP update
mpD.PF = [PF 1-PF]                                                        ;% -
mpD.n0 = ones(mpD.n,1,typeD).*n0                                          ;% porosity
mpD.l0 = ones(mpD.n,2,typeD).*(meD.h(1)/ni)./2                            ;% reference domain dimension
mpD.l  = mpD.l0                                                           ;% current domain dimension
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
mpD.dD = zeros(mpD.n,4,typeD)                                             ;% incremental deformation gradient
mpD.F  = zeros(mpD.n,4,typeD)+repmat([1 0 0 1],mpD.n,1)                   ;% deformation gradient
mpD.b  = zeros(mpD.n,4,typeD)+repmat([1 0 0 1],mpD.n,1)                   ;% left cauchy green deformation gradient
mpD.dU = zeros(mpD.n,meD.nDF(1),typeD)                                    ;% stretch deformation gradient
mpD.e  = zeros(nstr,mpD.n,typeD)                                          ;% strain
mpD.s  = zeros(nstr,mpD.n,typeD)                                          ;% stress tensor
mpD.sn = zeros(nstr,mpD.n,typeD)                                          ;% stress tensor
mpD.tau = zeros(nstr,mpD.n,typeD)                                         ;% deviatoric stess
mpD.ep = zeros(nstr,mpD.n,typeD)                                          ;% plastic strain
mpD.devep = zeros(nstr,mpD.n,typeD)                                       ;% deviatoric plastic strain
% ADDITIONAL QUANTITES
mpD.N   = zeros(mpD.n,meD.nNe,typeD)                                      ;% basis function
mpD.dNx = zeros(size(mpD.N),typeD)                                        ;% basis function x-gradient  
mpD.dNy = zeros(size(mpD.N),typeD)                                        ;% basis function y-gradient  
mpD.U   = zeros(mpD.n,meD.nDF(1),typeD)                                   ;% velocity vector 
mpD.B   = zeros(nstr,meD.nDF(1),mpD.n,typeD)                              ;% B matrix 
% CONNECTIVITY ARRAY
p2.e    = zeros(mpD.n,1      ,'uint32')                                   ;
p2.N    = zeros(mpD.n,meD.nNe,'uint32')                                   ;

end

