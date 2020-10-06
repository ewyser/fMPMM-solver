function [mpD] = DefUpdate_loop(meD,mpD,dt,l2g,p2n,nDoF)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
uGIMP = 0; cpGIMP = 1;
%% UPDATE VORTICITY & STRAIN-RATE
ix    = l2g(1:mpD.n    ,:)                                                ;% pre-computation of x index
iy    = l2g(mpD.n+1:end,:)                                                ;% pre-computation of y index
% INITIALIZE STRAIN-RATE
mpD.U               = zeros(mpD.n,size(mpD.B,2))                          ;%
mpD.U(:,1:meD.DoF:end) = meD.u(ix)                                           ;% nodal x-displacement
mpD.U(:,2:meD.DoF:end) = meD.u(iy)                                           ;% nodal y-displacement
for mp = 1:mpD.n
    iDx = meD.DoF*p2n(mp,:)-1                                                ;% x-component index
    iDy = iDx+1                                                           ;% y-component index
    % STRAIN-RATE
    mpD.e(:,mp) = mpD.B(:,:,mp)*mpD.U(mp,:)';
     %SPIN
    mpD.w(mp,3)= 0.5*(mpD.dSy(mp,:)*meD.u(iDx)-mpD.dSx(mp,:)*meD.u(iDy))  ;% mp spin tensor xy
    % INCREMENTAL DEFORMATION UPDATE
    mpD.dD(mp,1) = 1.0+(mpD.dSx(mp,:)*meD.u(iDx))                         ;% mp incremental deformation tensor xx
    mpD.dD(mp,2) = 0.0+(mpD.dSy(mp,:)*meD.u(iDx))                         ;% mp incremental deformation tensor xy
    mpD.dD(mp,3) = 0.0+(mpD.dSx(mp,:)*meD.u(iDy))                         ;% mp incremental deformation tensor yx
    mpD.dD(mp,4) = 1.0+(mpD.dSy(mp,:)*meD.u(iDy))                         ;% mp incremental deformation tensor yy
end
% DEFORMATION GRADIENT UPDATE
mpD.F(:,1)  = mpD.dD(:,1).*mpD.F(:,1)+mpD.dD(:,2).*mpD.F(:,3)             ;
mpD.F(:,2)  = mpD.dD(:,1).*mpD.F(:,2)+mpD.dD(:,2).*mpD.F(:,4)             ;
mpD.F(:,3)  = mpD.dD(:,3).*mpD.F(:,1)+mpD.dD(:,4).*mpD.F(:,3)             ;
mpD.F(:,4)  = mpD.dD(:,3).*mpD.F(:,2)+mpD.dD(:,4).*mpD.F(:,4)             ;
% MATERIAL DOMAIN, VOLUME & POROSITY UPDATE
mpD.J = mpD.dD(:,1).*mpD.dD(:,4)-mpD.dD(:,3).*mpD.dD(:,2)                 ;% determinant of the incremental deformation tensor
mpD.V = mpD.J.*mpD.V                                                      ;% volume of mp
if(uGIMP==1)
    mpD.l = mpD.l0                                                        ;%
elseif(cpGIMP==1)
    %mpD.l = mpD.J.*mpD.l                                                  ;% determinant
    %mpD.l = mpD.dD(:,[1 4]).*mpD.l                                        ;% principal component of dF
    mpD.dU= sqrt([(mpD.dD(:,1).^2+mpD.dD(:,3).^2),...
        (mpD.dD(:,2).^2+mpD.dD(:,4).^2)])                                 ;% stretch part of dF, see Charlton etal, 2017
    mpD.l = mpD.dU.*mpD.l                                                 ;% half-length of mp
end
end

