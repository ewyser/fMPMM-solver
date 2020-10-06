function [mpD] = DefUpdate(meD,mpD,l2g,c2N)
% f_DEFUPDATE : function which updates deformation-related information at
% the material point level, i.e., vorticity, strain-rate, incremental and
% current deformation gradient and material point properties (dimensions)
%   Detailed explanation goes here 
%% UPDATE VORTICITY & STRAIN-RATE
ix    = l2g(1:mpD.n    ,:)                                                ;% pre-computation of x index
iy    = l2g(mpD.n+1:end,:)                                                ;% pre-computation of y index
% INITIALIZE STRAIN-RATE
mpD.U               = zeros(mpD.n,size(mpD.B,2))                          ;%
mpD.U(:,1:meD.DoF:end) = meD.u(ix)                                           ;% nodal x-displacement
mpD.U(:,2:meD.DoF:end) = meD.u(iy)                                           ;% nodal y-displacement
for mp = 1:mpD.n
    iDx = meD.DoF*c2N(mp,:)-1                                                ;% x-component index
    iDy = iDx+1                                                           ;% y-component index
    % STRAIN-RATE
    mpD.e(:,mp) = mpD.B(:,:,mp)*mpD.U(mp,:)';
     %SPIN
    mpD.w(mp,3)= 0.5*(mpD.dSy(mp,:)*meD.u(iDx)-mpD.dSx(mp,:)*meD.u(iDy))  ;% mp spin tensor xy
    % INCREMENTAL DEFORMATION UPDATE
    mpD.dF(mp,1) = 1.0+(mpD.dSx(mp,:)*meD.u(iDx))                         ;% mp incremental deformation tensor xx
    mpD.dF(mp,2) = 0.0+(mpD.dSy(mp,:)*meD.u(iDx))                         ;% mp incremental deformation tensor xy
    mpD.dF(mp,3) = 0.0+(mpD.dSx(mp,:)*meD.u(iDy))                         ;% mp incremental deformation tensor yx
    mpD.dF(mp,4) = 1.0+(mpD.dSy(mp,:)*meD.u(iDy))                         ;% mp incremental deformation tensor yy
end

%% DEFORMATION GRADIENT UPDATE
mpD.F(:,1)  = mpD.dF(:,1).*mpD.F(:,1)+mpD.dF(:,2).*mpD.F(:,3)             ;% deformation gradient xx
mpD.F(:,2)  = mpD.dF(:,1).*mpD.F(:,2)+mpD.dF(:,2).*mpD.F(:,4)             ;% deformation gradient xy
mpD.F(:,3)  = mpD.dF(:,3).*mpD.F(:,1)+mpD.dF(:,4).*mpD.F(:,3)             ;% deformation gradient yx
mpD.F(:,4)  = mpD.dF(:,3).*mpD.F(:,2)+mpD.dF(:,4).*mpD.F(:,4)             ;% deformation gradient yy
%--------------------------------------------------------------------------%

%% UPDATE MATERIAL POINT DOMAIN
mpD.r1(:,1) = mpD.F(:,1).*mpD.r01(:,1)+mpD.F(:,2).*mpD.r01(:,2)           ;% update domain vector 1
mpD.r1(:,2) = mpD.F(:,3).*mpD.r01(:,1)+mpD.F(:,4).*mpD.r01(:,2)           ;% update domain vector 1
mpD.r2(:,1) = mpD.F(:,1).*mpD.r02(:,1)+mpD.F(:,2).*mpD.r02(:,2)           ;% update domain vector 2
mpD.r2(:,2) = mpD.F(:,3).*mpD.r02(:,1)+mpD.F(:,4).*mpD.r02(:,2)           ;% update domain vector 3
r1          = [-0.5 0.5 0.5 -0.5]                                         ;%
r2          = [-0.5 -0.5 0.5 0.5]                                         ;%
mpD.xc      = repmat(mpD.x(:,1),1,4)+mpD.r1(:,1).*r1+mpD.r2(:,1).*r2      ;% update corners x coordinates
mpD.yc      = repmat(mpD.x(:,2),1,4)+mpD.r1(:,2).*r1+mpD.r2(:,2).*r2      ;% update corners y coordinates
%--------------------------------------------------------------------------%

%% MATERIAL DOMAIN, VOLUME & POROSITY UPDATE
V     = cross([mpD.r1 zeros(mpD.n,1)],[mpD.r2 zeros(mpD.n,1)])            ;% update mp volume
mpD.V = V(:,end)                                                          ;% save new volume
mpD.J = mpD.V./mpD.V0                                                     ;% determinant of the deformation gradient tensor 
mpD.n0= 1-(1-mpD.n0)./mpD.J                                               ;% upadte porosity

end

