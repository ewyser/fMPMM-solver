function [mpD] = DefUpdate(meD,mpD,l2g)
% f_DEFUPDATE : function which updates deformation-related information at
% the material point level, i.e., vorticity, strain-rate, incremental and
% current deformation gradient and material point properties (dimensions)
%   Detailed explanation goes here 
%% PREPROCESSING
iDx    = l2g(1:mpD.n    ,:)                                               ;% x component global node index
iDy    = l2g(mpD.n+1:end,:)                                               ;% y component global node index
nstr   = size(mpD.B,1)                                                    ;% number of stress component
%--------------------------------------------------------------------------%

%% SPIN UPDATE
mpD.w(:,3)= 0.5*sum((mpD.dSy.*meD.u(iDx)-mpD.dSx.*meD.u(iDy)),2)          ;% spin in (x,y) plane
%--------------------------------------------------------------------------%

%% STRAIN CALCULATION
mpD.U(:)               = 0.0                                              ;% initialize nodal displacement vector
mpD.U(:,1:meD.DoF:end) = meD.u(iDx)                                       ;% nodal x-displacement
mpD.U(:,2:meD.DoF:end) = meD.u(iDy)                                       ;% nodal y-displacement
mpD.e = permute(sum(mpD.B.*repmat(permute(mpD.U,[3 2 1]),nstr,1),2),[1 3 2]);% mp strain
%--------------------------------------------------------------------------%

%% INCREMENTAL DEFORMATION GRADIENT UPDATE deformation gradient increment
mpD.dF(:,1) = 1.0+sum(mpD.dSx.*meD.u(iDx),2)                              ;% incremental deformation gradient xx
mpD.dF(:,2) = 0.0+sum(mpD.dSy.*meD.u(iDx),2)                              ;% incremental deformation gradient xy
mpD.dF(:,3) = 0.0+sum(mpD.dSx.*meD.u(iDy),2)                              ;% incremental deformation gradient yx
mpD.dF(:,4) = 1.0+sum(mpD.dSy.*meD.u(iDy),2)                              ;% incremental deformation gradient yy
%--------------------------------------------------------------------------%

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

