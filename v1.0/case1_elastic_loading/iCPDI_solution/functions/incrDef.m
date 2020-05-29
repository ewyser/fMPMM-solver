function [mpD] = incrDef(mpD,uvw,F,DoF,iDx,iDy,nstr)
% f_DEFUPDATE : function which updates deformation-related information at
% the material point level, i.e., vorticity, strain-rate, incremental and
% current deformation gradient and material point properties (dimensions)
%   Detailed explanation goes here
%% SPIN UPDATE
mpD.w(:,3)= 0.5*sum((mpD.dSy.*uvw(iDx)-mpD.dSx.*uvw(iDy)),2)          ;% spin in (x,y) plane
%--------------------------------------------------------------------------%
%% STRAIN UPDATE
mpD.U(:)           = 0.0                                              ;% initialize nodal displacement vector
mpD.U(:,1:DoF:end) = uvw(iDx)                                       ;% nodal x-displacement
mpD.U(:,2:DoF:end) = uvw(iDy)                                       ;% nodal y-displacement
mpD.e = permute(sum(mpD.B.*repmat(permute(mpD.U,[3 2 1]),nstr,1),2),[1 3 2]);% mp strain
%--------------------------------------------------------------------------%
%% INCREMENTAL DEFORMATION GRADIENT UPDATE deformation gradient increment
mpD.dF(:,1) = 1.0+sum(mpD.dSx.*uvw(iDx),2)                              ;% incremental deformation gradient xx
mpD.dF(:,2) = 0.0+sum(mpD.dSy.*uvw(iDx),2)                              ;% incremental deformation gradient xy
mpD.dF(:,3) = 0.0+sum(mpD.dSx.*uvw(iDy),2)                              ;% incremental deformation gradient yx
mpD.dF(:,4) = 1.0+sum(mpD.dSy.*uvw(iDy),2)                              ;% incremental deformation gradient yy
%--------------------------------------------------------------------------%

%% DEFORMATION GRADIENT UPDATE
mpD.F(:,1)  = mpD.dF(:,1).*F(:,1)+mpD.dF(:,2).*F(:,3)                     ;% deformation gradient xx
mpD.F(:,2)  = mpD.dF(:,1).*F(:,2)+mpD.dF(:,2).*F(:,4)                     ;% deformation gradient xy
mpD.F(:,3)  = mpD.dF(:,3).*F(:,1)+mpD.dF(:,4).*F(:,3)                     ;% deformation gradient yx
mpD.F(:,4)  = mpD.dF(:,3).*F(:,2)+mpD.dF(:,4).*F(:,4)                     ;% deformation gradient yy
%--------------------------------------------------------------------------%

%% MATERIAL DOMAIN, VOLUME & POROSITY UPDATE
mpD.J   = mpD.dF(:,1).*mpD.dF(:,4)-mpD.dF(:,3).*mpD.dF(:,2)                 ;% determinant of the deformation gradient increment
clear iDx iDy                                                             ;% clear temporary variables
%--------------------------------------------------------------------------%

end

