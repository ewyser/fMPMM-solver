function [mpD] = DefUpdate(meD,mpD,l2g)
% f_DEFUPDATE : function which updates deformation-related information at
% the material point level, i.e., vorticity, strain-rate, incremental and
% current deformation gradient and material point properties (dimensions)
%   Detailed explanation goes here
uGIMP = 0; cpGIMP = 1;
%% PREPROCESSING
iDx    = l2g(1:mpD.n    ,:)                                               ;% x component global node index
iDy    = l2g(mpD.n+1:end,:)                                               ;% y component global node index
nstr   = size(mpD.B,1)                                                    ;% number of stress component
%--------------------------------------------------------------------------%

%% SPIN UPDATE
mpD.w(:,3)= 0.5*sum((mpD.dNy.*meD.u(iDx)-mpD.dNx.*meD.u(iDy)),2)          ;% spin in (x,y) plane
%--------------------------------------------------------------------------%

%% STRAIN CALCULATION
mpD.U(:)               = 0.0                                              ;% initialize nodal displacement vector
mpD.U(:,1:meD.DoF:end) = meD.u(iDx)                                       ;% nodal x-displacement
mpD.U(:,2:meD.DoF:end) = meD.u(iDy)                                       ;% nodal y-displacement
mpD.e = permute(sum(mpD.B.*repmat(permute(mpD.U,[3 2 1]),nstr,1),2),[1 3 2]);% mp strain
%--------------------------------------------------------------------------%

%% INCREMENTAL DEFORMATION GRADIENT UPDATE deformation gradient increment
mpD.dD(:,1) = 1.0+sum(mpD.dNx.*meD.u(iDx),2)                              ;% deformation gradient increment xx
mpD.dD(:,2) = 0.0+sum(mpD.dNy.*meD.u(iDx),2)                              ;% deformation gradient increment xy
mpD.dD(:,3) = 0.0+sum(mpD.dNx.*meD.u(iDy),2)                              ;% deformation gradient increment yx
mpD.dD(:,4) = 1.0+sum(mpD.dNy.*meD.u(iDy),2)                              ;% deformation gradient increment yy
%--------------------------------------------------------------------------%

%% DEFORMATION GRADIENT UPDATE
mpD.F(:,1)  = mpD.dD(:,1).*mpD.F(:,1)+mpD.dD(:,2).*mpD.F(:,3)             ;% deformation gradient xx
mpD.F(:,2)  = mpD.dD(:,1).*mpD.F(:,2)+mpD.dD(:,2).*mpD.F(:,4)             ;% deformation gradient xy
mpD.F(:,3)  = mpD.dD(:,3).*mpD.F(:,1)+mpD.dD(:,4).*mpD.F(:,3)             ;% deformation gradient yx
mpD.F(:,4)  = mpD.dD(:,3).*mpD.F(:,2)+mpD.dD(:,4).*mpD.F(:,4)             ;% deformation gradient yy
%--------------------------------------------------------------------------%

%% MATERIAL DOMAIN, VOLUME & POROSITY UPDATE
mpD.J = mpD.dD(:,1).*mpD.dD(:,4)-mpD.dD(:,3).*mpD.dD(:,2)                 ;% determinant of the deformation gradient increment
mpD.V = mpD.J.*mpD.V                                                      ;% volume
if(uGIMP==1)
    %mpD.l = repmat(mpD.J,1,size(mpD.l,2)).*mpD.l                          ;% half-length
elseif(cpGIMP==1)
    mpD.dU= sqrt([(mpD.dD(:,1).^2+mpD.dD(:,3).^2),...
        (mpD.dD(:,2).^2+mpD.dD(:,4).^2)])                                 ;% diagonal components of the incremental stretch deformation, see Charlton etal, 2017
    mpD.l = mpD.dU.*mpD.l                                                 ;% -
end
clear iDx iDy                                                             ;% clear temporary variables
%--------------------------------------------------------------------------%

end

