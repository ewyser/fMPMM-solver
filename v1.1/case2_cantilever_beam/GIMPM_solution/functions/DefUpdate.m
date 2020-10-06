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

%% MATERIAL DOMAIN, VOLUME & POROSITY UPDATE
mpD.J = mpD.F(:,1).*mpD.F(:,4)-mpD.F(:,3).*mpD.F(:,2)                     ;% determinant of the deformation gradient increment
mpD.V = mpD.J.*mpD.V0                                                     ;% volume
if(uGIMP==1)
    mpD.l = mpD.l0                                                        ;% no update method
elseif(cpGIMP==1)
%    mpD.l = repmat(mpD.J.^(meD.DoF),1,size(mpD.l,2)).*mpD.l0               ;% update based on the determinant of dF
     mpD.dU= sqrt([(mpD.F(:,1).^2+mpD.F(:,3).^2),...
         (mpD.F(:,2).^2+mpD.F(:,4).^2)])                                 ;% update based on the diagonal components of the incremental stretch deformation dU, see Charlton etal, 2017
     mpD.l = mpD.dU.*mpD.l0                                                 ;% -   
%      mpD.l = mpD.F(:,[1 4]).*mpD.l0                                        ;% update based on the diagonal components of the incremental deformation dF, see Charlton etal, 2017
end
clear iDx iDy                                                             ;% clear temporary variables
%--------------------------------------------------------------------------%

end

