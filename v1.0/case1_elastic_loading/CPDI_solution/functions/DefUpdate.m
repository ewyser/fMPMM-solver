function [mpD] = DefUpdate(meD,mpD,mp,no)
% f_DEFUPDATE : function which updates deformation-related information at
% the material point level, i.e., vorticity, strain-rate, incremental and
% current deformation gradient and material point properties (dimensions)
%   Detailed explanation goes here 
%% PREPROCESSING
nX    = meD.DoF*no-1                                                      ;%
nY    = nX+1                                                              ;%
%--------------------------------------------------------------------------%
%% SPIN UPDATE
mpD.w(:,3)= 0.5*accumarray(mp,mpD.dSy.*meD.u(nX)-mpD.dSx.*meD.u(nY),[mpD.n 1]);% mp spin tensor xy
%--------------------------------------------------------------------------%

%% STRAIN CALCULATION
mpD.e(1,:)=accumarray(mp,mpD.dSx.*meD.u(nX)                   ,[mpD.n 1]) ;%
mpD.e(2,:)=accumarray(mp,mpD.dSy.*meD.u(nY)                   ,[mpD.n 1]) ;%
mpD.e(4,:)=accumarray(mp,mpD.dSy.*meD.u(nX)+mpD.dSx.*meD.u(nY),[mpD.n 1]) ;%
%--------------------------------------------------------------------------%

%% INCREMENTAL DEFORMATION GRADIENT UPDATE deformation gradient increment
mpD.dD(:,1) = 1.0+accumarray(mp,mpD.dSx.*meD.u(nX),[mpD.n 1])             ;% incremental deformation gradient tensor xx
mpD.dD(:,2) = 0.0+accumarray(mp,mpD.dSy.*meD.u(nX),[mpD.n 1])             ;% incremental deformation gradient tensor xy
mpD.dD(:,3) = 0.0+accumarray(mp,mpD.dSx.*meD.u(nY),[mpD.n 1])             ;% incremental deformation gradient tensor yx
mpD.dD(:,4) = 1.0+accumarray(mp,mpD.dSy.*meD.u(nY),[mpD.n 1])             ;% incremental deformation gradient tensor yy
%--------------------------------------------------------------------------%

%% DEFORMATION GRADIENT UPDATE
mpD.F(:,1)  = mpD.dD(:,1).*mpD.F(:,1)+mpD.dD(:,2).*mpD.F(:,3)             ;% deformation gradient xx
mpD.F(:,2)  = mpD.dD(:,1).*mpD.F(:,2)+mpD.dD(:,2).*mpD.F(:,4)             ;% deformation gradient xy
mpD.F(:,3)  = mpD.dD(:,3).*mpD.F(:,1)+mpD.dD(:,4).*mpD.F(:,3)             ;% deformation gradient yx
mpD.F(:,4)  = mpD.dD(:,3).*mpD.F(:,2)+mpD.dD(:,4).*mpD.F(:,4)             ;% deformation gradient yy
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
mpD.J = mpD.dD(:,1).*mpD.dD(:,4)-mpD.dD(:,3).*mpD.dD(:,2)                 ;% determinant of the deformation gradient tensor
V     = cross([mpD.r1 zeros(mpD.n,1)],[mpD.r2 zeros(mpD.n,1)])            ;% update mp volume
mpD.V = V(:,end)                                                          ;% save new volume
mpD.n0= 1-(1-mpD.n0)./mpD.J                                               ;% upadte porosity

end

