function [mpD] = DefUpdate(meD,mpD,N,mp,no,c2N)
% f_DEFUPDATE : function which updates deformation-related information at
% the material point level, i.e., vorticity, strain-rate, incremental and
% current deformation gradient and material point properties (dimensions)
%   Detailed explanation goes here 
%% PREPROCESSING
iDx    = meD.DoF*no-1                                                     ;%
iDy    = iDx+1                                                            ;%
%--------------------------------------------------------------------------%
%% SPIN UPDATE
mpD.w(:,3)= 0.5*accumarray(mp,mpD.dSy.*meD.u(iDx)-mpD.dSx.*meD.u(iDy),[mpD.n 1]);% mp spin tensor xy
%--------------------------------------------------------------------------%

%% STRAIN CALCULATION
mpD.e(1,:)=accumarray(mp,mpD.dSx.*meD.u(iDx)                    ,[mpD.n 1]);%
mpD.e(2,:)=accumarray(mp,mpD.dSy.*meD.u(iDy)                    ,[mpD.n 1]);%
mpD.e(4,:)=accumarray(mp,mpD.dSy.*meD.u(iDx)+mpD.dSx.*meD.u(iDy),[mpD.n 1]);%
%--------------------------------------------------------------------------%

%% INCREMENTAL DEFORMATION GRADIENT UPDATE deformation gradient increment
mpD.dF(:,1) = 1.0+accumarray(mp,mpD.dSx.*meD.u(iDx),[mpD.n 1])            ;% incremental deformation gradient xx
mpD.dF(:,2) = 0.0+accumarray(mp,mpD.dSx.*meD.u(iDy),[mpD.n 1])            ;% incremental deformation gradient xy
mpD.dF(:,3) = 0.0+accumarray(mp,mpD.dSy.*meD.u(iDx),[mpD.n 1])            ;% incremental deformation gradient yx
mpD.dF(:,4) = 1.0+accumarray(mp,mpD.dSy.*meD.u(iDy),[mpD.n 1])            ;% incremental deformation gradient yy
%--------------------------------------------------------------------------%

%% DEFORMATION GRADIENT UPDATE
mpD.F(:,1)  = mpD.dF(:,1).*mpD.F(:,1)+mpD.dF(:,2).*mpD.F(:,3)             ;% deformation gradient xx
mpD.F(:,2)  = mpD.dF(:,1).*mpD.F(:,2)+mpD.dF(:,2).*mpD.F(:,4)             ;% deformation gradient xy
mpD.F(:,3)  = mpD.dF(:,3).*mpD.F(:,1)+mpD.dF(:,4).*mpD.F(:,3)             ;% deformation gradient yx
mpD.F(:,4)  = mpD.dF(:,3).*mpD.F(:,2)+mpD.dF(:,4).*mpD.F(:,4)             ;% deformation gradient yy
mpD.J       = mpD.F(:,1) .*mpD.F(:,4)-mpD.F(:,3) .*mpD.F(:,2)             ;% determinant of the deformation gradient increment
%--------------------------------------------------------------------------%

%% UPDATE MATERIAL POINT DOMAIN AND COORDINATE
C=[1:4;5:8;9:12;13:16];
% CORNERS COORDINATE UPDATE
for c=1:4
    mpD.xc(:,c) = mpD.xc(:,c) + sum(N(:,C(c,:)).*meD.u(meD.DoF.*c2N(:,C(c,:))-1),2);
    mpD.yc(:,c) = mpD.yc(:,c) + sum(N(:,C(c,:)).*meD.u(meD.DoF.*c2N(:,C(c,:))  ),2);
end
% MP'S COORDINATE UPDATE
mpD.x = [mean(mpD.xc,2) mean(mpD.yc,2)];
% MP'S VOLUME UPDATE
mpD.V = ((mpD.xc(:,1).*mpD.yc(:,2)-mpD.xc(:,2).*mpD.yc(:,1))+...
         (mpD.xc(:,2).*mpD.yc(:,3)-mpD.xc(:,3).*mpD.yc(:,2))+...
         (mpD.xc(:,3).*mpD.yc(:,4)-mpD.xc(:,4).*mpD.yc(:,3))+...
         (mpD.xc(:,4).*mpD.yc(:,1)-mpD.xc(:,1).*mpD.yc(:,4))).*0.5;
%--------------------------------------------------------------------------%

end

