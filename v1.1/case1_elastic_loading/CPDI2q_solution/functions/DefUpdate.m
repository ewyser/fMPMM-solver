function [mpD] = DefUpdate(meD,mpD,N,l2g,c2N)
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
mpD.U(:,1:meD.DoF:end) = meD.u(iDx)                                       ;% nodal x-displacement
mpD.U(:,2:meD.DoF:end) = meD.u(iDy)                                       ;% nodal y-displacement
mpD.e = permute(sum(mpD.B.*repmat(permute(mpD.U,[3 2 1]),nstr,1),2),[1 3 2]);% mp strain
%--------------------------------------------------------------------------%

% % %% INCREMENTAL DEFORMATION GRADIENT UPDATE deformation gradient increment
% % mpD.dF(:,1) = 1.0+sum(mpD.dSx.*meD.u(iDx),2)                              ;% incremental deformation gradient xx
% % mpD.dF(:,2) = 0.0+sum(mpD.dSy.*meD.u(iDx),2)                              ;% incremental deformation gradient xy
% % mpD.dF(:,3) = 0.0+sum(mpD.dSx.*meD.u(iDy),2)                              ;% incremental deformation gradient yx
% % mpD.dF(:,4) = 1.0+sum(mpD.dSy.*meD.u(iDy),2)                              ;% incremental deformation gradient yy
% % %--------------------------------------------------------------------------%
% % 
% % %% DEFORMATION GRADIENT UPDATE
% % mpD.F(:,1)  = mpD.dF(:,1).*mpD.F(:,1)+mpD.dF(:,2).*mpD.F(:,3)             ;% deformation gradient xx
% % mpD.F(:,2)  = mpD.dF(:,1).*mpD.F(:,2)+mpD.dF(:,2).*mpD.F(:,4)             ;% deformation gradient xy
% % mpD.F(:,3)  = mpD.dF(:,3).*mpD.F(:,1)+mpD.dF(:,4).*mpD.F(:,3)             ;% deformation gradient yx
% % mpD.F(:,4)  = mpD.dF(:,3).*mpD.F(:,2)+mpD.dF(:,4).*mpD.F(:,4)             ;% deformation gradient yy
% % %--------------------------------------------------------------------------%

%% UPDATE MATERIAL POINT DOMAIN AND COORDINATE
C=1:4;
% CORNERS COORDINATE UPDATE
for c=1:4
    mpD.xc(:,c) = mpD.xc(:,c) + sum(N(:,C).*meD.u(meD.DoF.*c2N(:,C)-1),2);
    mpD.yc(:,c) = mpD.yc(:,c) + sum(N(:,C).*meD.u(meD.DoF.*c2N(:,C)  ),2);
    C = C+4;
end
% MP'S COORDINATE UPDATE
mpD.x = [mean(mpD.xc,2) mean(mpD.yc,2)];
% MP'S VOLUME UPDATE
mpD.V = ((mpD.xc(:,1).*mpD.yc(:,2)-mpD.xc(:,2).*mpD.yc(:,1))+...
         (mpD.xc(:,2).*mpD.yc(:,3)-mpD.xc(:,3).*mpD.yc(:,2))+...
         (mpD.xc(:,3).*mpD.yc(:,4)-mpD.xc(:,4).*mpD.yc(:,3))+...
         (mpD.xc(:,4).*mpD.yc(:,1)-mpD.xc(:,1).*mpD.yc(:,4))).*0.5;
mpD.J = mpD.V./mpD.V0;
%--------------------------------------------------------------------------%

end

