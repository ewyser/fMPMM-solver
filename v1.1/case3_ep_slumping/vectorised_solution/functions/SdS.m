function [mpD] = SdS(mpD,meD,c2N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% COMPUTE (X,Y)-BASIS FUNCTION
xC   = [repmat(mpD.xc(:,1),1,4),repmat(mpD.xc(:,2),1,4),...
        repmat(mpD.xc(:,3),1,4),repmat(mpD.xc(:,4),1,4)]                  ;% array of x coordinate of corners
yC   = [repmat(mpD.yc(:,1),1,4),repmat(mpD.yc(:,2),1,4),...
        repmat(mpD.yc(:,3),1,4),repmat(mpD.yc(:,4),1,4)]                  ;% array of y coordinate of corners
D    = xC-meD.x(c2N)                                                      ;% x distance from corner to node
[Nx] = NdN(D,meD.h(1))                                                    ;% - see function
D    = yC-meD.y(c2N)                                                      ;% y distance from corner to node
[Ny] = NdN(D,meD.h(2))                                                    ;% - see function
%% CONVOLUTION OF TEMPORARY BASIS FUNCTIONS
N    =  Nx.* Ny                                                           ;% FEM basis function
%% CPDI BASIS FUNCTION RECONSTRUCTION
x1   =  mpD.r1(:,2)-mpD.r2(:,2); x2 =  mpD.r1(:,2)+mpD.r2(:,2)            ;%
y1   =  mpD.r2(:,1)-mpD.r1(:,1); y2 = -mpD.r1(:,1)-mpD.r2(:,1)            ;%
wxC  = [repmat(x1,1,4) repmat(x2,1,4) -repmat(x1,1,4) -repmat(x2,1,4)]    ;%
wxC  =  repmat(1./(2*mpD.V),1,meD.nNp).*wxC                               ;%
wyC  = [repmat(y1,1,4) repmat(y2,1,4) -repmat(y1,1,4) -repmat(y2,1,4)]    ;%
wyC  =  repmat(1./(2*mpD.V),1,meD.nNp).*wyC                               ;%

mpD.S        = N.*0.25                                                    ;%
mpD.dSx      = N.*wxC                                                     ;%
mpD.dSy      = N.*wyC                                                     ;%

%% B MATRIX ASSEMBLY
iDx            = 1:meD.DoF:meD.nDoF(1)-1                                  ;% x component global node index
iDy            = iDx+1                                                    ;% y component global node index
mpD.B(1,iDx,:) = mpD.dSx'                                                 ;% -
mpD.B(2,iDy,:) = mpD.dSy'                                                 ;% -
mpD.B(3,iDx,:) = mpD.dSy'                                                 ;% -
mpD.B(3,iDy,:) = mpD.dSx'                                                 ;% -
%--------------------------------------------------------------------------%

end
function [N]=NdN(dXi,hi)
%% COMPUTE BASIS FUNCTIONS
c1 = (abs(dXi)<=hi)                                                       ;%
% BASIS FUNCTION
N=c1.*(1-abs(dXi)./hi)                                                    ;%
end
