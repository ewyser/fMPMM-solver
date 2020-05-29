function [mpD,N,p,n] = SdS(mpD,meD,c2N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% COMPUTE (X,Y)-BASIS FUNCTION
xC   = [repmat(mpD.xc(:,1),1,4),repmat(mpD.xc(:,2),1,4),...
        repmat(mpD.xc(:,3),1,4) repmat(mpD.xc(:,4),1,4)]                  ;% array of x coordinate of corners
yC   = [repmat(mpD.yc(:,1),1,4) repmat(mpD.yc(:,2),1,4),...
        repmat(mpD.yc(:,3),1,4) repmat(mpD.yc(:,4),1,4)]                  ;% array of y coordinate of corners
D    = xC-meD.x(c2N)                                                      ;% x distance from corner to node
[Nx] = NdN(D,meD.h(1))                                                    ;% - see function
D    = yC-meD.y(c2N)                                                      ;% y distance from corner to node
[Ny] = NdN(D,meD.h(2))                                                    ;% - see function
%--------------------------------------------------------------------------%

%% CONVOLUTION OF TEMPORARY BASIS FUNCTIONS
N    =  Nx.* Ny                                                           ;% FE basis function
%--------------------------------------------------------------------------%

%% CPDI2 BASIS FUNCTION RECONSTRUCTION
[wC,wxC,wyC] = weightCPDI2q(mpD.xc,mpD.yc,mpD.V,meD.nNp)                  ;%
i            = repmat((1:mpD.n)',1,meD.nNp)                               ;%
j            = c2N                                                        ;%
[p,n,mpD.S ] = find(sparse(i,j,N.*wC ,mpD.n,meD.nN))                      ;%
[~,~,mpD.dSx]= find(sparse(i,j,N.*wxC,mpD.n,meD.nN))                      ;%
[~,~,mpD.dSy]= find(sparse(i,j,N.*wyC,mpD.n,meD.nN))                      ;%
%--------------------------------------------------------------------------%

end
function [N]=NdN(dXi,hi)
%% COMPUTE BASIS FUNCTIONS
c1 = abs(dXi)<=hi                                                         ;%
% BASIS FUNCTION
N=c1.*(1-abs(dXi)./hi)                                                    ;%
end
function [wC,wxC,wyC]=weightCPDI2q(xc,yc,Vp,nNp)
%% COMPUTE WEIGHT AND GRADIENT WEIGHT
% COMPUTE WEIGHT
a  = ((xc(:,4)-xc(:,1)).*(yc(:,2)-yc(:,3)))-((xc(:,2)-xc(:,3)).*(yc(:,4)-yc(:,1)));%
b  = ((xc(:,3)-xc(:,4)).*(yc(:,1)-yc(:,2)))-((xc(:,1)-xc(:,2)).*(yc(:,3)-yc(:,4)));%
c1 = 6*Vp-a-b                                                             ;%
c2 = 6*Vp-a+b                                                             ;%
c3 = 6*Vp+a+b                                                             ;%
c4 = 6*Vp+a-b                                                             ;%
wC = [repmat(c1,1,4) repmat(c2,1,4) repmat(c3,1,4) repmat(c4,1,4)]        ;%
wC =  repmat(1./(24*Vp),1,nNp).*wC                                        ;%

% COMPUTE GRADIENT WEIGHT
dx1 = yc(:,2)-yc(:,4)                                                     ;%
dx2 = yc(:,3)-yc(:,1)                                                     ;%
dx3 = yc(:,4)-yc(:,2)                                                     ;%
dx4 = yc(:,1)-yc(:,3)                                                     ;%
dy1 = xc(:,4)-xc(:,2)                                                     ;%
dy2 = xc(:,1)-xc(:,3)                                                     ;%
dy3 = xc(:,2)-xc(:,4)                                                     ;%
dy4 = xc(:,3)-xc(:,1)                                                     ;%
wxC = [repmat(dx1,1,4) repmat(dx2,1,4) repmat(dx3,1,4) repmat(dx4,1,4)]   ;%
wxC =  repmat(1./(2*Vp),1,nNp).*wxC                                       ;%
wyC = [repmat(dy1,1,4) repmat(dy2,1,4) repmat(dy3,1,4) repmat(dy4,1,4)]   ;%
wyC =  repmat(1./(2*Vp),1,nNp).*wyC                                       ;%
end