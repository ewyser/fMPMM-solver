function [mpD] = SdS(mpD,meD,c2N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for mp = 1:mpD.n
    x1  =  mpD.r1(mp,2)-mpD.r2(mp,2)                                      ;% 
    x2  =  mpD.r1(mp,2)+mpD.r2(mp,2)                                      ;%
    y1  =  mpD.r2(mp,1)-mpD.r1(mp,1)                                      ;% 
    y2  = -mpD.r1(mp,1)-mpD.r2(mp,1)                                      ;%
    wxC = [x1 x2 -x1 -x2]                                                 ;%
    wxC = (1./(2*mpD.V(mp))).*wxC                                         ;%
    wyC = [y1 y2 -y1 -y2]                                                 ;%
    wyC = (1./(2*mpD.V(mp))).*wyC                                         ;%
    for c=1:4
        for nn = 1:4
            iD = (nn+c*4)-4;
            dX = mpD.xc(mp,c) - meD.x(c2N(mp,iD))                         ;% distance x from material point to node
            dY = mpD.yc(mp,c) - meD.y(c2N(mp,iD))                         ;% distance y from material point to node
            Nx = NdN(dX,meD.h(1))                                         ;% - see function
            Ny = NdN(dY,meD.h(2))                                         ;% - see function
            N  = Nx.* Ny                                                  ;% FEM basis function
            
            mpD.S(mp,iD)        = N.*0.25                                                    ;%
            mpD.dSx(mp,iD)      = N.*wxC(c)                                                     ;%
            mpD.dSy(mp,iD)      = N.*wyC(c)                                                     ;%
        end
    end
end
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
