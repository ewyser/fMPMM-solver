function [mpD] = SdS_loop(meD,mpD,p2n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% PREPROCESSING
for mp=1:mpD.n
    %% COMPUTE (X,Y)-BASIS FUNCTION
    for nn=1:meD.nNe
        dX       = mpD.x(mp,1) - meD.x(p2n(mp,nn))                        ;% distance x from material point to node
        dY       = mpD.x(mp,2) - meD.y(p2n(mp,nn))                        ;% distance y from material point to node
        [Sx,dSx] = NdN(dX,meD.h(1),mpD.l(mp,1))                           ;%
        [Sy,dSy] = NdN(dY,meD.h(2),mpD.l(mp,2))                           ;%
        %% CONVOLUTION OF BASIS FUNCTIONS
        mpD.S(mp,nn)    =  Sx.* Sy                                        ;% basis function
        mpD.dSx(mp,nn)  = dSx.* Sy                                        ;% x-gradient basis function
        mpD.dSy(mp,nn)  =  Sx.*dSy                                        ;% y-gradient basis function
    end
    %% B MATRIX INITIALIZATION
    mpD.B(:,:,mp) = 0.0;
    iDx       = 1:meD.DoF:meD.nDoF(1)                                     ;% x indexing
    iDy       = iDx+1                                                     ;% y indexing
    % B MATRIX ASSEMBLY
    mpD.B(1,iDx,mp)= mpD.dSx(mp,:)'                                       ;% strain-displacement matrix
    mpD.B(2,iDy,mp)= mpD.dSy(mp,:)'                                       ;% strain-displacement matrix
    mpD.B(3,iDx,mp)= mpD.dSy(mp,:)'                                       ;% strain-displacement matrix
    mpD.B(3,iDy,mp)= mpD.dSx(mp,:)'                                       ;% strain-displacement matrix
end
end
function [N,dN]=NdN(dX,h,lp)
%% COMPUTE BASIS FUNCTIONS (see Steffen etal, 2008)
lp  = 2.*lp;
if(abs(dX)<0.5*lp)
    N = 1-((4.*dX.^2+lp.^2)./(4.*h.*lp))                                  ;%
    dN= -((8.*dX)./(4.*h.*lp))                                            ;%
elseif(abs(dX)>=0.5*lp & abs(dX)<=(h-0.5*lp))
    N = 1-abs(dX)./h                                                      ;%
    dN= sign(dX).*(-1./h)                                                 ;%
elseif(abs(dX)>=(h-0.5*lp) & abs(dX)<(h+0.5*lp))
    N = ((h+0.5.*lp-abs(dX)).^2)./(2.*h.*lp)                              ;%
    dN= (-sign(dX).*(h+0.5.*lp-abs(dX)))./(h.*lp)                         ;%
else
    N = 0.0;
    dN= 0.0;
end
end

%% REFERENCE
% Steffen, M., Wallstedt, P.C., Guilkey, J.E., Berzins, M. (2008)
% Examination and Analysis of Implementation Choices within the Material
% Point Method (MPM). CMES. 31. 107-127.
%-------------------------------------------------------------------------%
% Coombs, W.M., Augarde, C.E. (2020) AMPLE: A Material Point Learning
% Environment. Adv. Eng. Soft. 139. 102748.
