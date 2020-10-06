function [meD,mpD] = mapN2p(meD,mpD,dt,l2g,c2N,bc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
for mp = 1:mpD.n
    iDx = meD.DoF*c2N(mp,:)-1                                             ;% x-component index
    iDy = iDx+1                                                           ;% y-component index
    mpD.v(mp,:) = mpD.v(mp,:)+dt*[mpD.S(mp,:)*meD.a(iDx)...
                                  mpD.S(mp,:)*meD.a(iDy)]                 ;% update mp velocity                         
    % MP'S MOMENTUM UPDATE
    mpD.p(mp,:) = mpD.v(mp,:).*repmat(mpD.m(mp),1,meD.DoF)                ;% update mp momentum
    % MP'S COORDINATE UPDATE
    mpD.x(mp,:) = mpD.x(mp,:)+dt*[mpD.S(mp,:)*meD.v(iDx)...
                                  mpD.S(mp,:)*meD.v(iDy)]                 ;% update mp position
end
%% UPDATE NODAL MOMENTUM WITH UPDATED MP MOMENTUM
meD.p(:) = 0.0                                                            ;% initialized momentum
for mp=1:mpD.n % for any material point
    for N=1:meD.nNp % for any neighboring node of any material point
        iDy          = meD.DoF*c2N(mp,N)                                  ;%
        meD.p(iDy-1) = meD.p(iDy-1)+mpD.S(mp,N).*mpD.p(mp,1)              ;%
        meD.p(iDy  ) = meD.p(iDy  )+mpD.S(mp,N).*mpD.p(mp,2)              ;%
    end
end
meD.u = dt*meD.p./meD.mr                                                  ;% incremental displacement 
meD.u(isnan(meD.u))=0.0                                                   ;% -
% BOUNDARY CONDITIONS: FIX DIRICHLET BOUNDARY CONDITIONS
meD.u(bc.x(:,1))=bc.x(:,2)                                                ;%
meD.u(bc.y(:,1))=bc.y(:,2)                                                ;% clear temporary variables
meD.u(bc.y(:,1)-1)=bc.y(:,2)                                                ;%
%% DISPLACEMENT UPDATE                                                                 
for mp = 1:mpD.n
    iDx = meD.DoF*c2N(mp,:)-1                                             ;% x-component index
    iDy = iDx+1                                                           ;% y-component index
    % MP'S DISPLACEMENT UPDATE
    mpD.u(mp,:) = mpD.u(mp,:)+[mpD.S(mp,:)*meD.u(iDx) mpD.S(mp,:)*meD.u(iDy)]      ;% update mp velocity
end
%--------------------------------------------------------------------------%

end

