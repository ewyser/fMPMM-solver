function [meD,mpD] = mapN2p(meD,mpD,dt,l2g,c2N,bc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% INTERPOLATE SOLUTIONS N to p
iDx        = meD.DoF*c2N-1                                                ;% x-component index
iDy        = iDx+1                                                        ;% y-component index
% VELOCITY UPDATE                                                                 
mpD.v      = mpD.v+dt*[sum(mpD.S.*meD.a(iDx),2) sum(mpD.S.*meD.a(iDy),2)] ;% update mp velocity
% MOMENTUM UPDATE
mpD.p      = mpD.v.*repmat(mpD.m,1,meD.DoF)                               ;% update mp momentum
% MP'S COORDINATE UPDATE
mpD.x      = mpD.x+dt*[sum(mpD.S.*meD.v(iDx),2) sum(mpD.S.*meD.v(iDy),2)] ;% update mp velocity
%--------------------------------------------------------------------------%
%% UPDATE NODAL MOMENTUM WITH UPDATED MP MOMENTUM (MUSL OR DOUBLE MAPPING)
meD.p(:) = 0.0                                                            ;% initialized momentum
p        = reshape([mpD.S.*repmat(mpD.p(:,1),1,meD.nNp) ;...
                    mpD.S.*repmat(mpD.p(:,2),1,meD.nNp)],...
                    mpD.n*meD.nDoF(1),1                  )                ;% preprocesing of nodal momentum 
meD.p    = accumarray(l2g(:),p,[meD.nDoF(2) 1])                           ;% nodal momentum 
meD.u    = dt*(meD.p./meD.mr)                                             ;% incremental displacement 
iD       = meD.mr==0                                                      ;% index of zero nodal mass
meD.u(iD)= 0.0                                                            ;% regularize incremental displacement when nodal mass is null
%--------------------------------------------------------------------------%

%% BOUNDARY CONDITIONS: FIX DIRICHLET BOUNDARY CONDITIONS
meD.u(bc.x(:,1))=bc.x(:,2)                                                ;%
meD.u(bc.y(:,1))=bc.y(:,2)                                                ;%
%--------------------------------------------------------------------------%

%% DISPLACEMENT UPDATE                                                                 
mpD.u    = mpD.u+[sum(mpD.S.*meD.u(iDx),2) sum(mpD.S.*meD.u(iDy),2)]      ;% update mp velocity
%--------------------------------------------------------------------------%

end

