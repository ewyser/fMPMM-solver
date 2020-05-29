function [meD,mpD] = mapN2p(meD,mpD,dt,mp,no,bc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% PREPROCESSING
iDx    = meD.DoF*no-1                                                     ;% x-component index
iDy    = iDx+1                                                            ;% y-component index
%% INTERPOLATE SOLUTIONS N to p
% VELOCITY UPDATE
mpD.v(:,1) = mpD.v(:,1)+dt.*accumarray(mp,mpD.S.*meD.a(iDx),[mpD.n 1])    ;% velocity update
mpD.v(:,2) = mpD.v(:,2)+dt.*accumarray(mp,mpD.S.*meD.a(iDy),[mpD.n 1])    ;% velocity update
% MOMENTUM UPDATE
mpD.p = mpD.v.*repmat(mpD.m,1,meD.DoF)                                    ;% momentum update
%--------------------------------------------------------------------------%

%% UPDATE NODAL MOMENTUM WITH UPDATED MP MOMENTUM (MUSL OR DOUBLE MAPPING)
meD.p(:) = 0.0                                                            ;% initialized momentum
meD.p    = accumarray([iDx;iDy],[mpD.S.*mpD.p(mp,1);...
                                 mpD.S.*mpD.p(mp,2)],[meD.nDF(2) 1])      ;% updated nodal momentum
meD.u    = dt*(meD.p./meD.mr)                                             ;% incremental displacement 
iD       = meD.mr==0                                                      ;% index of zero nodal mass
meD.u(iD)= 0.0                                                            ;% regularize incremental displacement when nodal mass is null
%--------------------------------------------------------------------------%

%% BOUNDARY CONDITIONS: FIX DIRICHLET BOUNDARY CONDITIONS
meD.u(bc.x(:,1))=meD.u(bc.x(:,1)).*bc.x(:,2)                              ;% 
meD.u(bc.y(:,1))=meD.u(bc.y(:,1)).*bc.y(:,2)                              ;% 
%--------------------------------------------------------------------------%

end

