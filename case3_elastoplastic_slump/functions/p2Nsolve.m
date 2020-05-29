function [meD] = p2Nsolve(meD,mpD,g,dt,mp,no,bc)
%f_P2NINT function that integrate every contributing material points to
%their neighboring nodes, which result in nodal mass, momentum and forces.
%The latter are further used to explicitly calculate a solution of the
%momentum balance equation.
%   Detailed explanation goes here
%% INITIALIZATION
% NODAL VECTOR INITIALIZATION
meD.m(:) = 0.0 ; meD.mr(:) = 0.0 ; meD.f(:) = 0.0 ; meD.d(:) = 0.0        ;% mass | repmated mass | force balance | damping force
meD.a(:) = 0.0 ;  meD.p(:) = 0.0 ; meD.v(:) = 0.0 ; meD.u(:) = 0.0        ;% acceleration | momentum | velocity | incremental displacement
%--------------------------------------------------------------------------%

%% CONTRIBUTION TO NODES
iDx   = meD.DoF*no-1                                                      ;% 
iDy   = iDx+1                                                             ;%

meD.m = accumarray(no,mpD.m(mp).*mpD.S,[meD.nN 1])                        ;% mass vector
meD.p = accumarray([iDx;iDy],[mpD.S.*mpD.p(mp,1);mpD.S.*mpD.p(mp,2)],[meD.nDF(2) 1]);% momentum vector
meD.f = accumarray(     iDy , mpD.S.*mpD.m(mp).*(-g),[meD.nDF(2) 1])      ;% external force vector
meD.f = meD.f-accumarray([iDx;iDy],mpD.V([mp;mp]).*[(mpD.dSx.*mpD.s(1,mp)'+mpD.dSy.*mpD.s(4,mp)');...
        (mpD.dSx.*mpD.s(4,mp)'+mpD.dSy.*mpD.s(2,mp)')],[meD.nDF(2) 1])    ;% out-of-balance force vector
%--------------------------------------------------------------------------%

%% SOLVE EXPLICIT MOMENTUM BALANCE EQUATION
% COMPUTE GLOBAL FORCE VECTOR                                                 
meD.d(iDx) = sqrt(meD.f(iDx).^2+meD.f(iDy).^2)                            ;% x component damping force vector
meD.d(iDy) = meD.d(iDx)                                                   ;% y component damping force vector
meD.f      = meD.f - meD.vd*meD.d.*sign(meD.p)                            ;% force vector - damping force vector
% UPDATE GLOBAL MOMENTUM VECTOR
meD.p      = meD.p + dt*meD.f                                             ;% forward euler momentum vector
% COMPUTE GLOBAL ACCELERATION AND VELOCITY VECTORS
meD.mr     = reshape(repmat(meD.m',meD.DoF,1),meD.nDF(2),1)               ;% repmat mass vector
meD.a      = meD.f./meD.mr                                                ;% compute acceleration vector
meD.v      = meD.p./meD.mr                                                ;% compute velocity vector
iD         = meD.mr==0                                                    ;% null nodal mass index
meD.a(iD)  = 0.0                                                          ;% zeroed acceleration vector
meD.v(iD)  = 0.0                                                          ;% zeroed  velocity vector
%--------------------------------------------------------------------------%

%% BOUNDARY CONDITIONS: FIX DIRICHLET BOUNDARY CONDITIONS
meD.a(bc.x(:,1))=meD.a(bc.x(:,1)).*bc.x(:,2)                              ;%
meD.a(bc.y(:,1))=meD.a(bc.y(:,1)).*bc.y(:,2)                              ;%
meD.v(bc.x(:,1))=meD.v(bc.x(:,1)).*bc.x(:,2)                              ;%
meD.v(bc.y(:,1))=meD.v(bc.y(:,1)).*bc.y(:,2)                              ;%
%-------------------------------------------------------------------------%

end