function [meD] = p2Nsolve(meD,mpD,g,dt,mp,no,BC)
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
meD.m = accumarray(no,mpD.m(mp).*mpD.S,[meD.nN 1]);
no    = meD.DoF*no-1;
meD.p = accumarray([no;no+1],[mpD.S.*mpD.p(mp,1);mpD.S.*mpD.p(mp,2)],[meD.nDF(2) 1]);
meD.f = accumarray([   no+1], mpD.S.*mpD.m(mp).*(-g),[meD.nDF(2) 1]);
meD.f = meD.f-accumarray([no;no+1],mpD.V([mp;mp]).*[(mpD.dSx.*mpD.s(1,mp)'+mpD.dSy.*mpD.s(4,mp)');...
    (mpD.dSx.*mpD.s(4,mp)'+mpD.dSy.*mpD.s(2,mp)')],[meD.nDF(2) 1]) ;
%--------------------------------------------------------------------------%

%% SOLVE EXPLICIT MOMENTUM BALANCE EQUATION
% UPDATE GLOBAL NODAL INFORMATIONS
iDx        = 1:meD.DoF:meD.nDF(2)-1                                       ;% x component of global vector node
iDy        = iDx+1                                                        ;% y component of global vector node
% COMPUTE GLOBAL NODAL FORCE                                                  
meD.d(iDx) = sqrt(meD.f(iDx).^2+meD.f(iDy).^2)                            ;% x component damping force
meD.d(iDy) = meD.d(iDx)                                                   ;% y component damping force
meD.f      = meD.f - meD.vd*meD.d.*sign(meD.p)                            ;% nodal force - damping force
% UPDATE GLOBAL NODAL MOMENTUM
meD.p      = meD.p + dt*meD.f                                             ;% forward euler momentum
% COMPUTE GLOBAL NODAL ACCELERATION AND VELOCITY
meD.mr     = reshape(repmat(meD.m',meD.DoF,1),meD.nDF(2),1)               ;% repmat nodal mass
iD         = meD.mr==0                                                    ;% null nodal mass index
meD.a      = meD.f./meD.mr                                                ;% compute nodal acceleration
meD.v      = meD.p./meD.mr                                                ;% compute nodal velocity
meD.a(iD)  = 0.0                                                          ;% zeroed nodal acceleration
meD.v(iD)  = 0.0                                                          ;% zeroed nodal velocity

% BOUNDARY CONDITIONS: FIX DIRICHLET BOUNDARY CONDITIONS
BCx =  [BC.xi;BC.xs]                                                      ;% boundary nodes index
BCy =  [BC.yi      ]                                                      ;% boundary nodes index
meD.a([meD.DoF*BCx-1;meD.DoF*BCy])=0.0                                    ;% fix BC's
meD.v([meD.DoF*BCx-1;meD.DoF*BCy])=0.0                                    ;% fix BC's
meD.a(meD.DoF*(BCy+1))=0.0                                                ;% fix BC's
meD.v(meD.DoF*(BCy+1))=0.0                                                ;% fix BC's
clear m p f fi iDx iDy iD BCx BCy                                         ;% clear temporary variables
%-------------------------------------------------------------------------%

end