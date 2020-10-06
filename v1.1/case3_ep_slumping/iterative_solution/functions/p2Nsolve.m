function [meD] = p2Nsolve(meD,mpD,g,dt,l2g,c2N,bc)
%f_P2NINT function that integrate every contributing material points to
%their neighboring nodes, which result in nodal mass, momentum and forces.
%The latter are further used to explicitly calculate a solution of the
%momentum balance equation.
%   Detailed explanation goes here
%% INITIALIZATION
% NODAL VECTOR INITIALIZATION
meD.m(:) = 0.0 ; meD.mr(:) = 0.0 ; meD.f(:) = 0.0 ; meD.d(:) = 0.0        ;% mass | repmated mass | force balance | damping force
meD.a(:) = 0.0 ;  meD.p(:) = 0.0 ; meD.v(:) = 0.0 ; meD.u(:) = 0.0        ;% acceleration | momentum | velocity | incremental displacement
meD.fi(:)= 0.0 ;
%--------------------------------------------------------------------------%
%% CONTRIBUTION TO NODES
for mp=1:mpD.n % for any material point
    for N=1:meD.nNp % for any neighboring node of any material pointv
        l             = c2N(mp,N)                                         ;%
        meD.m(l)      = meD.m(l)+mpD.S(mp,N).*mpD.m(mp)                   ;%
        iDy           = meD.DoF*l                                         ;%
        meD.p(iDy-1)  = meD.p(iDy-1)+mpD.S(mp,N).*mpD.p(mp,1)             ;%
        meD.p(iDy  )  = meD.p(iDy  )+mpD.S(mp,N).*mpD.p(mp,2)             ;%
        meD.f(iDy  )  = meD.f(iDy  )-mpD.S(mp,N).*mpD.m(mp).*g            ;%
        meD.fi(iDy-1) = meD.fi(iDy-1)+mpD.V(mp).*(mpD.B(:,2*N-1,mp)'*mpD.s(:,mp))        ;%
        meD.fi(iDy  ) = meD.fi(iDy  )+mpD.V(mp).*(mpD.B(:,2*N  ,mp)'*mpD.s(:,mp))        ;%                 
    end
end
%% SOLVE EXPLICIT MOMENTUM BALANCE EQUATION
meD.f      = meD.f-meD.fi                                                 ;%
% UPDATE NODAL QUANTITITES
iDx        = 1:meD.DoF:meD.nDoF(2)-1                                      ;%
iDy        = iDx+1                                                        ;%
% COMPUTE NODAL FORCE
meD.d(iDx) = sqrt(meD.f(iDx).^2+meD.f(iDy).^2)                            ;%
meD.d(iDy) = meD.d(iDx)                                                   ;%
meD.f      = meD.f - meD.vd.*meD.d.*sign(meD.p)                           ;%
% UPDATE NODAL MOMENTUM
meD.p      = meD.p + dt.*meD.f                                            ;%
% COMPUTE NODAL ACCELERATION AND VELOCITY
meD.mr = reshape(repmat(meD.m',meD.DoF,1),meD.nDoF(2),1)                  ;%
meD.a      = meD.f./meD.mr                                                ;%
meD.v      = meD.p./meD.mr                                                ;%
meD.a(isnan(meD.a))=0.0                                                   ;%
meD.v(isnan(meD.v))=0.0                                                   ;%
% BOUNDARY CONDITIONS: FIX DIRICHLET BOUNDARY CONDITIONS
meD.a(bc.x(:,1))=bc.x(:,2)                                                ;%
meD.a(bc.y(:,1))=bc.y(:,2)                                                ;%
meD.a(bc.y(:,1)-1)=bc.y(:,2)                                                ;%
meD.v(bc.x(:,1))=bc.x(:,2)                                                ;%
meD.v(bc.y(:,1))=bc.y(:,2)                                                ;%
meD.v(bc.y(:,1)-1)=bc.y(:,2)                                                ;%
end