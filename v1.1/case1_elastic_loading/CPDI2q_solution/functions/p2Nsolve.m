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
%--------------------------------------------------------------------------%

%% CONTRIBUTION TO NODES
% PREPROCESSING
m = reshape( mpD.S.*repmat(mpD.m,1,meD.nNp)      ,mpD.n*meD.nNp   ,1)     ;% preprocesing of mass global vector
p = reshape([mpD.S.*repmat(mpD.p(:,1),1,meD.nNp);...
             mpD.S.*repmat(mpD.p(:,2),1,meD.nNp)],mpD.n*meD.nDoF(1),1)    ;% preprocesing of momentum global vector
f = reshape([mpD.S.*0.0                         ;...
             mpD.S.*repmat(mpD.m,1,meD.nNp).*-g ],mpD.n*meD.nDoF(1),1)    ;% preprocesing of external force global vector  
fi= squeeze(sum(mpD.B.*repmat(reshape(mpD.s,size(mpD.s,1),1,mpD.n),1,meD.nDoF(1)),1)).*repmat(mpD.V',32,1); % 
% CONTRIBUTION FROM p TO N
meD.m = accumarray(c2N(:),m,[meD.nN     1])                               ;% mass global vector
meD.p = accumarray(l2g(:),p,[meD.nDoF(2) 1])                              ;% momentum global vector
meD.f = accumarray(l2g(:),f,[meD.nDoF(2) 1])                              ;% external force global vector
for n = 1:meD.nNp                                                          % BEGIN ITERATION OVER meD.nNe NEIGHBORING NODES FOR ALL MP                                                               
    l = [(meD.DoF*c2N(:,n)-1);(meD.DoF*c2N(:,n))]                         ;% local to global index
    meD.f = meD.f - accumarray(l,[fi(n*meD.DoF-1,:)';...
                                  fi(n*meD.DoF  ,:)'],[meD.nDoF(2) 1])    ;% external force global vector - nodal internal force global vector
end                                                                        % END ITERATION OVER meD.nNe NEIGHBORING NODES FOR ALL MP
%--------------------------------------------------------------------------%

%% SOLVE EXPLICIT MOMENTUM BALANCE EQUATION
iDx        = 1:meD.DoF:meD.nDoF(2)-1                                      ;% x component of global vector node
iDy        = iDx+1                                                        ;% y component of global vector node
% COMPUTE GLOBAL FORCE VECTOR                                                 
meD.d(iDx) = sqrt(meD.f(iDx).^2+meD.f(iDy).^2)                            ;% x component damping force vector
meD.d(iDy) = meD.d(iDx)                                                   ;% y component damping force vector
meD.f      = meD.f - meD.vd*meD.d.*sign(meD.p)                            ;% force vector - damping force vector
% COMPUTE GLOBAL ACCELERATION VECTOR
meD.mr     = reshape(repmat(meD.m',meD.DoF,1),meD.nDoF(2),1)              ;% repmat mass vector
meD.a      = meD.f./meD.mr                                                ;% compute acceleration vector
iD         = meD.mr==0                                                    ;% null nodal mass index
meD.a(iD)  = 0.0                                                          ;% zeroed acceleration vector
%--------------------------------------------------------------------------%

%% BOUNDARY CONDITIONS: FIX DIRICHLET BOUNDARY CONDITIONS
meD.a(bc.x(:,1))=bc.x(:,2)                                                ;%
meD.a(bc.y(:,1))=bc.y(:,2)                                                ;%
%-------------------------------------------------------------------------%

end