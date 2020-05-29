function [meD,mpD] = mapN2p(meD,mpD,dt,l2g,p2N,BC)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% INTERPOLATE SOLUTIONS N to p
iDx        = meD.DoF*p2N-1                                                ;% x-component index
iDy        = iDx+1                                                        ;% y-component index
% VELOCITY UPDATE                                                          
mpD.v(:,1) = mpD.PF(1)*(mpD.v(:,1)+dt*sum(mpD.N.*meD.a(iDx),2))+...
             mpD.PF(2)*(sum(mpD.N.*meD.v(iDx),2))                         ;% PIC-FLIP update
mpD.v(:,2) = mpD.PF(1)*(mpD.v(:,2)+dt*sum(mpD.N.*meD.a(iDy),2))+...
             mpD.PF(2)*(sum(mpD.N.*meD.v(iDy),2))                         ;% PIC-FLIP update
% MOMENTUM UPDATE
mpD.p      = mpD.v.*repmat(mpD.m,1,meD.DoF)                               ;% update mp momentum
% COORDINATE UPDATE
mpD.x      = mpD.x+dt*[sum(mpD.N.*meD.v(iDx),2) sum(mpD.N.*meD.v(iDy),2)] ;% update mp position
%--------------------------------------------------------------------------%

%% UPDATE NODAL MOMENTUM WITH UPDATED MP MOMENTUM (MUSL OR DOUBLE MAPPING)
meD.p(:) = 0.0                                                            ;% initialized momentum
P        = reshape([mpD.N.*repmat(mpD.p(:,1),1,meD.nNe) ;...
                    mpD.N.*repmat(mpD.p(:,2),1,meD.nNe)],...
                    mpD.n*meD.nDF(1),1                  )                 ;% preprocesing of nodal momentum 
meD.p    = accumarray(l2g(:),P,[meD.nDF(2) 1])                            ;% nodal momentum 
meD.u    = dt*(meD.p./meD.mr)                                             ;% incremental displacement 
iD       = meD.mr==0                                                      ;% index of zero nodal mass
meD.u(iD)= 0.0                                                            ;% regularize incremental displacement when nodal mass is null
%--------------------------------------------------------------------------%

%% BOUNDARY CONDITIONS: FIX DIRICHLET BOUNDARY CONDITIONS
BCx = [BC.xi        ;BC.xs      ]                                         ;% boundary nodes index
BCy = [BC.yi                    ]                                         ;% boundary nodes index
meD.u([meD.DoF*BCx-1;meD.DoF*BCy])=0.0                                    ;% fix BC's for displacement
meD.u( meD.DoF*(BCy+1)           )=0.0                                    ;% fix BC's for displacement
clear iDx iDy iD BCx BCy                                                  ;% clear temporary variables
%--------------------------------------------------------------------------%

end

