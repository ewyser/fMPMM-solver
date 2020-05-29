function [meD,mpD] = mapN2p(meD,mpD,dt,mp,no,BC)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% PREPROCESSING
no    = meD.DoF*no-1                                                      ;%
%% INTERPOLATE SOLUTIONS N to p
% VELOCITY UPDATE                                                          
mpD.v(:,1) =  mpD.PF(1)   .*(mpD.v(:,1)+dt.*accumarray(mp,mpD.S.*meD.a(no),[mpD.n 1]))+...
    mpD.PF(2).*(accumarray(mp,mpD.S.*meD.v(no),[mpD.n 1]))                      ;% PIC-FLIP update
mpD.v(:,2) = mpD.PF(1)   .*(mpD.v(:,2)+dt.*accumarray(mp,mpD.S.*meD.a(no+1),[mpD.n 1]))+...
    mpD.PF(2).*(accumarray(mp,mpD.S.*meD.v(no+1),[mpD.n 1]))                    ;% PIC-FLIP update
%------------------------------------------------------------------%

% MP'S MOMENTUM UPDATE
mpD.p = mpD.v.*repmat(mpD.m,1,meD.DoF);
%------------------------------------------------------------------%

% MP'S COORDINATE UPDATE
mpD.x = mpD.x+dt*[accumarray(mp,mpD.S.*meD.v(no  ),[mpD.n 1]),...
                  accumarray(mp,mpD.S.*meD.v(no+1),[mpD.n 1])];
%% UPDATE NODAL MOMENTUM WITH UPDATED MP MOMENTUM (MUSL OR DOUBLE MAPPING)
meD.p(:) = 0.0                                                            ;% initialized momentum
meD.p    = accumarray([no;no+1],[mpD.S.*mpD.p(mp,1);...
                                 mpD.S.*mpD.p(mp,2)],[meD.nDF(2) 1])      ;% updated nodal momentum
meD.u    = dt*(meD.p./meD.mr)                                             ;% incremental displacement 
iD       = meD.mr==0                                                      ;% index of zero nodal mass
meD.u(iD)= 0.0                                                            ;% regularize incremental displacement when nodal mass is null
%--------------------------------------------------------------------------%

%% BOUNDARY CONDITIONS: FIX DIRICHLET BOUNDARY CONDITIONS
BCx = [BC.xi        ;BC.xs      ]                                         ;% boundary nodes index
BCy = [BC.yi                    ]                                         ;% boundary nodes index
meD.u([meD.DoF*BCx-1;meD.DoF*BCx ;meD.DoF*BCy])=0.0                       ;% fix BC's for displacement
meD.u( meD.DoF*(BCy+1)           )=0.0                                    ;% fix BC's for displacement
clear iD BCx BCy                                                          ;% clear temporary variables
%--------------------------------------------------------------------------%

end

