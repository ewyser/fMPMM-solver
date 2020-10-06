function [bc,BC,xB,yB] = setBCs(meD,du)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% BOUNDARY CONDITIONS
xB      = [0 meD.L(1)]                                                    ;%
yB      = [0 meD.L(2)]                                                    ;%                                                        ;%
%--------------------------------------------------------------------------% 
[row]   = find(meD.y<=yB(1))                                              ;%
BC      = row                                                             ;%
bc_num  = [row*meD.DoF]                                                     ;%
bc      = [bc_num,-du*ones(size(bc_num))]                                 ;%
%--------------------------------------------------------------------------% 
% [row]   = find(meD.y>=yB(2))                                              ;%
% BC      = [BC;row]                                                        ;%
% bc_num  = row*meD.DoF                                                     ;%
% bc      = [bc;[bc_num,du*ones(size(bc_num))]]                             ;%
%--------------------------------------------------------------------------% 
[row]   = find(meD.x<=xB(1))                                              ;%
BC      = [BC;row]                                                        ;%
bc_num  = row*meD.DoF-1                                                   ;%
bc      = [bc;bc_num,du*ones(size(bc_num))]                               ;%
%--------------------------------------------------------------------------% 
[row]   = find(meD.x>=xB(2))                                              ;%
BC      = [BC;row]                                                        ;%
bc_num  = row*meD.DoF-1                                                   ;%
bc      = [bc;[bc_num,-du*ones(size(bc_num))]]                            ;%
end

