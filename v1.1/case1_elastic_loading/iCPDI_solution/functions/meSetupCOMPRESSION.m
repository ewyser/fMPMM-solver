function [meD] = meSetupCOMPRESSION(nEy,typeD)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% MESH INITIALIZATION
Ly      = 10                                                              ;% number of element in x direction
nEx     = 1                                                             ;% number of element in y direction
meD.L   = [1 Ly]
meD.h   = [1 meD.L(2)/nEy]                                     ;% [dx dy]
meD.h(1)= meD.h(2);
meD.L(1)= meD.h(1);

[xn,yn] = meshgrid(0.0:meD.h(1):meD.L(1),...
                             0.0:meD.h(2):1.0*meD.L(2)) ;%
xn      = flip(xn)                                                        ;%
yn      = flip(yn)                                                        ;%
               
               
meD.nNx  = size(xn,2)                                                      ;% number of nodes along x
meD.nNy  = size(yn,1)                                                      ;% number of nodes along y
meD.nN   = meD.nNx*meD.nNy                                                 ;% number of nodes
meD.nEx  = meD.nNx-1;
meD.nEy  = meD.nNy-1;
meD.nNe  = 16                                                               ;% number of node per element
meD.DoF  = 2                                                               ;% degree of freedom
meD.nDoF = meD.DoF.*[meD.nNe meD.nN]                                       ;% local and global number of degree of freedom
meD.x    = xn(:)                                                           ;% x coordinate
meD.y    = yn(:)                                                           ;% y coordinate
%% ELEMENT-NODE CONNECTIVITY
[meD.e2N] = e2N(meD.nNy,meD.nNx,meD.nEx,meD.nEy,4)                        ;% element to node topology
%-------------------------------------------------------------------------%

end

