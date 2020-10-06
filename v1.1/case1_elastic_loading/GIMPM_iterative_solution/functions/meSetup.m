function [meD,bc] = meSetup(nEx,typeD,l0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% MESH INITIALIZATION
meD.vd  = 0.0                                                            ;% viscous damping coefficient
Lx      = l0                                                         ;% mesh dimension (x,y)
Ly      = l0                                                              ;% number of element in x direction
nEy     = nEx                                                             ;% number of element in y direction
meD.L   = [Lx Ly]                                                         ;% mesh length in x            [m]
meD.h   = [meD.L(2)/nEy meD.L(2)/nEy]                                     ;% [dx dy]
Lx      = meD.h(2)
meD.L   = [Lx Ly]                                                         ;% mesh length in x  
[xn,yn] = meshgrid(0.0-meD.h(1):meD.h(1):meD.L(1)+meD.h(1),0.0-meD.h(2):meD.h(2):1.0*meD.L(2)+meD.h(2)) ;%
xn      = flip(xn)                                                        ;%
yn      = flip(yn)                                                        ;%
               
               
meD.nNx = size(xn,2)                                                      ;% number of nodes along x
meD.nNy = size(yn,1)                                                      ;% number of nodes along y
meD.nN  = meD.nNx*meD.nNy                                                 ;% number of nodes
meD.nEx = meD.nNx-1;
meD.nEy = meD.nNy-1;
meD.nNe = 16                                                              ;% number of node per element
meD.DoF = 2                                                               ;% degree of freedom
meD.nDoF = meD.DoF.*[meD.nNe meD.nN]                                       ;% local and global number of degree of freedom
meD.x   = xn(:)                                                           ;% x coordinate
meD.y   = yn(:)                                                           ;% y coordinate
% NODAL VECTOR INITIALIZATION                                              
meD.m   = zeros(meD.nN    ,1,typeD)                                       ;% mass
meD.mr  = zeros(meD.nDoF(2),1,typeD)                                       ;% repmated mass
meD.f   = zeros(meD.nDoF(2),1,typeD)                                       ;% force balance                                                       
meD.fi  = zeros(meD.nDoF(2),1,typeD)                                       ;% internal force
meD.d   = zeros(meD.nDoF(2),1,typeD)                                       ;% damping force
meD.a   = zeros(meD.nDoF(2),1,typeD)                                       ;% acceleration
meD.p   = zeros(meD.nDoF(2),1,typeD)                                       ;% momentum
meD.v   = zeros(meD.nDoF(2),1,typeD)                                       ;% velocity
meD.u   = zeros(meD.nDoF(2),1,typeD)                                       ;% displacement
%-------------------------------------------------------------------------%

%% ELEMENT-NODE CONNECTIVITY
[meD.e2N] = e2N(meD.nNy,meD.nNx,meD.nEx,meD.nEy,meD.nNe)                        ;% element to node topology
%-------------------------------------------------------------------------%

% BOUNDARY CONDITIONS
yinf   = 0.0                                                              ;%
meD.xB = [-meD.L(1)/2 meD.L(1)/2]                                         ;%
meD.xB = [0.0 meD.L(1)]                                         ;%
[row]  = find(meD.y<=yinf)                                                ;%
bc.y   = row                                                              ;%
[row]  = find(meD.x<=meD.xB(1))                                           ;%
bc.x   = row                                                              ;%
[row]  = find(meD.x>=meD.xB(2))                                           ;%
bc.x   = [bc.x;row]                                                       ;%

bc.y   = meD.DoF*bc.y(:,1);
bc.y   = [bc.y zeros(size(bc.y))];

bc.x   = meD.DoF*bc.x(:,1)-1 ;  
bc.x   = [bc.x zeros(size(bc.x))];

end

