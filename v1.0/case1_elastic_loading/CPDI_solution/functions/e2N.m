function [g_num] = e2N(nny,nnx,nelx,nely,nnpe)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
gnumbers= flip(reshape(1:(nnx*nny),nny ,nnx ))                            ;%
g_num   = zeros(nelx*nely,4)                                              ;%
iel     = 1                                                               ;%
disp('------------------------')                                          ;%
disp('4 nodes per quadrilateral element')                                 ;%
disp('------------------------')                                          ;%
for i = 1:nelx
    for j = 1:nely
        if(i>1 && i<nelx && j>1 && j<nely)
            g_num(iel,1) = gnumbers(j  ,i  )                              ;%
            g_num(iel,2) = gnumbers(j+1,i  )                              ;%
            g_num(iel,3) = gnumbers(j  ,i+1)                              ;%
            g_num(iel,4) = gnumbers(j+1,i+1)                              ;%
        end
        iel = iel+1;
    end
end


end

