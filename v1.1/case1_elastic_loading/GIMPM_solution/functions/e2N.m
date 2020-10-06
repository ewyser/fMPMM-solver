function [g_num] = e2N(nny,nnx,nelx,nely,nnpe)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
gnumbers= flip(reshape(1:(nnx*nny),nny ,nnx ))                            ;%
g_num   = zeros(nelx*nely,nnpe,'int32')                                   ;%
iel     = 1                                                               ;%
if(nnpe==4)
    disp('------------------------')                                      ;%
    disp('4 nodes per quadrilaterial element')                            ;%
    disp('------------------------')                                      ;%
    for i = 1:nelx
        for j = 1:nely
            if(i>1 && i<nelx && j>1 && j<nely)
                g_num(iel,1) = gnumbers(j  ,i  )                          ;%
                g_num(iel,2) = gnumbers(j+1,i  )                          ;%
                g_num(iel,3) = gnumbers(j  ,i+1)                          ;%
                g_num(iel,4) = gnumbers(j+1,i+1)                          ;%
            end
            iel = iel+1;
        end
    end
elseif(nnpe==16)
    disp('------------------------')                                      ;%
    disp('16 nodes per quadrilaterial element')                           ;%
    disp('------------------------')                                      ;%
    for i = 1:nelx
        for j = 1:nely
            if(i>1 && i<nelx && j>1 && j<nely)
                g_num(iel,1 ) = gnumbers(j-1,i-1)                         ;%
                g_num(iel,2 ) = gnumbers(j-0,i-1)                         ;%
                g_num(iel,3 ) = gnumbers(j+1,i-1)                         ;%
                g_num(iel,4 ) = gnumbers(j+2,i-1)                         ;%
                
                g_num(iel,5 ) = gnumbers(j-1,i  )                         ;%
                g_num(iel,6 ) = gnumbers(j-0,i  )                         ;%
                g_num(iel,7 ) = gnumbers(j+1,i  )                         ;%
                g_num(iel,8 ) = gnumbers(j+2,i  )                         ;%
                
                g_num(iel,9 ) = gnumbers(j-1,i+1)                         ;%
                g_num(iel,10) = gnumbers(j-0,i+1)                         ;%
                g_num(iel,11) = gnumbers(j+1,i+1)                         ;%
                g_num(iel,12) = gnumbers(j+2,i+1)                         ;%
                
                g_num(iel,13) = gnumbers(j-1,i+2)                         ;%
                g_num(iel,14) = gnumbers(j-0,i+2)                         ;%
                g_num(iel,15) = gnumbers(j+1,i+2)                         ;%
                g_num(iel,16) = gnumbers(j+2,i+2)                         ;%
            end
            iel = iel+1;
        end
    end
elseif(nnpe==36)
    disp('------------------------')                                      ;%
    disp('36 nodes per quadrilaterial element')                           ;%
    disp('------------------------')                                      ;%
    for i = 1:nelx
        for j = 1:nely
            if(i>2 && i<nelx-1 && j>2 && j<nely-1)
                g_num(iel,1) = gnumbers(j-2,i-2)                          ;%
                g_num(iel,2) = gnumbers(j-1,i-2)                          ;%
                g_num(iel,3) = gnumbers(j  ,i-2)                          ;%
                g_num(iel,4) = gnumbers(j+1,i-2)                          ;%
                g_num(iel,5) = gnumbers(j+2,i-2)                          ;%
                g_num(iel,6) = gnumbers(j+3,i-2)                          ;%
                
                g_num(iel,7) = gnumbers(j-2,i-1)                          ;%
                g_num(iel,8) = gnumbers(j-1,i-1)                          ;%
                g_num(iel,9) = gnumbers(j-0,i-1)                          ;%
                g_num(iel,10)= gnumbers(j+1,i-1)                          ;%
                g_num(iel,11)= gnumbers(j+2,i-1)                          ;%
                g_num(iel,12)= gnumbers(j+3,i-1)                          ;%
                
                g_num(iel,13)= gnumbers(j-2,i  )                          ;%
                g_num(iel,14)= gnumbers(j-1,i  )                          ;%
                g_num(iel,15)= gnumbers(j-0,i  )                          ;%
                g_num(iel,16)= gnumbers(j+1,i  )                          ;%
                g_num(iel,17)= gnumbers(j+2,i  )                          ;%
                g_num(iel,18)= gnumbers(j+3,i  )                          ;%
                
                g_num(iel,19)= gnumbers(j-2,i+1)                          ;%
                g_num(iel,20)= gnumbers(j-1,i+1)                          ;%
                g_num(iel,21)= gnumbers(j-0,i+1)                          ;%
                g_num(iel,22)= gnumbers(j+1,i+1)                          ;%
                g_num(iel,23)= gnumbers(j+2,i+1)                          ;%
                g_num(iel,24)= gnumbers(j+3,i+1)                          ;%
                
                g_num(iel,25)= gnumbers(j-2,i+2)                          ;%
                g_num(iel,26)= gnumbers(j-1,i+2)                          ;%
                g_num(iel,27)= gnumbers(j-0,i+2)                          ;%
                g_num(iel,28)= gnumbers(j+1,i+2)                          ;%
                g_num(iel,29)= gnumbers(j+2,i+2)                          ;%
                g_num(iel,30)= gnumbers(j+3,i+2)                          ;%
                
                g_num(iel,31)= gnumbers(j-2,i+3)                          ;%
                g_num(iel,32)= gnumbers(j-1,i+3)                          ;%
                g_num(iel,33)= gnumbers(j-0,i+3)                          ;%
                g_num(iel,34)= gnumbers(j+1,i+3)                          ;%
                g_num(iel,35)= gnumbers(j+2,i+3)                          ;%
                g_num(iel,36)= gnumbers(j+3,i+3)                          ;%
            end
            iel = iel+1;
        end
    end
else
end

end

