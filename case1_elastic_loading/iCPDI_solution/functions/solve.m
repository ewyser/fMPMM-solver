function [uvw,duvw,fErr] = solve(K,oobf,uvw,duvw,bc,fd,NRit)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
K(bc(:,1),:)=0.0;
D_K         =spdiags(K,0);
D_K(bc(:,1))=1.0;
K           =spdiags(D_K,0,K);
if(NRit==0)
    oobf(bc(:,1))= bc(:,2);
else
    oobf(bc(:,1))= 0.0;
end
duvw(fd)=K(fd,fd)\oobf(fd)                                    ;
uvw     =uvw+duvw                                             ;% update displacements

%fErr = norm(oobf);;
fErr = max(abs(duvw))/(max(abs(uvw)));
end

