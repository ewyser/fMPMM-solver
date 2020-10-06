function [K,fint] = getStiffness(meD,mpD,sig,vp,c2N,Del)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% UPDATE MP'S VOLUME AND B MATRIX IN THE UPDATED FRAME
vp   = vp.*mpD.J;
Bp   = zeros(4,meD.nDoF(1),mpD.n);
for mp = 1:mpD.n
    [B]   = getB(mpD.dSx(mp,:),mpD.dSy(mp,:),mpD.dF(mp,:));
    Bp(:,:,mp) = B;
end
%% INTERNAL FORCE CONTRIBUTION
% INITIALIZATION
fint = zeros(meD.nDoF(2),1)                                                ;%
% CONTRIBUTION FROM p TO N
fi= squeeze(sum(Bp.*repmat(reshape(sig,size(sig,1),1,mpD.n),1,meD.nDoF(1)),1)).*repmat(vp',32,1); % 
for n = 1:meD.nNe                                                          % BEGIN ITERATION OVER meD.nNe NEIGHBORING NODES FOR ALL MP
    l = [(meD.DoF*c2N(:,n)-1);(meD.DoF*c2N(:,n))]                         ;% local to global index
    fint = fint + accumarray(l,[fi(n*meD.DoF-1,:)';...
                                fi(n*meD.DoF  ,:)'],[meD.nDoF(2) 1])       ;% external force global vector - nodal internal force global vector
end
%% MATERIAL POINT STIFFNESS MATRIX
kprow = zeros(mpD.n*meD.nDoF(1)^2,1);
kpcol = kprow;
kpval = kpcol;
for mp = 1:mpD.n
    no    = reshape([meD.DoF.*c2N(mp,:)-1;meD.DoF.*c2N(mp,:)],meD.nDoF(1),1)';
    kp    = vp(mp)*(Bp(:,:,mp)'*Del*Bp(:,:,mp));
    iD          = (mp*length(no)^2+1)-length(no)^2:1:mp*length(no)^2;
    kprow(iD,1) = repmat(no,1,length(no))';
    kpcol(iD,1) = reshape(repmat(no,length(no),1),length(no)^2,1);
    kpval(iD,1) = kp(:);
end
%% GLOBAL STIFFNESS MATRIX ASSEMBLY
K = sparse(kprow,kpcol,kpval,meD.nDoF(2),meD.nDoF(2));

end
function [B] = getB(dNx,dNy,dF)
%% INITIALIZATION
B  = zeros(4,length(dNx)*2);
dF = reshape(dF,2,2)';
dN = dF\[dNx;dNy];
%% CALCULATE B MATRIX
B([1;4],1:2:end)=dN;
B([4;2],2:2:end)=dN;
end
