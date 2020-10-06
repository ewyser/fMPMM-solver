function [mpD] = DefUpdate_loop(meD,mpD,p2n,p2e)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
uGIMP = 0; cpGIMP = 1;
%% UPDATE VORTICITY & STRAIN-RATE
v0 = zeros(meD.nEx*meD.nEy,1);
vt = v0;
I  = eye(meD.DoF,meD.DoF);
for mp = 1:mpD.n
    pie     = p2e(mp)                                                     ;% active element
    iDx     = meD.DoF*p2n(mp,:)-1                                         ;% x-component index
    iDy     = iDx+1                                                       ;% y-component index
    uv      = [meD.u(iDx)';meD.u(iDy)']                                   ;% displacement matrix   
    % INCREMENTAL DEFORMATION  
    dF      = I+[mpD.dSx(mp,:)*uv';mpD.dSy(mp,:)*uv']'                    ;% incremental deformation gradient 
    F       = dF*reshape(mpD.F(mp,:),2,2)'                                ;% deformation gradient
    v0(pie) = v0(pie)+mpD.V0(mp)                                          ;% element initial volume
    vt(pie) = vt(pie)+mpD.V0(mp)*det(F)                                   ;% element update volume
    %----------------------------------------------------------------------%
    %% UPDATE MP DEFORMATION GRADIENT
    mpD.F(mp,:) = reshape(F',4,1)                                         ;% -
    mpD.dF(mp,:)= reshape(dF',4,1)                                        ;% -
    %----------------------------------------------------------------------%
end
cellJ=vt./v0;
[r]=find(vt==0);
cellJ(r)=0.0;

for mp = 1:mpD.n
    pie         = p2e(mp)                                                 ;% active element
    dF          = reshape(mpD.dF(mp,:),2,2)'                              ;% incremental deformation
    F           = reshape(mpD.F(mp,:) ,2,2)'                              ;% deformation gradient
    scale       = (cellJ(pie)./det(F)).^(1/3)                             ;% -
    dF          = scale.*dF                                               ;% scaled incremental deformation gradient
    F           = scale.*F                                                ;% scaled deformation gradient
    eps         = 0.5.*(dF+dF')-I                                         ;% incremental small strains
    ome         = 0.5.*(dF-dF')                                           ;% incremental vorticity
    mpD.e(:,mp) = [eps(1,1);eps(2,2);2*eps(1,2)]                      ;% -
    mpD.w(:,3 ) = [ome(1,2)]                                              ;% -
    %----------------------------------------------------------------------%
    %% MATERIAL DOMAIN, VOLUME & POROSITY UPDATE
    mpD.J(mp) = det(F)                                                    ;% determinant of the deformation gradient
    mpD.V(mp) = mpD.J(mp).*mpD.V0(mp)                                     ;% updated volume
    if(uGIMP==1)
        mpD.l(mp,:) = mpD.l0(mp,:)                                        ;% no update method
    elseif(cpGIMP==1)
        mpD.l(mp,:) = mpD.J(mp).*mpD.l0(mp,:)                             ;% update based on the determinant of F
        mpD.l(mp,:) = mpD.l0(mp,:).*diag(sqrt(F'*F))'                     ;% update based on the diagonal components of U
        %mpD.l(mp,:) = diag(F)                                             ;% update based on the diagonal components of F
    end
    mpD.F(mp,:) = reshape(F',4,1)                                         ;% -
end

end

