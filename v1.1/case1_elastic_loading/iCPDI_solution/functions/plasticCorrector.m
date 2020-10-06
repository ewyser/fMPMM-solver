function [sigN,mpD] = plasticCorrector(sig0,mpD,Hp,cohr,Del)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
c            = mpD.coh+(Hp*mpD.epII ); c(c<cohr)=cohr                     ;% hardening if Hp \neq 0 otherwise perfectly plastic material
%c            = mpD.coh                                                    ;
%% MOHR-COULOMB YIELD FUNCTION
ds           = sig0(1,:)-sig0(2,:)                                      ;% preprocessing: sig_xx - sig_yy
tau          = sqrt(0.25*(ds).^2 + sig0(4,:).^2)                         ;% deviatoric stress
sig          = 0.5*(sig0(1,:)+sig0(2,:))                                ;% mohr-circle center: pressure
f            = tau+sig.*sin(mpD.phi)-c.*cos(mpD.phi)                      ;% yield function f<=0, if f>0 then return mapping
sigN         = sig0                                                      ;% store stress 
%--------------------------------------------------------------------------%

%% PLASTIC CORRECTION: ONE STEP EXACT RETURN MAPPING
beta         = abs(c.*cos(mpD.phi)-sig.*sin(mpD.phi))./tau                ;%
dsigA        = 0.5*beta.*(ds)                                             ;%
dsigB        = c./tan(mpD.phi)                                            ;%
iDA          = (sig<=dsigB & f>0)                                         ;% logical indexing
sigN(1,iDA)   = sig(iDA)+dsigA(iDA)                                        ;%
sigN(2,iDA)   = sig(iDA)-dsigA(iDA)                                        ;%
sigN(4,iDA)   = beta(iDA).*sig0(4,iDA)                                    ;%
iDB           = (sig> dsigB & f>0)                                        ;% logical indexing
sigN(1,iDB) = dsigB(iDB)                                                ;%
sigN(2,iDB) = dsigB(iDB)                                                ;%
sigN(4,iDB) = 0.0                                                       ;%
%--------------------------------------------------------------------------%

%% COMPUTE NEW STRESS, VOLUMETRIC, DEVIATORIC AND SECOND DEVIATORIC INVARIANT
dsig         = sigN-sig0                                                  ;%
ep           = Del\dsig                                                   ;%
epV          = ep'*[1;1;1;0]                                              ;%
devep        = ep-repmat(epV',4,1).*repmat([1;1;1;0],1,mpD.n)             ;% 
mpD.depII    = sqrt(2/3*(ep(1,:).^2+ep(2,:).^2+...
               ep(3,:).^2+2.*ep(4,:).^2))                                 ;%
%--------------------------------------------------------------------------%

end

