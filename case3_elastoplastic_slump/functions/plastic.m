function [mpD,D] = plastic(mpD,Del)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
c            = mpD.coh                                                    ;% hardening if Hp \neq 0 otherwise perfectly plastic material
%% MOHR-COULOMB YIELD FUNCTION
ds           = mpD.s(1,:)-mpD.s(2,:)                                      ;% preprocessing: sig_xx - sig_yy
tau          = sqrt(0.25*(ds).^2 + mpD.s(4,:).^2)                         ;% deviatoric stress
sig          = 0.5*(mpD.s(1,:)+mpD.s(2,:))                                ;% mohr-circle center: pressure
f            = tau+sig.*sin(mpD.phi)-c.*cos(mpD.phi)                      ;% yield function f<=0, if f>0 then return mapping
mpD.sn       = mpD.s                                                      ;% store stress 
%--------------------------------------------------------------------------%

%% PLASTIC CORRECTION: ONE STEP EXACT RETURN MAPPING
beta         = abs(c.*cos(mpD.phi)-sig.*sin(mpD.phi))./tau                ;%
dsigA        = 0.5*beta.*(ds)                                             ;%
dsigB        = c./tan(mpD.phi)                                            ;%
iDA          = (sig<=dsigB & f>0)                                         ;% logical indexing
mpD.sn(1,iDA)= sig(iDA)+dsigA(iDA)                                        ;%
mpD.sn(2,iDA)= sig(iDA)-dsigA(iDA)                                        ;%
mpD.sn(4,iDA)= beta(iDA).*mpD.s(4,iDA)                                    ;%
iDB           = (sig> dsigB & f>0)                                        ;% logical indexing
mpD.sn(1,iDB) = dsigB(iDB)                                                ;%
mpD.sn(2,iDB) = dsigB(iDB)                                                ;%
mpD.sn(4,iDB) = 0.0                                                       ;%
%--------------------------------------------------------------------------%

%% COMPUTE NEW STRESS, VOLUMETRIC, DEVIATORIC AND SECOND DEVIATORIC INVARIANT
dsig         = mpD.sn-mpD.s                                               ;%
mpD.s        = mpD.sn                                                     ;%
mpD.ep       = Del\dsig                                                   ;%
mpD.epV      = mpD.ep'*[1;1;1;0]                                          ;%
mpD.devep    = mpD.ep-repmat(mpD.epV',4,1).*repmat([1;1;1;0],1,mpD.n)     ;%
mpD.epII     = mpD.epII+sqrt(2/3*(mpD.ep(1,:).^2+mpD.ep(2,:).^2+...
               mpD.ep(3,:).^2+2.*mpD.ep(4,:).^2))                         ;%
%--------------------------------------------------------------------------%

%% SAVE           
D            = [mpD.epII c]                                               ;%
clear c ds tau sig f beta dsigA dsigB iDA iDB dsig                        ;% clear temporary variables
%--------------------------------------------------------------------------%

end

