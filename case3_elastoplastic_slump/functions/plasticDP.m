function [mpD,D] = plasticDP(mpD,Gc,Kc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tensile  = 27.48e3;
c        = mpD.coh;
psi      = 0.0*pi/180;
mpD.P    = (mpD.s'*[1;1;1;0])./3;
mpD.devs = mpD.s-(mpD.P'.*[1;1;1;0]);
mpD.J2   = 0.5.*(mpD.devs(1,:).^2+mpD.devs(2,:).^2+mpD.devs(3,:).^2)+mpD.devs(4,:).^2;
tau      = sqrt(mpD.J2)';

qphi     = ((6.*sin(mpD.phi))./(sqrt(3).*(3+sin(mpD.phi))))';
qpsi     = ((6.*sin(psi))./(sqrt(3).*(3+sin(psi))))';
kphi     = ((6.*c.*cos(mpD.phi))./(sqrt(3).*(3+sin(mpD.phi))))';

sigmat   = min(tensile,(kphi./qphi));
sigmat   = tensile;
fs    = tau+qphi.*mpD.P-kphi; % yield function considering shear failure
ft    = mpD.P-sigmat;% yield function considering tensile failure
tauP  = kphi-qphi.*sigmat;
alpP  = sqrt(1+qphi.^2)-qphi;
h     = tau-tauP-alpP.*(mpD.P-sigmat);

IDX = ((fs>0)&(mpD.P<sigmat))|((h>0.0)&(mpD.P>=sigmat));
dlam= fs(IDX)./(Gc+Kc.*qphi(IDX).*qpsi);
Pn  = mpD.P(IDX)-Kc.*qpsi.*dlam;
taun= kphi(IDX)-qphi(IDX).*Pn;
mpD.s(:,IDX)= mpD.devs(:,IDX).*repmat(taun./tau(IDX),1,4)'+(Pn*[1 1 1 0])';
mpD.epII(IDX)= mpD.epII(IDX)+(dlam.*sqrt(1/3+(2/9).*qpsi.^2))';

IDX = ((h<=0.0)&(mpD.P>=sigmat));
dlam = (mpD.P(IDX)-sigmat)./Kc;
mpD.s(:,IDX) = mpD.s(:,IDX) + ((sigmat-mpD.P(IDX))*[1 1 1 0])';
mpD.epII(IDX) = mpD.epII(IDX) + (sqrt(2).*dlam./3)';

%% SAVE           
D            = [mpD.epII]                                                 ;%
clear c ds tau sig f beta dsigA dsigB iDA iDB dsig                        ;% clear temporary variables
%--------------------------------------------------------------------------%

end

