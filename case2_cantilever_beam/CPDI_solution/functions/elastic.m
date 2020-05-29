function [mpD] = elastic(mpD,Del,lame1,lame2)
%f_elastic : function to calculate 1) a rotated stress to ensure objectivity
%of the stress tensor and, 2) an elastic predictor
%   Detailed explanation goes here
%% PREPROCESSING
omega      = -mpD.w(:,3)'                                                 ;% incremental rotation in (x,y) plane
cos2       = cos(omega).^2                                                ;% -
sin2       = sin(omega).^2                                                ;% -
sincos     = sin(omega).*cos(omega)                                       ;% -
%--------------------------------------------------------------------------%

%% OBJECTIVE STRESS ROTATION DUE TO NON-OBJECTIVE STESS RATE ds_{ij}/dt
mpD.s(1,:) = mpD.s(1,:).*cos2 + ...
             mpD.s(2,:).*sin2 - 2.*mpD.s(4,:).*sincos                     ;% primed sig_xx
mpD.s(2,:) = mpD.s(1,:).*sin2 + ...
             mpD.s(2,:).*cos2 + 2.*mpD.s(4,:).*sincos                     ;% primed sig_yy
mpD.s(4,:) =(mpD.s(1,:)-mpD.s(2,:)).*sincos + ...
             mpD.s(4,:).*(cos2-sin2)                                      ;% primed sig_xy/sig_yx       
%--------------------------------------------------------------------------%

%% STRESS-STRAIN CONSTITUTIVE RELATION
mpD.s      = mpD.s+(Del*mpD.e)                                            ;% elastic constitutive relation

% P     = (lame2.*log(mpD.J))./mpD.J;
% dev   = (lame1./mpD.J).*(mpD.b-[1 1 0 0]);
% mpD.s = [P P P zeros(size(P))]'+dev';

clear omega cos2 sin2 sincos                                              ;% clear temporary variables  
%--------------------------------------------------------------------------%

end

